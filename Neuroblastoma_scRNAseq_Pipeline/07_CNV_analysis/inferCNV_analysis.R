#!/usr/bin/env Rscript
## =============================================================================
## inferCNV Analysis for NB Patient Samples
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Detect copy number variations using inferCNV with cluster 10 as reference
## Input:   Seurat object with clustering and annotations
## Output:  CNV scores, heatmaps, visualizations
## =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(infercnv)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
})

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")
REFERENCES_DIR <- file.path(RESULTS_DIR, "references")
OUTPUT_DIR <- file.path(RESULTS_DIR, "inferCNV_cluster10_NBAtlas_method")

SEURAT_FILE <- file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds")
GTF_FILE <- file.path(REFERENCES_DIR, "gencode.v43.annotation.gtf.gz")

# Reference cluster (immune cells for CNV normalization)
REFERENCE_CLUSTERS <- "10"

# inferCNV parameters (matching NBAtlas)
CUTOFF <- 0.1
DENOISE <- TRUE
HMM <- FALSE
N_THREADS <- 4
UMAP_REDUCTION <- "umap"
CHROMOSOME_ARMS <- "chr7,chr17q,chr1p,chr2p,chr3p,chr4p,chr11q,chr14q"

# Create output directories
out_dir <- normalizePath(OUTPUT_DIR, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

message("============================================================")
message("inferCNV (NBAtlas Method) Pipeline Starting")
message("============================================================")
message("Seurat object:      ", SEURAT_FILE)
message("Reference clusters: ", REFERENCE_CLUSTERS)
message("GTF annotation:     ", GTF_FILE)
message("Output directory:   ", OUTPUT_DIR)
message("============================================================")

ref_clusters <- trimws(strsplit(REFERENCE_CLUSTERS, ",")[[1]])
chr_arms_to_analyze <- trimws(strsplit(CHROMOSOME_ARMS, ",")[[1]])

## =============================================================================
## STEP 1: Load Seurat object
## =============================================================================

message("\n[1/9] Loading Seurat object...")

nb <- readRDS(SEURAT_FILE)

if (!inherits(nb, "Seurat")) {
  stop("Loaded object is not a Seurat object.")
}

if (!(UMAP_REDUCTION %in% names(nb@reductions))) {
  stop("UMAP reduction '", UMAP_REDUCTION, "' not found in Seurat object.")
}

cluster_ids <- as.character(nb$seurat_clusters)

message("  Cells per cluster (first 10):")
print(head(sort(table(cluster_ids), decreasing = TRUE), 10))

## =============================================================================
## STEP 2: Prepare reference annotations
## =============================================================================

message("\n[2/9] Preparing reference annotations...")

annotation_df <- data.frame(
  cell = colnames(nb),
  cluster = cluster_ids,
  infercnv_group = ifelse(cluster_ids %in% ref_clusters, "Reference", "Observed"),
  stringsAsFactors = FALSE
)

write.table(
  annotation_df[, c("cell", "cluster")],
  file = file.path(out_dir, "cell_annotations.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("  Reference cells: ", sum(annotation_df$infercnv_group == "Reference"))
message("  Observed cells:  ", sum(annotation_df$infercnv_group == "Observed"))

## =============================================================================
## STEP 3: Extract raw counts
## =============================================================================

message("\n[3/9] Extracting raw counts...")

DefaultAssay(nb) <- "RNA"

if (packageVersion("SeuratObject") >= "5.0.0") {
  counts <- GetAssayData(nb, layer = "counts")
} else {
  counts <- GetAssayData(nb, slot = "counts")
}

if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

message("  Matrix dimensions: ", nrow(counts), " genes x ", ncol(counts), " cells")

## =============================================================================
## STEP 4: Build gene ordering table from GTF
## =============================================================================

message("\n[4/9] Building gene ordering table from GTF...")

gtf_dt <- data.table::fread(
  GTF_FILE,
  sep = "\t",
  header = FALSE,
  skip = "#",
  col.names = c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
  )
)

genes_gtf <- gtf_dt[feature == "gene", .(
  gene = sub(".*gene_name \"([^\"]+)\".*", "\\1", attribute),
  chr = seqname,
  start = start,
  end = end
)]

genes_gtf$chr <- gsub("^chr", "", genes_gtf$chr)
genes_gtf <- genes_gtf[gene %in% rownames(counts)]
genes_gtf <- genes_gtf[!grepl("^(M|MT)$", chr)]

genes_gtf$chr <- factor(genes_gtf$chr, levels = c(1:22, "X", "Y"))
genes_gtf <- genes_gtf[order(chr, start)]
genes_gtf <- genes_gtf[!duplicated(gene)]

if (nrow(genes_gtf) == 0) {
  stop("No genes from the GTF matched the Seurat object.")
}

genes_gtf$chr <- paste0("chr", genes_gtf$chr)

gene_order_path <- file.path(out_dir, "gene_ordering.tsv")

write.table(
  genes_gtf,
  file = gene_order_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("  Genes included: ", nrow(genes_gtf))

## =============================================================================
## STEP 5: Create inferCNV object
## =============================================================================

message("\n[5/9] Creating inferCNV object...")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = file.path(out_dir, "cell_annotations.tsv"),
  delim = "\t",
  gene_order_file = gene_order_path,
  ref_group_names = ref_clusters
)

## =============================================================================
## STEP 6: Run inferCNV
## =============================================================================

message("\n[6/9] Running inferCNV...")

options(scipen = 100)

n_cells <- ncol(infercnv_obj@expr.data)
n_genes <- nrow(infercnv_obj@expr.data)

message("  Dataset size: ", n_cells, " cells x ", n_genes, " genes")

enable_plotting <- n_cells < 50000

if (!enable_plotting) {
  message("  WARNING: Dataset too large (", n_cells, " cells) - disabling heatmap generation")
  message("  CNV scores will still be calculated")
}

infercnv_res <- infercnv::run(
  infercnv_obj,
  cutoff = CUTOFF,
  out_dir = out_dir,
  cluster_by_groups = TRUE,
  denoise = DENOISE,
  HMM = HMM,
  num_threads = N_THREADS,
  plot_steps = enable_plotting,
  scale_data = TRUE
)

## =============================================================================
## STEP 7: Calculate Overall CNV Scores (NBAtlas Method)
## =============================================================================

message("\n[7/9] Calculating overall CNV scores (NBAtlas method)...")

expr_matrix <- infercnv_res@expr.data
cnv_scores <- colMeans(abs(expr_matrix))

cnv_scores_df <- data.frame(
  cell = names(cnv_scores),
  cnv_score = as.numeric(cnv_scores),
  stringsAsFactors = FALSE
)

cnv_scores_df <- cnv_scores_df %>%
  left_join(annotation_df, by = "cell")

write.table(
  cnv_scores_df,
  file = file.path(table_dir, "cnv_scores_per_cell.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cluster_summary <- cnv_scores_df %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_cnv = mean(cnv_score),
    median_cnv = median(cnv_score),
    sd_cnv = sd(cnv_score),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cnv))

write.table(
  cluster_summary,
  file = file.path(table_dir, "cnv_scores_by_cluster.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("  Mean CNV score across all cells: ", round(mean(cnv_scores), 4))

## =============================================================================
## STEP 8: Calculate Chromosome-Specific CNV Scores (NBAtlas Method)
## =============================================================================

message("\n[8/9] Calculating chromosome-specific CNV scores...")

gene_info <- genes_gtf
gene_info$gene_idx <- match(gene_info$gene, rownames(expr_matrix))
gene_info <- gene_info[!is.na(gene_info$gene_idx), ]

calculate_chr_score <- function(chr_name, expr_matrix, gene_info) {

  if (!grepl("[pq]$", chr_name)) {
    gene_indices <- gene_info$gene_idx[gene_info$chr == chr_name]
  } else {
    chr_base <- gsub("[pq]$", "", chr_name)
    arm <- substr(chr_name, nchar(chr_name), nchar(chr_name))

    chr_genes <- gene_info[gene_info$chr == chr_base, ]
    if (nrow(chr_genes) == 0) return(NULL)

    centromere_pos <- median(chr_genes$start)

    if (arm == "p") {
      gene_indices <- chr_genes$gene_idx[chr_genes$start < centromere_pos]
    } else {
      gene_indices <- chr_genes$gene_idx[chr_genes$start >= centromere_pos]
    }
  }

  if (length(gene_indices) < 10) {
    warning("  Skipping ", chr_name, ": fewer than 10 genes found")
    return(NULL)
  }

  chr_expr <- expr_matrix[gene_indices, , drop = FALSE]
  chr_scores <- colMeans(abs(chr_expr))

  message("  ", chr_name, ": ", length(gene_indices), " genes, mean score = ",
          round(mean(chr_scores), 4))

  return(chr_scores)
}

chr_scores_list <- list()

for (chr_arm in chr_arms_to_analyze) {
  scores <- calculate_chr_score(chr_arm, expr_matrix, gene_info)
  if (!is.null(scores)) {
    chr_scores_list[[chr_arm]] <- scores
  }
}

chr_scores_combined <- cnv_scores_df[, c("cell", "cluster", "infercnv_group")]

for (chr_arm in names(chr_scores_list)) {
  col_name <- paste0(chr_arm, "_score")
  chr_scores_combined[[col_name]] <- chr_scores_list[[chr_arm]][chr_scores_combined$cell]
}

write.table(
  chr_scores_combined,
  file = file.path(table_dir, "cnv_scores_per_chromosome.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## =============================================================================
## STEP 9: Add CNV Scores to Seurat Object & Generate UMAP Plots
## =============================================================================

message("\n[9/9] Adding CNV scores to Seurat object and generating UMAP plots...")

nb$infercnv_score <- NA
nb@meta.data[cnv_scores_df$cell, "infercnv_score"] <- cnv_scores_df$cnv_score

for (chr_arm in names(chr_scores_list)) {
  col_name <- paste0(chr_arm, "_score")
  nb@meta.data[[col_name]] <- NA
  nb@meta.data[names(chr_scores_list[[chr_arm]]), col_name] <- chr_scores_list[[chr_arm]]
}

message("\n  Generating UMAP plots (NBAtlas style)...")

# Overall CNV Score - RdBu color scheme
p1 <- FeaturePlot(
  nb,
  reduction = UMAP_REDUCTION,
  features = "infercnv_score",
  order = TRUE,
  min.cutoff = "q2",
  max.cutoff = "q98",
  raster = FALSE
) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle("InferCNV Score (Overall)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(
  filename = file.path(fig_dir, "InferCNV_score_UMAP_RdBu_HD.png"),
  plot = p1,
  width = 12,
  height = 10,
  dpi = 300
)

# Overall CNV Score - PuOr color scheme
p2 <- FeaturePlot(
  nb,
  reduction = UMAP_REDUCTION,
  features = "infercnv_score",
  order = TRUE,
  min.cutoff = "q2",
  max.cutoff = "q98",
  raster = FALSE
) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr"))) +
  ggtitle("InferCNV Score (Overall)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(
  filename = file.path(fig_dir, "InferCNV_score_UMAP_PuOr_HD.png"),
  plot = p2,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(fig_dir, "InferCNV_score_UMAP_PuOr.png"),
  plot = p2,
  width = 9,
  height = 7,
  dpi = 300
)

# Chromosome-specific scores
for (chr_arm in names(chr_scores_list)) {
  col_name <- paste0(chr_arm, "_score")

  p <- FeaturePlot(
    nb,
    reduction = UMAP_REDUCTION,
    features = col_name,
    order = TRUE,
    min.cutoff = "q2",
    max.cutoff = "q98",
    raster = FALSE
  ) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr"))) +
    ggtitle(paste0("CNV Score: ", chr_arm)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    filename = file.path(fig_dir, paste0("CNV_", chr_arm, "_UMAP_PuOr_HD.png")),
    plot = p,
    width = 12,
    height = 10,
    dpi = 300
  )

  ggsave(
    filename = file.path(fig_dir, paste0("CNV_", chr_arm, "_UMAP_PuOr.png")),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )

  message("  Generated: CNV_", chr_arm, "_UMAP plots")
}

# Cluster identity plot
p_clusters <- DimPlot(
  nb,
  reduction = UMAP_REDUCTION,
  group.by = "seurat_clusters",
  label = TRUE,
  raster = FALSE
) +
  ggtitle("Clusters (Reference cells highlighted)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(
  filename = file.path(fig_dir, "Clusters_UMAP.png"),
  plot = p_clusters,
  width = 10,
  height = 8,
  dpi = 300
)

# Save updated Seurat object with CNV scores
updated_seurat_path <- file.path(out_dir, "seurat_with_cnv_scores.rds")
saveRDS(nb, file = updated_seurat_path)

message("  Saved updated Seurat object: ", updated_seurat_path)

## =============================================================================
## Summary Report
## =============================================================================

message("\n============================================================")
message("inferCNV Pipeline Completed (NBAtlas Method)")
message("============================================================")
message("Results directory:  ", out_dir)
message("")
message("Output files:")
message("  - CNV scores:     ", file.path(table_dir, "cnv_scores_per_cell.tsv"))
message("  - Chr scores:     ", file.path(table_dir, "cnv_scores_per_chromosome.tsv"))
message("  - Cluster summary:", file.path(table_dir, "cnv_scores_by_cluster.tsv"))
message("  - UMAP plots:     ", fig_dir)
message("  - Updated Seurat: ", updated_seurat_path)
message("")
message("Chromosome arms analyzed:")
for (chr_arm in names(chr_scores_list)) {
  message("  - ", chr_arm)
}
message("")
message("Top 5 clusters by CNV score:")
print(head(cluster_summary[, c("cluster", "mean_cnv", "n_cells")], 5))
message("============================================================")
