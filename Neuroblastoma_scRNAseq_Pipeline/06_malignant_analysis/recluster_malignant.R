#!/usr/bin/env Rscript
## =============================================================================
## Re-cluster Malignant Cells
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Subset to malignant clusters and re-cluster for finer resolution
## Input:   Seurat object with initial clustering (jansky annotations)
## Output:  Malignant-only Seurat object with new clusters, markers, signatures
## =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readxl)
library(readr)

set.seed(12345)

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")
OUTPUT_DIR <- file.path(RESULTS_DIR, "malignant_reclustering")
REFERENCES_DIR <- file.path(RESULTS_DIR, "references")

# Input: Use corrected Seurat object
SEURAT_FILE <- file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds")

# Malignant clusters to keep
MALIGNANT_CLUSTERS <- c(0, 1, 2, 3, 4, 5, 6, 12, 16)

# Re-clustering parameters
N_PCS <- 30
CLUSTER_RESOLUTIONS <- c(0.3, 0.5, 0.8, 1.0)
DEFAULT_RESOLUTION <- 0.5

# Gene signature files
VANGRONINGEN_FILE <- file.path(REFERENCES_DIR, "vanGroningen_2017.xlsx")
DYER_FILE <- file.path(REFERENCES_DIR, "Dyer_Sig.xlsx")
INT_SIGNATURE_CSV <- file.path(REFERENCES_DIR, "Int_signature.csv")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

message("============================================================")
message("RE-CLUSTERING MALIGNANT CELLS")
message("============================================================")
message("Input: ", SEURAT_FILE)
message("Output: ", OUTPUT_DIR)
message("Malignant clusters: ", paste(MALIGNANT_CLUSTERS, collapse = ", "))
message("============================================================\n")

## =============================================================================
## STEP 1: LOAD OBJECT
## =============================================================================

message("[1/6] Loading Seurat object...")

if (!file.exists(SEURAT_FILE)) {
  stop("ERROR: Seurat file not found: ", SEURAT_FILE)
}

nb <- readRDS(SEURAT_FILE)
message("  Loaded: ", ncol(nb), " cells x ", nrow(nb), " genes")

if (!"seurat_clusters" %in% colnames(nb@meta.data)) {
  stop("ERROR: 'seurat_clusters' column not found in metadata")
}

message("\n  Original cluster distribution:")
cluster_table <- table(nb$seurat_clusters)
for (clust in names(sort(cluster_table, decreasing = TRUE))) {
  n_cells <- cluster_table[clust]
  is_malignant <- as.numeric(clust) %in% MALIGNANT_CLUSTERS
  marker <- if (is_malignant) " [MALIGNANT]" else ""
  message("    Cluster ", clust, ": ", n_cells, " cells", marker)
}

## =============================================================================
## STEP 2: SUBSET TO MALIGNANT CLUSTERS
## =============================================================================

message("\n[2/6] Subsetting to malignant clusters...")

nb$original_cluster <- nb$seurat_clusters
cluster_numeric <- as.numeric(as.character(nb$seurat_clusters))
malignant_cells <- cluster_numeric %in% MALIGNANT_CLUSTERS
nb_malignant <- subset(nb, cells = colnames(nb)[malignant_cells])

n_original <- ncol(nb)
n_malignant <- ncol(nb_malignant)
pct_malignant <- round(100 * n_malignant / n_original, 1)

message("  Original cells: ", n_original)
message("  Malignant cells: ", n_malignant, " (", pct_malignant, "%)")
message("  Excluded cells: ", n_original - n_malignant)

message("\n  Malignant subset cluster distribution:")
cluster_table_mal <- table(nb_malignant$original_cluster)
for (clust in names(sort(cluster_table_mal, decreasing = TRUE))) {
  n_cells <- cluster_table_mal[clust]
  message("    Cluster ", clust, ": ", n_cells, " cells")
}

## =============================================================================
## STEP 3: RE-RUN DIMENSIONALITY REDUCTION
## =============================================================================

message("\n[3/6] Re-running dimensionality reduction on malignant cells...")

DefaultAssay(nb_malignant) <- "RNA"
message("  Switching to RNA assay for fresh normalization on subset...")

message("  Normalizing data...")
nb_malignant <- NormalizeData(nb_malignant, verbose = FALSE)

message("  Finding variable features for malignant subset...")
nb_malignant <- FindVariableFeatures(nb_malignant,
                                      selection.method = "vst",
                                      nfeatures = 3000,
                                      verbose = FALSE)
message("    Found ", length(VariableFeatures(nb_malignant)), " variable features")

message("  Scaling data...")
nb_malignant <- ScaleData(nb_malignant, verbose = FALSE)

message("  Running PCA on malignant-specific variable features...")
nb_malignant <- RunPCA(nb_malignant, npcs = 50, verbose = FALSE)

message("  Running Harmony batch correction...")

if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("ERROR: harmony package is required for batch correction. Please install it.")
}

library(harmony)

batch_col <- NULL
if ("Patient" %in% colnames(nb_malignant@meta.data)) {
  batch_col <- "Patient"
} else if ("patient" %in% colnames(nb_malignant@meta.data)) {
  batch_col <- "patient"
} else if ("Batch" %in% colnames(nb_malignant@meta.data)) {
  batch_col <- "Batch"
}

if (!is.null(batch_col)) {
  message("    Using '", batch_col, "' for batch correction")
  message("    Patients: ", paste(unique(nb_malignant@meta.data[[batch_col]]), collapse = ", "))

  nb_malignant <- RunHarmony(nb_malignant,
                              group.by.vars = batch_col,
                              reduction = "pca",
                              reduction.save = "harmony",
                              theta = 2,
                              verbose = FALSE)
  reduction_to_use <- "harmony"
  message("    ✓ Harmony batch correction completed")
} else {
  warning("No batch column (Patient/patient/Batch) found! Using PCA without batch correction.")
  reduction_to_use <- "pca"
}

message("  Running UMAP...")
nb_malignant <- RunUMAP(nb_malignant,
                         reduction = reduction_to_use,
                         dims = 1:N_PCS,
                         verbose = FALSE)

message("  ✓ Dimensionality reduction complete")

## =============================================================================
## STEP 4: RE-CLUSTER AT MULTIPLE RESOLUTIONS
## =============================================================================

message("\n[4/6] Re-clustering at multiple resolutions...")

message("  Building neighbor graph...")
nb_malignant <- FindNeighbors(nb_malignant,
                               reduction = reduction_to_use,
                               dims = 1:N_PCS,
                               verbose = FALSE)

for (res in CLUSTER_RESOLUTIONS) {
  message("  Clustering at resolution ", res, "...")
  nb_malignant <- FindClusters(nb_malignant,
                                resolution = res,
                                verbose = FALSE)

  col_name <- paste0("malignant_clusters_res", res)
  nb_malignant@meta.data[[col_name]] <- nb_malignant$seurat_clusters

  n_clusters <- length(unique(nb_malignant@meta.data[[col_name]]))
  message("    Found ", n_clusters, " clusters at resolution ", res)
}

default_col <- paste0("malignant_clusters_res", DEFAULT_RESOLUTION)
nb_malignant$malignant_clusters <- nb_malignant@meta.data[[default_col]]
nb_malignant$seurat_clusters <- nb_malignant$malignant_clusters

message("\n  ✓ Clustering complete")
message("  Default resolution: ", DEFAULT_RESOLUTION)
message("  Default clusters: ", length(unique(nb_malignant$malignant_clusters)))

## =============================================================================
## STEP 5: GENERATE VISUALIZATIONS
## =============================================================================

message("\n[5/6] Generating visualizations...")

plot_dir <- file.path(OUTPUT_DIR, "plots")

message("  Creating original cluster UMAP...")
p1 <- DimPlot(nb_malignant,
              reduction = "umap",
              group.by = "original_cluster",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.5) +
  ggtitle("Original Clusters on Malignant UMAP") +
  theme(legend.position = "right")

ggsave(file.path(plot_dir, "01_original_clusters_on_malignant_UMAP.pdf"),
       p1, width = 10, height = 8)

message("  Creating new cluster UMAP (resolution ", DEFAULT_RESOLUTION, ")...")
p2 <- DimPlot(nb_malignant,
              reduction = "umap",
              group.by = "malignant_clusters",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.5) +
  ggtitle(paste0("Re-clustered Malignant Cells (res=", DEFAULT_RESOLUTION, ")")) +
  theme(legend.position = "right")

ggsave(file.path(plot_dir, "02_reclustered_malignant_UMAP.pdf"),
       p2, width = 10, height = 8)

message("  Creating comparison plot...")
p_compare <- p1 | p2
ggsave(file.path(plot_dir, "03_original_vs_reclustered_comparison.pdf"),
       p_compare, width = 18, height = 8)

message("  Creating multi-resolution comparison...")
resolution_plots <- list()
for (res in CLUSTER_RESOLUTIONS) {
  col_name <- paste0("malignant_clusters_res", res)
  p <- DimPlot(nb_malignant,
               reduction = "umap",
               group.by = col_name,
               label = TRUE,
               pt.size = 0.3) +
    ggtitle(paste0("Resolution ", res)) +
    theme(legend.position = "none")
  resolution_plots[[as.character(res)]] <- p
}

p_multi_res <- wrap_plots(resolution_plots, ncol = 2)
ggsave(file.path(plot_dir, "04_multi_resolution_comparison.pdf"),
       p_multi_res, width = 14, height = 14)

if ("Patient" %in% colnames(nb_malignant@meta.data)) {
  message("  Creating patient-colored UMAP...")
  p_patient <- DimPlot(nb_malignant,
                       reduction = "umap",
                       group.by = "Patient",
                       pt.size = 0.5) +
    ggtitle("Malignant Cells by Patient") +
    theme(legend.position = "right")

  ggsave(file.path(plot_dir, "05_malignant_by_patient.pdf"),
         p_patient, width = 10, height = 8)
}

if ("Fraction" %in% colnames(nb_malignant@meta.data)) {
  message("  Creating fraction-colored UMAP...")
  p_fraction <- DimPlot(nb_malignant,
                        reduction = "umap",
                        group.by = "Fraction",
                        cols = c("F1" = "blue", "F2" = "red"),
                        pt.size = 0.5) +
    ggtitle("Malignant Cells by Fraction (F1 vs F2)") +
    theme(legend.position = "right")

  ggsave(file.path(plot_dir, "06_malignant_by_fraction.pdf"),
         p_fraction, width = 10, height = 8)
}

jansky_cols <- c("jansky_transfer_label", "predicted.id", "jansky_singleR_label")
available_jansky <- jansky_cols[jansky_cols %in% colnames(nb_malignant@meta.data)]

if (length(available_jansky) > 0) {
  jansky_col <- available_jansky[1]
  message("  Creating Jansky annotation UMAP...")
  p_jansky <- DimPlot(nb_malignant,
                      reduction = "umap",
                      group.by = jansky_col,
                      label = TRUE,
                      repel = TRUE,
                      pt.size = 0.5) +
    ggtitle(paste0("Malignant Cells - Jansky Annotations (", jansky_col, ")")) +
    theme(legend.position = "right",
          legend.text = element_text(size = 8))

  ggsave(file.path(plot_dir, "07_malignant_jansky_annotations.pdf"),
         p_jansky, width = 12, height = 8)
}

message("  ✓ Visualizations complete")

## =============================================================================
## STEP 6: SCORE GENE SIGNATURES
## =============================================================================

message("\n[5b/6] Scoring gene signatures on malignant cells...")

vanGroningen_sigs <- list()
if (file.exists(VANGRONINGEN_FILE)) {
  vg_df <- read_xlsx(VANGRONINGEN_FILE, sheet = 1)
  if ("Gene" %in% names(vg_df) && "Cell_State" %in% names(vg_df)) {
    adrn_genes <- vg_df %>% filter(Cell_State == "ADRN") %>% pull(Gene) %>% na.omit() %>% as.character()
    mes_genes <- vg_df %>% filter(Cell_State == "MES") %>% pull(Gene) %>% na.omit() %>% as.character()
    if (length(adrn_genes) > 0) vanGroningen_sigs[["ADRN"]] <- adrn_genes
    if (length(mes_genes) > 0) vanGroningen_sigs[["MES"]] <- mes_genes
  }
  message("  Loaded van Groningen: ADRN (", length(vanGroningen_sigs$ADRN), "), MES (", length(vanGroningen_sigs$MES), ")")
} else {
  message("  Warning: van Groningen file not found: ", VANGRONINGEN_FILE)
}

Dyer_sigs <- list()
if (file.exists(DYER_FILE)) {
  dyer_df <- read_xlsx(DYER_FILE, sheet = 1)
  if ("Gene" %in% names(dyer_df) && "Cell_State" %in% names(dyer_df)) {
    adrn_genes <- dyer_df %>% filter(Cell_State == "ADRN") %>% pull(Gene) %>% na.omit() %>% as.character()
    mes_genes <- dyer_df %>% filter(Cell_State == "MES") %>% pull(Gene) %>% na.omit() %>% as.character()
    if (length(adrn_genes) > 0) Dyer_sigs[["ADRN"]] <- adrn_genes
    if (length(mes_genes) > 0) Dyer_sigs[["MES"]] <- mes_genes
  }
  message("  Loaded Dyer signatures")
} else {
  message("  Warning: Dyer file not found: ", DYER_FILE)
}

Int_sigs <- list()
if (file.exists(INT_SIGNATURE_CSV)) {
  int_df <- read_csv(INT_SIGNATURE_CSV, show_col_types = FALSE)
  gene_col <- if ("Gene" %in% colnames(int_df)) "Gene" else colnames(int_df)[1]
  int_genes <- unique(na.omit(as.character(int_df[[gene_col]])))
  int_genes <- intersect(int_genes, rownames(nb_malignant))
  if (length(int_genes) > 0) Int_sigs[["Int"]] <- int_genes
  message("  Loaded Int signature: ", length(int_genes), " genes")
} else {
  message("  Warning: Int signature file not found: ", INT_SIGNATURE_CSV)
}

all_sigs <- list()
for (sig_name in names(vanGroningen_sigs)) {
  all_sigs[[paste0("vanGroningen_", sig_name)]] <- vanGroningen_sigs[[sig_name]]
}
for (sig_name in names(Dyer_sigs)) {
  all_sigs[[paste0("Dyer_", sig_name)]] <- Dyer_sigs[[sig_name]]
}
for (sig_name in names(Int_sigs)) {
  all_sigs[[sig_name]] <- Int_sigs[[sig_name]]
}

if (length(all_sigs) > 0) {
  for (sig_name in names(all_sigs)) {
    message("  Scoring: ", sig_name, " (", length(all_sigs[[sig_name]]), " genes)")
    nb_malignant <- AddModuleScore(
      object = nb_malignant,
      features = list(all_sigs[[sig_name]]),
      name = paste0(sig_name, "_Score"),
      assay = "SCT",
      ctrl = 100,
      seed = 12345
    )
    old_name <- paste0(sig_name, "_Score1")
    if (old_name %in% colnames(nb_malignant@meta.data)) {
      colnames(nb_malignant@meta.data)[colnames(nb_malignant@meta.data) == old_name] <- sig_name
    }
  }

  sig_cols <- names(all_sigs)
  message("  ✓ Scored ", length(sig_cols), " signatures (using SCT assay)")
}

## =============================================================================
## STEP 7: SAVE FINAL RESULTS
## =============================================================================

message("\n[6/6] Saving final results...")

output_file <- file.path(DATA_DIR, "NB_patients_malignant_reclustered.rds")
saveRDS(nb_malignant, output_file)
message("  Saved: ", basename(output_file))

output_file_copy <- file.path(OUTPUT_DIR, "NB_patients_malignant_complete.rds")
saveRDS(nb_malignant, output_file_copy)
message("  Saved: ", basename(output_file_copy), " (complete object with all metadata)")

cluster_mapping <- nb_malignant@meta.data %>%
  select(starts_with("original_cluster"), starts_with("malignant_clusters")) %>%
  distinct() %>%
  arrange(original_cluster)

write.csv(cluster_mapping,
          file.path(OUTPUT_DIR, "tables", "cluster_mapping.csv"),
          row.names = TRUE)
message("  Saved: cluster_mapping.csv")

write.csv(nb_malignant@meta.data,
          file.path(OUTPUT_DIR, "tables", "malignant_cell_metadata_complete.csv"),
          row.names = TRUE)
message("  Saved: malignant_cell_metadata_complete.csv")

## =============================================================================
## SUMMARY
## =============================================================================

message("\n============================================================")
message("RE-CLUSTERING & ANALYSIS COMPLETE!")
message("============================================================")
message("\nSummary:")
message("  Original cells: ", n_original)
message("  Malignant cells: ", n_malignant, " (", pct_malignant, "%)")
message("  Original malignant clusters: ", paste(MALIGNANT_CLUSTERS, collapse = ", "))
message("  New clusters (res=", DEFAULT_RESOLUTION, "): ",
        length(unique(nb_malignant$malignant_clusters)))

message("\nCluster counts at each resolution:")
for (res in CLUSTER_RESOLUTIONS) {
  col_name <- paste0("malignant_clusters_res", res)
  n_clust <- length(unique(nb_malignant@meta.data[[col_name]]))
  message("  Resolution ", res, ": ", n_clust, " clusters")
}

message("\nOutput files:")
message("  Main Seurat object: ", output_file)
message("  Complete object (for sharing): ", output_file_copy)
message("  Plots: ", file.path(OUTPUT_DIR, "plots"))
message("  Tables: ", file.path(OUTPUT_DIR, "tables"))

message("\n✓ Malignant cell analysis complete!")
message("============================================================\n")
