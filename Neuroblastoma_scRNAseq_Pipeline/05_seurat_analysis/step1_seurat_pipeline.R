#!/usr/bin/env Rscript
## =============================================================================
## Neuroblastoma Patient Pipeline - Step 1: Seurat Analysis Pipeline
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Complete pipeline from normalized SCE to publication-ready objects
## Input:   Normalized SCE object from normalization step
## Output:  Seurat objects, QC visualizations, cluster assignments
## =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(scuttle)
  library(Seurat)
  library(Matrix)
  library(irlba)
  library(BiocParallel)
  library(biomaRt)
  library(EnsDb.Hsapiens.v86)
  library(patchwork)
  library(harmony)
  library(ggplot2)
})

set.seed(12345)

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

# Project and data directories
PROJECT_DIR    <- "/path/to/NB_patient_analysis"
RESULTS_DIR    <- file.path(PROJECT_DIR, "results")
DATA_DIR       <- file.path(RESULTS_DIR, "data")
PLOT_DIR       <- file.path(RESULTS_DIR, "plots")

# Input files
SCE_PATH       <- file.path(RESULTS_DIR, "NB_patients_SCE-norm_relaxed_cleaned.RDS")
BIOMART_PATH   <- file.path(PROJECT_DIR, "src", "ensembl_biomart.csv")

# Output prefix
OUTPUT_PREFIX  <- "NB_patients"

# Analysis parameters
RUN_HARMONY    <- TRUE    # Set FALSE to skip Harmony batch correction
CLUSTER_RES    <- 0.2     # Clustering resolution
N_PCS          <- 30      # Number of PCs for UMAP/clustering
N_VARIABLE     <- 1000    # Number of variable features for SCTransform

# QC thresholds
MIN_GENES      <- 200     # Minimum genes per cell
MAX_MT_PERCENT <- 10      # Maximum mitochondrial percentage
MIN_RIBO       <- 0.1     # Minimum ribosomal percentage

# Patient ID mapping (Batch number -> Patient ID)
BATCH_TO_PATIENT <- c(
  "1" = "Patient_770",
  "2" = "Patient_637",
  "3" = "Patient_644",
  "4" = "Patient_649",
  "5" = "Patient_1087"
)

## =============================================================================
## SETUP
## =============================================================================

# Create output directories
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

message("=============================================================================")
message("Neuroblastoma Patient Analysis Pipeline")
message("=============================================================================")
message("Configuration:")
message("  Input SCE: ", SCE_PATH)
message("  Output directory: ", RESULTS_DIR)
message("  Plots directory: ", PLOT_DIR)
message("  Data directory: ", DATA_DIR)
message("  Run Harmony: ", RUN_HARMONY)
message("  Clustering resolution: ", CLUSTER_RES)
message("=============================================================================\n")

# Check input file exists
if (!file.exists(SCE_PATH)) {
  stop("ERROR: SCE file not found at: ", SCE_PATH)
}

## =============================================================================
## STEP 1: LOAD DATA
## =============================================================================

message("[1/12] Loading normalized SCE from: ", SCE_PATH)
nb.sce <- readRDS(SCE_PATH)
message("  Dimensions: ", nrow(nb.sce), " genes x ", ncol(nb.sce), " cells")
message("  Metadata columns: ", paste(colnames(colData(nb.sce)), collapse = ", "))

## =============================================================================
## STEP 2: GENE ANNOTATION
## =============================================================================

message("\n[2/12] Loading gene annotations...")

# Load or download biomart annotation
if (file.exists(BIOMART_PATH)) {
  message("  Using existing biomart annotation: ", BIOMART_PATH)
  ens.bm <- read.csv(BIOMART_PATH, stringsAsFactors = FALSE)
} else {
  message("  Downloading from Ensembl (hsapiens_gene_ensembl)...")
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ens.bm <- getBM(
    attributes = c("ensembl_gene_id", "chromosome_name", "external_gene_name", "description"),
    mart = mart
  )
  write.csv(ens.bm, BIOMART_PATH, row.names = FALSE)
  message("  Saved annotation to: ", BIOMART_PATH)
}

# Get lincRNA genes from EnsDb
ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes   <- ens.86.genes$gene_id[ens.86.genes$gene_biotype == "lincRNA"]

# Filter biomart annotation
ens.filt <- ens.bm
ens.filt$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt$description)
ens.filt <- ens.filt[ens.filt$chromosome_name %in% c(1:22, "X", "Y"), ]
ens.filt <- ens.filt[ens.filt$external_gene_name != "", ]
ens.filt <- ens.filt[!duplicated(ens.filt$external_gene_name), ]

# Convert gene IDs to symbols if needed
gene_ids <- rownames(nb.sce)
is_ensembl <- grepl("^ENSG", gene_ids[1])

if (is_ensembl) {
  message("  Converting Ensembl IDs to gene symbols...")
  match_idx <- match(gene_ids, ens.filt$ensembl_gene_id)
  keep_gene <- !is.na(match_idx)
  if (!any(keep_gene)) {
    stop("ERROR: No Ensembl gene IDs matched the biomart reference")
  }
  nb.sce <- nb.sce[keep_gene, ]
  match_idx <- match_idx[keep_gene]
  ens.filt <- ens.filt[match_idx, ]
  rownames(nb.sce) <- ens.filt$external_gene_name
  message("  Converted ", sum(keep_gene), " genes to symbols")
} else {
  message("  Gene symbols detected, matching to biomart...")
  match_idx <- match(gene_ids, ens.filt$external_gene_name)
  keep_gene <- !is.na(match_idx)
  if (!any(keep_gene)) {
    stop("ERROR: Gene symbols could not be matched to biomart reference")
  }
  nb.sce <- nb.sce[keep_gene, ]
  match_idx <- match_idx[keep_gene]
  ens.filt <- ens.filt[match_idx, ]
  message("  Matched ", sum(keep_gene), " genes")
}

rownames(ens.filt) <- rownames(nb.sce)

# Define QC gene sets
mt.genes   <- grep("^MT-",  rownames(nb.sce), value = TRUE)
ribo.genes <- grep("^RP[SL]", rownames(nb.sce), value = TRUE)
linc.genes <- intersect(linc.genes, ens.bm$ensembl_gene_id)
linc.genes <- ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% linc.genes]
linc.genes <- intersect(linc.genes, rownames(nb.sce))

message("  QC gene sets: ", length(mt.genes), " MT, ",
        length(ribo.genes), " Ribo, ", length(linc.genes), " lincRNA")

## =============================================================================
## STEP 3: CALCULATE QC METRICS
## =============================================================================

message("\n[3/12] Calculating per-cell QC metrics...")
nb.sce <- addPerCellQCMetrics(
  nb.sce,
  flatten = TRUE,
  subsets = list(
    mt = mt.genes,
    ribo = ribo.genes,
    linc = linc.genes
  )
)

# Add gene annotation to rowData
rowData(nb.sce) <- ens.filt

# Save pre-QC object
pre_qc_path <- file.path(DATA_DIR, paste0(OUTPUT_PREFIX, "_sce_preQC.rds"))
message("  Saving pre-QC SCE: ", basename(pre_qc_path))
saveRDS(nb.sce, pre_qc_path)

## =============================================================================
## STEP 4: VERIFY PATIENT AND SAMPLE METADATA
## =============================================================================

message("\n[4/12] Verifying patient and sample metadata...")

# Cleaned object already has Patient_ID - rename to Patient for consistency
if ("Patient_ID" %in% colnames(colData(nb.sce))) {
  nb.sce$Patient <- nb.sce$Patient_ID
  message("  Patient column verified (from Patient_ID)")
} else if (!"Patient" %in% colnames(colData(nb.sce))) {
  # Fallback: map from Batch if needed
  nb.sce$Patient <- BATCH_TO_PATIENT[as.character(nb.sce$Batch)]
  message("  Created Patient column from Batch mapping")
}
message("  Found ", length(unique(nb.sce$Patient)), " patients")

# Verify Sample column exists
if (!"Sample" %in% colnames(colData(nb.sce))) {
  stop("ERROR: Sample column missing from cleaned object!")
}
message("  Found ", length(unique(nb.sce$Sample)), " unique samples")

# Print sample distribution
sample_counts <- table(nb.sce$Sample)
message("  Sample distribution:")
for (s in names(sample_counts)) {
  message("    ", s, ": ", sample_counts[s], " cells")
}

## =============================================================================
## STEP 5: CREATE SEURAT OBJECT
## =============================================================================

message("\n[5/12] Creating Seurat object with full metadata...")

# Select metadata columns to keep
meta_keep <- intersect(
  colnames(colData(nb.sce)),
  c("Sample", "Barcode", "Batch", "Patient", "Patient_ID",
    "Fraction", "Replicate", "Tumor_Type", "Sequencing_Run",
    "scDblFinder.class", "scDblFinder.score", "scDblFinder.weighted", "scDblFinder.cxds_score",
    "detected", "total", "subsets_mt_percent", "subsets_ribo_percent", "subsets_linc_percent")
)
metadata.df <- as.data.frame(colData(nb.sce))[, meta_keep, drop = FALSE]

# Create Seurat object
nb.seurat <- CreateSeuratObject(
  counts = assay(nb.sce, "counts"),
  assay = "RNA",
  meta.data = metadata.df
)
message("  Created Seurat object: ", ncol(nb.seurat), " cells, ", nrow(nb.seurat), " genes")

## =============================================================================
## STEP 6: CELL CYCLE SCORING
## =============================================================================

message("\n[6/12] Scoring cell cycle phases...")

# Normalize for cell cycle scoring
nb.seurat <- NormalizeData(nb.seurat, verbose = FALSE)

# Score cell cycle
cc.genes <- Seurat::cc.genes.updated.2019
nb.seurat <- CellCycleScoring(
  nb.seurat,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

# Add cell cycle scores back to SCE
nb.sce$Seurat.Phase <- nb.seurat$Phase
nb.sce$Seurat.S     <- nb.seurat$S.Score
nb.sce$Seurat.G2M   <- nb.seurat$G2M.Score

# Print cell cycle distribution
phase_counts <- table(nb.seurat$Phase)
message("  Cell cycle distribution:")
for (p in names(phase_counts)) {
  message("    ", p, ": ", phase_counts[p], " cells (",
          round(100 * phase_counts[p] / sum(phase_counts), 1), "%)")
}

## =============================================================================
## STEP 7: QUALITY CONTROL FILTERING
## =============================================================================

message("\n[7/12] Applying QC filters...")
message("  Thresholds:")
message("    Min genes/cell: ", MIN_GENES)
message("    Max MT%: ", MAX_MT_PERCENT)
message("    Min Ribo%: ", MIN_RIBO)

# Apply filters
detected.filt  <- colnames(nb.sce)[nb.sce$detected >= MIN_GENES]
rowcounts.filt <- rownames(nb.sce)[Matrix::rowSums(counts(nb.sce)) > 5]
mt.genes.filt  <- grep("^MT-", rownames(nb.sce), invert = TRUE, value = TRUE)

genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

# Filter SCE
nb.filt.sce <- nb.sce[genes.filt, detected.filt]
mito.filt <- colnames(nb.filt.sce)[nb.filt.sce$subsets_mt_percent < MAX_MT_PERCENT]
ribo.filt <- colnames(nb.filt.sce)[nb.filt.sce$subsets_ribo_percent > MIN_RIBO]
qc.filt <- intersect(mito.filt, ribo.filt)
nb.filt.sce <- nb.filt.sce[, qc.filt]

# Filter Seurat object (metadata automatically preserved)
nb.filt.seurat <- subset(nb.seurat, cells = qc.filt, features = genes.filt)

# IMPORTANT: Verify cell cycle scores were transferred
if (!all(c("Phase", "S.Score", "G2M.Score") %in% colnames(nb.filt.seurat@meta.data))) {
  stop("ERROR: Cell cycle scores not found in filtered Seurat object!")
}

message("  Cell cycle scores verified in filtered object")

# Report filtering results
message("  Before filtering: ", ncol(nb.sce), " cells, ", nrow(nb.sce), " genes")
message("  After filtering:  ", ncol(nb.filt.sce), " cells, ", nrow(nb.filt.sce), " genes")
message("  Cells removed: ", ncol(nb.sce) - ncol(nb.filt.sce),
        " (", round(100 * (ncol(nb.sce) - ncol(nb.filt.sce)) / ncol(nb.sce), 1), "%)")
message("  Genes removed: ", nrow(nb.sce) - nrow(nb.filt.sce))

# Save filtered objects
sce_filtered_path    <- file.path(DATA_DIR, paste0(OUTPUT_PREFIX, "_sce_filtered.rds"))
seurat_filtered_path <- file.path(DATA_DIR, paste0(OUTPUT_PREFIX, "_seurat_filtered.rds"))
message("  Saving filtered objects:")
message("    ", basename(sce_filtered_path))
message("    ", basename(seurat_filtered_path))
saveRDS(nb.filt.sce,    sce_filtered_path)
saveRDS(nb.filt.seurat, seurat_filtered_path)

## =============================================================================
## STEP 8: SCTRANSFORM NORMALIZATION
## =============================================================================

message("\n[8/12] Running SCTransform normalization...")
message("  Method: glmGamPoi")
message("  Regressing out: Cell cycle score, MT%")
message("  Variable features: ", N_VARIABLE)

DefaultAssay(nb.filt.seurat) <- "RNA"
nb.filt.seurat$Seurat.Cycle.Score <- nb.filt.seurat$S.Score - nb.filt.seurat$G2M.Score

nb.filt.seurat <- SCTransform(
  nb.filt.seurat,
  method = "glmGamPoi",
  vst.flavor = "v2",
  vars.to.regress = c("Seurat.Cycle.Score", "subsets_mt_percent"),
  variable.features.n = N_VARIABLE,
  do.scale = TRUE,
  do.center = TRUE,
  verbose = TRUE
)

# Set SCT as default assay after SCTransform
DefaultAssay(nb.filt.seurat) <- "SCT"
message("  SCTransform complete (SCT set as default assay)")

## =============================================================================
## STEP 9: PCA + UMAP (NO HARMONY)
## =============================================================================

message("\n[9/12] Running PCA and UMAP (no batch correction)...")
message("  Computing 100 PCs...")

nb.noharm <- nb.filt.seurat
nb.noharm <- RunPCA(
  nb.noharm,
  npcs = 100,
  verbose = FALSE,
  features = nb.noharm@assays$SCT@var.features
)

message("  Running UMAP with ", N_PCS, " PCs...")
nb.noharm <- RunUMAP(nb.noharm, reduction = "pca", dims = 1:N_PCS)
nb.noharm <- FindNeighbors(nb.noharm, reduction = "pca", dims = 1:N_PCS)
nb.noharm <- FindClusters(nb.noharm, resolution = CLUSTER_RES)

# Ensure SCT is default assay
DefaultAssay(nb.noharm) <- "SCT"

n_clusters <- length(unique(nb.noharm$seurat_clusters))
message("  Found ", n_clusters, " clusters at resolution ", CLUSTER_RES)

# Save no-Harmony object
noharmony_path <- file.path(DATA_DIR, paste0(OUTPUT_PREFIX, "_seurat_noHarmony.rds"))
message("  Saving: ", basename(noharmony_path))
saveRDS(nb.noharm, noharmony_path)

# Save individual no-Harmony plots
message("  Creating individual no-Harmony plots...")

p1 <- DimPlot(nb.noharm, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("No Harmony - Clusters")
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_clusters.pdf")),
       plot = p1, width = 9, height = 7)

p2 <- DimPlot(nb.noharm, group.by = "Patient") +
  ggtitle("No Harmony - By Patient")
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_patient.pdf")),
       plot = p2, width = 9, height = 7)

p3 <- DimPlot(nb.noharm, group.by = "Sample", label = FALSE) +
  ggtitle(paste0("No Harmony - By Sample (n=", length(unique(nb.noharm$Sample)), ")")) +
  theme(legend.text = element_text(size = 8))
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_samples.pdf")),
       plot = p3, width = 12, height = 7)

p4 <- FeaturePlot(nb.noharm, features = "detected", order = TRUE) +
  ggtitle("Detected features (no Harmony)")
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_detected.pdf")),
       plot = p4, width = 9, height = 7)

p5 <- FeaturePlot(nb.noharm, features = "subsets_mt_percent", order = TRUE) +
  ggtitle("% MT genes (no Harmony)")
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_mt.pdf")),
       plot = p5, width = 9, height = 7)

p6 <- FeaturePlot(nb.noharm, features = "S.Score", order = TRUE) +
  ggtitle("S-phase score (no Harmony)")
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_S_score.pdf")),
       plot = p6, width = 9, height = 7)

p7 <- FeaturePlot(nb.noharm, features = "G2M.Score", order = TRUE) +
  ggtitle("G2M-phase score (no Harmony)")
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_G2M_score.pdf")),
       plot = p7, width = 9, height = 7)

p8 <- DimPlot(nb.noharm, group.by = "Phase") +
  ggtitle("Cell Cycle Phases (no Harmony)")
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_noHarmony_phases.pdf")),
       plot = p8, width = 9, height = 7)

## =============================================================================
## STEP 10: HARMONY BATCH CORRECTION (OPTIONAL)
## =============================================================================

if (RUN_HARMONY) {
  message("\n[10/12] Running Harmony batch correction...")
  message("  Batch variable: Batch (Patient ID)")
  message("  Theta: 2")
  message("  Max iterations: 50")

  # First run PCA on the filtered object
  message("  Running PCA (100 PCs) before Harmony...")
  nb.harm <- RunPCA(
    nb.filt.seurat,
    npcs = 100,
    verbose = FALSE,
    features = nb.filt.seurat@assays$SCT@var.features
  )

  # Now run Harmony on the PCA embeddings
  message("  Running Harmony integration...")
  nb.harm <- RunHarmony(
    nb.harm,
    group.by.vars = "Batch",
    assay.use = "SCT",
    reduction.save = "Harmony",
    theta = 2,
    max.iter.harmony = 50,
    plot_convergence = FALSE
  )

  message("  Running UMAP on Harmony embeddings...")
  nb.harm <- RunUMAP(nb.harm, reduction = "Harmony", dims = 1:N_PCS)
  nb.harm <- FindNeighbors(nb.harm, reduction = "Harmony", dims = 1:N_PCS)
  nb.harm <- FindClusters(nb.harm, resolution = CLUSTER_RES)

  # Ensure SCT is default assay
  DefaultAssay(nb.harm) <- "SCT"

  n_clusters_harm <- length(unique(nb.harm$seurat_clusters))
  message("  Found ", n_clusters_harm, " clusters at resolution ", CLUSTER_RES)

  # Save Harmony object
  harmony_path <- file.path(DATA_DIR, paste0(OUTPUT_PREFIX, "_seurat_harmony.rds"))
  message("  Saving: ", basename(harmony_path))
  saveRDS(nb.harm, harmony_path)

  # Create multi-page PDF with Harmony plots
  message("  Creating Harmony plots...")
  pdf(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony.pdf")),
      width = 9, height = 7)
  print(DimPlot(nb.harm, group.by = "seurat_clusters", label = TRUE) +
          ggtitle("Clusters (Harmony)"))
  print(DimPlot(nb.harm, group.by = "Patient") +
          ggtitle("By Patient (Harmony)"))
  print(DimPlot(nb.harm, group.by = "Sample") +
          ggtitle("By Sample (Harmony)") + NoLegend())
  print(FeaturePlot(nb.harm, features = "detected", order = TRUE) +
          ggtitle("Detected genes"))
  print(FeaturePlot(nb.harm, features = "subsets_mt_percent", order = TRUE) +
          ggtitle("% MT genes"))
  dev.off()

  # Save individual Harmony UMAP plots
  message("  Saving individual Harmony plots...")

  p1 <- DimPlot(nb.harm, group.by = "seurat_clusters", label = TRUE) +
    ggtitle("Harmony - Clusters")
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_clusters.pdf")),
         plot = p1, width = 9, height = 7)

  p2 <- DimPlot(nb.harm, group.by = "Patient") +
    ggtitle("Harmony - By Patient")
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_patient.pdf")),
         plot = p2, width = 9, height = 7)

  p3 <- DimPlot(nb.harm, group.by = "Sample", label = FALSE) +
    ggtitle(paste0("Harmony - By Sample (n=", length(unique(nb.harm$Sample)), ")")) +
    theme(legend.text = element_text(size = 8))
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_samples.pdf")),
         plot = p3, width = 12, height = 7)

  p4 <- FeaturePlot(nb.harm, features = "detected", order = TRUE) +
    ggtitle("Detected features")
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_detected.pdf")),
         plot = p4, width = 9, height = 7)

  p5 <- FeaturePlot(nb.harm, features = "subsets_mt_percent", order = TRUE) +
    ggtitle("% MT genes")
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_mt.pdf")),
         plot = p5, width = 9, height = 7)

  p6 <- FeaturePlot(nb.harm, features = "S.Score", order = TRUE) +
    ggtitle("S-phase score")
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_S_score.pdf")),
         plot = p6, width = 9, height = 7)

  p7 <- FeaturePlot(nb.harm, features = "G2M.Score", order = TRUE) +
    ggtitle("G2M-phase score")
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_G2M_score.pdf")),
         plot = p7, width = 9, height = 7)

  p8 <- DimPlot(nb.harm, group.by = "Phase") +
    ggtitle("Cell Cycle Phases")
  ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_UMAP_Harmony_phases.pdf")),
         plot = p8, width = 9, height = 7)

} else {
  message("\n[10/12] Skipping Harmony integration (RUN_HARMONY = FALSE)")
}

## =============================================================================
## STEP 11: QC VISUALIZATION
## =============================================================================

message("\n[11/12] Creating QC visualizations...")

# Violin plots by sample
message("  Creating violin plots...")
qc_theme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

violin_pdf <- file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_QC_violin.pdf"))
pdf(violin_pdf, width = 10, height = 8)

# Detected genes
print(ggplot(as.data.frame(colData(nb.filt.sce)),
             aes(x = Sample, y = detected, fill = Sample)) +
        geom_violin(trim = FALSE) +
        stat_summary(fun = median, geom = "crossbar", width = 0.7, linewidth = 0.6) +
        labs(y = "Genes / cell", title = "Detected genes per cell") +
        qc_theme)

# UMI counts
print(ggplot(as.data.frame(colData(nb.filt.sce)),
             aes(x = Sample, y = total, fill = Sample)) +
        geom_violin(trim = FALSE) +
        stat_summary(fun = median, geom = "crossbar", width = 0.7, linewidth = 0.6) +
        scale_y_log10(labels = scales::comma) +
        labs(y = "UMI / cell", title = "Detected UMI per cell") +
        qc_theme)

# Mitochondrial percentage
print(ggplot(as.data.frame(colData(nb.filt.sce)),
             aes(x = Sample, y = subsets_mt_percent, fill = Sample)) +
        geom_violin(trim = FALSE) +
        stat_summary(fun = median, geom = "crossbar", width = 0.7, linewidth = 0.6) +
        labs(y = "Mitochondrial %", title = "Proportion of mitochondrial genes") +
        qc_theme)

# Ribosomal percentage
print(ggplot(as.data.frame(colData(nb.filt.sce)),
             aes(x = Sample, y = subsets_ribo_percent, fill = Sample)) +
        geom_violin(trim = FALSE) +
        stat_summary(fun = median, geom = "crossbar", width = 0.7, linewidth = 0.6) +
        labs(y = "Ribosomal %", title = "Proportion of ribosomal genes") +
        qc_theme)

dev.off()

# Feature plots on UMAP
message("  Creating feature plots on UMAP...")

umap_obj <- if (RUN_HARMONY) nb.harm else nb.noharm
suffix <- if (RUN_HARMONY) "Harmony" else "noHarmony"

feature_pdf <- file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_QC_features_", suffix, ".pdf"))
pdf(feature_pdf, width = 12, height = 8)

print(FeaturePlot(umap_obj, features = "detected", order = TRUE) +
        ggtitle("Detected features"))
print(FeaturePlot(umap_obj, features = "subsets_mt_percent", order = TRUE) +
        ggtitle("% mitochondrial genes"))
print(FeaturePlot(umap_obj, features = "subsets_ribo_percent", order = TRUE) +
        ggtitle("% ribosomal genes"))
print(FeaturePlot(umap_obj, features = "S.Score", order = TRUE) +
        ggtitle("Seurat Cycle S-phase"))
print(FeaturePlot(umap_obj, features = "G2M.Score", order = TRUE) +
        ggtitle("Seurat Cycle G2M-phase"))
print(DimPlot(umap_obj, group.by = "Phase") +
        ggtitle("Seurat Cycle Phases"))

dev.off()

# Individual QC violin plots
message("  Creating individual QC violin plots by sample...")

qc_data <- as.data.frame(colData(nb.filt.sce))
sample_order <- names(sort(table(qc_data$Sample), decreasing = TRUE))
qc_data$Sample <- factor(qc_data$Sample, levels = sample_order)

library(RColorBrewer)
n_samples <- length(unique(qc_data$Sample))
sample_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_samples)

qc_violin_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

# Detected genes per cell
p_genes <- ggplot(qc_data, aes(x = Sample, y = detected, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_hline(yintercept = MIN_GENES, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.5) +
  scale_fill_manual(values = sample_colors) +
  labs(x = "", y = "Genes / cell", title = "Detected genes per cell") +
  qc_violin_theme +
  annotate("text", x = 1, y = MIN_GENES, label = paste0(MIN_GENES, " genes"),
           vjust = -0.5, hjust = 0, size = 3)

ggsave(
  file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_QC_violin_detected_genes.pdf")),
  plot = p_genes,
  width = max(10, n_samples * 0.6),
  height = 6
)

# Mitochondrial percentage
p_mito <- ggplot(qc_data, aes(x = Sample, y = subsets_mt_percent, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_hline(yintercept = MAX_MT_PERCENT, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 15, linetype = "dotted", color = "red", linewidth = 0.5, alpha = 0.5) +
  scale_fill_manual(values = sample_colors) +
  labs(x = "", y = "Mitochondrial\nexpression % of total",
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  qc_violin_theme +
  coord_cartesian(ylim = c(0, max(40, max(qc_data$subsets_mt_percent, na.rm = TRUE)))) +
  annotate("text", x = 1, y = MAX_MT_PERCENT, label = paste0(MAX_MT_PERCENT, "%"),
           vjust = -0.5, hjust = 0, size = 3)

ggsave(
  file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_QC_violin_mitochondrial.pdf")),
  plot = p_mito,
  width = max(10, n_samples * 0.6),
  height = 6
)

# Ribosomal percentage
p_ribo <- ggplot(qc_data, aes(x = Sample, y = subsets_ribo_percent, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_hline(yintercept = MIN_RIBO, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.5) +
  scale_fill_manual(values = sample_colors) +
  labs(x = "", y = "Ribosomal\nexpression % of total",
       title = "Proportion of ribosomal genes\nin total cell expression") +
  qc_violin_theme +
  coord_cartesian(ylim = c(0, max(8, max(qc_data$subsets_ribo_percent, na.rm = TRUE)))) +
  annotate("text", x = 1, y = MIN_RIBO, label = paste0(MIN_RIBO, "%"),
           vjust = -0.5, hjust = 0, size = 3)

ggsave(
  file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_QC_violin_ribosomal.pdf")),
  plot = p_ribo,
  width = max(10, n_samples * 0.6),
  height = 6
)

message("  Saved 3 individual QC violin plots")

## =============================================================================
## STEP 12: SUMMARY
## =============================================================================

message("\n[12/12] Pipeline complete!")
message("\n=============================================================================")
message("SUMMARY")
message("=============================================================================")
message("Input:")
message("  ", ncol(nb.sce), " cells, ", nrow(nb.sce), " genes (pre-QC)")
message("\nOutput:")
message("  ", ncol(nb.filt.sce), " cells, ", nrow(nb.filt.sce), " genes (post-QC)")
message("  ", length(unique(nb.filt.seurat$Patient)), " patients")
message("  ", length(unique(nb.filt.seurat$Sample)), " samples")
if (RUN_HARMONY) {
  message("  ", length(unique(nb.harm$seurat_clusters)), " clusters (Harmony, res=", CLUSTER_RES, ")")
} else {
  message("  ", length(unique(nb.noharm$seurat_clusters)), " clusters (no Harmony, res=", CLUSTER_RES, ")")
}
message("\nSaved objects in: ", DATA_DIR)
message("  ", basename(pre_qc_path))
message("  ", basename(sce_filtered_path))
message("  ", basename(seurat_filtered_path))
message("  ", basename(noharmony_path))
if (RUN_HARMONY) {
  message("  ", basename(harmony_path))
}
message("\nSaved plots in: ", PLOT_DIR)
message("  (", length(list.files(PLOT_DIR, pattern = "*.pdf")), " PDF files)")
message("\n=============================================================================")
message("Next steps:")
message("  1. Load Harmony object: readRDS('", harmony_path, "')")
message("  2. Find marker genes (step 1b)")
message("  3. Annotate clusters (step 1c)")
message("  4. Score gene signatures (step 1d)")
message("  5. Jansky annotation (step 1e)")
message("=============================================================================\n")

# Save session info
sessionInfo_path <- file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_sessionInfo.txt"))
writeLines(capture.output(sessionInfo()), sessionInfo_path)
message("Session info saved to: ", basename(sessionInfo_path))

message("\nDone!")
