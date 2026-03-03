#!/usr/bin/env Rscript
## =============================================================================
## Visualize InferCNV Results
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Create publication-quality UMAP visualizations of inferCNV results
## Input:   Seurat object with inferCNV scores from main analysis
## Output:  High-resolution UMAP plots with multiple color schemes
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(dplyr)
})

set.seed(123)

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")

# Input: Updated Seurat object with CNV scores
SEURAT_FILE <- file.path(DATA_DIR, "NB_patients_seurat_with_CNV_scores.rds")
INFERCNV_RESULTS_DIR <- file.path(RESULTS_DIR, "inferCNV_cluster10_NBAtlas_method")

OUTPUT_DIR <- file.path(RESULTS_DIR, "inferCNV_visualizations")

# Features to visualize (must exist in Seurat object)
CNV_SCORE_FEATURE <- "infercnv_score"

# UMAP reduction
UMAP_REDUCTION <- "umap"

# Additional chromosome-specific features (optional, comment out if not available)
CHROMOSOME_FEATURES <- c(
  "chr7_score", "chr17q_score", "chr1p_score", "chr2p_score",
  "chr3p_score", "chr4p_score", "chr11q_score", "chr14q_score"
)

# Output settings
DPI_STANDARD <- 300
DPI_HD <- 300

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(OUTPUT_DIR, "figures")
dir.create(fig_dir, showWarnings = FALSE)

message("============================================================")
message("INFERCNV RESULT VISUALIZATION")
message("============================================================")
message("Input Seurat:         ", SEURAT_FILE)
message("InferCNV results:     ", INFERCNV_RESULTS_DIR)
message("Output directory:     ", OUTPUT_DIR)
message("UMAP reduction:       ", UMAP_REDUCTION)
message("============================================================\n")

## =============================================================================
## STEP 1: LOAD SEURAT OBJECT WITH CNV SCORES
## =============================================================================

message("[1/5] Loading Seurat object with CNV scores...")

if (!file.exists(SEURAT_FILE)) {
  stop("ERROR: Seurat file not found: ", SEURAT_FILE)
}

seurat_obj <- readRDS(SEURAT_FILE)
message("  Loaded: ", ncol(seurat_obj), " cells")

if (!(UMAP_REDUCTION %in% names(seurat_obj@reductions))) {
  available_reductions <- paste(names(seurat_obj@reductions), collapse = ", ")
  stop("ERROR: ", UMAP_REDUCTION, " not found. Available: ", available_reductions)
}

if (!(CNV_SCORE_FEATURE %in% colnames(seurat_obj@meta.data))) {
  available_features <- paste(colnames(seurat_obj@meta.data), collapse = ", ")
  stop("ERROR: ", CNV_SCORE_FEATURE, " not found in metadata. Available: ", available_features)
}

message("  ✓ Seurat object validated")
message("  Metadata columns: ", paste(colnames(seurat_obj@meta.data)[1:5], collapse = ", "), ", ...")

## =============================================================================
## STEP 2: CREATE HELPER FUNCTIONS
## =============================================================================

message("\n[2/5] Defining visualization functions...")

create_cnv_umap <- function(seurat_obj, feature, reduction, title, output_path,
                            color_scheme = "PuOr", dpi = 300, hd = FALSE) {

  if (!(feature %in% colnames(seurat_obj@meta.data))) {
    warning("Feature ", feature, " not found in metadata. Skipping.")
    return(NULL)
  }

  feature_data <- seurat_obj@meta.data[[feature]]

  if (all(is.na(feature_data))) {
    warning("Feature ", feature, " has all NA values. Skipping.")
    return(NULL)
  }

  valid_cutoff_q2 <- quantile(feature_data, 0.02, na.rm = TRUE)
  valid_cutoff_q98 <- quantile(feature_data, 0.98, na.rm = TRUE)

  p <- FeaturePlot(
    seurat_obj,
    reduction = reduction,
    features = feature,
    order = TRUE,
    min.cutoff = valid_cutoff_q2,
    max.cutoff = valid_cutoff_q98,
    raster = !hd
  ) +
    scale_colour_gradientn(
      colours = rev(brewer.pal(n = 11, name = color_scheme)),
      name = feature,
      labels = function(x) round(x, 3)
    ) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      legend.text = element_text(size = 9)
    ) +
    guides(colour = guide_colorbar(title.position = "top", barwidth = 8, barheight = 0.5))

  ggsave(output_path, p, width = if(hd) 15 else 10, height = if(hd) 12 else 8, dpi = dpi)

  return(p)
}

message("  ✓ Functions defined")

## =============================================================================
## STEP 3: OVERALL INFERCNV SCORE - MULTIPLE COLOR SCHEMES
## =============================================================================

message("\n[3/5] Creating overall inferCNV score visualizations...")

message("  Creating RdBu color scheme (standard)...")
create_cnv_umap(
  seurat_obj,
  feature = CNV_SCORE_FEATURE,
  reduction = UMAP_REDUCTION,
  title = "InferCNV Score (Overall)",
  output_path = file.path(fig_dir, "InferCNV_score_UMAP_RdBu.png"),
  color_scheme = "RdBu",
  dpi = DPI_STANDARD
)

message("  Creating RdBu color scheme (HD)...")
create_cnv_umap(
  seurat_obj,
  feature = CNV_SCORE_FEATURE,
  reduction = UMAP_REDUCTION,
  title = "InferCNV Score (Overall)",
  output_path = file.path(fig_dir, "InferCNV_score_UMAP_RdBu_HD.png"),
  color_scheme = "RdBu",
  dpi = DPI_HD,
  hd = TRUE
)

message("  Creating PuOr color scheme (standard)...")
create_cnv_umap(
  seurat_obj,
  feature = CNV_SCORE_FEATURE,
  reduction = UMAP_REDUCTION,
  title = "InferCNV Score (Overall)",
  output_path = file.path(fig_dir, "InferCNV_score_UMAP_PuOr.png"),
  color_scheme = "PuOr",
  dpi = DPI_STANDARD
)

message("  Creating PuOr color scheme (HD)...")
create_cnv_umap(
  seurat_obj,
  feature = CNV_SCORE_FEATURE,
  reduction = UMAP_REDUCTION,
  title = "InferCNV Score (Overall)",
  output_path = file.path(fig_dir, "InferCNV_score_UMAP_PuOr_HD.png"),
  color_scheme = "PuOr",
  dpi = DPI_HD,
  hd = TRUE
)

message("  ✓ Overall score visualizations complete")

## =============================================================================
## STEP 4: CHROMOSOME-SPECIFIC SCORES
## =============================================================================

message("\n[4/5] Creating chromosome-specific CNV score visualizations...")

chr_features_present <- CHROMOSOME_FEATURES[CHROMOSOME_FEATURES %in% colnames(seurat_obj@meta.data)]

if (length(chr_features_present) == 0) {
  message("  WARNING: No chromosome-specific features found in metadata.")
  message("  Available features: ", paste(colnames(seurat_obj@meta.data)[1:10], collapse = ", "))
} else {
  message("  Found ", length(chr_features_present), " chromosome features")

  for (chr_feature in chr_features_present) {
    chr_name <- gsub("_score", "", chr_feature)
    message("    Processing: ", chr_feature)

    create_cnv_umap(
      seurat_obj,
      feature = chr_feature,
      reduction = UMAP_REDUCTION,
      title = paste0("CNV Score: ", chr_name),
      output_path = file.path(fig_dir, paste0("CNV_", chr_name, "_UMAP_PuOr.png")),
      color_scheme = "PuOr",
      dpi = DPI_STANDARD
    )

    create_cnv_umap(
      seurat_obj,
      feature = chr_feature,
      reduction = UMAP_REDUCTION,
      title = paste0("CNV Score: ", chr_name),
      output_path = file.path(fig_dir, paste0("CNV_", chr_name, "_UMAP_PuOr_HD.png")),
      color_scheme = "PuOr",
      dpi = DPI_HD,
      hd = TRUE
    )
  }

  message("  ✓ Chromosome-specific visualizations complete")
}

## =============================================================================
## STEP 5: CLUSTER IDENTITY REFERENCE PLOT
## =============================================================================

message("\n[5/5] Creating cluster identity reference plot...")

if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  p_clusters <- DimPlot(
    seurat_obj,
    reduction = UMAP_REDUCTION,
    group.by = "seurat_clusters",
    label = TRUE,
    label.size = 4,
    repel = TRUE,
    pt.size = 0.5
  ) +
    ggtitle("Cluster Identity (Reference for CNV patterns)") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  ggsave(
    file.path(fig_dir, "Clusters_identity_UMAP.png"),
    p_clusters,
    width = 12,
    height = 10,
    dpi = DPI_STANDARD
  )

  message("  Saved: Clusters_identity_UMAP.png")
}

if ("Patient" %in% colnames(seurat_obj@meta.data)) {
  p_patient <- DimPlot(
    seurat_obj,
    reduction = UMAP_REDUCTION,
    group.by = "Patient",
    pt.size = 0.5
  ) +
    ggtitle("Patient Identity (Reference)") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  ggsave(
    file.path(fig_dir, "Patient_identity_UMAP.png"),
    p_patient,
    width = 12,
    height = 10,
    dpi = DPI_STANDARD
  )

  message("  Saved: Patient_identity_UMAP.png")
}

if ("Fraction" %in% colnames(seurat_obj@meta.data)) {
  p_fraction <- DimPlot(
    seurat_obj,
    reduction = UMAP_REDUCTION,
    group.by = "Fraction",
    cols = c("F1" = "blue", "F2" = "red"),
    pt.size = 0.5
  ) +
    ggtitle("Fraction Identity (F1 pre-therapy vs F2 post-therapy)") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  ggsave(
    file.path(fig_dir, "Fraction_identity_UMAP.png"),
    p_fraction,
    width = 12,
    height = 10,
    dpi = DPI_STANDARD
  )

  message("  Saved: Fraction_identity_UMAP.png")
}

message("  ✓ Reference plots complete")

## =============================================================================
## CREATE README
## =============================================================================

readme_content <- "# InferCNV Visualization Results

## Overview
This directory contains publication-quality UMAP visualizations of inferCNV CNV detection results.

## File Descriptions

### Overall InferCNV Score
- **InferCNV_score_UMAP_RdBu.png** - Red-Blue color scheme (standard resolution)
- **InferCNV_score_UMAP_RdBu_HD.png** - Red-Blue color scheme (high resolution)
- **InferCNV_score_UMAP_PuOr.png** - Purple-Orange color scheme (standard resolution)
- **InferCNV_score_UMAP_PuOr_HD.png** - Purple-Orange color scheme (high resolution)

Higher CNV scores (red/purple) indicate more aneuploidy; lower scores (blue/orange) indicate diploid regions.

### Chromosome-Specific CNV Scores
- **CNV_chr*.png** - InferCNV scores for specific chromosomes/chromosome arms
  - chr7 (MYCN amplification)
  - chr17q (frequent deletion)
  - chr1p (frequent deletion)
  - chr2p (frequent amplification)
  - chr3p, chr4p, chr11q, chr14q

### Reference Plots
- **Clusters_identity_UMAP.png** - Seurat clustering identity (for reference)
- **Patient_identity_UMAP.png** - Patient identity (for comparing CNVs across patients)
- **Fraction_identity_UMAP.png** - Fraction identity (F1 pre-therapy vs F2 post-therapy)

## Technical Details

- Color schemes:
  - RdBu (Red-Blue): Traditional diverging scale, good for presentations
  - PuOr (Purple-Orange): Colorblind-friendly diverging scale

- Resolution:
  - Standard (300 dpi): For documents, presentations
  - HD (300 dpi, larger size): For detailed inspection, figure panels

## Interpretation

InferCNV uses expression patterns to infer copy number variation:
- Positive scores indicate likely aneuploidy (duplications/gains)
- Negative scores indicate likely loss (deletions)
- Scores near 0 indicate diploid status

## Citation

These visualizations are based on inferCNV (https://github.com/broadinstitute/inferCNV) analysis
results, using the NBAtlas-recommended parameters with Cluster 10 as the reference.

"

writeLines(readme_content, file.path(fig_dir, "README.md"))
message("  Created: README.md")

## =============================================================================
## SUMMARY REPORT
## =============================================================================

message("\n============================================================")
message("VISUALIZATION COMPLETE")
message("============================================================")

message("\nFiles generated:")
message("  Standard resolution (300 dpi):")
message("    - InferCNV_score_UMAP_RdBu.png")
message("    - InferCNV_score_UMAP_PuOr.png")
if (length(chr_features_present) > 0) {
  message("    - CNV_chr*_UMAP_PuOr.png (", length(chr_features_present), " files)")
}

message("\n  High resolution (300 dpi, larger):")
message("    - InferCNV_score_UMAP_RdBu_HD.png")
message("    - InferCNV_score_UMAP_PuOr_HD.png")
if (length(chr_features_present) > 0) {
  message("    - CNV_chr*_UMAP_PuOr_HD.png (", length(chr_features_present), " files)")
}

message("\n  Reference identity plots:")
message("    - Clusters_identity_UMAP.png")
message("    - Patient_identity_UMAP.png")
message("    - Fraction_identity_UMAP.png")
message("    - README.md")

message("\nOutput directory: ", fig_dir)

message("\n============================================================")
message("✓ InferCNV visualization complete!")
message("============================================================\n")
