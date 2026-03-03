#!/usr/bin/env Rscript
## =============================================================================
## Neuroblastoma Patient Pipeline - Step 1c: Manual Cluster Annotation
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Apply manual cell-type annotations based on marker analysis
## Input:   Harmony-integrated Seurat object, cluster_annotations.csv
## Output:  Annotated Seurat object with CellType_Manual metadata
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
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
TABLE_DIR      <- file.path(RESULTS_DIR, "tables_step1c")
PLOT_DIR       <- file.path(RESULTS_DIR, "plots_step1c")

# Input annotation file (user-created with cluster -> cell_type mapping)
ANNOTATION_CSV <- file.path(PROJECT_DIR, "cluster_annotations.csv")

# Auto-detect most recent Seurat object
# Priority: Jansky > Signatures > Harmony
SEURAT_CANDIDATES <- c(
  file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds"),
  file.path(DATA_DIR, "NB_patients_seurat_with_signatures.rds"),
  file.path(DATA_DIR, "NB_patients_seurat_harmony.rds")
)

SEURAT_PATH <- NULL
for (cand in SEURAT_CANDIDATES) {
  if (file.exists(cand)) {
    SEURAT_PATH <- cand
    break
  }
}

if (is.null(SEURAT_PATH)) {
  stop("No Seurat object found. Expected one of:\n",
       paste("  -", SEURAT_CANDIDATES, collapse = "\n"))
}

# Output to same location as input
ANNOTATED_OUT <- SEURAT_PATH
METADATA_OUT  <- file.path(TABLE_DIR, "manual_annotations_metadata.csv")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

ensure_dir(TABLE_DIR)
ensure_dir(PLOT_DIR)

## =============================================================================
## LOAD SEURAT OBJECT
## =============================================================================

message("Loading Seurat object from: ", SEURAT_PATH)
nb.seurat <- readRDS(SEURAT_PATH)
message("  Cells: ", ncol(nb.seurat), "  Genes: ", nrow(nb.seurat))

# Check what data is already present
message("\n  Checking existing metadata...")
if ("CellType_Manual" %in% colnames(nb.seurat@meta.data)) {
  message("    Manual annotations already present (will be updated)")
}
signature_cols <- c("vanGroningen_ADRN_score", "Dyer_ADRN_score", "Dyer_MES_score")
if (any(signature_cols %in% colnames(nb.seurat@meta.data))) {
  message("    Gene signatures detected (will be preserved)")
}
if ("jansky_singleR_label" %in% colnames(nb.seurat@meta.data)) {
  message("    Jansky annotations detected (will be preserved)")
}

clusters <- levels(nb.seurat$seurat_clusters)
message("\n  Detected clusters: ", paste(clusters, collapse = ", "))

## =============================================================================
## LOAD ANNOTATION TABLE
## =============================================================================

if (!file.exists(ANNOTATION_CSV)) {
  # Create template file for the user to fill in
  message("\nAnnotation file not found. Creating template at: ", ANNOTATION_CSV)
  template <- tibble(cluster = clusters, cell_type = NA_character_)
  write_csv(template, ANNOTATION_CSV)

  message("\n=============================================================================")
  message("TEMPLATE CREATED")
  message("=============================================================================")
  message("Please fill in the 'cell_type' column in:")
  message("  ", ANNOTATION_CSV)
  message("\nUse marker results from Step 1b to guide your annotations:")
  message("  ", file.path(RESULTS_DIR, "tables_step1b/FindAllMarkers_by_cluster.xlsx"))
  message("\nThen rerun this script.")
  message("=============================================================================\n")

  stop("Annotation template created. Please fill it in and rerun.")
}

message("\nLoading annotation table from: ", ANNOTATION_CSV)
anno_tbl <- read_csv(ANNOTATION_CSV, show_col_types = FALSE) %>%
  mutate(cluster = as.character(cluster))

if (!all(c("cluster", "cell_type") %in% names(anno_tbl))) {
  stop("Annotation file must contain columns 'cluster' and 'cell_type'.")
}

missing <- setdiff(clusters, anno_tbl$cluster)
if (length(missing) > 0) {
  stop("Annotation file is missing cluster IDs: ", paste(missing, collapse = ", "))
}

# Check for NA cell types
na_clusters <- anno_tbl$cluster[is.na(anno_tbl$cell_type)]
if (length(na_clusters) > 0) {
  stop("The following clusters have NA cell_type: ", paste(na_clusters, collapse = ", "),
       "\nPlease assign cell types to all clusters.")
}

anno_tbl <- anno_tbl %>% filter(cluster %in% clusters)

## =============================================================================
## ATTACH ANNOTATIONS TO SEURAT OBJECT
## =============================================================================

message("\nApplying manual annotations...")
cluster_to_type <- setNames(anno_tbl$cell_type, anno_tbl$cluster)
nb.seurat$CellType_Manual <- cluster_to_type[as.character(nb.seurat$seurat_clusters)]

if (any(is.na(nb.seurat$CellType_Manual))) {
  warn_clusters <- unique(nb.seurat$seurat_clusters[is.na(nb.seurat$CellType_Manual)])
  warning("Some cells remain unannotated: ", paste(warn_clusters, collapse = ", "))
}

# Print annotation summary
annot_summary_stats <- table(nb.seurat$CellType_Manual)
message("\n=============================================================================")
message("ANNOTATION SUMMARY")
message("=============================================================================")
for (ct in names(annot_summary_stats)) {
  message(sprintf("%-20s: %6d cells", ct, annot_summary_stats[ct]))
}
message("=============================================================================\n")

## =============================================================================
## VISUALIZATIONS
## =============================================================================

message("Creating annotation visualizations...")

# UMAP by cell type
p_umap <- DimPlot(nb.seurat, group.by = "CellType_Manual", label = TRUE,
                  repel = TRUE, label.size = 4) +
  ggtitle("Manual cell type annotations") +
  theme(legend.position = "right")

ggsave(file.path(PLOT_DIR, "UMAP_annotated_cell_types.pdf"),
       plot = p_umap, width = 12, height = 8)

# UMAP by cell type (no labels)
p_umap_nolab <- DimPlot(nb.seurat, group.by = "CellType_Manual") +
  ggtitle("Manual cell type annotations")

ggsave(file.path(PLOT_DIR, "UMAP_annotated_cell_types_no_labels.pdf"),
       plot = p_umap_nolab, width = 10, height = 8)

# Bar plot of cell counts per type
counts_df <- as.data.frame(table(nb.seurat$CellType_Manual))
names(counts_df) <- c("CellType", "Cells")

p_bar <- ggplot(counts_df, aes(x = reorder(CellType, -Cells), y = Cells, fill = CellType)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  guides(fill = "none") +
  labs(x = "Cell type", y = "Number of cells",
       title = "Cell counts by manual annotation") +
  theme(axis.text = element_text(size = 10))

ggsave(file.path(PLOT_DIR, "barplot_cell_type_counts.pdf"),
       plot = p_bar, width = 8, height = 6)

# Proportion by patient
prop_by_patient <- nb.seurat@meta.data %>%
  group_by(Patient, CellType_Manual) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Patient) %>%
  mutate(prop = n / sum(n))

p_prop <- ggplot(prop_by_patient, aes(x = Patient, y = prop, fill = CellType_Manual)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Proportion", fill = "Cell Type",
       title = "Cell type composition by patient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PLOT_DIR, "barplot_cell_type_by_patient.pdf"),
       plot = p_prop, width = 10, height = 6)

## =============================================================================
## SAVE RESULTS
## =============================================================================

message("\nSaving results...")

# Save annotated Seurat object
message("  Annotated Seurat object: ", ANNOTATED_OUT)
saveRDS(nb.seurat, ANNOTATED_OUT)

# Verify all data is preserved
message("\n  Verifying data preservation...")
all_metadata <- colnames(nb.seurat@meta.data)
signature_present <- sum(signature_cols %in% all_metadata)
jansky_present <- "jansky_singleR_label" %in% all_metadata
manual_present <- "CellType_Manual" %in% all_metadata

if (signature_present > 0) {
  message("    Gene signatures preserved (", signature_present, " columns)")
}
if (jansky_present) {
  message("    Jansky annotations preserved")
}
if (manual_present) {
  message("    Manual annotations present (CellType_Manual)")
}
message("    Total metadata columns: ", length(all_metadata))

# Save metadata table
message("  Metadata table: ", METADATA_OUT)
meta_out <- nb.seurat@meta.data %>%
  mutate(Cell = colnames(nb.seurat)) %>%
  relocate(Cell, Patient, Sample, Batch, seurat_clusters, CellType_Manual) %>%
  arrange(CellType_Manual, seurat_clusters)

write_csv(meta_out, METADATA_OUT)
message("    Exported ", ncol(meta_out), " metadata columns (preserving all previous data)")

# Save annotation mapping table
annot_map <- data.frame(
  cluster = names(cluster_to_type),
  cell_type = cluster_to_type
) %>% arrange(cluster)

write_csv(annot_map, file.path(TABLE_DIR, "cluster_to_celltype_mapping.csv"))

## =============================================================================
## COMPLETE
## =============================================================================

message("\n=============================================================================")
message("Manual cluster annotation complete!")
message("=============================================================================")
message("Updated Seurat object saved to:")
message("  ", ANNOTATED_OUT)
message("\n  All previous data preserved (signatures, Jansky annotations, etc.)")
message("  Manual annotations added/updated (CellType_Manual)")
message("\nPlots saved to:")
message("  ", PLOT_DIR)
message("\nMetadata exported to:")
message("  ", METADATA_OUT)
message("\nNext steps:")
message("  - Review annotation plots")
message("  - Continue with downstream analysis")
message("=============================================================================\n")

message("Done!")
