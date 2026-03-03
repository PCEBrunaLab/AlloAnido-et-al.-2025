#!/usr/bin/env Rscript
## =============================================================================
## Neuroblastoma Patient Pipeline - Step 1b: Find Cluster Markers
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Identify marker genes for each cluster using Seurat's FindAllMarkers
## Input:   Harmony-integrated Seurat object from Step 1
## Output:  Marker gene tables, visualizations, and dotplots
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(openxlsx)
  library(tidyr)
})

set.seed(12345)

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

# Project and data directories
PROJECT_DIR    <- "/path/to/NB_patient_analysis"
RESULTS_DIR    <- file.path(PROJECT_DIR, "results")
DATA_DIR       <- file.path(RESULTS_DIR, "data")
TABLE_DIR      <- file.path(RESULTS_DIR, "tables_step1b")
PLOT_DIR       <- file.path(RESULTS_DIR, "plots_step1b")

# Input Seurat object
SEURAT_PATH    <- file.path(DATA_DIR, "NB_patients_seurat_harmony.rds")

## Helper function
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

ensure_dir(TABLE_DIR)
ensure_dir(PLOT_DIR)

## =============================================================================
## LOAD SEURAT OBJECT
## =============================================================================

if (!file.exists(SEURAT_PATH)) {
  stop("Harmony Seurat object not found: ", SEURAT_PATH,
       "\nRun step1_seurat_pipeline.R before finding markers.")
}

message("Loading Seurat object from: ", SEURAT_PATH)
nb.seurat <- readRDS(SEURAT_PATH)
message("  Cells: ", ncol(nb.seurat), "  Genes: ", nrow(nb.seurat))

clusters <- levels(nb.seurat$seurat_clusters)
message("  Detected ", length(clusters), " clusters: ", paste(clusters, collapse = ", "))

## =============================================================================
## SET DEFAULT ASSAY
## =============================================================================

# Set default assay to SCT for all plotting and analysis
if ("SCT" %in% names(nb.seurat@assays)) {
  DefaultAssay(nb.seurat) <- "SCT"
  message("Using SCT assay for marker analysis")
} else {
  DefaultAssay(nb.seurat) <- "RNA"
  message("Warning: SCT assay not found, using RNA assay")
}

## =============================================================================
## CANONICAL MARKER FEATURE PLOTS
## =============================================================================

message("\nGenerating canonical marker feature plots...")

marker_panels <- list(
  "Neuroendocrine" = c("PHOX2B", "MYCN", "NXPH1", "SYT1", "DBH"),
  "Endothelial"    = c("EGFL7", "EMCN", "PLVAP", "PECAM1", "VWF"),
  "Schwann"        = c("CDH19", "PLP1", "PTPRZ1", "MPZ", "S100B"),
  "Fibroblast"     = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  "T_Cell"         = c("CD3D", "CD3E", "CD2", "CD8A", "CD4"),
  "NK_Cell"        = c("KLRF1", "KLRC1", "KLRD1", "NKG7", "GNLY"),
  "B_Cell"         = c("MS4A1", "CD79A", "CD79B", "VPREB3"),
  "Plasma"         = c("IGHG1", "IGHG2", "IGHG3", "JCHAIN"),
  "Myeloid"        = c("LYZ", "IL1B", "C1QC", "CD68", "CD14"),
  "pDC"            = c("LILRA4", "SCT", "PTCRA", "IRF7"),
  "RBC"            = c("HBA1", "HBA2", "HBB", "HBG1")
)

# Create feature plots for each panel
for (panel_name in names(marker_panels)) {
  genes <- marker_panels[[panel_name]]
  existing_genes <- genes[genes %in% rownames(nb.seurat)]

  if (length(existing_genes) == 0) {
    warning("No marker genes found for panel: ", panel_name)
    next
  }

  message("  Plotting markers for: ", panel_name)
  p <- FeaturePlot(nb.seurat, features = existing_genes,
                   min.cutoff = "q2", max.cutoff = "q98",
                   order = TRUE, ncol = 3)

  ggsave(file.path(PLOT_DIR, paste0("canonical_markers_", panel_name, ".pdf")),
         plot = p, width = 12, height = ceiling(length(existing_genes)/3) * 4)
}

## =============================================================================
## DOTPLOT OF CANONICAL MARKERS
## =============================================================================

message("\nCreating dotplot of canonical markers...")

all_markers <- unique(unlist(marker_panels))
existing_markers <- all_markers[all_markers %in% rownames(nb.seurat)]

if (length(existing_markers) > 0) {
  p_dot <- DotPlot(nb.seurat, features = existing_markers,
                   group.by = "seurat_clusters") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Canonical markers by cluster")

  ggsave(file.path(PLOT_DIR, "canonical_markers_dotplot.pdf"),
         plot = p_dot, width = 16, height = 8)
}

## =============================================================================
## FIND ALL MARKERS
## =============================================================================

message("\nRunning FindAllMarkers...")
message("  This may take some time depending on dataset size...")

# Set default assay to SCT for marker finding
if ("SCT" %in% names(nb.seurat@assays)) {
  DefaultAssay(nb.seurat) <- "SCT"
  message("  Using SCT assay for marker finding")
} else {
  DefaultAssay(nb.seurat) <- "RNA"
  message("  Warning: SCT assay not found, using RNA assay")
}

Idents(nb.seurat) <- "seurat_clusters"

FindAllMarkers_results <- FindAllMarkers(
  nb.seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE,
  assay = DefaultAssay(nb.seurat)
)

message("  Found ", nrow(FindAllMarkers_results), " marker genes across all clusters")

# Add a score column (pct.1 / pct.2 * avg_log2FC)
FindAllMarkers_results$score <-
  FindAllMarkers_results$pct.1 / (FindAllMarkers_results$pct.2 + 0.01) *
  FindAllMarkers_results$avg_log2FC

FindAllMarkers_results <- FindAllMarkers_results %>%
  arrange(cluster, desc(score))

## =============================================================================
## SAVE RESULTS
## =============================================================================

message("\nSaving marker results...")

# Save as RDS
markers_rds <- file.path(TABLE_DIR, "FindAllMarkers_results.rds")
saveRDS(FindAllMarkers_results, markers_rds)
message("  RDS: ", markers_rds)

# Save as CSV
markers_csv <- file.path(TABLE_DIR, "FindAllMarkers_results.csv")
write_csv(FindAllMarkers_results, markers_csv)
message("  CSV: ", markers_csv)

# Save as Excel (one sheet per cluster)
markers_by_cluster <- split(FindAllMarkers_results, FindAllMarkers_results$cluster)
markers_xlsx <- file.path(TABLE_DIR, "FindAllMarkers_by_cluster.xlsx")
write.xlsx(markers_by_cluster, markers_xlsx, overwrite = TRUE)
message("  Excel: ", markers_xlsx)

## =============================================================================
## TOP MARKERS SUMMARY
## =============================================================================

message("\nCreating top marker summaries...")

# Top 10 markers per cluster
top_markers <- FindAllMarkers_results %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 10) %>%
  ungroup()

top_markers_csv <- file.path(TABLE_DIR, "top10_markers_per_cluster.csv")
write_csv(top_markers, top_markers_csv)
message("  Top 10 markers per cluster: ", top_markers_csv)

# Create summary table
summary_table <- FindAllMarkers_results %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    top3_genes = paste(head(gene, 3), collapse = ", "),
    .groups = "drop"
  )

summary_csv <- file.path(TABLE_DIR, "cluster_marker_summary.csv")
write_csv(summary_table, summary_csv)
message("  Cluster summary: ", summary_csv)

# Print summary to console
message("\n=============================================================================")
message("MARKER SUMMARY BY CLUSTER")
message("=============================================================================")
for (i in seq_len(nrow(summary_table))) {
  message(sprintf("Cluster %s: %d markers | Top genes: %s",
                  summary_table$cluster[i],
                  summary_table$n_markers[i],
                  summary_table$top3_genes[i]))
}
message("=============================================================================\n")

## =============================================================================
## DOTPLOT OF TOP MARKERS
## =============================================================================

message("Creating dotplot of top marker genes per cluster...")

top_genes_for_plot <- FindAllMarkers_results %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 3) %>%
  pull(gene) %>%
  unique()

if (length(top_genes_for_plot) > 0) {
  p_top <- DotPlot(nb.seurat, features = top_genes_for_plot,
                   group.by = "seurat_clusters") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Top 3 marker genes per cluster")

  ggsave(file.path(PLOT_DIR, "top_markers_dotplot.pdf"),
         plot = p_top, width = max(12, length(top_genes_for_plot) * 0.4), height = 8)
}

## =============================================================================
## HEATMAP OF TOP MARKERS (SUBSAMPLED)
## =============================================================================

message("Creating heatmap of top markers...")

# Subsample cells for heatmap
n_cells_per_cluster <- 100
sampled_cells <- c()

for (clust in clusters) {
  cells_in_cluster <- Cells(nb.seurat)[nb.seurat$seurat_clusters == clust]
  if (length(cells_in_cluster) <= n_cells_per_cluster) {
    sampled_cells <- c(sampled_cells, cells_in_cluster)
  } else {
    sampled_cells <- c(sampled_cells, sample(cells_in_cluster, n_cells_per_cluster))
  }
}

nb.subset <- subset(nb.seurat, cells = sampled_cells)

# Use SCT assay for heatmap (already normalized)
if ("SCT" %in% names(nb.subset@assays)) {
  DefaultAssay(nb.subset) <- "SCT"
} else {
  DefaultAssay(nb.subset) <- "RNA"
  nb.subset <- NormalizeData(nb.subset, verbose = FALSE)
  nb.subset <- ScaleData(nb.subset, features = top_genes_for_plot, verbose = FALSE)
}

if (length(top_genes_for_plot) > 0) {
  p_heatmap <- DoHeatmap(nb.subset, features = top_genes_for_plot,
                         group.by = "seurat_clusters", size = 3,
                         assay = DefaultAssay(nb.subset)) +
    scale_fill_gradientn(colors = c("blue", "white", "red"))

  ggsave(file.path(PLOT_DIR, "top_markers_heatmap.pdf"),
         plot = p_heatmap, width = 12, height = max(8, length(top_genes_for_plot) * 0.2))
}

## =============================================================================
## COMPREHENSIVE DOTPLOT BY CELL TYPE
## =============================================================================

message("\nCreating comprehensive cell type marker dotplot...")

cell_type_markers <- list(
  "NE" = c("PHOX2B", "MYCN", "NXPH1", "SYT1", "DBH"),
  "Endothelial" = c("EGFL7", "EMCN", "PLVAP", "PECAM1", "VWF"),
  "Schwann" = c("CDH19", "PLP1", "PTPRZ1", "MPZ", "S100B"),
  "Fibroblast" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  "T cell" = c("CD3D", "CD3E", "CD2", "CD8A", "CD4"),
  "NK cell" = c("KLRF1", "KLRC1", "KLRD1", "NKG7", "GNLY"),
  "B cell" = c("MS4A1", "CD79A", "CD79B", "VPREB3"),
  "Plasma" = c("IGHG1", "IGHG2", "IGHG3", "JCHAIN"),
  "Myeloid" = c("LYZ", "IL1B", "C1QC", "CD68", "CD14"),
  "pDC" = c("LILRA4", "SCT", "PTCRA", "IRF7"),
  "RBCs" = c("HBA1", "HBA2", "HBB", "HBG1")
)

# Flatten the markers list
all_markers_ordered <- unlist(cell_type_markers, use.names = FALSE)
existing_markers_ordered <- all_markers_ordered[all_markers_ordered %in% rownames(nb.seurat)]

if (length(existing_markers_ordered) > 0) {
  # Create dotplot with all markers
  p_comprehensive <- DotPlot(
    nb.seurat,
    features = existing_markers_ordered,
    group.by = "seurat_clusters",
    cols = c("lightgrey", "red"),
    dot.scale = 8
  ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Marker",
      y = "Cluster",
      title = "Canonical Cell Type Markers by Cluster"
    ) +
    scale_size_continuous(name = "Percent\nExpressed", range = c(0, 10)) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2, title = "Average\nExpression")
    )

  # Save comprehensive dotplot
  ggsave(
    file.path(PLOT_DIR, "comprehensive_markers_dotplot.pdf"),
    plot = p_comprehensive,
    width = max(16, length(existing_markers_ordered) * 0.35),
    height = max(8, length(clusters) * 0.4)
  )

  message("  Saved comprehensive marker dotplot")

  # Extract data for cell-type grouped plot
  plot_data_list <- list()

  for (marker in existing_markers_ordered) {
    if (marker %in% rownames(nb.seurat)) {
      expr_data <- FetchData(nb.seurat, vars = c(marker, "seurat_clusters"))
      colnames(expr_data) <- c("expression", "cluster")

      summary_data <- expr_data %>%
        group_by(cluster) %>%
        summarise(
          pct_expressed = 100 * mean(expression > 0),
          avg_expression = mean(expression),
          .groups = "drop"
        ) %>%
        mutate(gene = marker)

      plot_data_list[[marker]] <- summary_data
    }
  }

  # Combine all data
  plot_data <- bind_rows(plot_data_list)

  # Add cell type annotations
  plot_data$cell_type <- NA
  for (ct_name in names(cell_type_markers)) {
    ct_genes <- cell_type_markers[[ct_name]]
    plot_data$cell_type[plot_data$gene %in% ct_genes] <- ct_name
  }

  # Create the plot with cell types as facets
  plot_data$gene <- factor(plot_data$gene, levels = existing_markers_ordered)
  plot_data$cell_type <- factor(plot_data$cell_type, levels = names(cell_type_markers))

  p_by_celltype <- ggplot(plot_data, aes(x = gene, y = cluster)) +
    geom_point(aes(size = pct_expressed, color = avg_expression)) +
    scale_size_continuous(name = "Percent\nExpressed", range = c(0, 8)) +
    scale_color_gradient(name = "Average\nExpression", low = "lightgrey", high = "red") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey95")
    ) +
    labs(
      x = "Marker",
      y = "Cluster",
      title = "Cell Type Markers by Cluster (Grouped by Cell Type)"
    ) +
    facet_grid(. ~ cell_type, scales = "free_x", space = "free_x")

  ggsave(
    file.path(PLOT_DIR, "markers_by_celltype_dotplot.pdf"),
    plot = p_by_celltype,
    width = max(20, length(existing_markers_ordered) * 0.4),
    height = max(8, length(clusters) * 0.5)
  )

  message("  Saved cell-type-grouped marker dotplot")
}

## =============================================================================
## COMPLETE
## =============================================================================

message("\n=============================================================================")
message("FindAllMarkers analysis complete!")
message("=============================================================================")
message("Next steps:")
message("  1. Review marker tables in: ", TABLE_DIR)
message("  2. Inspect plots in: ", PLOT_DIR)
message("  3. Create cluster_annotations.csv with your cell type labels")
message("  4. Run step1c_annotate_clusters.R to apply annotations")
message("=============================================================================\n")

message("Done!")
