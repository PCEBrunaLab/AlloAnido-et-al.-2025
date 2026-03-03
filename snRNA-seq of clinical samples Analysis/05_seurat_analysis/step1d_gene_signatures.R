#!/usr/bin/env Rscript
## =============================================================================
## Neuroblastoma Patient Pipeline - Step 1d: Gene Signature Scoring
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Score neuroblastoma gene signatures (van Groningen, Dyer, Int)
## Input:   Annotated Seurat object with manual cell type labels
## Output:  Seurat object with signature scores, visualizations, tables
## =============================================================================

library(Seurat)
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(patchwork)

set.seed(12345)

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

# Project and data directories
PROJECT_DIR    <- "/path/to/NB_patient_analysis"
RESULTS_DIR    <- file.path(PROJECT_DIR, "results")
DATA_DIR       <- file.path(RESULTS_DIR, "data")
TABLE_DIR      <- file.path(RESULTS_DIR, "tables_step1d")
PLOT_DIR       <- file.path(RESULTS_DIR, "plots_step1d")
REF_DIR        <- file.path(PROJECT_DIR, "references")

# Input Seurat objects
SEURAT_ANNOTATED <- file.path(DATA_DIR, "NB_patients_seurat_annotated.rds")
SEURAT_HARMONY   <- file.path(DATA_DIR, "NB_patients_seurat_harmony.rds")
SEURAT_PATH      <- if (file.exists(SEURAT_ANNOTATED)) SEURAT_ANNOTATED else SEURAT_HARMONY
OUTPUT_PATH      <- file.path(DATA_DIR, "NB_patients_seurat_with_signatures.rds")

# Reference gene signatures
VANGRONINGEN_FILE <- file.path(REF_DIR, "vanGroningen_2017.xlsx")
DYER_FILE        <- file.path(REF_DIR, "Dyer_Sig.xlsx")
INT_SIGNATURE_CSV <- file.path(REF_DIR, "Int_signature.csv")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## LOAD SEURAT OBJECT
## =============================================================================

if (!file.exists(SEURAT_PATH)) {
  stop("Seurat object not found: ", SEURAT_PATH)
}

message("Loading Seurat object...")
nb.seurat <- readRDS(SEURAT_PATH)
message("  Cells: ", ncol(nb.seurat), ", Genes: ", nrow(nb.seurat))

has_manual_annot <- "CellType_Manual" %in% colnames(nb.seurat@meta.data)
group_var <- if (has_manual_annot) "CellType_Manual" else "seurat_clusters"

## =============================================================================
## LOAD VAN GRONINGEN SIGNATURES
## =============================================================================

message("\nLoading van Groningen signatures...")
vanGroningen_sigs <- list()

if (file.exists(VANGRONINGEN_FILE)) {
  vg_df <- read_xlsx(VANGRONINGEN_FILE, sheet = 1)

  if ("Gene" %in% names(vg_df) && "Cell_State" %in% names(vg_df)) {
    adrn_genes <- vg_df %>% filter(Cell_State == "ADRN") %>% pull(Gene) %>% na.omit() %>% as.character()
    mes_genes <- vg_df %>% filter(Cell_State == "MES") %>% pull(Gene) %>% na.omit() %>% as.character()

    if (length(adrn_genes) > 0) {
      vanGroningen_sigs[["ADRN"]] <- adrn_genes
      message("  Loaded ADRN: ", length(adrn_genes), " genes")
    }
    if (length(mes_genes) > 0) {
      vanGroningen_sigs[["MES"]] <- mes_genes
      message("  Loaded MES: ", length(mes_genes), " genes")
    }
  }
}

## =============================================================================
## LOAD DYER SIGNATURES
## =============================================================================

message("\nLoading Dyer signatures...")
Dyer_sigs <- list()

if (file.exists(DYER_FILE)) {
  dyer_df <- read_xlsx(DYER_FILE, sheet = 1)

  if ("Gene" %in% names(dyer_df) && "Cell_State" %in% names(dyer_df)) {
    adrn_genes <- dyer_df %>% filter(Cell_State == "ADRN") %>% pull(Gene) %>% na.omit() %>% as.character()
    mes_genes <- dyer_df %>% filter(Cell_State == "MES") %>% pull(Gene) %>% na.omit() %>% as.character()

    # Try to find Sympathoblast (multiple possible names)
    symp_patterns <- c("Sympathoblast", "SYMP", "Sympathetic", "sympathoblast")
    symp_genes <- NULL
    for (pattern in symp_patterns) {
      temp <- dyer_df %>% filter(Cell_State == pattern) %>% pull(Gene) %>% na.omit() %>% as.character()
      if (length(temp) > 0) {
        symp_genes <- temp
        break
      }
    }

    if (length(adrn_genes) > 0) {
      Dyer_sigs[["ADRN"]] <- adrn_genes
      message("  Loaded ADRN: ", length(adrn_genes), " genes")
    }
    if (length(mes_genes) > 0) {
      Dyer_sigs[["MES"]] <- mes_genes
      message("  Loaded MES: ", length(mes_genes), " genes")
    }
    if (!is.null(symp_genes) && length(symp_genes) > 0) {
      Dyer_sigs[["Sympathoblast"]] <- symp_genes
      message("  Loaded Sympathoblast: ", length(symp_genes), " genes")
    }
  }
}

## =============================================================================
## LOAD INT SIGNATURE
## =============================================================================

message("\nLoading Int signature...")
Int_sigs <- list()

if (file.exists(INT_SIGNATURE_CSV)) {
  int_df <- read_csv(INT_SIGNATURE_CSV, show_col_types = FALSE)
  gene_col <- if ("Gene" %in% colnames(int_df)) "Gene" else colnames(int_df)[1]
  int_genes <- unique(na.omit(as.character(int_df[[gene_col]])))

  # Filter to genes present in Seurat object
  int_genes <- intersect(int_genes, rownames(nb.seurat))

  if (length(int_genes) > 0) {
    Int_sigs[["Int"]] <- int_genes
    message("  Loaded Int: ", length(int_genes), " genes")
  }
}

## =============================================================================
## SCORE SIGNATURES
## =============================================================================

message("\nScoring signatures...")

# Set default assay to SCT
if ("SCT" %in% names(nb.seurat@assays)) {
  DefaultAssay(nb.seurat) <- "SCT"
  message("  Using SCT assay")
} else {
  DefaultAssay(nb.seurat) <- "RNA"
  message("  Warning: Using RNA assay (SCT not found)")
}

# Combine all signatures
all_sigs <- list()
if (length(vanGroningen_sigs) > 0) {
  for (sig_name in names(vanGroningen_sigs)) {
    all_sigs[[paste0("vanGroningen_", sig_name)]] <- vanGroningen_sigs[[sig_name]]
  }
}
if (length(Dyer_sigs) > 0) {
  for (sig_name in names(Dyer_sigs)) {
    all_sigs[[paste0("Dyer_", sig_name)]] <- Dyer_sigs[[sig_name]]
  }
}
if (length(Int_sigs) > 0) {
  for (sig_name in names(Int_sigs)) {
    all_sigs[[sig_name]] <- Int_sigs[[sig_name]]
  }
}

if (length(all_sigs) == 0) {
  stop("No signatures loaded. Check signature files in: ", REF_DIR)
}

# Score each signature
for (sig_name in names(all_sigs)) {
  message("  Scoring: ", sig_name, " (", length(all_sigs[[sig_name]]), " genes)")

  nb.seurat <- AddModuleScore(
    object = nb.seurat,
    features = list(all_sigs[[sig_name]]),
    name = paste0(sig_name, "_Score"),
    ctrl = 100,
    seed = 12345
  )

  # Rename column (AddModuleScore adds "1" suffix)
  old_name <- paste0(sig_name, "_Score1")
  if (old_name %in% colnames(nb.seurat@meta.data)) {
    colnames(nb.seurat@meta.data)[colnames(nb.seurat@meta.data) == old_name] <- sig_name
  }
}

sig_cols <- names(all_sigs)
message("  Scored ", length(sig_cols), " signatures")

## =============================================================================
## CREATE PLOTS
## =============================================================================

message("\nCreating plots...")

umap_reduction <- if ("umap" %in% names(nb.seurat@reductions)) "umap" else names(nb.seurat@reductions)[1]

# Plot van Groningen signatures
vg_names <- grep("^vanGroningen_", sig_cols, value = TRUE)
if (length(vg_names) > 0) {
  plot_list <- list()
  for (sig in vg_names) {
    color_scale <- if (grepl("ADRN", sig)) {
      scale_color_gradient(low = "lightgrey", high = "purple")
    } else {
      scale_color_gradient(low = "lightgrey", high = "orange")
    }
    plot_list[[sig]] <- FeaturePlot(nb.seurat, features = sig, reduction = umap_reduction, order = TRUE) +
      color_scale + ggtitle(sig)
  }
  ggsave(file.path(PLOT_DIR, "vanGroningen_signatures_UMAP.pdf"),
         wrap_plots(plot_list, ncol = 2), width = 12, height = 6 * ceiling(length(vg_names)/2))

  p_vln <- VlnPlot(nb.seurat, features = vg_names, group.by = group_var, pt.size = 0, ncol = length(vg_names))
  ggsave(file.path(PLOT_DIR, "vanGroningen_signatures_violin.pdf"), p_vln, width = 5 * length(vg_names), height = 6)
}

# Plot Dyer signatures
dy_names <- grep("^Dyer_", sig_cols, value = TRUE)
if (length(dy_names) > 0) {
  plot_list <- list()
  for (sig in dy_names) {
    color_scale <- if (grepl("ADRN", sig)) {
      scale_color_gradient(low = "lightgrey", high = "purple")
    } else if (grepl("MES", sig)) {
      scale_color_gradient(low = "lightgrey", high = "orange")
    } else {
      scale_color_gradient(low = "lightgrey", high = "darkgreen")
    }
    plot_list[[sig]] <- FeaturePlot(nb.seurat, features = sig, reduction = umap_reduction, order = TRUE) +
      color_scale + ggtitle(sig)
  }
  ggsave(file.path(PLOT_DIR, "Dyer_signatures_UMAP.pdf"),
         wrap_plots(plot_list, ncol = 2), width = 12, height = 6 * ceiling(length(dy_names)/2))

  p_vln <- VlnPlot(nb.seurat, features = dy_names, group.by = group_var, pt.size = 0, ncol = length(dy_names))
  ggsave(file.path(PLOT_DIR, "Dyer_signatures_violin.pdf"), p_vln, width = 5 * length(dy_names), height = 6)
}

# Plot Int signature
if ("Int" %in% sig_cols) {
  p_umap <- FeaturePlot(nb.seurat, features = "Int", reduction = umap_reduction, order = TRUE) +
    scale_color_gradient(low = "lightgrey", high = "darkblue") + ggtitle("Int Signature")
  ggsave(file.path(PLOT_DIR, "Int_signature_UMAP.pdf"), p_umap, width = 10, height = 8)

  p_vln <- VlnPlot(nb.seurat, features = "Int", group.by = group_var, pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Int Signature")
  ggsave(file.path(PLOT_DIR, "Int_signature_violin.pdf"), p_vln, width = 12, height = 6)
}

## =============================================================================
## SAVE TABLES
## =============================================================================

message("\nSaving tables...")

sig_scores <- nb.seurat@meta.data %>% mutate(Cell = colnames(nb.seurat))
base_cols <- c("Cell", "Patient", "Sample", group_var)
available_cols <- intersect(c(base_cols, sig_cols), colnames(sig_scores))
write_csv(sig_scores %>% select(all_of(available_cols)),
          file.path(TABLE_DIR, "signature_scores_per_cell.csv"))

sig_summary <- nb.seurat@meta.data %>%
  group_by(across(all_of(group_var))) %>%
  summarise(across(all_of(sig_cols), mean, na.rm = TRUE), .groups = "drop")
summary_filename <- if (has_manual_annot) "signature_scores_by_celltype.csv" else "signature_scores_by_cluster.csv"
write_csv(sig_summary, file.path(TABLE_DIR, summary_filename))

message("  Saved tables")

## =============================================================================
## SAVE OBJECT
## =============================================================================

message("\nSaving Seurat object...")
saveRDS(nb.seurat, OUTPUT_PATH)
message("  Saved: ", OUTPUT_PATH)

## =============================================================================
## SUMMARY
## =============================================================================

message("\n============================================================")
message("Signature scoring complete!")
message("============================================================")
message("Signatures scored: ", length(sig_cols))
for (sig_name in sig_cols) {
  message("  - ", sig_name)
}
message("\nOutputs:")
message("  - Seurat object: ", basename(OUTPUT_PATH))
message("  - Plots: ", PLOT_DIR)
message("  - Tables: ", TABLE_DIR)
message("============================================================\n")
