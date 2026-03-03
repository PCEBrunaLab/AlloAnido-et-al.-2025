#!/usr/bin/env Rscript
## =============================================================================
## Verify Final Seurat Object - Check All Pipeline Information Is Present
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Validate final Seurat object contains all required QC metrics,
##          annotations, and metadata for downstream analysis
## Input:   Final Seurat object (RDS file)
## Output:  Console report with comprehensive verification checklist
## =============================================================================

## =============================================================================
## CONFIGURATION
## =============================================================================

DATA_DIR <- "/path/to/NB_patient_analysis/results/data"
OBJECT_NAME_1 <- "NB_patients_seurat_with_jansky.rds"
OBJECT_NAME_2 <- "NB_patients_seurat_with_signatures.rds"

## =============================================================================

library(Seurat)

## Path to final object
FINAL_OBJECT <- file.path(DATA_DIR, OBJECT_NAME_1)

if (!file.exists(FINAL_OBJECT)) {
  FINAL_OBJECT <- file.path(DATA_DIR, OBJECT_NAME_2)
}

if (!file.exists(FINAL_OBJECT)) {
  stop("Final object not found! Expected:\n",
       "  - ", file.path(DATA_DIR, OBJECT_NAME_1), "\n",
       "  - ", file.path(DATA_DIR, OBJECT_NAME_2))
}

cat("============================================================\n")
cat("VERIFYING FINAL SEURAT OBJECT\n")
cat("============================================================\n")
cat("Loading: ", FINAL_OBJECT, "\n\n")

nb <- readRDS(FINAL_OBJECT)

## Basic info
cat("BASIC INFO:\n")
cat("  Cells: ", ncol(nb), "\n")
cat("  Genes: ", nrow(nb), "\n")
cat("  Assays: ", paste(names(nb@assays), collapse = ", "), "\n")
cat("  Reductions: ", paste(names(nb@reductions), collapse = ", "), "\n")
cat("  Metadata columns: ", ncol(nb@meta.data), "\n\n")

## Check metadata columns
cat("METADATA CHECK:\n")
meta_cols <- colnames(nb@meta.data)

checks <- list(
  "Patient info" = any(c("Patient", "Patient_ID") %in% meta_cols),
  "Sample info" = "Sample" %in% meta_cols,
  "Fraction/Replicate" = all(c("Fraction", "Replicate") %in% meta_cols),
  "Seurat clusters" = "seurat_clusters" %in% meta_cols,
  "QC metrics" = all(c("detected", "total", "subsets_mt_percent") %in% meta_cols)
)

checks[["Manual annotations"]] <- "CellType_Manual" %in% meta_cols

signature_cols <- c("vanGroningen_ADRN_score", "vanGroningen_MES_score",
                    "Dyer_ADRN_score", "Dyer_MES_score", "Dyer_Sympathoblast_score")
checks[["Gene signatures"]] <- sum(signature_cols %in% meta_cols) >= 3

jansky_cols <- c("jansky_singleR_label", "jansky_transfer_label", "predicted.id")
checks[["Jansky annotations"]] <- any(jansky_cols %in% meta_cols)

for (check_name in names(checks)) {
  status <- if (checks[[check_name]]) "✓" else "✗"
  cat("  ", status, " ", check_name, "\n")
}

## Check reductions
cat("\nREDUCTIONS CHECK:\n")
reductions_needed <- c("harmony", "umap", "pca")
for (red in reductions_needed) {
  found <- any(grepl(red, names(nb@reductions), ignore.case = TRUE))
  status <- if (found) "✓" else "✗"
  cat("  ", status, " ", toupper(red), "\n")
}

## Summary
all_ok <- all(unlist(checks)) &&
          any(grepl("harmony", names(nb@reductions), ignore.case = TRUE)) &&
          any(grepl("umap", names(nb@reductions), ignore.case = TRUE))

cat("\n============================================================\n")
if (all_ok) {
  cat("✓ ALL CHECKS PASSED - Object is ready for downstream analysis!\n")
} else {
  cat("⚠ SOME CHECKS FAILED - Review missing components above\n")
}
cat("============================================================\n")

cat("\nOBJECT LOCATION:\n")
cat("  ", FINAL_OBJECT, "\n\n")

cat("METADATA COLUMNS SUMMARY:\n")
cat("  Total columns: ", length(meta_cols), "\n")
if ("CellType_Manual" %in% meta_cols) {
  cat("  Cell types: ", length(unique(nb$CellType_Manual)), "\n")
  cat("    ", paste(unique(nb$CellType_Manual), collapse = ", "), "\n")
}
if (any(jansky_cols %in% meta_cols)) {
  jansky_col <- jansky_cols[jansky_cols %in% meta_cols][1]
  cat("  Jansky stages: ", length(unique(nb@meta.data[[jansky_col]])), "\n")
  cat("    ", paste(unique(nb@meta.data[[jansky_col]]), collapse = ", "), "\n")
}

cat("\nDone!\n")
