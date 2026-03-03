#!/usr/bin/env Rscript
## =============================================================================
## Inspect Normalized SingleCellExperiment Object
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Verify normalized SCE object has all required components and metadata
## Input:   Normalized SCE object (RDS file)
## Output:  Console report with quality metrics and metadata verification
## =============================================================================

## =============================================================================
## CONFIGURATION
## =============================================================================

SCE_PATH <- "/path/to/NB_patient_analysis/results/NB_patients_SCE-norm_relaxed.RDS"

## =============================================================================

library(SingleCellExperiment)

## Load the object
sce <- readRDS(SCE_PATH)

## =============================================================================
## 1. BASIC DIMENSIONS
## =============================================================================
cat("\n=== DIMENSIONS ===\n")
cat("Cells: ", ncol(sce), "\n")
cat("Genes: ", nrow(sce), "\n")

## =============================================================================
## 2. ASSAYS (Raw and Normalized Counts)
## =============================================================================
cat("\n=== ASSAYS ===\n")
print(assayNames(sce))
cat("Expected: 'counts' and 'logcounts'\n")

## =============================================================================
## 3. CELL METADATA (colData)
## =============================================================================
cat("\n=== CELL METADATA COLUMNS ===\n")
print(colnames(colData(sce)))

cat("\n=== SAMPLE DISTRIBUTION ===\n")
if ("Sample" %in% colnames(colData(sce))) {
    print(table(sce$Sample))
}

cat("\n=== PATIENT DISTRIBUTION ===\n")
if ("Patient_ID" %in% colnames(colData(sce))) {
    print(table(sce$Patient_ID))
}

cat("\n=== FRACTION DISTRIBUTION ===\n")
if ("Fraction" %in% colnames(colData(sce))) {
    print(table(sce$Fraction))
}

cat("\n=== REPLICATE DISTRIBUTION ===\n")
if ("Replicate" %in% colnames(colData(sce))) {
    print(table(sce$Replicate))
}

## =============================================================================
## 4. CHECK SPECIFIC METADATA COLUMNS
## =============================================================================
cat("\n=== ESSENTIAL METADATA CHECK ===\n")
essential_meta <- c("Sample", "Patient_ID", "Tumor_Type", "Fraction",
                    "Replicate", "Sequencing_Run", "Barcode", "Batch", "sizeFactor")

for (col in essential_meta) {
    if (col %in% colnames(colData(sce))) {
        cat("✓ ", col, "\n", sep = "")
    } else {
        cat("✗ MISSING: ", col, "\n", sep = "")
    }
}

## =============================================================================
## 5. EXAMPLE METADATA (First 3 cells)
## =============================================================================
cat("\n=== EXAMPLE METADATA (first 3 cells) ===\n")
meta_cols <- c("Sample", "Patient_ID", "Fraction", "Replicate", "Barcode")
meta_cols_present <- meta_cols[meta_cols %in% colnames(colData(sce))]
if (length(meta_cols_present) > 0) {
    print(as.data.frame(colData(sce)[1:3, meta_cols_present]))
}

## =============================================================================
## 6. GENE METADATA (rowData)
## =============================================================================
cat("\n=== GENE METADATA COLUMNS ===\n")
if (ncol(rowData(sce)) > 0) {
    print(colnames(rowData(sce)))
    cat("Number of gene annotation columns: ", ncol(rowData(sce)), "\n")
} else {
    cat("No gene metadata (rowData is empty)\n")
}

cat("\n=== GENE ROWNAMES (first 10) ===\n")
print(head(rownames(sce), 10))

## =============================================================================
## 7. SIZE FACTORS
## =============================================================================
cat("\n=== SIZE FACTORS ===\n")
if (!is.null(sizeFactors(sce))) {
    cat("Size factors present: YES\n")
    cat("Summary:\n")
    print(summary(sizeFactors(sce)))
} else {
    cat("Size factors present: NO\n")
}

## =============================================================================
## 8. COUNTS SPARSITY
## =============================================================================
cat("\n=== DATA QUALITY ===\n")
cat("Raw counts sparsity: ",
    round(100 * sum(counts(sce) == 0) / length(counts(sce)), 2), "%\n")

## =============================================================================
## 9. SUMMARY TABLE: Cells per Patient × Fraction × Replicate
## =============================================================================
cat("\n=== CELLS PER PATIENT × FRACTION × REPLICATE ===\n")
if (all(c("Patient_ID", "Fraction", "Replicate") %in% colnames(colData(sce)))) {
    summary_table <- table(sce$Patient_ID, sce$Fraction, sce$Replicate)
    print(summary_table)
}

## =============================================================================
## 10. CHECK IF READY FOR SHARING
## =============================================================================
cat("\n=== READY FOR SHARING? ===\n")

ready <- TRUE

if (!all(c("counts", "logcounts") %in% assayNames(sce))) {
    cat("✗ Missing assays\n")
    ready <- FALSE
} else {
    cat("✓ Both raw and normalized counts present\n")
}

essential <- c("Sample", "Patient_ID", "Fraction", "Replicate")
missing_meta <- essential[!essential %in% colnames(colData(sce))]
if (length(missing_meta) > 0) {
    cat("✗ Missing metadata:", paste(missing_meta, collapse = ", "), "\n")
    ready <- FALSE
} else {
    cat("✓ Essential metadata present\n")
}

if (is.null(sizeFactors(sce))) {
    cat("⚠ Size factors missing (might be okay)\n")
} else {
    cat("✓ Size factors present\n")
}

if (ready) {
    cat("\n✅ OBJECT IS READY FOR SHARING!\n")
} else {
    cat("\n⚠️ OBJECT MAY BE MISSING SOME COMPONENTS\n")
}

## =============================================================================
## FINAL SUMMARY
## =============================================================================
cat("\n=== FINAL SUMMARY ===\n")
cat("Total cells: ", ncol(sce), "\n")
cat("Total genes: ", nrow(sce), "\n")
cat("Assays: ", paste(assayNames(sce), collapse = ", "), "\n")
cat("Metadata columns: ", ncol(colData(sce)), "\n")
