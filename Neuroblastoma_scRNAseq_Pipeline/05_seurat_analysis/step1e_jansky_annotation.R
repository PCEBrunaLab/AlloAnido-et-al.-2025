#!/usr/bin/env Rscript
## =============================================================================
## Neuroblastoma Patient Pipeline - Step 1e: Jansky Reference Annotation
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Annotate cells with Jansky developmental stages using SingleR
## Input:   Seurat object with signatures, Jansky reference data
## Output:  Annotated Seurat object with developmental stage predictions
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(RColorBrewer)
  library(patchwork)
})

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

# Project and data directories
PROJECT_DIR      <- "/path/to/NB_patient_analysis"
RESULTS_DIR      <- file.path(PROJECT_DIR, "results")
DATA_DIR         <- file.path(RESULTS_DIR, "data")
REF_DIR          <- file.path(PROJECT_DIR, "references")
OUTPUT_DIR       <- file.path(RESULTS_DIR, "step1e_jansky")

# Input Seurat object (from step 1d with signatures)
SEURAT_INPUT     <- file.path(DATA_DIR, "NB_patients_seurat_with_signatures.rds")

# Output subdirectories
PLOT_DIR         <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR        <- file.path(OUTPUT_DIR, "tables")
QC_DIR           <- file.path(OUTPUT_DIR, "qc_metrics")

# Analysis parameters
N_CORES                <- 8          # Number of cores for SingleR
DIMS                   <- 1:30       # Dimensions for label transfer
REDUCTION              <- "umap"     # Reduction for visualization
LOW_CONF_THRESHOLD     <- 0.1        # Delta score threshold for low-confidence

# Validate inputs
if (!file.exists(SEURAT_INPUT)) {
  stop(
    "Seurat object with signatures not found!\n",
    "Expected: ", SEURAT_INPUT, "\n",
    "Please run step1d_gene_signatures.R first."
  )
}

if (!dir.exists(REF_DIR)) {
  stop("Reference directory not found: ", REF_DIR)
}

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)

message("============================================================")
message("STEP 1e: Jansky developmental stage annotation")
message(" Input object  : ", SEURAT_INPUT)
message(" Reference dir : ", REF_DIR)
message(" Output dir    : ", OUTPUT_DIR)
message(" Cores         : ", N_CORES)
message(" Dimensions    : ", paste(range(DIMS), collapse = ":"))
message("============================================================")

## =============================================================================
## LOAD REFERENCE DATA
## =============================================================================

message("\n>>> STEP 1: Loading Jansky reference data...")

load_reference_data <- function(data_dir) {
  cat("  Loading reference files...\n")
  ref_file <- file.path(data_dir, "adrenal_medulla_Seurat.RDS")
  annot_file <- file.path(data_dir, "adrenal_medulla_annot.RDS")
  counts_file <- file.path(data_dir, "adrenal_medulla_counts.RDS")

  if (!file.exists(ref_file)) stop("Missing reference Seurat object at ", ref_file)
  if (!file.exists(annot_file)) stop("Missing reference annotation at ", annot_file)
  if (!file.exists(counts_file)) stop("Missing reference counts at ", counts_file)

  reference <- tryCatch(readRDS(ref_file), error = function(e) {
    message("WARNING: Unable to load reference Seurat object directly")
    NULL
  })
  annotation <- readRDS(annot_file)
  counts_ref <- readRDS(counts_file)

  cat("  Reference genes:", nrow(counts_ref), "\n")
  cell_n <- if (!is.null(reference@meta.data)) nrow(reference@meta.data) else ncol(counts_ref)
  cat("  Reference cells:", cell_n, "\n")

  list(reference = reference, annotation = annotation, counts = counts_ref)
}

ref_data <- load_reference_data(REF_DIR)

# Construct compatible Seurat reference object
counts_ref <- ref_data$counts
jansky_ref <- tryCatch({
  obj <- ref_data$reference
  if (!is.null(obj)) {
    suppressWarnings(
      CreateSeuratObject(
        counts = GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts"),
        project = "Jansky_reference"
      )
    )
  } else {
    NULL
  }
}, error = function(e) NULL)

if (is.null(jansky_ref)) {
  jansky_ref <- CreateSeuratObject(counts = counts_ref, project = "Jansky_reference")
}

annotation_df <- ref_data$annotation

# Extract reference labels
extract_reference_labels <- function(annotation_df, cell_names, preferred_column = NULL) {
  labels <- NULL
  column_used <- NULL

  if (is.vector(annotation_df) && !is.list(annotation_df)) {
    if (!is.null(names(annotation_df))) {
      labels <- annotation_df[cell_names]
      column_used <- "vector"
    } else if (length(annotation_df) == length(cell_names)) {
      labels <- annotation_df
      names(labels) <- cell_names
      column_used <- "vector"
    }
  } else if (is.data.frame(annotation_df)) {
    df <- annotation_df
    if (!is.null(df$cell)) rownames(df) <- df$cell
    if (!is.null(df$Cell)) rownames(df) <- df$Cell
    if (!is.null(df$barcode)) rownames(df) <- df$barcode
    if (!is.null(df$Barcode)) rownames(df) <- df$Barcode

    candidate_cols <- setdiff(colnames(df), c("cell", "Cell", "barcode", "Barcode"))
    if (length(candidate_cols) > 0) {
      choose_col <- if (!is.null(preferred_column) && preferred_column %in% candidate_cols) {
        preferred_column
      } else {
        candidate_cols[1]
      }
      column_used <- choose_col
      df <- df[intersect(rownames(df), cell_names), , drop = FALSE]
      labels <- df[cell_names, column_used, drop = TRUE]
    }
  }

  list(labels = labels, column = column_used)
}

# Extract labels from reference
if (!is.null(ref_data$reference) && !is.null(ref_data$reference@active.ident)) {
  ref_labels <- as.character(ref_data$reference@active.ident)
  names(ref_labels) <- names(ref_data$reference@active.ident)
  ref_labels <- ref_labels[colnames(jansky_ref)]
  message("  Using labels from reference@active.ident")
} else {
  label_info <- extract_reference_labels(annotation_df, colnames(jansky_ref), NULL)
  ref_labels <- label_info$labels
}

if (is.null(ref_labels) || all(is.na(ref_labels))) {
  stop("Unable to derive reference labels; please check reference structure.")
}

ref_labels <- as.character(ref_labels)
jansky_ref$reference_annotation <- ref_labels
Idents(jansky_ref) <- "reference_annotation"

cat("  Reference developmental stages loaded\n")

## =============================================================================
## LOAD AND PREPARE QUERY DATA
## =============================================================================

message("\n>>> STEP 2: Loading query data...")

nb <- readRDS(SEURAT_INPUT)

# Set SCT as default assay
if ("SCT" %in% names(nb@assays)) {
  DefaultAssay(nb) <- "SCT"
  cat("  Using SCT assay as default\n")
} else if ("RNA" %in% names(nb@assays)) {
  DefaultAssay(nb) <- "RNA"
  cat("  Using RNA assay as default\n")
} else {
  stop("Seurat object does not contain SCT or RNA assay.")
}

if (!"RNA" %in% names(nb@assays)) {
  stop("Seurat object does not contain RNA assay required for SingleR.")
}

cat("  Query cells:", ncol(nb), "\n")
cat("  Query genes:", nrow(nb), "\n")

## =============================================================================
## PREPROCESSING
## =============================================================================

message("\n>>> STEP 3: Preprocessing reference and query...")

# Prepare reference
if (!"RNA" %in% names(jansky_ref@assays)) {
  DefaultAssay(jansky_ref) <- names(jansky_ref@assays)[1]
} else {
  DefaultAssay(jansky_ref) <- "RNA"
}

cat("  Normalizing reference...\n")
jansky_ref <- NormalizeData(jansky_ref)
jansky_ref <- FindVariableFeatures(jansky_ref, selection.method = "vst", nfeatures = 3000)
jansky_ref <- ScaleData(jansky_ref, features = VariableFeatures(jansky_ref), verbose = FALSE)
jansky_ref <- RunPCA(jansky_ref, features = VariableFeatures(jansky_ref), verbose = FALSE)

# Prepare query
cat("  Checking query normalization status...\n")
current_assay <- DefaultAssay(nb)
DefaultAssay(nb) <- "RNA"

if (!"data" %in% slotNames(nb@assays$RNA) || length(nb@assays$RNA@data) == 0) {
  cat("    Normalizing RNA assay for SingleR...\n")
  nb <- NormalizeData(nb, assay = "RNA", verbose = FALSE)
} else {
  cat("    RNA assay already normalized\n")
}

if (length(VariableFeatures(nb, assay = "RNA")) == 0) {
  cat("    Finding variable features in RNA assay...\n")
  nb <- FindVariableFeatures(nb, assay = "RNA", selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

if (!"pca" %in% names(nb@reductions)) {
  cat("    Running PCA on SCT assay...\n")
  DefaultAssay(nb) <- "SCT"
  if (length(VariableFeatures(nb, assay = "SCT")) == 0) {
    nb <- FindVariableFeatures(nb, assay = "SCT", selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  }
  nb <- ScaleData(nb, assay = "SCT", features = VariableFeatures(nb, assay = "SCT"), verbose = FALSE)
  nb <- RunPCA(nb, assay = "SCT", features = VariableFeatures(nb, assay = "SCT"), verbose = FALSE)
}

# Set SCT as default
if ("SCT" %in% names(nb@assays)) {
  DefaultAssay(nb) <- "SCT"
  cat("  SCT set as default assay\n")
} else {
  DefaultAssay(nb) <- "RNA"
}

bp <- MulticoreParam(workers = N_CORES, progressbar = TRUE)

cat("  Preprocessing complete\n")

## =============================================================================
## RUN SINGLER ANNOTATION
## =============================================================================

message("\n>>> STEP 4: Running SingleR annotation...")

pred_singleR <- SingleR(
  test = GetAssayData(nb, assay = "RNA", slot = "data"),
  assay.type.test = "logcounts",
  ref = GetAssayData(jansky_ref, assay = "RNA", slot = "data"),
  assay.type.ref = "logcounts",
  labels = jansky_ref$reference_annotation,
  BPPARAM = bp
)

# Add predictions to Seurat object
singleR_labels <- pred_singleR$labels[match(colnames(nb), rownames(pred_singleR))]
singleR_scores <- apply(pred_singleR$scores, 1, max)[match(colnames(nb), rownames(pred_singleR))]
nb$jansky_singleR_label <- singleR_labels
nb$jansky_singleR_score <- singleR_scores

# Add pruned labels
singleR_pruned <- pred_singleR$pruned.labels[match(colnames(nb), rownames(pred_singleR))]
nb$jansky_singleR_pruned <- singleR_pruned

# Set SCT back as default
if ("SCT" %in% names(nb@assays)) {
  DefaultAssay(nb) <- "SCT"
}

message("  SingleR label distribution:")
print(table(nb$jansky_singleR_label, useNA = "ifany"))

# Calculate delta scores
pred_singleR_scores_sorted <- t(apply(pred_singleR$scores, 1, function(x) sort(x, decreasing = TRUE)))
delta_scores <- pred_singleR_scores_sorted[,1] - pred_singleR_scores_sorted[,2]
nb$jansky_singleR_delta <- delta_scores[match(colnames(nb), rownames(pred_singleR))]

# Flag low-confidence cells
nb$jansky_low_confidence <- nb$jansky_singleR_delta < LOW_CONF_THRESHOLD
n_low_conf <- sum(nb$jansky_low_confidence, na.rm = TRUE)
pct_low_conf <- n_low_conf / length(nb$jansky_low_confidence) * 100

cat("  SingleR annotation complete\n")

## =============================================================================
## RUN SEURAT LABEL TRANSFER
## =============================================================================

message("\n>>> STEP 5: Running Seurat label transfer...")

# Find common features
if ("SCT" %in% names(nb@assays) && length(VariableFeatures(nb, assay = "SCT")) > 0) {
  nb_var_features <- VariableFeatures(nb, assay = "SCT")
} else {
  nb_var_features <- VariableFeatures(nb, assay = "RNA")
}

common_features <- intersect(VariableFeatures(jansky_ref), nb_var_features)
if (length(common_features) < 200) {
  common_features <- intersect(rownames(jansky_ref), rownames(nb))
}

cat("  Using", length(common_features), "common features\n")

# Label transfer
transfer_assay <- if ("SCT" %in% names(nb@assays)) "SCT" else "RNA"
DefaultAssay(nb) <- transfer_assay

anchors <- FindTransferAnchors(
  reference = jansky_ref,
  query = nb,
  features = common_features,
  dims = DIMS,
  reference.reduction = "pca",
  normalization.method = if (transfer_assay == "SCT") "SCT" else "LogNormalize"
)

transfer_pred <- TransferData(
  anchorset = anchors,
  refdata = jansky_ref$reference_annotation,
  dims = DIMS
)

nb <- AddMetaData(nb, metadata = transfer_pred)
nb$jansky_transfer_label <- nb$predicted.id

# Set SCT as default
if ("SCT" %in% names(nb@assays)) {
  DefaultAssay(nb) <- "SCT"
}

cat("  Label transfer complete\n")

## =============================================================================
## METHOD COMPARISON
## =============================================================================

message("\n>>> STEP 6: Comparing SingleR vs Seurat Transfer...")

agreement <- sum(nb$jansky_singleR_label == nb$predicted.id, na.rm = TRUE)
total_compared <- sum(!is.na(nb$jansky_singleR_label) & !is.na(nb$predicted.id))
pct_agree <- if (total_compared > 0) agreement / total_compared * 100 else 0

message(sprintf("  Agreement between methods: %.1f%% (%d / %d cells)",
                pct_agree, agreement, total_compared))

cat("  Method comparison complete\n")

## =============================================================================
## SAVE RESULTS
## =============================================================================

message("\n>>> STEP 7: Saving results...")

# Ensure SCT is default
if ("SCT" %in% names(nb@assays)) {
  DefaultAssay(nb) <- "SCT"
}

annotated_rds <- file.path(OUTPUT_DIR, "NB_patients_seurat_with_jansky.rds")
main_output <- file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds")

saveRDS(nb, annotated_rds)
saveRDS(nb, main_output)

# Export metadata
meta_cols_to_export <- c(
  "Patient", "Patient_ID", "Sample", "Batch", "Fraction", "Replicate",
  "seurat_clusters", "CellType_Manual", "manual_annotation",
  "jansky_singleR_label", "jansky_singleR_score", "jansky_singleR_pruned",
  "jansky_singleR_delta", "jansky_low_confidence",
  "predicted.id", "prediction.score.max", "jansky_transfer_label",
  "vanGroningen_ADRN", "vanGroningen_MES",
  "Dyer_ADRN", "Dyer_MES", "Dyer_Sympathoblast",
  "detected", "total", "subsets_mt_percent", "subsets_ribo_percent"
)

meta_cols_present <- intersect(meta_cols_to_export, colnames(nb@meta.data))
meta_out <- nb@meta.data %>% select(all_of(meta_cols_present))

write.csv(
  meta_out,
  file.path(TABLE_DIR, "NB_patients_jansky_annotations.csv"),
  row.names = TRUE
)

cat("  Results saved\n")

## =============================================================================
## GENERATE VISUALIZATIONS
## =============================================================================

message("\n>>> STEP 8: Generating visualizations...")

# Bar plots
pdf(file.path(PLOT_DIR, "01_barplot_singleR_distribution.pdf"), width = 10, height = 6)
print(
  ggplot(as.data.frame(table(nb$jansky_singleR_label)), aes(y = reorder(Var1, Freq), x = Freq, fill = Var1)) +
    geom_col() +
    labs(x = "Number of cells", y = "Developmental Stage",
         title = "SingleR predicted Jansky developmental stages") +
    theme_minimal() +
    theme(legend.position = "none")
)
dev.off()

# UMAP plots
if (REDUCTION %in% names(nb@reductions)) {
  pdf(file.path(PLOT_DIR, "06_umap_singleR_labels.pdf"), width = 12, height = 8)
  print(
    DimPlot(nb, reduction = REDUCTION, group.by = "jansky_singleR_label",
            label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3) +
      ggtitle("Jansky developmental stages (SingleR)") +
      theme(legend.text = element_text(size = 8))
  )
  dev.off()

  pdf(file.path(PLOT_DIR, "09_umap_transfer_labels.pdf"), width = 12, height = 8)
  print(
    DimPlot(nb, reduction = REDUCTION, group.by = "predicted.id",
            label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3) +
      ggtitle("Jansky developmental stages (Seurat Transfer)") +
      theme(legend.text = element_text(size = 8))
  )
  dev.off()
}

cat("  Visualizations complete\n")

## =============================================================================
## SUMMARY STATISTICS
## =============================================================================

message("\n>>> STEP 9: Generating summary statistics...")

summary_stats <- data.frame(
  Cell_Type = names(table(pred_singleR$labels)),
  Count = as.numeric(table(pred_singleR$labels)),
  Percentage = as.numeric(prop.table(table(pred_singleR$labels)) * 100)
) %>% arrange(desc(Count))

write.csv(summary_stats, file.path(TABLE_DIR, "singleR_prediction_summary.csv"), row.names = FALSE)

cat("  Summary statistics generated\n")

## =============================================================================
## FINAL SUMMARY
## =============================================================================

message("\n============================================================")
message("ANNOTATION COMPLETE!")
message("============================================================")
cat(sprintf("Total cells analysed: %d\n", ncol(nb)))
cat(sprintf("Unique SingleR labels: %d\n", length(unique(na.omit(nb$jansky_singleR_label)))))
cat(sprintf("Unique transfer labels: %d\n", length(unique(na.omit(nb$predicted.id)))))
cat(sprintf("Low confidence cells: %.1f%%\n", pct_low_conf))
cat(sprintf("Method agreement: %.1f%%\n", pct_agree))

cat("\nAnnotated object saved to:\n")
cat("  ", main_output, "\n")
cat("  ^ Use this for downstream analysis\n")

message("============================================================\n")
