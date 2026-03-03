#!/usr/bin/env Rscript
## Longitudinal Analysis for NB Patient Samples (F1 vs F2)
## Author: Ayeh Sadr | Pipeline: Neuroblastoma Patient scRNA-seq Analysis

library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(ComplexHeatmap); library(circlize)
set.seed(12345)

## CONFIGURATION
PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")
OUTPUT_DIR <- file.path(RESULTS_DIR, "longitudinal_analysis")
SEURAT_FILE <- file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds")

MIN_CELLS_PER_GROUP <- 50
LOGFC_THRESHOLD <- 0.25
PADJ_THRESHOLD <- 0.05

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "DE_genes"), recursive = TRUE, showWarnings = FALSE)

message("LONGITUDINAL ANALYSIS - NB PATIENT SAMPLES (F1 vs F2)")
message("Input: ", SEURAT_FILE); message("Output: ", OUTPUT_DIR); message("")

## LOAD DATA
message("[1/8] Loading Seurat object...")
seurat_obj <- readRDS(SEURAT_FILE)
message("Loaded: ", ncol(seurat_obj), " cells x ", nrow(seurat_obj), " genes")

if ("Patient_ID" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$Patient_ID <- NULL
  message("Removed incorrect 'Patient_ID' column")
}

## IDENTIFY LONGITUDINAL PATIENTS
message("\n[2/8] Identifying patients with longitudinal samples...")
longitudinal_summary <- seurat_obj@meta.data %>%
  group_by(Patient, Fraction) %>% summarise(N_cells = n(), .groups = "drop") %>%
  group_by(Patient) %>%
  summarise(has_F1 = "F1" %in% Fraction, has_F2 = "F2" %in% Fraction, n_timepoints = n(), .groups = "drop") %>%
  filter(has_F1 & has_F2)

longitudinal_patients <- longitudinal_summary$Patient
seurat_long <- subset(seurat_obj, subset = Patient %in% longitudinal_patients)
seurat_long$Patient_Fraction <- paste(seurat_long$Patient, seurat_long$Fraction, sep = "_")
message("Found ", length(longitudinal_patients), " patients with longitudinal data")

## CELL COMPOSITION
message("\n[3/8] Analyzing cell composition changes...")
composition <- seurat_long@meta.data %>%
  group_by(Patient, Fraction, seurat_clusters) %>% summarise(N_cells = n(), .groups = "drop") %>%
  group_by(Patient, Fraction) %>% mutate(Proportion = N_cells / sum(N_cells)) %>% ungroup()
write.csv(composition, file.path(OUTPUT_DIR, "tables", "cell_composition_by_timepoint.csv"), row.names = FALSE)

p1 <- ggplot(composition, aes(x = Fraction, y = Proportion, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "fill") + facet_wrap(~Patient, nrow = 1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
  labs(title = "Cell Type Composition by Timepoint", x = "Timepoint (Fraction)", y = "Proportion", fill = "Cluster") +
  scale_y_continuous(labels = scales::percent)
ggsave(file.path(OUTPUT_DIR, "plots", "01_cell_composition_by_timepoint.pdf"), p1, width = 12, height = 6)
message("Saved: 01_cell_composition_by_timepoint.pdf")

## UMAP VISUALIZATIONS
message("\n[4/8] Generating UMAP visualizations...")
for (patient in longitudinal_patients) {
  patient_obj <- subset(seurat_long, subset = Patient == patient)
  if (ncol(patient_obj) == 0) { message("Patient ", patient, ": No cells, skipping"); next }
  message("Patient ", patient, ": Creating UMAP plots...")
  
  p_timepoint <- DimPlot(patient_obj, group.by = "Fraction", reduction = "umap", pt.size = 0.5,
                         cols = c("F1" = "blue", "F2" = "red")) +
    ggtitle(paste0("Patient ", patient, ": F1 vs F2")) + theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  ggsave(file.path(OUTPUT_DIR, "plots", paste0("umap_timepoint_patient_", patient, ".pdf")), p_timepoint, width = 8, height = 7)
  
  p_clusters <- DimPlot(patient_obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE, pt.size = 0.5) +
    ggtitle(paste0("Patient ", patient, ": Clusters")) + theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))
  ggsave(file.path(OUTPUT_DIR, "plots", paste0("umap_clusters_patient_", patient, ".pdf")), p_clusters, width = 8, height = 7)
  
  p_split <- DimPlot(patient_obj, group.by = "seurat_clusters", split.by = "Fraction", reduction = "umap", label = TRUE, pt.size = 0.4) +
    ggtitle(paste0("Patient ", patient, ": Clusters by Timepoint")) + theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))
  ggsave(file.path(OUTPUT_DIR, "plots", paste0("umap_split_timepoint_patient_", patient, ".pdf")), p_split, width = 14, height = 6)
  
  message("  ✓ Saved 3 UMAP plots for patient ", patient)
}

## DIFFERENTIAL EXPRESSION
message("\n[5/8] Running differential expression analysis (F1 vs F2)...")
run_DE_F1_vs_F2 <- function(seurat_obj, patient_id) {
  patient_obj <- subset(seurat_obj, subset = Patient == patient_id)
  n_f1 <- sum(patient_obj$Fraction == "F1")
  n_f2 <- sum(patient_obj$Fraction == "F2")
  if (n_f1 < MIN_CELLS_PER_GROUP | n_f2 < MIN_CELLS_PER_GROUP) {
    message("Patient ", patient_id, ": Skipping (F1=", n_f1, ", F2=", n_f2, ")")
    return(NULL)
  }
  message("Patient ", patient_id, ": Running DE (F1=", n_f1, ", F2=", n_f2, " cells)")
  if ("SCT" %in% names(patient_obj@assays)) {
    DefaultAssay(patient_obj) <- "SCT"
  } else if ("RNA" %in% names(patient_obj@assays)) {
    DefaultAssay(patient_obj) <- "RNA"
  } else {
    return(NULL)
  }
  Idents(patient_obj) <- "Fraction"
  de_genes <- tryCatch({
    FindMarkers(patient_obj, ident.1 = "F2", ident.2 = "F1", test.use = "wilcox",
                logfc.threshold = 0, min.pct = 0.1, verbose = FALSE)
  }, error = function(e) { return(NULL) })
  if (is.null(de_genes) || nrow(de_genes) == 0) return(NULL)
  de_genes$gene <- rownames(de_genes)
  de_genes$Patient <- patient_id
  de_genes$Significant <- de_genes$p_val_adj < PADJ_THRESHOLD & abs(de_genes$avg_log2FC) > LOGFC_THRESHOLD
  return(de_genes)
}

all_de_results <- list()
for (patient in longitudinal_patients) {
  de_result <- run_DE_F1_vs_F2(seurat_long, patient)
  if (!is.null(de_result)) {
    all_de_results[[patient]] <- de_result
    write.csv(de_result, file.path(OUTPUT_DIR, "DE_genes", paste0("DE_F2_vs_F1_patient_", patient, ".csv")), row.names = FALSE)
  }
}
message("Completed DE analysis for ", length(all_de_results), " patients")

if (length(all_de_results) > 0) {
  all_de_combined <- do.call(rbind, all_de_results)
  write.csv(all_de_combined, file.path(OUTPUT_DIR, "tables", "DE_F2_vs_F1_all_patients.csv"), row.names = FALSE)
}

## HEATMAPS
message("\n[6/8] Generating heatmaps...")
n_heatmaps_generated <- 0
if (length(all_de_results) > 0) {
  expr_data <- tryCatch({ LayerData(seurat_long, assay = "SCT", layer = "data") },
                        error = function(e) { GetAssayData(seurat_long, assay = "SCT", slot = "data") })
  message("Expression data loaded: ", nrow(expr_data), " genes x ", ncol(expr_data), " cells")
  
  for (patient in names(all_de_results)) {
    de_data <- all_de_results[[patient]]
    if (is.null(de_data) || nrow(de_data) == 0) next
    
    top_genes <- de_data %>% filter(Significant == TRUE) %>% arrange(p_val_adj, desc(abs(avg_log2FC))) %>% head(50) %>% pull(gene)
    if (length(top_genes) < 5) next
    
    top_genes <- intersect(top_genes, rownames(expr_data))
    if (length(top_genes) < 5) next
    
    message("Patient ", patient, ": Creating heatmap with ", length(top_genes), " genes...")
    f1_cells <- which(seurat_long$Patient == patient & seurat_long$Fraction == "F1")
    f2_cells <- which(seurat_long$Patient == patient & seurat_long$Fraction == "F2")
    
    expr_matrix <- cbind(rowMeans(expr_data[top_genes, f1_cells, drop = FALSE]),
                        rowMeans(expr_data[top_genes, f2_cells, drop = FALSE]))
    colnames(expr_matrix) <- c("F1", "F2")
    rownames(expr_matrix) <- top_genes
    
    pdf_file <- file.path(OUTPUT_DIR, "plots", paste0("heatmap_patient_", patient, ".pdf"))
    pdf(pdf_file, width = 6, height = 10)
    col_fun <- colorRamp2(c(min(expr_matrix), median(expr_matrix), max(expr_matrix)), c("blue", "white", "red"))
    print(Heatmap(expr_matrix, name = "Expression", col = col_fun, cluster_rows = TRUE, cluster_columns = FALSE,
                  show_row_names = TRUE, show_column_names = TRUE,
                  column_title = paste0("Patient ", patient, ": Top DE Genes (F1 vs F2)"), row_names_gp = gpar(fontsize = 8)))
    dev.off()
    n_heatmaps_generated <- n_heatmaps_generated + 1
    message("  ✓ Saved: ", pdf_file)
  }
  message("Total heatmaps generated: ", n_heatmaps_generated)
}

## VOLCANO PLOTS
message("\n[7/8] Generating volcano plots...")
n_volcano_generated <- 0
if (length(all_de_results) > 0) {
  for (patient in names(all_de_results)) {
    de_data <- all_de_results[[patient]]
    if (is.null(de_data) || nrow(de_data) == 0) next
    
    n_sig <- sum(de_data$Significant, na.rm = TRUE)
    n_up <- sum(de_data$Significant & de_data$avg_log2FC > 0, na.rm = TRUE)
    n_down <- sum(de_data$Significant & de_data$avg_log2FC < 0, na.rm = TRUE)
    message("Patient ", patient, ": Creating volcano plot (", n_sig, " significant)")
    
    de_data$Change <- ifelse(de_data$Significant & de_data$avg_log2FC > 0, "Up in F2",
                            ifelse(de_data$Significant & de_data$avg_log2FC < 0, "Down in F2", "Not significant"))
    
    p <- ggplot(de_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Change)) +
      geom_point(alpha = 0.5, size = 1.5) + geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD), linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("Up in F2" = "red", "Down in F2" = "blue", "Not significant" = "gray80")) + theme_classic() +
      labs(title = paste0("Patient ", patient, ": F2 vs F1"), x = "Log2 Fold Change (F2/F1)", y = "-log10(Adjusted p-value)", color = "") +
      theme(legend.position = "top", plot.title = element_text(face = "bold", size = 14))
    
    pdf_file <- file.path(OUTPUT_DIR, "plots", paste0("volcano_patient_", patient, ".pdf"))
    ggsave(pdf_file, p, width = 8, height = 7)
    n_volcano_generated <- n_volcano_generated + 1
    message("  ✓ Saved: ", pdf_file)
  }
  message("Total volcano plots generated: ", n_volcano_generated)
}

## SUMMARY
message("\n============================================================")
message("LONGITUDINAL ANALYSIS COMPLETE")
message("============================================================")
message("Longitudinal patients: ", length(longitudinal_patients))
message("Total cells analyzed: ", ncol(seurat_long))
if (length(all_de_results) > 0) {
  message("Patients with DE results: ", length(all_de_results))
  message("Heatmaps generated: ", n_heatmaps_generated, " / ", length(all_de_results))
  message("Volcano plots generated: ", n_volcano_generated, " / ", length(all_de_results))
}
message("Output files: ", OUTPUT_DIR)
message("============================================================\n")
