#!/usr/bin/env Rscript
## Technical Replicate Analysis - Assess variation between replicates (i vs ii)
## Author: Ayeh Sadr | Pipeline: Neuroblastoma Patient scRNA-seq Analysis

library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(corrplot)
set.seed(12345)

## CONFIGURATION
PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")
OUTPUT_DIR <- file.path(RESULTS_DIR, "replicate_analysis")
SEURAT_FILE <- file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

message("TECHNICAL REPLICATE ANALYSIS - NB PATIENT SAMPLES")
message("Input: ", SEURAT_FILE); message("Output: ", OUTPUT_DIR); message("")

## LOAD DATA
message("[1/6] Loading Seurat object...")
seurat_obj <- readRDS(SEURAT_FILE)
message("Loaded: ", ncol(seurat_obj), " cells")

required_cols <- c("Patient", "Fraction", "Replicate")
missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
if (length(missing_cols) > 0) stop("ERROR: Missing columns: ", paste(missing_cols, collapse = ", "))

if ("Patient_ID" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$Patient_ID <- NULL
  message("Removed incorrect 'Patient_ID' column")
}

## IDENTIFY SAMPLES WITH REPLICATES
message("\n[2/6] Identifying samples with technical replicates...")
replicate_values <- unique(seurat_obj@meta.data$Replicate)
message("Unique Replicate values: ", paste(replicate_values, collapse = ", "))

sample_summary <- seurat_obj@meta.data %>%
  group_by(Patient, Fraction, Replicate) %>% summarise(Sample = first(Sample), N_cells = n(), .groups = "drop") %>%
  arrange(Patient, Fraction, Replicate)

message("\nSample summary:"); print(sample_summary)

replicate_pairs <- sample_summary %>%
  group_by(Patient, Fraction) %>% summarise(has_i = "i" %in% Replicate, has_ii = "ii" %in% Replicate,
                                           n_replicates = n(), .groups = "drop") %>%
  filter(has_i & has_ii) %>% mutate(pair_id = paste(Patient, Fraction, sep = "_"))

message("Found ", nrow(replicate_pairs), " patient-fraction pairs with both replicates")
if (nrow(replicate_pairs) == 0) stop("ERROR: No replicate pairs found.")

## CALCULATE REPLICATE CORRELATIONS
message("\n[3/6] Calculating pseudo-bulk replicate correlations...")
calculate_replicate_correlation <- function(seurat_obj, patient_id, fraction) {
  cells_rep_i <- which(seurat_obj$Patient == patient_id & seurat_obj$Fraction == fraction & seurat_obj$Replicate == "i")
  cells_rep_ii <- which(seurat_obj$Patient == patient_id & seurat_obj$Fraction == fraction & seurat_obj$Replicate == "ii")
  if (length(cells_rep_i) == 0 | length(cells_rep_ii) == 0) return(list(correlation = NA, n_cells_i = 0, n_cells_ii = 0))
  
  expr_data <- tryCatch({ LayerData(seurat_obj, assay = "SCT", layer = "data") },
                        error = function(e) { GetAssayData(seurat_obj, assay = "SCT", slot = "data") })
  pseudobulk_i <- rowMeans(expr_data[, cells_rep_i, drop = FALSE])
  pseudobulk_ii <- rowMeans(expr_data[, cells_rep_ii, drop = FALSE])
  cor_value <- cor(pseudobulk_i, pseudobulk_ii, method = "pearson")
  return(list(correlation = cor_value, n_cells_i = length(cells_rep_i), n_cells_ii = length(cells_rep_ii)))
}

replicate_cors <- data.frame()
for (i in 1:nrow(replicate_pairs)) {
  patient <- replicate_pairs$Patient[i]
  fraction <- replicate_pairs$Fraction[i]
  result <- calculate_replicate_correlation(seurat_obj, patient, fraction)
  replicate_cors <- rbind(replicate_cors, data.frame(Patient = patient, Fraction = fraction,
                                                     Correlation = result$correlation,
                                                     N_cells_rep_i = result$n_cells_i,
                                                     N_cells_rep_ii = result$n_cells_ii,
                                                     Total_cells = result$n_cells_i + result$n_cells_ii))
}

replicate_cors$Status <- ifelse(is.na(replicate_cors$Correlation), "Missing/Failed",
                                ifelse(replicate_cors$Correlation > 0.95, "Excellent",
                                       ifelse(replicate_cors$Correlation > 0.90, "Good", "Adequate")))

message("\nReplicate Correlation Summary:"); print(replicate_cors)
write.csv(replicate_cors, file.path(OUTPUT_DIR, "tables", "replicate_correlations.csv"), row.names = FALSE)

## VISUALIZE CORRELATIONS
message("\n[4/6] Creating visualizations...")
p1 <- ggplot(replicate_cors[!is.na(replicate_cors$Correlation),],
             aes(x = paste(Patient, Fraction, sep = "_"), y = Correlation, fill = Status)) +
  geom_bar(stat = "identity") + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "orange", linewidth = 0.8) +
  geom_text(aes(label = round(Correlation, 3)), vjust = -0.5, size = 3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), plot.title = element_text(face = "bold", size = 14)) +
  ylim(0, 1.05) + labs(title = "Technical Replicate Correlation (Pseudo-bulk)",
                       subtitle = "Red line: >0.95, Orange line: >0.90", x = "Patient_Fraction", y = "Pearson Correlation") +
  scale_fill_manual(values = c("Excellent" = "#2E7D32", "Good" = "#FFA726", "Adequate" = "#FFD54F"))

ggsave(file.path(OUTPUT_DIR, "plots", "01_replicate_correlation_all.pdf"), p1, width = 12, height = 7)
ggsave(file.path(OUTPUT_DIR, "plots", "01_replicate_correlation_all.png"), p1, width = 12, height = 7, dpi = 300)
message("Saved: 01_replicate_correlation_all.pdf/png")

## TECHNICAL VS BIOLOGICAL VARIATION
message("\n[5/6] Comparing technical vs biological variation...")
longitudinal_patients <- replicate_pairs %>%
  group_by(Patient) %>% summarise(has_F1 = "F1" %in% Fraction, has_F2 = "F2" %in% Fraction, .groups = "drop") %>%
  filter(has_F1 & has_F2) %>% pull(Patient)

message("Found ", length(longitudinal_patients), " patients with longitudinal data: ", paste(longitudinal_patients, collapse = ", "))

if (length(longitudinal_patients) > 0) {
  compare_technical_biological <- function(seurat_obj, patient_id) {
    expr_data <- tryCatch({ LayerData(seurat_obj, assay = "SCT", layer = "data") },
                          error = function(e) { GetAssayData(seurat_obj, assay = "SCT", slot = "data") })
    f1_i <- which(seurat_obj$Patient == patient_id & seurat_obj$Fraction == "F1" & seurat_obj$Replicate == "i")
    f1_ii <- which(seurat_obj$Patient == patient_id & seurat_obj$Fraction == "F1" & seurat_obj$Replicate == "ii")
    f2_i <- which(seurat_obj$Patient == patient_id & seurat_obj$Fraction == "F2" & seurat_obj$Replicate == "i")
    f2_ii <- which(seurat_obj$Patient == patient_id & seurat_obj$Fraction == "F2" & seurat_obj$Replicate == "ii")
    
    f1_rep_i <- rowMeans(expr_data[, f1_i, drop = FALSE])
    f1_rep_ii <- rowMeans(expr_data[, f1_ii, drop = FALSE])
    f2_rep_i <- rowMeans(expr_data[, f2_i, drop = FALSE])
    f2_rep_ii <- rowMeans(expr_data[, f2_ii, drop = FALSE])
    f1_avg <- (f1_rep_i + f1_rep_ii) / 2
    f2_avg <- (f2_rep_i + f2_rep_ii) / 2
    
    cor_f1_tech <- cor(f1_rep_i, f1_rep_ii, method = "pearson")
    cor_f2_tech <- cor(f2_rep_i, f2_rep_ii, method = "pearson")
    cor_f1_vs_f2 <- cor(f1_avg, f2_avg, method = "pearson")
    
    return(data.frame(Patient = patient_id, Comparison = c("F1: Rep-i vs Rep-ii (Technical)", "F2: Rep-i vs Rep-ii (Technical)", "F1 vs F2 (Biological)"),
                     Correlation = c(cor_f1_tech, cor_f2_tech, cor_f1_vs_f2), Type = c("Technical", "Technical", "Biological")))
  }
  
  all_comparisons <- data.frame()
  for (patient in longitudinal_patients) {
    comp <- compare_technical_biological(seurat_obj, patient)
    all_comparisons <- rbind(all_comparisons, comp)
  }
  
  message("\nTechnical vs Biological Variation:"); print(all_comparisons)
  write.csv(all_comparisons, file.path(OUTPUT_DIR, "tables", "technical_vs_biological_variation.csv"), row.names = FALSE)
  
  p3 <- ggplot(all_comparisons, aes(x = Comparison, y = Correlation, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    geom_text(aes(label = round(Correlation, 3)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
    ylim(0, 1.05) + labs(title = "Technical vs Biological Variation",
                         subtitle = "Technical correlation should be HIGH, F1 vs F2 may be LOWER", x = "", y = "Pearson Correlation") +
    facet_wrap(~Patient, ncol = length(longitudinal_patients)) + scale_fill_manual(values = c("Technical" = "#4CAF50", "Biological" = "#2196F3"))
  
  ggsave(file.path(OUTPUT_DIR, "plots", "03_technical_vs_biological_variation.pdf"), p3, width = 10, height = 6)
  message("Saved: 03_technical_vs_biological_variation.pdf")
}

## SUMMARY
message("\n[6/6] Summary...")
message("\n============================================================")
message("REPLICATE ANALYSIS COMPLETE")
message("============================================================")
message("Total replicate pairs analyzed: ", nrow(replicate_cors))
message("Mean correlation (all pairs): ", round(mean(replicate_cors$Correlation, na.rm = TRUE), 3))
message("Longitudinal patients: ", length(longitudinal_patients))

excellent <- sum(replicate_cors$Status == "Excellent", na.rm = TRUE)
good <- sum(replicate_cors$Status == "Good", na.rm = TRUE)
adequate <- sum(replicate_cors$Status == "Adequate", na.rm = TRUE)

message("\nQuality Assessment:")
message("  Excellent (>0.95): ", excellent, " pairs")
message("  Good (>0.90): ", good, " pairs")
message("  Adequate (<0.90): ", adequate, " pairs")

message("\nOutput files: ", OUTPUT_DIR)
message("============================================================\n")
