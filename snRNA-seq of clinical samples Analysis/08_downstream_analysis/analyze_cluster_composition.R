#!/usr/bin/env Rscript
## Cluster Composition Analysis - See original for full documentation
## Author: Ayeh Sadr | Pipeline: Neuroblastoma Patient scRNA-seq Analysis

library(Seurat); library(ggplot2); library(dplyr); library(tidyr); 
library(patchwork); library(ggalluvial); library(scales)
set.seed(12345)

## CONFIGURATION
PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")
OUTPUT_DIR <- file.path(RESULTS_DIR, "composition_analysis")
SEURAT_FILE <- file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds")

MIN_CHANGE_PCT <- 5.0
MIN_CHANGE_ABSOLUTE <- 100
FISHER_PADJ_THRESHOLD <- 0.05

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

message("CLUSTER COMPOSITION ANALYSIS - LONGITUDINAL NB PATIENTS")
message("Input: ", SEURAT_FILE); message("Output: ", OUTPUT_DIR); message("")

## LOAD & IDENTIFY LONGITUDINAL PATIENTS
message("[1/9] Loading data and identifying longitudinal patients...")
seurat_obj <- readRDS(SEURAT_FILE)
if ("SCT" %in% names(seurat_obj@assays)) {
  DefaultAssay(seurat_obj) <- "SCT"
  message("Using SCT assay")
} else {
  DefaultAssay(seurat_obj) <- "RNA"
  message("Using RNA assay")
}
message("Loaded: ", ncol(seurat_obj), " cells")

longitudinal_summary <- seurat_obj@meta.data %>%
  group_by(Patient, Fraction) %>% summarise(N_cells = n(), .groups = "drop") %>%
  group_by(Patient) %>%
  summarise(has_F1 = "F1" %in% Fraction, has_F2 = "F2" %in% Fraction, n_timepoints = n(), .groups = "drop") %>%
  filter(has_F1 & has_F2)

longitudinal_patients <- longitudinal_summary$Patient
seurat_long <- subset(seurat_obj, subset = Patient %in% longitudinal_patients)
message("Analyzing ", ncol(seurat_long), " cells from ", length(longitudinal_patients), " patients")

## EXTRACT CLUSTER COLORS
message("\n[2/9] Extracting cluster colors...")
all_clusters <- sort(unique(as.numeric(as.character(seurat_long@meta.data$seurat_clusters))))
p_temp <- DimPlot(seurat_long, group.by = "seurat_clusters", reduction = "umap", label = FALSE)
gg_build <- ggplot_build(p_temp); plot_data <- gg_build$data[[1]]
if (nrow(plot_data) > 0 && "colour" %in% colnames(plot_data) && "group" %in% colnames(plot_data)) {
  color_map <- plot_data[, c("group", "colour")][!duplicated(plot_data$group), ]
  cluster_factor <- factor(seurat_long@meta.data$seurat_clusters)
  cluster_levels <- levels(cluster_factor)
  CLUSTER_COLORS <- setNames(color_map$colour, as.character(cluster_levels[color_map$group]))
  missing <- setdiff(as.character(all_clusters), names(CLUSTER_COLORS))
  if (length(missing) > 0) {
    default_colors <- Seurat::DiscretePalette(length(all_clusters), palette = "alphabet2")
    names(default_colors) <- as.character(all_clusters)
    CLUSTER_COLORS <- c(CLUSTER_COLORS, default_colors[missing])
  }
  CLUSTER_COLORS <- CLUSTER_COLORS[as.character(all_clusters)]
} else {
  CLUSTER_COLORS <- Seurat::DiscretePalette(length(all_clusters), palette = "alphabet2")
  names(CLUSTER_COLORS) <- as.character(all_clusters)
}
message("Extracted colors for ", length(CLUSTER_COLORS), " clusters")

## CALCULATE COMPOSITION
message("\n[3/9] Calculating cluster compositions...")
composition <- seurat_long@meta.data %>%
  group_by(Patient, Fraction, seurat_clusters) %>% summarise(Count = n(), .groups = "drop") %>%
  group_by(Patient, Fraction) %>%
  mutate(Total = sum(Count), Proportion = Count / Total, Percentage = Proportion * 100) %>% ungroup()
message("Composition calculated for ", nrow(composition), " patient-timepoint-cluster combinations")
write.csv(composition, file.path(OUTPUT_DIR, "tables", "00_cluster_composition_raw.csv"), row.names = FALSE)

## VISUALIZATIONS: STACKED BARPLOTS
message("\n[4/9] Creating stacked barplots...")
p_stacked <- ggplot(composition, aes(x = paste(Patient, Fraction, sep = "_"), y = Percentage, fill = factor(seurat_clusters))) +
  geom_bar(stat = "identity", color = "white", size = 0.3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), legend.position = "right") +
  labs(title = "Cluster Composition: Longitudinal Analysis (F1 → F2)", x = "Patient_Timepoint", y = "Percentage of Cells (%)", fill = "Cluster") +
  scale_fill_manual(values = CLUSTER_COLORS) +
  geom_text(data = composition %>% filter(Percentage > 3), aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 2.5, color = "black")
ggsave(file.path(OUTPUT_DIR, "plots", "01_composition_stacked_all.pdf"), p_stacked, width = 12, height = 7)

p_sidebyside <- ggplot(composition, aes(x = Fraction, y = Percentage, fill = factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
  facet_wrap(~Patient, ncol = length(longitudinal_patients)) + theme_classic() +
  theme(strip.background = element_rect(fill = "gray90"), strip.text = element_text(size = 12, face = "bold")) +
  labs(title = "Cluster Composition Changes: F1 → F2", x = "Timepoint", y = "Percentage of Cells (%)", fill = "Cluster") +
  scale_fill_manual(values = CLUSTER_COLORS)
ggsave(file.path(OUTPUT_DIR, "plots", "02_composition_sidebyside.pdf"), p_sidebyside, width = 4 * length(longitudinal_patients), height = 6)
message("Saved stacked barplots")

## CALCULATE & VISUALIZE CHANGES
message("\n[5/9] Calculating cluster changes (F2 - F1)...")
composition_wide <- composition %>%
  select(Patient, Fraction, seurat_clusters, Percentage, Count) %>%
  pivot_wider(names_from = Fraction, values_from = c(Percentage, Count), values_fill = list(Percentage = 0, Count = 0)) %>%
  mutate(Change_Percentage = Percentage_F2 - Percentage_F1, Change_Count = Count_F2 - Count_F1,
         Fold_Change = ifelse(Percentage_F1 == 0, NA, Percentage_F2 / Percentage_F1),
         Status = case_when(abs(Change_Percentage) > 10 ~ "Major Change (>10%)", abs(Change_Percentage) > 5 ~ "Moderate Change (>5%)",
                           abs(Change_Percentage) > 2 ~ "Minor Change (>2%)", TRUE ~ "Stable"))
write.csv(composition_wide, file.path(OUTPUT_DIR, "tables", "01_cluster_changes_F1_to_F2.csv"), row.names = FALSE)

p_changes <- ggplot(composition_wide, aes(x = factor(seurat_clusters), y = Change_Percentage, fill = Status)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) + geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  geom_hline(yintercept = c(-10, -5, 5, 10), linetype = "dashed", color = "gray50", alpha = 0.5) +
  facet_wrap(~Patient, ncol = 1) + theme_classic() +
  theme(strip.background = element_rect(fill = "gray90"), strip.text = element_text(size = 12, face = "bold"), legend.position = "bottom") +
  labs(title = "Cluster Proportion Changes: F2 - F1", subtitle = "Positive = Expansion in F2, Negative = Contraction in F2",
       x = "Cluster", y = "Change in Percentage (%)", fill = "Change Category") +
  scale_fill_manual(values = c("Major Change (>10%)" = "#E53935", "Moderate Change (>5%)" = "#FFA726",
                               "Minor Change (>2%)" = "#FFD54F", "Stable" = "#9E9E9E"))
ggsave(file.path(OUTPUT_DIR, "plots", "03_cluster_changes_barplot.pdf"), p_changes, width = 12, height = 4 * length(longitudinal_patients))

## STATISTICAL TESTING (FISHER'S EXACT TEST)
message("\n[6/9] Running statistical tests...")
test_cluster_changes <- function(composition_df, patient_id) {
  patient_data <- composition_df %>% filter(Patient == patient_id)
  totals <- patient_data %>% group_by(Fraction) %>% summarize(Total = sum(Count), .groups = "drop")
  f1_total <- totals$Total[totals$Fraction == "F1"]
  f2_total <- totals$Total[totals$Fraction == "F2"]
  test_results <- data.frame()
  for (clust in unique(patient_data$seurat_clusters)) {
    f1_data <- patient_data %>% filter(Fraction == "F1", seurat_clusters == clust)
    f2_data <- patient_data %>% filter(Fraction == "F2", seurat_clusters == clust)
    f1_count <- ifelse(nrow(f1_data) > 0, f1_data$Count, 0)
    f2_count <- ifelse(nrow(f2_data) > 0, f2_data$Count, 0)
    cont_table <- matrix(c(f1_count, f1_total - f1_count, f2_count, f2_total - f2_count), nrow = 2, byrow = TRUE)
    fisher_result <- fisher.test(cont_table)
    test_results <- rbind(test_results, data.frame(Patient = patient_id, Cluster = clust, F1_Count = f1_count, F2_Count = f2_count,
                                                   F1_Prop = f1_count / f1_total, F2_Prop = f2_count / f2_total,
                                                   Change_Prop = (f2_count / f2_total) - (f1_count / f1_total),
                                                   Fold_Change = ifelse(f1_count == 0, NA, (f2_count / f2_total) / (f1_count / f1_total)),
                                                   P_Value = fisher_result$p.value, Odds_Ratio = fisher_result$estimate))
  }
  test_results$P_Adjusted <- p.adjust(test_results$P_Value, method = "BH")
  test_results$Significant <- test_results$P_Adjusted < FISHER_PADJ_THRESHOLD
  return(test_results)
}

all_test_results <- data.frame()
for (patient in longitudinal_patients) {
  test_result <- test_cluster_changes(composition, patient)
  all_test_results <- rbind(all_test_results, test_result)
}
write.csv(all_test_results, file.path(OUTPUT_DIR, "tables", "02_statistical_tests_fisher.csv"), row.names = FALSE)

sig_changes <- all_test_results %>% filter(Significant) %>% arrange(Patient, P_Adjusted)
if (nrow(sig_changes) > 0) {
  message("  Significant changes found: ", nrow(sig_changes))
} else {
  message("  No significant changes found")
}

## VOLCANO PLOTS
message("\n[7/9] Creating volcano plots...")
volcano_data <- all_test_results %>%
  mutate(Change_Percentage = Change_Prop * 100, NegLog10P = ifelse(P_Adjusted > 0, -log10(P_Adjusted), -log10(.Machine$double.xmin)),
         Label = ifelse(Significant & abs(Change_Percentage) > MIN_CHANGE_PCT, as.character(Cluster), NA),
         Direction = case_when(Significant & Change_Percentage > MIN_CHANGE_PCT ~ "Expanded",
                              Significant & Change_Percentage < -MIN_CHANGE_PCT ~ "Contracted", TRUE ~ "Not Significant"))

p_volcano <- ggplot(volcano_data, aes(x = Change_Percentage, y = NegLog10P, color = Direction, size = abs(Change_Percentage))) +
  geom_point(alpha = 0.7) + geom_text(aes(label = Label), vjust = -0.5, size = 3, color = "black", na.rm = TRUE, check_overlap = TRUE) +
  geom_vline(xintercept = c(-MIN_CHANGE_PCT, MIN_CHANGE_PCT), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(FISHER_PADJ_THRESHOLD), linetype = "dashed", color = "gray50") +
  facet_wrap(~Patient, ncol = length(longitudinal_patients)) + theme_classic() +
  theme(strip.background = element_rect(fill = "gray90"), strip.text = element_text(size = 12, face = "bold")) +
  labs(title = "Volcano Plot: Cluster Composition Changes", x = "Change in Proportion (F2 - F1, %)", y = "-log10(Adjusted P-value)",
       color = "Status", size = "|Change (%)|") +
  scale_color_manual(values = c("Expanded" = "#4CAF50", "Contracted" = "#F44336", "Not Significant" = "gray70")) +
  scale_size_continuous(range = c(1, 6))
ggsave(file.path(OUTPUT_DIR, "plots", "04_volcano_composition_changes.pdf"), p_volcano, width = 5 * length(longitudinal_patients), height = 6)

## ALLUVIAL DIAGRAMS
message("\n[8/9] Creating alluvial flow diagrams...")
for (patient in longitudinal_patients) {
  patient_comp <- composition %>% filter(Patient == patient)
  flow_data <- patient_comp %>% select(Fraction, seurat_clusters, Count) %>%
    mutate(Fraction = factor(Fraction, levels = c("F1", "F2")))
  
  p_alluvial <- ggplot(flow_data, aes(x = Fraction, stratum = seurat_clusters, alluvium = seurat_clusters, y = Count,
                                      fill = factor(seurat_clusters), label = seurat_clusters)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray", alpha = 0.6, width = 0.3) +
    geom_stratum(alpha = 0.8, color = "white", width = 0.3) + geom_text(stat = "stratum", size = 3) + theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14)) +
    labs(title = paste0("Patient ", patient, ": Cluster Flow (F1 → F2)"), subtitle = "Width = number of cells", x = "Timepoint", y = "Number of Cells") +
    scale_fill_manual(values = CLUSTER_COLORS)
  
  ggsave(file.path(OUTPUT_DIR, "plots", paste0("05_alluvial_patient_", patient, ".pdf")), p_alluvial, width = 8, height = 10)
}

## UMAP VISUALIZATIONS
message("\n[9/9] Creating UMAP visualizations...")
for (patient in longitudinal_patients) {
  patient_obj <- subset(seurat_long, subset = Patient == patient)
  if (ncol(patient_obj) == 0) next
  
  p_umap_timepoint <- DimPlot(patient_obj, group.by = "Fraction", reduction = "umap", pt.size = 0.5,
                              cols = c("F1" = "blue", "F2" = "red")) + ggtitle(paste0("Patient ", patient, ": F1 vs F2")) + theme_minimal()
  ggsave(file.path(OUTPUT_DIR, "plots", paste0("06_umap_timepoint_patient_", patient, ".pdf")), p_umap_timepoint, width = 8, height = 7)
  
  p_umap_clusters <- DimPlot(patient_obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE, pt.size = 0.5) +
    ggtitle(paste0("Patient ", patient, ": Clusters")) + theme_minimal()
  ggsave(file.path(OUTPUT_DIR, "plots", paste0("07_umap_clusters_patient_", patient, ".pdf")), p_umap_clusters, width = 8, height = 7)
}

message("\n============================================================")
message("CLUSTER COMPOSITION ANALYSIS COMPLETE")
message("Output: ", OUTPUT_DIR)
message("============================================================\n")
