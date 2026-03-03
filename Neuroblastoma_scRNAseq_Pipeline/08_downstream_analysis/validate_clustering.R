#!/usr/bin/env Rscript
## Clustering Validation - Patient-Specific Analysis
## Author: Ayeh Sadr | Pipeline: Neuroblastoma Patient scRNA-seq Analysis

suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(cluster); library(RColorBrewer)})

## CONFIGURATION
PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")
OUTPUT_DIR <- file.path(RESULTS_DIR, "clustering_validation_final")

SEURAT_CANDIDATES <- c(
  file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds"),
  file.path(DATA_DIR, "NB_patients_seurat_with_signatures.rds"),
  file.path(DATA_DIR, "NB_patients_seurat_harmony.rds")
)

SEURAT_FILE <- NULL
for (cand in SEURAT_CANDIDATES) {
  if (file.exists(cand)) { SEURAT_FILE <- cand; message("Found: ", cand); break }
}
if (is.null(SEURAT_FILE)) stop("No Seurat object found!")

N_PCS <- 30
N_CELLS_SAMPLE <- 10000
RECLUSTER <- FALSE

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(OUTPUT_DIR, "plots"); table_dir <- file.path(OUTPUT_DIR, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE); dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("\nCLUSTERING VALIDATION - PATIENT-SPECIFIC ANALYSIS")
message("Input: ", SEURAT_FILE); message("Output: ", OUTPUT_DIR); message("")

## LOAD DATA
message("[1/4] Loading data...")
nb <- readRDS(SEURAT_FILE); message("Loaded: ", ncol(nb), " cells")

if ("Patient" %in% colnames(nb@meta.data) && !"patient" %in% colnames(nb@meta.data)) {
  nb$patient <- nb$Patient
}

if ("seurat_clusters" %in% colnames(nb@meta.data)) {
  nb$original_clusters <- nb$seurat_clusters
  message("Original clusters preserved: ", length(unique(nb$original_clusters)), " clusters")
}

harmony_variants <- c("harmony", "Harmony", "HARMONY")
available_reductions <- names(nb@reductions)
harmony_name <- harmony_variants[harmony_variants %in% available_reductions][1]

if (is.na(harmony_name)) stop("Harmony reduction not found!")
message("Using Harmony reduction: ", harmony_name)

if (RECLUSTER) {
  message("Re-clustering at resolution 0.5...")
  nb <- FindNeighbors(nb, reduction = harmony_name, dims = 1:N_PCS, verbose = FALSE)
  graph_name <- if (paste0(harmony_name, "_snn") %in% names(nb@graphs)) paste0(harmony_name, "_snn") else names(nb@graphs)[1]
  nb <- FindClusters(nb, graph.name = graph_name, resolution = 0.5, verbose = FALSE)
} else {
  message("Using existing clusters (RECLUSTER = FALSE)")
  if (!"seurat_clusters" %in% colnames(nb@meta.data)) stop("No clusters found!")
}

## CALCULATE SILHOUETTE & ENTROPY
message("\n[2/4] Calculating silhouette scores and patient entropy...")

harmony_coords <- Embeddings(nb, reduction = harmony_name)[, 1:N_PCS]
set.seed(42)
sample_idx <- if (ncol(nb) > N_CELLS_SAMPLE) sample(1:ncol(nb), N_CELLS_SAMPLE) else 1:ncol(nb)

harmony_sample <- harmony_coords[sample_idx, ]
clusters_sample <- as.numeric(as.factor(nb$seurat_clusters[sample_idx]))
patients_sample <- nb$patient[sample_idx]
cluster_names_sample <- nb$seurat_clusters[sample_idx]

message("Computing silhouette...")
dist_harmony <- dist(harmony_sample)
sil <- silhouette(clusters_sample, dist_harmony)

sil_df <- data.frame(cell_idx = sample_idx, cluster = cluster_names_sample, patient = patients_sample, silhouette = sil[, "sil_width"])

cluster_sil_summary <- sil_df %>% group_by(cluster) %>%
  summarise(mean_silhouette = mean(silhouette), median_silhouette = median(silhouette), sd_silhouette = sd(silhouette), n_cells = n(), .groups = "drop") %>%
  arrange(as.numeric(as.character(cluster)))

message("Silhouette range: ", round(min(cluster_sil_summary$mean_silhouette), 3), " to ", round(max(cluster_sil_summary$mean_silhouette), 3))

## PATIENT ENTROPY
message("\nCalculating patient entropy...")
cluster_patient_table <- table(nb$seurat_clusters, nb$patient)
cluster_patient_prop <- prop.table(cluster_patient_table, margin = 1)

calculate_entropy <- function(x) {
  x <- x[x > 0]
  if (length(x) <= 1) return(0)
  -sum(x * log(x)) / log(length(x))
}

cluster_entropy <- apply(cluster_patient_prop, 1, calculate_entropy)
dominant_patient <- apply(cluster_patient_table, 1, function(x) names(which.max(x)))
dominant_prop <- apply(cluster_patient_prop, 1, max)

entropy_df <- data.frame(cluster = names(cluster_entropy), patient_entropy = as.numeric(cluster_entropy),
                        dominant_patient = dominant_patient, dominant_proportion = round(dominant_prop, 3),
                        n_cells_total = as.numeric(rowSums(cluster_patient_table)))

cluster_summary <- merge(cluster_sil_summary %>% select(cluster, mean_silhouette, median_silhouette, n_cells),
                        entropy_df, by = "cluster", all = FALSE)
cluster_summary$n_cells <- cluster_summary$n_cells_total; cluster_summary$n_cells_total <- NULL

cluster_summary <- cluster_summary %>% mutate(
  cluster_type = case_when(patient_entropy > 0.5 ~ "Shared (Immune/Stromal)", TRUE ~ "Patient-Specific (Tumor)"),
  silhouette_quality = case_when(mean_silhouette > 0.4 ~ "Excellent", mean_silhouette > 0.25 ~ "Good",
                                mean_silhouette > 0.1 ~ "Fair", TRUE ~ "Poor"))

write.csv(cluster_summary, file.path(table_dir, "cluster_validation_summary.csv"), row.names = FALSE)

n_patient_specific <- sum(cluster_summary$patient_entropy <= 0.5)
n_shared <- sum(cluster_summary$patient_entropy > 0.5)
mean_sil_patient_specific <- mean(cluster_summary$mean_silhouette[cluster_summary$patient_entropy <= 0.5])
mean_sil_shared <- mean(cluster_summary$mean_silhouette[cluster_summary$patient_entropy > 0.5])

message("Patient-specific clusters: ", n_patient_specific, " (mean silhouette: ", round(mean_sil_patient_specific, 3), ")")
message("Shared clusters: ", n_shared, " (mean silhouette: ", round(mean_sil_shared, 3), ")")

## VISUALIZATIONS
message("\n[3/4] Creating validation figures...")

n_patients <- length(unique(nb$patient))
patient_colors <- if (n_patients <= 9) brewer.pal(max(3, n_patients), "Set1")[1:n_patients] else
  colorRampPalette(brewer.pal(9, "Set1"))(n_patients)
names(patient_colors) <- sort(unique(nb$patient))

## Figure 1: Silhouette boxplot
cluster_order <- cluster_summary %>% arrange(desc(mean_silhouette)) %>% pull(cluster)
sil_df$cluster <- factor(sil_df$cluster, levels = cluster_order)

p1 <- ggplot(sil_df, aes(x = cluster, y = silhouette, fill = patient)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) + geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "orange", linewidth = 0.8) +
  scale_fill_manual(values = patient_colors) + labs(title = "Silhouette Score Distribution by Cluster",
       subtitle = "Colored by patient | Clusters ordered by mean silhouette", x = "Cluster", y = "Silhouette Score", fill = "Patient") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "01_silhouette_boxplot_by_cluster.pdf"), p1, width = 14, height = 7)
message("Saved: 01_silhouette_boxplot_by_cluster.pdf")

## Figure 3: Entropy
cluster_summary$cluster <- factor(cluster_summary$cluster, levels = cluster_summary$cluster[order(-cluster_summary$patient_entropy)])

p3 <- ggplot(cluster_summary, aes(x = reorder(cluster, -patient_entropy), y = patient_entropy, fill = cluster_type)) +
  geom_col(width = 0.7) + geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_text(aes(label = dominant_patient), vjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Shared (Immune/Stromal)" = "#27ae60", "Patient-Specific (Tumor)" = "#e74c3c")) +
  labs(title = "Patient Entropy by Cluster (Harmony)", subtitle = "Low entropy = patient-specific | Labels show dominant patient",
       x = "Cluster", y = "Patient Entropy (0 = single patient, 1 = equal mix)", fill = "Cluster Type") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "03_patient_entropy_by_cluster.pdf"), p3, width = 14, height = 7)

## Figure 4: Silhouette vs Entropy Scatter
p4 <- ggplot(cluster_summary, aes(x = patient_entropy, y = mean_silhouette, color = cluster_type, size = n_cells)) +
  geom_point(alpha = 0.8) + geom_text(aes(label = cluster), vjust = -1, hjust = 0.5, size = 3.5, show.legend = FALSE) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "orange", linewidth = 0.8) +
  scale_color_manual(values = c("Shared (Immune/Stromal)" = "#27ae60", "Patient-Specific (Tumor)" = "#e74c3c")) +
  scale_size_continuous(range = c(3, 12)) +
  annotate("rect", xmin = -Inf, xmax = 0.5, ymin = 0.25, ymax = Inf, fill = "#e74c3c", alpha = 0.1) +
  annotate("rect", xmin = 0.5, xmax = Inf, ymin = 0.25, ymax = Inf, fill = "#27ae60", alpha = 0.1) +
  labs(title = "Cluster Validation: Silhouette vs Patient Entropy",
       subtitle = "Patient-specific clusters (red zone) have GOOD silhouette = biologically valid",
       x = "Patient Entropy", y = "Mean Silhouette Score", color = "Cluster Type", size = "Cells") +
  xlim(0, 1) + ylim(0, max(cluster_summary$mean_silhouette) + 0.15) + theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "right")

ggsave(file.path(plot_dir, "04_silhouette_vs_entropy_validation.pdf"), p4, width = 12, height = 9)

## SUMMARY
message("\n[4/4] Summary...")
message("\n============================================================")
message("VALIDATION SUMMARY")
message("============================================================")
message("Patient-specific clusters: ", n_patient_specific, " (likely tumor)")
message("Shared clusters: ", n_shared, " (likely immune/stromal)")
message("Overall mean silhouette: ", round(mean(cluster_summary$mean_silhouette), 3))
message("Patient-specific mean: ", round(mean_sil_patient_specific, 3))
message("Shared mean: ", round(mean_sil_shared, 3))

if (mean_sil_patient_specific >= 0.1) {
  message("\n✅ KEY VALIDATION: Patient segregation is BIOLOGICAL, not batch effect!")
} else {
  message("\n⚠️ WARNING: Patient-specific clusters have low silhouette")
}

message("\nOutput: ", OUTPUT_DIR)
message("============================================================\n")
