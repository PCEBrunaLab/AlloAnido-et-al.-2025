#!/usr/bin/env Rscript
## Correlation Visualization: Technical Replicates vs Biological Timepoints
## Author: Ayeh Sadr | Pipeline: Neuroblastoma Patient scRNA-seq Analysis

library(Seurat); library(ggplot2); library(dplyr); library(patchwork); library(viridis)
set.seed(12345)

## CONFIGURATION
PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR <- file.path(RESULTS_DIR, "data")
OUTPUT_DIR <- file.path(RESULTS_DIR, "correlation_analysis")
SEURAT_FILE <- file.path(DATA_DIR, "NB_patients_seurat_with_jansky.rds")

PATIENTS <- c("649", "1087")
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

message("Loading Seurat object...")
seurat_obj <- readRDS(SEURAT_FILE)

## HELPER FUNCTIONS
get_pseudobulk <- function(seurat_obj, patient, fraction, replicate = NULL) {
  cells <- seurat_obj$Patient == patient & seurat_obj$Fraction == fraction
  if (!is.null(replicate)) cells <- cells & seurat_obj$Replicate == replicate
  expr <- tryCatch(LayerData(seurat_obj, assay = "SCT", layer = "data"),
                  error = function(e) GetAssayData(seurat_obj, assay = "SCT", slot = "data"))
  rowMeans(expr[, cells, drop = FALSE])
}

density_scatter <- function(x, y, title, xlab, ylab) {
  df <- data.frame(x = x, y = y)
  r <- cor(x, y, method = "pearson")
  ggplot(df, aes(x, y)) + geom_hex(bins = 100) + scale_fill_viridis(option = "plasma", trans = "log10", name = "Count") +
    geom_abline(slope = 1, intercept = 0, color = "white", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, label = sprintf("r = %.4f", r), hjust = 1.1, vjust = -1, size = 5, fontface = "bold", color = "white") +
    labs(title = title, x = xlab, y = ylab) + theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 11), panel.grid.minor = element_blank(), legend.position = "right") +
    coord_fixed()
}

message("Calculating pseudo-bulk profiles...")
pseudobulk <- list()
for (pt in PATIENTS) {
  pseudobulk[[pt]] <- list(
    F1_i = get_pseudobulk(seurat_obj, pt, "F1", "i"),
    F1_ii = get_pseudobulk(seurat_obj, pt, "F1", "ii"),
    F2_i = get_pseudobulk(seurat_obj, pt, "F2", "i"),
    F2_ii = get_pseudobulk(seurat_obj, pt, "F2", "ii"),
    F1 = get_pseudobulk(seurat_obj, pt, "F1"),
    F2 = get_pseudobulk(seurat_obj, pt, "F2")
  )
}

## FIGURE 1: DENSITY SCATTER PLOTS
message("Generating density scatter plots...")
plots <- list()
for (pt in PATIENTS) {
  pb <- pseudobulk[[pt]]
  plots[[paste0(pt, "_F1_tech")]] <- density_scatter(pb$F1_i, pb$F1_ii,
    paste0("Patient ", pt, " - F1 Replicates (Technical)"), "F1 Rep-i", "F1 Rep-ii")
  plots[[paste0(pt, "_F2_tech")]] <- density_scatter(pb$F2_i, pb$F2_ii,
    paste0("Patient ", pt, " - F2 Replicates (Technical)"), "F2 Rep-i", "F2 Rep-ii")
  plots[[paste0(pt, "_bio")]] <- density_scatter(pb$F1, pb$F2,
    paste0("Patient ", pt, " - F1 vs F2 (Biological)"), "F1 (pre-therapy)", "F2 (post-therapy)")
}

combined <- (plots$`649_F1_tech` | plots$`649_F2_tech` | plots$`649_bio`) /
            (plots$`1087_F1_tech` | plots$`1087_F2_tech` | plots$`1087_bio`) +
  plot_annotation(title = "Technical vs Biological Correlation",
                 subtitle = "Technical replicates (high r expected) vs Timepoints (lower r = biological change)",
                 theme = theme(plot.title = element_text(face = "bold", size = 16)))

ggsave(file.path(OUTPUT_DIR, "plots", "01_density_scatter_all.pdf"), combined, width = 15, height = 10, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "plots", "01_density_scatter_all.png"), combined, width = 15, height = 10, dpi = 300)
message("Saved: 01_density_scatter_all.pdf/png")

## FIGURE 2: CORRELATION HEATMAP
message("Generating correlation heatmap...")
plot_cor_heatmap <- function(pb, patient) {
  samples <- c("F1_i", "F1_ii", "F2_i", "F2_ii")
  mat <- sapply(samples, function(s) pb[[s]])
  cor_mat <- cor(mat, method = "pearson")
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Sample1", "Sample2", "Correlation")
  cor_df$Type <- ifelse((grepl("F1", cor_df$Sample1) & grepl("F1", cor_df$Sample2)) |
                         (grepl("F2", cor_df$Sample1) & grepl("F2", cor_df$Sample2)), "Technical", "Biological")
  cor_df$Type[cor_df$Sample1 == cor_df$Sample2] <- "Self"
  
  ggplot(cor_df, aes(Sample1, Sample2, fill = Correlation)) + geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Correlation)), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0.9, limits = c(0.7, 1)) +
    labs(title = paste("Patient", patient), x = "", y = "") + theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) + coord_fixed()
}

heatmap_649 <- plot_cor_heatmap(pseudobulk[["649"]], "649 (Responder)")
heatmap_1087 <- plot_cor_heatmap(pseudobulk[["1087"]], "1087 (Resistant)")

heatmap_combined <- heatmap_649 | heatmap_1087 +
  plot_annotation(title = "Pairwise Correlation Matrix",
                 subtitle = "Diagonal blocks = Technical replicates | Off-diagonal = Biological (F1 vs F2)",
                 theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUTPUT_DIR, "plots", "02_correlation_heatmap.pdf"), heatmap_combined, width = 12, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "plots", "02_correlation_heatmap.png"), heatmap_combined, width = 12, height = 5, dpi = 300)
message("Saved: 02_correlation_heatmap.pdf/png")

## FIGURE 3: BAR PLOT SUMMARY
message("Generating summary bar plot...")
calc_cors <- function(pb, patient) {
  data.frame(Patient = patient, Comparison = c("F1: i vs ii", "F2: i vs ii", "F1 vs F2"), Type = c("Technical", "Technical", "Biological"),
            Correlation = c(cor(pb$F1_i, pb$F1_ii), cor(pb$F2_i, pb$F2_ii), cor(pb$F1, pb$F2)))
}

cor_summary <- rbind(calc_cors(pseudobulk[["649"]], "649 (Responder)"), calc_cors(pseudobulk[["1087"]], "1087 (Resistant)"))

bar_plot <- ggplot(cor_summary, aes(x = Comparison, y = Correlation, fill = Type)) +
  geom_col(width = 0.7) + geom_hline(yintercept = 0.95, linetype = "dashed", color = "#D32F2F", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), vjust = -0.5, size = 4, fontface = "bold") +
  facet_wrap(~Patient) + scale_fill_manual(values = c("Technical" = "#4CAF50", "Biological" = "#2196F3")) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.1)) +
  labs(title = "Technical vs Biological Correlation Summary", subtitle = "Red dashed line = 0.95 threshold",
       x = "", y = "Pearson Correlation") + theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), strip.text = element_text(face = "bold", size = 12),
        legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "plots", "03_correlation_summary_barplot.pdf"), bar_plot, width = 10, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "plots", "03_correlation_summary_barplot.png"), bar_plot, width = 10, height = 6, dpi = 300)
message("Saved: 03_correlation_summary_barplot.pdf/png")

## SAVE TABLE
write.csv(cor_summary, file.path(OUTPUT_DIR, "correlation_summary_649_1087.csv"), row.names = FALSE)

## SUMMARY
message("\n", strrep("=", 60))
message("CORRELATION ANALYSIS COMPLETE")
message(strrep("=", 60))
print(cor_summary)
message("\nInterpretation:")
message("  • Technical correlations >0.95 = excellent reproducibility")
message("  • Lower F1 vs F2 = biological change (therapy response)")
message("Output: ", OUTPUT_DIR)
message(strrep("=", 60), "\n")
