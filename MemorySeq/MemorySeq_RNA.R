### MemorySeq data analysis ###
## RNAseq Data ##
#date: Feb 2025
# Following paper: "Memory Sequencing Reveals Heritable Single-Cell Gene Expression Programs Associated with Distinct Cellular Behaviors"

# Initialize renv in this directory
# Install and initialize renv
install.packages("renv")
renv::init()

# Create main project directories
dir.create("Rnaseq/data", showWarnings = FALSE)      # For raw and processed data
dir.create("Rnaseq/plots", showWarnings = FALSE)     # For all plots/figures
dir.create("Rnaseq/tables", showWarnings = FALSE)    # For result tables and gene lists

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "GenomicFeatures", "biomart"))

# Install CRAN packages using renv
renv::install(c(
  "dplyr",
  "tidyr",
  "data.table",
  "ggplot2"
))

# Load all packages
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(DESeq2)
library(GenomicFeatures)
library(stringr)
library(biomaRt)

# Save snapshot of environment
renv::snapshot()

## Preprocessing count matrix ##

# Rename the sample column in metadata for easier handling
colnames(metadata)[1] <- "Sample"

# Keep essential columns from counts
counts_clean <- counts[, c("Geneid", "Length")]
colnames(counts_clean)[1] <- "gene_id"

# Get sample columns
sample_cols <- grep("^X\\.data", colnames(counts), value = TRUE)

# Clean up sample names
clean_names <- sapply(sample_cols, function(x) {
  sample_name <- str_extract(x, "Mix\\.[0-9]+|P[0-9]+\\.[A-Z]+[0-9]+")
  sample_name <- str_replace(sample_name, "Mix\\.", "Mix_")
  return(sample_name)
})

# Add count data with clean names
counts_clean <- cbind(counts_clean, counts[, sample_cols])
colnames(counts_clean)[3:ncol(counts_clean)] <- clean_names

# Check for sample matching between metadata and counts
missing_samples <- setdiff(metadata$Sample, clean_names)
extra_samples <- setdiff(clean_names, metadata$Sample)

if(length(missing_samples) > 0) {
  message("Samples in metadata but missing from counts: ", 
          paste(missing_samples, collapse = ", "))
}

if(length(extra_samples) > 0) {
  message("Samples in counts but missing from metadata: ",
          paste(extra_samples, collapse = ", "))
}

# Find common samples between metadata and counts
common_samples <- intersect(metadata$Sample, clean_names)

# Filter metadata and counts_clean to keep only common samples
metadata_filtered <- metadata[metadata$Sample %in% common_samples,]
counts_clean_filtered <- counts_clean[, c("gene_id", "Length", common_samples)]

# Order counts_clean columns to match metadata sample order
sample_order <- metadata_filtered$Sample
final_cols <- c("gene_id", "Length", sample_order)

# Ensure all columns in final_cols are actually present in counts_clean_filtered
if(!all(final_cols %in% colnames(counts_clean_filtered))) {
  missing_cols <- final_cols[!final_cols %in% colnames(counts_clean_filtered)]
  stop("The following columns are missing from counts_clean_filtered: ", paste(missing_cols, collapse = ", "))
}

# Reorder counts_clean_filtered
counts_clean_final <- counts_clean_filtered[, final_cols]

# Print the first few rows of the final counts data frame
head(counts_clean_final)

# Print summary
message("\nProcessing complete!")
message("Number of genes: ", nrow(counts_clean))
message("Number of samples: ", ncol(counts_clean) - 2)

# Show the first few rows and columns of the processed data
message("\nPreview of processed data:")
print(head(counts_clean[, 1:10]))

## Annotating Genes with BioMart ##

# 1. Load the Ensembl dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 2. Define the attributes to retrieve
attributes <- c("ensembl_gene_id", "hgnc_symbol")

# 3. Extract Ensembl IDs from counts_clean_final
ensembl_ids <- counts_clean_final$gene_id

# 4. Query BioMart
gene_info <- getBM(
  attributes = attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# 5. Rename the Ensembl gene ID column
gene_info <- gene_info %>%
  dplyr::rename(gene_id = ensembl_gene_id)


# 6. Join the gene information to your counts data
counts_annotated <- left_join(counts_clean_final, gene_info, by = "gene_id")

# 7. Save the annotated data
output_file_annotated <- "RNAseq/data/annotated_counts.csv"
write.csv(counts_annotated, file = output_file_annotated, row.names = FALSE, quote = FALSE)

## Merging Sample Information ##

# 1. Creat sample information
sample_info <- data.frame(
  sample_name = colnames(counts_clean)[-(1:2)],
  is_control = grepl("Mix_", colnames(counts_clean)[-(1:2)])
)

# 2. Merge metadata with sample_info
final_metadata <- left_join(metadata_filtered, sample_info, by = c("Sample" = "sample_name"))

# 3. Save the final metadata
output_file_metadata <- "RNAseq/data/final_metadata.csv"
write.csv(final_metadata, file = output_file_metadata, row.names = FALSE, quote = FALSE)

################################################################################
# Step 1: Initial RNA-seq Data Processing
################################################################################
# Load required libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(GenomicFeatures)
library(IRanges)

# Step 1a: Load and prepare data
counts <- read.csv("RNAseq/data/annotated_counts.csv")
metadata <- read.csv("RNAseq/data/final_metadata.csv")

# Assign control vs sorted sample labels based on `is_control`
metadata$controls <- ifelse(metadata$is_control == TRUE, "Controls", "SortedSample")

# Melt count data to long format
melted_counts <- reshape2::melt(counts, 
                                id.vars = c("gene_id", "hgnc_symbol"), 
                                variable.name = "sampleID",
                                value.name = "counts")

# Merge with metadata
data <- merge(melted_counts, 
              metadata[, c("Sample", "controls", "Sortedcell")],
              by.x = "sampleID", 
              by.y = "Sample")

# Step 1b: Calculate TPM values
message("Calculating TPM values...")

# Load gene lengths from counts file
length_table <- counts[, c("gene_id", "Length")]
colnames(length_table) <- c("gene_id", "length")

# Add gene lengths and calculate TPM
data <- left_join(data, length_table, by = "gene_id")

data <- data %>%
  group_by(sampleID) %>%
  mutate(
    totalMappedReads = sum(counts, na.rm = TRUE),
    rpm = 1e6 * counts / totalMappedReads,
    rpk = counts / (length / 1000),  # Normalize by gene length in kb
    rpkScalePerMillion = sum(rpk, na.rm = TRUE) / 1e6,
    tpm = rpk / rpkScalePerMillion
  ) %>%
  ungroup()

# Step 1c: Filter low-quality samples and low-expression genes
message("Filtering samples and genes...")

# Remove samples with less than 500,000 reads
data_filtered <- data %>%
  group_by(sampleID) %>%
  mutate(totalMappedReads = sum(counts, na.rm = TRUE)) %>%
  filter(totalMappedReads > 0.5e6) %>%
  ungroup()

# Calculate mean expression per gene
gene_stats <- data_filtered %>%
  group_by(gene_id) %>%
  summarize(meanRPM = mean(rpm, na.rm = TRUE)) %>%
  mutate(keep = meanRPM > 1)  # Keep genes with mean RPM > 1

# Filter genes using dplyr syntax
gene_stats_filtered <- filter(gene_stats, keep)

# Ensure selection works without errors
data_filtered <- inner_join(data_filtered, gene_stats_filtered, by = "gene_id") %>%
  dplyr::select(-meanRPM, -keep)

# Step 1d: Save processed data
message("Saving processed data...")
dir.create("results", showWarnings = FALSE)
write.csv(data_filtered, "RNAseq/data/processed_counts.csv", row.names = FALSE, quote = FALSE)

# Step 1e: Generate QC plots
message("Generating QC plots...")

# Load ggplot2
library(ggplot2)

# Library size distribution 
p1 <- ggplot(data_filtered, aes(x = sampleID, y = totalMappedReads / 1e6, fill = controls)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("SortedSample" = "blue", "Controls" = "pink")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Million mapped reads", title = "Library Size Distribution")

# Save Library Size Distribution plot
ggsave("RNAseq/plots/library_size_distribution.pdf", p1, width = 12, height = 6)
ggsave("RNAseq/plots/library_size_distribution.png", p1, width = 12, height = 6, dpi = 300)

# Epression distribution 
p2 <- ggplot(data_filtered, aes(x = log2(tpm + 1), color = controls)) +
  geom_density(size = 1) +
  scale_color_manual(values = c("SortedSample" = "royalblue3", "Controls" = "pink2")) +
  theme_minimal() +
  labs(x = "log2(Average TPM + 1)", y = "Density", title = "Expression Distribution")

# Save Expression Distribution plot
ggsave("RNAseq/plots/expression_distribution.pdf", p2, width = 8, height = 6)
ggsave("RNAseq/plots/expression_distribution.png", p2, width = 8, height = 6, dpi = 300)

# Summary statistics
message("\nProcessing complete! Summary statistics:")
message("Number of samples: ", length(unique(data_filtered$sampleID)))
message("Number of genes after filtering: ", length(unique(data_filtered$gene_id))) # 14035
message("Number of control samples: ", length(unique(filter(data_filtered, controls == "Controls")$sampleID)))
message("Number of sorted samples: ", length(unique(filter(data_filtered, controls == "SortedSample")$sampleID)))

# Step 1f: Save entire environment as an RDS file
save.image(file = "RNAseq/data/processed_data_step1.RData")
################################################################################
# Step 2: Batch Effect Analysis and Correction
################################################################################
# Step 2: Batch Effect Analysis and Correction using limma
library(dplyr)
library(ggplot2)
library(limma)
library(DESeq2)
library(tibble) 
library(tidyr)


# Load processed data if not in environment
if(!exists("data_filtered")) {
  load("data/processed/step1_processed_data.RData")
}

# 2a. Prepare data for batch analysis
# Method 1: Using base R to create expression matrix
# First let's check the duplicates
duplicates <- data_filtered %>%
  dplyr::group_by(gene_id, sampleID) %>%
  dplyr::summarise(n = n(), .groups = 'drop') %>%
  dplyr::filter(n > 1)

print("Number of duplicate combinations:")
print(nrow(duplicates))

# Create expression matrix with mean values for duplicates
expr_matrix <- data_filtered %>%
  dplyr::group_by(gene_id, sampleID) %>%
  dplyr::summarise(tpm = mean(tpm), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = sampleID,
                     values_from = tpm)

# Convert to matrix format
rownames(expr_matrix) <- expr_matrix$gene_id
expr_matrix <- expr_matrix[,-1]
expr_matrix <- as.matrix(expr_matrix)

# Check the result
print("Dimensions of expression matrix:")
print(dim(expr_matrix))
print("First few rows and columns:")
print(expr_matrix[1:5, 1:5])

# Let's also check if we still have any NA values
print("Number of NA values:")
print(sum(is.na(expr_matrix)))

# Save the expression matrix for the next steps
save(expr_matrix, file = "RNAseq/data/expression_matrix.RData")

# Convert to matrix
expr_matrix <- as.matrix(expr_matrix)

# Log transform TPM values
log_expr <- log2(expr_matrix + 1)

# 2b. Principal Component Analysis before correction
pca_before <- prcomp(t(log_expr))
pca_df_before <- data.frame(
  PC1 = pca_before$x[,1],
  PC2 = pca_before$x[,2],
  Batch = metadata$Batch,
  SampleType = metadata$controls
)

# Create PCA plot before correction
p1 <- ggplot(pca_df_before, 
             aes(x = PC1, y = PC2, color = factor(Batch), shape = SampleType)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA Plot Before Batch Correction") +
  theme(plot.title = element_text(hjust = 0.5))

# 2c. Test for batch effects
batch_test <- aov(PC1 ~ factor(Batch), data = pca_df_before)
batch_pvalue <- summary(batch_test)[[1]]["Pr(>F)"][1,1]

# 2d. Perform batch correction using limma's removeBatchEffect
# Create design matrix for sample groups
design <- model.matrix(~metadata$controls)

# Perform batch correction
corrected_matrix <- removeBatchEffect(
  log_expr,
  batch = metadata$Batch,
  design = design
)

# Convert back to TPM scale
corrected_tpm <- 2^corrected_matrix - 1

# 2e. PCA after correction
pca_after <- prcomp(t(corrected_matrix))
pca_df_after <- data.frame(
  PC1 = pca_after$x[,1],
  PC2 = pca_after$x[,2],
  Batch = metadata$Batch,
  SampleType = metadata$controls
)

# Create PCA plot after correction
p2 <- ggplot(pca_df_after, 
             aes(x = PC1, y = PC2, color = factor(Batch), shape = SampleType)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA Plot After Batch Correction") +
  theme(plot.title = element_text(hjust = 0.5))

# Save plots
dir.create("plots/batch_effects", recursive = TRUE, showWarnings = FALSE)
ggsave("RNAseq/plots/batch_effects/pca_before_correction.pdf", p1, width = 8, height = 6)
ggsave("RNAseq/plots/batch_effects/pca_before_correction.png", p1, width = 8, height = 6,dpi = 300)

ggsave("RNAseqplots/batch_effects/pca_after_correction.pdf", p2, width = 8, height = 6)
ggsave("RNAseq/plots/batch_effects/pca_after_correction.png", p2, width = 8, height = 6, dpi = 300)

################################################################################
# Step 3: Separate Variability Analysis for Controls and Sorted Samples
################################################################################

# Complete Variability Analysis following paper exactly
library(dplyr)
library(ggplot2)
library(moments)  # for skewness and kurtosis
library(RColorBrewer)

# Calculate ALL metrics mentioned in paper for both groups
metrics_complete <- data_filtered %>%
  group_by(gene_id, controls, hgnc_symbol) %>%
  summarize(
    CV_tpm = sd(tpm) / mean(tpm),
    Mean_tpm = mean(tpm),
    skewness = moments::skewness(tpm),
    kurtosis = moments::kurtosis(tpm),
    .groups = 'drop'
  )

# Separate Controls and Sorted Samples
control_metrics <- filter(metrics_complete, controls == "Controls")
sorted_metrics <- filter(metrics_complete, controls == "SortedSample")

# Fit Poisson regression models
# model_controls <- glm(CV_tpm ~ log2(Mean_tpm), 
#                       family = poisson(link="log"), 
#                       data = control_metrics)

model_sorted <- glm(CV_tpm ~ log2(Mean_tpm), 
                    family = poisson(link = "log"), 
                    data = sorted_metrics,
                    na.action = na.exclude)

# Calculate residuals with all metrics
fit_data_controls <- data.frame(
  gene_id = control_metrics$gene_id,
  hgnc_symbol = control_metrics$hgnc_symbol,
  group = "Controls",
  Mean_tpm = control_metrics$Mean_tpm,
  CV_tpm = control_metrics$CV_tpm,
  skewness = control_metrics$skewness,
  kurtosis = control_metrics$kurtosis,
  fitted_y = fitted(model_controls),
  resid = resid(model_controls)
)

fit_data_sorted <- data.frame(
  gene_id = sorted_metrics$gene_id,
  hgnc_symbol = sorted_metrics$hgnc_symbol,
  group = "SortedSample",
  Mean_tpm = sorted_metrics$Mean_tpm,
  CV_tpm = sorted_metrics$CV_tpm,
  skewness = sorted_metrics$skewness,
  kurtosis = sorted_metrics$kurtosis,
  fitted_y = fitted(model_sorted),
  resid = resid(model_sorted)
)

# Apply paper's exact cutoffs
mean_tpm_cutoff <- 1.5  # From paper
percentile_cutoff <- 0.98  # From paper

# Calculate cutoffs
cutoff_controls <- quantile(
  filter(fit_data_controls, Mean_tpm > mean_tpm_cutoff)$resid, 
  percentile_cutoff
)

cutoff_sorted <- quantile(
  filter(fit_data_sorted, Mean_tpm > mean_tpm_cutoff)$resid, 
  percentile_cutoff
)

# Identify heritable genes
heritable_genes_controls <- fit_data_controls %>%
  filter(Mean_tpm > mean_tpm_cutoff & resid > cutoff_controls)

heritable_genes_sorted <- fit_data_sorted %>%
  filter(Mean_tpm > mean_tpm_cutoff & resid > cutoff_sorted)

# Calculate Cook's distance for pairwise correlations
calculate_cooks_distance <- function(gene1_expr, gene2_expr) {
  model <- lm(gene2_expr ~ gene1_expr)
  return(cooks.distance(model))
}

# Get expression data for heritable genes
# heritable_expr <- data_filtered %>%
#   filter(gene_id %in% heritable_genes_sorted$gene_id,
#          controls == "SortedSample") %>%
#   select(gene_id, sampleID, tpm) %>%
#   pivot_wider(names_from = sampleID, values_from = tpm)

heritable_expr <- data_filtered %>%
  filter(gene_id %in% heritable_genes_sorted$gene_id,
         controls == "SortedSample") %>%
  dplyr::select(gene_id, sampleID, tpm) %>%
  pivot_wider(names_from = sampleID, values_from = tpm)


# Calculate correlation matrix and Cook's distance for pairs
cor_matrix <- cor(t(heritable_expr[,-1]), method = "pearson")

# Create plots
# 1. CV vs Mean with both controls and sorted samples
p1 <- ggplot() +
  stat_bin2d(data = rbind(
    filter(fit_data_controls, !(gene_id %in% heritable_genes_controls$gene_id)),
    filter(fit_data_sorted, !(gene_id %in% heritable_genes_sorted$gene_id))),
    aes(x = log2(Mean_tpm), y = CV_tpm),
    bins = 50) +
  scale_fill_gradientn(colors = colorRampPalette(c("black", "#56B4E9"))(100)) +
  geom_point(data = rbind(heritable_genes_controls, heritable_genes_sorted),
             aes(x = log2(Mean_tpm), y = CV_tpm, color = group),
             size = 0.5, alpha = 0.5) +
  theme_classic() +
  facet_wrap(~group) +
  labs(x = "log2(Mean TPM)", y = "Coefficient of Variation")

# 2. Skewness vs CV
p2 <- ggplot(rbind(fit_data_controls, fit_data_sorted),
             aes(x = CV_tpm, y = skewness, color = group)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "CV", y = "Skewness")

# 3. Distribution of residuals
p3 <- ggplot() +
  geom_histogram(data = fit_data_controls, aes(x = resid, fill = "Controls"), 
                 alpha = 0.5, bins = 100) +
  geom_histogram(data = fit_data_sorted, aes(x = resid, fill = "Sorted"), 
                 alpha = 0.5, bins = 100) +
  geom_vline(xintercept = c(cutoff_controls, cutoff_sorted), 
             linetype = "dashed", color = c("blue", "red")) +
  theme_classic() +
  labs(x = "Residuals", y = "Count")

# Save plots
dir.create("RNAseq/plots/variability", recursive = TRUE, showWarnings = FALSE)
ggsave("RNAseq/plots/variability/cv_vs_mean_complete.pdf", p1, width = 10, height = 5)
ggsave("RNAseq/plots/variability/skewness_vs_cv.pdf", p2, width = 6, height = 6)
ggsave("RNAseq/plots/variability/residuals_distribution_combined.pdf", p3, width = 8, height = 6)

# Save results
dir.create("RNAseq/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(heritable_genes_sorted %>% 
            arrange(desc(resid)), 
          "RNAseq/tables/heritable_genes_sorted.csv", 
          row.names = FALSE)
write.csv(heritable_genes_controls %>% 
            arrange(desc(resid)), 
          "RNAseq/tables/heritable_genes_controls.csv", 
          row.names = FALSE)

# Print summary statistics
message("\nComplete Variability Analysis Summary:")
message("Controls:")
message("  Number of heritable genes: ", nrow(heritable_genes_controls))
message("  Residual cutoff value: ", round(cutoff_controls, 3))
message("\nSorted Samples:")
message("  Number of heritable genes: ", nrow(heritable_genes_sorted))
message("  Residual cutoff value: ", round(cutoff_sorted, 3))

# Save complete analysis results
variability_complete <- list(
  metrics = metrics_complete,
  controls = list(
    model = model_controls,
    heritable_genes = heritable_genes_controls,
    cutoff = cutoff_controls
  ),
  sorted = list(
    model = model_sorted,
    heritable_genes = heritable_genes_sorted,
    cutoff = cutoff_sorted
  ),
  correlation = list(
    matrix = cor_matrix
  )
)

save(variability_complete, 
     file = "RNAseq/data/variability_complete_results.RData")

################################################################################
# Generate all plots following paper's figures
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)

# Create directories for different plot types
plot_dirs <- c("main_figures", "supplementary", "qc", "correlation")
for(dir in plot_dirs) {
  dir.create(file.path("RNAseq/plots", dir), recursive = TRUE, showWarnings = FALSE)
}

###############################################
# Figure 1: Initial QC and Variability Plots
###############################################

# 1A. Library Size Distribution
p1a <- ggplot(data_filtered, 
              aes(x = sampleID, y = totalMappedReads/1e6, fill = controls)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Controls" = "pink", "SortedSample" = "blue")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Million mapped reads", 
       title = "Library Size Distribution")

# 1B. Expression Distribution
p1b <- ggplot(data_filtered, 
              aes(x = log2(tpm + 1), color = controls)) +
  geom_density() +
  scale_color_manual(values = c("Controls" = "pink", "SortedSample" = "blue")) +
  theme_classic() +
  labs(x = "log2(TPM + 1)", y = "Density",
       title = "Expression Distribution")

# 1C. Expression distributions for top variable genes
top_genes <- head(heritable_genes_sorted$hgnc_symbol, 10)
p1c <- data_filtered %>%
  filter(hgnc_symbol %in% top_genes) %>%
  ggplot(aes(x = log2(tpm + 1), color = controls)) +
  geom_density() +
  geom_rug() +
  facet_wrap(~hgnc_symbol, scales = "free_y", ncol = 1) +
  theme_classic() +
  labs(x = "log2(TPM + 1)", y = "Density")

# 1D. CV vs Mean Expression
p1d <- ggplot() +
  stat_bin2d(data = rbind(
    filter(fit_data_controls, !(gene_id %in% heritable_genes_controls$gene_id)),
    filter(fit_data_sorted, !(gene_id %in% heritable_genes_sorted$gene_id))),
    aes(x = log2(Mean_tpm), y = CV_tpm),
    bins = 50) +
  scale_fill_gradientn(colors = colorRampPalette(c("black", "#56B4E9"))(100)) +
  geom_point(data = rbind(heritable_genes_controls, heritable_genes_sorted),
             aes(x = log2(Mean_tpm), y = CV_tpm),
             color = 'green', size = 0.5, alpha = 0.5) +
  theme_classic() +
  facet_wrap(~group) +
  labs(x = "log2(Mean TPM)", y = "CV")

# 1E. Residuals Distribution
p1e <- ggplot() +
  geom_histogram(data = fit_data_controls, 
                 aes(x = resid, fill = "Controls"), 
                 alpha = 0.5, bins = 100) +
  geom_histogram(data = fit_data_sorted, 
                 aes(x = resid, fill = "Sorted"), 
                 alpha = 0.5, bins = 100) +
  geom_vline(xintercept = c(cutoff_controls, cutoff_sorted), 
             linetype = "dashed", 
             color = c("blue", "red")) +
  theme_classic() +
  labs(x = "Residuals", y = "Count")

###############################################
# Supplementary Figure S1
###############################################

# S1A. Skewness vs CV
ps1a <- ggplot(rbind(fit_data_controls, fit_data_sorted),
               aes(x = CV_tpm, y = skewness, color = group)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "CV", y = "Skewness")

# S1B. Kurtosis vs CV
ps1b <- ggplot(rbind(fit_data_controls, fit_data_sorted),
               aes(x = CV_tpm, y = kurtosis, color = group)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "CV", y = "Kurtosis")

# S1C. Mean vs Variance
ps1c <- ggplot(rbind(fit_data_controls, fit_data_sorted),
               aes(x = log2(Mean_tpm), 
                   y = log2((CV_tpm * Mean_tpm)^2), 
                   color = group)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "log2(Mean TPM)", y = "log2(Variance)")

###############################################
# Figure 5: Correlation Analysis
###############################################

# 5A. Correlation Heatmap
# Calculate correlation matrix for heritable genes
heritable_expr <- data_filtered %>%
  filter(gene_id %in% heritable_genes_sorted$gene_id,
         controls == "SortedSample") %>%
  select(gene_id, sampleID, tpm) %>%
  pivot_wider(names_from = sampleID, values_from = tpm)

cor_matrix <- cor(t(heritable_expr[,-1]), method = "pearson")

# Create correlation heatmap
pdf("RNAseq/plots/correlation/correlation_heatmap.pdf", width = 12, height = 12)
pheatmap(cor_matrix,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Gene Expression Correlation")
dev.off()

# Save all main figure plots
ggsave("RNAseq/plots/main_figures/Fig1A_library_size.pdf", p1a, width = 10, height = 6)
ggsave("RNAseq/plots/main_figures/Fig1B_expression_dist.pdf", p1b, width = 8, height = 6)
ggsave("RNAseq/plots/main_figures/Fig1C_top_genes_dist.pdf", p1c, width = 8, height = 15)
ggsave("RNAseq/plots/main_figures/Fig1D_cv_mean.pdf", p1d, width = 12, height = 6)
ggsave("RNAseq/plots/main_figures/Fig1E_residuals.pdf", p1e, width = 8, height = 6)

# Save supplementary figure plots
ggsave("RNAseq/plots/supplementary/FigS1A_skewness.pdf", ps1a, width = 8, height = 6)
ggsave("RNAseq/plots/supplementary/FigS1B_kurtosis.pdf", ps1b, width = 8, height = 6)
ggsave("RNAseq/plots/supplementary/FigS1C_mean_variance.pdf", ps1c, width = 8, height = 6)

# Create combined figure panels
# Figure 1 panel
fig1_panel <- grid.arrange(p1a, p1b, p1c, p1d, p1e,
                           layout_matrix = rbind(c(1,2),
                                                 c(3,4),
                                                 c(5,5)),
                           heights = c(1, 2, 1))
ggsave("RNAseq/plots/main_figures/Figure1_complete.pdf", 
       fig1_panel, width = 15, height = 20)

# Supplementary Figure S1 panel
figS1_panel <- grid.arrange(ps1a, ps1b, ps1c,
                            ncol = 2)
ggsave("RNAseq/plots/supplementary/FigureS1_complete.pdf", 
       figS1_panel, width = 12, height = 8)

# Save plot objects for future use
save(list = ls(pattern = "^p"), 
     file = "RNAseq/plots/plot_objects.RData")

################################################################################
#Correlation Analysis
################################################################################

# Correlation Analysis following paper's methodology
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dendextend)

# Load variability results if not in environment
if(!exists("variability_results")) {
  load("RNAseq/data/variability_analysis_results.RData")
}

# Create expression matrix for variable genes
# For Sorted Samples
sorted_expr_matrix <- data_filtered %>%
  filter(controls == "SortedSample",
         gene_id %in% variability_results$sorted$variable_genes$gene_id) %>%
  select(gene_id, sampleID, tpm) %>%
  pivot_wider(names_from = sampleID, 
              values_from = tpm) %>%
  column_to_rownames("gene_id")

# Calculate correlation matrix
cor_matrix_sorted <- cor(t(sorted_expr_matrix), method = "pearson")

# Add gene symbols for labeling
gene_symbols <- variability_results$sorted$variable_genes$hgnc_symbol
names(gene_symbols) <- variability_results$sorted$variable_genes$gene_id
rownames(cor_matrix_sorted) <- gene_symbols[rownames(cor_matrix_sorted)]
colnames(cor_matrix_sorted) <- gene_symbols[colnames(cor_matrix_sorted)]

# Create correlation heatmap (Figure 5-style)
# Define colors
col_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Perform hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix_sorted), method = "complete")

# Create heatmap
pdf("RNAseq/plots/correlation/correlation_heatmap.pdf", width = 12, height = 12)
pheatmap(cor_matrix_sorted,
         color = col_palette,
         breaks = seq(-1, 1, length.out = 101),
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Gene Expression Correlation in Sorted Samples")
dev.off()

# Create clustered gene groups
# Cut tree to get clusters (adjust h or k based on your data)
clusters <- cutree(hc, k = 3)  # or use h parameter for height-based cutting

# Add cluster information to gene list
gene_clusters <- data.frame(
  gene_id = names(clusters),
  cluster = clusters
) %>%
  left_join(variability_results$sorted$variable_genes, by = "gene_id")

# Save cluster assignments
write.csv(gene_clusters %>% 
            arrange(cluster, desc(resid)), 
          "RNAseq/tables/gene_clusters.csv", 
          row.names = FALSE)

# Create scatter plots for highly correlated gene pairs
# Find top correlated pairs
cor_long <- cor_matrix_sorted %>%
  as.data.frame() %>%
  rownames_to_column("gene1") %>%
  pivot_longer(-gene1, names_to = "gene2", values_to = "correlation") %>%
  filter(gene1 < gene2) %>%  # Remove duplicates and self-correlations
  arrange(desc(abs(correlation))) %>%
  head(10)  # Top 10 correlations

# Create scatter plots for top pairs
for(i in 1:nrow(cor_long)) {
  pair_data <- data_filtered %>%
    filter(controls == "SortedSample",
           hgnc_symbol %in% c(cor_long$gene1[i], cor_long$gene2[i])) %>%
    select(sampleID, hgnc_symbol, tpm) %>%
    pivot_wider(names_from = hgnc_symbol, values_from = tpm)
  
  p <- ggplot(pair_data, 
              aes_string(x = paste0("`", cor_long$gene1[i], "`"), 
                         y = paste0("`", cor_long$gene2[i], "`"))) +
    geom_point() +
    theme_classic() +
    labs(title = sprintf("r = %.3f", cor_long$correlation[i]))
  
  ggsave(sprintf("RNAseq/plots/correlation/gene_pair_%d.pdf", i), 
         p, width = 6, height = 6)
}

# Print summary
message("\nCorrelation Analysis Summary:")
message("Number of variable genes analyzed: ", nrow(sorted_expr_matrix))
message("Number of clusters identified: ", length(unique(clusters)))
message("\nTop correlations:")
print(head(cor_long, 5))

# Save correlation results
correlation_results <- list(
  correlation_matrix = cor_matrix_sorted,
  clustering = list(
    hclust_obj = hc,
    clusters = clusters
  ),
  top_correlations = cor_long
)

save(correlation_results, 
     file = "RNAseq/data/correlation_analysis_results.RData")



