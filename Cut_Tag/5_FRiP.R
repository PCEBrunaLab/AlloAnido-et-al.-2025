#!/usr/bin/env Rscript

# Peak Analysis from Summary File
# Analyzes peak number, width, and FRiP for R85 CUT&Tag data

library(tidyverse)
library(GenomicRanges)
library(Rsamtools)
library(ggpubr)

# Set paths
projPath <- "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise"
peakPath <- paste0(projPath, "/peakCalling/MACS2")
bamPath <- paste0(projPath, "/bam_files_dedup")

# Read peak summary
peak_summary <- read.table(paste0(peakPath, "/peak_summary_separate_reps.txt"), 
                           header = TRUE, sep = "\t")

# Filter for default q=0.1 results only
peak_summary_filtered <- peak_summary %>%
  filter(Q_value == 0.1, Mode == "default") %>%
  distinct(Replicate, Comparison, Target, .keep_all = TRUE)

# Your color scheme
target_colors <- c("H3K4me3"="#E41A1C", "KDM5A"="#4DAF4A", "KDM5B"="#984EA3")

# =====================================================
# 1. PEAK NUMBER ANALYSIS
# =====================================================

print("=== PEAK NUMBER ANALYSIS ===")

# Summary statistics for peak numbers
peak_number_summary <- peak_summary_filtered %>%
  group_by(Target, Comparison) %>%
  summarise(
    Mean_Peaks = mean(Peak_Count),
    SD_Peaks = sd(Peak_Count),
    Rep1_Peaks = Peak_Count[Replicate == "Rep1"],
    Rep2_Peaks = Peak_Count[Replicate == "Rep2"],
    CV = (SD_Peaks / Mean_Peaks) * 100  # Coefficient of variation
  ) %>%
  mutate(Treatment = sub("_vs_.*", "", Comparison)) %>%
  mutate(Control = sub(".*_vs_", "", Comparison))

print(peak_number_summary)

# Plot peak numbers
# Plot peak numbers - CORRECTED VERSION
p_peak_number <- peak_summary_filtered %>%
  mutate(Treatment = sub("_vs_.*", "", Comparison)) %>%
  mutate(Treatment = factor(Treatment, levels = c("D7_Cis", "D14_Cis", "D14_C70"))) %>%
  ggplot(aes(x = Target, y = Peak_Count, fill = Target)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_point(aes(shape = as.factor(Replicate)), 
             position = position_dodge(width = 0.3), size = 3) +
  facet_wrap(~Treatment, scales = "free_y") +
  scale_fill_manual(values = target_colors) +
  scale_shape_manual(values = c(16, 17), name = "Replicate") +
  theme_bw(base_size = 10) +
  labs(title = "Peak Numbers by Treatment and Target",
       y = "Number of Peaks",
       x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom")

ggsave("peak_numbers_analysis.pdf", p_peak_number, width = 10, height = 6)

# =====================================================
# 2. PEAK WIDTH ANALYSIS
# =====================================================

print("\n=== PEAK WIDTH ANALYSIS ===")

# Function to read peak widths
get_peak_widths <- function(sample_id, condition, target, rep) {
  # Map to file path based on your directory structure
  if (condition == "D7_Cis" || condition == "D14_Cis" || condition == "D14_C70") {
    peak_file <- paste0(peakPath, "/human/Untreated_control/", condition, "/", 
                        target, "/rep", rep, "/macs2_peak_q0.1_peaks.narrowPeak.bed")
  } else {
    return(NULL)  # Skip if not in expected conditions
  }
  
  if (file.exists(peak_file)) {
    peaks <- read.table(peak_file, header = FALSE, sep = "\t")
    widths <- abs(peaks$V3 - peaks$V2)
    return(data.frame(
      Sample_ID = sample_id,
      Condition = condition,
      Target = target,
      Replicate = rep,
      Width = widths
    ))
  } else {
    return(NULL)
  }
}

# Collect all peak widths
all_widths <- list()
for (i in 1:nrow(peak_summary_filtered)) {
  row <- peak_summary_filtered[i,]
  condition <- strsplit(as.character(row$Comparison), "_vs_")[[1]][1]
  rep_num <- as.numeric(gsub("Rep", "", row$Replicate))
  sample_id <- strsplit(as.character(row$Sample_IDs), "_vs_")[[1]][1]
  
  widths <- get_peak_widths(sample_id, condition, row$Target, rep_num)
  if (!is.null(widths)) {
    all_widths[[length(all_widths) + 1]] <- widths
  }
}

# Combine all width data
if (length(all_widths) > 0) {
  peak_width_data <- bind_rows(all_widths)
  
  # Summary statistics
  width_summary <- peak_width_data %>%
    group_by(Target, Condition) %>%
    summarise(
      Median_Width = median(Width),
      Mean_Width = mean(Width),
      Min_Width = min(Width),
      Max_Width = max(Width),
      Q25 = quantile(Width, 0.25),
      Q75 = quantile(Width, 0.75)
    )
  
  print(width_summary)
  
  # Plot peak width distributions
  p_peak_width <- peak_width_data %>%
    ggplot(aes(x = Target, y = Width, fill = Target)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    facet_wrap(~Condition) +
    scale_y_log10(breaks = c(100, 500, 1000, 5000, 10000),
                  labels = c("100", "500", "1K", "5K", "10K")) +
    scale_fill_manual(values = target_colors) +
    theme_bw(base_size = 10) +
    labs(title = "Peak Width Distribution",
         y = "Peak Width (bp, log scale)",
         x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(hjust = 0.5, size = 12),
          legend.position = "bottom")
  
  ggsave("peak_width_distribution.pdf", p_peak_width, width = 10, height = 6)
}

# =====================================================
# 3. FRiP CALCULATION
# =====================================================

print("\n=== FRiP CALCULATION ===")

# Read QC summary for total mapped reads
qc_data <- read.table(paste0(projPath, "/QC/qc_summary.txt"), 
                      header = TRUE, sep = "\t")


# Function to calculate FRiP 
calculate_frip <- function(peak_file, bam_file, total_mapped_reads) {
  if (!file.exists(peak_file) || !file.exists(bam_file)) {
    return(NA)
  }
  
  # Read peaks
  peaks <- read.table(peak_file, header = FALSE, sep = "\t")
  peak_gr <- GRanges(seqnames = peaks$V1,
                     ranges = IRanges(start = peaks$V2, end = peaks$V3))
  
  # Count reads in peaks using Rsamtools
  param <- ScanBamParam(which = peak_gr, 
                        flag = scanBamFlag(isUnmappedQuery = FALSE))
  
  # Get total count across all peaks - FIXED
  bam_file_handle <- BamFile(bam_file)
  reads_in_peaks <- sum(countBam(bam_file_handle, param = param)$records)
  
  # Calculate FRiP
  frip <- (reads_in_peaks / total_mapped_reads) * 100
  return(frip)
}

# First, define the corrected function
calculate_frip <- function(peak_file, bam_file, total_mapped_reads) {
  if (!file.exists(peak_file) || !file.exists(bam_file)) {
    return(NA)
  }
  
  # Read peaks
  peaks <- read.table(peak_file, header = FALSE, sep = "\t")
  peak_gr <- GRanges(seqnames = peaks$V1,
                     ranges = IRanges(start = peaks$V2, end = peaks$V3))
  
  # Count reads in peaks using Rsamtools
  param <- ScanBamParam(which = peak_gr, 
                        flag = scanBamFlag(isUnmappedQuery = FALSE))
  
  # Get total count across all peaks - FIXED
  bam_file_handle <- BamFile(bam_file)
  reads_in_peaks <- sum(countBam(bam_file_handle, param = param)$records)
  
  # Calculate FRiP
  frip <- (reads_in_peaks / total_mapped_reads) * 100
  return(frip)
}

# Now calculate FRiP for each sample
frip_results <- peak_summary_filtered %>%
  rowwise() %>%
  mutate(
    Sample_ID = strsplit(Sample_IDs, "_vs_")[[1]][1],
    Condition = strsplit(Comparison, "_vs_")[[1]][1],
    Rep_Num = as.numeric(gsub("Rep", "", Replicate)),
    Peak_File = paste0(peakPath, "/human/Untreated_control/", 
                       Condition, "/", Target, "/rep", Rep_Num, 
                       "/macs2_peak_q0.1_peaks.narrowPeak.bed"),
    BAM_File = paste0(bamPath, "/", Sample_ID, "_dedup.bam")
  ) %>%
  left_join(qc_data %>% select(Sample, Mapped_Reads), 
            by = c("Sample_ID" = "Sample")) %>%
  rowwise() %>%
  mutate(
    FRiP = calculate_frip(Peak_File, BAM_File, Mapped_Reads)
  )

# Then calculate the summary
frip_summary <- frip_results %>%
  group_by(Target, Condition) %>%
  summarise(
    Mean_FRiP = mean(FRiP, na.rm = TRUE),
    SD_FRiP = sd(FRiP, na.rm = TRUE),
    Rep1_FRiP = FRiP[Replicate == "Rep1"],
    Rep2_FRiP = FRiP[Replicate == "Rep2"]
  )

print(frip_summary)



# Plot FRiP
p_frip <- frip_results %>%
  ggplot(aes(x = Target, y = FRiP, fill = Target)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_point(aes(shape = as.factor(Replicate)), 
             position = position_dodge(width = 0.3), size = 3) +
  facet_wrap(~Condition) +
  scale_fill_manual(values = target_colors) +
  scale_shape_manual(values = c(16, 17), name = "Replicate") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "orange", alpha = 0.5) +
  theme_bw(base_size = 10) +
  labs(title = "Fraction of Reads in Peaks (FRiP)",
       y = "FRiP (%)",
       x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom") +
  annotate("text", x = 3.5, y = 5, label = "IgG threshold", size = 3, color = "red") +
  annotate("text", x = 3.5, y = 20, label = "Good signal", size = 3, color = "orange")

ggsave("frip_analysis.pdf", p_frip, width = 10, height = 6)

# =====================================================
# 4. COMBINED SUMMARY REPORT
# =====================================================

print("\n=== COMBINED SUMMARY REPORT ===")

# Create combined summary
combined_summary <- peak_number_summary %>%
  select(Target, Treatment, Mean_Peaks, CV) %>%
  left_join(
    frip_summary %>% select(Target, Condition, Mean_FRiP),
    by = c("Target" = "Target", "Treatment" = "Condition")
  )

print(combined_summary)

# Quality assessment
quality_assessment <- combined_summary %>%
  mutate(
    Peak_Quality = case_when(
      Target == "H3K4me3" & Mean_Peaks > 100 ~ "Good",
      Target %in% c("KDM5A", "KDM5B") & Mean_Peaks > 50 ~ "Good",
      TRUE ~ "Check"
    ),
    FRiP_Quality = case_when(
      Target == "H3K4me3" & Mean_FRiP > 30 ~ "Excellent",
      Target == "H3K4me3" & Mean_FRiP > 20 ~ "Good",
      Target %in% c("KDM5A", "KDM5B") & Mean_FRiP > 10 ~ "Good",
      Mean_FRiP < 5 ~ "Poor",
      TRUE ~ "Moderate"
    ),
    Replicate_Consistency = case_when(
      CV < 20 ~ "Excellent",
      CV < 50 ~ "Good",
      TRUE ~ "Poor"
    )
  )

print("\n=== QUALITY ASSESSMENT ===")
print(quality_assessment)

# Save all summaries
write.table(peak_number_summary, "peak_number_summary.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(frip_summary, "frip_summary.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(quality_assessment, "quality_assessment.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Create combined figure
combined_plot <- ggarrange(p_peak_number, p_frip,
                           ncol = 1, nrow = 2,
                           common.legend = TRUE,
                           legend = "bottom")

ggsave("peak_analysis_combined.pdf", combined_plot, width = 12, height = 10)

print("\n=== Analysis Complete ===")
print("Output files:")
print("  - peak_numbers_analysis.pdf")
print("  - peak_width_distribution.pdf") 
print("  - frip_analysis.pdf")
print("  - peak_analysis_combined.pdf")
print("  - peak_number_summary.txt")
print("  - frip_summary.txt")
print("  - quality_assessment.txt")