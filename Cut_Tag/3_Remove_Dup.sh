#!/bin/bash
#SBATCH --job-name=CUTTag_RemoveDup_Compare_R85
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%j.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%j.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL

# Load modules
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc

# Activate mamba environment
conda activate r_env

# Load additional modules
module load SAMtools/1.11
module load picard-tools/2.23.8
module load R/4.1.0

# Define directories - Updated for R85 project
BAM_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/bam_files
OUT_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/bam_files_dedup
METRICS_DIR=${OUT_DIR}/metrics
LOG_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/4_Remove_Dup
QC_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/QC_final
COMP_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/QC_dedup_comparison

# Create output directories
mkdir -p $OUT_DIR
mkdir -p $METRICS_DIR
mkdir -p $LOG_DIR
mkdir -p $COMP_DIR

# Start fresh debug log
echo "=== Starting R85 Duplicate Removal ===" > ${LOG_DIR}/debug.log
date >> ${LOG_DIR}/debug.log

# Get sample list
SAMPLE_LIST=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/sample_list.txt

# Clean up the sample list to remove any trailing spaces or hidden characters
TMP_SAMPLE_LIST=${LOG_DIR}/clean_sample_list.txt

# Remove any carriage returns and trailing spaces
tr -d '\r' < "$SAMPLE_LIST" | sed 's/[[:space:]]*$//' > "$TMP_SAMPLE_LIST"

# Debug: Print cleaned sample list
echo "=== Cleaned sample list contents ===" >> ${LOG_DIR}/debug.log
cat "$TMP_SAMPLE_LIST" >> ${LOG_DIR}/debug.log
echo "=== End of sample list ===" >> ${LOG_DIR}/debug.log

# Initialize post-dedup summary file
echo -e "Sample\tTotal_Reads_PostDedup\tMapped_PostDedup\tDuplicates_Removed\tPercent_Dup_Removed\tUnique_Fragments_Final\tSeq_Depth_PostDedup_M\tCondition\tTarget" > ${COMP_DIR}/qc_summary_postdedup.txt

# Process each sample
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
    # Remove any trailing spaces
    SAMPLE=$(echo "$SAMPLE" | tr -d '[:space:]')
    
    echo "=== Processing sample: '$SAMPLE' ===" >> ${LOG_DIR}/debug.log
    
    # Input and output files
    INPUT_BAM=${BAM_DIR}/${SAMPLE}_sorted.bam
    OUTPUT_BAM=${OUT_DIR}/${SAMPLE}_dedup.bam
    METRICS_FILE=${METRICS_DIR}/${SAMPLE}.dedup_metrics.txt
    
    echo "Looking for BAM file: $INPUT_BAM" >> ${LOG_DIR}/debug.log
    
    # Check if input BAM exists
    if [ ! -f "$INPUT_BAM" ]; then
        echo "ERROR: Input BAM file not found for sample $SAMPLE" >> ${LOG_DIR}/debug.log
        echo "Attempted path: $INPUT_BAM" >> ${LOG_DIR}/debug.log
        echo "SAMPLE SKIPPED!" >> ${LOG_DIR}/debug.log
        continue
    fi
    
    echo "BAM file found successfully" >> ${LOG_DIR}/debug.log
    
    # Check if output already exists
    if [ -f "$OUTPUT_BAM" ]; then
        echo "Output already exists for $SAMPLE, skipping deduplication..." >> ${LOG_DIR}/debug.log
    else
        echo "Removing duplicates from $SAMPLE..."
        echo "Removing duplicates from $SAMPLE..." >> ${LOG_DIR}/debug.log
        
        # Remove duplicates using Picard
        picard MarkDuplicates \
            I=$INPUT_BAM \
            O=$OUTPUT_BAM \
            M=$METRICS_FILE \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR=/tmp
        
        # Check if deduplication was successful
        if [ $? -ne 0 ]; then
            echo "Error: Deduplication failed for $SAMPLE" >> ${LOG_DIR}/debug.log
            continue
        fi
        
        # Index the deduplicated BAM file
        echo "Indexing BAM file for $SAMPLE..." >> ${LOG_DIR}/debug.log
        samtools index $OUTPUT_BAM
    fi
    
    # Generate statistics for the deduplicated BAM
    echo "Generating statistics for $SAMPLE deduplicated BAM..." >> ${LOG_DIR}/debug.log
    samtools flagstat $OUTPUT_BAM > ${LOG_DIR}/${SAMPLE}_dedup.flagstat.txt
    
    # Collect post-dedup metrics
    TOTAL_POST=$(samtools view -c $OUTPUT_BAM)
    MAPPED_POST=$(samtools view -c -F 4 $OUTPUT_BAM)
    SEQ_DEPTH_POST=$(echo "scale=2; $TOTAL_POST/1000000" | bc)
    UNIQUE_FINAL=$(samtools view -c -f 2 $OUTPUT_BAM | awk '{print int($1/2)}')
    
    # Parse metrics file for duplicate information
    if [ -f "$METRICS_FILE" ]; then
        DUPS_REMOVED=$(grep -A 1 "LIBRARY" $METRICS_FILE | tail -1 | cut -f 8)
        PERCENT_DUP=$(grep -A 1 "LIBRARY" $METRICS_FILE | tail -1 | cut -f 9 | awk '{printf "%.1f", $1*100}')
    else
        DUPS_REMOVED="NA"
        PERCENT_DUP="NA"
    fi
    
    # Determine metadata from sample ID
    SAMPLE_ID=$(echo $SAMPLE | sed 's/R85_0*//')
    SAMPLE_NUM=$(echo $SAMPLE_ID | bc)
    
    # Assign Target
    if [ $SAMPLE_NUM -ge 1 ] && [ $SAMPLE_NUM -le 12 ]; then
        TARGET="H3K4me3"
    elif [ $SAMPLE_NUM -ge 13 ] && [ $SAMPLE_NUM -le 24 ]; then
        TARGET="KDM5A"
    elif [ $SAMPLE_NUM -ge 25 ] && [ $SAMPLE_NUM -le 36 ]; then
        TARGET="KDM5B"
    elif [ $SAMPLE_NUM -ge 37 ] && [ $SAMPLE_NUM -le 48 ]; then
        TARGET="IgG"
    else
        TARGET="Unknown"
    fi
    
    # Determine Condition
    POSITION=$(( ($SAMPLE_NUM - 1) % 6 + 1 ))
    case $POSITION in
        1) CONDITION="MES" ;;
        2) CONDITION="ADRN" ;;
        3) CONDITION="Untreated" ;;
        4) CONDITION="D7_Cis" ;;
        5) CONDITION="D14_Cis" ;;
        6) CONDITION="D14_C70" ;;
        *) CONDITION="Unknown" ;;
    esac
    
    # Add to post-dedup summary
    echo -e "${SAMPLE}\t${TOTAL_POST}\t${MAPPED_POST}\t${DUPS_REMOVED}\t${PERCENT_DUP}\t${UNIQUE_FINAL}\t${SEQ_DEPTH_POST}\t${CONDITION}\t${TARGET}" >> ${COMP_DIR}/qc_summary_postdedup.txt
    
    echo "Completed processing $SAMPLE" >> ${LOG_DIR}/debug.log
    echo "Completed processing $SAMPLE"
    
done < "$TMP_SAMPLE_LIST"

echo "Deduplication Complete"
echo "Deduplication Complete" >> ${LOG_DIR}/debug.log

# Clean up temporary file
rm -f "$TMP_SAMPLE_LIST"

# Copy original QC summary if it exists
if [ -f "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/QC/qc_summary.txt" ]; then
    cp /data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/QC/qc_summary.txt ${COMP_DIR}/qc_summary_original.txt
    echo "Original QC summary copied to comparison directory" >> ${LOG_DIR}/debug.log
else
    echo "WARNING: Original QC summary not found at /data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/QC/qc_summary.txt" >> ${LOG_DIR}/debug.log
fi

# Generate comparison plots
echo "Generating comparison plots..." >> ${LOG_DIR}/debug.log

# Create R script for comparison plots
cat > ${COMP_DIR}/plot_dedup_comparison.R << 'EOF'
#!/usr/bin/env Rscript

# Force R to use only conda environment libraries
.libPaths(.libPaths()[1])
cat("Using library path:", .libPaths()[1], "\n")

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(scales)

cat("All packages loaded successfully!\n")

# Read data with error checking
if (!file.exists("qc_summary_original.txt")) {
    stop("Original QC summary file not found: qc_summary_original.txt")
}

if (!file.exists("qc_summary_postdedup.txt")) {
    stop("Post-deduplication QC summary file not found: qc_summary_postdedup.txt")
}

data_orig <- read.table("qc_summary_original.txt", header=TRUE, sep="\t")
data_post <- read.table("qc_summary_postdedup.txt", header=TRUE, sep="\t")

# Check if data frames are empty
if (nrow(data_orig) == 0) {
    stop("Original QC summary file is empty")
}

if (nrow(data_post) == 0) {
    stop("Post-deduplication QC summary file is empty")
}

cat("Original data rows:", nrow(data_orig), "\n")
cat("Post-dedup data rows:", nrow(data_post), "\n")

# Add stage identifier
data_orig$Stage <- "Before Dedup"
data_post$Stage <- "After Dedup"

# Check column names and prepare post-dedup data
cat("Original data columns:", colnames(data_orig), "\n")
cat("Post-dedup data columns:", colnames(data_post), "\n")

# Prepare post-dedup data to match original format
if ("Percent_Dup_Removed" %in% colnames(data_post)) {
    data_post$Duplication_Rate <- as.numeric(data_post$Percent_Dup_Removed)
} else {
    cat("Warning: Percent_Dup_Removed column not found in post-dedup data\n")
    data_post$Duplication_Rate <- 0  # Default value
}

if ("Seq_Depth_PostDedup_M" %in% colnames(data_post)) {
    data_post$Seq_Depth_M <- data_post$Seq_Depth_PostDedup_M
} else {
    cat("Warning: Seq_Depth_PostDedup_M column not found in post-dedup data\n")
    data_post$Seq_Depth_M <- data_post$Seq_Depth_M  # Use existing column
}

if ("Unique_Fragments_Final" %in% colnames(data_post)) {
    data_post$Unique_Fragments <- data_post$Unique_Fragments_Final
} else {
    cat("Warning: Unique_Fragments_Final column not found in post-dedup data\n")
    data_post$Unique_Fragments <- data_post$Unique_Fragments  # Use existing column
}

# Combine relevant columns with error checking
tryCatch({
    data_combined <- bind_rows(
        data_orig %>% select(Sample, Condition, Target, Duplication_Rate, Seq_Depth_M, Unique_Fragments, Stage),
        data_post %>% select(Sample, Condition, Target, Duplication_Rate, Seq_Depth_M, Unique_Fragments, Stage)
    )
    cat("Data combination successful. Combined rows:", nrow(data_combined), "\n")
}, error = function(e) {
    cat("Error combining data:", e$message, "\n")
    cat("Trying alternative approach...\n")
    
    # Alternative: use only original data if combination fails
    data_combined <- data_orig %>% 
        select(Sample, Condition, Target, Duplication_Rate, Seq_Depth_M, Unique_Fragments, Stage)
    cat("Using original data only. Rows:", nrow(data_combined), "\n")
})

# Set factor orders
condition_order <- c("ADRN", "D14_C70", "D14_Cis", "D7_Cis", "MES", "Untreated")
data_combined$Condition <- factor(data_combined$Condition, levels=condition_order)
data_combined$Target <- factor(data_combined$Target, levels=c("H3K4me3", "IgG", "KDM5A", "KDM5B"))
data_combined$Stage <- factor(data_combined$Stage, levels=c("Before Dedup", "After Dedup"))

# Colors
target_colors <- c("H3K4me3"="#E41A1C", "IgG"="#377EB8", "KDM5A"="#4DAF4A", "KDM5B"="#984EA3")

# Create comparison data for arrows
comparison_data <- data_orig %>%
    select(Sample, Condition, Target, 
           Depth_Before = Seq_Depth_M, 
           Unique_Before = Unique_Fragments) %>%
    left_join(
        data_post %>% select(Sample, 
                           Depth_After = Seq_Depth_M, 
                           Unique_After = Unique_Fragments),
        by = "Sample"
    ) %>%
    mutate(
        Depth_Reduction = Depth_Before - Depth_After,
        Depth_Reduction_Pct = (Depth_Reduction / Depth_Before) * 100,
        Unique_Change = Unique_After - Unique_Before,
        Unique_Change_Pct = (Unique_Change / Unique_Before) * 100
    )

# Plot 1: Depth Reduction by Target and Condition
p1 <- ggplot(comparison_data, aes(x=interaction(Condition, Target, sep="."), 
                                  y=Depth_Reduction_Pct, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +
    geom_hline(yintercept=0, linetype="solid", color="gray50") +
    labs(title="A. Sequencing Depth Reduction (%) After Deduplication - R85 CUT&Tag",
         y="Depth Reduction (%)",
         x="Condition.Target") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          axis.text.y=element_text(size=9),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="right",
          legend.text=element_text(size=9),
          legend.title=element_text(size=10),
          panel.grid.minor=element_blank())

# Plot 2: Sequencing Depth After Deduplication (like script 2)
p2 <- ggplot(comparison_data, aes(x=interaction(Condition, Target, sep="."), 
                                  y=Depth_After/1000000, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +
    labs(title="B. Sequencing Depth After Deduplication - R85 CUT&Tag",
         y="Sequencing Depth per Million",
         x="Condition.Target") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          axis.text.y=element_text(size=9),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="right",
          legend.text=element_text(size=9),
          legend.title=element_text(size=10),
          panel.grid.minor=element_blank())

# Save individual plots
ggsave("Depth_Reduction_After_Dedup.pdf", p1, width=10, height=6, dpi=300)
ggsave("Sequencing_Depth_After_Dedup.pdf", p2, width=10, height=6, dpi=300)

# Create combined figure with both plots
library(gridExtra)
pdf("Deduplication_Analysis_Combined.pdf", width=12, height=10)
grid.arrange(p1, p2, ncol=1)
dev.off()

# Generate summary statistics for depth reduction
cat("\n=== DEPTH REDUCTION SUMMARY ===\n")

# Summary by target
summary_by_target <- comparison_data %>%
    group_by(Target) %>%
    summarise(
        Mean_Depth_Before = mean(Depth_Before, na.rm=TRUE),
        Mean_Depth_After = mean(Depth_After, na.rm=TRUE),
        Mean_Depth_Reduction_Pct = mean(Depth_Reduction_Pct, na.rm=TRUE),
        Median_Depth_Reduction_Pct = median(Depth_Reduction_Pct, na.rm=TRUE),
        Min_Depth_Reduction_Pct = min(Depth_Reduction_Pct, na.rm=TRUE),
        Max_Depth_Reduction_Pct = max(Depth_Reduction_Pct, na.rm=TRUE)
    ) %>%
    mutate(across(where(is.numeric), ~round(., 2)))

print(summary_by_target)

# Summary by condition
summary_by_condition <- comparison_data %>%
    group_by(Condition) %>%
    summarise(
        Mean_Depth_Reduction_Pct = mean(Depth_Reduction_Pct, na.rm=TRUE),
        Median_Depth_Reduction_Pct = median(Depth_Reduction_Pct, na.rm=TRUE)
    ) %>%
    mutate(across(where(is.numeric), ~round(., 2)))

cat("\n=== DEPTH REDUCTION BY CONDITION ===\n")
print(summary_by_condition)

# Save summary tables
write.table(summary_by_target, "depth_reduction_by_target.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(summary_by_condition, "depth_reduction_by_condition.txt", sep="\t", row.names=FALSE, quote=FALSE)

cat("\n=== Analysis Complete ===\n")
cat("Results saved in current directory\n")
cat("Individual plots:\n")
cat("  - Depth_Reduction_After_Dedup.pdf\n")
cat("  - Sequencing_Depth_After_Dedup.pdf\n")
cat("Combined plot: Deduplication_Analysis_Combined.pdf\n")
EOF

# Run R script with proper library path
cd $COMP_DIR

# Set R library path to use conda environment packages
export R_LIBS_USER=""
export R_LIBS_SITE=""

# Temporarily rename problematic personal R library to avoid conflicts
if [ -d "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3" ]; then
    echo "Temporarily renaming problematic R library..."
    mv "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3" "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3_backup"
fi

# Run R script
Rscript plot_dedup_comparison.R

# Restore the personal R library
if [ -d "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3_backup" ]; then
    echo "Restoring personal R library..."
    mv "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3_backup" "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3"
fi

echo "=== Deduplication and Comparison Complete ===" >> ${LOG_DIR}/debug.log
date >> ${LOG_DIR}/debug.log

echo "Deduplication and comparison analysis complete!"
echo "Results saved in: $COMP_DIR"
echo "Main report: ${COMP_DIR}/Dedup_Comparison_Analysis.pdf"