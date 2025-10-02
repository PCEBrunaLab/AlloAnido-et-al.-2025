#!/bin/bash
#SBATCH --job-name=CUTTag_QC_R85
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
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

# Define directories - Following R85 project structure
BAM_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/bam_files
DEDUP_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/bam_files_dedup
QC_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/QC_final
LOG_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/6_QC_Final

# Create output directories
mkdir -p $QC_DIR
mkdir -p $LOG_DIR

# Get sample list
SAMPLE_LIST=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/sample_list.txt

# Clean up the sample list
TMP_SAMPLE_LIST=${LOG_DIR}/clean_sample_list.txt
tr -d '\r' < "$SAMPLE_LIST" | sed 's/[[:space:]]*$//' > "$TMP_SAMPLE_LIST"

# Start QC log
echo "=== Starting QC Analysis ===" > ${LOG_DIR}/qc_analysis.log
date >> ${LOG_DIR}/qc_analysis.log

# Create summary header
echo -e "Sample\tTotal_Reads\tMapped_Reads\tAlignment_Rate\tDuplication_Rate\tUnique_Fragments\tSeq_Depth_M\tMapped_Frags_M\tCondition\tTarget" > ${QC_DIR}/qc_summary.txt

# Process each sample
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
    # Remove any trailing spaces
    SAMPLE=$(echo "$SAMPLE" | tr -d '[:space:]')
    
    echo "Processing QC for $SAMPLE..." >> ${LOG_DIR}/qc_analysis.log
    
    # Define BAM files
    ORIGINAL_BAM=${BAM_DIR}/${SAMPLE}_sorted.bam
    DEDUP_BAM=${DEDUP_DIR}/${SAMPLE}_dedup.bam
    
    # Check if files exist
    if [ ! -f "$ORIGINAL_BAM" ]; then
        echo "Original BAM not found for $SAMPLE, skipping..." >> ${LOG_DIR}/qc_analysis.log
        continue
    fi
    
    # Calculate metrics from ORIGINAL BAM
    TOTAL_READS=$(samtools view -c $ORIGINAL_BAM)
    MAPPED_READS=$(samtools view -c -F 4 $ORIGINAL_BAM)
    ALIGNMENT_RATE=$(echo "scale=1; ($MAPPED_READS * 100)/$TOTAL_READS" | bc)
    SEQ_DEPTH=$(echo "scale=2; $TOTAL_READS/1000000" | bc)
    MAPPED_FRAGS=$(echo "scale=2; ($MAPPED_READS/2)/1000000" | bc)
    
    # Run Picard MarkDuplicates to get duplication metrics (without removing)
    echo "Calculating duplication metrics for $SAMPLE..."
    picard MarkDuplicates \
        I=$ORIGINAL_BAM \
        O=${QC_DIR}/${SAMPLE}_marked.bam \
        M=${QC_DIR}/${SAMPLE}_dup_metrics.txt \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=/tmp \
        QUIET=true
    
    # Extract duplication rate from Picard metrics
    if [ -f "${QC_DIR}/${SAMPLE}_dup_metrics.txt" ]; then
        DUP_RATE=$(grep -A 1 "LIBRARY" ${QC_DIR}/${SAMPLE}_dup_metrics.txt | tail -1 | cut -f9 | awk '{printf "%.1f", $1*100}')
        if [ -z "$DUP_RATE" ]; then
            DUP_RATE="0"
        fi
    else
        DUP_RATE="0"
    fi
    
    # Count unique fragments (non-duplicate paired reads / 2)
    if [ -f "${QC_DIR}/${SAMPLE}_marked.bam" ]; then
        UNIQUE_FRAGS=$(samtools view -c -f 2 -F 1024 ${QC_DIR}/${SAMPLE}_marked.bam | awk '{print int($1/2)}')
    else
        UNIQUE_FRAGS="0"
    fi
    
    # Determine metadata from sample ID
    SAMPLE_ID=$(echo $SAMPLE | sed 's/R85_0*//')  # Remove R85_ and leading zeros
    SAMPLE_NUM=$(echo $SAMPLE_ID | bc)  # Convert to number
    
    # Assign Target based on ID ranges
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
    
    # Determine Condition based on position within each group of 6
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
    
    # Add to summary
    echo -e "${SAMPLE}\t${TOTAL_READS}\t${MAPPED_READS}\t${ALIGNMENT_RATE}\t${DUP_RATE}\t${UNIQUE_FRAGS}\t${SEQ_DEPTH}\t${MAPPED_FRAGS}\t${CONDITION}\t${TARGET}" >> ${QC_DIR}/qc_summary.txt
    
    # Clean up temporary marked BAM
    rm -f ${QC_DIR}/${SAMPLE}_marked.bam
    rm -f ${QC_DIR}/${SAMPLE}_marked.bam.bai
    
    echo "Completed QC for $SAMPLE" >> ${LOG_DIR}/qc_analysis.log
    
done < "$TMP_SAMPLE_LIST"

# Create R script for plotting
cat > ${QC_DIR}/plot_qc.R << 'EOF'
#!/usr/bin/env Rscript

# Force R to use only conda environment libraries
.libPaths(.libPaths()[1])  # Keep only the first (conda) library path
cat("Using library path:", .libPaths()[1], "\n")

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)

cat("All packages loaded successfully from conda environment!\n")

# Read data
data <- read.table("qc_summary.txt", header=TRUE, sep="\t")

# Set factor orders - ALPHABETICAL for consistency
condition_order <- c("ADRN", "D14_C70", "D14_Cis", "D7_Cis", "MES", "Untreated")
data$Condition <- factor(data$Condition, levels=condition_order)
data$Target <- factor(data$Target, levels=c("H3K4me3", "IgG", "KDM5A", "KDM5B"))

# Colors - match your PDF exactly
target_colors <- c("H3K4me3"="#E41A1C", "IgG"="#377EB8", "KDM5A"="#4DAF4A", "KDM5B"="#984EA3")

# Create interaction factor for x-axis grouping
data$Condition.Mark <- interaction(data$Condition, data$Target, sep=".")

# Plot with exact formatting
p1 <- ggplot(data, aes(x=Condition.Mark, y=Seq_Depth_M, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +  # Changed legend title
    labs(title="A. Sequencing Depth - R85 CUT&Tag", 
         y="Sequencing Depth per Million", 
         x="Condition.Mark") +
    theme_bw(base_size=10) +  # Smaller base size
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),  # Smaller x labels
          axis.text.y=element_text(size=9),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),  # Centered title
          legend.position="right",
          legend.text=element_text(size=9),
          legend.title=element_text(size=10),
          panel.grid.minor=element_blank())  # Remove minor gridlines

# Plot 2: Duplication Rate (NEW - same style)
p2 <- ggplot(data, aes(x=Condition.Mark, y=Duplication_Rate, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +
    geom_hline(yintercept=40, linetype="dashed", color="orange", alpha=0.5) +  # Warning line
    geom_hline(yintercept=60, linetype="dashed", color="red", alpha=0.5) +     # Critical line
    labs(title="B. Duplication Rate - R85 CUT&Tag", 
         y="Duplication Rate (%)", 
         x="Condition.Mark") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          axis.text.y=element_text(size=9),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="right",
          legend.text=element_text(size=9),
          legend.title=element_text(size=10),
          panel.grid.minor=element_blank())

# Plot 3: Unique Library Size (same style)
p3 <- ggplot(data, aes(x=Condition.Mark, y=Unique_Fragments/1000, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +
    labs(title="C. Unique Library Size - R85 CUT&Tag", 
         y="Unique Fragments (Thousands)", 
         x="Condition.Mark") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          axis.text.y=element_text(size=9),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="right",
          legend.text=element_text(size=9),
          legend.title=element_text(size=10),
          panel.grid.minor=element_blank())

# Save individual plots with consistent dimensions
ggsave("QC_01_Sequencing_Depth.pdf", p1, width=10, height=5, dpi=300)
ggsave("QC_02_Duplication_Rate.pdf", p2, width=10, height=5, dpi=300)
ggsave("QC_03_Unique_Library.pdf", p3, width=10, height=5, dpi=300)

# Create combined figure with all three plots
library(gridExtra)
pdf("CUT_Tag_QC_Combined.pdf", width=12, height=12)
grid.arrange(p1, p2, p3, ncol=1)
dev.off()

# Print summary
cat("\n=== QC SUMMARY STATISTICS ===\n")
summary_by_target <- aggregate(cbind(Alignment_Rate, Duplication_Rate, Seq_Depth_M) ~ Target, 
                               data, function(x) round(mean(x), 1))
print(summary_by_target)

# Flag problematic samples
problems <- data[data$Alignment_Rate < 80 | data$Duplication_Rate > 70,]
if(nrow(problems) > 0) {
    cat("\n=== WARNING: Samples with QC issues ===\n")
    print(problems[,c("Sample","Alignment_Rate","Duplication_Rate","Target")])
    write.table(problems, "flagged_samples.txt", sep="\t", row.names=FALSE, quote=FALSE)
}

cat("\nQC plots saved successfully!\n")
EOF

# Run R script with proper library path
cd $QC_DIR

# Set R library path to use conda environment packages
export R_LIBS_USER=""
export R_LIBS_SITE=""

# Temporarily rename problematic personal R library to avoid conflicts
if [ -d "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3" ]; then
    echo "Temporarily renaming problematic R library..."
    mv "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3" "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3_backup"
fi

# Run R script
Rscript plot_qc.R

# Restore the personal R library
if [ -d "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3_backup" ]; then
    echo "Restoring personal R library..."
    mv "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3_backup" "/data/scratch/DMP/DUDMP/PAEDCANC/asadr/R/x86_64-redhat-linux-gnu-library/4.3"
fi

# Clean up
rm -f "$TMP_SAMPLE_LIST"

echo "=== QC Analysis Complete ===" >> ${LOG_DIR}/qc_analysis.log
date >> ${LOG_DIR}/qc_analysis.log

echo "QC analysis complete!"
echo "Results saved in: $QC_DIR"
echo "Main report: ${QC_DIR}/CUT_Tag_QC_Report.pdf"