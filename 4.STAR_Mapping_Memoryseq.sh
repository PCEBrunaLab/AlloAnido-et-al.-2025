#!/bin/bash
#SBATCH --job-name=STAR_Mapping_Memoryseq
#SBATCH --partition=smp
#SBATCH --array=1-24%3    # Run 3 jobs simultaneously
#SBATCH --cpus-per-task=8 # 8 cores per task
#SBATCH --time=120:00:00  # 5 days
#SBATCH --mem=128G        # 16GB * 8 cores
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Logs/%x_%A_%a.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Logs/%x_%A_%a.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --export=ALL

# Load necessary modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

# Activate STAR environment
conda activate star2.7.6a

# Define directories
TRIM_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Data_RNAseq/trimmed_fastq
OUT_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/BAM_Files
GENOME_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Genome/STAR_index
SAMPLE_LIST=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/sample_list.txt

# Create output directory
mkdir -p $OUT_DIR

# Read samples into an array
mapfile -t SAMPLES < "$SAMPLE_LIST"

# Calculate which samples to process in this array task
start_idx=$(( ($SLURM_ARRAY_TASK_ID-1)*3 ))
TASK_SAMPLES=("${SAMPLES[@]:$start_idx:3}")

# Process samples assigned to this array task
for sample in "${TASK_SAMPLES[@]}"; do
    if [ ! -z "$sample" ]; then
        echo "Processing sample: $sample"
        STAR \
            --runThreadN 8 \
            --genomeDir $GENOME_DIR \
            --readFilesIn ${TRIM_DIR}/${sample}_R1_001_val_1.fq.gz ${TRIM_DIR}/${sample}_R2_001_val_2.fq.gz \
            --readFilesCommand zcat \
            --outFileNamePrefix ${OUT_DIR}/${sample}_ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --quantMode GeneCounts
    fi
done