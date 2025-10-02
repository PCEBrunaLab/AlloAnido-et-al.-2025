#!/bin/bash
#SBATCH --job-name=Mapping_CUT_Tag_R85
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --array=1-48%8
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%A_%a.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%A_%a.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda activate bowtie2.4.2
module load SAMtools/1.11

# Define directories and variables
GENOME_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Genome
BOWTIE2_INDEX=${GENOME_DIR}/bowtie2_index/GRCh38

# Directory containing fastq files
FASTQ_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/R85_run1883

# Output directories
OUT_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/bam_files
LOG_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs
QC_DIR=${OUT_DIR}/qc

# Number of threads to use
PPN=$SLURM_CPUS_PER_TASK

# Create output directories
mkdir -p $OUT_DIR
mkdir -p $LOG_DIR
mkdir -p $QC_DIR

# Get sample list
SAMPLE_LIST=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/sample_list.txt

# Get the specific sample for this array job
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)

# Check if we got a valid sample
if [ -z "$SAMPLE" ]; then
    echo "Error: No sample found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Define input files - using wildcard to match Illumina naming
R1=$(ls ${FASTQ_DIR}/${SAMPLE}_S*_L*_R1_001.fastq.gz 2>/dev/null | head -n1)
R2=$(ls ${FASTQ_DIR}/${SAMPLE}_S*_L*_R2_001.fastq.gz 2>/dev/null | head -n1)

# Check if input files exist
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Error: Input files not found for sample $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2"
    exit 1
fi

# Check if output already exists
if [ -f "${OUT_DIR}/${SAMPLE}_sorted.bam" ]; then
    echo "[$(date)] Skipping $SAMPLE - output already exists"
    exit 0
fi

echo "Processing sample: $SAMPLE (Array job $SLURM_ARRAY_TASK_ID)"

# Mapping with CUT&Tag specific parameters
echo "[$(date)] Running bowtie2 for $SAMPLE..."
bowtie2 --local  \
    --very-sensitive \
    --no-unal \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    -p $PPN \
    -x $BOWTIE2_INDEX \
    -1 $R1 \
    -2 $R2 \
    2> ${LOG_DIR}/${SAMPLE}.bowtie2.log | \
    samtools sort -@ $PPN -O bam -o ${OUT_DIR}/${SAMPLE}_sorted.bam

# Check if mapping was successful
if [ $? -ne 0 ]; then
    echo "Error: bowtie2 mapping failed for $SAMPLE"
    exit 1
fi

# Index BAM
echo "[$(date)] Indexing BAM file for $SAMPLE..."
samtools index -@ $PPN ${OUT_DIR}/${SAMPLE}_sorted.bam

# Generate mapping statistics
echo "[$(date)] Generating mapping statistics for $SAMPLE..."
samtools flagstat ${OUT_DIR}/${SAMPLE}_sorted.bam > ${QC_DIR}/${SAMPLE}.flagstat.txt

# Log completion
echo "[$(date)] Completed mapping for sample: $SAMPLE"