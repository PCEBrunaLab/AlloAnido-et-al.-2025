#!/bin/bash
#SBATCH --job-name=TrimGalore_Memoryseq
#SBATCH --partition=compute
#SBATCH --array=1-74%3
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL

module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda activate trim-galore0.6.6

FASTQ_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Data_RNAseq/00_fastq
OUT_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Data_RNAseq/trimmed_fastq
SAMPLE_LIST=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/sample_list.txt

mkdir -p $OUT_DIR

# Read sample from the list based on array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)

if [ ! -f "${OUT_DIR}/${SAMPLE}_R1_001_val_1.fq.gz" ] || [ ! -f "${OUT_DIR}/${SAMPLE}_R2_001_val_2.fq.gz" ]; then
   echo "Processing sample: $SAMPLE"
   trim_galore --paired --illumina --cores 8 \
       -o $OUT_DIR \
       ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
       ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz
else
   echo "Skipping $SAMPLE - already processed"
fi