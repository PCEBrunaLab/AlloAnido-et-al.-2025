#!/bin/bash
#SBATCH --job-name=QC_Memoryseq
#SBATCH --partition=compute
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=8000
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --export=ALL

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

SCRATCH_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq

# Change to the job directory
cd $SCRATCH_DIR

# Changing directories to where the fastq files are located
cd /data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Data/00_fastq

# Running FASTQC
module load FastQC/0.11.9
mkdir -p fastqc
# Find all of the fastq files
find /data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Data/00_fastq -name "*.gz" -exec fastqc -o fastqc/ {} \;

# Run multiqc to compile individual fastqc files, this helps visualization of fastqc reports
conda activate multiqc1.9
mkdir -p multiqc
multiqc /data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Data/00_fastq/fastqc/ -o multiqc
