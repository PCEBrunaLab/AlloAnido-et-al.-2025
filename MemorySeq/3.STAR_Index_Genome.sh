#!/bin/bash
#SBATCH --job-name=STAR_Index
#SBATCH --partition=smp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda activate star2.7.6a

# Define directories
GENOME_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Genome
INDEX_DIR=${GENOME_DIR}/STAR_index

# Create index directory
mkdir -p $INDEX_DIR

# Run STAR indexing
STAR --runMode genomeGenerate \
     --runThreadN 24 \
     --genomeDir $INDEX_DIR \
     --genomeFastaFiles ${GENOME_DIR}/Homo_sapiens.GRCh38.dna.toplevel.fa \
     --sjdbGTFfile ${GENOME_DIR}/Homo_sapiens.GRCh38.113.gtf \
     --sjdbOverhang 149