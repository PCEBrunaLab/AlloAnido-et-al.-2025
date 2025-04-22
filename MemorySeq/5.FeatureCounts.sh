#!/bin/bash
#SBATCH --job-name=featureCounts_index
#SBATCH --partition=smp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL

# Load modules
module load anaconda/3 SAMtools/1.11
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda activate memoryseq_env

# Define directories
BAM_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/BAM_Files
GTF_FILE=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Genome/Homo_sapiens.GRCh38.113.gtf
OUT_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/MemorySeq/Counts

mkdir -p $OUT_DIR

# Index BAM files
for bam in ${BAM_DIR}/*_Aligned.sortedByCoord.out.bam; do
    samtools index $bam
done

# Run featureCounts
featureCounts \
    -T 24 \
    -p \
    -t exon \
    -g gene_id \
    -a $GTF_FILE \
    -o ${OUT_DIR}/all_samples_counts.txt \
    ${BAM_DIR}/*_Aligned.sortedByCoord.out.bam