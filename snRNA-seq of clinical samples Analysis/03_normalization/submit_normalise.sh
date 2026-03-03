#!/bin/bash
## =============================================================================
## Neuroblastoma Patient scRNA-seq - Normalization
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
## =============================================================================
## Step 3 of scRNA-seq pipeline: Multi-batch normalization
## Reads all 18 NB patient SCE files (post-emptyDrops), combines them,
## and performs multi-batch normalization.
## Output: NB_patients_SCE-norm.RDS
## =============================================================================

#SBATCH --job-name=NB_normalisation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=USER_EMAIL
#SBATCH --mail-type=ALL

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

BASE_DIR="/path/to/NB_patient_analysis"
SCE_DIR="${BASE_DIR}/data/SCE_objects/relaxed_filters"
SRC_DIR="${BASE_DIR}/src"
OUTPUT_DIR="${BASE_DIR}/results"
LOGS_DIR="${BASE_DIR}/logs"

# Create output directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOGS_DIR}

# Output file
OUTPUT="${OUTPUT_DIR}/NB_patients_SCE-norm_relaxed.RDS"

## =============================================================================
## ACTIVATE CONDA ENVIRONMENT
## =============================================================================

module use /opt/software/easybuild/modules/all/
module load Mamba

source ~/.bashrc
conda activate hotspot

## =============================================================================
## RUN NORMALIZATION
## =============================================================================

cd ${BASE_DIR}

echo "====================================================="
echo "Starting Normalization for NB Patient Samples"
echo "====================================================="
echo "SCE directory:  $SCE_DIR"
echo "Output:         $OUTPUT"
echo "Cores:          ${SLURM_CPUS_PER_TASK}"
echo "Start time:     $(date)"
echo "====================================================="

srun Rscript ${SRC_DIR}/norm_counts.R \
    --scedir=${SCE_DIR} \
    --breaks \
    --output=${OUTPUT}

## =============================================================================
## CHECK COMPLETION
## =============================================================================

if [ $? -eq 0 ] && [ -f "$OUTPUT" ]; then
    echo ""
    echo "====================================================="
    echo "Normalization COMPLETED SUCCESSFULLY"
    echo "====================================================="
    echo "Output file: $OUTPUT"
    echo "File size:   $(du -h ${OUTPUT} | cut -f1)"
    echo "End time:    $(date)"
    echo "====================================================="
else
    echo ""
    echo "====================================================="
    echo "ERROR: Normalization FAILED"
    echo "====================================================="
    echo "Check logs: ${LOGS_DIR}/"
    echo "====================================================="
    exit 1
fi
