#!/bin/bash
## =============================================================================
## Neuroblastoma Patient scRNA-seq - inferCNV Visualization
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
## =============================================================================
## Creates UMAP visualizations of inferCNV results
## Projects CNV scores onto cell embeddings for interpretation
## =============================================================================

#SBATCH --job-name=NB_viz_inferCNV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
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
LOG_DIR="${BASE_DIR}/logs"
SRC_DIR="${BASE_DIR}/src"
R_SCRIPT="${SRC_DIR}/visualize_inferCNV_results.R"

# Create directories if they don't exist
mkdir -p "${LOG_DIR}"

## =============================================================================
## ACTIVATE CONDA ENVIRONMENT
## =============================================================================

module use /opt/software/easybuild/modules/all/
module load Mamba

source ~/.bashrc
conda activate infercnv_env

## =============================================================================
## RUN VISUALIZATION
## =============================================================================

cd "${BASE_DIR}"

echo "============================================================"
echo "inferCNV UMAP Visualization"
echo "============================================================"
echo "Started: $(date)"
echo "R script: ${R_SCRIPT}"
echo "============================================================"

srun Rscript "${R_SCRIPT}"

## =============================================================================
## CHECK COMPLETION
## =============================================================================

if [ $? -eq 0 ]; then
    echo ""
    echo "============================================================"
    echo "Visualization COMPLETED SUCCESSFULLY"
    echo "============================================================"
    echo "Completed: $(date)"
    echo "============================================================"
else
    echo ""
    echo "============================================================"
    echo "ERROR: Visualization FAILED"
    echo "============================================================"
    echo "Check logs: ${LOG_DIR}/"
    echo "============================================================"
    exit 1
fi
