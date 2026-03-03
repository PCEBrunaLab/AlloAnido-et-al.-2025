#!/bin/bash
## =============================================================================
## Neuroblastoma Patient scRNA-seq - inferCNV Analysis
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
## =============================================================================
## CNV detection using inferCNV with Cluster 10 as reference
## Cluster 10: immune/stromal cells (diploid reference)
## =============================================================================

#SBATCH --job-name=NB_inferCNV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
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
PROJECT_ROOT="/path/to/Patient_Samples"
RESULTS_DIR="${BASE_DIR}/results"
REFERENCE_DIR="${BASE_DIR}/references"
LOG_DIR="${BASE_DIR}/logs"
SRC_DIR="${BASE_DIR}/src"

# Create directories if they don't exist
mkdir -p "${LOG_DIR}" "${RESULTS_DIR}/inferCNV_cluster10"

# Input Seurat object (prefer Jansky-annotated version if available)
SEURAT_RDS="${RESULTS_DIR}/data/NB_patients_seurat_with_jansky.rds"

# Gene annotation file
GTF_FILE="${REFERENCE_DIR}/gencode.v43.annotation.gtf.gz"

# Output directory
OUTPUT_DIR="${RESULTS_DIR}/inferCNV_cluster10"

# R script for inferCNV analysis
R_SCRIPT="${SRC_DIR}/inferCNV_cluster10_reference.R"

## =============================================================================
## ACTIVATE CONDA ENVIRONMENT
## =============================================================================

module use /opt/software/easybuild/modules/all/
module load Mamba
module load OpenBLAS/0.3.20-GCC-11.3.0

source ~/.bashrc
conda activate infercnv_env

## =============================================================================
## SET LIBRARY PATHS
## =============================================================================

# Set LD_LIBRARY_PATH to include OpenBLAS so system R libraries can find it
# This is needed because system-installed R packages (like ape) were compiled against OpenBLAS
if [[ -n "${EBROOTOPENBLAS}" ]]; then
    export LD_LIBRARY_PATH="${EBROOTOPENBLAS}/lib:${LD_LIBRARY_PATH}"
    echo "OpenBLAS module loaded: ${EBROOTOPENBLAS}"
    echo "LD_LIBRARY_PATH includes: ${EBROOTOPENBLAS}/lib"
fi

# Also add conda library path if available
if [[ -n "${CONDA_PREFIX}" ]]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"
    echo "CONDA_PREFIX: ${CONDA_PREFIX}"
fi

echo "Final LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"

## =============================================================================
## CHECK INPUTS
## =============================================================================

echo "====================================================="
echo "Running inferCNV for NB patients"
echo "====================================================="
echo "Seurat object: ${SEURAT_RDS}"
echo "GTF annotation: ${GTF_FILE}"
echo "Output dir:    ${OUTPUT_DIR}"
echo "Start time:    $(date)"
echo "====================================================="

if [[ ! -f "${SEURAT_RDS}" ]]; then
    echo "ERROR: Seurat RDS not found at ${SEURAT_RDS}" >&2
    exit 1
fi

if [[ ! -f "${GTF_FILE}" ]]; then
    echo "ERROR: GTF file not found at ${GTF_FILE}" >&2
    exit 1
fi

echo "✓ Input files verified"

## =============================================================================
## RUN INFERCNV
## =============================================================================

cd "${BASE_DIR}"

THREADS=${SLURM_CPUS_PER_TASK:-8}

srun Rscript "${R_SCRIPT}" \
    --seurat_rds "${SEURAT_RDS}" \
    --cluster_column seurat_clusters \
    --reference_clusters 10 \
    --gtf "${GTF_FILE}" \
    --output_dir "${OUTPUT_DIR}" \
    --threads "${THREADS}"

## =============================================================================
## CHECK COMPLETION
## =============================================================================

if [ $? -eq 0 ]; then
    echo ""
    echo "====================================================="
    echo "inferCNV COMPLETED SUCCESSFULLY"
    echo "====================================================="
    echo "Output directory: ${OUTPUT_DIR}"
    echo "End time: $(date)"
    echo "====================================================="
else
    echo ""
    echo "====================================================="
    echo "ERROR: inferCNV FAILED"
    echo "====================================================="
    echo "Check logs: ${LOG_DIR}/"
    echo "====================================================="
    exit 1
fi
