#!/bin/bash
## =============================================================================
## Neuroblastoma Patient scRNA-seq - EmptyDrops Summary Report
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
## =============================================================================
## Generates a summary table of cell calling results
## Run this AFTER all EmptyDrops jobs have completed
## =============================================================================

#SBATCH --job-name=emptyDrops_report_NB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --partition=compute
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=USER_EMAIL
#SBATCH --mail-type=ALL

set -e  # Exit on error

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

BASE_DIR="/path/to/NB_patient_analysis"
CELLRANGER_DIR="${BASE_DIR}/cellranger_output"
DATA_DIR="${BASE_DIR}/data/SCE_objects/relaxed_filters"
REPORT_DIR="${BASE_DIR}/reports"
SCRIPT_DIR="${BASE_DIR}/src"

# Create report directory
mkdir -p ${REPORT_DIR}

OUTPUT_FILE="${REPORT_DIR}/emptyDrops_summary_NB_relaxed_filters.tsv"

## =============================================================================
## ACTIVATE CONDA ENVIRONMENT
## =============================================================================

module use /opt/software/easybuild/modules/all/
module load OpenBLAS/0.3.20-GCC-11.3.0
module load Mamba

source ~/.bashrc
conda activate hotspot

## =============================================================================
## GENERATE REPORT
## =============================================================================

echo "================================================================================"
echo "EmptyDrops Summary Report - NB Samples"
echo "================================================================================"
echo "Data directory:         ${DATA_DIR}"
echo "Cell Ranger directory:  ${CELLRANGER_DIR}"
echo "Output file:            ${OUTPUT_FILE}"
echo "Start time:             $(date)"
echo "================================================================================"

Rscript ${SCRIPT_DIR}/emptyDrops_report.R \
    --datadir=${DATA_DIR} \
    --cellrangerdir=${CELLRANGER_DIR} \
    --output=${OUTPUT_FILE} \
    --prefix=NB

## =============================================================================
## CHECK COMPLETION
## =============================================================================

if [ $? -eq 0 ] && [ -f "$OUTPUT_FILE" ]; then
    echo ""
    echo "================================================================================"
    echo "REPORT GENERATION COMPLETED SUCCESSFULLY"
    echo "================================================================================"
    echo "Output file: ${OUTPUT_FILE}"
    echo "File size:   $(du -h ${OUTPUT_FILE} | cut -f1)"
    echo "End time:    $(date)"
    echo "================================================================================"
    echo ""
    echo "View the report with:"
    echo "  cat ${OUTPUT_FILE}"
    echo "  or"
    echo "  column -t -s $'\\t' ${OUTPUT_FILE} | less -S"
    echo "================================================================================"
else
    echo ""
    echo "================================================================================"
    echo "ERROR: Report generation FAILED"
    echo "================================================================================"
    exit 1
fi
