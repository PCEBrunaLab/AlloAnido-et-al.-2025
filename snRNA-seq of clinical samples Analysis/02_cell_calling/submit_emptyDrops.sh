#!/bin/bash
## =============================================================================
## Neuroblastoma Patient scRNA-seq - EmptyDrops Cell Calling
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
## =============================================================================
## Step 2 of scRNA-seq pipeline (after Cell Ranger)
## This script calls true cells vs empty droplets using EmptyDrops
## =============================================================================

#SBATCH --job-name=emptyDrops_NB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --array=1-18
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=USER_EMAIL
#SBATCH --mail-type=ALL

set -e  # Exit on error
set -u  # Exit on undefined variable

## =============================================================================
## CONFIGURATION - Edit these paths for your environment
## =============================================================================

BASE_DIR="/path/to/NB_patient_analysis"
CELLRANGER_DIR="${BASE_DIR}/cellranger_output"
OUTPUT_DIR="${BASE_DIR}/data/SCE_objects"
RELAXED_DIR="${OUTPUT_DIR}/relaxed_filters"
SCRIPT_DIR="${BASE_DIR}/src"

# Create output directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${RELAXED_DIR}

## =============================================================================
## SAMPLE MAPPING - NB Patients (18 samples)
## =============================================================================
## Format: "SAMPLE_ID:PATIENT_ID"
## All 18 NB samples treated identically

declare -a SAMPLES=(
    "S35_0010:770-F2-i"       # Array 1
    "S35_0011:770-F1-i"       # Array 2
    "S35_0012:637-F1-i"       # Array 3
    "S35_0013:649-F2-i"       # Array 4
    "S35_0017:644-F1-i"       # Array 5
    "S35_0018:1087-F1-i"      # Array 6
    "S35_0019:637-F2-i"       # Array 7
    "S35_0020:649-F1-i"       # Array 8
    "S35_0021:1087-F2-i"      # Array 9
    "S35_0025:637-F2-ii"      # Array 10
    "S35_0026:644-F1-ii"      # Array 11
    "S35_0027:1087-F1-ii"     # Array 12
    "S35_0028:770-F1-ii"      # Array 13 (resequenced)
    "S35_0029:770-F2-ii"      # Array 14
    "S35_0033:637-F1-ii"      # Array 15
    "S35_0034:1087-F2-ii"     # Array 16
    "S35_0035:649-F2-ii"      # Array 17
    "S35_0036:649-F1-ii"      # Array 18
)

## =============================================================================
## ACTIVATE CONDA ENVIRONMENT
## =============================================================================

module use /opt/software/easybuild/modules/all/
module load Mamba

source ~/.bashrc
mamba activate hotspot

## =============================================================================
## GET SAMPLE INFO FOR THIS ARRAY TASK
## =============================================================================

TASK_ID=${SLURM_ARRAY_TASK_ID}

# Get sample info
SAMPLE_INFO="${SAMPLES[$TASK_ID-1]}"
IFS=':' read -r SAMPLE_ID PATIENT_ID <<< "$SAMPLE_INFO"

echo "================================================================================"
echo "EmptyDrops Cell Calling - Neuroblastoma Patient"
echo "================================================================================"
echo "Array Task:    $TASK_ID"
echo "Sample ID:     $SAMPLE_ID"
echo "Patient ID:    $PATIENT_ID"
echo "Start Time:    $(date)"
echo "================================================================================"

## =============================================================================
## CHECK INPUT FILES
## =============================================================================

# Path to Cell Ranger output
H5_FILE="${CELLRANGER_DIR}/${SAMPLE_ID}/outs/raw_feature_bc_matrix.h5"

if [ ! -f "$H5_FILE" ]; then
    echo "ERROR: H5 file not found: $H5_FILE"
    echo "Make sure Cell Ranger has completed for this sample"
    exit 1
fi

echo "Input H5 file: $H5_FILE"
echo "File size:     $(du -h $H5_FILE | cut -f1)"

## =============================================================================
## RUN EMPTYDROPS
## =============================================================================

OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_SCE.RDS"

echo ""
echo "================================================================================"
echo "Running EmptyDrops"
echo "================================================================================"
echo "Sample:        $SAMPLE_ID"
echo "Input:         $H5_FILE"
echo "Output:        $RELAXED_DIR/${SAMPLE_ID}_SCE.RDS"
echo "Cores:         ${SLURM_CPUS_PER_TASK}"
echo "UMI Threshold: 500"
echo "Doublet Rate:  0.04"
echo "================================================================================"

Rscript ${SCRIPT_DIR}/emptyDrops_cell_calling.R \
    --h5file=${H5_FILE} \
    --id=${SAMPLE_ID} \
    --ncores=${SLURM_CPUS_PER_TASK} \
    --umithreshold=500 \
    --output=${OUTPUT_FILE} \
    --doubletrate=0.04

## =============================================================================
## CHECK COMPLETION
## =============================================================================

RELAXED_OUTPUT="${RELAXED_DIR}/${SAMPLE_ID}_SCE.RDS"
DIAG_FILE="${RELAXED_DIR}/${SAMPLE_ID}_SCE_emptyDrops_diagnostics.tsv"

if [ $? -eq 0 ] && [ -f "$RELAXED_OUTPUT" ]; then
    echo ""
    echo "================================================================================"
    echo "EmptyDrops COMPLETED SUCCESSFULLY"
    echo "================================================================================"
    echo "Sample:     $SAMPLE_ID ($PATIENT_ID)"
    echo "Output:     $RELAXED_OUTPUT"
    echo "File size:  $(du -h $RELAXED_OUTPUT | cut -f1)"
    if [ -f "$DIAG_FILE" ]; then
        echo "Diagnostics: $DIAG_FILE"
        echo "Diag size:   $(du -h $DIAG_FILE | cut -f1)"
    fi
    echo "End Time:   $(date)"
    echo "================================================================================"
else
    echo ""
    echo "================================================================================"
    echo "ERROR: EmptyDrops FAILED"
    echo "================================================================================"
    echo "Sample:     $SAMPLE_ID ($PATIENT_ID)"
    echo "Check logs: ${BASE_DIR}/logs/"
    echo "================================================================================"
    exit 1
fi
