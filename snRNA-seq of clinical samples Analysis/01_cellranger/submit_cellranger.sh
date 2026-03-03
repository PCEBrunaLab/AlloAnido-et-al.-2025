#!/bin/bash
## =============================================================================
## Neuroblastoma Patient scRNA-seq - Cell Ranger Count Submission
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
## =============================================================================
## This script processes 18 NB patient samples through Cell Ranger count.
## Each sample has 5 lanes that Cell Ranger will automatically merge.
## =============================================================================

#SBATCH --job-name=cellranger_NB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=smp
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
FASTQ_DIR="/path/to/S35_GEX_run1896"
OUTPUT_DIR="${BASE_DIR}/cellranger_output"
TRANSCRIPTOME="/path/to/refdata-gex-GRCh38-2020-A"
CELLRANGER_PATH="/path/to/cellranger-8.0.0"
LOGS_DIR="${BASE_DIR}/logs"

# Create directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOGS_DIR}

## =============================================================================
## SAMPLE MAPPING - NB Patients (18 samples)
## =============================================================================
## Format: "GF_ID:SAMPLE_NAME:PATIENT_ID:EXPECTED_CELLS"
## All 18 NB samples treated identically

declare -a SAMPLES=(
    "S35_0010:S35_0010:770-F2-i:5000"       # Array 1  - Patient 770, F2, replicate i
    "S35_0011:S35_0011:770-F1-i:5000"       # Array 2  - Patient 770, F1, replicate i
    "S35_0012:S35_0012:637-F1-i:5000"       # Array 3  - Patient 637, F1, replicate i
    "S35_0013:S35_0013:649-F2-i:5000"       # Array 4  - Patient 649, F2, replicate i
    "S35_0017:S35_0017:644-F1-i:5000"       # Array 5  - Patient 644, F1, replicate i
    "S35_0018:S35_0018:1087-F1-i:5000"      # Array 6  - Patient 1087, F1, replicate i
    "S35_0019:S35_0019:637-F2-i:5000"       # Array 7  - Patient 637, F2, replicate i
    "S35_0020:S35_0020:649-F1-i:5000"       # Array 8  - Patient 649, F1, replicate i
    "S35_0021:S35_0021:1087-F2-i:5000"      # Array 9  - Patient 1087, F2, replicate i
    "S35_0025:S35_0025:637-F2-ii:5000"      # Array 10 - Patient 637, F2, replicate ii
    "S35_0026:S35_0026:644-F1-ii:5000"      # Array 11 - Patient 644, F1, replicate ii
    "S35_0027:S35_0027:1087-F1-ii:5000"     # Array 12 - Patient 1087, F1, replicate ii
    "S35_0028:S35_0028:770-F1-ii:5000"      # Array 13 - Patient 770, F1, replicate ii (resequenced)
    "S35_0029:S35_0029:770-F2-ii:5000"      # Array 14 - Patient 770, F2, replicate ii
    "S35_0033:S35_0033:637-F1-ii:5000"      # Array 15 - Patient 637, F1, replicate ii
    "S35_0034:S35_0034:1087-F2-ii:5000"     # Array 16 - Patient 1087, F2, replicate ii
    "S35_0035:S35_0035:649-F2-ii:5000"      # Array 17 - Patient 649, F2, replicate ii
    "S35_0036:S35_0036:649-F1-ii:5000"      # Array 18 - Patient 649, F1, replicate ii
)

## =============================================================================
## SETUP CELL RANGER
## =============================================================================

export PATH=${CELLRANGER_PATH}:$PATH

# Verify Cell Ranger is available
if ! command -v cellranger &> /dev/null; then
    echo "ERROR: Cell Ranger not found in PATH"
    echo "PATH: $PATH"
    exit 1
fi

echo "Cell Ranger version:"
cellranger --version

## =============================================================================
## GET SAMPLE INFO FOR THIS ARRAY TASK
## =============================================================================

TASK_ID=${SLURM_ARRAY_TASK_ID}

# Get sample info for this task
SAMPLE_INFO="${SAMPLES[$TASK_ID-1]}"  # Arrays are 0-indexed
IFS=':' read -r GF_ID SAMPLE_NAME PATIENT_ID EXPECTED_CELLS <<< "$SAMPLE_INFO"

echo "==============================================================================="
echo "NEUROBLASTOMA PATIENT - Cell Ranger Count"
echo "==============================================================================="
echo "Array Task ID:     $TASK_ID"
echo "GF_ID:             $GF_ID"
echo "Sample Name:       $SAMPLE_NAME"
echo "Patient ID:        $PATIENT_ID"
echo "Expected Cells:    $EXPECTED_CELLS"
echo "==============================================================================="

## =============================================================================
## VERIFY FASTQ FILES EXIST
## =============================================================================

echo "Checking for FASTQ files..."
FASTQ_COUNT=$(ls ${FASTQ_DIR}/${SAMPLE_NAME}_*.fastq.gz 2>/dev/null | wc -l)

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTQ files found for sample ${SAMPLE_NAME}"
    echo "Searched in: ${FASTQ_DIR}/${SAMPLE_NAME}_*.fastq.gz"
    exit 1
fi

echo "Found $FASTQ_COUNT FASTQ files for ${SAMPLE_NAME}"
echo "Expected: 10 files (5 lanes × 2 reads)"

if [ "$FASTQ_COUNT" -ne 10 ]; then
    echo "WARNING: Expected 10 files but found $FASTQ_COUNT"
fi

# Show first few files
echo "Sample files:"
ls ${FASTQ_DIR}/${SAMPLE_NAME}_*.fastq.gz | head -5
echo "..."

## =============================================================================
## RUN CELL RANGER COUNT
## =============================================================================

echo ""
echo "==============================================================================="
echo "Starting Cell Ranger Count"
echo "==============================================================================="
echo "Output directory: ${OUTPUT_DIR}/${GF_ID}"
echo "Start time: $(date)"
echo "==============================================================================="

cd ${OUTPUT_DIR}

cellranger count \
    --id=${GF_ID} \
    --transcriptome=${TRANSCRIPTOME} \
    --fastqs=${FASTQ_DIR} \
    --sample=${SAMPLE_NAME} \
    --expect-cells=${EXPECTED_CELLS} \
    --create-bam=true \
    --localcores=${SLURM_CPUS_PER_TASK} \
    --localmem=120

## =============================================================================
## CHECK COMPLETION
## =============================================================================

if [ $? -eq 0 ]; then
    echo ""
    echo "==============================================================================="
    echo "Cell Ranger Count COMPLETED SUCCESSFULLY"
    echo "==============================================================================="
    echo "Sample:     ${SAMPLE_NAME} (${PATIENT_ID})"
    echo "Output:     ${OUTPUT_DIR}/${GF_ID}"
    echo "End time:   $(date)"

    # Check for output files
    if [ -f "${OUTPUT_DIR}/${GF_ID}/outs/web_summary.html" ]; then
        echo "QC Report:  ${OUTPUT_DIR}/${GF_ID}/outs/web_summary.html"
    fi

    if [ -d "${OUTPUT_DIR}/${GF_ID}/outs/filtered_feature_bc_matrix" ]; then
        BARCODE_COUNT=$(zcat ${OUTPUT_DIR}/${GF_ID}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l)
        echo "Cells:      $BARCODE_COUNT"
    fi

    echo "==============================================================================="
else
    echo ""
    echo "==============================================================================="
    echo "ERROR: Cell Ranger Count FAILED"
    echo "==============================================================================="
    echo "Sample:     ${SAMPLE_NAME} (${PATIENT_ID})"
    echo "Check logs: ${LOGS_DIR}/"
    echo "==============================================================================="
    exit 1
fi
