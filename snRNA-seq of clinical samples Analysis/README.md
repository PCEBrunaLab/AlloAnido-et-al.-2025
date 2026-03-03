# Neuroblastoma Patient scRNA-seq Analysis Pipeline

Single-cell RNA sequencing analysis pipeline for neuroblastoma (NB) patient tumour biopsies, from raw FASTQ processing through malignant cell characterisation, gene program discovery, and CNV detection.

**Author:** Ayeh Sadr

---

## Dataset Overview

Six matched samples at diagnosis and post-treatment were obtained from VIVO biobank.

- **Tumour type:** Neuroblastoma (NB)
- **Patients:** 5 (770, 637, 644, 649, 1087)
- **Timepoints:** F1 (diagnosis) and F2 (post-treatment)
- **Replicates:** Technical replicates (i, ii) per timepoint
- **Total samples:** 18
- **Sequencing:** 10x Genomics 3' GEX, sequencing run S35_GEX_run1896
- **Reference genome:** GRCh38 (refdata-gex-GRCh38-2020-A)

### Sample Summary

| Patient | F1 rep-i | F1 rep-ii | F2 rep-i | F2 rep-ii |
|---------|----------|-----------|----------|-----------|
| 770     | S35_0011 | S35_0028  | S35_0010 | S35_0029  |
| 637     | S35_0012 | S35_0033  | S35_0019 | S35_0025  |
| 644     | S35_0017 | S35_0026  | --       | --        |
| 649     | S35_0020 | S35_0036  | S35_0013 | S35_0035  |
| 1087    | S35_0018 | S35_0027  | S35_0021 | S35_0034  |

---

## Repository Structure

```
.
├── README.md
├── sample_metadata.csv
├── 01_cellranger/
│   └── submit_cellranger.sh
├── 02_cell_calling/
│   ├── submit_emptyDrops.sh
│   ├── emptyDrops_cell_calling.R
│   ├── submit_emptyDrops_report.sh
│   └── emptyDrops_report.R
├── 03_normalization/
│   ├── submit_normalise.sh
│   └── norm_counts.R
├── 04_qc_inspection/
│   ├── inspect_normalized_SCE.R
│   └── verify_final_object.R
├── 05_seurat_analysis/
│   ├── step1_seurat_pipeline.R
│   ├── step1b_find_markers.R
│   ├── step1c_annotate_clusters.R
│   ├── step1d_gene_signatures.R
│   └── step1e_jansky_annotation.R
├── 06_malignant_analysis/
│   └── recluster_malignant.R
├── 07_CNV_analysis/
│   ├── submit_inferCNV.sh
│   ├── inferCNV_analysis.R
│   ├── submit_copykat.sh
│   ├── copykat_analysis.R
│   ├── submit_visualize_inferCNV.sh
│   └── visualize_inferCNV.R
└── 08_downstream_analysis/
    ├── analyze_cluster_composition.R
    ├── analyze_longitudinal.R
    ├── analyze_replicates.R
    ├── validate_clustering.R
    └── correlation_analysis.R
```

---

## Pipeline Workflow

### Step 1: Cell Ranger Count
```bash
sbatch 01_cellranger/submit_cellranger.sh
```
Runs Cell Ranger 8.0.0 as a SLURM array job across all 18 samples.

### Step 2: EmptyDrops Cell Calling
```bash
sbatch 02_cell_calling/submit_emptyDrops.sh
sbatch 02_cell_calling/submit_emptyDrops_report.sh   # after all EmptyDrops complete
```
Relaxed cell-calling filters (NBAtlas pipeline): 0.1% FDR, >=200 genes, >=500 UMIs, adaptive MT filtering, scDblFinder doublet removal.

### Step 3: Normalization
```bash
sbatch 03_normalization/submit_normalise.sh
```
Merges per-sample SCE objects with patient metadata, deconvolution size factors, multi-batch normalisation.

### Step 4: QC Inspection
```bash
Rscript 04_qc_inspection/inspect_normalized_SCE.R
Rscript 04_qc_inspection/verify_final_object.R
```

### Step 5: Seurat Analysis Pipeline
```bash
Rscript 05_seurat_analysis/step1_seurat_pipeline.R     # SCE -> Seurat, SCTransform, Harmony
Rscript 05_seurat_analysis/step1b_find_markers.R       # FindAllMarkers
Rscript 05_seurat_analysis/step1c_annotate_clusters.R  # manual annotation
Rscript 05_seurat_analysis/step1d_gene_signatures.R    # ADRN/MES scoring
Rscript 05_seurat_analysis/step1e_jansky_annotation.R  # developmental stages
```

### Step 6: Malignant Cell Analysis
```bash
Rscript 06_malignant_analysis/recluster_malignant.R
```

### Step 7: CNV Detection
```bash
sbatch 07_CNV_analysis/submit_inferCNV.sh
sbatch 07_CNV_analysis/submit_visualize_inferCNV.sh
```
inferCNV using cluster 10 as diploid reference.

### Step 8: Downstream Analyses
```bash
Rscript 08_downstream_analysis/analyze_cluster_composition.R
Rscript 08_downstream_analysis/analyze_longitudinal.R
Rscript 08_downstream_analysis/analyze_replicates.R
Rscript 08_downstream_analysis/validate_clustering.R
Rscript 08_downstream_analysis/correlation_analysis.R
```

---

## Configuration

All scripts use config variables at the top. Edit the **CONFIGURATION** section before running.

### Shell Scripts
```bash
## CONFIGURATION
BASE_DIR="/path/to/NB_patient_analysis"
FASTQ_DIR="/path/to/S35_GEX_run1896"
TRANSCRIPTOME="/path/to/refdata-gex-GRCh38-2020-A"
USER_EMAIL="your.email@institution.ac.uk"
```

### R Scripts
```r
## CONFIGURATION
PROJECT_DIR <- "/path/to/NB_patient_analysis"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR    <- file.path(RESULTS_DIR, "data")
```

---

## Environment Setup

### Conda Environments
```bash
# Main analysis environment
mamba create -n hotspot
mamba activate hotspot
mamba install -c conda-forge r-base r-seurat r-ggplot2 r-patchwork r-optparse
mamba install -c bioconda bioconductor-singlecellexperiment bioconductor-dropletutils \
  bioconductor-scran bioconductor-scater bioconductor-batchelor \
  bioconductor-biomart bioconductor-scdblfinder bioconductor-harmony

# inferCNV (separate environment)
mamba create -n infercnv_env
mamba install -c bioconda bioconductor-infercnv
```

### Software Versions
- Cell Ranger 8.0.0
- R >= 4.2, Seurat v5
- Reference: GRCh38 (refdata-gex-GRCh38-2020-A)

---

## Required Reference Files

Place in `references/` directory:

- `refdata-gex-GRCh38-2020-A/` - Cell Ranger transcriptome
- `gencode.v43.annotation.gtf.gz` - Gene annotation (inferCNV)
- `adrenal_medulla_Seurat.RDS`, `adrenal_medulla_annot.RDS`, `adrenal_medulla_counts.RDS` - Jansky reference
- `vanGroningen_2017.xlsx` - ADRN/MES signatures
- `Dyer_Sig.xlsx` - Dyer signatures
- `Int_signature.csv` - Intermediate signature

---

## Key Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| EmptyDrops FDR | 0.001 (0.1%) | NBAtlas |
| Min detected genes | 200 | NBAtlas |
| Min UMI counts | 500 | NBAtlas |
| Doublet rate | 0.04 | scDblFinder |
| SCTransform features | 1000 | Seurat |
| Harmony PCs | 30 | Standard |
| Clustering resolution | 0.2 | Empirical |
| inferCNV cutoff | 0.1 | NBAtlas |

---

## Key Outputs

### Seurat Objects
- `NB_patients_seurat_harmony.rds` - Base with Harmony
- `NB_patients_seurat_with_signatures.rds` - With ADRN/MES scores
- `NB_patients_seurat_with_jansky.rds` - With developmental stages
- `NB_patients_malignant_reclustered.rds` - Malignant re-clustered

### Tables
- `malignant_cell_metadata_complete.csv`
- `malignant_cluster_markers.csv`
- `emptyDrops_summary_NB_relaxed_filters.tsv`


