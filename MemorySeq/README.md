# MemorySeq RNA-seq Analysis Pipeline

This repository contains scripts and analysis code for processing and analyzing RNA-seq data from a MemorySeq experiment. The analysis follows the methodology described in the paper "Memory Sequencing Reveals Heritable Single-Cell Gene Expression Programs Associated with Distinct Cellular Behaviors".

## Overview

The pipeline processes RNA-seq data from MemorySeq experiments, performs quality control, data normalization, and identifies heritable gene expression patterns. The analysis workflow includes batch effect correction, variability analysis, and correlation analysis to identify genes with heritable expression patterns.

## Pipeline Components

### 1. Pre-processing Scripts (SLURM)

The pipeline uses SLURM for high-performance computing:

- `1.QC_Memoryseq.sh`: Quality control using FastQC and MultiQC
- `2.Trim_Memoryseq.sh`: Adapter trimming using Trim Galore
- `3.STAR_Index_Genome.sh`: Indexing reference genome for STAR alignment
- `4.STAR_Mapping_Memoryseq.sh`: RNA-seq read alignment using STAR
- `5.FeatureCounts.sh`: Counting reads in genomic features

### 2. R Analysis Pipeline

The main analysis is performed in R, following these key steps:

#### Initial Setup and Data Loading
- Environment setup with renv for reproducibility
- Data preprocessing and annotation
- TPM calculation and filtering of low-quality samples/genes

#### Batch Effect Analysis and Correction
- PCA analysis to identify batch effects
- Batch correction using limma's removeBatchEffect
- Visual assessment of correction effectiveness

#### Variability Analysis
- Coefficient of variation (CV) calculation
- Poisson regression modeling to identify highly variable genes
- Identification of heritable genes based on residual values

#### Correlation Analysis
- Analysis of co-expression patterns among heritable genes
- Hierarchical clustering to identify gene modules
- Visualization of gene clusters and correlations



## Analysis Workflow

1. **Quality Control**: Assessment of raw sequencing data quality
2. **Read Trimming**: Removal of adapters and low-quality sequences
3. **Alignment**: Mapping reads to the human reference genome
4. **Quantification**: Counting reads per gene
5. **Data Preprocessing**: Filtering, normalization, and annotation
6. **Batch Effect Correction**: Removing technical variability
7. **Variability Analysis**: Identifying genes with heritable expression patterns
8. **Correlation Analysis**: Exploring co-expression patterns

## Requirements

### SLURM Scripts
- FastQC 0.11.9
- MultiQC 1.9
- Trim Galore 0.6.6
- STAR 2.7.6a
- SAMtools 1.11
- featureCounts (Subread package)

### R Analysis
- R 4.0+
- renv (for environment management)
- Required R packages:
  - DESeq2
  - GenomicFeatures
  - biomaRt
  - dplyr
  - tidyr
  - data.table
  - ggplot2
  - limma
  - moments
  - pheatmap
  - RColorBrewer
  - dendextend

## Usage

### 1. SLURM Scripts
Run the SLURM scripts in sequence:

```bash
sbatch 1.QC_Memoryseq.sh
sbatch 2.Trim_Memoryseq.sh
sbatch 3.STAR_Index_Genome.sh
sbatch 4.STAR_Mapping_Memoryseq.sh
sbatch 5.FeatureCounts.sh
```

### 2. R Analysis

The R script `MemorySeq_RNA.R` contains the complete analysis workflow:

```R
# Run the full analysis script
source("MemorySeq_RNA.R")
```
