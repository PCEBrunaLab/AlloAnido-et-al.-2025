# CUT&Tag Analysis Pipeline

A comprehensive bioinformatics pipeline for analyzing CUT&Tag (Cleavage Under Targets and Tagmentation) data, from raw sequencing reads to peak calling, quality control, and visualization.

## Overview

This pipeline processes CUT&Tag data through mapping, quality control, duplicate removal, peak calling, and visualization. It is optimized for analyzing histone modifications (H3K4me3) and histone dem ethylases (KDM5A, KDM5B) across multiple experimental conditions.

## Table of Contents

1. [Requirements](#requirements)
2. [Pipeline Steps](#pipeline-steps)
3. [Quick Start](#quick-start)
4. [Detailed Workflow](#detailed-workflow)
5. [Output Files](#output-files)
6. [Customization](#customization)

## Requirements

### Software Dependencies

- **Bowtie2** (v2.4.2 or higher) - Read alignment
- **SAMtools** (v1.11 or higher) - BAM file manipulation
- **Picard Tools** (v2.23.8 or higher) - Duplicate marking and removal
- **MACS2** (v2.2.7.1 or higher) - Peak calling
- **deepTools** (v3.5.6 or higher) - BigWig generation and visualization
- **bedtools** (v2.29.2 or higher) - Genomic interval operations
- **R** (v4.1.0 or higher) with packages:
  - tidyverse
  - GenomicRanges
  - Rsamtools
  - ggplot2
  - ggpubr

### Reference Files

- **Genome Index**: Bowtie2 index for GRCh38 human genome
- **Gene Annotations**: BED format file with gene coordinates (`genes.bed`)
- **Gene Lists** (optional): CSV files with gene sets for condition-specific analysis

### Computational Resources

- HPC cluster with SLURM workload manager (scripts can be adapted for other systems)
- Recommended: 8-128 GB RAM depending on step
- Recommended: 8 CPU cores per job

## Pipeline Steps

### Step 1: Read Mapping
**Script**: `1_Mapping_CUT_Tag_R85.sh`

Maps paired-end FASTQ files to the reference genome using Bowtie2 with CUT&Tag-specific parameters.

**Parameters**:
- `--local`: Local alignment mode for CUT&Tag
- `--very-sensitive`: High sensitivity alignment
- `-I 10 -X 700`: Insert size range suitable for CUT&Tag

**Input**: Paired-end FASTQ files (`.fastq.gz`)  
**Output**: Sorted BAM files (`*_sorted.bam`), alignment statistics

### Step 2: Quality Control
**Script**: `2_CUTTag_QC_R85.sh`

Generates comprehensive QC metrics and visualizations for all samples.

**Metrics Calculated**:
- Total and mapped reads
- Alignment rates
- Duplication rates
- Unique library complexity
- Sequencing depth

**Output**: QC summary table, QC plots (PDF)

### Step 3: Duplicate Removal
**Script**: `3_Remove_Dup.sh`

Removes PCR duplicates using Picard MarkDuplicates and generates before/after comparison plots.

**Input**: Sorted BAM files  
**Output**: Deduplicated BAM files (`*_dedup.bam`), comparison metrics

### Step 4: Peak Calling
**Script**: `4_MACS2_R85_Untreated_Control.sh`

Calls peaks using MACS2 with multiple parameter combinations for optimization.

**Features**:
- Separate peak calling for each replicate
- Multiple q-value thresholds tested (0.1, 0.001, 0.00001)
- Parameter testing (default, nomodel, broad, nomodel+broad)
- Uses untreated samples as control

**Input**: Deduplicated BAM files  
**Output**: Peak files (narrowPeak/broadPeak format), peak summaries

### Step 5: Peak Quality Analysis
**Script**: `5_FRiP.R`

Calculates Fraction of Reads in Peaks (FRiP) and generates peak quality metrics.

**Metrics**:
- Peak numbers per condition
- Peak width distributions
- FRiP scores
- Replicate consistency

**Output**: Peak quality reports, FRiP plots (PDF)

### Step 6: BigWig Generation
**Script**: `6_Format_BigWig_R85.sh`

Creates normalized BigWig files for genome browser visualization.

**Normalization**: RPGC (Reads Per Genomic Content)  
**Features**:
- Individual replicate BigWigs
- Merged replicate BigWigs
- Binsize: 10 bp

**Input**: Deduplicated BAM files  
**Output**: Normalized BigWig files (`.bw`)

### Step 7: Genome-wide Heatmaps
**Script**: `7_Heatmap_R85_CUTTag.sh`

Generates genome-wide heatmaps showing signal at TSS and peak regions across time points.

**Features**:
- TSS-centered heatmaps (±3 kb)
- Peak-centered heatmaps (±3 kb)
- Multiple Homer peak merging distances (5kb, 7.5kb, 10kb)
- Time series visualization

**Output**: Heatmap PDFs, profile plots, DeepTools matrices

### Step 8: Gene-Specific Analysis
**Script**: `8_GeneSpecific_Peak_Analysis.sh`

Analyzes specific gene sets defined in CSV files.

**Features**:
- TSS-centered analysis for gene lists
- Peak assignment within 5 kb of genes
- Active gene identification (H3K4me3-based)
- KDM5A/B presence at active genes

**Input**: Gene list CSV, gene annotations (BED)  
**Output**: Gene-specific heatmaps, peak assignments, summaries

### Step 9: Cell State Analysis
**Script**: `9_CellState_GeneList_Analysis.sh`

Analyzes cell state-specific gene signatures (ADRN, MES, Intermediate).

**Features**:
- TSS-centered analysis
- Gene body-scaled analysis (full gene coverage)
- Peak-centered analysis
- Active gene identification per cell state

**Input**: Cell state gene lists (CSV)  
**Output**: Cell state-specific heatmaps, summaries


## Detailed Workflow

### Data Organization

```
project/
├── raw_fastq/              # Raw FASTQ files
├── bam_files/              # Aligned BAM files
├── bam_files_dedup/        # Deduplicated BAM files
├── QC_final/               # QC metrics and plots
├── peakCalling/
│   └── MACS2/              # Peak files
├── alignment/
│   └── bigwig/             # Normalized coverage tracks
├── heatmap/
│   ├── final_timeseries/   # Genome-wide heatmaps
│   ├── gene_specific_peaks/# Gene-specific analysis
│   └── cellstate_genelist/ # Cell state analysis
├── Reference/              # Gene annotations and lists
└── logs/                   # Job logs and errors
```


