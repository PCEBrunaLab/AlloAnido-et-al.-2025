# BASiCS pipeline

## Overview

This workflow was used for the data preparation and analysis of single-cell RNA-sequencing data to quantify levels of transcriptioanl noise. This analysis was performed either between cell states in untreated cells (neuroblastoma and hepatoblastoma models) or between phases of a cisplatin treatment course.

For additional information and instructions see the BASiCS GitHub page (https://github.com/catavallejos/BASiCS) and relevant publications (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004333).


## Structure

1. **Data preparation** 
   - Load pre-processed scRNA-seq data (see https://github.com/PCEBrunaLab/Roux-et-al.-2024/tree/main/Single_Cell_RNA_Sequencing), add replicate information, remove lowly expressed genes and eparate data into individual objects.

2. **Generation of Markov chains**
   - Estimation of gene expression parameters (including mean and over-disperssion) with BASiCS_MCMC.

3. **Compare transcriptional noise levels between groups**
   - Calculate differential expression and over-dispersion between cell states using BASiCS_TestDE, plot results and extract gene information.
     
4. **Identify highly variable genes within individual groups**
   - Use BASiCS_DetectHVG to identify genes with hihg variability within phases of cisplatin treatment course.
