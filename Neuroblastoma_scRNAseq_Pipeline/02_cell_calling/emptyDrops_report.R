#!/usr/bin/env Rscript
## =============================================================================
## EmptyDrops Summary Report
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Generate summary table of cell calling results across all samples
## Input:   SCE objects and filter info files
## Output:  TSV summary table with cell counts and filtering statistics
## =============================================================================

library(SingleCellExperiment)
library(DropletUtils)
library(optparse)

## =============================================================================
## Parse Command Line Arguments
## =============================================================================

parser <- OptionParser()
parser <- add_option(parser, c("-d", "--datadir"), type="character",
                     help="Directory containing SCE objects")
parser <- add_option(parser, c("-c", "--cellrangerdir"), type="character",
                     help="Directory containing Cell Ranger outputs")
parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file for summary table (TSV)")
parser <- add_option(parser, c("-p", "--prefix"), type="character",
                     default="NB", help="Sample prefix for report header")

opt <- parse_args(parser)

## =============================================================================
## Find All Sample Files
## =============================================================================

message("Searching for SCE objects in: ", opt$datadir)
sce_files <- list.files(opt$datadir, pattern="_SCE\\.RDS$", full.names=TRUE)
message("Found ", length(sce_files), " SCE objects")

## =============================================================================
## Extract Summary Statistics for Each Sample
## =============================================================================

results <- data.frame()

for (sce_file in sce_files) {

    sample_id <- gsub("_SCE\\.RDS$", "", basename(sce_file))
    message("\nProcessing: ", sample_id)

    sce <- readRDS(sce_file)
    final_cells <- ncol(sce)
    message("  Final cells after QC: ", final_cells)

    filter_file <- gsub("_SCE\\.RDS$", "_filterInfo.tsv", sce_file)
    if (file.exists(filter_file)) {
        filter_info <- read.table(filter_file, header=TRUE, sep="\t", check.names=FALSE)
        if ("Gt1kUMIs" %in% colnames(filter_info)) {
            cells_called <- filter_info$N.NotEmpty
            cells_post_gene <- filter_info$Gt1kUMIs
            removed_low_complexity <- filter_info$N.NotEmpty - filter_info$Gt1kUMIs
            removed_high_mt <- filter_info$Gt1kUMIs - final_cells
            removed_doublets <- if ("Removed_Doublets" %in% colnames(filter_info)) filter_info$Removed_Doublets else 0
        } else if ("Gt500Genes" %in% colnames(filter_info)) {
            cells_called <- filter_info$N.NotEmpty
            cells_post_gene <- filter_info$Gt500Genes
            removed_low_complexity <- filter_info$Removed_Low_Genes
            removed_high_mt <- filter_info$Removed_High_MT
            removed_doublets <- if ("Removed_Doublets" %in% colnames(filter_info)) filter_info$Removed_Doublets else 0
        } else if ("Gt200Genes" %in% colnames(filter_info)) {
            cells_called <- filter_info$N.NotEmpty
            cells_post_gene <- filter_info$Gt200Genes
            removed_low_complexity <- filter_info$Removed_Low_Genes
            removed_high_mt <- filter_info$Removed_High_MT
            removed_doublets <- if ("Removed_Doublets" %in% colnames(filter_info)) filter_info$Removed_Doublets else 0
        } else {
            warning("  Unrecognised filterInfo columns in ", filter_file)
            cells_called <- filter_info$N.NotEmpty
            cells_post_gene <- NA
            removed_low_complexity <- NA
            removed_high_mt <- NA
            removed_doublets <- if ("Removed_Doublets" %in% colnames(filter_info)) filter_info$Removed_Doublets else NA
        }
    } else {
        message("  Warning: Filter info file not found")
        cells_called <- NA
        cells_post_gene <- NA
        removed_low_complexity <- NA
        removed_high_mt <- NA
        removed_doublets <- NA
    }

    h5_file <- file.path(opt$cellrangerdir, sample_id, "outs", "raw_feature_bc_matrix.h5")
    if (file.exists(h5_file)) {
        raw_sce <- read10xCounts(h5_file, sample.names=sample_id)
        total_droplets <- ncol(raw_sce)
        empty_droplets <- total_droplets - cells_called
        message("  Total droplets: ", total_droplets)
        message("  Empty droplets: ", empty_droplets)
        message("  Called cells: ", cells_called)
        rm(raw_sce)
    } else {
        message("  Warning: H5 file not found")
        total_droplets <- NA
        empty_droplets <- NA
    }

    if (is.na(removed_low_complexity)) {
        cells_low_complexity <- if (!is.na(cells_called) && !is.na(cells_post_gene)) cells_called - cells_post_gene else NA
    } else {
        cells_low_complexity <- removed_low_complexity
    }
    if (is.na(removed_high_mt)) {
        cells_high_mt <- if (!is.na(cells_post_gene)) cells_post_gene - final_cells else NA
    } else {
        cells_high_mt <- removed_high_mt
    }

    results <- rbind(results, data.frame(
        Sample_ID = sample_id,
        Total_Droplets = total_droplets,
        Empty_Droplets = empty_droplets,
        Cells_Called_EmptyDrops = cells_called,
        Removed_Low_Complexity = cells_low_complexity,
        Removed_High_MT = cells_high_mt,
        Removed_Doublets = removed_doublets,
        Final_Cells = final_cells,
        Percent_Cells_Retained = ifelse(is.na(total_droplets) || total_droplets == 0,
                                        NA,
                                        round(100 * final_cells / total_droplets, 2))
    ))
}

if (nrow(results) == 0) {
    warning("No SCE files found, writing empty report.")
    write.table(data.frame(), file=opt$output, sep="\t", quote=FALSE, row.names=FALSE)
    quit(save="no", status=0)
}

## =============================================================================
## Add Summary Row
## =============================================================================

summary_row <- data.frame(
    Sample_ID = "TOTAL",
    Total_Droplets = sum(results$Total_Droplets, na.rm=TRUE),
    Empty_Droplets = sum(results$Empty_Droplets, na.rm=TRUE),
    Cells_Called_EmptyDrops = sum(results$Cells_Called_EmptyDrops, na.rm=TRUE),
    Removed_Low_Complexity = sum(results$Removed_Low_Complexity, na.rm=TRUE),
    Removed_High_MT = sum(results$Removed_High_MT, na.rm=TRUE),
    Removed_Doublets = sum(results$Removed_Doublets, na.rm=TRUE),
    Final_Cells = sum(results$Final_Cells, na.rm=TRUE),
    Percent_Cells_Retained = {
        total_droplets_sum <- sum(results$Total_Droplets, na.rm=TRUE)
        if (total_droplets_sum == 0) NA else round(100 * sum(results$Final_Cells, na.rm=TRUE) / total_droplets_sum, 2)
    }
)

results <- rbind(results, summary_row)

## =============================================================================
## Save Summary Table
## =============================================================================

message("\nSaving summary table to: ", opt$output)
write.table(results, file=opt$output, sep="\t", quote=FALSE, row.names=FALSE)

## =============================================================================
## Print Summary to Console
## =============================================================================

message("\n================================================================================")
message("EMPTYDROPS SUMMARY - ", opt$prefix, " SAMPLES")
message("================================================================================")
print(results)
message("================================================================================")
message("Summary table saved to: ", opt$output)
message("================================================================================")
