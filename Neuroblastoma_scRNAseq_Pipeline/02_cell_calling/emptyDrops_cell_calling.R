#!/usr/bin/env Rscript
## =============================================================================
## EmptyDrops Cell Calling (Relaxed Filters for Patient Biopsies)
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Cell calling using emptyDrops with relaxed quality thresholds
##          suitable for patient biopsy samples (0.1% FDR, >=200 genes,
##          >=500 UMIs, adaptive MT filtering, scDblFinder doublet removal)
## Input:   10X raw_feature_bc_matrix.h5 file
## Output:  SingleCellExperiment RDS object and filtering summary
## =============================================================================

# Required packages
required.packages <- c(
  "devtools", "SingleCellExperiment", "DropletUtils", "scran", "Matrix",
  "MatrixGenerics", "BiocParallel", "biomaRt", "optparse", "scater",
  "Cairo", "scDblFinder"
)

# Load libraries
library(devtools)
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(scater)
library(Matrix)
library(MatrixGenerics)
library(BiocParallel)
library(biomaRt)
library(optparse)
library(scDblFinder)

start <- Sys.time()

## =============================================================================
## Parse Command Line Arguments
## =============================================================================

parser <- OptionParser()
parser <- add_option(parser, c("-a", "--h5file"), type = "character",
                     help = "Path to 10X raw_feature_bc_matrix.h5 file")
parser <- add_option(parser, c("-i", "--id"), type = "character",
                     help = "Sample identifier")
parser <- add_option(parser, c("-n", "--ncores"), type = "numeric",
                     help = "Number of cores for parallel steps")
parser <- add_option(parser, c("-u", "--umithreshold"), type = "numeric",
                     help = "UMI threshold defining ambient background for emptyDrops")
parser <- add_option(parser, c("-o", "--output"), type = "character",
                     help = "Desired output path for SCE object (filename only)")
parser <- add_option(parser, c("-d", "--doubletrate"), type = "numeric", default = NA,
                     help = "Expected doublet rate (e.g. 0.04 for 4%%). Leave NA for auto-estimate")

opt <- parse_args(parser)

if (is.null(opt$h5file) || is.null(opt$id) || is.null(opt$ncores) ||
    is.null(opt$umithreshold) || is.null(opt$output)) {
  stop("Missing required arguments. Please provide --h5file, --id, --ncores, --umithreshold, and --output.")
}

if (!file.exists(opt$h5file)) {
  stop("Input HDF5 file not found: ", opt$h5file)
}

## =============================================================================
## Prepare Output Location (relaxed_filters subdirectory)
## =============================================================================

output_dir <- dirname(opt$output)
if (output_dir == "") {
  output_dir <- getwd()
}
output_basename <- basename(opt$output)
relaxed_dir <- file.path(output_dir, "relaxed_filters")
if (!dir.exists(relaxed_dir)) {
  dir.create(relaxed_dir, recursive = TRUE, showWarnings = FALSE)
}
output_path <- file.path(relaxed_dir, output_basename)
filterinfo_path <- gsub("_SCE\\.RDS$", "_filterInfo.tsv", output_path)

message("Prepared relaxed output path: ", output_path)

## =============================================================================
## Parallel backend
## =============================================================================

mcparam <- MulticoreParam(workers = opt$ncores)
register(mcparam)

## =============================================================================
## Load raw counts
## =============================================================================

message("Reading 10X matrix: ", opt$h5file)
in.sce <- read10xCounts(samples = opt$h5file, sample.names = opt$id, BPPARAM = mcparam)

## =============================================================================
## EmptyDrops with 0.1% FDR
## =============================================================================

message("Running emptyDrops (FDR ≤ 0.001)")
intest_calls <- emptyDrops(in.sce, assay.type = "counts", niters = 20000,
                           ignore = 4999, BPPARAM = mcparam,
                           lower = opt$umithreshold, retain = Inf)

sig_cells <- intest_calls$FDR <= 0.001 & !is.na(intest_calls$FDR)
message("Keeping ", sum(sig_cells), " barcodes passing FDR 0.001")
sce.drops <- in.sce[, sig_cells]

cells_called <- sum(sig_cells)

## =============================================================================
## Filter by detected genes (≥ 200) and UMI counts (≥ 500)
## =============================================================================

nonzero.genes <- MatrixGenerics::colSums2(counts(sce.drops) > 0)
total.counts <- MatrixGenerics::colSums2(counts(sce.drops))

message("Filtering cells with ≥ 200 detected genes AND ≥ 500 UMI counts")
sce.drops <- sce.drops[, (nonzero.genes >= 200) & (total.counts >= 500)]

cells_after_gene <- ncol(sce.drops)

## =============================================================================
## Remove high mitochondrial fraction (Adaptive filtering)
## =============================================================================

message("Querying BioMart for mitochondrial gene IDs")
biomaRt.connection <- useMart("ensembl", "hsapiens_gene_ensembl")
gene.df <- getBM(attributes = c("ensembl_gene_id", "chromosome_name"),
                 filters = "ensembl_gene_id",
                 values = rownames(sce.drops),
                 mart = biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id

mt.genes <- gene.df$ensembl_gene_id[gene.df$chromosome_name == "MT"]
mt.counts <- counts(sce.drops)[rownames(sce.drops) %in% mt.genes, , drop = FALSE]
mt.fraction <- MatrixGenerics::colSums2(mt.counts) / MatrixGenerics::colSums2(counts(sce.drops))

message("Filtering cells with high mitochondrial content (adaptive)")
mt.p <- pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction) * 5, lower.tail = FALSE)
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.1)])
message("  MT fraction threshold inferred: ", signif(mt.lim, 3), " (", signif(mt.lim * 100, 3), "%)")

sce.drops <- sce.drops[, mt.fraction < mt.lim]

cells_after_mt <- ncol(sce.drops)

## =============================================================================
## Doublet detection (scDblFinder)
## =============================================================================

set.seed(12345)
message("Running scDblFinder for doublet detection")
sce.drops <- scDblFinder(sce.drops,
                         dbr = if (is.na(opt$doubletrate)) NULL else opt$doubletrate,
                         BPPARAM = mcparam,
                         verbose = FALSE)

doublet_labels <- colData(sce.drops)$scDblFinder.class
sce.drops <- sce.drops[, doublet_labels == "singlet"]

cells_final <- ncol(sce.drops)
removed_doublets <- sum(doublet_labels != "singlet")

message("Removed ", removed_doublets, " predicted doublets")

## =============================================================================
## Save outputs
## =============================================================================

message("Writing relaxed-filter SCE to: ", output_path)
saveRDS(sce.drops, file = output_path)

message("Saving filtering summary to: ", filterinfo_path)

filt.df <- data.frame(
  Sample_ID = opt$id,
  N.NotEmpty = cells_called,
  Removed_Low_Genes = cells_called - cells_after_gene,
  Gt200Genes = cells_after_gene,
  Removed_High_MT = cells_after_gene - cells_after_mt,
  Post_MT = cells_after_mt,
  Removed_Doublets = removed_doublets,
  Final_Cells = cells_final
)

write.table(filt.df, file = filterinfo_path, sep = "\t", quote = FALSE, row.names = FALSE)

diagnostics_path <- gsub("_SCE\\.RDS$", "_emptyDrops_diagnostics.tsv", output_path)
message("Saving emptyDrops diagnostics to: ", diagnostics_path)

ed_diagnostics <- data.frame(
  Barcode = colnames(in.sce),
  Total_UMI = intest_calls$Total,
  LogProb = intest_calls$LogProb,
  PValue = intest_calls$PValue,
  FDR = intest_calls$FDR,
  Limited = intest_calls$Limited,
  stringsAsFactors = FALSE
)

write.table(ed_diagnostics, file = diagnostics_path, sep = "\t",
            quote = FALSE, row.names = FALSE)

message("All done. Time taken: ", Sys.time() - start)
