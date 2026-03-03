#!/usr/bin/env Rscript
## =============================================================================
## Normalization for Non-Multiplexed Patient Samples
## =============================================================================
## Author: Ayeh Sadr
## Pipeline: Neuroblastoma Patient scRNA-seq Analysis
##
## Purpose: Library size normalization and log-transformation with optional
##          batch processing for large datasets
## Input:   SCE .RDS files from cell calling step
## Output:  Combined normalized SingleCellExperiment object
## =============================================================================

library(DropletUtils)
library(SingleCellExperiment)
library(scran)
library(scater)
library(batchelor)
library(irlba)
library(optparse)
library(Matrix)
library(umap)
library(BiocParallel)

## =============================================================================
## Parse Command Line Arguments
## =============================================================================

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--scedir"), type="character",
                     help="Directory containing SCE .RDS files")

parser <- add_option(parser, c("-b", "--breaks"), action="store_true", default=FALSE,
                     help="Break data into chunks for size factor estimation (matches original pipeline logic)")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file path for combined normalized SCE object")

opt <- parse_args(parser)

if(is.null(opt$scedir)){
  stop("Missing --scedir argument specifying the directory containing SCE .RDS files")
}

if(!dir.exists(opt$scedir)){
  stop("The provided --scedir path does not exist: ", opt$scedir)
}

if(is.null(opt$output)){
  stop("Missing --output argument for the combined normalized SCE file")
}

## =============================================================================
## Read Sample Metadata for Paired Sample Tracking
## =============================================================================

data_dir <- dirname(dirname(opt$scedir))
metadata_file <- file.path(data_dir, "sample_metadata.csv")

message("Looking for metadata at: ", metadata_file)

if(file.exists(metadata_file)){
  message("Metadata file found!")
  message("Loading sample metadata from: ", metadata_file)

  sample_meta <- read.csv(metadata_file, stringsAsFactors=FALSE)

  if(!"GF_ID" %in% colnames(sample_meta)){
    stop("ERROR: Metadata file found but GF_ID column is missing!")
  }

  rownames(sample_meta) <- sample_meta$GF_ID

  message("Successfully loaded metadata for ", nrow(sample_meta), " samples")
  message("  Metadata columns: ", paste(colnames(sample_meta), collapse=", "))
  message("  Sample IDs: ", paste(head(sample_meta$GF_ID, 6), collapse=", "),
          if(nrow(sample_meta) > 6) "..." else "")

  required_cols <- c("GF_ID", "Patient_ID", "Tumor_Type", "Fraction", "Replicate", "Sequencing_Run")
  missing_cols <- required_cols[!required_cols %in% colnames(sample_meta)]
  if(length(missing_cols) > 0){
    stop("ERROR: Metadata is missing required columns: ", paste(missing_cols, collapse=", "))
  }
  message("All required metadata columns present")

} else {
  stop("ERROR: Metadata file not found at: ", metadata_file,
       "\nPlease ensure sample_metadata.csv exists in the data/ directory.",
       "\nSCE directory: ", opt$scedir,
       "\nExpected metadata: ", metadata_file)
}

## =============================================================================
## Read All Patient SCE Files
## =============================================================================

message("Looking for SCE files in: ", opt$scedir)
setwd(opt$scedir)

sce.file.list <- list.files(pattern="*_SCE.RDS")
message("Found ", length(sce.file.list), " SCE files")
print(sce.file.list)

sce.names <- gsub(sce.file.list, pattern="(.*)_SCE\\.RDS", replacement="\\1")
print(sce.names)
names(sce.file.list) <- sce.names

## =============================================================================
## Load SCE Objects and Combine
## =============================================================================

sce.list <- list()
for(x in seq_along(sce.names)){
  x.nom <- sce.names[x]
  message("Loading ", x.nom)
  x.sce <- readRDS(sce.file.list[[x.nom]])

  scdblfinder_cols <- grep("^scDblFinder", colnames(rowData(x.sce)), value=TRUE)
  if(length(scdblfinder_cols) > 0){
    message("  Removing ", length(scdblfinder_cols), " scDblFinder metadata columns from gene data")
    rowData(x.sce) <- rowData(x.sce)[, !colnames(rowData(x.sce)) %in% scdblfinder_cols, drop=FALSE]
  }

  colnames(x.sce) <- paste(x.nom, colData(x.sce)$Barcode, sep="_")

  colData(x.sce)$Sample <- x.nom

  if(!is.null(sample_meta) && x.nom %in% rownames(sample_meta)){
    colData(x.sce)$Patient_ID <- sample_meta[x.nom, "Patient_ID"]
    colData(x.sce)$Tumor_Type <- sample_meta[x.nom, "Tumor_Type"]
    colData(x.sce)$Fraction <- sample_meta[x.nom, "Fraction"]
    colData(x.sce)$Replicate <- sample_meta[x.nom, "Replicate"]
    colData(x.sce)$Sequencing_Run <- sample_meta[x.nom, "Sequencing_Run"]
    message("  Patient ", sample_meta[x.nom, "Patient_ID"],
            " ", sample_meta[x.nom, "Fraction"],
            " (replicate ", sample_meta[x.nom, "Replicate"], ")")
  } else {
    colData(x.sce)$Patient_ID <- NA
    colData(x.sce)$Tumor_Type <- NA
    colData(x.sce)$Fraction <- NA
    colData(x.sce)$Replicate <- NA
    colData(x.sce)$Sequencing_Run <- NA
    message("  Metadata not found for ", x.nom)
  }

  sce.list[[x.nom]] <- x.sce
}

message("Combining ", length(sce.list), " SCE files together")
big.sce <- do.call(cbind, sce.list)
counts(big.sce) <- as(counts(big.sce), "dgCMatrix")

message("Combined SCE dimensions: ", nrow(big.sce), " genes x ", ncol(big.sce), " cells")

## =============================================================================
## Verify Metadata Was Successfully Added
## =============================================================================

message("\n=== METADATA VERIFICATION ===")
if(!is.null(sample_meta)){
  message("Cells with Patient_ID: ", sum(!is.na(big.sce$Patient_ID)), " / ", ncol(big.sce))
  message("Cells with Fraction info: ", sum(!is.na(big.sce$Fraction)), " / ", ncol(big.sce))

  if(sum(!is.na(big.sce$Patient_ID)) == 0){
    stop("ERROR: No cells have Patient_ID metadata! Something went wrong during metadata assignment.")
  }

  message("\nCells per patient:")
  patient_table <- table(big.sce$Patient_ID)
  for(p in names(patient_table)){
    message("  ", p, ": ", patient_table[p], " cells")
  }

  message("\nCells per timepoint:")
  fraction_table <- table(big.sce$Fraction)
  for(f in names(fraction_table)){
    message("  ", f, ": ", fraction_table[f], " cells")
  }

  message("\nMetadata successfully added to all cells")
} else {
  warning("No metadata was loaded - cells will only have Sample and Barcode information")
}

## =============================================================================
## Gene Filtering for Normalization
## =============================================================================

n.cells <- ncol(big.sce)
n.genes <- nrow(big.sce)
message(paste0("Computing gene expression sparsity over ", n.cells, " cells"))

keep_genes <- rowMeans(counts(big.sce)) > 0.01
genes <- rownames(big.sce)

message(paste0("Using ", sum(keep_genes), " genes for size factor estimation"))

## =============================================================================
## Size Factor Estimation and Normalization
## =============================================================================

if(isTRUE(opt$breaks)){
  poss.ints <- c(3:50)
  mods <- ncol(big.sce) %% poss.ints
  if(sum(mods == 0) < 1){
    warning("No integer values <50 that creates an integer divisor - setting to 5 chunks")
    n.chunks <- 5
  } else{
    n.chunks <- min(poss.ints[which(mods == 0)])
  }

  max.n <- floor(ncol(big.sce)/n.chunks)
  d <- seq_along(c(1:ncol(big.sce)))
  max.points <- unlist(lapply(split(d, ceiling(d/max.n)), max))

  sce.list <- list()
  for(q in seq_along(max.points)){
    q.max <- max.points[q]
    if(q == 1){
      x.sce <- big.sce[, 1:q.max]
    } else{
      x.sce <- big.sce[, (old.max+1):q.max]
    }
    old.max <- q.max

    cluster.size <- ceiling(ncol(x.sce) * 0.1)
    message(paste0("Chunk ", q, ": Cluster size set to ", cluster.size))

    clusters <- quickCluster(x.sce, min.size=cluster.size,
                             subset.row=keep_genes,
                             method="igraph")
    max.size <- floor(cluster.size/2)

    size.inc <- ceiling(max.size * 0.5)
    message(paste0("Estimating size factors using ", size.inc, " cell increments"))

    x.sce <- computeSumFactors(x.sce,
                                max.cluster.size=max.size,
                                positive=TRUE,
                                subset.row=keep_genes,
                                BPPARAM=MulticoreParam(workers=2),
                                assay.type='counts',
                                clusters=clusters)

    neg.sf <- sum(sizeFactors(x.sce) < 0)
    if(neg.sf > 0){
      message(paste0(neg.sf, " negative size factors estimated - consider higher sparsity threshold"))
    }

    message("Normalising single-cell expression values for chunk ", q)
    x.sce <- logNormCounts(x.sce)
    colData(x.sce)$Batch <- as.character(q)

    batch.file <- gsub(opt$output, pattern="\\.RDS",
                       replacement=paste0("_Batch", q, ".RDS"))
    saveRDS(x.sce, file=batch.file)

    sce.list[[paste0(q)]] <- x.sce

    sink(file="/dev/null")
    rm(list=c("x.sce", "clusters"))
    gc()
    sink(file=NULL)
  }

  big.sce <- do.call(cbind, sce.list)

  message(paste0("Performing multi-batch normalisation across ",
                 length(unique(big.sce$Batch)), " batches"))
  big.sce <- multiBatchNorm(big.sce, batch=big.sce$Batch,
                            assay.type="counts",
                            min.mean=0.01,
                            normalize.all=TRUE,
                            preserve.single=TRUE)
} else {
  cluster.size <- ceiling(ncol(big.sce) * 0.01)
  message(paste0("Cluster size set to ", cluster.size))

  clusters <- quickCluster(big.sce, min.size=cluster.size,
                           subset.row=keep_genes,
                           BPPARAM=MulticoreParam(workers=12),
                           method="igraph")
  max.size <- floor(cluster.size/2)

  size.inc <- ceiling(max.size * 0.5)
  message(paste0("Estimating size factors using ", size.inc, " cell increments"))

  big.sce <- computeSumFactors(big.sce,
                                max.cluster.size=max.size,
                                positive=TRUE,
                                subset.row=keep_genes,
                                assay.type='counts',
                                clusters=clusters)

  neg.sf <- sum(sizeFactors(big.sce) < 0)
  if(neg.sf > 0){
    message(paste0(neg.sf, " negative size factors estimated - consider higher sparsity threshold"))
  }

  message("Normalising single-cell expression values")
  big.sce <- logNormCounts(big.sce)
}

## =============================================================================
## Clean Up and Save
## =============================================================================

big.sce <- big.sce[!is.na(rownames(rowData(big.sce))), ]
rownames(big.sce) <- genes[!is.na(rownames(rowData(big.sce)))]

message("Saving normalized SCE to: ", opt$output)
saveRDS(big.sce, file=opt$output)

message("All done!")
message("Final dimensions: ", nrow(big.sce), " genes x ", ncol(big.sce), " cells")
