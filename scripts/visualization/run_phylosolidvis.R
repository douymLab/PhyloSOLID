#!/usr/bin/env Rscript
# Run PhyloSOLIDvis::run_all on a PhyloSOLID phylo/ directory.
# Args: inputpath outputpath [annotation_file]
# Exit 2 if PhyloSOLIDvis is not installed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: run_phylosolidvis.R <inputpath> <outputpath> [annotation_file]\n")
  quit(status = 1)
}

inputpath <- args[[1]]
outputpath <- args[[2]]
annotation_file <- if (length(args) >= 3) args[[3]] else ""

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  cat("ERROR: ggplot2 is not installed in this R environment\n")
  quit(status = 1)
}
if (!requireNamespace("PhyloSOLIDvis", quietly = TRUE)) {
  cat("PhyloSOLIDvis is not installed in this R environment\n")
  quit(status = 2)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(PhyloSOLIDvis)
})

dir.create(outputpath, recursive = TRUE, showWarnings = FALSE)

kwargs <- list(
  inputpath = inputpath,
  outputpath = outputpath,
  verbose = TRUE
)

use_annotation <- nzchar(annotation_file) &&
  !annotation_file %in% c("None", "none", "NA") &&
  file.exists(annotation_file)
if (use_annotation) {
  kwargs$annotation_file <- annotation_file
  cat("Using annotation file:", annotation_file, "\n")
} else {
  cat("No annotation file; drawing tree/heatmap without annotation layers\n")
}

do.call(run_all, kwargs)
cat("PhyloSOLIDvis finished:", outputpath, "\n")
