#!/usr/bin/env Rscript

library(converTree)
library(ape)
library(treeio)
library(purrr)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript convert_tree.R <input_file> <output_file>")
}

tree_file <- args[1]
output_file <- args[2]

# Main conversion function
convert_tree <- function(input_file, output_file) {
  # Read data and convert
  treedata_truth <- cf2treedata(input_file)
  phylo_tree <- treedata_truth %>% ape::as.phylo()
  
  # Index-label helper
  index_label <- function(node) {
    treedata_truth[treedata_truth$node == node, ]$label
  }
  
  # Update tip labels
  phylo_tree[["tip.label"]] <- map(phylo_tree[["tip.label"]], index_label) %>% unlist()
  
  # Generate Newick-format tree
  nwk_tree <- treeio::write.tree(phylo_tree)
  
  # Remove internal-node numbers
  nwk_tree <- gsub(")[0-9]+", ")", nwk_tree)
  
  # Save to file
  writeLines(nwk_tree, output_file)
  
  cat("Conversion completed. Results saved to:", output_file, "\n")
  cat("Generated tree:", nwk_tree, "\n")
}

# Run conversion
convert_tree(tree_file, output_file)