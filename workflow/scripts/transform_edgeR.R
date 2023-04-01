#!/usr/bin/env Rscript

set.seed(19760620)

# borrowed from L425-457 and L572-577 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate singe otu_count file from
# another otu_count file
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

raw_otu_file <- args[1]
method <- args[2] #c("upperquartile", "TMM", "RLE")

edger_otu_file <- gsub("RDS", paste0(method, ".RDS"), raw_otu_file)

raw_otu_data <- readRDS(raw_otu_file)

edger_data <- function(raw_data, method) {

  if (!taxa_are_rows(raw_data)) {
    raw_data <- t(raw_data)
  }

  x <- as(otu_table(raw_data), "matrix")

  # See if adding a single observation, 1, 
  # everywhere (so not zeros) prevents errors
  # without needing to borrow and modify 
  # calcNormFactors (and its dependent functions)
  # It did. This fixed all problems. 
  # Can the 1 be reduced to something smaller and still work?
  x <- x + 1

  # Now turn into a DGEList
  y <- edgeR::DGEList(counts = x, remove.zeros = TRUE)

  # Perform edgeR-encoded normalization, using the specified method (...)
  z <- edgeR::calcNormFactors(y, method = method)

  return(z)
}

edger_otu_data <- lapply(raw_otu_data, edger_data, method = method)

saveRDS(edger_otu_data, edger_otu_file)
