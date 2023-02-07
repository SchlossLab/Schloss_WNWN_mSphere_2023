#!/usr/bin/env Rscript

# borrowed from L392-401 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate singe otu_count file from
# another otu_count file
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)

args <- commandArgs(trailingOnly = TRUE)

raw_otu_file <- args[1]
normalize_otu_file <- gsub("RDS", "proportion.RDS", raw_otu_file)

raw_otu_data <- readRDS(raw_otu_file)

# Normalize total sequences represented
# PDS: each sample takes on the total abundance of the most abundant sample
#      in the dataset; the correlation between the values from this approach
#      and normalizing to the smallest sample or to a total of 1 is 1.0
normf <- function(x, tot = max(phyloseq::sample_sums(raw_otu_data))) {

  tot * x / sum(x)

}

normalize_otu_data <- phyloseq::transform_sample_counts(raw_otu_data, normf)

saveRDS(normalize_otu_data, normalize_otu_file)
