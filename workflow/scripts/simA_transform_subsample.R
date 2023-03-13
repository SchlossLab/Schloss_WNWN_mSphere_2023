#!/usr/bin/env Rscript

set.seed(19760620)

# borrowed from L409-416 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate singe otu_count file from
# another otu_count file
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)

args <- commandArgs(trailingOnly = TRUE)

raw_otu_file <- args[1]
smalltrim <- as.numeric(args[2]) #value between 0 and 1
iters <- as.numeric(args[3])

tag <- gsub("0\\.", "", format(smalltrim, nsmall = 2L))
method <- ifelse(iters == 1, "subsample", "rarefaction")

subsample_otu_file <- gsub("RDS",
                           paste0(method, tag, ".RDS"),
                           raw_otu_file)

raw_otu_data <- readRDS(raw_otu_file)

# Set the minimum value as the smallest library quantile, `smalltrim`
samplemin <- quantile(phyloseq::sample_sums(raw_otu_data), smalltrim)

subsample_otu_data <- replicate(
                        iters,
                        phyloseq::rarefy_even_depth(raw_otu_data,
                                                    sample.size = samplemin,
                                                    rngseed = FALSE,
                                                    replace = FALSE,
                                                    trimOTUs = TRUE,
                                                    verbose = FALSE)
                        )

saveRDS(subsample_otu_data, subsample_otu_file)
