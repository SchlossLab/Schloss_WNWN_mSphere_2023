#!/usr/bin/env Rscript --vanilla

# borrowed from L366-386 and XXXXX of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate single otu_count file from
# a non-filtered version
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

no_filter_file <- args[1]
filter_file <- gsub("nofilter", "filter", no_filter_file)

no_filter <- readRDS(no_filter_file)

#remove OTUs that have total abundance across samples less than 3
filter <- phyloseq::prune_taxa(phyloseq::taxa_sums(no_filter) > 2.5,
                                  no_filter)

#remove samples that have total abundance less than 3(?)
#doesn't appear to do anything...
filter <- phyloseq::prune_samples(phyloseq::sample_sums(filter) > 2.5,
                                     filter)

# Remove OTUs not appearing in at least 3 samples
if (taxa_are_rows(filter)) { #should be false
  y <- otu_table(filter)
  wh1 <- apply(aaply(y, 1, function(x){x >= 1}), MARGIN=1, sum) >= 3
} else {
  y <- as(t(otu_table(filter)), "matrix")
  wh1 <- apply(aaply(y, 1, function(x){x >= 1}), MARGIN=1, sum) >= 3
}
filter <- prune_taxa(wh1, filter)

saveRDS(filter, filter_file)
