#!/usr/bin/env Rscript --vanilla

# borrowed from L754-837 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate file with list of distances
# from a phyloseq formatted otu_count file
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)

args <- commandArgs(trailingOnly = TRUE)

transformed_otu_file <- args[1]
pretty_method <- args[2]

stopifnot(pretty_method %in% c("bray", "euclidean", "wunifrac", "uunifrac"))

distance_file <- gsub("RDS",
                      paste0(pretty_method, ".RDS"),
                      transformed_otu_file)

method <- ifelse(pretty_method == "uunifrac", "unifrac", pretty_method)

# wnwn didn't use these with the edgeR approaches. not sure why not. will use
# the counts element from the list. suspect won't use / using wrong with these
# distance calculators. if any of the factor values are NA/NaN will return an
# NA

mean_distance <- NA
skip <- FALSE
transformed_otu <- NULL

if (grepl("upperquartile|RLE|TMM", transformed_otu_file)) {

  library(edgeR)
  data(GlobalPatterns)

  edger_data <- readRDS(transformed_otu_file)
  counts <- as.data.frame(edger_data$count)
  factors <- edger_data$samples

  skip <- TRUE

  if (!any(is.na(factors$norm.factors))) {
    skip <- FALSE

    samples <- data.frame(sample = gsub("\\d*", "", colnames(counts)))
    rownames(samples) <- colnames(counts)

    transformed_otu <- list(phyloseq(otu_table(counts, taxa_are_rows = TRUE),
                                sample_data(samples),
                                phy_tree(GlobalPatterns)))
  }
} else {
  transformed_otu <- readRDS(transformed_otu_file)

  if (!is.list(transformed_otu)) {
    transformed_otu <- list(transformed_otu)
  }
}

if (!skip) {
  # get the mean distance from across the elements of the list. this is
  # necessary when using rarefied datasets. in previous experience, the sd
  # is miniscule relative to the mean
  distances <- lapply(transformed_otu, phyloseq::distance, method = method)
  mean_distance <- Reduce("+", distances) / length(distances)
}

saveRDS(mean_distance, distance_file)
