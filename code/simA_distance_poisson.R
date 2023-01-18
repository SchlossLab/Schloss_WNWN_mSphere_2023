#!/usr/bin/env Rscript --vanilla

# borrowed from L730-739 and L800-806 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate file with list of distances
# from a phyloseq formatted otu_count file
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)
library(PoiClaClu)

args <- commandArgs(trailingOnly = TRUE)

transformed_otu_file <- args[1]
distance_file <- gsub("RDS", "poisson.RDS", transformed_otu_file)


poisson_dist <- function(physeq, ...){
  # PoiClaClu expects samples as rows and taxa/variables as columns
  if (phyloseq::taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  x <- as(phyloseq::otu_table(physeq), "matrix")
  dd <- PoiClaClu::PoissonDistance(x, ...)$dd
  attributes(dd)$Labels <- phyloseq::sample_names(physeq)

  return(dd)
}

# wnwn didn't use these with the edgeR approaches. not sure why not. will use
# the counts element from the list. suspect won't use / using wrong with these
# distance calculators. if any of the factor values are NA/NaN will return an
# NA

mean_distance <- NA
skip <- FALSE
transformed_otu <- NULL

if (grepl("upperquartile|RLE|TMM", transformed_otu_file)) {

  library(edgeR)

  edger_data <- readRDS(transformed_otu_file)
  counts <- as.data.frame(edger_data$count)
  factors <- edger_data$samples

  skip <- TRUE

  if (!any(is.na(factors$norm.factors))) {
    skip <- FALSE

    samples <- data.frame(sample = gsub("\\d*", "", colnames(counts)))
    rownames(samples) <- colnames(counts)

    transformed_otu <- list(phyloseq(otu_table(counts, taxa_are_rows = TRUE),
                                sample_data(samples)))
  }
} else {
  transformed_otu <- readRDS(transformed_otu_file)

  transformed_otu <- ifelse(is.list(transformed_otu),
                          transformed_otu,
                          list(transformed_otu))
}

if (!skip) {
  # get the mean distance from across the elements of the list. this is
  # necessary when using rarefied datasets. in previous experience, the sd
  # is miniscule relative to the mean
  distances <- lapply(transformed_otu, poisson_dist)
  mean_distance <- Reduce("+", distances) / length(distances)
}

saveRDS(mean_distance, distance_file)
