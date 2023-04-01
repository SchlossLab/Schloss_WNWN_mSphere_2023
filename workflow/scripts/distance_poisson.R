#!/usr/bin/env Rscript

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

edger_to_phyloseq <- function(edger) {

  counts <- as.data.frame(edger$count)
  factors <- edger$samples

  if (!any(is.na(factors$norm.factors))) {

    samples <- data.frame(sample = gsub("\\d*", "", colnames(counts)))
    rownames(samples) <- colnames(counts)

    transformed_otu <- phyloseq(otu_table(counts, taxa_are_rows = TRUE),
                                sample_data(samples),
                                phy_tree(GlobalPatterns))
  } else {
    transformed_otu <- NA
  }

  transformed_otu
}


get_distances <- function(otu_counts) {

  d <- NA

  # get the mean distance from across the elements of the list. this is
  # necessary when using rarefied datasets. in previous experience, the sd
  # is miniscule relative to the mean
  if (is.list(otu_counts)) {

    distances <- lapply(otu_counts, poisson_dist)
    d <- Reduce("+", distances) / length(distances)

  } else {

    d <- poisson_dist(otu_counts)

  }

  return(d)

}

# pds: wnwn didn't use these with the edgeR approaches. not sure why not. will
# use the counts element from the list. suspect won't use / using wrong with
# these distance calculators. if any of the factor values are NA/NaN will return
# an NA

mean_distance <- NA
transformed_otu_list <- readRDS(transformed_otu_file)

if (grepl("upperquartile|RLE|TMM", transformed_otu_file)) {

  library(edgeR)
  data(GlobalPatterns)

  transformed_otu_list <- lapply(transformed_otu_list, edger_to_phyloseq)

}

if (any(is.na(transformed_otu_list))) {
  mean_distance <- rep(NA, length(transformed_otu_list))

} else {
  mean_distance <- lapply(transformed_otu_list, get_distances)
}

saveRDS(mean_distance, distance_file)
