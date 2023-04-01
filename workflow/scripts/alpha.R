#!/usr/bin/env Rscript

# borrowed from L839-909 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate file with list of distances
# from a phyloseq formatted otu_count file

# be sure to first run: conda activate nr-s1

library(phyloseq)

args <- commandArgs(trailingOnly = TRUE)

phyloseq_list_file <- args[1]
output_file <- gsub("RDS", "alpha.tsv", phyloseq_list_file)

get_sobs <- function(x) {

  sobs <- NA

  if (all(x >= 0)) {
    sobs <- sum(x > 0)
  }

  return(sobs)
}

get_shannon <- function(x) {

  x <- x[x != 0]

  shannon <- NA

  if (all(x > 0)){
    n <- sum(x)

    p <- x / n

    shannon <- -1 * sum(p * log(p))
  }

  return(shannon)
}

get_alpha_subsamples <- function(x, f) {

  if (class(x) == "DGEList") {
    otu_count_table <- t(x$counts)
  } else {
    otu_count_table <- phyloseq::otu_table(x)

    if (phyloseq::taxa_are_rows(x)) {
      otu_count_table <- t(otu_count_table)
    }
  }

  f_outupt <- apply(otu_count_table, 1, f)

  return(f_outupt)

}

run_test <- function(otu_data) {

  if (is.list(otu_data) & class(otu_data) != "DGEList") {

    sobs_list <- lapply(otu_data, get_alpha_subsamples, get_sobs)
    sobs_df <- as.data.frame(do.call(rbind, sobs_list))
    sobs <- apply(sobs_df, 2, mean)

    shannon_list <- lapply(otu_data, get_alpha_subsamples, get_shannon)
    shannon_df <- as.data.frame(do.call(rbind, shannon_list))
    shannon <- apply(shannon_df, 2, mean)

  } else {

    sobs <- get_alpha_subsamples(otu_data, get_sobs)
    shannon <- get_alpha_subsamples(otu_data, get_shannon)

  }

  treatment_group <- gsub("\\d*", "", names(shannon))

  sobs_p <- NA
  sobs_feces <- NA
  sobs_ocean <- NA

  if (!any(is.na(sobs))) {
    sobs_p <- wilcox.test(sobs ~ treatment_group, exact = FALSE)$p.value
    sobs_feces <- median(sobs[treatment_group == "Feces"])
    sobs_ocean <- median(sobs[treatment_group == "Ocean"])
  }

  shannon_p <- NA
  shannon_feces <- NA
  shannon_ocean <- NA

  if (!any(is.na(shannon))) {
    shannon_p <- wilcox.test(shannon ~ treatment_group, exact = FALSE)$p.value
    shannon_feces <- median(shannon[treatment_group == "Feces"])
    shannon_ocean <- median(shannon[treatment_group == "Ocean"])
  }

  data.frame(sobs_ocean = sobs_ocean,
            sobs_feces = sobs_feces,
            sobs_p = ifelse(is.na(sobs_p), NA, sobs_p),
            shannon_ocean = shannon_ocean,
            shannon_feces = shannon_feces,
            shannon_p = ifelse(is.na(shannon_p), NA, shannon_p)
          )
}

phyloseq_list <- readRDS(phyloseq_list_file)

output_data <- do.call(rbind, lapply(phyloseq_list, run_test))
output_data$conditions <- gsub("data/.*/(.*)\\.RDS", "\\1", phyloseq_list_file)
output_data$reps <- seq_along(phyloseq_list)

write.table(output_data,
            file = output_file,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
