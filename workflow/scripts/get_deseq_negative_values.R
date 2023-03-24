#!/usr/bin/env Rscript

# nr-s1

library(phyloseq)

fraction <- c(1, 1.15, 1.25, 1.5, 1.75, 2, 2.5, 3.5)
filter <- c("filter", "nofilter")
#sim <- c("sim_a", "sim_log")
n_seqs <- c(1000, 2000, 5000, 10000, 50000)


conditions <- expand.grid(fraction = fraction,
                          filter = filter,
                          n_seqs = n_seqs,
                          rep = 1:100)

conditions$sim_a <- paste0("data/sim_a/", conditions$fraction, "_",
                          conditions$n_seqs, "_", conditions$rep, ".",
                          conditions$filter, ".deseq.RDS")
conditions$sim_log <- paste0("data/sim_log/l",
                          conditions$fraction, "_", conditions$n_seqs, "_",
                          conditions$rep, ".", conditions$filter, ".deseq.RDS")


get_neg_fraction <- function(rds) {
  counts <- phyloseq::otu_table(readRDS(rds))
  (sum(counts < 0)) / length(counts)
}

is_neg_matrix <- function(rds) {
  any(phyloseq::otu_table(readRDS(rds)) < 0)
}

conditions$sim_a_fraction <- sapply(conditions$sim_a, get_neg_fraction)
conditions$sim_log_fraction <- sapply(conditions$sim_log, get_neg_fraction)

conditions$sim_a_negative <- sapply(conditions$sim_a, is_neg_matrix)
conditions$sim_log_negative <- sapply(conditions$sim_log, is_neg_matrix)

conditions <- conditions[, !names(conditions) %in% c("sim_a", "sim_log")]

write.table(conditions, "data/deseq_negative_values.tsv.gz", sep = "\t",
            quote = FALSE, row.names = FALSE)
