#!/usr/bin/env Rscript

library(phyloseq)
data("GlobalPatterns")

sample_counts <- data.frame(samples = phyloseq::sample_names(GlobalPatterns),
                            n_seqs = phyloseq::sample_sums(GlobalPatterns))

write.table(sample_counts,
            file = "data/globalpatterns_distribution.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
