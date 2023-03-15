#!/usr/bin/env Rscript

# nr-s1

library(phyloseq)
data("GlobalPatterns")

gp_counts <- data.frame(samples = phyloseq::sample_names(GlobalPatterns),
                            n_seqs = phyloseq::sample_sums(GlobalPatterns),
                            simulation = "GlobalPatterns")
rownames(gp_counts) <- NULL

gp_min <- min(gp_counts$n_seqs)
gp_max <- max(gp_counts$n_seqs)

log_counts <- data.frame(samples = as.character(1:80),
                        n_seqs = floor(exp(seq(log(gp_min),
                                              log(gp_max),
                                              length.out = 80))),
                        simulation = "Log-scaled")

sample_counts <- rbind(gp_counts, log_counts)

write.table(sample_counts,
            file = "data/simulated_nseqs_distros.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
