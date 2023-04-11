#!/usr/bin/env Rscript

# nr-s1

library(phyloseq)

model <- c("gp", "log")
distribution <- c("random", "skew")
effect_size <- c(1, 1.15, 1.25, 1.5, 1.75, 2, 2.5, 3.5)
filter <- c("filter", "nofilter")
n_seqs <- c(1000, 2000, 5000, 10000, 50000)

paths <- expand.grid("data", model, distribution)
paths <- paste(paths$Var1, paths$Var2, paths$Var3, sep = "/")

conditions <- expand.grid(path = paths,
                          effect_size = effect_size,
                          filter = filter,
                          n_seqs = n_seqs)

conditions$rds <- paste0(conditions$path, "/", conditions$effect_size, "_",
                    conditions$n_seqs, ".", conditions$filter, ".deseq.RDS")


get_neg_fraction_each_rep <- function(ps_object) {

  counts <- phyloseq::otu_table(ps_object)
  (sum(counts < 0)) / length(counts)

}

is_neg_matrix <- function(ps_object) {

  any(phyloseq::otu_table(ps_object) < 0)

}

driver <- function(rds) {

  phyloseq_objects <- readRDS(rds)
  fractions <- sapply(phyloseq_objects, get_neg_fraction_each_rep)
  negatives <- mean(sapply(phyloseq_objects, is_neg_matrix))

  list(rds = rds,
    mean = mean(fractions),
    lci = quantile(fractions, prob = 0.025),
    uci = quantile(fractions, prob = 0.975),
    has_negative = negatives
  )

}

results <- do.call(rbind, lapply(conditions$rds, driver))

data <- merge(conditions, results, by = "rds")
data$model <- gsub("data/(.*)/.*", "\\1", data$path)
data$distribution <- gsub("data/.*/(.*)", "\\1", data$path)
data <- data[, !(colnames(data) %in% c("path", "rds"))]
data$mean <- as.numeric(data$mean)
data$lci <- as.numeric(data$lci)
data$uci <- as.numeric(data$uci)
data$has_negative <- as.numeric(data$has_negative)

write.table(data, file = "data/deseq_negative_values.tsv.gz", sep = "\t",
            quote = FALSE, row.names = FALSE)
