#!/usr/bin/env Rscript

# be sure to first run: conda activate nr-modern

library(tidyverse)

read_files <- function(x) {

      readr::read_tsv(x,
               col_types = readr::cols(method = readr::col_character(),
                                       .default = readr::col_double())
      ) %>%
      dplyr::mutate(conditions = x, .before = tidyselect::everything())

}

tsv_files <- list.files(path = "data/sim_a",
                        pattern = "clusters.tsv",
                        full.names = TRUE)

map_dfr(tsv_files, read_files) %>%
  rename(subset = fracCorrectPred,
        all = fracCorrect,
        cluster = method) %>%
  mutate(conditions = str_replace(conditions, "data/sim_a/", ""),
        model = str_replace(conditions,
                             "\\.n?o?filter.*", ""),
         processing = str_replace(conditions,
                                  ".*\\.(n?o?filter.*).clusters.tsv",
                                  "\\1")) %>%
  select(-conditions) %>%
  separate(model,
           into = c("fraction", "n_seqs", "rep"),
           sep = "_") %>%
  separate(processing,
           into = c("filter", "transform", "distance"),
           sep = "\\.") %>%
  select(fraction, n_seqs, rep, filter, transform, distance,
         cluster, subset, all) %>%
  write_tsv("data/simulation_cluster_accuracy.tsv")
