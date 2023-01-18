#!/usr/bin/env Rscript --vanilla

# be sure to first run: conda activate nr-modern

library(tidyverse)

adonis_files <- list.files(path = "data/sim_a",
                           pattern = "adonis.tsv",
                           full.names = TRUE)

read_adonis <- function(x) {
  
    read_tsv(x,
             col_types = cols(conditions = col_character(),
                              .default = col_double())) %>%
    mutate(filter_transform = str_replace(conditions,
                                          ".*\\.(n?o?filter.*)",
                                          "\\1"),
          conditions = str_replace(conditions,
                                   "(.*)\\.n?o?filter.*",
                                   "\\1")) %>%
    separate(filter_transform,
             into = c("filter", "transform", "distance"),
             sep = "\\.") %>%
    separate(conditions, into = c("fraction", "n_seqs", "rep"), sep = "_") %>%
    mutate(n_seqs = as.numeric(n_seqs),
           rep = as.numeric(rep)) %>%
    select(fraction, n_seqs, rep, filter, transform, distance, r_squared, p_value)

}

map_dfr(adonis_files, read_adonis) %>%
  write_tsv("data/simulation_adonis.tsv")
