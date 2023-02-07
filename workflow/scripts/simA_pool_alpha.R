#!/usr/bin/env Rscript

# be sure to first run: conda activate nr-modern

library(tidyverse)

alpha_files <- list.files(path = "data/sim_a",
                          pattern = "alpha.tsv",
                          full.names = TRUE)

read_alpha <- function(x) {
  read_delim(x, delim = " ",
             col_types = cols(conditions = col_character(),
                              .default = col_double())) %>%
    mutate(filter_transform = str_replace(conditions,
                                          ".*\\.(n?o?filter.*)",
                                          "\\1"),
          conditions = str_replace(conditions,
                                   "\\.n?o?filter.*",
                                   "")) %>%
    separate(filter_transform, into = c("filter", "transform"), sep = "\\.") %>%
    separate(conditions, into = c("fraction", "n_seqs", "rep"), sep = "_") %>%
    mutate(n_seqs = as.numeric(n_seqs),
           rep = as.numeric(rep)) %>%
    select(fraction, n_seqs, rep, filter, transform, everything())
}

map_dfr(alpha_files, read_alpha) %>%
  write_tsv("data/simulation_alpha.tsv")
