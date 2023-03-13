#!/usr/bin/env Rscript

# be sure to first run: conda activate nr-modern

library(tidyverse)

not_all_na <- function(x) any(!is.na(x))

args <- commandArgs(trailingOnly = TRUE)

output_filename <- args[1]
#output_filename <- "data/simulation_adonis.tsv.gz"

pattern <- str_replace(output_filename, ".*_(.*\\.tsv).*", "pool\\.\\*\\1")

input_filenames <- list.files(path = "data/sim_a",
                           pattern = pattern,
                           full.names = TRUE)


read_delim(input_filenames, show_col_types = FALSE) %>%
  mutate(conditions = str_replace(conditions, "data/sim_a/", ""),
         model = str_replace(conditions, "\\.n?o?filter.*", ""),
         processing = str_replace(conditions, ".*\\.(n?o?filter.*)", "\\1")) %>%
  select(-conditions) %>%
  separate(model,
           into = c("fraction", "n_seqs", "rep"),
           sep = "_") %>%
  separate(processing,
           into = c("filter", "transform", "distance"),
           sep = "\\.", fill = "right") %>%
  select(fraction, n_seqs, rep, filter, transform, distance, everything()) %>%
  select(where(not_all_na)) %>%
  write_tsv(output_filename)
