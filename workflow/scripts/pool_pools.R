#!/usr/bin/env Rscript

# be sure to first run: conda activate nr-modern

library(tidyverse)

not_all_na <- function(x) any(!is.na(x))

read_pool_files <- function(input_filenames) {

  readr::read_delim(input_filenames, show_col_types = FALSE, id = "path") %>% 
    dplyr::mutate(
         conditions = str_replace(conditions, "l?(\\d)", "\\1"),
         model = str_replace(conditions, "\\.n?o?filter.*", ""),
         processing = str_replace(conditions, ".*\\.(n?o?filter.*)", "\\1")) %>%
   select(-conditions) %>%
   separate(model,
            into = c("fraction", "n_seqs"),
            sep = "_") %>%
   separate(processing,
            into = c("filter", "transform", "distance"),
            sep = "\\.", fill = "right") %>% 
   separate(path,
            into = c("parent", "model", "distribution", "file"),
            sep = "\\/") %>% 
   select(-file, -parent) %>%
   select(model, distribution, fraction, n_seqs, reps, filter,
      transform, distance, everything()) %>%
   select(where(not_all_na))

}

args <- commandArgs(trailingOnly = TRUE)

output_filename <- args[1]
#output_filename <- "data/simulation_clusters.tsv.gz"

pattern <- str_replace(output_filename, ".*_(.*\\.tsv).*", "\\1")

input_filenames <- list.files(path = "data",
                              pattern = pattern,
                              recursive = TRUE,
                              full.names = TRUE)

read_pool_files(input_filenames) %>%
  write_tsv(output_filename)
