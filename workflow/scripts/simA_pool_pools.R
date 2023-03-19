#!/usr/bin/env Rscript

# be sure to first run: conda activate nr-modern

library(tidyverse)

not_all_na <- function(x) any(!is.na(x))

read_pool_files <- function(input_filenames) {

  readr::read_delim(input_filenames, show_col_types = FALSE) %>%
    dplyr::mutate(
         simulation = if_else(str_detect(conditions, "^s?l"), "sim_log", "sim_a"),
         conditions = str_replace(conditions, "l?(\\d)", "\\1"),
         model = str_replace(conditions, "\\.n?o?filter.*", ""),
         processing = str_replace(conditions, ".*\\.(n?o?filter.*)", "\\1")) %>%
   select(-conditions) %>%
   separate(model,
            into = c("fraction", "n_seqs", "rep"),
            sep = "_") %>%
   separate(processing,
            into = c("filter", "transform", "distance"),
            sep = "\\.", fill = "right") %>%
   select(simulation, fraction, n_seqs, rep, filter,
      transform, distance, everything()) %>%
   select(where(not_all_na))

}

args <- commandArgs(trailingOnly = TRUE)

output_filename <- args[1]
#output_filename <- "data/simulation_adonis.tsv.gz"

pattern <- str_replace(output_filename, ".*_(.*\\.tsv).*", "pool\\.\\*\\1")

input_filenames_a <- list.files(path = "data/sim_a",
                           pattern = pattern,
                           full.names = TRUE)

input_filenames_log <- list.files(path = "data/sim_log",
                           pattern = pattern,
                           full.names = TRUE)

sim_a_pools <- read_pool_files(input_filenames_a)
sim_log_pools <- read_pool_files(input_filenames_log)

bind_rows(sim_a_pools, sim_log_pools) %>%
  write_tsv(output_filename)
