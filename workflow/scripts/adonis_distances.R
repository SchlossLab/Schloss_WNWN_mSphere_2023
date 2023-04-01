#!/usr/bin/env Rscript

library(vegan)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

## Borrowed from code/simA_cluster_distances.R...
# NA distances. So far, it looks like it is very rare for any distance
# matrices to have `NA` values, but if it does happen it makes `pam` fail. The
# reason for this happening is unclear right now, but appears constrained to
# the `PoissonDistance` results. My first guess is that these simulated samples
# have nothing in common, and so `PoissonDistance` gives them a value of `NA`
# rather than infinity (or 1, if it is normalized distance). Instead of removing
# these distance matrices (which is one option for this simulation), since this
# is a rare event, I will instead give these entries the max value that
# otherwise occurred in the distance matrix.

# PDS: also double check for NA distances from when upperquartile, TMM, and RLE
#      produced non-sensical scaling factors - nothing had non finite values

# Define function to replace infinite or NA values in distance matrices
# with 1.25x the max value - rewrote for speed and to fix infinite values (there
# weren't actually any infinite values that I could find)

fixbad_dists <- function(dd) {

  not_finite <- !is.finite(dd)
  n_not_finite <- sum(not_finite)
  n_negative <- sum(dd < 0, na.rm = TRUE)

  if (n_negative != 0) {
    dd <- NA
  } else if (n_not_finite != 0) {

    if (n_not_finite == length(not_finite)) {

      dd <- NA

    } else if (n_not_finite != 0) {

      mx <- max(dd[!not_finite])
      mx <- 1.25 * mx
      dd[not_finite] <- mx
    }

  }

  return(dd)
}


distance_file <- args[1]
conditions <- str_replace(distance_file, ".*\\/(.*)\\.RDS", "\\1")
output_file <- str_replace(distance_file, "RDS", "adonis.tsv")

distances <- readRDS(distance_file)

run_adonis <- function(distance_matrix) {

  # this is needed for poisson distances
  distance_matrix <- fixbad_dists(distance_matrix)

  test_result <- tibble(r_squared = NA,
                        p_value = NA)

  if (!is.na(distance_matrix[1])) {
    sample_ids <- labels(distance_matrix)
    group_ids <- str_replace(sample_ids, "\\d*", "")

    test <- adonis2(as.dist(distance_matrix) ~ group_ids)

    r_squared <- test[["R2"]][1]
    p_value <- test[["Pr(>F)"]][1]

    test_result <- tibble(r_squared = r_squared,
                          p_value = p_value)

  }

  return(test_result)
}

do.call(rbind, lapply(distances, run_adonis)) %>%
  mutate(conditions = conditions,
        reps = seq_along(distances)) %>%
  write_tsv(output_file)
