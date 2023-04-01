#!/usr/bin/env Rscript

# borrowed from L839-909 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate file with list of distances
# from a phyloseq formatted otu_count file
#
# be sure to first run: conda activate nr-s1

library(cluster)

args <- commandArgs(trailingOnly = TRUE)

distance_file <- args[1]

pretty_distance_file <- gsub(".*/(.*)\\.RDS", "\\1", distance_file)
output_file <- gsub("RDS", "clusters.tsv", distance_file)

distances <- readRDS(distance_file)
samples_per_type <- 40L  # need to hardcode

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

  if (n_not_finite != 0) {

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

# Evaluate the results of each clustering.
binaryclusteval <- function(x, tot) {

# Assumes argument, x, is a named cluster vector with two components,
# and that the "true" clustering can be inferred from the name of each
# element, namely, the non-digit, non-delimiter word in each name.

 clres <- c(fracCorrectPred = NA, fracCorrect = NA)

  if (!is.na(x[1])) {

    truth <- as(factor(gsub("[[:digit:]]{1,}", "", names(x))), "numeric")
    xinv <- ifelse(x == 1, 2, 1)
    correct <- max(sum(x == truth), sum(xinv == truth))
    clres <- c(fracCorrectPred = (correct / length(x)),
              fracCorrect = (correct / tot))

  }

  return(clres)
}


cluster_distance_matrix <- function(distance_matrix) {
  # this is needed for poisson distances
  distance_matrix <- fixbad_dists(distance_matrix)


  # get this to generate all three types of clusters in one go...

  # PDS: the following comment is wrong- pam comes from cluster package
  # Binary clustering using `limma::pam` (partitioning around medoids) with
  # `k==2`, and alternatively `stats::kmeans` and `stats::hclust`.

  pam <- NA
  kmeans <- NA
  hclust <- NA

  if (!is.na(distance_matrix[1])) {

    pam <- cluster::pam(distance_matrix, k = 2, cluster.only = TRUE)

    cmdscale_results <- stats::cmdscale(distance_matrix, k = 2)

    if (ncol(cmdscale_results) != 0) {

      kmeans <- stats::kmeans(cmdscale_results,
                                centers = 2)$cluster
    }

    hclust <- stats::cutree(stats::hclust(distance_matrix), k = 2)

  }

  pam_clres <- binaryclusteval(pam, 2 * samples_per_type)
  kmeans_clres <- binaryclusteval(kmeans, 2 * samples_per_type)
  hclust_clres <- binaryclusteval(hclust, 2 * samples_per_type)

  clres <- rbind(pam_clres, kmeans_clres, hclust_clres)

  clres <- data.frame(method = c("pam", "kmeans", "hclust"),
                      clres,
                      row.names = NULL)

  return(clres)

}

output_data <- do.call(rbind, lapply(distances, cluster_distance_matrix))
output_data$conditions <- pretty_distance_file
output_data$reps <- rep(seq_along(distances), each = 3)

write.table(output_data,
            file = output_file,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
