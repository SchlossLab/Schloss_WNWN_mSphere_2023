#!/usr/bin/env Rscript

set.seed(20140206)

library(phyloseq)
data("GlobalPatterns")

sampletypes <- c("Feces", "Ocean")

# Number of OTUs in simulation.
n_otus <- 2000L

# Minimum number of reads to consider an OTU "observed" in a sample
minobs <- 1L

# Samples per simulation
samples_per_type <- 40L

sampsums <- phyloseq::sample_sums(GlobalPatterns)
keepsamples <- phyloseq::sample_data(GlobalPatterns)$SampleType %in% sampletypes
template <- phyloseq::prune_samples(keepsamples, GlobalPatterns)

samobs <- apply(phyloseq::otu_table(template),
                1,
                function(x, m) sum(x > m), m = minobs)

otudf <- data.frame(prev = samobs, sums = taxa_sums(template))
otudf <- otudf[order(-otudf$prev, -otudf$sums), ]

template <- prune_taxa(rownames(otudf)[1:n_otus], template)
template1 <- phyloseq::subset_samples(template, SampleType == sampletypes[1])
template1 <- phyloseq::merge_samples(template1, "SampleType")

template2 <- phyloseq::subset_samples(template, SampleType == sampletypes[2])
template2 <- phyloseq::merge_samples(template2, "SampleType")

distros <- otu_table(rbind(otu_table(template1),
                          otu_table(template2)),
                    taxa_are_rows = FALSE)

min_size <- min(apply(distros, 1, sum))

get_shannon <- function(x) {

  x <- as.numeric(x)
  x <- x[x!=0]
  ra <- x / sum(x)
  -sum(ra * log(ra))

}

subsample_alpha_beta <- function(size) {

  subsample <- phyloseq::rarefy_even_depth(distros,
                                          sample.size = size,
                                          rngseed = FALSE,
                                          replace = FALSE,
                                          trimOTUs = TRUE,
                                          verbose = FALSE)

  richness <- apply(subsample > 0, 1, sum)

  shannon <- c(hocean = get_shannon(subsample["Ocean", ]),
              hfeces = get_shannon(subsample["Feces", ])
              )

  shared <- sum(subsample["Ocean", ] > 0 & subsample["Feces", ] > 0)
  jaccard <- 1 - shared / (sum(richness) - shared)

  bray <- as.numeric(phyloseq::distance(subsample, method = "bray"))

  return(data.frame(depth = size,
                    rocean = richness[["Ocean"]], rfeces = richness[["Feces"]],
                    hocean = shannon[["hocean"]], hfeces = shannon[["hfeces"]],
                    jaccard = jaccard, bray = bray))
}

get_subsamples <- function(size){

  subsamples <- replicate(100,
                    subsample_alpha_beta(size),
                    simplify = FALSE)
  subsamples <- do.call(rbind, subsamples)

  return(subsamples)
}

depth_subsamples <- sapply(c(1000, 2000, 5000, 10000, 50000, min_size),
      get_subsamples, simplify = FALSE)
depth_subsamples <- do.call(rbind, depth_subsamples)



richness <- apply(distros > 0, 1, sum)

shannon <- c(hocean = get_shannon(distros["Ocean", ]),
            hfeces = get_shannon(distros["Feces", ])
            )

shared <- sum(distros["Ocean", ] > 0 & distros["Feces", ] > 0)
jaccard <- 1 - (shared / (sum(richness) - shared))

bray <- as.numeric(phyloseq::distance(distros, method = "bray"))

full_dataset <- data.frame(
                  depth = 0,
                  rocean = richness[["Ocean"]], rfeces = richness[["Feces"]],
                  hocean = shannon[["hocean"]], hfeces = shannon[["hfeces"]],
                  jaccard = jaccard, bray = bray
                )

depth_subsamples <- rbind(full_dataset, depth_subsamples)

write.table(depth_subsamples, file = "data/rarefy_ocean_feces.tsv.gz",
            sep = "\t", row.names = FALSE, quote = FALSE)
