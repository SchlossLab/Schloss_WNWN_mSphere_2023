#!/usr/bin/env Rscript

# borrowed from L57-346 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate singe otu_count file from
# user specified values of mixfacs, n_l, and rep
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)

args <- commandArgs(trailingOnly = TRUE)


# Effect Size. The artificial mix fraction.
mixfac <- as.numeric(args[1]) #c(1, 1.15, 1.25, 1.5, 1.75, 2, 2.5, 3.5)

# Mean sampling depth
n_l <- as.numeric(args[2]) #c(1000, 2000, 5000, 10000, 50000)

# Vector of the replicate numbers to repeat for
# each comb of simulation parameters (n_l, etc)
rep <- as.numeric(args[3]) #1:5

rel_path <- "data/sim_a/"

dir.create(rel_path, showWarnings = FALSE)

file_name <- paste0(rel_path,
                    paste(mixfac, n_l, rep, sep = "_"),
                    ".nofilter.RDS")

set.seed(20140206 + rep)
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


# Define function for mixing
mix_microbes <- function(template1, template2, unmixfac) {

  # Check that the number of taxa are equal
  if (!identical(phyloseq::ntaxa(template1), ntaxa(template2))) {
    stop("number of OTUs b/w template1 and template2 must be equal")
  }

  # Expects templates to be a 1-sample dataset.
  if (!identical(phyloseq::nsamples(template1), 1L)) {
    stop("`template1` should have only 1 sample")
  }
  if (!identical(phyloseq::nsamples(template2), 1L)){
    stop("`template2` should have only 1 sample")
  }

  # Enforce taxa_are_rows to FALSE
  if (phyloseq::taxa_are_rows(template1)) template1 <- t(template1)
  if (phyloseq::taxa_are_rows(template2)) template2 <- t(template2)

  # Define a vector version of each template for subsampling
  x1 <- as(phyloseq::otu_table(template1), "numeric")
  x2 <- as(phyloseq::otu_table(template2), "numeric")

  # Create 'dirty' multinomial. Defined artificial mixing. Create mixed
  # multinomial by adding counts from the other in precise proportion, a total
  # of Library Size / Effect Size
  add_to_template1 <- round((sum(x1) * x2 / (sum(x2) * unmixfac)), 0)

  # Add them together to create "dirty" multinomial
  # This adds the fractional subsampled abundances to the template
  mat1 <- matrix((add_to_template1 + x1), nrow = 1)
  rownames(mat1) <- phyloseq::sample_names(template1)
  colnames(mat1) <- phyloseq::taxa_names(template1)
  phyloseq::otu_table(template1) <- phyloseq::otu_table(mat1,
                                                        taxa_are_rows = FALSE)

  return(template1)
}

# rescale the sum of reads in the original raw(-ish) template data
# to the expected library size being requested here
sumsim <- function(n_l, sumtemplate, samples_per_type) {
# `n_l` - expected size target
# `sumtemplate` - the template vector of library sizes observed in template
# `samples_per_type` - The number of sample sizes to return

  scaled_sums <- round(n_l * (sumtemplate / median(sumtemplate)))
  rand_depths <- sample(scaled_sums, size = samples_per_type, replace = TRUE)

  return(rand_depths)
}

microbesim <- function(postfix = "sim", template, templatex, unmixfac,
                       samples_per_type, n_l) {
# Generate `samples_per_type` simulated microbiomes with `n_l` total reads each
# (all the same, or n_l has length equal to the value of `samples_per_type`),
# with subsamples drawn from `template`.
# `postfix` is a dummy idenitifer added to help distinguish
# simulated samples in downstream code.

# Perform the mixing here, so each replicate is a separate random mix
  template <- mix_microbes(template, templatex, unmixfac)

# call the proporitions vector `pi`, similar to nomenclature from DMN
  pi <- phyloseq::taxa_sums(template)

# n_l must be a scalar (recycled as the number of reads for every simulation) or
# it can be vector of length equal to J, the number of samples being simulated.
  if (length(samples_per_type) != 1) {
    stop("Length of samples_per_type should be 1.")
  }
  if (length(samples_per_type) != 1 & length(n_l) != samples_per_type) {
    stop("n_l should be length 1, or length samples_per_type.")
  }

# Actually create the simulated abundance table
  simat <- mapply(
    function(i, x, sample_size){
      if (FALSE) print(i)  # i is a dummy iterator
      phyloseq:::rarefaction_subsample(x, sample_size)
    },
    i = 1:samples_per_type,
    sample_size = n_l,
    MoreArgs = list(x = pi),
    SIMPLIFY = TRUE
  )
  simat <- t(simat)

# Add the OTU names to the OTU (column) indices
  colnames(simat) <- names(pi)

# Add new simulated sample_names to the row (sample) indices
  rownames(simat) <- paste0(seq_len(nrow(simat)), postfix)

# Put simulated abundances together with metadata as a phyloseq object
  otu <- phyloseq::otu_table(simat, taxa_are_rows = FALSE)

# Define data.frame that will become sample_data
  sdf <- data.frame(sample = phyloseq::sample_names(otu),
                    type = "simulated",
                    postfix = postfix,
                    stringsAsFactors = FALSE)
  rownames(sdf) <- phyloseq::sample_names(otu)
  
  sdf <- phyloseq::sample_data(sdf)

# Return a phyloseq object
  return(phyloseq::phyloseq(otu, sdf))
}



# Rarely a simulation has a weird value and fails.
# Catch these with `try`, and repeat the simulation call
# if error (it will be a new seed)
try_again <- TRUE
infiniteloopcounter <- 1

while (try_again && infiniteloopcounter < 5) {

  n1   <- sumsim(n_l, sampsums, samples_per_type)
  n2   <- sumsim(n_l, sampsums, samples_per_type)

  sim1 <- microbesim(sampletypes[1], template1, template2,
                     mixfac, samples_per_type, n1)
  sim2 <- microbesim(sampletypes[2], template2, template1,
                     mixfac, samples_per_type, n2)

  if (is.null(sim1) || is.null(sim2) || is.null(n1) || is.null(n2) ||
        inherits(sim1, "try-error") || inherits(sim2, "try-error")) {
    try_again <- TRUE
    infiniteloopcounter <- infiniteloopcounter + 1
  } else {
    try_again <- FALSE
  }
}

if (infiniteloopcounter >= 5) {
  stop("Consistent error found during simulation. Need to investigate cause.")
}


# Merge the two simulated datasets together into one phyloseq object
# and add back tree. then remove any otus that have zero abundance across all
# samples
sim <- phyloseq::merge_phyloseq(sim1, sim2)
sim <- phyloseq::merge_phyloseq(sim,
                                phyloseq::tax_table(GlobalPatterns),
                                phyloseq::phy_tree(GlobalPatterns))
sim <- prune_taxa(taxa_sums(sim) > 0, sim)


saveRDS(sim, file_name)