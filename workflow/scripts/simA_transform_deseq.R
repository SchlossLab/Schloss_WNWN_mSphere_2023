#!/usr/bin/env Rscript

set.seed(19760620)

# borrowed from L478-504 and L581 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate singe otu_count file from
# another otu_count file
#
# be sure to first run: conda activate nr-s1 

library(phyloseq)
library(DESeq)

args <- commandArgs(trailingOnly = TRUE)

raw_otu_file <- args[1]

deseq_otu_file <- gsub("RDS", "deseq.RDS", raw_otu_file)
raw_otu_data <- readRDS(raw_otu_file)
deseq_otu_data <- raw_otu_data

# this was the default value in the deseq_varstab function and it is not altered
# elsewhere in the code
sample_conditions <- rep("A", nsamples(raw_otu_data))

# Enforce orientation.
if (!phyloseq::taxa_are_rows(raw_otu_data)) {
  raw_otu_data <- t(raw_otu_data)
}

x <- as(phyloseq::otu_table(raw_otu_data), "matrix")

# The same tweak as for edgeR to avoid NaN problems
# that cause the workflow to stall/crash.
x <- x + 1

# Create annotated data.frame with the taxonomy table, in case it is useful
# later
tax_adf <- as(data.frame(as(phyloseq::tax_table(raw_otu_data), "matrix"),
                         stringsAsFactors = FALSE), "AnnotatedDataFrame")

cds <- DESeq::newCountDataSet(x, sample_conditions, featureData = tax_adf)

# First estimate library size factors
cds <- DESeq::estimateSizeFactors(cds)

# Variance estimation, passing along additional options
# these were the options that were passed along at L581
cds <- DESeq::estimateDispersions(cds,
                                  method = "blind",
                                  sharingMode = "maximum",
                                  fitType = "local")

# Determine which column(s) have the dispersion estimates
dispcol <- grep("disp\\_", colnames(fData(cds)))

# Enforce that there are no infinite values in the dispersion estimates
if (any(!is.finite(fData(cds)[, dispcol]))) {
  fData(cds)[which(!is.finite(fData(cds)[, dispcol])), dispcol] <- 0.0
}

vsmat <- exprs(DESeq::varianceStabilizingTransformation(cds))
otu_table(deseq_otu_data) <- otu_table(vsmat, taxa_are_rows = TRUE)

saveRDS(deseq_otu_data, deseq_otu_file)
