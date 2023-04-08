#!/usr/bin/env Rscript

# borrowed from L604-722 and L776-788 of simulation-cluster-accuracy-server.Rmd
# code edited for clarity and style and to generate file with list of distances
# from a phyloseq formatted otu_count file
#
# be sure to first run: conda activate nr-s1 

################################################################################
# Distance function adapted for this analysis
# Expects a DGEList as main input.
# Can calculate sample-wise distances for both "biological coefficient of variation" 
# or the lead log fold change ("logFC") distance using "pairwise"" or "common"
# `gene.selection` methods.
# From edgeR::plotMDS.DGEList:
# "Distances on the plot represent coefficient of variation of expression 
# between samples for the top genes that best distinguish the samples"
# root mean squared difference for top features.
# (of the logFC transformed counts, via edgeR::cpm)
# Root mean squared difference of the logFC transformed counts of the top OTUs only.
################################################################################
# Adapted from edgeR::plotMDS.DGEList
#	Multidimensional scaling plot of digital gene expression profiles
#	Yunshun Chen, Mark Robinson and Gordon Smyth
#	23 May 2011.  Last modified 7 Feb 2013.
################################################################################
edgeRdist <- function(x, top = 500, method = "logFC",
                      gene.selection = "pairwise", prior.count = 2) {
	require("edgeR")
  require("phyloseq")
  if (inherits(x, "phyloseq")) {
    # Convert phyloseq object to an otherwise untransformed 
    # abundance matrix in DGEList form
    x = edgeRnorm(x, method="none")

		if(is.na(x)) {return(NA)}
    if( !inherits(x, "DGEList") ){
      stop("x not converted to DGEList for some reason. Dispatch now flawed.")
    }
  }

	#	Remove rows with missing or Inf values
	ok <- is.finite(x$counts)
	if(!all(ok)) {
		x <- x[rowSums(ok)>0, ]
		x$samples$lib.size <- rowSums(x$counts)
	}
	nprobes <- nrow(x)
	nsamples <- ncol(x)
  # Check that `top` exceeds `nprobes` by at least a small margin.
  # Give warning if you have to change it.
  if( top > (0.8*nprobes) ){
    warning("The value of `top` was ", top, " for a dataset with ", nprobes, 
            "probes. \n Decreasing `top` to ", round(0.8*nprobes), " probes.")
    top = 0.8*nprobes
  }
  
	# Initialize distance object, dd
	cn = colnames(x)
	dd = matrix(0, nrow=nsamples, ncol=nsamples, dimnames=list(cn, cn))
	########################################
	#	Default method is to convert to moderated logCPM and call limma plotMDS
	########################################
	method <- match.arg(method, c("logFC", "bcv", "BCV"))
	if(method=="logFC"){
		y <- edgeR::cpm(x, log=TRUE, prior.count=prior.count)
		# Here we're essentially dispatching to 
		# code derived from the limma::plotMDS function,
		# but only the dist-calculating portion.

		#	Check top. `top` is handled differently between logFC and bcv
		top <- min(top, nprobes)
		#	Check gene.selection
		gene.selection <- match.arg(gene.selection, c("pairwise", "common"))
		#	Distance matrix from pairwise leading fold changes
		dd <- matrix(0, nrow=nsamples, ncol=nsamples, dimnames=list(cn, cn))
		topindex <- nprobes - top + 1
		if(gene.selection=="pairwise") {
			# Distance measure is mean of top squared deviations for each pair of arrays
			for (i in 2:(nsamples))
				for (j in 1:(i-1)){
					dd[i, j] <- sqrt(mean(
						sort.int((y[, i]-y[, j])^2, partial=topindex)[topindex:nprobes]
					))
				}
		} else {
			# gene.selection == "common".
			#	Same genes used for all comparisons
			s <- rowMeans((y-rowMeans(y))^2)
			q <- quantile(s, p=(topindex-1.5)/(nprobes-1))
			y <- y[s >= q, ]
			for (i in 2:(nsamples)){
				dd[i, 1:(i-1)] <- sqrt(colMeans((y[, i]-y[, 1:(i-1), drop=FALSE])^2))
			}
		}
		return(as.dist(dd))
	}
	########################################
	########################################
	#	Code past this point is for method="bcv"
	########################################
	x$samples$group <- factor(rep.int(1, nsamples))
	#	Check value for top
	if (top < nprobes) { 
		twd <- edgeR::estimateTagwiseDisp(edgeR::estimateCommonDisp(x), grid.length=100) 
		o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:top]
		subdata <- x$counts[o,,drop=FALSE]
	} else {
		subdata<-x$counts
	}

	lib.size <- x$samples$lib.size * x$samples$norm.factors
	myFun <- function(delta, y, der=1){
		return(sum(condLogLikDerDelta(y, delta, der)))
	}

	for(i in 2:(nsamples)){
		for(j in 1:(i - 1)){
			mm <- subdata[,c(i,j)]
			rs5 <- rowSums(mm) > 5
			lib <- lib.size[c(i, j)]
			norm <- t(t(mm)/lib) * exp(mean(log(lib)))
			delta <- optimize(myFun, interval=c(0.0001, .99), tol=0.000001,
												maximum=TRUE, y=norm[rs5, ], der=0)
			dd[i, j] = sqrt( delta$maximum / (1-delta$maximum) )
		}
	}
	return(as.dist(dd))
}

edgeRnorm = function(physeq, ...){
	require("edgeR")
  require("phyloseq")
  #physeq = simlist[["55000_1e-04"]]
  #z0 = simlisttmm[["55000_1e-04"]]
  #physeq = simlist[["1000_0.2"]]
  #z0 = simlisttmm[["1000_0.2"]]
	# Enforce orientation.
	if( !taxa_are_rows(physeq) ){
		physeq <- t(physeq)
	}
	x = as(otu_table(physeq), "matrix")


#methods cannot take negative values
	if(any(x < 0)){
		return(NA)
	}


  # See if adding a single observation, 1, 
  # everywhere (so not zeros) prevents errors
  # without needing to borrow and modify 
  # calcNormFactors (and its dependent functions)
  # It did. This fixed all problems. 
  # Can the 1 be reduced to something smaller and still work?
  x = x + 1
  # Now turn into a DGEList
	y = edgeR::DGEList(counts=x, remove.zeros=TRUE)
	# Perform edgeR-encoded normalization, using the specified method (...)
	z = edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
	# Don't need the following additional steps, which are also
	# built-in to some of the downstream distance methods.
	#z1 = estimateCommonDisp(z)
	#z2 = estimateTagwiseDisp(z1)
	return(z)
}

################################################################################

library(phyloseq)
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

transformed_otu_file <- args[1]
method <- args[2] # c("logFC", "bcv")

distance_file <- gsub("RDS", paste0(method, ".RDS"), transformed_otu_file)

mean_distance <- NA

edger_to_phyloseq <- function(edger) {

  counts <- as.data.frame(edger$count)
  factors <- edger$samples

  if (!any(is.na(factors$norm.factors))) {

    samples <- data.frame(sample = gsub("\\d*", "", colnames(counts)))
    rownames(samples) <- colnames(counts)

    transformed_otu <- phyloseq(otu_table(counts, taxa_are_rows = TRUE),
                                sample_data(samples),
                                phy_tree(GlobalPatterns))
  } else {
    transformed_otu <- NA
  }

  transformed_otu
}


get_distances <- function(otu_counts) {

  d <- NA

  # get the mean distance from across the elements of the list. this is
  # necessary when using rarefied datasets. in previous experience, the sd
  # is miniscule relative to the mean
  if (is.list(otu_counts)) {

		distances <- lapply(otu_counts,
												edgeRdist,
												method = method, top = 2000,
												gene.selection = "pairwise")
    d <- Reduce("+", distances) / length(distances)

  } else {

    d <- edgeRdist(otu_counts, method = method, top = 2000,
									gene.selection = "pairwise")

  }

  return(d)

}

# pds: wnwn didn't use these with the edgeR approaches. not sure why not. will
# use the counts element from the list. suspect won't use / using wrong with
# these distance calculators. if any of the factor values are NA/NaN will return
# an NA

mean_distance <- NA
transformed_otu_list <- readRDS(transformed_otu_file)

if (grepl("upperquartile|RLE|TMM", transformed_otu_file)) {

  library(edgeR)
  data(GlobalPatterns)

  transformed_otu_list <- lapply(transformed_otu_list, edger_to_phyloseq)

}

if (any(is.na(transformed_otu_list))) {
  mean_distance <- rep(NA, length(transformed_otu_list))

} else {
  mean_distance <- lapply(transformed_otu_list, get_distances)
}


saveRDS(mean_distance, distance_file)
