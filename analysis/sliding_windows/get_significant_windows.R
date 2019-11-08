#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Extract significant sliding windows from a p-value matrix


options(scipen=100000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Generic helper function to load numeric BED-style matrix
load.bed <- function(infile){
  x <- read.table(infile, header=T, comment.char="", sep="\t")
  colnames(x)[1] <- "chr"
  x[, -1] <- apply(x[, -1], 2, as.numeric)
  return(x)
}


# Load p-value matrix
load.pvalues <- function(pvalues.in, bed, p.is.phred){
  pvalues <- load.bed(pvalues.in)
  if(p.is.phred==T){
    pvalues[, -c(1:3)] <- apply(pvalues[, -c(1:3)], 2, function(p){10^-p})
  }
  pvalues <- merge(pvalues, bed, all.x=F, all.y=T, 
                   by=c("chr", "start", "end"), sort=F)
  return(pvalues)
}


# Load odds ratio matrix
load.ors <- function(ors.in, bed, or.is.ln){
  ors <- load.bed(ors.in)
  if(or.is.ln==T){
    ors[, -c(1:3)] <- apply(ors[, -c(1:3)], 2, exp)
  }
  ors <- merge(ors, bed, all.x=F, all.y=T, 
               by=c("chr", "start", "end"), sort=F)
  return(ors)
}


# Load nominal cohort count matrix
load.nomsig <- function(nomsig.in, bed){
  nomsig <- load.bed(nomsig.in)
  nomsig <- merge(nomsig, bed, all.x=F, all.y=T, 
               by=c("chr", "start", "end"), sort=F)
  return(nomsig)
}


# Identify significant windows based on a numeric matrix
get.sig.idxs <- function(x, cutoff, direction){
  apply(x[, -c(1:3)], 2, function(vals){
    if(direction=="ge"){
      which(vals>=cutoff)
    }else if(direction=="le"){
      which(vals<=cutoff)
    }else{
      stop(paste("get.sig.windows: direction", direction, "is not recognized."))
    }
  })
}


# Merge lists of significant window indexes
merge.sig <- function(bed, sig.pvals, sig.ors, sig.nom){
  if(is.null(sig.pvals)){
    stop("merge.sig requires p-values to be supplied.")
  }
  sig.all <- lapply(1:length(sig.pvals), function(i){
    sig <- 1:nrow(bed)
    if(!is.null(sig.pvals)){
      sig <- intersect(sig, sig.pvals[[i]])
    }
    if(!is.null(sig.ors)){
      sig <- intersect(sig, sig.ors[[i]])
    }
    if(!is.null(sig.nom)){
      sig <- intersect(sig, sig.nom[[i]])
    }
  })
  names(sig.all) <- names(sig.pvals)
  return(sig.all)
}


# Annotate matrix of bins by significance
mark.sig.windows <- function(bed, sig.merged){
  nvals.all <- sapply(1:length(sig.merged), function(i){
    nvals <- rep(FALSE, nrow(bed))
    if(length(sig.merged[[i]]) > 0){
      nvals[sig.merged[[i]]] <- TRUE
    }
    return(nvals)
  })
  colnames(nvals.all) <- names(sig.merged)
  cbind(bed, nvals.all)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--pvalues"), default=NULL, 
              help="matrix of p-values per window per phenotype. [default %default]"),
  make_option(c("--p-is-phred"), type="logical", default=F, action="store_true",
              help="supplied P-values are Phred-scaled (-log10[P]). [default %default]"),
  make_option(c("--max-p"), type="numeric", default=0.05, 
              help="maximum P-value to consider significant. [default %default]"),
  make_option(c("--odds-ratios"), default=NULL, 
              help="matrix of odds ratios (or lower bounds) per window per phenotype. [default %default]"),
  make_option(c("--or-is-ln"), type="logical", default=F, action="store_true",
              help="supplied odds ratios are natural log-scaled. [default %default]"),
  make_option(c("--min-or"), type="numeric", default=1, 
              help="minimum odds ratio to consider significant. [default %default]"),
  make_option(c("--nominal-counts"), default=NULL, 
              help="matrix of nominally significant cohorts per window per phenotype. [default %default]"),
  make_option(c("--min-nominal"), type="numeric", default=1, 
              help="minimum number of nominal cohorts to consider significant. [default %default]"),
  make_option(c("-o", "--out-prefix"), type="character", default="./significant_windows.",
              help="prefix for writing out all results. [default %default]", metavar="path")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog windows",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 1){
  stop("Must supply unannotated windows as only positional argument.\n")
}

# Writes args & opts to vars
bed.in <- args$args[1]
pvalues.in <- opts$pvalues
p.is.phred <- opts$`p-is-phred`
max.p <- opts$`max-p`
ors.in <- opts$`odds-ratios`
or.is.ln <- opts$`or-is-ln`
min.or <- opts$`min-or`
nomsig.in <- opts$`nominal-counts`
min.nom <- opts$`min-nominal`
out.prefix <- opts$`out-prefix`

# # DEV PARAMETERS:
# bed.in <- "~/scratch/GRCh37.200kb_bins_10kb_steps.raw.bed.gz"
# pvalues.in <- "~/scratch/DEL.pval_matrix.bed.gz"
# p.is.phred <- T
# max.p <- 0.000000629506182857198
# ors.in <- "~/scratch/DEL.lnOR_lower_matrix.bed.gz"
# or.is.ln <- T
# min.or <- 2
# nomsig.in <- "~/scratch/DEL.nominal_cohort_counts.bed.gz"
# min.nom <- 2
# out.prefix <- "~/scratch/sig_bins.test."

# Read windows
bed <- load.bed(bed.in)

# Load p-values and determine significant windows, if provided
if(!is.null(pvalues.in)){
  pvalues <- load.pvalues(pvalues.in, bed, p.is.phred)
  sig.pvals <- get.sig.idxs(pvalues, max.p, "le")
}else{
  sig.pvals <- NULL
}

# Load odds ratios and determine significant windows, if provided
if(!is.null(ors.in)){
  ors <- load.ors(ors.in, bed, or.is.ln)
  sig.ors <- get.sig.idxs(ors, min.or, "ge")
}else{
  sig.ors <- NULL
}

# Load count of nominally significant cohorts and determine significant windows, if provided
if(!is.null(nomsig.in)){
  nomsig <- load.nomsig(nomsig.in, bed)
  sig.nom <- get.sig.idxs(nomsig, min.nom, "ge")
}else{
  sig.nom <- NULL
}

# Determine final significant windows
sig.merged <- merge.sig(bed, sig.pvals, sig.ors, sig.nom)
bed.sigAnno <- mark.sig.windows(bed, sig.merged)
bed.sigOnly <- bed[apply(bed.sigAnno[, -c(1:3)], 1, function(vals){any(vals)}), 1:3]

# Write out
colnames(bed.sigAnno)[1] <- "#chr"
write.table(bed.sigAnno, paste(out.prefix, "all_windows_labeled.bed", sep=""), 
            sep="\t", quote=F, row.names=F, col.names=T)
colnames(bed.sigOnly)[1] <- "#chr"
write.table(bed.sigOnly, paste(out.prefix, "significant_windows.bed", sep=""), 
            sep="\t", quote=F, row.names=F, col.names=T)

