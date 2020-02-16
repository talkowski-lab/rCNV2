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


# Load .tsv of p-value cutoff matrix per HPO
load.p.cutoffs <- function(p.cutoffs.in){
  p.cutoffs <- read.table(p.cutoffs.in, header=T, sep="\t", comment.char="")
  colnames(p.cutoffs) <- c("hpo", "max.p")
  p.cutoffs$max.p <- as.numeric(p.cutoffs$max.p)
  return(p.cutoffs)
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


# Identify significant windows based on p-value
get.sig.pval.idxs <- function(pvalues, p.cutoffs){
  vals <- lapply(4:ncol(pvalues), function(i){
    hpo <- unlist(strsplit(colnames(pvalues)[i], split=".", fixed=T))[1]
    max.p <- p.cutoffs$max.p[which(p.cutoffs$hpo == hpo)]
    which(pvalues[, i] <= max.p)
  })
  names(vals) <- colnames(pvalues)[-c(1:3)]
  return(vals)
}


# Identify significant windows based on a numeric matrix
get.sig.idxs <- function(x, cutoff, direction){
  lapply(4:ncol(x), function(i){
    vals <- x[, i]
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
merge.sig <- function(bed, sig.pvals, sig.secondary.pvals, sig.ors, sig.nom, secondary.or.nom=F){
  if(is.null(sig.pvals)){
    stop("merge.sig requires p-values to be supplied.")
  }
  if(secondary.or.nom==T){
    sig.sec.or.nom <- lapply(1:length(sig.secondary.pvals), function(x){
      sort(unique(c(sig.secondary.pvals[[x]], sig.nom[[x]])))
    })
    names(sig.sec.or.nom) <- names(sig.pvals)
  }
  sig.all <- lapply(1:length(sig.pvals), function(i){
    sig <- 1:nrow(bed)
    if(!is.null(sig.pvals)){
      sig <- intersect(sig, sig.pvals[[i]])
    }
    if(!is.null(sig.ors)){
      sig <- intersect(sig, sig.ors[[i]])
    }
    if(secondary.or.nom==F){
      if(!is.null(sig.secondary.pvals)){
        sig <- intersect(sig, sig.secondary.pvals[[i]])
      }
      if(!is.null(sig.nom)){
        sig <- intersect(sig, sig.nom[[i]])
      }
    }else{
      sig <- intersect(sig, sig.sec.or.nom[[i]])
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
  make_option(c("--secondary-pvalues"), default=NULL, 
              help="matrix of secondary p-values per window per phenotype. [default %default]"),
  make_option(c("--p-is-phred"), type="logical", default=F, action="store_true",
              help="supplied P-values are Phred-scaled (-log10[P]). [default %default]"),
  make_option(c("--p-cutoffs"), default=NULL, 
              help="tsv of maximum P-value to consider significant per phenotype. [default: 0.05 for all]"),
  make_option(c("--odds-ratios"), default=NULL, 
              help="matrix of odds ratios (or lower bounds) per window per phenotype. [default %default]"),
  make_option(c("--or-is-ln"), type="logical", default=F, action="store_true",
              help="supplied odds ratios are natural log-scaled. [default %default]"),
  make_option(c("--min-secondary-p"), type="numeric", default=1, 
              help="minimum secondary P-value to consider significant. [default %default]"),
  make_option(c("--min-or"), type="numeric", default=1, 
              help="minimum odds ratio to consider significant. Supply as unscaled OR (will be transformed if needed). [default %default]"),
  make_option(c("--nominal-counts"), default=NULL, 
              help="matrix of nominally significant cohorts per window per phenotype. [default %default]"),
  make_option(c("--min-nominal"), type="numeric", default=1, 
              help="minimum number of nominal cohorts to consider significant. [default %default]"),
  make_option(c("--secondary-or-nom"), type="logical", default=F, action="store_true",
              help="include windows that meet either secondary P-value or min. nominal significant cohort count. [default %default]"),
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
secondary.pvalues.in <- opts$`secondary-pvalues`
p.is.phred <- opts$`p-is-phred`
p.cutoffs.in <- opts$`p-cutoffs`
ors.in <- opts$`odds-ratios`
or.is.ln <- opts$`or-is-ln`
min.secondary.p <- opts$`min-secondary-p`
min.or <- opts$`min-or`
nomsig.in <- opts$`nominal-counts`
min.nom <- opts$`min-nominal`
secondary.or.nom <- opts$`secondary-or-nom`
out.prefix <- opts$`out-prefix`

# # DEV PARAMETERS:
# bed.in <- "~/scratch/GRCh37.200kb_bins_10kb_steps.raw.bed.gz"
# pvalues.in <- "~/scratch/DUP.pval_matrix.bed.gz"
# secondary.pvalues.in <- "~/scratch/DUP.secondary_pval_matrix.bed.gz"
# p.cutoffs.in <- "~/scratch/sliding_window.rCNV.DUP.empirical_genome_wide_pval.hpo_cutoffs.tsv"
# p.is.phred <- T
# ors.in <- "~/scratch/DUP.lnOR_lower_matrix.bed.gz"
# or.is.ln <- T
# min.secondary.p <- 0.05
# min.or <- 1
# nomsig.in <- "~/scratch/DUP.nominal_cohort_counts.bed.gz"
# min.nom <- 2
# secondary.or.nom <- T
# out.prefix <- "~/scratch/sig_bins.test."

# Read windows
bed <- load.bed(bed.in)

# Load p-values and determine significant windows, if provided
if(!is.null(pvalues.in)){
  pvalues <- load.pvalues(pvalues.in, bed, p.is.phred)
  p.cutoffs <- load.p.cutoffs(p.cutoffs.in)
  sig.pvals <- get.sig.pval.idxs(pvalues, p.cutoffs)
}else{
  sig.pvals <- NULL
}

# Load secondary p-values and determine significant windows, if provided
if(!is.null(pvalues.in)){
  secondary.pvalues <- load.pvalues(secondary.pvalues.in, bed, p.is.phred)
  secondary.p.cutoffs <- p.cutoffs
  secondary.p.cutoffs$max.p <- min.secondary.p
  sig.secondary.pvals <- get.sig.pval.idxs(secondary.pvalues, secondary.p.cutoffs)
}else{
  sig.secondary.pvals <- NULL
}

# Load odds ratios and determine significant windows, if provided
if(!is.null(ors.in)){
  if(or.is.ln==TRUE){
    min.or <- log(min.or)
  }
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
sig.merged <- merge.sig(bed, sig.pvals, sig.secondary.pvals, sig.ors, sig.nom, secondary.or.nom)
bed.sigAnno <- mark.sig.windows(bed, sig.merged)
bed.sigOnly <- bed[apply(bed.sigAnno[, -c(1:3)], 1, function(vals){any(vals)}), 1:3]

# Write out
colnames(bed.sigAnno)[1] <- "#chr"
write.table(bed.sigAnno, paste(out.prefix, "all_windows_labeled.bed", sep=""), 
            sep="\t", quote=F, row.names=F, col.names=T)
colnames(bed.sigOnly)[1] <- "#chr"
write.table(bed.sigOnly, paste(out.prefix, "significant_windows.bed", sep=""), 
            sep="\t", quote=F, row.names=F, col.names=T)

