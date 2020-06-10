#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot FDR permutation results from sliding window meta-analysis for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000, check.names=F)


######################
### DATA FUNCTIONS ###
######################
# Load HPO sample size information
load.hpos <- function(hpos.in){
  hpos <- read.table(hpos.in, sep="\t", header=T, comment.char="")
  data.frame("hpo"=gsub(":", "", hpos[, 1], fixed=T),
             "n"=as.numeric(hpos[, 3]))
}

# Calculate FDR across a range of cutoffs for a vector of p-values
calc.fdr <- function(pv, cutoffs){
  pv <- pv[which(!is.na(pv))]
  n.pv <- length(pv)
  sapply(cutoffs, function(x){
    length(which(pv>=x))/n.pv
  })
}

# Read p-values and calculate FDR CDF
load.fdrs <- function(pvals.in, cutoffs, cnvtype=NULL){
  pvals <- read.table(pvals.in, header=T, sep="\t", comment.char="")
  if(!is.null(cnvtype)){
    pvals <- pvals[, grep(cnvtype, colnames(pvals), fixed=T)]
  }
  fdr.mat <- apply(pvals, 2, calc.fdr, cutoffs=cutoffs)
  rownames(fdr.mat) <- cutoffs
  return(fdr.mat)
}

# Calculate minimum P-value cutoff corresponding to target FDR for a vector of p-values
hit.fdr.target <- function(fdrv, cutoffs, target){
  cutoffs[min(c(head(which(fdrv<=target), 1), length(fdrv)))]
}

# Calculate P-value cutoffs per permutation corresponding to target FDRs for a matrix of p-values
get.fdr.cutoffs <- function(fdr.mat, fdr.target){
  fdr.cutoffs <- as.data.frame(apply(fdr.mat, 2, hit.fdr.target, cutoffs=cutoffs, target=fdr.target))
  colnames(fdr.cutoffs) <- "fdr.cutoff"
  return(fdr.cutoffs)
}

# Create cutoff table of all permutations
make.cutoff.mat <- function(fdr.mat, hpos, fdr.target){
  perm.hpos <- unlist(lapply(strsplit(colnames(fdr.mat), split=".", fixed=T), function(vals){vals[1]}))
  perm.n <- as.numeric(sapply(perm.hpos, function(hpo){hpos$n[which(hpos$hpo==hpo)]}))
  case.frac <- perm.n / (perm.n + hpos$n[which(hpos$hpo=="HEALTHY_CONTROL")])
  perm.idx <- unlist(lapply(strsplit(colnames(fdr.mat), split=".", fixed=T), function(vals){vals[3]}))
  fdr.res <- get.fdr.cutoffs(fdr.mat, fdr.target)
  cutoff.mat <- data.frame("hpo"=perm.hpos,
                           "perm"=perm.idx,
                           "n.cases"=perm.n,
                           "case.frac"=case.frac,
                           fdr.res)
  rownames(cutoff.mat) <- NULL
  return(cutoff.mat)
}

# Calculate a statistic for permuted cutoffs for each phenotype
get.cutoff.stat <- function(cutoff.mat, stat){
  stat.df <- as.data.frame(do.call("rbind", lapply(unique(cutoff.mat$hpo), function(hpo){
    if(stat=="max"){
      x <- max(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], na.rm=T)
    }else if(stat=="mean"){
      x <- mean(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], na.rm=T)
    }else if(stat=="median"){
      x <- median(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], na.rm=T)
    }else if(stat=="q3"){
      x <- quantile(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], probs=0.75, na.rm=T)
    }
    c(hpo,
      head(cutoff.mat[which(cutoff.mat$hpo==hpo), 3:4], 1),
      unlist(x))
  })))
  colnames(stat.df) <- c("hpo", "n.cases", "case.frac", "fdr.cutoff")
  stat.df[, -1] <- apply(stat.df[, -1], 2, as.numeric)
  return(stat.df)
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)
require(funr, quietly=T)
require(MASS, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--max-cutoff"), type="numeric", default=20,
              help="max P-value cutoff to evaluate (Phred-scaled) [default %default]"),
  make_option(c("--cutoff-step"), type="numeric", default=0.05,
              help="P-value increments to evaluate (Phred-scaled) [default %default]"),
  make_option(c("--fdr-target"), type="numeric", default=0.01,
              help="FDR target [default %default]"))

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog del.tsv dup.tsv hpo_table.tsv out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

if(length(args$args) != 4){
  stop("Incorrect number of positional arguments specified.")
}

del.pvals.in <- args$args[1]
dup.pvals.in <- args$args[2]
hpos.in <- args$args[3]
out.prefix <- args$args[4]
max.cutoff <- opts$`max-cutoff`
cutoff.step <- opts$`cutoff-step`
fdr.target <- opts$`fdr-target`

# # DEV PARAMETERS
# del.pvals.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.sliding_window.meta_analysis.stats.permuted_p_values.tsv.gz"
# dup.pvals.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DUP.sliding_window.meta_analysis.stats.permuted_p_values.tsv.gz"
# hpos.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/sliding_windows_meta_fdr_perm_test"
# max.cutoff <- 20
# cutoff.step <- 0.05
# fdr.target <- 0.000003715428

# Load sample size info
hpos <- load.hpos(hpos.in)

# Set FDR cutoffs
cutoffs <- seq(0, max.cutoff, cutoff.step)

# Calculate FDR CDFs for each permutation
del.fdr.mat <- load.fdrs(del.pvals.in, cutoffs, "DEL")
dup.fdr.mat <- load.fdrs(dup.pvals.in, cutoffs, "DUP")

# Calculate p-value cutoffs for each permutation for each target
del.cutoff.mat <- make.cutoff.mat(del.fdr.mat, hpos, fdr.target)
dup.cutoff.mat <- make.cutoff.mat(dup.fdr.mat, hpos, fdr.target)

# Compute median cutoffs
del.median.cutoffs <- get.cutoff.stat(del.cutoff.mat, "median")
dup.median.cutoffs <- get.cutoff.stat(dup.cutoff.mat, "median")


