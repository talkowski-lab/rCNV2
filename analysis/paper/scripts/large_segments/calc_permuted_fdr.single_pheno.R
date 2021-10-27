#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute permuted FDR for a single phenotype from the sliding window meta-analyses


options(stringsAsFactors=F, scipen=1000, check.names=F)


######################
### DATA FUNCTIONS ###
######################
# Load meta-analysis summary stats & extract vector of primary P-values
load.pvals <- function(sumstats.in, p.colname="meta_neg_log10_p", p.is.neg.log10=T){
  sumstats <- read.table(sumstats.in, header=T, sep="\t", comment.char="", check.names=F)
  pvals <- as.numeric(as.vector(sumstats[, which(colnames(sumstats) == p.colname)]))
  if(p.is.neg.log10==T){
    pvals <- 10^-pvals
  }
  return(pvals)
}

# Compute empirical FDR
calc.fdr <- function(pvals, fdr.target){
  n <- length(pvals)
  n.true <- ceiling(fdr.target * n)
  sort(pvals)[n.true]
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--fdr-target"), type="numeric", default=10e-8,
              help="FDR target [default %default]"))

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog sumstats.bed outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

if(length(args$args) != 2){
  stop("Incorrect number of positional arguments specified.")
}

sumstats.in <- args$args[1]
outfile <- args$args[2]
fdr.target <- opts$`fdr-target`

# # DEV PARAMETERS
# sumstats.in <- "~/scratch/HP0012759.rCNV.DEL.sliding_window.meta_analysis.stats.perm_1.bed.gz"
# outfile <- "~/scratch/sliding_windows_meta_fdr_perm_test.tsv"
# fdr.target <- 0.000003715428

# Read sumstats & extract P-values
pvals <- load.pvals(sumstats.in)

# Compute empirical cutoff for equivalent FDR
fdr <- calc.fdr(pvals, fdr.target)

# Write to outfile
write.table(data.frame(fdr), outfile, col.names=F, row.names=F, sep="\t", quote=F)

