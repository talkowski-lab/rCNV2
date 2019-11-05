#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Extract significant sliding windows from a p-value matrix


options(scipen=1000, stringsAsFactors=F)


# Load P-value matrix
load.windows <- function(bed.in, p.is.phred){
  bed <- read.table(bed.in, header=T, comment.char="", sep="\t")
  colnames(bed)[1] <- "chr"
  bed[, -1] <- apply(bed[, -1], 2, as.numeric)
  if(p.is.phred==T){
    bed[, -c(1:3)] <- apply(bed[, -c(1:3)], 2, function(p){10^-p})
  }
  return(bed)
}


# Identify significant windows
get.sig.windows <- function(bed, cutoff){
  keep.idx <- which(sapply(1:nrow(bed), function(i){
    p.i <- as.numeric(bed[i, -c(1:3)])
    p.i <- p.i[which(!is.na(p.i))]
    any(p.i <= cutoff)
  }))
  bed[keep.idx, 1:3]
}


# Pad windows
pad.windows <- function(bed, pad){
  bed[, 2] <- sapply(bed[, 2], function(x){max(c(0, x-pad))})
  bed[, 3] <- bed[, 3]+pad
  return(bed)
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--p-is-phred"), type="logical", default=F, action="store_true",
              help="supplied P-values are Phred-scaled (-log10[P]) [default %default]"),
  make_option(c("--cutoff"), type="numeric", default=10^-8, 
              help="P-value of significance threshold [default %default]",
              metavar="numeric"),
  make_option(c("--pad-windows"), type="numeric", default=0, 
              help="amount to pad left and right coordinates of significant windows [default %default]",
              metavar="numeric"),
  make_option(c("-o", "--outfile"), type="character", default="stdout",
              help="output file [default %default]", metavar="path")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog windows",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 1){
  stop("Incorrect number of required positional arguments\n")
}

# Writes args & opts to vars
bed.in <- args$args[1]
p.is.phred <- opts$`p-is-phred`
cutoff <- opts$cutoff
pad <- opts$`pad-windows`
outfile <- opts$outfile

# # DEV PARAMETERS:
# bed.in <- "~/scratch/DEL.pval_matrix.bed.gz"
# p.is.phred <- T
# cutoff <- 0.0000192931
# pad <- 1000000
# outfile <- "stdout"

# Read windows
bed <- load.windows(bed.in, p.is.phred)

# Find significant windows
bed.sig <- get.sig.windows(bed, cutoff)

# Pad significant windows, if optioned
if(pad > 0){
  bed.sig <- pad.windows(bed.sig, pad)
}

# Write out
colnames(bed.sig)[1] <- "#chr"
if(outfile %in% c("stdout", "/dev/stdout", "-")){
  outfile <- stdout()
}
write.table(bed.sig, outfile, sep="\t", quote=F,
            row.names=F, col.names=T)

