#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Split genes based on meta-analysis odds ratios


options(stringsAsFactors=F, scipen=1000)


#################
### FUNCTIONS ###
#################
# Load stats .tsv
load.stats <- function(path, use.zscore=FALSE){
  x <- read.table(path, header=T, sep="\t", comment.char="")[, -c(1:3)]
  colnames(x)[1] <- "gene"
  if(use.zscore==TRUE){
    if(any(is.na(x$meta_z))){
      x <- x[-which(is.na(x$meta_z)), ]
    }
  }else{
    if(any(is.na(x$meta_lnOR))){
      x <- x[-which(is.na(x$meta_lnOR)), ]
    }
  }
  return(x)
}

# Extract top/bottom genes based on percentile
get.tails <- function(stats, pct, use.zscore){
  if(use.zscore==TRUE){
    x <- stats$meta_z
  }else{
    x <- stats$meta_lnOR
  }
  cutoffs <- quantile(x, probs=c(pct, 1-pct))
  bottom <- sort(stats$gene[which(x<=cutoffs[1])])
  top <- sort(stats$gene[which(x>=cutoffs[2])])
  tails <- list("top"=top, "bottom"=bottom)
  return(tails)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--pct"), default=0.5,
              help="Percentile cutoff for top & bottom"),
  make_option(c("--use-zscore"), action="store_true", default=FALSE,
              help="Threshold based on Z-score (rather than odds ratio)")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog stats.bed out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Two positional arguments: stats.bed and out_prefix\n")
}

# Writes args & opts to vars
stats.in <- args$args[1]
out.prefix <- args$args[2]
pct <- opts$pct
zscore <- opts$`use-zscore`

# # DEV PARAMS
# stats.in <- "~/scratch/HP0000118.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
# out.prefix <- "~/scratch/split.test"
# pct <- 0.5
# zscore <- T

# Load stats
stats <- load.stats(stats.in, zscore)

# Determine cutoffs & extract genes in tails
tails <- get.tails(stats, pct, zscore)

# Write to file
write.table(tails$top, paste(out.prefix, "top.genes.list", sep="."),
            col.names=F, quote=F, row.names=F)
write.table(tails$bottom, paste(out.prefix, "bottom.genes.list", sep="."),
            col.names=F, quote=F, row.names=F)

