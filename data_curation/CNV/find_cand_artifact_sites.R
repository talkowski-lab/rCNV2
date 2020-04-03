#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Identify windows with CNVs found at a much higher frequency in one cohort than any other


options(stringsAsFactors=F, scipen=1000)


#################
### FUNCTIONS ###
#################

# Read frequency matrix & filter to bins with at least one cohort above min.top
read.freqs <- function(freq.in, min.top=0.01, max.mean=0.05){
  x <- read.table(freq.in, sep="\t", header=T)
  cat(paste(prettyNum(nrow(x), big.mark=","), "total bins"))
  # Find bins with top > min.top
  rowmax <- apply(x[, -c(1:3)], 1, max, na.rm=T)
  pass <- which(rowmax>=min.top)
  cat(paste(prettyNum(length(pass), big.mark=","), "/",
              prettyNum(nrow(x), big.mark=","), " bins ",
              "(", round(100*length(pass)/nrow(x), 1), 
              "%) with at least one platform above ",
              round(100*min.top, 1), "% frequency", sep=""))
  x <- x[pass, ]
  # Find bins with mean < max.mean
  rowmean <- apply(x[, -c(1:3)], 1, mean, na.rm=T)
  pass2 <- which(rowmean<max.mean)
  cat(paste(prettyNum(length(pass2), big.mark=","), "/",
              prettyNum(nrow(x), big.mark=","), " bins ",
              "(", round(100*length(pass2)/nrow(x), 1), 
              "%) with mean platform frequency below ",
              round(100*max.mean, 1), "% frequency", sep=""))
  x[pass2, ]
}

# Find candidate artifacts
find.arts <- function(x, min.diff=0.01){
  rowmax <- apply(x[, -c(1:3)], 1, max, na.rm=T)
  rowmean <- apply(x[, -c(1:3)], 1, function(vals){
    best <- head(which(vals==max(vals, na.rm=T)), 1)
    max(vals[-c(best)], na.rm=T)
  })
  rowdiff <- rowmax - rowmean
  pass3 <- which(rowdiff>min.diff)
  cat(paste(prettyNum(length(pass3), big.mark=","), "/",
              prettyNum(nrow(x), big.mark=","), " bins ",
              "(", round(100*length(pass3)/nrow(x), 1), 
              "%) where absolute difference between top platform and mean platform freq is above ",
              round(100*min.diff, 1), "% frequency", sep=""))
  x[pass3, 1:3]
}



#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--min-diff"), type="numeric", default=0.01, 
              help="Maxmimum difference in CNV frequency between highest and mean-highest to flag. [default %default]"),
  make_option(c("--min-top"), type="numeric", default=0.01, 
              help="Minimum frequency in top cohort to flag. [default %default]"),
  make_option(c("--max-mean"), type="numeric", default=0.05, 
              help="Maximum frequency in mean cohort to flag. [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog freq_matrix.bed.gz outfile.bed",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Two positional arguments: freq_matrix.bed.gz and outfile.bed\n")
}

# Writes args & opts to vars
freq.in <- args$args[1]
outbed <- args$args[2]
min.diff <- opts$`min-diff`
min.top <- opts$`min-top`
max.mean <- opts$`max-mean`

# # DEV PARAMTERS
# setwd("~/scratch")
# freq.in <- "DUP.freq_matrix.bed.gz"
# outbed <- "DUP.candidate_artifacts.bed"
# min.diff <- 0.01
# min.top <- 0.01
# max.mean <- 0.05

# Read frequency matrix
x <- read.freqs(freq.in, min.top, max.mean)

# Identify candidate artifacts
arts <- find.arts(x, min.diff)

# Write list of potential artifacts
colnames(arts) <- c("#chr", "start", "end")
write.table(arts, outbed, col.names=T, row.names=F, quote=F, sep="\t")
