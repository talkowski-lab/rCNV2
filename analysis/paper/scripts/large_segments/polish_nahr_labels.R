#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot grid summarizing large segment association across phenotypes for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load unmodified segment table
load.raw.segs <- function(segs.in){
  segs <- read.table(segs.in, header=T, sep="\t", comment.char="")
  colnames(segs)[1] <- gsub("X.", "", colnames(segs)[1], fixed=T)
  return(segs)
}

# Load true NAHR labels
load.labels <- function(labels.in){
  labels <- read.table(labels.in, header=T, sep="\t", comment.char="")
  labels.v <- labels[, 2]
  names(labels.v) <- labels[, 1]
  return(labels.v)
}

# Overwrite NAHR labels based on manually curated true classifications
polish.nahr.labels <- function(segs, labels){
  segs$nahr[sapply(names(labels), function(id){which(segs$region_id == id)})] <- as.numeric(labels)
  return(segs)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog segs.tsv nahr_labels.tsv outfile", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: segs.tsv, nahr_labels.tsv, and outfile\n", sep=" "))
}

# Writes args & opts to vars
segs.in <- args$args[1]
labels.in <- args$args[2]
outfile <- args$args[3]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# labels.in <- "~/scratch/rCNV2_analysis_d1.correct_nahr_labels.tsv"
# outfile <- "~/scratch/rCNV2_analysis_d1.master_segments.polished_labels.bed"

# Load (raw) segment table
segs <- load.raw.segs(segs.in)

# Load correct NAHR labels
labels <- load.labels(labels.in)

# Polish NAHR labels
segs <- polish.nahr.labels(segs, labels)

# Write to outfile
colnames(segs)[1] <- paste("#", colnames(segs)[1], sep="")
write.table(segs, outfile, col.names=T, row.names=F, sep="\t", quote=F)
