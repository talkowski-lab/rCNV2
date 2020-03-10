#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Join two or more gene metadata BED files


options(scipen=100000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Load metadata
load.meta <- function(mdat.in){
  x <- read.table(mdat.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(x)[1] <- gsub("^#", "", colnames(x)[1])
  return(x)
}

# Join list of 2+ metadata files
join.mdats <- function(mdat.list){
  merged <- mdat.list[[1]]
  for(i in 2:length(mdat.list)){
    merged <- merge(merged, mdat.list[[2]],
                    by=c("chr", "start", "end", "gene"),
                    all=F, sort=F)
  }
  return(merged)
}


#####################
### RSCRIPT BLOCK ###
#####################
# Use standard positional args to accept any number of metadata files
mdat.list.in <- commandArgs(trailingOnly=T)

# Load metadata files
mdat.list <- lapply(mdat.list.in, load.meta)

# Merge metadata files
mdat <- join.mdats(mdat.list)

# Write out merged metadata file
colnames(mdat)[1] <- paste("#", colnames(mdat)[1], sep="")
write.table(mdat, stdout(), col.names=T, row.names=F,
            sep="\t", quote=F)

