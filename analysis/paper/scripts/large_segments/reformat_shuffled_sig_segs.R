#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Reformat permuted significant large segments to match loci-style BED file


options(scipen=1000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Load table of permuted loci
load.perms <- function(perms.in){
  perm <- read.table(perm.in, header=T, sep="\t", comment.char="", check.names=F)
  perm$perm_region_id <- perm$region_id
  perm$region_id <- sapply(perm$region_id, function(id){
    gsub("perm[0-9]+_", "", id, fixed=F)
  })
  return(perm)
}

# Read table of original loci and restrict to columns to keep
# Also returns character vector of names of columns to keep
load.loci <- function(loci.in){
  loci <- read.table(loci.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(loci)[1] <- gsub("#", "", colnames(loci)[1], fixed=T)
  orig.columns <- colnames(loci)
  cols.to.drop <- c("chr", "start_min", "end_max", "cnv", "cred_interval_coords",
                    "n_genes", "genes")
  loci <- loci[, -which(colnames(loci) %in% cols.to.drop)]
  return(list(loci, orig.columns))
}

# Merge old association stats into permuted loci, and reorder to match original locus table
merge.tables <- function(perm, loci, orig.columns){
  df <- merge(perm, loci, by="region_id", all.x=T, all.y=F, sort=F, 
              suffixes=c(".perm", ".original"))

  df$region_id <- df$perm_region_id
  df$perm_region_id <- NULL
  colnames(df)[which(colnames(df) == "start")] <- "start_min"
  colnames(df)[which(colnames(df) == "end")] <- "end_max"
  
  df <- df[with(df, order(perm, chr, start_min, end_max)), ]
  df[orig.columns]
}


################
###RSCRIPT BLOCK
################
# Load required libraries
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("-z", "--bgzip"), action="store_true",  
              help="bgzip output BED [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog perm.bed loci.bed out.bed",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
perm.in <- args$args[1]
loci.in <- args$args[2]
outfile <- args$args[3]
bgzip <- opts$bgzip

# # Dev parameters
# setwd("~/scratch")
# perm.in <- "~/scratch/test_perms.w_genes.bed.gz"
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# outfile <- "~/scratch/reformatted_test_perms.bed"
# bgzip <- T

# Load data
perm <- load.perms(perm.in)
loci.list <- load.loci(loci.in)
loci <- loci.list[[1]]
orig.columns <- loci.list[[2]]

# Merge & reorder permuted & original loci
merged <- merge.tables(perm, loci, orig.columns)
colnames(merged)[1] <- paste("#", colnames(merged)[1], sep="")

# Write to outfile
if(tools::file_ext(outfile) %in% c("gz", "bgz", "gzip", "bgzip")){
  outfile <- tools::file_path_sans_ext(outfile)
}
write.table(merged, outfile, sep="\t", quote=F, col.names=T, row.names=F)
if(bgzip==T){
  system(paste("bgzip -f", outfile), wait=T, intern=F)
}