#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Extract list of genes with insufficient CNVs for gene score model training


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############

# Read an input file of association statistics
read.stats <- function(stats.in, prefix){
  # Read data & subset to necessary columns
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")
  colnames(stats)[1] <- "chr"
  cols.to.keep <- c("chr", "start", "end", "gene", "case_alt", "case_ref",
                    "control_alt", "control_ref")
  stats <- stats[, which(colnames(stats) %in% cols.to.keep)]
  colnames(stats)[-(1:4)] <- paste(prefix, colnames(stats)[-(1:4)], sep=".")
  return(stats)
}

# Merge a list of association statistics
combine.stats <- function(stats.list){
  # Merge all cohorts
  merged <- stats.list[[1]]
  for(i in 2:length(stats.list)){
    merged <- merge(merged, stats.list[[i]], 
                    by=c("chr", "start", "end", "gene"),
                    all=F, sort=F)
  }
  merged[, -c(1:4)] <- apply(merged[, -c(1:4)], 2, as.numeric)
  merged$total_alt <- apply(merged[, grep("_alt", colnames(merged), fixed=T)], 1, sum, na.rm=T)
  return(merged)
}


################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--min-cnvs"), default=5,
              help="Minimum number of total CNVs required per gene"),
  make_option(c("--gene-counts-out"), help="Output tsv with CNV counts per gene.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog infile outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
infile <- args$args[1]
outfile <- args$args[2]
min.cnvs <- opts$`min-cnvs`
counts.out <- opts$`gene-counts-out`

# # Dev parameters
# setwd("~/scratch/")
# infile <- "~/scratch/gene_meta_in.test.tsv"
# outfile <- "~/scratch/underpowered_genes.test.list"
# min.cnvs <- 5
# counts.out <- "~/scratch/counts.tsv"

# Read list of cohorts to meta-analyze
cohort.info <- read.table(infile, header=F, sep="\t")
ncohorts <- nrow(cohort.info)
stats.list <- lapply(1:ncohorts, function(i){read.stats(cohort.info[i, 2], 
                                                        cohort.info[i, 1])})
names(stats.list) <- cohort.info[, 1]

# Merge stats
stats.merged <- combine.stats(stats.list)

# Extract BED of genes with total_alt < min.cnvs and write to table
out.bed <- stats.merged[which(stats.merged$total_alt < min.cnvs), 1:4]
colnames(out.bed)[1] <- paste("#", colnames(out.bed)[1], sep="")
write.table(out.bed, outfile, col.names=T, row.names=F, quote=F, sep="\t")

# If optioned, write table of CNV counts per gene to file
if(!is.null(counts.out)){
  counts <- data.frame("gene"=stats.merged$gene, "cnvs"=stats.merged$total_alt)
  colnames(counts)[1] <- paste("#", colnames(counts)[1], sep="")
  write.table(counts, counts.out, col.names=T, row.names=F, quote=F, sep="\t")
}
