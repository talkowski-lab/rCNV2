#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Extract significant windows from meta-analysis summary stats


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load meta-analysis summary stats
load.sumstats <- function(sumstats.in){
  sumstats <- read.table(sumstats.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(sumstats)[1] <- gsub("#", "", colnames(sumstats)[1], fixed=T)
  keep.cols <- c("chr", "start", "end", "n_nominal_cohorts", "meta_neg_log10_p", "meta_neg_log10_p_secondary", "meta_neg_log10_fdr_q")
  sumstats[, which(colnames(sumstats) %in% keep.cols)]
}

# Filter summary stats based on P-value cutoffs
filter.sumstats <- function(sumstats, primary.cutoff=6, secondary.cutoff=-log10(0.05),
                            min.nominal=2, secondary.or.nominal=FALSE,
                            fdr.cutoff=NULL){
  primary.hits <- which(sumstats$meta_neg_log10_p >= primary.cutoff)
  secondary.hits <- which(sumstats$meta_neg_log10_p_secondary >= secondary.cutoff)
  nominal.hits <- which(sumstats$n_nominal_cohorts >= min.nominal)
  if(secondary.or.nominal){
    hits <- intersect(primary.hits, unique(c(secondary.hits, nominal.hits)))
  }else{
    hits <- intersect(primary.hits, intersect(secondary.hits, nominal.hits))
  }
  if(!is.null(fdr.cutoff)){
    fdr.hits <- which(sumstats$meta_neg_log10_fdr_q >= fdr.cutoff)
    hits <- unique(sort(c(hits, fdr.hits)))
  }
  sumstats[hits, ]
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--primary-p-cutoff"), default=10e-6, help="Significant primary P-value cutoff."),
  make_option(c("--secondary-p-cutoff"), default=0.05, help="Significant secondary P-value cutoff."),
  make_option(c("--min-nominal"), default=2, help="Minimum number of nominally significant individual cohorts."),
  make_option(c("--secondary-or-nominal"), action="store_true", help="Require either secondary P cutoff and/or min nominal."),
  make_option(c("--fdr-q-cutoff", help="Supplement with FDR Q-value cutoff. [default: do not consider FDR]"))
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog sumstats.bed sig_windows_out.bed", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Seven positional arguments required: sumstats.bed and sign_windows_out.bed\n"))
}

# Writes args & opts to vars
sumstats.in <- args$args[1]
outfile <- args$args[2]
primary.cutoff <- -log10(as.numeric(opts$`primary-p-cutoff`))
secondary.cutoff <- -log10(as.numeric(opts$`secondary-p-cutoff`))
min.nominal <- as.numeric(opts$`min-nominal`)
secondary.or.nominal <- opts$`secondary-or-nominal`
fdr.cutoff <- -log10(as.numeric(opts$`fdr-q-cutoff`))

# # DEV PARAMETERS
# sumstats.in <- "~/scratch/HP0012759.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz"
# outfile <- "~/scratch/test.bed"
# primary.cutoff <- -log10(10e-6)
# secondary.cutoff <- -log10(0.05)
# min.nominal <- 2
# secondary.or.nominal <- T
# fdr.cutoff <- -log10(0.01)

# Load summary stats
sumstats <- load.sumstats(sumstats.in)

# Filter sumstats
sumstats.filtered <- filter.sumstats(sumstats, primary.cutoff, secondary.cutoff,
                                     min.nominal, secondary.or.nominal,
                                     fdr.cutoff)

# Write output bed
colnames(sumstats.filtered)[1] <- paste("#", colnames(sumstats.filtered)[1], sep="")
write.table(sumstats.filtered[, 1:3], outfile, col.names=T, row.names=F, sep="\t", quote=F)
