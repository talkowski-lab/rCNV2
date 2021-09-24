#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute BED of best P-value for each window from sliding window meta-analysis
# (Helper function to speed up segment permtuation)


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog in.bed out.bed", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: in.bed and out.bed \n", sep=" "))
}

# Writes args & opts to vars
pvals.in <- args$args[1]
outfile <- args$args[2]

# # DEV PARAMETERS
# pvals.in <- "~/scratch/rCNV2_analysis_d1.DEL.meta_neg_log10_p.all_hpos.bed.gz"
# outfile <- "~/scratch/test_max_pval_per_window.DEL.bed"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Load p-value matrix
pvals <- load.pval.matrix(pvals.in, has.coords=T, p.is.neg.log10=T)

# Compute best P-value across any phenotype per window, and write to output file
df.out <- pvals$coords
df.out$best_p <- -log10(apply(pvals$pvals, 1, min, na.rm=T))
df.out$best_p[which(is.infinite(df.out$best_p) | is.na(df.out$best_p) | is.nan(df.out$best_p))] <- 0
colnames(df.out)[1] <- paste("#", colnames(df.out)[1], sep="")
write.table(df.out, outfile, col.names=T, row.names=F, sep="\t", quote=F)
