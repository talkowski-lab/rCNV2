#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Extract list of genes with insufficient CNVs for gene score model training


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############


################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--max-se"), default=1,
              help="Maximum standard error to tolerate")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog meta.stats out.bed",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
infile <- args$args[1]
outfile <- args$args[2]
max.se <- opts$`max-se`

# # Dev parameters
# setwd("~/scratch/")
# infile <- "~/scratch/rCNV2_analysis_d2.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
# outfile <- "~/scratch/underpowered_genes.test.list"
# max.se <- 1

# Load meta-analysis stats
meta.stats <- read.table(infile, header=T, sep="\t", comment.char="", check.names=F)

# Compute SE for each gene (assuming 95% CI)
meta.stats$se <- (meta.stats$meta_lnOR_upper - meta.stats$meta_lnOR) / qnorm(0.975)

# Extract BED of genes with SE > max.se and write to table
out.bed <- meta.stats[which(meta.stats$se > max.se), 1:4]
write.table(out.bed, outfile, col.names=T, row.names=F, quote=F, sep="\t")
