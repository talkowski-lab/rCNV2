#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Determine best gene meta-analysis parameters based on distribution of lambda ratios

# Note: assumes mixture of parameters considered during Feb 2020 methods development
#       at a later date, some of these parameters/frequencies/models may be obsolete

# Requires output from compare_gene_burden_models.sh as input


options(scipen=1000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Load tsv of lambda ratios
load.tsv <- function(tsv.in){
  x <- read.table(tsv.in, header=T, sep="\t")
  
  return(x)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# Read arguments
args <- commandArgs(trailingOnly=T)
tsv.in <- as.character(args[1])
out.prefix <- as.character(args[2])

# DEV args
tsv.in <- "~/scratch/rCNV2.gene_burden_meta_analysis.lambda_ratios.tsv"

# Load data
x <- load.tsv(tsv.in)



