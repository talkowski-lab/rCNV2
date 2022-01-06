#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute performance of pHaplo vs pTriplo scores for specific positive and negative gene sets for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(flux, quietly=T)
require(Hmisc, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv true.genes false.genes out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, true.genes, false.genes, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
true.genes.in <- args$args[2]
false.genes.in <- args$args[3]
out.prefix <- args$args[4]

# # DEV PARAMETERS
# setwd("~/scratch/")
# scores.in <- "rCNV.gene_scores.tsv.gz"
# true.genes.in <- "gold_standard.haploinsufficient.genes.list"
# false.genes.in <- "gold_standard.haplosufficient.genes.list"
# out.prefix <- "phi_vs_pts_test"

# Load gene lists
true.genes <- as.character(read.table(true.genes.in, header=F)[, 1])
false.genes <- as.character(read.table(false.genes.in, header=F)[, 1])

# Load & evaluate all scores in input list
scores <- load.scores(scores.in)
evals <- lapply(c("pHaplo", "pTriplo"), evaluate.score, stats=scores,
                true.genes=true.genes, false.genes=false.genes)

# Plot ROC
pdf(paste(out.prefix, "pHaplo_vs_pTriplo.roc.pdf", sep="."),
    height=2.75, width=2.75)
plot.roc(evals, colors=cnv.colors[1:2], auc.text.colors=rep("white", 2), grid.col=NA)
dev.off()

# Plot PRC
pdf(paste(out.prefix, "pHaplo_vs_pTriplo.prc.pdf", sep="."),
    height=2.75, width=2.75)
plot.prc(evals, colors=cnv.colors[1:2], auc.text.colors=rep("white", 2), grid.col=NA)
dev.off()

