#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute performance of pHI vs pTS scores for specific positive and negative gene sets for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(flux, quietly=T)
require(Hmisc, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

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
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# setwd("~/scratch/")
# scores.in <- "rCNV.gene_scores.tsv.gz"
# true.genes.in <- "gold_standard.haploinsufficient.genes.list"
# false.genes.in <- "gold_standard.haplosufficient.genes.list"
# out.prefix <- "phi_vs_pts_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load gene lists
true.genes <- as.character(read.table(true.genes.in, header=F)[, 1])
false.genes <- as.character(read.table(false.genes.in, header=F)[, 1])

# Load & evaluate all scores in input list
scores <- load.scores(scores.in)
evals <- lapply(c("pHI", "pTS"), evaluate.score, stats=scores,
                true.genes=true.genes, false.genes=false.genes)

# Plot ROC
pdf(paste(out.prefix, "pHI_vs_pTS.roc.pdf", sep="."),
    height=2.75, width=2.75)
plot.roc(evals, colors=cnv.colors[1:2], auc.text.colors=rep("white", 2))
dev.off()

# Plot PRC
pdf(paste(out.prefix, "pHI_vs_pTS.prc.pdf", sep="."),
    height=2.75, width=2.75)
plot.prc(evals, colors=cnv.colors[1:2], auc.text.colors=rep("white", 2))
dev.off()

