#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute performance of various ML models trained for gene scoring


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load a single raw score (output directly from gene scoring procedure)
load.raw.score <- function(path){
  x <- read.table(path, header=T, sep="\t", comment.char="")[, c(1, 3)]
  colnames(x) <- c("gene", "score")
  x[, 2] <- as.numeric(x[, 2])
  return(x)
}
  
# Load all scores from an input .tsv of scores and names
load.all.scores <- function(scores.in){
  score.list <- read.table(scores.in, header=F, sep="\t")
  score.list <- score.list[order(score.list[, 1]), ]
  scores <- lapply(score.list[, 2], load.raw.score)
  names(scores) <- score.list[, 1]
  return(scores)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot simple legend
simple.legend <- function(labels, colors){
  par(mar=rep(0.25, 4), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, length(labels)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  points(x=rep(0.1, length(labels)), y=(1:length(labels))-0.5, pch=22, 
         col=blueblack, bg=colors, cex=1.8)
  text(x=rep(0.1, length(labels)), y=(1:length(labels))-0.58, pos=4, labels=labels, xpd=T)
}

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
args <- parse_args(OptionParser(usage=paste("%prog scores_list.tsv true.genes false.genes out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores_list.tsv, true.genes, false.genes, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
true.genes.in <- args$args[2]
false.genes.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# setwd("~/scratch/")
# scores.in <- "evaluation.DEL.input.tsv"
# cnv.type <- "DEL"
# true.genes.in <- "gold_standard.haploinsufficient.genes.list"
# false.genes.in <- "gold_standard.haplosufficient.genes.list"
# out.prefix <- "ml_model_eval_test"
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
scores <- load.all.scores(scores.in)
evals <- lapply(scores, evaluate.score, score="score", 
                true.genes=true.genes, false.genes=false.genes)

# Plot ROC
pdf(paste(out.prefix, "model_eval.roc.pdf", sep="."),
    height=2.75, width=2.75)
plot.roc(evals)
dev.off()

# Plot PRC
pdf(paste(out.prefix, "model_eval.prc.pdf", sep="."),
    height=2.75, width=2.75)
plot.prc(evals)
dev.off()

# Plot legend
pdf(paste(out.prefix, "model_eval.legend.pdf", sep="."),
    height=2, width=2)
simple.legend(labels=rev(ml.model.abbrevs[names(evals)]),
              colors=viridis(length(evals)))
dev.off()

