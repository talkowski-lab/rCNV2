#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot ROC analysis of various gene fine-mapping approaches


options(scipen=100000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Load truth set
load.truth <- function(truth.in){
  x <- read.table(truth.in, sep="\t", header=T, comment.char="")
  colnames(x)[1] <- "HPO"
  return(x)
}

# Compute ROC
roc <- function(stats, steps=seq(1, 0, -0.001)){
  roc_res <- as.data.frame(t(sapply(steps, function(minPIP){
    idxs <- which(stats$PIP > minPIP)
    ftrue <- length(which(stats$true[idxs])) / length(which(stats$true))
    fother <- length(which(!stats$true[idxs])) / length(which(!stats$true))
    fall <- length(idxs) / nrow(stats)
    return(c(minPIP, fall, fother, ftrue))
  })))
  colnames(roc_res) <- c("maxPIP", "frac_all", "frac_other", "frac_true")
  return(roc_res)
}

# Load a single dataset and annotate vs truth sets
load.data.single <- function(path, truth){
  x <- read.table(path, header=T, sep="\t", comment.char="")
  colnames(x)[1] <- "HPO"
  truth$true <- TRUE
  x <- merge(x, truth, all.x=T, all.y=F, by=c("HPO", "gene"))
  x$true[which(is.na(x$true))] <- FALSE
  return(x)
}

# Wrapper to load all datasets
load.datasets <- function(data.in, truth){
  datlist <- read.table(data.in, header=F, sep="\t")
  colnames(datlist) <- c("name", "color", "path")
  data <- lapply(1:nrow(datlist), function(i){
    stats <- load.data.single(datlist$path[i], truth)
    roc.res <- roc(stats)
    roc.opt <- optimize.roc(roc)
    return(list("stats"=stats,
                "roc"=roc.res,
                "roc.opt"=roc.opt,
                "color"=datlist$color[i]))
  })
  names(data) <- datlist$name
  return(data)
}

# ROC plot
plot.roc <- function(data){
  par(mar=c(3, 3, 1, 1))
  plot(x=c(0, 1), y=c(0, 1), type="n",
       xaxs="i", yaxs="i", xlab="", ylab="")
  abline(0, 1, col="gray70")
  lapply(data, function(x){
    points(x$roc$frac_other, x$roc$frac_true,
           type="l", col=x$color, lwd=3)
  })
  legend("bottomright", legend=names(data), lwd=5, cex=0.75, 
         col=sapply(data, function(x){x$color}), bty="n")
  mtext(1, line=2, text="Fraction of other associations retained")
  mtext(2, line=2, text="Fraction of true positive associations retained")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog roc_input.tsv truth_genes.tsv",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Three positional arguments: data.tsv, truth_genes.tsv, and path to plot.pdf\n")
}

# Writes args & opts to vars
data.in <- args$args[1]
truth.in <- args$args[2]
plot.out <- args$args[3]

# DEV PARAMTERS
setwd("~/scratch")
data.in <- "finemap_roc_input.tsv"
truth.in <- "constrained_truth_set.tsv"
plot.out <- "finemap_roc.pdf"

# Read truth genes
truth <- load.truth(truth.in)

# Load each input, annotate vs. truth, and store as list
data <- load.datasets(data.in, truth)

# Plot ROC
pdf(plot.out, height=4, width=4)
plot.roc(data)
dev.off()
