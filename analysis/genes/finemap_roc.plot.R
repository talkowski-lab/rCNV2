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
  colnames(roc_res) <- c("minPIP", "frac_all", "frac_other", "frac_true")
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

# Find ROC-optimal point
optimize.roc <- function(roc.res){
  # Compute Euclidean distance from (0, 1)
  x2 <- (roc.res$frac_other - 0) ^ 2
  y2 <- (roc.res$frac_true - 1) ^ 2
  d <- sqrt(x2 + y2)
  d.best <- min(d)
  best.idx <- sort(which(d == d.best))[1]
  return(roc.res[best.idx, ])
}

# Wrapper to load all datasets
load.datasets <- function(data.in, truth){
  datlist <- read.table(data.in, header=F, sep="\t")
  colnames(datlist) <- c("name", "color", "lty", "path")
  data <- lapply(1:nrow(datlist), function(i){
    stats <- load.data.single(datlist$path[i], truth)
    roc.res <- roc(stats)
    roc.opt <- optimize.roc(roc.res)
    roc.auc <- flux::auc(roc.res$frac_other, roc.res$frac_true)
    return(list("stats"=stats,
                "roc"=roc.res,
                "roc.opt"=roc.opt,
                "midPIP"=roc.res[which(roc.res$minPIP == 0.5), ],
                "auc"=roc.auc,
                "color"=datlist$color[i],
                "lty"=datlist$lty[i]))
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
  lorder <- order(-sapply(data, function(x){x$auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$roc$frac_other, x$roc$frac_true,
           type="l", col=x$color, lwd=2, lty=x$lty)
  })
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$roc.opt$frac_other, x$roc.opt$frac_true,
           pch=19, col=x$color)
    points(x$midPIP$frac_other, x$midPIP$frac_true,
           pch=23, bg=x$color, col="black")
  })
  legend("bottomright", lwd=5, cex=0.75, bty="n",
         col=sapply(data, function(x){x$color})[lorder], 
         legend=paste(names(data), " (AUC=",
                      sapply(data, function(x){format(round(x$auc, 2), nsmall=2)}),
                      ")", sep="")[lorder])
  legend("topleft", cex=0.75, bty="n", pch=c(19, 23),
         legend=c("ROC-optimal cutoff", "PIP > 0.5"),
         pt.bg=c(NA, "gray50"))
  mtext(1, line=2, text="Fraction of other associations retained")
  mtext(2, line=2, text="Fraction of true positive associations retained")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(flux, quietly=T)

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
truth.in <- "clingen_truth_set.tsv"
plot.out <- "finemap_roc.pdf"

# Read truth genes
truth <- load.truth(truth.in)

# Load each input, annotate vs. truth, and store as list
data <- load.datasets(data.in, truth)

# Plot ROC
pdf(plot.out, height=4, width=4)
plot.roc(data)
dev.off()
