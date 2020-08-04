#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot quality assessment analyses of gene scores


options(stringsAsFactors=F, scipen=1000)
colors <- c("#D43925", "#2376B2", "#7459B2")


#################
### FUNCTIONS ###
#################
# Load gene score stats .tsv
load.stats <- function(path){
  x <- read.table(path, header=T, sep="\t", comment.char="")[, 1:3]
  colnames(x)[1] <- "gene"
  x[, 2:3] <- apply(x[, 2:3], 2, as.numeric)
  return(x)
}

# Load truth sets
load.truth <- function(truth.in){
  truth.list <- read.table(truth.in, header=F, sep="\t", comment.char="")
  truth.sets <- lapply(1:nrow(truth.list), function(i){
    x <- read.table(truth.list[i, 2], sep="\t", header=F, comment.char="")[, 1]
    return(list("truth.genes"=x, "name"=truth.list[i, 1], 
                "path"=truth.list[i, 2], "color"=truth.list[i, 3]))
  })
  names(truth.sets) <- truth.list[, 1]
  return(truth.sets)
}

# Compute ROC
roc <- function(stats, score, truth.genes, neg.genes, steps=seq(1, 0, -0.001)){
  x <- data.frame("score" = stats[, which(colnames(stats) == score)],
                  "true" = stats$gene %in% truth.genes,
                  "neg" = stats$gene %in% neg.genes)
  roc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score > k)
    ftrue <- length(which(x$true[idxs])) / length(which(x$true))
    fneg <- length(which(x$neg[idxs])) / length(which(x$neg))
    fother <- length(which(!x$true[idxs])) / length(which(!x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, fother, ftrue, fneg))
  })))
  colnames(roc_res) <- c("min_score", "frac_all", "frac_other", "frac_true", "frac_neg")
  return(roc_res)
}

# Compute PRC
prc <- function(stats, score, truth.genes, neg.genes, steps=seq(1, 0, -0.001)){
  x <- data.frame("score" = stats[, which(colnames(stats) == score)],
                  "true" = stats$gene %in% truth.genes,
                  "neg" = stats$gene %in% neg.genes)
  prc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score > k)
    prec <- length(which(x$true[idxs])) / (length(which(x$true[idxs])) + length(which(x$neg[idxs])))
    recall <- length(which(x$true[idxs])) / length(which(x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, prec, recall))
  })))
  colnames(prc_res) <- c("minPIP", "frac_all", "precision", "recall")
  return(prc_res)
}

# Wrapper to calculate all plotting data
calc.plot.data <- function(stats, score, truth, neg.genes){
  data <- lapply(1:length(truth), function(i){
    roc.res <- roc(stats, score, truth[[i]]$truth.genes, neg.genes)
    roc.auc <- flux::auc(roc.res$frac_neg, roc.res$frac_true)
    prc.res <- prc(stats, score, truth[[i]]$truth.genes, neg.genes)
    prc.auc <- flux::auc(prc.res$recall, prc.res$precision)
    enrich.res <- data.frame("min_score"=roc.res$min_score,
                             "enrichment"=roc.res$frac_true / roc.res$frac_all)
    return(list("roc"=roc.res,
           "roc.auc"=roc.auc,
           "prc"=prc.res,
           "prc.auc"=prc.auc,
           "enrich"=enrich.res,
           "color"=truth[[i]]$color))
  })
  names(data) <- names(truth)
  return(data)
}

# ROC plot
plot.roc <- function(data, title="Receiver Operating Characteristic"){
  par(mar=c(3, 3, 1.5, 1))
  plot(x=c(0, 1), y=c(0, 1), type="n",
       xaxs="i", yaxs="i", xlab="", ylab="")
  abline(0, 1, col="gray70", lty=2)
  lorder <- order(-sapply(data, function(x){x$roc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$roc$frac_neg, x$roc$frac_true,
           type="l", col=x$color, lwd=3)
  })
  legend("bottomright", pch=19, pt.cex=1.5, cex=0.75, bty="n",
         col=sapply(data, function(x){x$color})[lorder], 
         legend=paste(names(data), " (AUC=",
                      sapply(data, function(x){format(round(x$roc.auc, 2), nsmall=2)}),
                      ")", sep="")[lorder])
  mtext(1, line=2, text="False positive rate")
  mtext(2, line=2, text="True positive rate")
  mtext(3, line=0.1, text=title)
}


# PRC plot
plot.prc <- function(data, title="Precision/Recall"){
  par(mar=c(3, 3, 1.5, 1))
  prec.max <- max(sapply(data, function(l){l$prc$precision}), na.rm=T)
  prec.baselines <- sapply(data, function(x){tail(x$prc$precision, 1)})
  plot(x=c(0, 1), y=c(0, prec.max), type="n",
       xaxs="i", yaxs="i", xlab="", ylab="")
  lorder <- order(-sapply(data, function(x){x$prc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    abline(h=prec.baselines[i], col=x$color, lty=2)
    points(x$prc$recall, x$prc$precision,
           type="l", col=x$color, lwd=2)
  })
  legend("topright", pch=19, pt.cex=1.5, cex=0.75, bty="n",
         col=sapply(data, function(x){x$color})[lorder], 
         legend=paste(names(data), " (AUC=",
                      sapply(data, function(x){format(round(x$prc.auc, 2), nsmall=2)}),
                      ")", sep="")[lorder])
  mtext(1, line=2, text="Recall")
  mtext(2, line=2, text="Precision")
  mtext(3, line=0.1, text=title)
}

# Enrichment plot
plot.enrichment <- function(data, title="Gene Set Enrichment"){
  ymax <- max(sapply(data, function(x){max(x$enrich$enrichment, na.rm=T)}), na.rm=T)
  par(mar=c(2.75, 4, 1.5, 1))
  plot(x=c(0, 1), y=c(0, ymax), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  abline(h=1, lty=2, col="gray80")
  sapply(1:length(data), function(i){
    points(data[[i]]$enrich, col=data[[i]]$color, type="l", lwd=2)
  })
  axis(1, at=seq(0, 1, 0.2), labels=NA)
  axis(1, at=seq(0, 1, 0.2), tick=F, line=-0.5)
  mtext(1, line=1.75, text="Min Score")
  axis(2, las=2)
  mtext(2, line=3, text="Fold-Enrichment")
  mtext(3, line=0.2, text=title, font=2)
  legend("topleft", pch=19, pt.cex=1.5, cex=0.75, bty="n",
         col=sapply(data, function(x){x$color}), 
         legend=names(data))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(flux, quietly=T)
require(Hmisc, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog scores.tsv del_truth.tsv dup_truth.tsv del_false.tsv dup_false.tsv out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop("Six positional arguments: scores.tsv, del_truth.tsv, dup_truth.tsv, del_false.tsv, dup_false.tsv, and out_prefix\n")
}

# Writes args & opts to vars
stats.in <- args$args[1]
del_truth.in <- args$args[2]
dup_truth.in <- args$args[3]
del_false.in <- args$args[4]
dup_false.in <- args$args[5]
out.prefix <- args$args[6]

# # DEV PARAMTERS
# setwd("~/scratch")
# stats.in <- "rCNV.gene_scores.tsv.gz"
# del_truth.in <- "DEL.roc_truth_sets.tsv"
# dup_truth.in <- "DUP.roc_truth_sets.tsv"
# del_false.in <- "gold_standard.haplosufficient.genes.list"
# dup_false.in <- "gold_standard.triplosensitive.genes.list"
# out.prefix <- "rCNV_gene_scoring_qc"

# Read gene score stats
stats <- load.stats(stats.in)

# Read truth sets
del.truth <- load.truth(del_truth.in)
dup.truth <- load.truth(dup_truth.in)
del.false <- read.table(del_false.in, header=F)[, 1]
dup.false <- read.table(dup_false.in, header=F)[, 1]

# Compute plot stats for DEL & DUP
del.data <- calc.plot.data(stats, "pHI", del.truth, del.false)
dup.data <- calc.plot.data(stats, "pTS", dup.truth, dup.false)

# Generate all DEL/pHI plots
pdf(paste(out.prefix, "pHI.roc.pdf", sep="."), height=4, width=4)
plot.roc(del.data)
dev.off()
pdf(paste(out.prefix, "pHI.prc.pdf", sep="."), height=4, width=4)
plot.prc(del.data)
dev.off()
pdf(paste(out.prefix, "pHI.enrich.pdf", sep="."), height=4, width=4)
plot.enrichment(del.data)
dev.off()

# Generate all DUP/pTS plots
pdf(paste(out.prefix, "pTS.roc.pdf", sep="."), height=4, width=4)
plot.roc(dup.data)
dev.off()
pdf(paste(out.prefix, "pTS.prc.pdf", sep="."), height=4, width=4)
plot.prc(dup.data)
dev.off()
pdf(paste(out.prefix, "pTS.enrich.pdf", sep="."), height=4, width=4)
plot.enrichment(dup.data)
dev.off()

