#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute RMSE, AUROC, and AUPRC for various gene scoring models


options(stringsAsFactors=F, scipen=1000)


#################
### FUNCTIONS ###
#################
# Load & clean original BFDPs
load.orig.bfdps <- function(bfdps.in){
  x <- read.table(bfdps.in, header=T, sep="\t", comment.char="")
  x[, which(colnames(x) %in% c("gene", "bfdp"))]
}

# Load a single gene score .tsv
load.stats <- function(path){
  x <- read.table(path, header=T, sep="\t", comment.char="")[, c(1, 3)]
  colnames(x)[1] <- "gene"
  x[, 2] <- as.numeric(x[, 2])
  return(x)
}

# Compute ROC
roc <- function(stats, score, true.genes, false.genes, steps=seq(1, 0, -0.001)){
  x <- data.frame("score" = stats[, which(colnames(stats) == score)],
                  "true" = stats$gene %in% true.genes,
                  "false" = stats$gene %in% false.genes)
  roc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score > k)
    ftrue <- length(which(x$true[idxs])) / length(which(x$true))
    ffalse <- length(which(x$false[idxs])) / length(which(x$false))
    fother <- length(which(!x$true[idxs])) / length(which(!x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, fother, ftrue, ffalse))
  })))
  colnames(roc_res) <- c("min_score", "frac_all", "frac_other", "frac_true", "frac_false")
  return(roc_res)
}

# Compute PRC
prc <- function(stats, score, true.genes, false.genes, steps=seq(1, 0, -0.001)){
  x <- data.frame("score" = stats[, which(colnames(stats) == score)],
                  "true" = stats$gene %in% true.genes,
                  "false" = stats$gene %in% false.genes)
  prc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score > k)
    prec <- length(which(x$true[idxs])) / (length(which(x$true[idxs])) + length(which(x$false[idxs])))
    recall <- length(which(x$true[idxs])) / length(which(x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, prec, recall))
  })))
  colnames(prc_res) <- c("minPIP", "frac_all", "precision", "recall")
  return(prc_res)
}

# Compute RMSE vs original BFDPs
calc.rmse <- function(stats, score, orig.bfdps){
  x <- merge(stats, orig.bfdps, all=F, sort=F, by="gene")
  keep.idx <- which(apply(x[, 2:3], 1, function(vals){all(!is.na(vals))}))
  sqrt(sum((x[keep.idx, 2] - x[keep.idx, 3])^2, na.rm=T) / length(keep.idx))
}

# Wrapper to calculate all plotting data for a single model
calc.plot.data <- function(stats, score, true.genes, false.genes){
  roc.res <- roc(stats, score, true.genes, false.genes)
  roc.auc <- flux::auc(roc.res$frac_false, roc.res$frac_true)
  prc.res <- prc(stats, score, true.genes, false.genes)
  prc.auc <- flux::auc(prc.res$recall, prc.res$precision)
  rmse <- calc.rmse(stats, score, orig.bfdps)
  return(list("roc"=roc.res,
              "roc.auc"=roc.auc,
              "prc"=prc.res,
              "prc.auc"=prc.auc,
              "rmse"=rmse))
}

# Wrapper to calculate all plotting data for all models
load.all.models <- function(inlist.in, orig.bdfps, true.genes, false.genes){
  inlist <- read.table(inlist.in, header=F, comment.char="", sep="\t")
  data <- lapply(1:nrow(inlist), function(i){
    stats <- load.stats(inlist[i, 2])
    calc.plot.data(stats, "score", true.genes, false.genes)
  })
  names(data) <- inlist[, 1]
  return(data)
}

# Compile summary table of models
get.sum.table <- function(data){
  x <- do.call("rbind", lapply(data, function(l){c(l$rmse, l$roc.auc, l$prc.auc)}))
  x <- as.data.frame(cbind(names(data), x))
  colnames(x) <- c("model", "mean_rmse", "auroc", "auprc")
  rownames(x) <- NULL
  return(x)
}


# ROC plot
plot.roc <- function(data, title="Receiver Operating Characteristic"){
  colors <- viridis(length(data))
  par(mar=c(3, 3, 1.5, 1))
  plot(x=c(0, 1), y=c(0, 1), type="n",
       xaxs="i", yaxs="i", xlab="", ylab="")
  abline(0, 1, col="gray70", lty=2)
  lorder <- order(-sapply(data, function(x){x$roc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$roc$frac_false, x$roc$frac_true,
           type="l", lwd=3, col=colors[i])
  })
  legend("bottomright", pch=19, pt.cex=1.5, cex=0.75, bty="n",
         col=colors[lorder], 
         legend=paste(names(data), " (AUC=",
                      sapply(data, function(x){format(round(x$roc.auc, 2), nsmall=2)}),
                      ")", sep="")[lorder])
  mtext(1, line=2, text="False positive rate")
  mtext(2, line=2, text="True positive rate")
  mtext(3, line=0.1, text=title)
}

# PRC plot
plot.prc <- function(data, title="Precision/Recall"){
  colors <- viridis(length(data))
  par(mar=c(3, 3, 1.5, 1))
  prec.max <- max(sapply(data, function(l){l$prc$precision}), na.rm=T)
  prec.baseline <- mean(sapply(data, function(x){tail(x$prc$precision, 1)}), na.rm=T)
  plot(x=c(0, 1), y=c(0, prec.max), type="n",
       xaxs="i", yaxs="i", xlab="", ylab="")
  abline(h=prec.baseline, lty=2)
  lorder <- order(-sapply(data, function(x){x$prc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$prc$recall, x$prc$precision,
           type="l", col=colors[i], lwd=2)
  })
  legend("topright", pch=19, pt.cex=1.5, cex=0.75, bty="n",
         col=colors[lorder], 
         legend=paste(names(data), " (AUC=",
                      sapply(data, function(x){format(round(x$prc.auc, 2), nsmall=2)}),
                      ")", sep="")[lorder])
  mtext(1, line=2, text="Recall")
  mtext(2, line=2, text="Precision")
  mtext(3, line=0.1, text=title)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(flux, quietly=T)
require(Hmisc, quietly=T)
require(viridisLite, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog inputs.tsv orig_bfdps.tsv true_genes.txt false_genes.txt out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop("Three positional arguments: inputs.tsv, orig_bfdps.tsv, true_genes.txt, false_genes.txt, and out_prefix\n")
}

# Writes args & opts to vars
inlist.in <- args$args[1]
bfdps.in <- args$args[2]
true_genes.in <- args$args[3]
false_genes.in <- args$args[4]
out.prefix <- args$args[5]

# # DEV PARAMTERS
# setwd("~/scratch")
# inlist.in <- "rCNV.DEL.model_evaluation.input.tsv"
# bfdps.in <- "rCNV.DEL.gene_abfs.tsv"
# true_genes.in <- "gold_standard.haploinsufficient.genes.list"
# false_genes.in <- "gold_standard.haplosufficient.genes.list"
# out.prefix <- "rCNV_gene_scoring_model_comparison"

# Load original BFDPs
orig.bfdps <- load.orig.bfdps(bfdps.in)

# Read gene lists
true.genes <- read.table(true_genes.in, header=F)[, 1]
false.genes <- read.table(false_genes.in, header=F)[, 1]

# Load data
data <- load.all.models(inlist.in, orig.bdfps, true.genes, false.genes)

# Write table with stats
write.table(get.sum.table(data),
            paste(out.prefix, "summary_table.tsv", sep="."), 
            col.names=T, row.names=F, sep="\t", quote=F)

# Plot ROC
pdf(paste(out.prefix, "roc.pdf", sep="."), height=4, width=4)
plot.roc(data)
dev.off()

# Plot PRC
pdf(paste(out.prefix, "prc.pdf", sep="."), height=4, width=4)
plot.prc(data)
dev.off()

