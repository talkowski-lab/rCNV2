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
    idxs <- which(x$score >= k)
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
    idxs <- which(x$score >= k)
    prec <- length(which(x$true[idxs])) / (length(which(x$true[idxs])) + length(which(x$false[idxs])))
    recall <- length(which(x$true[idxs])) / length(which(x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, prec, recall))
  })))
  colnames(prc_res) <- c("minPIP", "frac_all", "precision", "recall")
  return(prc_res)
}

# # Compute RMSE vs original BFDPs
# calc.rmse <- function(stats, score, orig.bfdps){
#   x <- merge(stats, orig.bfdps, all=F, sort=F, by="gene")
#   keep.idx <- which(apply(x[, 2:3], 1, function(vals){all(!is.na(vals))}))
#   sqrt(sum((x[keep.idx, 2] - x[keep.idx, 3])^2, na.rm=T) / length(keep.idx))
# }

# Wrapper to calculate all plotting data for a single model
calc.plot.data <- function(stats, score, true.genes, false.genes){
  roc.res <- roc(stats, score, true.genes, false.genes)
  roc.auc <- flux::auc(roc.res$frac_false, roc.res$frac_true)
  prc.res <- prc(stats, score, true.genes, false.genes)
  prc.auc <- flux::auc(prc.res$recall, prc.res$precision)
  # rmse <- calc.rmse(stats, score, orig.bfdps)
  return(list("roc"=roc.res,
              "roc.auc"=roc.auc,
              "prc"=prc.res,
              "prc.auc"=prc.auc))
}

# Wrapper to calculate all plotting data for all models
load.all.models <- function(inlist.in, true.genes, false.genes){
  inlist <- read.table(inlist.in, header=F, comment.char="", sep="\t")
  data <- lapply(1:nrow(inlist), function(i){
    stats <- load.stats(inlist[i, 2])
    calc.plot.data(stats, "score", true.genes, false.genes)
  })
  names(data) <- inlist[, 1]
  return(data)
}

# Compile summary table of models
get.sum.table <- function(del, dup){
  sumdat <- t(sapply(names(del), function(model){
    del.auroc <- del[[which(names(del)==model)]]$roc.auc
    del.auprc <- del[[which(names(del)==model)]]$prc.auc
    dup.auroc <- dup[[which(names(dup)==model)]]$roc.auc
    dup.auprc <- dup[[which(names(dup)==model)]]$prc.auc
    mean.auc <- harmonic.mean(c(del.auroc, del.auprc, dup.auroc, dup.auprc))
    c(del.auroc, del.auprc, dup.auroc, dup.auprc, mean.auc)
  }))
  
  x <- as.data.frame(cbind(names(del), sumdat))
  x[, 2:ncol(x)] <- apply(x[, 2:ncol(x)], 2, as.numeric)
  colnames(x) <- c("model", "del.auroc", "del.auprc", "dup.auroc", "dup.auprc", "mean.auc")
  rownames(x) <- NULL
  return(x[order(-x$mean.auc), ])
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
  legend("bottomleft", pch=19, pt.cex=1.5, cex=0.75, bty="n",
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
require(psych, quietly=T)
require(viridisLite, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog del.tsv dup.tsv del.true.genes dup.true.genes del.false.genes dup.false.genes out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 7){
  stop("Seven positional arguments: del.tsv, dup.tsv, del.true.genes, dup.true.genes, del.false.genes, dup.false.genes, and out.prefix\n")
}

# Writes args & opts to vars
del.in <- args$args[1]
dup.in <- args$args[2]
del.true_genes.in <- args$args[3]
dup.true_genes.in <- args$args[4]
del.false_genes.in <- args$args[5]
dup.false_genes.in <- args$args[6]
out.prefix <- args$args[7]

# # DEV PARAMTERS
# setwd("~/scratch")
# del.in <- "rCNV.DEL.model_evaluation.input.tsv"
# dup.in <- "rCNV.DUP.model_evaluation.input.tsv"
# del.true_genes.in <- "gold_standard.haploinsufficient.genes.list"
# dup.true_genes.in <- "gold_standard.triplosensitive.genes.list"
# del.false_genes.in <- "gold_standard.haplosufficient.genes.list"
# dup.false_genes.in <- "gold_standard.triploinsensitive.genes.list"
# out.prefix <- "rCNV_gene_scoring_model_comparison"

# Read gene lists
del.true.genes <- read.table(del.true_genes.in, header=F)[, 1]
dup.true.genes <- read.table(dup.true_genes.in, header=F)[, 1]
del.false.genes <- read.table(del.false_genes.in, header=F)[, 1]
dup.false.genes <- read.table(dup.false_genes.in, header=F)[, 1]

# Load data
del <- load.all.models(del.in, del.true.genes, del.false.genes)
dup <- load.all.models(dup.in, dup.true.genes, dup.false.genes)

# Compute table with stats
write.table(get.sum.table(del, dup),
            paste(out.prefix, "summary_table.tsv", sep="."), 
            col.names=T, row.names=F, sep="\t", quote=F)

# Plot ROCs
pdf(paste(out.prefix, "del.roc.pdf", sep="."), height=4, width=4)
plot.roc(del, title="Deletion ROC")
dev.off()
pdf(paste(out.prefix, "dup.roc.pdf", sep="."), height=4, width=4)
plot.roc(dup, title="Duplication ROC")
dev.off()

# Plot PRCs
pdf(paste(out.prefix, "del.prc.pdf", sep="."), height=4, width=4)
plot.prc(del, title="Deletion PRC")
dev.off()
pdf(paste(out.prefix, "dup.prc.pdf", sep="."), height=4, width=4)
plot.prc(dup, title="Duplication PRC")
dev.off()

