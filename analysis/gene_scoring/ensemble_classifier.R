#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Ensemble classification of gene scores


options(stringsAsFactors=F, scipen=1000)


#################
### FUNCTIONS ###
#################
# Load a single gene score stats .tsv
load.scores.single <- function(path, prefix=NULL){
  x <- read.table(path, header=T, sep="\t", comment.char="")[, 1:3]
  colnames(x)[1] <- "gene"
  if(!is.null(prefix)){
    colnames(x)[-1] <- paste(prefix, colnames(x)[-1], sep=".")
  }
  x[, 2:3] <- apply(x[, 2:3], 2, as.numeric)
  return(x)
}

# Load all gene score stats from list
load.scores <- function(scores_file.in){
  sl <- read.table(scores_file.in)[, 1]
  models <- as.vector(sapply(sl, function(path){
    parts <- unlist(strsplit(path, split=".", fixed=T))
    as.character(parts[length(parts)-1])
  }))
  scores <- load.scores.single(sl[1], models[1])
  for(i in 2:length(sl)){
    scores <- merge(scores,
                    load.scores.single(sl[i], models[i]),
                    by="gene", all=T, sort=F)
  }
  return(list("scores"=scores, "models"=models))
}

# Compute ROC
roc <- function(stats, score, truth.genes, neg.genes, steps=seq(1, 0, -0.001), eq="gt"){
  x <- data.frame("score" = stats[, which(colnames(stats) == score)],
                  "true" = stats$gene %in% truth.genes,
                  "neg" = stats$gene %in% neg.genes)
  roc_res <- as.data.frame(t(sapply(steps, function(k){
    if(eq=="gt"){
      idxs <- which(x$score > k)
    }else{
      idxs <- which(x$score < k)
    }
    n.true <- length(which(x$true[idxs]))
    ftrue <- n.true / length(which(x$true))
    n.neg <- length(which(x$neg[idxs]))
    fneg <- n.neg / length(which(x$neg))
    fother <- length(which(!x$true[idxs])) / length(which(!x$true))
    fall <- length(idxs) / nrow(x)
    acc <- n.true / (n.true + n.neg)
    d.opt <- sqrt((0-fneg)^2 + (1-ftrue)^2)
    return(c(k, fall, fother, ftrue, fneg, acc, d.opt))
  })))
  colnames(roc_res) <- c("cutoff", "frac_all", "frac_other", "frac_true", "frac_neg", "accuracy", "d.opt")
  roc.opt.idx <- head(which(roc_res$d.opt == min(roc_res$d.opt, na.rm=T)), 1)
  return(list("roc.res"=roc_res, 
              "roc.opt"=c("cutoff"=roc_res$cutoff[roc.opt.idx], 
                          "accuracy"=roc_res$accuracy[roc.opt.idx])))
}
  
# Calculate weights for each model
get.weights <- function(scores, models, pos_genes, neg_genes, eval.metric="pred_bfdp", eq="lt"){
  weights <- sapply(models, function(model){
    roc.res <- roc(scores, score=paste(model, eval.metric, sep="."),
                   pos_genes, neg_genes, eq=eq)
    as.numeric(roc.res$roc.opt[2])
  })
  weights <- weights / sum(weights)
  names(weights) <- models
  return(weights)
}

# Compute weighted average of BFDPs for ensemble classifier
get.ensemble.bfdps <- function(scores, models, weights){
  weighted.bfdps <- sapply(models, function(model){
    w <- weights[which(names(weights)==model)]
    bfdps <- scores[, which(colnames(scores) == paste(model, "pred_bfdp", sep="."))]
    w * bfdps
  })
  apply(weighted.bfdps, 1, sum)
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
args <- parse_args(OptionParser(usage="%prog scores.tsv truth.genes false.genes out.tsv",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Four positional arguments: scores.tsv, truth.genes, false.genes, and out.tsv\n")
}

# Writes args & opts to vars
scores_file.in <- args$args[1]
pos_genes.in <- args$args[2]
neg_genes.in <- args$args[3]
outfile <- args$args[4]

# # DEV PARAMTERS
# setwd("~/scratch")
# scores_file.in <- "~/scratch/ensemble_input.tsv"
# pos_genes.in <- "~/scratch/gold_standard.haploinsufficient.genes.list"
# neg_genes.in <- "~/scratch/gold_standard.haplosufficient.genes.list"
# outfile <- "~/scratch/ensemble_scores.test.tsv"

# Read gene scores
scores <- load.scores(scores_file.in)
models <- scores$models
scores <- scores$scores

# Read truth sets
pos_genes <- read.table(pos_genes.in)[, 1]
neg_genes <- read.table(neg_genes.in)[, 1]

# Determine weights for each model
weights <- get.weights(scores, models, pos_genes, neg_genes)

# Compute ensemble BFDPs, normalized scores, and percentile
scores$ensemble.pred_bfdp <- get.ensemble.bfdps(scores, models, weights)
scores$ensemble.score <- pnorm(scale(scores$ensemble.pred_bfdp, center=T, scale=T), lower.tail=F)
scores$ensemble.quantile <- 100 * order(scores$ensemble.pred_bfdp) / nrow(scores)

# Reformat output file and write out
out.df <- data.frame("gene"=scores$gene,
                     "pred_bfdp"=scores$ensemble.pred_bfdp,
                     "score"=scores$ensemble.score,
                     "quantile"=scores$ensemble.quantile)
write.table(out.df, outfile, col.names=T, row.names=F, sep="\t", quote=F)
