#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compare gene scores to population CNV data from gnomAD


options(stringsAsFactors=F, scipen=1000)
colors <- c("#D43925", "#2376B2", "#7459B2")


#################
### FUNCTIONS ###
#################
# Load gene score stats .tsv and merge with variation metadata
load.scores <- function(scores.in, vardat.in){
  scores <- read.table(scores.in, header=T, sep="\t", comment.char="")[, 1:3]
  colnames(scores)[1] <- "gene"
  scores[, 2:3] <- apply(scores[, 2:3], 2, as.numeric)
  vardat <- read.table(vardat.in, header=T, sep="\t", comment.char="")[, -c(1:3)]
  vardat[, -1] <- apply(vardat[, -1], 2, as.numeric)
  merge(scores, vardat, by="gene", all.x=T, all.y=F, sort=F)
}

# Load covariates
load.covariates <- function(cov.in){
  cov <- read.table(cov.in, header=T, sep="\t", comment.char="")[, -c(1:3)]
  cov[, -1] <- apply(cov[, -1], 2, as.numeric)
  return(cov)
}

# Predict number of SVs expected per gene in gnomAD-SV
calc.expected <- function(dat, cov){
  csqs <- colnames(dat)[grep("gnomad_sv", colnames(dat), fixed=T)]
  exp <- as.data.frame(sapply(csqs, function(csq){
    tdat <- merge(dat[, c("gene", csq)], cov, all.x=T, all.y=F, sort=F, by="gene")
    rownames(tdat) <- tdat$gene
    tdat$gene <- NULL
    colnames(tdat)[1] <- "sv"
    fit <- glm.nb(sv ~ ., data=tdat)
    predict.glm(fit, newdata=tdat[, -1], type="response")
  }))
  colnames(exp) <- paste("exp", csqs, sep=".")
  exp$gene <- rownames(exp)
  merge(dat, exp, by="gene", all.x=T, all.y=F, sort=F)
}


# Compute obs/exp SVs in gnomAD-SV across bins of rCNV score
gnomad.by.score.bin <- function(dat, score="pHI", var="gnomad_sv_lof_del", n.bins=10){
  x <- as.numeric(dat[, score])
  y.obs <- as.numeric(dat[, var])
  y.exp <- as.numeric(dat[, paste("exp", var, sep=".")])
  drop.idx <- which(is.na(x) | is.na(y.obs) | is.na(y.exp))
  if(length(drop.idx) > 0){
    x <- x[-drop.idx]
    y.obs <- y.obs[-drop.idx]
    y.exp <- y.exp[-drop.idx]
  }
  bins <- seq(0, 1, length.out=n.bins+1)
  mids <- (bins[1:n.bins] + bins[2:(n.bins+1)])/2
  oe <- sapply(1:n.bins, function(i){
    idxs <- which(x>=bins[i] & x<=bins[i+1])
    sum.obs <- sum(y.obs[idxs])
    sum.exp <- sum(y.exp[idxs])
    sum.obs / sum.exp
  })
  avgs <- sapply(1:n.bins, function(i){
    idxs <- which(x>=bins[i] & x<=bins[i+1])
    mean(y.obs[idxs])
  })
  return(data.frame("score"=mids, "mean"=avgs, "oe"=oe))
}

# Plot pTS & pHI for a single gnomAD-SV metric
plot.oe <- function(dat, var, metric="oe", n.bins=10, title=NULL){
  phi <- gnomad.by.score.bin(dat, "pHI", var, n.bins)
  pts <- gnomad.by.score.bin(dat, "pTS", var, n.bins)
  if(metric=="oe"){
    y.del <- phi$oe
    y.dup <- pts$oe
    ylab <- "Obs/Exp in gnomAD-SV"
  }else{
    y.del <- phi$mean
    y.dup <- pts$mean
    ylab <- "Mean SVs per Gene"
  }
  ylims <- range(c(y.del, y.dup))
  plot(NA, xlim=c(0, 1), ylim=ylims,
       xlab="Score Bin (pHI or pTS)", ylab=ylab, main=title)
  if(metric=="oe"){
    abline(h=1, lty=2, col="gray40")
  }
  points(x=phi$score-0.01, y=y.del, type="l", col=colors[1])
  points(x=pts$score+0.01, y=y.dup, type="l", col=colors[2])
  points(x=phi$score-0.01, y=y.del, pch=19, col=colors[1])
  points(x=pts$score+0.01, y=y.dup, pch=19, col=colors[2])
  legend("bottomleft", legend=c("pHI", "pTS"), fill=colors[1:2])
}




#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(MASS, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog scores.tsv variation_data.bed gene_covariants.bed out.pdf",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Four positional arguments: scores.tsv, variation_data.bed, gene_covariates.bed, and out.pdf\n")
}

# Writes args & opts to vars
scores.in <- args$args[1]
vardat.in <- args$args[2]
cov.in <- args$args[3]
outfile <- args$args[4]

# # DEV PARAMTERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# vardat.in <- "~/scratch/gencode.v19.canonical.pext_filtered.variation_features.bed.gz"
# cov.in <- "~/scratch/gencode.v19.canonical.pext_filtered.genomic_features.eigenfeatures.bed.gz"
# outfile <- "~/scratch/test_gnomad_cor.pdf"

# Read gene scores and merge with variation metadata
dat <- load.scores(scores.in, vardat.in)

# Load gene covariates and predict expected # of SVs per gene
cov <- load.covariates(cov.in)
dat <- calc.expected(dat, cov)

# Plot enrichments
pdf(outfile, height=6, width=8)
par(mfrow=c(2, 2))
plot.oe(dat, "gnomad_sv_lof_del", title="LoF Deletions in gnomAD-SV", n.bins=20)
plot.oe(dat, "gnomad_sv_cg", title="CG Duplications in gnomAD-SV", n.bins=20)
plot.oe(dat, "gnomad_sv_lof_del", title="LoF Deletions in gnomAD-SV", n.bins=20, metric="mean")
plot.oe(dat, "gnomad_sv_cg", title="CG Duplications in gnomAD-SV", n.bins=20, metric="mean")
dev.off()

