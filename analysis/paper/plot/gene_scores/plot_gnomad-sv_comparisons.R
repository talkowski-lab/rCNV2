#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compare rates of genic SVs in gnomAD-SV vs rCNV gene scores


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Predict number of SVs expected per gene in gnomAD-SV
calc.expected <- function(var.meta, genomic.meta, train.genes=NULL){
  if(is.null(train.genes)){
    train.genes <- unique(genomic.meta$gene)
  }
  csqs <- colnames(var.meta)[grep("gnomad_sv", colnames(var.meta), fixed=T)]
  exp <- as.data.frame(sapply(csqs, function(csq){
    all.dat <- merge(var.meta[, c("gene", csq)], genomic.meta, all.x=T, all.y=F, sort=F, by="gene")
    rownames(all.dat) <- all.dat$gene
    all.dat$gene <- NULL
    colnames(all.dat)[1] <- "sv"
    tdat <- all.dat[which(rownames(all.dat) %in% train.genes), ]
    fit <- glm.nb(sv ~ ., data=tdat)
    predict.glm(fit, newdata=all.dat[, -1], type="response")
  }))
  colnames(exp) <- paste("exp", csqs, sep=".")
  exp$gene <- rownames(exp)
  merge(var.meta, exp, by="gene", all.x=T, all.y=F, sort=F)
}

# Helper functions for boot (bootstrapping confidence intervals)
summed.oe <- function(obs.exp.df, indexes=NULL){
  if(is.null(indexes)){
    indexes <- 1:nrow(obs.exp.df)
  }
  sum.obs <- sum(obs.exp.df[indexes, 1])
  sum.exp <- sum(obs.exp.df[indexes, 2])
  sum.obs / sum.exp
}
mean.for.boot <- function(vals, indexes){mean(vals[indexes], na.rm=T)}

# Compute obs/exp SVs in gnomAD-SV across bins of rCNV score
gnomad.by.score.bin <- function(oe.dat, scores, score="pHI", 
                                var="gnomad_sv_lof_del", n.bins=10){
  dat <- merge(scores, oe.dat, all=F, sort=F, by="gene")
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
  oe.vals <- t(sapply(1:n.bins, function(i){
    idxs <- which(x>=bins[i] & x<=bins[i+1])
    obs.exp.df <- data.frame(y.obs[idxs], y.exp[idxs])
    oe <- summed.oe(obs.exp.df)
    ci <- boot.ci(boot(data=obs.exp.df, statistic=summed.oe, R=1000), 
                  conf=0.95, type="norm")$normal[, -1]
    as.numeric(c(oe, ci))
  }))
  colnames(oe.vals) <- c("oe", "oe.lower", "oe.upper")
  avgs <- t(sapply(1:n.bins, function(i){
    idxs <- which(x>=bins[i] & x<=bins[i+1])
    estimate <- mean(y.obs[idxs])
    ci <- boot.ci(boot(data=y.obs[idxs], statistic=mean.for.boot, R=1000), 
                  conf=0.95, type="norm")$normal[, -1]
    as.numeric(c(estimate, ci))
  }))
  colnames(avgs) <- c("mean", "mean.lower", "mean.upper")
  cbind(data.frame("score"=mids), avgs, oe.vals)
}

# Compute obs/exp SVs in gnomAD-SV, stratified by low/high pHI & pTS
gnomad.stratified <- function(oe.dat, scores, high.cutoff=0.8, low.cutoff=0.5, measure="oe"){
  dat <- merge(scores, oe.dat, all=F, sort=F, by="gene")
  dat <- dat[which(!is.na(dat$pHI) & !(is.na(dat$pTS))), ]
  low.low.idx <- which(dat$pHI<=low.cutoff & dat$pTS<=low.cutoff)
  high.low.idx <- which(dat$pHI>=high.cutoff & dat$pTS<=low.cutoff)
  low.high.idx <- which(dat$pHI<=low.cutoff & dat$pTS>=high.cutoff)
  high.high.idx <- which(dat$pHI>=high.cutoff & dat$pTS>=high.cutoff)
  idxs <- list(low.low.idx, high.low.idx, low.high.idx, high.high.idx)
  strat.vals <- lapply(c("gnomad_sv_lof_del", "gnomad_sv_cg"), function(var){
    as.data.frame(do.call("rbind", lapply(idxs, function(idxs.i){
      y.obs <- as.numeric(dat[idxs.i, var])
      y.exp <- as.numeric(dat[idxs.i, paste("exp", var, sep=".")])
      drop.idx <- which(is.na(y.obs) | is.na(y.exp))
      if(length(drop.idx) > 0){
        y.obs <- y.obs[-drop.idx]
        y.exp <- y.exp[-drop.idx]
      }
      if(measure=="oe"){
        obs.exp.df <- data.frame(y.obs, y.exp)
        estimate <- summed.oe(obs.exp.df)
        ci <- boot.ci(boot(data=obs.exp.df, statistic=summed.oe, R=1000), 
                      conf=0.95, type="norm")$normal[, -1]
      }else{
        estimate <- mean(y.obs)
        ci <- boot.ci(boot(data=y.obs, statistic=mean.for.boot, R=1000), 
                      conf=0.95, type="norm")$normal[, -1]
      }
      return(as.numeric(c(estimate, ci)))
    })))
  })
  names(strat.vals) <- c("LoF", "CG")
  return(strat.vals)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot a single comparison of rCNV score vs. gnomAD-SV counts
plot.oe <- function(oe.dat, scores, score, var, metric="oe", n.bins=10,
                    xlab=NULL, parmar=c(2.25, 4, 0.5, 0.5)){
  # Collect plot data
  plot.dat <- gnomad.by.score.bin(oe.dat, scores, score, var, n.bins)
  if(score=="pHI"){
    color <- cnv.colors[1]
    descrip <- "LoF Dels."
  }else if(score=="pTS"){
    color <- cnv.colors[2]
    descrip <- "CG Dups."
  }else{
    color <- blueblack
    descrip <- "SVs"
  }
  if(metric=="oe"){
    x <- plot.dat[, c(1, grep("oe", colnames(plot.dat), fixed=T))]
    ylab <- paste("Obs/Exp ", descrip, "\nin gnomAD-SV", sep="")
  }else{
    x <- plot.dat[, c(1, grep("mean", colnames(plot.dat), fixed=T))]
    ylab <- paste(descrip, "per Gene\nin gnomAD-SV")
  }
  ylims <- range(x[, -1], na.rm=T)
  if(is.null(xlab)){
    xlab <- paste("Genes Binned by", score)
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 1), ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, bty="n", col=bluewhite)
  abline(h=axTicks(2), v=axTicks(1), col="white")
  if(metric=="oe"){
    abline(h=1, col=blueblack, lty=2)
  }
  
  # Add points
  segments(x0=x[, 1], x1=x[, 1], y0=x[, 3], y1=x[, 4],
           col=color, lend="round", lwd=1.5)
  points(x=x[, 1], y=x[, 2], pch=19, col=color)
  
  # Add axes
  x.ax.at <- axTicks(1)
  axis(1, at=c(-100, 100), col=blueblack, tck=0, labels=NA)
  axis(1, at=x.ax.at, labels=NA, tck=-0.025, col=blueblack)
  sapply(x.ax.at, function(k){
    axis(1, at=k, tick=F, line=-0.75)
  })
  mtext(1, text=xlab, line=1.1)
  y.ax.at <- axTicks(2)
  axis(2, at=c(-100, 100), col=blueblack, tck=0, labels=NA)
  axis(2, at=y.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(2, at=y.ax.at, tick=F, line=-0.65, las=2)
  mtext(2, text=ylab, line=2)
}

# Plot rates of SVs in bins by CNV score (either percentile or absolute bin)
plot.oe.by.scorebin <- function(oe.dat, scores, score, var, metric="oe", n.bins=10,
                                x.label=NULL, parse.x.label=T, x.label.line=2,
                                x.ax.labels=NULL, parse.x.ax.labels=T, ylims=NULL,
                                cex.x.ax.labels=1, pt.color=blueblack, 
                                baseline.color=blueblack, null.color=bluewhite, 
                                parmar=c(3.5, 3.5, 0.5, 0.5)){
  # Collect plot data
  all.plot.dat <- gnomad.by.score.bin(oe.dat, scores, score, var, n.bins)
  if(metric == "oe"){
    baseline <- 1
    plot.dat <- all.plot.dat[, c("score", "oe", "oe.lower", "oe.upper")]
  }else{
    baseline <- NA
    plot.dat <- all.plot.dat[, c("score", "mean", "mean.lower", "mean.upper")]
  }
  colnames(plot.dat) <- c("score", "estimate", "lower", "upper")
  if(is.null(ylims)){
    ylims <- range(plot.dat[, -1], na.rm=T)
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, n.bins), ylim=ylims, xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, bty="n", col=bluewhite)
  y.ax.at <- sort(unique(round(axTicks(2), 1)))
  abline(h=y.ax.at, col="white")
  abline(h=c(0, baseline), lty=c(1, 2), col=c(null.color, baseline.color))
  text(x=par("usr")[1]-(0.05*(par("usr")[2]-par("usr")[1])),
       y=baseline+(0.035*(par("usr")[4]-par("usr")[3])), 
       labels="Expected", cex=0.85, font=3, col=baseline.color, pos=4)
  
  # Add points
  pt.x.at <- (1:n.bins)-0.5
  segments(x0=pt.x.at, x1=pt.x.at, y0=plot.dat$lower, y1=plot.dat$upper, lend="round", 
           lwd=2.5, col=pt.color)
  points(x=pt.x.at, y=plot.dat$estimate, pch=19, col=pt.color, cex=1.5)
  
  # Add X axis
  if(!is.null(x.ax.labels)){
    sapply(1:length(pt.x.at), function(i){
      if(parse.x.ax.labels==T){
        axis(1, at=pt.x.at[i], tick=F, line=-1, labels=parse(text=x.ax.labels[i]), 
             cex.axis=cex.x.ax.labels)
      }else{
        axis(1, at=pt.x.at[i], tick=F, line=-1, labels=x.ax.labels[i], 
             cex.axis=cex.x.ax.labels)
      }
    })
  }
  if(parse.x.label==T){
    mtext(1, text=parse(text=x.label), line=x.label.line)
  }else{
    mtext(1, text=x.label, line=x.label.line)
  }
  
  # Add Y axis
  axis(2, at=c(-100, 100), tck=0, col=blueblack, labels=NA)
  axis(2, at=y.ax.at, col=blueblack, labels=NA, tck=-0.025)
  axis(2, at=y.ax.at, tick=F, labels=y.ax.at, las=2, line=-0.6)
  if(metric == "oe"){
    mtext(2, text="Obs/Exp CNVs per Gene", line=1.75)
  }else{
    mtext(2, text="Mean CNVs per Gene", line=1.75)
  }
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(MASS, quietly=T)
require(boot, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv variation.meta.bed genomic.meta.eigen.bed train.genes out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: scores.tsv, variation.meta.bed, genomic.meta.eigen.bed, training.genes, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
var.meta.in <- args$args[2]
genomic.meta.in <- args$args[3]
mu.training.genes.in <- args$args[4]
out.prefix <- args$args[5]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# setwd("~/scratch/")
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# var.meta.in <- "gencode.v19.canonical.pext_filtered.variation_features.bed.gz"
# genomic.meta.in <- "gencode.v19.canonical.pext_filtered.genomic_features.eigenfeatures.bed.gz"
# mu.training.genes.in <- "gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list"
# out.prefix <- "gene_score_gnomad_comparisons"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load scores
scores <- load.scores(scores.in)

# Load gene metadata
var.meta <- load.gene.metadata(var.meta.in)
genomic.meta <- load.gene.metadata(genomic.meta.in)

# Train mutation rate model to predict expected number of variants per gene
train.genes <- as.character(read.table(mu.training.genes.in)[, 1])
oe.dat <- calc.expected(var.meta, genomic.meta, train.genes)

# Plot rCNV scores vs obs/exp ratios of gnomAD SVs
pdf(paste(out.prefix, "scores_vs_gnomAD-SV.pHI_oe.pdf", sep="."),
    height=2, width=2.6)
plot.oe(oe.dat, scores, "pHI", "gnomad_sv_lof_del", n.bins=10)
dev.off()
pdf(paste(out.prefix, "scores_vs_gnomAD-SV.pTS_oe.pdf", sep="."),
    height=2, width=2.6)
plot.oe(oe.dat, scores, "pTS", "gnomad_sv_cg", n.bins=10)
dev.off()

# Plot rCNV scores vs mean SVs per gene from gnomAD-SV
pdf(paste(out.prefix, "scores_vs_gnomAD-SV.pHI_raw.pdf", sep="."),
    height=2, width=2.6)
plot.oe(oe.dat, scores, "pHI", "gnomad_sv_lof_del", metric="raw", n.bins=10)
dev.off()
pdf(paste(out.prefix, "scores_vs_gnomAD-SV.pTS_raw.pdf", sep="."),
    height=2, width=2.6)
plot.oe(oe.dat, scores, "pTS", "gnomad_sv_cg", metric="raw", n.bins=10)
dev.off()

# Plot gnomAD-SV stratified by high/low pHI & pTS
strat.vals <- gnomad.stratified(oe.dat, scores, high.cutoff=0.9, low.cutoff=0.5, measure="oe")
pdf(paste(out.prefix, "scores_vs_gnomAD-SV.stratified.pdf", sep="."),
    height=2.25, width=3)
plot.stratified.metric(strat.vals, y.title="\"gnomAD-SV Obs/Exp\"")
dev.off()

# Plot large quintile-binned panel of pHI for main figure
pdf(paste(out.prefix, "scores_vs_gnomAD-SV.pHI_raw.quintiled_large.pdf", sep="."),
    height=2.75, width=2.1)
plot.oe.by.scorebin(oe.dat, scores, "pHI", "gnomad_sv_lof_del", metric="mean", n.bins=5,
                    x.label="All Autosomal Genes", 
                    parse.x.label=F, x.label.line=1.1,
                    x.ax.labels=paste("\"Q\"[", 1:5, "]", sep=""),
                    parse.x.ax.labels=T, cex.x.ax.labels=0.9, 
                    pt.color=cnv.colors[1], baseline.color=NA, 
                    null.color=blueblack, parmar=c(3.25, 2.75, 0.5, 0.5))
mtext(1, line=2.1, text="in Quintiles by pHI")
dev.off()

# Plot large quintile-binned panel of pHI for main figure
pdf(paste(out.prefix, "scores_vs_gnomAD-SV.pTS_raw.quintiled_large.pdf", sep="."),
    height=2.75, width=2.1)
plot.oe.by.scorebin(oe.dat, scores, "pTS", "gnomad_sv_cg", metric="mean", n.bins=5,
                    x.label="All Autosomal Genes", 
                    parse.x.label=F, x.label.line=1.1,
                    x.ax.labels=paste("\"Q\"[", 1:5, "]", sep=""),
                    parse.x.ax.labels=T, cex.x.ax.labels=0.9, 
                    pt.color=cnv.colors[2], baseline.color=NA,
                    null.color=blueblack, parmar=c(3.25, 2.75, 0.5, 0.5))
mtext(1, line=2.1, text="in Quintiles by pTS")
dev.off()
