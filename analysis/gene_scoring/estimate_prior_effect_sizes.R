#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute prior effect sizes for rCNV gene scoring


# Load libraries
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(vioplot, quietly=T)
require(beeswarm, quietly=T)


# Set global params
options(stringsAsFactors=F, scipen=1000)


#################
### FUNCTIONS ###
#################
# Compute table of empirical priors
calc.priors <- function(meta.dat, cnv, neg, constrained, gs, xlist=c()){
  neg.or <- gene.meta.otf(meta.dat, setdiff(neg, xlist), cnv)[1:3]
  constr.or <- gene.meta.otf(meta.dat, setdiff(constrained, xlist), cnv)[1:3]
  gs.or <- gene.meta.otf(meta.dat, setdiff(gs, xlist), cnv)[1:3]
  res <- as.data.frame(cbind(c("theta0", "constr", "theta1"),
                      t(data.frame(neg.or, constr.or, gs.or))))
  rownames(res) <- NULL
  colnames(res) <- c("geneset", "lnOR", "lnOR_lower", "lnOR_upper")
  res[, -1] <- apply(res[, -1], 2, as.numeric)
  return(res)
}

# Plot prior effect size by gene set
plot.prior.lnors <- function(dat, cnv, title=NULL, ylims=NULL){
  # Get plot data
  par(mar=c(3, 3.5, 2, 1))
  colors <- c("gray75", control.cnv.colors[cnv], cnv.colors[cnv])
  if(is.null(ylims)){
    ylims <- range(dat[, -1], na.rm=T)
  }

  # Plot effect sizes
  plot(NA, xlim=c(0, 3), ylim=ylims, ylab="", xlab="", xaxt="n", las=2)
  abline(h=axTicks(2), lty=2, col="gray80")
  mtext(3, font=2, line=0.25, text=title)
  xlabs <- c("Gold-Standard\nTolerant Genes", "Constrained\nGenes",
             "Gold-Standard\nIntolerant Genes")
  sapply(1:length(xlabs), function(x){
    axis(1, at=x-0.5, tick=F, labels=xlabs[x])
  })
  mtext(2, line=2.5, text="ln(Odds Ratio)")
  segments(x0=(1:3)-0.5, x1=(1:3)-0.5, y0=dat$lnOR_lower, y1=dat$lnOR_upper,
           lwd=3, col=colors)
  points(x=(1:3)-0.5, y=dat$lnOR, pch=18, cex=2, col=colors)
}


#####################
### RSCRIPT BLOCK ###
#####################
# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog meta_inputs.tsv exclude.bed constrained.genes.list ",
                                            "del.gs_true.genes.list del.gs_false.genes.list",
                                            "dup.gs_true.genes.list dup.gs_false.genes.list",
                                            "out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 8){
  stop(paste("Eight positional arguments: meta_inputs.tsv, exclude.bed,",
             "constrained.genes.list, del.gs_true.genes.list, del.gs_false.genes.list",
             "dup.gs_true.genes.list dup.gs_false.genes.list out.prefix\n", sep=" "))
}

# Writes args & opts to vars
meta.inputs.in <- args$args[1]
excludelist.in <- args$args[2]
constrained.in <- args$args[3]
del.gs_true.in <- args$args[4]
del.gs_false.in <- args$args[5]
dup.gs_true.in <- args$args[6]
dup.gs_false.in <- args$args[7]
out.prefix <- args$args[8]

# # DEV PARAMTERS
# setwd("~/scratch")
# meta.inputs.in <- "prior_estimation.meta_inputs.tsv"
# excludelist.in <- "~/scratch/rCNV.gene_scoring.training_gene_excludelist.bed.gz"
# constrained.in <- "~/scratch/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
# del.gs_true.in <- "~/scratch/gold_standard.haploinsufficient.genes.list"
# del.gs_false.in <- "~/scratch/gold_standard.haplosufficient.genes.list"
# dup.gs_true.in <- "~/scratch/gold_standard.triplosensitive.genes.list"
# dup.gs_false.in <- "~/scratch/gold_standard.triploinsensitive.genes.list"
# out.prefix <- "rCNV_prior_estimation"

# Read all gene lists
xlist <- read.table(excludelist.in, header=F)[, 4]
constrained <- read.table(constrained.in, header=F)[, 1]
del.gs <- read.table(del.gs_true.in, header=F)[, 1]
del.neg <- read.table(del.gs_false.in, header=F)[, 1]
dup.gs <- read.table(dup.gs_true.in, header=F)[, 1]
dup.neg <- read.table(dup.gs_false.in, header=F)[, 1]

# Load meta-analysis data
meta.dat <- load.otf.meta.dat(meta.inputs.in)

# Compute prior effectsizes
del.priors <- calc.priors(meta.dat, "DEL", del.neg, constrained, del.gs, xlist)
dup.priors <- calc.priors(meta.dat, "DUP", dup.neg, constrained, dup.gs, xlist)

# Plot effect sizes
ylims <- range(cbind(del.priors[, -1], dup.priors[, -1]), na.rm=T)
pdf(paste(out.prefix, "observed_effect_sizes.pdf", sep="."), height=4, width=10)
par(mfrow=c(1, 2))
plot.prior.lnors(del.priors, "DEL", title="Deletions", ylims=ylims)
plot.prior.lnors(dup.priors, "DUP", title="Duplications", ylims=ylims)
dev.off()

# Write table of priors
del.priors$cnv <- "DEL"
dup.priors$cnv <- "DUP"
out.df <- rbind(del.priors, dup.priors)[, c("geneset", "cnv", "lnOR")]
colnames(out.df) <- c("#parameter", "cnv", "value")
write.table(out.df[which(out.df[, 1] != "constr"), ],
            paste(out.prefix, "empirical_prior_estimates.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

