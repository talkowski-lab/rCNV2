#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute meta-analysis effect size as a function of number of CNVs


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############

# Load meta analysis stats, add CNV counts, and simplify
load.data <- function(stats.in, counts.in){
  st <- read.table(stats.in, header=T, sep="\t", comment.char="")
  ct <- read.table(counts.in, header=T, sep="\t", comment.char="")
  colnames(ct)[1] <- "gene"
  x <- merge(st, ct, all=T, by="gene", sort=F)
  x$meta_lnOR_ci_width <- x$meta_lnOR_upper - x$meta_lnOR_lower
  x[, c("gene", "meta_lnOR", "meta_lnOR_ci_width", "cnvs")]
}

# Compute sd of effect sizes as function of CNV count
sd.by.cnv <- function(dat, max.eval=50){
  x <- as.data.frame(t(sapply(0:max.eval, function(min.cnv){
    idxs <- which(dat$cnvs >= min.cnv)
    c(min.cnv, sd(dat$meta_lnOR[idxs], na.rm=T))
  })))
  colnames(x) <- c("min.cnv", "sd")
  return(x)
}

# Compute average CI width of effect sizes as function of CNV count
ci.width.by.cnv <- function(dat, max.eval=50){
  x <- as.data.frame(t(sapply(0:max.eval, function(min.cnv){
    idxs <- which(dat$cnvs >= min.cnv)
    c(min.cnv, mean(dat$meta_lnOR_ci_width[idxs], na.rm=T))
  })))
  colnames(x) <- c("min.cnv", "mean.ci.width")
}

# Plot SD versus min CNVs for DEL and DUP
plot.sds <- function(dels, dups, min.sd){
  # Get plot values
  xlims <- range(c(dels[, 1], dups[, 1]), na.rm=T)
  ylims <- c(0, max(c(dels[, 2], dups[, 2]), na.rm=T))
  
  # Plot
  par(mfrow=c(1, 2), mar=c(4, 4, 2, 1), bty="n")
  plot(dels, col="firebrick", pch=19, xlim=xlims, ylim=ylims,
       xlab="Min. # of CNVs", ylab="Std. Dev. of ln(OR)", main="Deletions",
       panel.first=c(abline(h=min.sd, lty=2)))
  plot(dups, col="dodgerblue3", pch=19, xlim=xlims, ylim=ylims,
       xlab="Min. # of CNVs", ylab="Std. Dev. of ln(OR)", main="Duplications",
       panel.first=c(abline(h=min.sd, lty=2)))
}

################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--max-sd"), default=1,
              help="Maximum standard deviation of effect size to tolerate")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog del.stats.bed dup.stats.bed del.counts.tsv dup.counts.tsv outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
del.stats.in <- args$args[1]
dup.stats.in <- args$args[2]
del.counts.in <- args$args[3]
dup.counts.in <- args$args[4]
outfile <- args$args[5]
min.sd <- opts$`min-sd`

# Dev parameters
setwd("~/scratch/")
del.stats.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
dup.stats.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz"
del.counts.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.counts_per_gene.tsv"
dup.counts.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DUP.counts_per_gene.tsv"
outfile <- "~/scratch/test.sd_vs_cnvs.pdf"
min.sd <- 1

# Load data
del <- load.data(del.stats.in, del.counts.in)
dup <- load.data(dup.stats.in, dup.counts.in)

# Compute sd vs number of CNVs
del.sd <- sd.by.cnv(del)
dup.sd <- sd.by.cnv(dup)

# Plot results
pdf(outfile, height=4, width=8)
plot.sds(del.sd, dup.sd, min.sd)
dev.off()
