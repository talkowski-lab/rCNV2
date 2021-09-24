#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Calculate mean effect size variance versus number of CNVs 


options(stringsAsFactors=F, scipen=1000)


#################
### FUNCTIONS ###
#################
# Compute standard error, standard deviation, and variance given 95% CI, mean, and sample size (number of CNVs)
calc.variance <- function(mean, ci.upper, cnvs){
  se <- (ci.upper - mean) / 1.96
  sd <- sqrt(cnvs) * se
  return(c(se, sd, sd ^ 2))
}

# Load meta-analysis stats, compute variance and merge with counts
load.data <- function(stats.in, counts.in){
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")[, -c(1:3)]
  counts <- read.table(counts.in, header=T, sep="\t", comment.char="")
  colnames(counts)[1] <- "gene"
  stats <- merge(stats, counts, by="gene", all=T, sort=F)
  vars <- as.data.frame(t(sapply(1:nrow(stats), function(i){
    calc.variance(stats$meta_lnOR[i], stats$meta_lnOR_upper[i], stats$cnvs[i])
  })))
  colnames(vars) <- c("se", "sd", "var")
  cbind(stats, vars)
}

# Bin genes by minimum CNV count and plot average binwise variance
plot.var.by.cnvs <- function(stats, metric="var", max.eval=30){
  x <- as.data.frame(t(sapply(1:max.eval, function(i){
    c(i, mean(stats[which(stats$cnvs>=i), metric], na.rm=T))
  })))
  ylims <- c(0, max(x[, 2], na.rm=T))
  if(metric=="var"){
    ylab <- "Mean variance"
  }else if(metric=="se"){
    ylab <- "Mean standard error"
  }else if(metric=="sd"){
    ylab <- "Mean standard deviation"
  }
  par(mar=c(4, 4, 1, 1))
  plot(x, pch=19, col="grey30", xlab="Min. CNVs per Gene", ylab=ylab, ylim=ylims)
  abline(h=1, lty=2)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog stats.bed counts.tsv out.pdf",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Three positional arguments: stats.bed, counts.tsv, and out.pdf\n")
}

# Writes args & opts to vars
stats.in <- args$args[1]
counts.in <- args$args[2]
outfile <- args$args[3]

# # DEV PARAMTERS
# setwd("~/scratch")
# stats.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
# counts.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.counts_per_gene.tsv"
# outfile <- "~/scratch/variance_vs_counts.test.pdf"

# Load data
stats <- load.data(stats.in, counts.in)

# Plot variance vs number of genes
pdf(outfile, height=3, width=4)
plot.var.by.cnvs(stats, metric="se")
dev.off()
