#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute prior effect sizes for rCNV gene scoring


# Load libraries
require(optparse, quietly=T)
require(vioplot, quietly=T)
require(beeswarm, quietly=T)


# Set global params
options(stringsAsFactors=F, scipen=1000)
colors <- c("#D43925", "#2376B2", "#7459B2")


#################
### FUNCTIONS ###
#################
# Read and merge all data
read.data <- function(path, blacklist, constrained, gs, neg){
  x <- read.table(path, header=T, sep="\t", comment.char="")
  prune <- which(x$gene %in% blacklist)
  if(length(prune) > 0){
    x <- x[-prune, ]
  }
  x$var <- (x$meta_lnOR/1.96)^2
  x$constrained <- x$gene %in% constrained
  x$gs <- x$gene %in% gs
  x$neg <- x$gene %in% neg
  return(x)
}

# Plot prior effect size by gene set
plot.prior.lnors <- function(dat, color, title, pct=1){
  par(mar=c(3, 3.5, 2, 1))
  cutoffs <- quantile(dat$meta_lnOR, probs=c(1-pct, pct), na.rm=T)
  
  # Plot effect sizes
  plot(x=c(-1, 3), y=range(dat$meta_lnOR, na.rm=T), type="n", 
       ylab="", xlab="", xaxt="n", las=2)
  abline(h=axTicks(2), lty=2, col="gray80")
  mtext(3, font=2, line=0.25, text=title)
  xlabs <- c("Gold-Standard\nTolerant Genes", "All\nGenes", 
             "Constrained\nGenes", "Gold-Standard\nIntolerant Genes")
  sapply(1:length(xlabs), function(x){
    axis(1, at=x-1.5, tick=F, labels=xlabs[x])
  })
  mtext(2, line=2, text="ln(Odds Ratio)")
  vioplot(dat$meta_lnOR[which(!is.na(dat$meta_lnOR) & dat$neg & dat$meta_lnOR<=cutoffs[2])],
          at=-0.5, add=T, h=0.25,
          col="gray75", drawRect=F)
  boxplot(dat$meta_lnOR[which(dat$neg & dat$meta_lnOR<=cutoffs[2])], at=-0.5, add=T, col=NA, outline=F, boxwex=1/3, 
          staplewex=1/3, lty=1, lwd=2, yaxt="n", col.border="white")
  vioplot(dat$meta_lnOR[which(!is.na(dat$meta_lnOR))], at=0.5, add=T, h=0.25,
          col="gray40", drawRect=F)
  boxplot(dat$meta_lnOR, at=0.5, add=T, col=NA, outline=F, boxwex=1/3, 
          staplewex=1/3, lty=1, lwd=2, yaxt="n", col.border="white")
  vioplot(dat$meta_lnOR[which(!is.na(dat$meta_lnOR) & dat$constrained)], 
          at=1.5, add=T, h=0.25, col=color, drawRect=F)
  boxplot(dat$meta_lnOR[which(!is.na(dat$meta_lnOR) & dat$constrained)], 
          at=1.5, add=T, col=NA, outline=F, boxwex=1/3, staplewex=1/3, lty=1, lwd=2, yaxt="n")
  beeswarm(dat$meta_lnOR[which(dat$gs & dat$meta_lnOR>=cutoffs[1])], add=T, at=2.5, pch=19, col=color,
           corral="wrap", corralwidth=0.8, cex=0.5)
  boxplot(dat$meta_lnOR[which(dat$gs & dat$meta_lnOR>=cutoffs[1])], at=2.5, add=T, col=NA, outline=F, 
          boxwex=1/3, staplewex=1/3, lty=1, lwd=2, yaxt="n")
  abline(h=median(dat$meta_lnOR, na.rm=T))
  
  # # Plot variance
  # plot(x=c(-1, 3), y=range(dat$var, na.rm=T), type="n", 
  #      ylab="", xlab="", xaxt="n", las=2)
  # abline(h=axTicks(2), lty=2, col="gray80")
  # mtext(3, font=2, line=0.25, text="Variance of ln(OR)")
  # axis(1, at=-0.5:2.5, tick=F, labels=c("Gold-Standard\nTolerant Genes", "All\nGenes", 
  #                                       "Constrained\nGenes", "Gold-Standard\nIntolerant Genes"))
  # mtext(2, line=2, text="Variance of ln(OR)")
  # vioplot(dat$var[which(!is.na(dat$var) & dat$neg)], at=-0.5, add=T, h=0.25,
  #         col="gray75", drawRect=F)
  # boxplot(dat$var[which(dat$neg)], at=-0.5, add=T, col=NA, outline=F, boxwex=1/3, 
  #         staplewex=1/3, lty=1, lwd=2, yaxt="n", col.border="white")
  # vioplot(dat$var[which(!is.na(dat$var))], at=0.5, add=T, h=0.25,
  #         col="gray40", drawRect=F)
  # boxplot(dat$var, at=0.5, add=T, col=NA, outline=F, boxwex=1/3, 
  #         staplewex=1/3, lty=1, lwd=2, yaxt="n", col.border="white")
  # vioplot(dat$var[which(!is.na(dat$var) & dat$constrained)], 
  #         at=1.5, add=T, h=0.25, col=color, drawRect=F)
  # boxplot(dat$var[which(!is.na(dat$var) & dat$constrained)], 
  #         at=1.5, add=T, col=NA, outline=F, boxwex=1/3, staplewex=1/3, lty=1, lwd=2, yaxt="n")
  # beeswarm(dat$var[which(dat$gs)], add=T, at=2.5, pch=19, col=color,
  #          corral="wrap", corralwidth=0.8, cex=0.5)
  # boxplot(dat$var[which(dat$gs)], at=2.5, add=T, col=NA, outline=F, 
  #         boxwex=1/3, staplewex=1/3, lty=1, lwd=2, yaxt="n")
  # abline(h=median(dat$var, na.rm=T))
}

# Compute table of empirical priors
calc.priors <- function(del, dup, pct=1){
  vals <- unlist(lapply(list(del, dup), function(df){
    cutoffs <- quantile(df$meta_lnOR, probs=c(1-pct, pct), na.rm=T)
    c(median(df$meta_lnOR[which(df$neg & df$meta_lnOR<=cutoffs[2])], na.rm=T),
      median(df$meta_lnOR[which(df$gs & df$meta_lnOR>=cutoffs[1])], na.rm=T),
      median(df$var[which(df$neg & df$meta_lnOR<=cutoffs[2])], na.rm=T),
      median(df$var[which(df$gs & df$meta_lnOR>=cutoffs[1])], na.rm=T))
  }))
  data.frame("parameter"=rep(c("theta0", "theta1", "var0", "var1"), 2), 
             "cnv"=c(rep("DEL", 4), rep("DUP", 4)), 
             "estimate"=vals)
}

# Plot H0/H1 curves
plot.example.priors <- function(del, dup, colors){
  # Get plot values  
  x <- seq(-10, 10, 0.05)
  var <- median(del$var[which(del$gs)], na.rm=T)
  theta0.del <- median(del$meta_lnOR[which(del$neg)], na.rm=T)
  del0 <- 1-sapply(x, function(k){pnorm(k, mean=theta0.del, sd=sqrt(var))})
  theta0.dup <- median(dup$meta_lnOR[which(dup$neg)], na.rm=T)
  dup0 <- 1-sapply(x, function(k){pnorm(k, mean=theta0.dup, sd=sqrt(var))})
  theta1 <- median(del$meta_lnOR[which(del$gs)], na.rm=T)
  alt <- sapply(x, function(k){pnorm(k, mean=theta1, sd=sqrt(var))})
  
  # Plot
  par(mar=c(3.5, 3.5, 2, 1), mfrow=c(1, 2))
  
  # DEL
  plot(x=range(del$meta_lnOR, na.rm=T), y=c(0, 1), type="n",
       ylab="", xlab="", xaxt="n", las=2)
  abline(h=axTicks(2), lty=2, col="gray80")
  mtext(3, font=2, line=0.25, text="Deletion H0/H1")
  axis(1)
  mtext(1, line=2.25, text="ln(Odds Ratio)")
  mtext(2, line=2.5, text="Probability")
  points(x, del0, type="l", col="gray50", lwd=3)
  points(x, alt, type="l", col=colors[1], lwd=3)
  legend("right", col=c("gray50", colors[1]), pch=NA, lwd=3, lty=1,
         legend=c("H0: Gene is not HI", "H1: Gene is HI"), cex=0.75,
         bg="white")
  
  # DUP
  plot(x=range(dup$meta_lnOR, na.rm=T), y=c(0, 1), type="n",
       ylab="", xlab="", xaxt="n", las=2)
  abline(h=axTicks(2), lty=2, col="gray80")
  mtext(3, font=2, line=0.25, text="Duplication H0/H1")
  axis(1)
  mtext(1, line=2.25, text="ln(Odds Ratio)")
  mtext(2, line=2.5, text="Probability")
  points(x, dup0, type="l", col="gray50", lwd=3)
  points(x, alt, type="l", col=colors[2], lwd=3)
  legend("right", col=c("gray50", colors[2]), pch=NA, lwd=3, lty=1,
         legend=c("H0: Gene is not TS", "H1: Gene is TS"), cex=0.75,
         bg="white")
}


#####################
### RSCRIPT BLOCK ###
#####################
# List of command-line options
option_list <- list(
  make_option(c("--pct"), default=0.75,
              help="Percentile cutoff for top & bottom")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog del.bed.gz dup.bed.gz del.exclude.list dup.exclude.list",
                                            "constrained.genes.list del.gs_true.genes.list del.gs_false.genes.list",
                                            "dup.gs_true.genes.list dup.gs_false.genes.list",
                                            "out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 10){
  stop(paste("Ten positional arguments: del.bed.gz, dup.bed.gz, del.exclude.list, dup.exclude.list,",
             "constrained.genes.list, del.gs_true.genes.list, del.gs_false.genes.list", 
             "dup.gs_true.genes.list dup.gs_false.genes.list out.prefix\n", sep=" "))
}

# Writes args & opts to vars
del.in <- args$args[1]
dup.in <- args$args[2]
del.blacklist.in <- args$args[3]
dup.blacklist.in <- args$args[4]
constrained.in <- args$args[5]
del.gs_true.in <- args$args[6]
del.gs_false.in <- args$args[7]
dup.gs_true.in <- args$args[8]
dup.gs_false.in <- args$args[9]
out.prefix <- args$args[10]
pct <- opts$pct

# # DEV PARAMTERS
# setwd("~/scratch")
# del.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
# dup.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz"
# del.blacklist.in <- "~/scratch/DEL.training_blacklist.genes.list"
# dup.blacklist.in <- "~/scratch/DUP.training_blacklist.genes.list"
# constrained.in <- "~/scratch/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
# del.gs_true.in <- "~/scratch/gold_standard.haploinsufficient.genes.list"
# del.gs_false.in <- "~/scratch/gold_standard.haplosufficient.genes.list"
# dup.gs_true.in <- "~/scratch/gold_standard.triplosensitive.genes.list"
# dup.gs_false.in <- "~/scratch/gold_standard.triploinsensitive.genes.list"
# out.prefix <- "rCNV_prior_estimation"
# pct <- 0.75

# Read all gene lists
del.blacklist <- read.table(del.blacklist.in, header=F)[, 1]
dup.blacklist <- read.table(del.blacklist.in, header=F)[, 1]
constrained <- read.table(constrained.in, header=F)[, 1]
del.gs <- read.table(del.gs_true.in, header=F)[, 1]
del.neg <- read.table(del.gs_false.in, header=F)[, 1]
dup.gs <- read.table(dup.gs_true.in, header=F)[, 1]
dup.neg <- read.table(dup.gs_false.in, header=F)[, 1]

# Read meta-analysis stats and annotate with gene lists
del <- read.data(del.in, del.blacklist, constrained, del.gs, del.neg)
dup <- read.data(dup.in, dup.blacklist, constrained, dup.gs, dup.neg)

# Plot effect sizes
pdf(paste(out.prefix, "observed_effect_sizes.pdf", sep="."), height=4, width=12)
par(mfrow=c(1, 2))
plot.prior.lnors(del, colors[1], title="Deletions", pct=pct)
plot.prior.lnors(dup, colors[2], title="Duplications", pct=pct)
dev.off()

# Compute table of priors
priors <- calc.priors(del, dup, pct)
colnames(priors)[1] <- paste("#", colnames(priors)[1], sep="")
write.table(priors, paste(out.prefix, "empirical_prior_estimates.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

# Plot mock H0/H1 cdf curves
pdf(paste(out.prefix, "mock_h0h1_curves.pdf", sep="."), height=4, width=10)
plot.example.priors(del, dup, colors)
dev.off()

