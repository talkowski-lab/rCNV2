#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Calculate empirical FDR from permuted sliding window data


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############

# Load HPO sample size information
load.hpos <- function(hpos.in){
  hpos <- read.table(hpos.in, sep="\t", header=T, comment.char="")
  data.frame("hpo"=gsub(":", "", hpos[, 1], fixed=T),
             "n"=as.numeric(hpos[, 3]))
}

# Calculate FDR across a range of cutoffs for a vector of p-values
calc.fdr <- function(pv, cutoffs){
  pv <- pv[which(!is.na(pv))]
  n.pv <- length(pv)
  sapply(cutoffs, function(x){
    length(which(pv>=x))/n.pv
  })
}

# Read p-values and calculate FDR CDF
load.fdrs <- function(pvals.in, cutoffs, cnvtype=NULL){
  pvals <- read.table(pvals.in, header=T, sep="\t", comment.char="")
  if(!is.null(cnvtype)){
    pvals <- pvals[, grep(cnvtype, colnames(pvals), fixed=T)]
  }
  fdr.mat <- apply(pvals, 2, calc.fdr, cutoffs=cutoffs)
  rownames(fdr.mat) <- cutoffs
  return(fdr.mat)
}

# Calculate minimum P-value cutoff corresponding to target FDR for a vector of p-values
hit.fdr.target <- function(fdrv, cutoffs, target){
  cutoffs[min(c(head(which(fdrv<=target), 1), length(fdrv)))]
}

# Calculate P-value cutoffs per permutation corresponding to target FDRs for a matrix of p-values
get.fdr.cutoffs <- function(fdr.mat, fdr.target){
  fdr.cutoffs <- as.data.frame(apply(fdr.mat, 2, hit.fdr.target, cutoffs=cutoffs, target=fdr.target))
  colnames(fdr.cutoffs) <- "fdr.cutoff"
  return(fdr.cutoffs)
}

# Create cutoff table of all permutations
make.cutoff.mat <- function(fdr.mat, hpos, fdr.target){
  perm.hpos <- unlist(lapply(strsplit(colnames(fdr.mat), split=".", fixed=T), function(vals){vals[1]}))
  perm.n <- as.numeric(sapply(perm.hpos, function(hpo){hpos$n[which(hpos$hpo==hpo)]}))
  case.frac <- perm.n / (perm.n + hpos$n[which(hpos$hpo=="HEALTHY_CONTROL")])
  perm.idx <- unlist(lapply(strsplit(colnames(fdr.mat), split=".", fixed=T), function(vals){vals[3]}))
  fdr.res <- get.fdr.cutoffs(fdr.mat, fdr.target)
  cutoff.mat <- data.frame("hpo"=perm.hpos,
                           "perm"=perm.idx,
                           "n.cases"=perm.n,
                           "case.frac"=case.frac,
                           fdr.res)
  rownames(cutoff.mat) <- NULL
  return(cutoff.mat)
}

# Calculate a statistic for permuted cutoffs for each phenotype
get.cutoff.stat <- function(cutoff.mat, stat){
  stat.df <- as.data.frame(do.call("rbind", lapply(unique(cutoff.mat$hpo), function(hpo){
    if(stat=="max"){
      x <- max(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], na.rm=T)
    }else if(stat=="mean"){
      x <- mean(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], na.rm=T)
    }else if(stat=="median"){
      x <- median(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], na.rm=T)
    }else if(stat=="Q3"){
      x <- quantile(cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)], probs=0.75, na.rm=T)
    }
    c(hpo,
      head(cutoff.mat[which(cutoff.mat$hpo==hpo), 3:4], 1),
      unlist(x))
  })))
  colnames(stat.df) <- c("hpo", "n.cases", "case.frac", "fdr.cutoff")
  stat.df[, -1] <- apply(stat.df[, -1], 2, as.numeric)
  return(stat.df)
}

# Fit exponential decay to a list of points
fit.exp.decay <- function(x, y){
  # Following example on https://rpubs.com/mengxu/exponential-model
  
  train.df <- data.frame("x"=as.numeric(x), 
                         "y"=as.numeric(y))
  
  # Estimate initial parameters
  theta.0 <- 0.5 * min(train.df$y)
  model.0 <- lm(log(y - theta.0) ~ x, data=train.df)
  alpha.0 <- exp(coef(model.0)[1])
  beta.0 <- coef(model.0)[2]
  start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
  
  # Re-fit the model with estimated starting parameters
  return(nls(y ~ alpha * exp(beta * x) + theta, start = start, data=train.df,
             control=nls.control(maxiter=1000, warnOnly=T)))
}

# Calculate adjusted p-value cutoffs per HPO term
get.adjusted.cutoffs <- function(cutoff.stat.df, pred.n=NULL){
  fit <- fit.exp.decay(cutoff.stat.df$n.cases, cutoff.stat.df$fdr.cutoff)
  if(is.null(pred.n)){
    pred.n <- cutoff.stat.df$n.cases
    pred.labs <- unlist(cutoff.stat.df$hpo)
    col1.name <- "#hpo"
  }else{
    pred.labs <- pred.n
    col1.name <- "#n_cases"
  }
  fit.df <- data.frame("x"=pred.n)
  fit.df$y <- predict(fit, newdata=fit.df)
  out.df <- as.data.frame(cbind(pred.labs, unlist(10^-fit.df$y)))
  out.df <- out.df[!duplicated(out.df), ]
  colnames(out.df) <- c(col1.name, "min_p")
  return(out.df)
}

# Plot FDR, annotated with means and fitted curve
plot.fdrs <- function(cutoff.mat, cutoff.stat.df, stat, fdr.target, title=NULL, floor=T){
  
  phred.target <- -log10(fdr.target)
  if(is.null(title)){
    title <- stat
  }
  
  # Prep plot area & add points
  par(mar=c(3, 3, 2, 0.5))
  plot(x=c(0, max(cutoff.mat$n.cases)),
       y=c(0, max(cutoff.mat$fdr.cutoff)),
           type="n", xaxt="n", yaxt="n", xlab="", ylab="")
  points(x=cutoff.mat$n.cases, y=cutoff.mat$fdr.cutoff, col="gray75")
  abline(h=phred.target, lty=2)
  points(x=cutoff.stat.df$n.cases, y=cutoff.stat.df$fdr.cutoff, pch=15)
  
  # Add trendline
  fit <- fit.exp.decay(cutoff.stat.df$n.cases, cutoff.stat.df$fdr.cutoff)
  fit.df <- data.frame("x"=round(quantile(0:par("usr")[2], probs=seq(0, 1, 0.005))))
  fit.df$y <- predict(fit, newdata=fit.df)
  if(floor==T){
    fit.df$y[which(fit.df$y < phred.target)] <- phred.target
  }
  lines(fit.df$x, fit.df$y, col="red", lwd=3)
  
  # Add axes
  axis(1, labels=NA)
  axis(1, at=axTicks(1), tick=F, line=-0.5, cex.axis=0.85, labels=axTicks(1)/1000)
  mtext(1, line=1.5, text="Case samples (thousands)")
  axis(2, labels=NA)
  axis(2, at=axTicks(2), tick=F, line=-0.5, las=2, cex.axis=0.85)
  mtext(2, line=1.5, text=bquote(-log[10](italic(P))))
  mtext(3, line=0.1, text=title, font=2)
  
  # Add legend
  legend("topright", pch=c(1, 15, NA, NA), lwd=c(1, 1, 3, 1), lty=c(NA, NA, 1, 2), 
         col=c("gray75", "black", "red", "black"), legend=c("Permutation", stat, "Fit", "FDR Target"), 
         cex=0.9)
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)
require(MASS, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--cnv"), type="character", default=NULL,
              help="restrict analysis to DEL or DUP [default: consider both DEL + DUP]"),
  make_option(c("--max-cutoff"), type="numeric", default=20,
              help="max P-value cutoff to evaluate (Phred-scaled) [default %default]"),
  make_option(c("--cutoff-step"), type="numeric", default=0.05,
              help="P-value increments to evaluate (Phred-scaled) [default %default]"),
  make_option(c("--fdr-target"), type="numeric", default=0.01,
              help="FDR target [default %default]"),
  make_option(c("--plot"), type="character", default=NULL,
              help="path to .png of FDR fit vs. permutations [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog pval_matrix hpo_table out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

if(length(args$args) != 3){
  stop("Incorrect number of positional arguments specified.")
}

pvals.in <- args$args[1]
hpos.in <- args$args[2]
out.prefix <- args$args[3]
cnvtype <- opts$cnv
max.cutoff <- opts$`max-cutoff`
cutoff.step <- opts$`cutoff-step`
fdr.target <- opts$`fdr-target`
plot.out <- opts$`plot`

# # DEV PARAMETERS
# pvals.in <- "~/scratch/rCNV.permuted_pval_matrix.txt.gz"
# hpos.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/meta_cutoffs.test"
# cnvtype <- "DEL"
# max.cutoff <- 20
# cutoff.step <- 0.05
# fdr.target <- 0.00000385862
# plot.out <- "~/scratch/meta_cutoffs.test.png"

# Load sample size info
hpos <- load.hpos(hpos.in)

# Set FDR cutoffs
cutoffs <- seq(0, max.cutoff, cutoff.step)

# Calculate FDR CDFs for each permutation
fdr.mat <- load.fdrs(pvals.in, cutoffs, cnvtype)

# Calculate p-value cutoffs for each permutation for each target
cutoff.mat <- make.cutoff.mat(fdr.mat, hpos, fdr.target)

# Calculate various stats of cutoffs for each phenotype
mean.cutoffs <- get.cutoff.stat(cutoff.mat, "mean")
median.cutoffs <- get.cutoff.stat(cutoff.mat, "median")
q3.cutoffs <- get.cutoff.stat(cutoff.mat, "Q3")
max.cutoffs <- get.cutoff.stat(cutoff.mat, "max")

# Format & write output file matched to HPOs
df.out.mean <- get.adjusted.cutoffs(mean.cutoffs)
outfile.hpos <- paste(out.prefix, ".hpo_cutoffs.tsv", sep="")
write.table(df.out.mean, outfile.hpos, sep="\t", quote=F,
            col.names=T, row.names=F)

# Format & write output file against arbitrary sample size steps
df.out.ladder <- get.adjusted.cutoffs(mean.cutoffs, 
                                      pred.n=as.numeric(sapply(3:5, function(x){seq(1, 9.9, 0.1) * 10 ^ x})))
outfile.ladder <- paste(out.prefix, ".ncase_cutoff_ladder.tsv", sep="")
write.table(df.out.ladder, outfile.ladder, sep="\t", quote=F,
            col.names=T, row.names=F)

# Plot FDR data, if optioned
if(!is.null(plot.out)){
  # png(plot.out, res=300, height=5*300, width=6*300)
  # par(mfrow=c(2, 2))
  # plot.fdrs(cutoff.mat, mean.cutoffs, "Mean", fdr.target, floor=F)
  # plot.fdrs(cutoff.mat, median.cutoffs, "Median", fdr.target, floor=F)
  # plot.fdrs(cutoff.mat, q3.cutoffs, "Third quartile", fdr.target, floor=F)
  # plot.fdrs(cutoff.mat, max.cutoffs, "Max", fdr.target, floor=F)
  # dev.off()
  png(plot.out, res=300, height=4*300, width=5*300)
  plot.fdrs(cutoff.mat, mean.cutoffs, "Mean", fdr.target, 
            title=paste(cnvtype, "Permutation Results for FDR =",
                        format(fdr.target, scientific=T, digits=3)),
            floor=F)
  dev.off()
}

