#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
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
    }else if(stat=="q3"){
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

# Calculate weighted mean
calc.wmean <- function(vals, weights){
  sum(weights * vals) / sum(weights)
}

# Calculate adjusted p-value cutoffs per HPO term
get.adjusted.cutoffs <- function(cutoff.mat, stat, pred.n=NULL, exclude.outliers=T, linear.fit=F){
  # Determine outlier permutations and phenotypes, if optioned
  if(exclude.outliers==T){
    # Find individual outlier permutations (matched on hpo)
    outlier.perms <- unlist(lapply(unique(cutoff.mat$hpo), function(hpo){
      vals <- cutoff.mat$fdr.cutoff[which(cutoff.mat$hpo==hpo)]
      vals.quants <- quantile(vals, probs=c(0.25, 0.75), na.rm=T)
      vals.iqr <- vals.quants[2] - vals.quants[1]
      outlier.cutoffs <- c(vals.quants[1] - (1.5 * vals.iqr),
                           vals.quants[2] + (1.5 * vals.iqr))
      which(cutoff.mat$hpo==hpo & (cutoff.mat$fdr.cutoff < outlier.cutoffs[1] | cutoff.mat$fdr.cutoff > outlier.cutoffs[2]))
    }))
    # Find outlier phenotypes
    cutoff.stat.df.noOutliers <- get.cutoff.stat(cutoff.mat[-outlier.perms, ], stat)
    stat.quants <- quantile(cutoff.stat.df.noOutliers$fdr.cutoff, probs=c(0.25, 0.75), na.rm=T)
    stat.iqr <- stat.quants[2] - stat.quants[1]
    outlier.cutoffs <- c(stat.quants[1] - (1.5 * stat.iqr),
                         stat.quants[2] + (1.5 * stat.iqr))
    outlier.hpos <- unlist(cutoff.stat.df.noOutliers$hpo[which(cutoff.stat.df.noOutliers$fdr.cutoff < outlier.cutoffs[1]
                                                        | cutoff.stat.df.noOutliers$fdr.cutoff > outlier.cutoffs[2])])
    # Combine indexes to exclude
    outlier.exclude.idxs <- unique(c(outlier.perms, which(cutoff.mat$hpo %in% outlier.hpos)))
    cutoff.stat.df <- get.cutoff.stat(cutoff.mat[-outlier.exclude.idxs, ], stat)
  }else{
    cutoff.stat.df <- get.cutoff.stat(cutoff.mat, stat)
  }

  if(is.null(pred.n)){
    pred.labs <- unlist(unique(cutoff.mat$hpo))
    pred.n <- sapply(pred.labs, function(hpo){head(cutoff.mat$n.cases[which(cutoff.mat$hpo==hpo)], 1)})
    col1.name <- "#hpo"
  }else{
    pred.labs <- pred.n
    col1.name <- "#n_cases"
  }

  if(linear.fit==T){
    weights <- sqrt(cutoff.stat.df$n.cases)
    wmean <- calc.wmean(cutoff.stat.df$fdr.cutoff, weights)
    fit.df <- data.frame("x"=pred.n, "y"=wmean)
    fit <- NULL
  }else{
    fit <- fit.exp.decay(cutoff.stat.df$n.cases, cutoff.stat.df$fdr.cutoff)
    fit.df <- data.frame("x"=pred.n)
    fit.df$y <- predict(fit, newdata=fit.df)
  }
  out.df <- as.data.frame(cbind(pred.labs, unlist(10^-fit.df$y)))
  out.df <- out.df[!duplicated(out.df), ]
  colnames(out.df) <- c(col1.name, "min_p")
  return(list("df"=out.df, "fit"=fit))
}

# Plot FDR, annotated with means and fitted curve
plot.fdrs <- function(cutoff.mat, cutoff.stat.df, stat, fdr.target,
                      model="exponential", title=NULL, linear.fit=F, plot.exp=T, floor=F){

  neglog10.target <- -log10(fdr.target)
  if(is.null(title)){
    title <- stat
  }

  # Prep plot area & add points for raw permutations
  par(mar=c(3, 3, 2, 0.5))
  plot(x=c(0, max(cutoff.mat$n.cases)),
       y=c(0, max(cutoff.mat$fdr.cutoff)),
           type="n", xaxt="n", yaxt="n", xlab="", ylab="")
  points(x=cutoff.mat$n.cases, y=cutoff.mat$fdr.cutoff, col="gray75")
  abline(h=neglog10.target, lty=2)

  # Add trendline
  if(linear.fit==T){
    abline(h=calc.wmean(cutoff.stat.df$fdr.cutoff, sqrt(cutoff.stat.df$n.cases)),
           lwd=3, col="red")
  }else if(model=="exponential"){
    fit <- get.adjusted.cutoffs(cutoff.mat, tolower(stat), exclude.outliers=T)$fit
    fit.df <- data.frame("x"=round(quantile(0:par("usr")[2], probs=seq(0, 1, 0.005))))
    fit.df$y <- predict(fit, newdata=fit.df)
    if(floor==T){
      fit.df$y[which(fit.df$y < neglog10.target)] <- neglog10.target
    }
    if(plot.exp==T){
      lines(fit.df$x, fit.df$y, col="red", lwd=3)
    }
  }

  # Add group means
  points(x=cutoff.stat.df$n.cases, y=cutoff.stat.df$fdr.cutoff, pch=15)

  # Add axes
  axis(1, labels=NA)
  axis(1, at=axTicks(1), tick=F, line=-0.5, cex.axis=0.85, labels=axTicks(1)/1000)
  mtext(1, line=1.5, text="Case samples (thousands)")
  axis(2, labels=NA)
  axis(2, at=axTicks(2), tick=F, line=-0.5, las=2, cex.axis=0.85)
  mtext(2, line=1.5, text=bquote(-log[10](italic(P))))
  mtext(3, line=0.1, text=title, font=2)

  # Add legend
  if(plot.exp==T){
    legend("topright", pch=c(1, 15, NA, NA), lwd=c(1, 1, 3, 1), lty=c(NA, NA, 1, 2),
           col=c("gray75", "black", "red", "black"), cex=0.9,
           legend=c("Permutation", stat, "Exponential Fit", "FDR Target"))
  }else if(linear.fit==T){
    legend("topright", pch=c(1, 15, NA, NA), lwd=c(1, 1, 3, 1), lty=c(NA, NA, 1, 2),
           col=c("gray75", "black", "red", "black"), cex=0.9,
           legend=c("Permutation", stat, "Weighted Mean", "FDR Target"))
  }else{
    legend("topright", pch=c(1, 15, NA), lwd=c(1, 1, 1), lty=c(NA, NA, 2),
           col=c("gray75", "black", "black"), cex=0.9,
           legend=c("Permutation", stat, "FDR Target"))
  }
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
  make_option(c("--linear-fit"), action="store_true", default=FALSE,
              help="Compute cross-phenotype weighted mean for universal cutoff instead of empirical mean per phenotype [default %default]"),
  make_option(c("--flat-ladder"), action="store_true", default=FALSE,
              help="Enforce true FDR target for ladder [default %default]"),
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
linear.fit <- opts$`linear-fit`
flat.ladder <- opts$`flat-ladder`
plot.out <- opts$`plot`

# # DEV PARAMETERS
# pvals.in <- "~/scratch/rCNV.DUP.permuted_pval_matrix.txt.gz"
# hpos.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/meta_cutoffs.test"
# cnvtype <- "DUP"
# max.cutoff <- 20
# cutoff.step <- 0.05
# fdr.target <- 3.097318E-6
# linear.fit <- T
# flat.ladder <- T
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
q3.cutoffs <- get.cutoff.stat(cutoff.mat, "q3")
max.cutoffs <- get.cutoff.stat(cutoff.mat, "max")

# Format & write output file matched to HPOs
# df.out.mean <- get.adjusted.cutoffs(cutoff.mat, "mean", exclude.outliers=F)$df
# df.out.mean.noOutliers <- get.adjusted.cutoffs(cutoff.mat, "mean", exclude.outliers=T)$df
if(linear.fit==T){
  # Compute flat weighted mean of all HPOs, if optioned
  df.out.mean <- get.adjusted.cutoffs(cutoff.mat, "mean", exclude.outliers=F,
                                      linear.fit=linear.fit)$df
}else{
  # Otherwise, compute empirical mean for each HPO independently
  df.out.mean <- data.frame(unlist(as.character(mean.cutoffs$hpo)),
                                      as.numeric(10^-mean.cutoffs$fdr.cutoff))
}
colnames(df.out.mean) <- c("#hpo", "min_p")
outfile.hpos <- paste(out.prefix, ".hpo_cutoffs.tsv", sep="")
write.table(df.out.mean, outfile.hpos, sep="\t", quote=F,
            col.names=T, row.names=F)

# Format & write output file against arbitrary sample size steps
if(flat.ladder == T){
  df.out.ladder <- data.frame("#n_cases"=as.numeric(sapply(3:5, function(x){seq(1, 9.9, 0.1) * 10 ^ x})),
                              "min_p"=fdr.target)
  colnames(df.out.ladder)[1] <- "#n_cases"
}else{
  df.out.ladder <- get.adjusted.cutoffs(cutoff.mat, "mean",
                                        pred.n=as.numeric(sapply(3:5, function(x){seq(1, 9.9, 0.1) * 10 ^ x})),
                                        exclude.outliers=T)$df
}
outfile.ladder <- paste(out.prefix, ".ncase_cutoff_ladder.tsv", sep="")
write.table(df.out.ladder, outfile.ladder, sep="\t", quote=F,
            col.names=T, row.names=F)

# Plot FDR data, if optioned
if(!is.null(plot.out)){
  plot.exp <- all(flat.ladder==F & linear.fit==F)
  png(plot.out, res=300, height=4*300, width=5*300)
  plot.fdrs(cutoff.mat, mean.cutoffs, "Mean", fdr.target,
            title=paste(cnvtype, "Permutation Results for FDR =",
                        format(fdr.target, scientific=T, digits=3)),
            linear.fit=linear.fit, plot.exp=plot.exp, floor=F)
  dev.off()
}

