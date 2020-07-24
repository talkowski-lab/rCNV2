#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Optimize number of genes per CNV for gene scoring


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############

# Process a single input tsv of genes with counts and compute CDFs
load.data.single <- function(path, max.eval){
  x <- read.table(path, sep="\t", header=T, comment.char="")
  case.idx <- which(x$pheno=="case")
  control.idx <- which(x$pheno=="control")
  cdfs <- as.data.frame(do.call("rbind", lapply(1:max.eval, function(k){
    kidx <- which(x$n_genes <= k)
    case.signal <- sum(x$n_pos[intersect(case.idx, kidx)])
    case.noise <- sum(x$n_neg[intersect(case.idx, kidx)])
    case.ratio <- case.signal/case.noise
    control.signal <- sum(x$n_pos[intersect(control.idx, kidx)])
    control.noise <- sum(x$n_neg[intersect(control.idx, kidx)])
    control.ratio <- control.signal/control.noise
    cc.ratio <- case.ratio / control.ratio
    cc.ratio.diff <- case.ratio - control.ratio
    c(k, 
      case.signal, case.noise, case.ratio,
      control.signal, control.noise, control.ratio,
      cc.ratio, cc.ratio.diff)
  })))
  cdfs <- as.data.frame(apply(cdfs, 2, as.numeric))
  colnames(cdfs) <- c("n_genes", "case_pos", "case_neg", "case_ratio",
                      "control_pos", "control_neg", "control_ratio",
                      "ratio_of_ratios", "cc_ratio_diff")
  return(cdfs)
}

# Load & process input tsvs for all cohorts
load.data <- function(infile, max.eval){
  x <- read.table(infile, sep="\t", header=F)
  dlist <- lapply(1:nrow(x), function(i){
    load.data.single(x[i, 2], max.eval)
  })
  names(dlist) <- as.character(x[, 1])
  return(dlist)
}

# Joint optimization of ratio differences 
joint.opt <- function(dlist, max.eval){
  # Compute joint ratio diff weighted by proportion of positive hits in cases
  wdiffs <- sapply(1:max.eval, function(k){
    case.hits <- sapply(dlist, function(df){df$case_pos[which(df$n_genes==k)]})
    weights <- case.hits / sum(case.hits)
    bests <- sapply(dlist, function(df){max(df$cc_ratio_diff, na.rm=T)})
    diffs <- sapply(dlist, function(df){df$cc_ratio_diff[which(df$n_genes==k)]})
    sum(weights * (diffs-bests))
  })
  dlist[[1]]$n_genes[which(wdiffs==max(wdiffs, na.rm=T))]
}

# Plot case/control ratio diffs for a single cohort
plot.diffs <- function(dlist, cohort, ymax=NULL, color="red"){
  x <- as.data.frame(dlist[[cohort]])
  xmax <- max(dlist[[1]]$n_genes)
  if(is.null(ymax)){
    ymax <- max(c(x$case_ratio, x$control_ratio), na.rm=T)
  }
  best.idx <- which(x$cc_ratio_diff==max(x$cc_ratio_diff, na.rm=T))
  par(mar=c(4, 4, 3, 0.5), bty="n")
  plot(x$n_genes, x$control_ratio, xlim=c(0, xmax), ylim=c(0, ymax), 
       xaxs="i", yaxs="i", xpd=T,
       lwd=3, type="l", col="gray50",
       xlab="Max. Genes per CNV", ylab="Signal-to-Noise Ratio", main=cohort,
       panel.first=c(polygon(x=c(x$n_genes, rev(x$n_genes)),
                             y=c(x$control_ratio, rev(x$case_ratio)),
                             col=adjustcolor(color, alpha=0.2), border=NA),
                     segments(x0=x$n_genes[best.idx], x1=x$n_genes[best.idx],
                              y0=x$case_ratio[best.idx], y1=x$control_ratio[best.idx])))
  points(x$n_genes, x$case_ratio, type="l", lwd=3, col=color)
}

# Wrapper for plotting all cohorts
plot.all <- function(del, dup){
  # Get plot data
  ymax.del <- max(sapply(del, function(df){max(c(df$case_ratio, df$control_ratio), na.rm=T)}), na.rm=T)
  ymax.dup <- max(sapply(dup, function(df){max(c(df$case_ratio, df$control_ratio), na.rm=T)}), na.rm=T)
  
  # Prep layout
  par(mfrow=c(2, length(del)))
  
  # Plot signal:noise ratio for all cohorts
  sapply(names(del), function(cohort){
    plot.diffs(del, cohort, ymax.del)
  })
  legend("topright", lwd=3, col=c("red", "grey50"), legend=c("Cases", "Controls"), cex=0.9, bty="n")
  sapply(names(dup), function(cohort){
    plot.diffs(dup, cohort, ymax.dup, "blue")
  })
  legend("topright", lwd=3, col=c("blue", "grey50"), legend=c("Cases", "Controls"), cex=0.9, bty="n")
}


################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--max-eval"), default=80,
              help="Maximum number of genes per CNVs to evaluate")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog del.in dup.in out.pdf",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
del.in <- args$args[1]
dup.in <- args$args[2]
out.pdf <- args$args[3]
max.eval <- opts$`max-eval`

# # Dev parameters
# del.in <- "~/scratch/optimize_genes_per_cnv.DEL.input.tsv"
# dup.in <- "~/scratch/optimize_genes_per_cnv.DUP.input.tsv"
# out.pdf <- "~/scratch/test.pdf"
# max.eval <- 80
# setwd("~/scratch/")

# Load data and compute summaries
del <- load.data(del.in, max.eval)
dup <- load.data(dup.in, max.eval)

# Run optimization
del.opt <- joint.opt(del, max.eval)
dup.opt <- joint.opt(dup, max.eval)
print(paste("Optimal cutoffs:", del.opt, "(DEL);", dup.opt, "(DUP)"))

# Plot optimization
pdf(out.pdf, height=4, width=8)
plot.all(del, dup)
dev.off()

