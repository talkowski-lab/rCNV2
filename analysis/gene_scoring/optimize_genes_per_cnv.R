#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
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
  x$true_only <- x$n_pos>0 & x$n_neg<1
  x$ambiguous <- x$n_pos>0 & x$n_neg>0
  # Approximate proportion of true informative CNVs added
  # for each step in max.eval
  if(any(x$true_only)){
    true.added <- dnorm(1:max.eval,
                        mean=mean(x$n_genes[x$true_only]),
                        sd=sd(x$n_genes[x$true_only]))
  }else{
    true.added <- rep(1/max.eval, len=max.eval)
  }
  # Approximate proportion of uninformative CNVs (noise) added
  # for each step in max.eval
  if(any(x$ambiguous)){
    ambig.added <- dnorm(1:max.eval,
                         mean=mean(x$n_genes[x$ambiguous]),
                         sd=sd(x$n_genes[x$ambiguous]))
  }else{
    ambig.added <- rep(1/max.eval, len=max.eval)
  }
  cdfs <- data.frame("true.added" = true.added,
                     "ambig.added" = ambig.added,
                     "d" = true.added - ambig.added)
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

# Optimize cutoffs for a list of cohorts
opt.genes.per.cnv <- function(dlist, max.eval){
  # Compute optimal cutoff per cohort
  opt <- sapply(dlist, function(df){
    # Optimal cutoff is the first value where the derivative <= 0
    # Subtract one to go to step immediately before this point
    min(which(df$d <= 0)) - 1
  })
  names(opt) <- names(dlist)

  # Average derivatives across cohorts
  d.avg <- apply(do.call("cbind", lapply(dlist, function(df){df$d})), 1, mean)
  joint.opt <- min(which(d.avg <= 0)) - 1
  return(list("per.cohort"=opt, "joint"=joint.opt))
}

# Plot optimization derivative for a single cohort
plot.diffs <- function(df, cohort, opt, ylims=NULL, color="black"){
  d <- df$d
  xmax <- nrow(df)
  if(is.null(ylims)){
    ylims <- max(d, na.rm=T)
  }
  par(mar=c(4, 4, 3, 0.5), bty="n")
  plot(1:xmax, d, xlim=c(0, xmax), ylim=ylims,
       xaxs="i", yaxs="i", xpd=T,
       lwd=3, type="l", col=color,
       xlab="Max. Genes / CNV", ylab="Optimization Derivative", main=cohort,
       panel.first=c(abline(h=0, col="gray70")))
  abline(v=opt)
}

# Wrapper for plotting all cohorts
plot.all <- function(del, dup, del.opt, dup.opt){
  # Get plot data
  ylims.del <- range(sapply(del, function(df){range(df$d, na.rm=T)}), na.rm=T)
  ylims.dup <- range(sapply(dup, function(df){range(df$d, na.rm=T)}), na.rm=T)

  # Prep layout
  par(mfrow=c(2, length(del)))

  # Plot signal:noise ratio for all cohorts
  sapply(names(del), function(cohort){
    plot.diffs(del[[cohort]], cohort, del.opt[cohort],
               ylims.del, cnv.colors[1])
  })
  sapply(names(dup), function(cohort){
    plot.diffs(dup[[cohort]], cohort, dup.opt[cohort],
               ylims.dup, cnv.colors[2])
  })
}


################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--max-eval"), default=100,
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
# max.eval <- 100
# setwd("~/scratch/")

# Load data and compute summaries
del <- load.data(del.in, max.eval)
dup <- load.data(dup.in, max.eval)

# Run optimization per CNV type
del.opt <- opt.genes.per.cnv(del, max.eval)
dup.opt <- opt.genes.per.cnv(dup, max.eval)
print(paste("Optimal cutoffs per CNV type:",
            floor(mean(del.opt$joint)), "(DEL);",
            floor(mean(dup.opt$joint)), "(DUP)"))

# Run joint optimization
joint.opt <- opt.genes.per.cnv(c(del, dup), max.eval)
print(paste("Optimal cutoffs when averaged across both CNV type:",
            floor(mean(joint.opt$joint))))

# Plot optimization
pdf(out.pdf, height=4, width=10)
plot.all(del, dup, del.opt$per.cohort, dup.opt$per.cohort)
dev.off()

