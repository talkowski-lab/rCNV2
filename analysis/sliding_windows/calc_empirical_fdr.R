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
get.fdr.cutoffs <- function(fdr.mat, fdr.targets){
  fdr.cutoffs <- sapply(fdr.targets, function(t){
    apply(fdr.mat, 2, hit.fdr.target, cutoffs=cutoffs, target=t)
  })
  colnames(fdr.cutoffs) <- paste("FDR", fdr.targets, sep="_")
  return(fdr.cutoffs)
}

# Create cutoff table of all permutations
make.cutoff.mat <- function(fdr.mat, hpos, fdr.targets){
  perm.hpos <- unlist(lapply(strsplit(colnames(fdr.mat), split=".", fixed=T), function(vals){vals[1]}))
  perm.n <- as.numeric(sapply(perm.hpos, function(hpo){hpos$n[which(hpos$hpo==hpo)]}))
  perm.idx <- unlist(lapply(strsplit(colnames(fdr.mat), split=".", fixed=T), function(vals){vals[3]}))
  fdr.res <- get.fdr.cutoffs(fdr.mat, fdr.targets)
  cutoff.mat <- data.frame("hpo"=perm.hpos,
                           "perm"=perm.idx,
                           "n.samples"=perm.n,
                           fdr.res)
  rownames(cutoff.mat) <- NULL
  return(cutoff.mat)
}


# Plot FDR traces, annotated with cutoffs
plot.fdrs <- function(fdr.mat, cutoffs, cutoff.table){
  # Get plot values
  xmax <- 1.5*max(cutoff.table[, 2])
  ymax <- 1.5*max(cutoff.table[, 1])
  
  # Prep plot area
  par(mar=c(3, 3, 0.5, 0.5))
  plot(x=c(0, xmax), y=c(0, ymax), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  
  # Add axes
  axis(1, labels=NA)
  axis(1, tick=F, line=-0.5, cex.axis=0.85)
  mtext(1, line=1.5, text=bquote(-log[10](italic(p))))
  axis(2, labels=NA)
  axis(2, at=axTicks(2), tick=F, line=-0.5, las=2, cex.axis=0.85,
       labels=paste(round(100*axTicks(2)), "%", sep=""))
  mtext(2, line=1.8, text="FDR")
  
  # Add FDR traces
  sapply(1:ncol(fdr.mat), function(k){
    points(x=cutoffs, y=fdr.mat[, k], type="l",
           col=adjustcolor("black", alpha=0.2))
  })
  
  # Add optimized FDR cutoffs
  textLab.ybuffer <- 0.05*(par("usr")[4] - par("usr")[3])
  sapply(1:nrow(cutoff.table), function(i){
    abline(h=cutoff.table[i, 1], col="red", lty=2)
    points(x=cutoff.table[i, 2], y=cutoff.table[i, 1],
           pch=19, col="red")
    text(x=cutoff.table[i, 2], y=cutoff.table[i, 1] + textLab.ybuffer,
         labels=paste("FDR = ", round(100*cutoff.table[i, 1], 2), "%\n",
                      "-log10(P) = ", round(cutoff.table[i, 2], 2), sep=""),
         col="red", pos=4)
  })
  
  
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
  make_option(c("--fdr-targets"), type="character", default="0.1,0.05,0.01",
              help="FDR targets to calculate (comma-delimited list) [default %default]"),
  make_option(c("--plot"), type="character", default=NULL,
              help="path to .png of FDR calculations [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog pval_matrix hpo_table outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

if(length(args$args) != 3){
  stop("Incorrect number of positional arguments specified.")
}

pvals.in <- args$args[1]
hpos.in <- args$args[2]
outfile <- args$args[3]
cnvtype <- opts$cnv
max.cutoff <- opts$`max-cutoff`
cutoff.step <- opts$`cutoff-step`
fdr.targets <- as.numeric(unlist(strsplit(opts$`fdr-targets`, split=",")))
plot.out <- opts$`plot`

# # DEV PARAMETERS
# pvals.in <- "~/scratch/rCNV.permuted_pval_matrix.txt.gz"
# hpos.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# outfile <- "~/scratch/meta_cutoffs.test.txt"
# cnvtype <- "DEL"
# max.cutoff <- 20
# cutoff.step <- 0.05
# fdr.targets <- c(0.1, 0.05, 0.01)
# plot.out <- "~/scratch/meta_cutoffs.test.png"

# Load sample size info
hpos <- load.hpos(hpos.in)

# Set FDR cutoffs
cutoffs <- seq(0, max.cutoff, cutoff.step)

# Calculate FDR CDFs for each permutation
fdr.mat <- load.fdrs(pvals.in, cutoffs, cnvtype)

# Calculate p-value cutoffs for each permutation for each target
cutoff.mat <- make.cutoff.mat(fdr.mat, hpos, fdr.targets)

# Determine 99th percentile of FDRs per P-value cutoff
fdr.cutoffs <- get.fdr.cutoffs(fdr.mat, fdr.targets)

# Format & write output file
df.out <- data.frame("#fdr"=fdr.targets,
                     "phred_p_cutoff"=fdr.cutoffs,
                     "p_cutoff"=10^-fdr.cutoffs)
write.table(df.out, outfile, sep="\t", quote=F, 
            col.names=T, row.names=F)

# Plot FDR traces, if optioned
if(!is.null(plot.out)){
  png(plot.out, res=300, height=4*300, width=6*300)
  plot.fdrs(fdr.mat, cutoffs, df.out)
  dev.off()
}




