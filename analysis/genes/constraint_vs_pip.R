#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Assess enrichments of constrained genes vs. fine-mapped PIPs for rapid QC


options(scipen=100000, stringsAsFactors=F)
require(rCNV2)


#################
### FUNCTIONS ###
#################
# Load PIPs from gene stats .tsv
load.stats <- function(stats.in){
  s <- read.table(stats.in, header=T, sep="\t", comment.char="")
  s <- s[which(!is.na(s$credible_set)), c("gene", "PIP_final")]
  r <- as.data.frame(do.call("rbind", lapply(unique(s$gene), function(g){
    c(g, max(s$PIP_final[which(s$gene == g)], na.rm=T))
  })))
  colnames(r) <- c("gene", "PIP")
  r$PIP <- as.numeric(r$PIP)
  r[rev(order(r$PIP)), ]
}

# Compute fraction of constrained genes as a function of PIP
get.frac.constr <- function(pips, constr, steps, equality="ge"){
  res <- as.data.frame(t(sapply(steps, function(step){
    if(equality == "ge"){
      idx <- which(pips$PIP >= step)
    }else if(equality == "lt"){
      idx <- which(pips$PIP < step)
    }
    hits <- length(intersect(idx, which(pips$gene %in% constr)))
    total <- length(idx)
    return(c(step, hits / total, hits, total))
  })))
  colnames(res) <- c("cutoff", "frac", "hits", "total")
  return(res)
}

# Plot constrained fractions vs. PIP
plot.constr.frac <- function(left.frac, right.frac, baseline){
  # Get plot parameters
  xlims <- c(1, 0)
  # ylims <- c(0, max(c(left.frac$frac, right.frac$frac), na.rm=T))
  ylims <- c(0, 0.5)

  # Prep plot area
  par(mar=c(3, 3, 0.25, 0.25), bty="n")
  plot(NA, xlim=xlims, ylim=ylims, xaxt="n", yaxt="n", xlab="", ylab="")

  # Add lines
  abline(h=baseline, lty=5, col=ns.color)
  points(x=left.frac$cutoff, y=left.frac$frac, lwd=2, col="#F95700FF", type="l")
  points(x=right.frac$cutoff, y=right.frac$frac, lwd=2, col="#00A4CCFF", type="l")

  # Add X-axis
  axis(1, labels=NA, tck=-0.025)
  axis(1, tick=F, line=-0.65, cex.axis=0.85)
  mtext(1, line=1.25, text="PIP")

  # Add Y-axis
  axis(2, at=c(0, 10e10), labels=NA, tck=0)
  axis(2, labels=NA, tck=-0.025)
  axis(2, tick=F, line=-0.65, cex.axis=0.85, las=2)
  mtext(2, line=1.75, text="Fraction Constrained")

  # Add legend
  legend("topright", col=c("#F95700FF", "#00A4CCFF", ns.color),
         cex=0.7, lty=c(1, 1, 5), lwd=c(2, 2, 1),
         legend=c("Genes above", "Genes below", "Baseline"))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog stats.tsv constrained.list out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Three positional arguments: stats.tsv, constrained.list, and out_prefix\n")
}

# Writes args & opts to vars
stats.in <- args$args[1]
constr.in <- args$args[2]
out.prefix <- args$args[3]
baseline <- 3052/18641 #Baseline rate of contsrained genes as a fraction of all genes

# # DEV PARAMTERS
# setwd("~/scratch")
# stats.in <- "rCNV.DEL.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv"
# constr.in <- "gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
# out.prefix <- "finemap_qc"

# Read data
pips <- load.stats(stats.in)
constr <- sort(unique(as.character(read.table(constr.in, header=F)[, 1])))

# Compute constrained proportions vs. PIP cutoffs
steps <- seq(1, 0, -0.01)
left.frac <- get.frac.constr(pips, constr, steps)
right.frac <- get.frac.constr(pips, constr, steps, equality="lt")

# Plot left & right fractions
pdf(paste(out.prefix, "pip_vs_constrained.pdf", sep="."),
    height=3.5, width=4)
plot.constr.frac(left.frac, right.frac, baseline)
dev.off()
