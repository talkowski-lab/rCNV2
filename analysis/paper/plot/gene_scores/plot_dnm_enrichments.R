#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Assess enrichment of de novo mutations from ASC & DDD exomes vs. gene-level intolerance scores


options(stringsAsFactors=F, scipen=1000)


# Constants
csqs <- c("lof", "mis", "syn")
dnm.cohorts <- c("ddd", "asc", "asc_unaffected")


######################
### DATA FUNCTIONS ###
######################
# Split genes into groups based on score bin
bin.genes <- function(scores, score, n.bins=10){
  breaks <- seq(0, 1, length.out=n.bins+1)
  lapply(2:length(breaks), function(i){
    scores$gene[which(scores[, score] >= breaks[i-1]
                      & scores[, score] <= breaks[i])]
  })
}

# Compute expected DNMs for a single consequence & cohort for all genes
calc.exp.dnms.single <- function(meta, cohort, csq, neutral.genes){
  keep.idxs <- which(meta$gene %in% neutral.genes)
  obs <- meta[, paste(cohort, "dn", csq, sep="_")]
  mu <- meta[, paste("gnomad_mu", csq, sep="_")]
  obs.sum <- sum(obs[keep.idxs], na.rm=T)
  mu.sum <- sum(mu[keep.idxs], na.rm=T)
  scalar <- obs.sum / mu.sum
  exp <- mu * scalar
  meta[, paste(cohort, "dn", csq, "exp", sep="_")] <- exp
  return(meta)
}

# Compute expected DNMs for all consequences & cohorts
calc.exp.dnms <- function(meta, cohorts, csqs, neutral.genes){
  for(cohort in cohorts){
    for(csq in csqs){
      meta <- calc.exp.dnms.single(meta, cohort, csq, neutral.genes)
    }
  }
  return(meta)
}

# Helper function for DNM O/e boot (bootstrapping confidence intervals)
summed.oe <- function(obs.exp.df, indexes=NULL){
  if(is.null(indexes)){
    indexes <- 1:nrow(obs.exp.df)
  }
  sum.obs <- sum(obs.exp.df[indexes, 1])
  sum.exp <- sum(obs.exp.df[indexes, 2])
  sum.obs / sum.exp
}

# Compute obs/exp DNMs for a single consequence & cohort for a set of genes
calc.dnm.oe.single <- function(genes, meta, cohort, csq){
  obs.exp.df <- meta[which(meta$gene %in% genes), 
                     grep(paste(cohort, "dn", csq, sep="_"), colnames(meta))]
  oe <- summed.oe(obs.exp.df)
  ci <- boot.ci(boot(data=obs.exp.df, statistic=summed.oe, R=1000), 
                conf=0.95, type="norm")$normal[, -1]
  as.numeric(c(oe, ci))
}

# Compute obs/exp DNMs for a single consequence & cohort for all sets of genes
calc.dnm.oe.multi <- function(gene.groups, meta, cohort, csq){
  oe.df <- as.data.frame(do.call("rbind", lapply(gene.groups, function(genes){
    calc.dnm.oe.single(genes, meta, cohort, csq)
  })))
  colnames(oe.df) <- c("oe", "lower", "upper")
  return(oe.df)
}

# Compute obs/exp DNMs for all consequences & genes for a single cohort
calc.dnm.oe.cohort <- function(gene.groups, meta, cohort, csqs){
  oe.data <- lapply(csqs, function(csq){
    calc.dnm.oe.multi(gene.groups, meta, cohort, csq)
  })
  names(oe.data) <- csqs
  return(oe.data)
}

# Helper function for BCA fraction disrupted boot (bootstrapping confidence intervals)
frac.disrupted <- function(bca.counts, indexes=NULL){
  if(is.null(indexes)){
    indexes <- 1:length(bca.counts)
  }
  length(which(bca.counts[indexes] > 0)) / length(indexes)
}

# Compute fraction of genes disrupted by at least one BCA for a set of genes
calc.frac.bcas <- function(gene.groups, meta){
  bca.counts <- meta$redin_any_bca
  baseline <- frac.disrupted(bca.counts)
  res <- data.frame(t(sapply(gene.groups, function(genes){
    idxs <- which(meta$gene %in% genes)
    frac <- frac.disrupted(bca.counts[idxs])
    ci <- boot.ci(boot(data=bca.counts[idxs], statistic=frac.disrupted, R=1000), 
                  conf=0.95, type="norm")$normal[, -1]
    as.numeric(c(frac, ci))
  })))
  colnames(res) <- c("frac", "lower", "upper")
  return(list("bins"=res, "baseline"=baseline))
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot rCNV score vs. DNM enrichment for a single cohort
plot.dnm.oe <- function(scores, score, meta, cohort, csqs, n.bins=10,
                        xlab=NULL, ylab=NULL, ymax=NULL, blue.bg=TRUE,
                        parmar=c(2.25, 2.6, 0.5, 0.5)){
  # Collect plot data
  gene.groups <- bin.genes(scores, score, n.bins)
  plot.dat <- calc.dnm.oe.cohort(gene.groups, meta, cohort, csqs)
  
  # Set parameters
  if(is.null(ymax)){
    ymax <- max(unlist(lapply(plot.dat, range, na.rm=T)), na.rm=T)
  }
  ylims <- c(0, ymax)
  if(is.null(xlab)){
    xlab <- paste("Genes Binned by", score)
  }
  if(is.null(ylab)){
    ylab <- "Fold-Enrichment"
  }
  if(blue.bg==TRUE){
    plot.bg <- bluewhite
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- "white"
  }else{
    plot.bg <- "white"
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- NA
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 1), ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  abline(h=axTicks(2), v=axTicks(1), col=grid.col)
  abline(h=1, col=blueblack, lty=2)
  
  # Add points
  x.breaks <- seq(0, 1, length.out=n.bins+1)
  x.at <- (x.breaks[-1] + x.breaks[-length(x.breaks)])/2
  x.mod <- 0.015
  sapply(length(plot.dat):1, function(i){
    segments(x0=x.at+(x.mod*(2-i)), x1=x.at+(x.mod*(2-i)), 
             y0=plot.dat[[i]][, 2], y1=plot.dat[[i]][, 3],
             col=snv.colors[i], lend="round", lwd=1.5)
    points(x=x.at+(x.mod*(2-i)), y=plot.dat[[i]][, 1], pch=19, col=snv.colors[i])
  })
  
  # Add axes
  x.ax.at <- axTicks(1)
  axis(1, at=c(-100, 100), col=blueblack, tck=0, labels=NA)
  axis(1, at=x.ax.at, labels=NA, tck=-0.025, col=blueblack)
  sapply(x.ax.at, function(k){
    axis(1, at=k, tick=F, line=-0.75)
  })
  mtext(1, text=xlab, line=1.1)
  y.ax.at <- axTicks(2)
  axis(2, at=c(-100, 100), col=blueblack, tck=0, labels=NA)
  axis(2, at=y.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(2, at=y.ax.at, tick=F, line=-0.65, las=2)
  mtext(2, text=ylab, line=1.6)
}

# Plot rCNV score vs. BCA disruptions
plot.bca.fracs <- function(scores, score, meta, n.bins=10,
                        xlab=NULL, ymax=NULL, blue.bg=TRUE,
                        parmar=c(2.25, 4, 0.5, 0.5)){
  # Collect plot data
  gene.groups <- bin.genes(scores, score, n.bins)
  all.plot.dat <- calc.frac.bcas(gene.groups, meta)
  plot.dat <- all.plot.dat$bins
  baseline <- all.plot.dat$baseline
  
  # Set parameters
  if(is.null(ymax)){
    ymax <- max(plot.dat, na.rm=T)
  }
  ylims <- c(0, ymax)
  if(is.null(xlab)){
    xlab <- paste("Genes Binned by", score)
  }
  if(score=="pHI"){
    color <- cnv.colors[1]
  }else if(score=="pTS"){
    color <- cnv.colors[2]
  }else{
    color <- blueblack
  }
  if(blue.bg==TRUE){
    plot.bg <- bluewhite
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- "white"
  }else{
    plot.bg <- "white"
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- NA
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 1), ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  abline(h=axTicks(2), v=axTicks(1), col=grid.col)
  abline(h=baseline, col=blueblack, lty=2)
  text(x=par("usr")[1]-(0.05*(par("usr")[2]-par("usr")[1])),
       y=baseline+(0.04*(par("usr")[4]-par("usr")[3])), 
       labels="Baseline", cex=0.85, font=3, pos=4, col=blueblack)
  
  # Add points
  x.breaks <- seq(0, 1, length.out=n.bins+1)
  x.at <- (x.breaks[-1] + x.breaks[-length(x.breaks)])/2
  segments(x0=x.at, x1=x.at, 
           y0=plot.dat[, 2], y1=plot.dat[, 3],
           col=color, lend="round", lwd=1.5)
  points(x=x.at, y=plot.dat[, 1], pch=19, col=color)
  
  # Add axes
  x.ax.at <- axTicks(1)
  axis(1, at=c(-100, 100), col=blueblack, tck=0, labels=NA)
  axis(1, at=x.ax.at, labels=NA, tck=-0.025, col=blueblack)
  sapply(x.ax.at, function(k){
    axis(1, at=k, tick=F, line=-0.75)
  })
  mtext(1, text=xlab, line=1.1)
  y.ax.at <- axTicks(2)
  axis(2, at=c(-100, 100), col=blueblack, tck=0, labels=NA)
  axis(2, at=y.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(2, at=y.ax.at, labels=100*y.ax.at, tick=F, line=-0.65, las=2)
  mtext(2, text="Genes Disrupted by", line=3)
  mtext(2, text=bquote(italic("De Novo") ~ "BCAs" ~ ("%")), line=1.8)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(boot, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog rCNV_scores.tsv gene.meta.bed neutral.genes out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, gene.meta.bed, neutral.genes, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
meta.in <- args$args[2]
neutral.genes.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# meta.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.bed.gz"
# neutral.genes.in <- "~/scratch/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list"
# out.prefix <- "~/scratch/test_gene_score_dnm_analyses"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load gene scores & divide genes by score
scores <- load.scores(scores.in)

# Load list of neutral genes (for calibrating mutation rates)
neutral.genes <- as.character(read.table(neutral.genes.in, header=F)[, 1])

# Load gene metadata & compute expected DNMs per cohort & consequence per gene
meta <- load.gene.metadata(meta.in)
meta <- meta[which(meta$gene %in% scores$gene), ]
meta <- calc.exp.dnms(meta, dnm.cohorts, csqs, neutral.genes)

# Plot DNM enrichment per score per exome cohort 
sapply(c("pHI", "pTS"), function(score){
  sapply(dnm.cohorts, function(cohort){
    pdf(paste(out.prefix, "dnm_enrichments", cohort, score, "pdf", sep="."),
        height=2, width=3)
    plot.dnm.oe(scores, score, meta, cohort, csqs, ymax=5, blue.bg=FALSE)
    dev.off()
  })
})

# Plot fraction of genes disrupted by de novo BCA
sapply(c("pHI", "pTS"), function(score){
  pdf(paste(out.prefix, "dnm_enrichments.BCA", score, "pdf", sep="."),
      height=2, width=2.6)
  plot.bca.fracs(scores, score, meta, blue.bg=FALSE)
  dev.off()
})

