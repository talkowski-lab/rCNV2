#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot FDR permutation results from sliding window meta-analysis for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000, check.names=F)


######################
### DATA FUNCTIONS ###
######################
# Load HPO sample size information
load.hpos <- function(hpos.in){
  hpos <- read.table(hpos.in, sep="\t", header=T, comment.char="")
  data.frame("hpo"=gsub(":", "", hpos[, 1], fixed=T),
             "n"=as.numeric(hpos[, 3]))
}

# Load & phred-scale precomputed FDR matrix
load.fdrs <- function(fdrs.in){
  fdrs <- read.table(fdrs.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(fdrs) <- sapply(colnames(fdrs), function(x){unlist(strsplit(x, split=".", fixed=T))[1]})
  fdrs <- as.data.frame(apply(fdrs, 2, function(vals){-log10(as.numeric(vals))}))
  return(fdrs)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot permutation results against targeted FDR
plot.fdrs <- function(fdrs, hpos, cnv, fdr.target,
                      xlims=NULL, ylims=NULL, blue.bg=TRUE,
                      parmar=c(2.2, 2.5, 0.25, 0.25)){
  # Get plot data
  if(is.null(xlims)){
    xlims <- range(log10(hpos$n), na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- c(0, max(as.matrix(fdrs), na.rm=T))
  }
  medians <- apply(fdrs, 2, median, na.rm=T)
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
  plot(NA, xlim=xlims, ylim=ylims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty=plot.bty, border=plot.border, col=plot.bg)
  abline(h=axTicks(2), v=log10(logscale.major), col=grid.col)
  
  # Add points
  sapply(hpos[, 1], function(hpo){
    n <- log10(hpos$n[which(hpos$hpo==hpo)])
    hpo.idx <- which(colnames(fdrs)==hpo)
    points(x=rep(n, nrow(fdrs)), y=fdrs[, hpo.idx],
           pch=19, cex=0.2, col=control.cnv.colors[cnv])
  })
  
  # Add boxes for medians
  sapply(hpos[, 1], function(hpo){
    n <- log10(hpos$n[which(hpos$hpo==hpo)])
    hpo.idx <- which(colnames(fdrs)==hpo)
    points(x=n, y=medians[hpo.idx], pch=22, 
           bg=cnv.colors[cnv], col=cnv.blacks[cnv])
  })
  
  # Add line for Bonferroni
  abline(h=fdr.target, col=blueblack, lty=2, lwd=1.5)
  
  # Add X-axis
  axis(1, at=log10(logscale.minor), tck=-0.015, col=blueblack, labels=NA)
  axis(1, at=log10(logscale.major), tck=-0.03, col=blueblack, labels=NA)
  axis(1, at=log10(c(1000, 10000, 100000, 1000000)), tick=F, line=-0.75,
       labels=c("1", "10", "100", "1,000"))
  mtext(1, text="Cases (Thousands)", line=1.1)
  
  # Add Y-axis
  axis(2, at=c(-10e10, 10e10), labels=NA, col=blueblack, tck=0)
  axis(2, tck=-0.03, col=blueblack, labels=NA)
  axis(2, tick=F, las=2, line=-0.65)
  mtext(2, text=bquote("Permuted" ~ -log[10](italic(P))), line=1.25)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
  make_option(c("--fdr-target"), type="numeric", default=10e-8,
              help="FDR target [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog del.tsv dup.tsv hpo_table.tsv out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

if(length(args$args) != 4){
  stop("Incorrect number of positional arguments specified.")
}

del.fdrs.in <- args$args[1]
dup.fdrs.in <- args$args[2]
hpos.in <- args$args[3]
out.prefix <- args$args[4]
fdr.target <- -log10(opts$`fdr-target`)
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# del.fdrs.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.sliding_window.meta_analysis.stats.permuted_fdrs.tsv.gz"
# dup.fdrs.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DUP.sliding_window.meta_analysis.stats.permuted_fdrs.tsv.gz"
# hpos.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/sliding_windows_meta_fdr_perm_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# fdr.target <- -log10(0.000003715428)
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load sample size info
hpos <- load.hpos(hpos.in)
hpos <- hpos[which(hpos$hpo != "HEALTHY_CONTROL"), ]

# Load precomputed FDRs
del.fdrs <- load.fdrs(del.fdrs.in)
dup.fdrs <- load.fdrs(dup.fdrs.in)
fdrs <- list("DEL" = del.fdrs, "DUP" = dup.fdrs)

# Make one plot each for deletions and duplications
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(out.prefix, "fdr_by_hpo", cnv, "pdf", sep="."),
      height=2.25, width=2.5)
  plot.fdrs(fdrs[[cnv]], hpos, cnv, fdr.target, blue.bg=FALSE)
  dev.off()
})

