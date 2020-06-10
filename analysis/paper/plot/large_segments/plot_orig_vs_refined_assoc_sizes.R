#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot comparison of original significant region vs. credble set for rCNV2 large segment final analyses


options(stringsAsFactors=F, scipen=1000, check.names=F)


######################
### DATA FUNCTIONS ###
######################
# Load table of new/old sizes
load.sizes <- function(tsv.in){
  x <- read.table(tsv.in, header=T, sep="\t", comment.char="")
  x$bg <- as.character(cnv.colors[x$cnv])
  x$border <- as.character(cnv.blacks[x$cnv])
  x$original_size <- log10(x$original_size)
  x$refined_size <- log10(x$refined_size)
  return(x)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Scatterplot of original vs refined sizes
sizes.scatter <- function(sizes, parmar=c(3.2, 3.2, 0.75, 0.75)){
  # Get plot data
  ax.lims <- range(c(sizes$original_size, sizes$refined_size), na.rm=T)
  ax.lims <- c(floor(ax.lims[1]), ceiling(ax.lims[2]))
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=ax.lims, ylim=ax.lims,
       xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty="n", border=NA, col=bluewhite)
  abline(h=log10(logscale.major.bp), v=log10(logscale.major.bp), col="white")
  abline(0, 1, col=blueblack)
  
  # Add points
  points(sizes$original_size, sizes$refined_size, pch=22, bg=sizes$bg, col=sizes$border)
  
  # Add axes
  sapply(1:2, function(i){
    axis(i, at=log10(logscale.minor), labels=NA, tck=-0.015, col=blueblack)
    axis(i, at=log10(logscale.major), labels=NA, tck=-0.03, col=blueblack)  
  })
  sapply(1:length(logscale.major.bp), function(i){
    axis(1, at=log10(logscale.major.bp)[i], labels=logscale.major.bp.labels[i],
         tick=F, line=-0.65)
    axis(2, at=log10(logscale.major.bp)[i], labels=logscale.major.bp.labels[i],
         tick=F, line=-0.65, las=2)
  })
  mtext(1, text="Original size", line=1.25)
  mtext(2, text="Credible interval size", line=2.25)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog sizes.tsv output_prefix", option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Two positional arguments required: sizes.tsv and output_prefix\n")
}

# Writes args & opts to vars
tsv.in <- args$args[1]
out.prefix <- args$args[2]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# tsv.in <- "~/scratch/rCNV2_analysis_d1.associations.old_vs_new_size.tsv.gz"
# out.prefix <- "~/scratch/new_old_sizes"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load table of sizes
sizes <- load.sizes(tsv.in)

# Plot correlation of sizes
pdf(paste(out.prefix, "original_vs_refined_sizes.pdf", sep="."),
    height=2.65, width=2.65)
sizes.scatter(sizes)
dev.off()
