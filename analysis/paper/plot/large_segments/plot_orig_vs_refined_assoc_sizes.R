#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot comparison of original significant region vs. credble set for rCNV2 large segment final analyses


options(stringsAsFactors=F, scipen=1000, check.names=F)


######################
### DATA FUNCTIONS ###
######################
# Append old size onto segments
load.sizes <- function(tsv.in, segs){
  x <- read.table(tsv.in, header=T, sep="\t", comment.char="")
  merge(segs, x[, c("region_id", "original_size")],
                by="region_id", all.x=T, all.y=F, sort=F)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Scatterplot of original vs refined sizes
sizes.scatter <- function(sizes, blue.bg=TRUE, parmar=c(3.2, 3.2, 0.75, 0.75)){
  # Get plot data
  ax.lims <- range(c(sizes$original_size, sizes$refined_size), na.rm=T)
  ax.lims <- c(floor(ax.lims[1]), ceiling(ax.lims[2]))
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
  plot(NA, xlim=ax.lims, ylim=ax.lims,
       xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty=plot.bty, border=plot.border, col=plot.bg)
  abline(h=log10(logscale.major.bp), v=log10(logscale.major.bp), col=grid.col)
  abline(0, 1, col=blueblack)

  # Add points
  points(sizes$original_size, sizes$refined_size,
         pch=sizes$pt.pch, bg=sizes$pt.bg, col=sizes$pt.border)

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
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog sizes.tsv segs.tsv output_prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Three positional arguments required: sizes.tsv, segs.tsv and output_prefix\n")
}

# Writes args & opts to vars
tsv.in <- args$args[1]
segs.in <- args$args[2]
out.prefix <- args$args[3]

# # DEV PARAMETERS
# tsv.in <- "~/scratch/rCNV2_analysis_d2.segments.old_vs_new_size.tsv.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# out.prefix <- "~/scratch/new_old_sizes"

# Load table of all segments & annotate with original size
segs <- load.segment.table(segs.in)
segs <- segs[which(segs$any_sig), ]
segs <- load.sizes(tsv.in, segs)

# Plot correlation of sizes
pdf(paste(out.prefix, "original_vs_refined_sizes.pdf", sep="."),
    height=2.5, width=2.8)
segs.scatter(segs, x=log10(segs$original_size), y=log10(segs$size),
             xlims=log10(c(100000, max(segs$size))),
             ylims=log10(c(100000, max(segs$size))),
             blue.bg=FALSE, add.cor=F, add.lm=F, abline.a=0, abline.b=1, pt.cex=0.65,
             x.at=log10(logscale.demi.bp), x.labs.at=log10(logscale.major.bp), x.labs=logscale.major.bp.labels,
             y.at=log10(logscale.demi.bp), y.labs=logscale.demi.bp.labels,
             xtitle=bquote(log[10]("Original Size")), x.title.line=1.55,
             ytitle=bquote(log[10]("Refined Size")), y.title.line=2.65,
             parmar=c(2.5, 3.9, 0.5, 1))
dev.off()
