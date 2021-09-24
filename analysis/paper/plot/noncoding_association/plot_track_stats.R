#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distributions of raw genome annotation tracks for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


##########################
### PLOTTING FUNCTIONS ###
##########################
# Simple horizontal barplot of tracks per family
family.barplot <- function(stats, parmar=c(0.5, 12, 2.25, 2)){
  # Get plot data
  pdat <- sort(table(stats$family[which(stats$cnv=="DEL")]), decreasing=TRUE)
  n.bars <- length(pdat)
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 1.1*max(pdat)), ylim=c(n.bars, 0),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  
  # Add top axis & gridlines
  x.ax.at <- 1000*round(axTicks(1)/1000, 0)
  # abline(v=x.ax.at, col=bluewhite)
  axis(3, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(3, at=x.ax.at, tck=-0.025, labels=NA, col=blueblack)
  sapply(x.ax.at, function(x){
    axis(3, at=x, tick=F, line=-0.8, labels=round(x/1000, 1))
  })
  mtext(3, text="Annotation Classes (x1,000)", line=1.25)
  
  # Add rectangles & labels
  rect(xleft=rep(0, n.bars), xright=pdat,
       ybottom=(1:n.bars)-0.15, ytop=(1:n.bars)-0.85,
       col=nc.anno.family.colors[names(pdat)],
       border=blueblack)
  text(x=pdat, y=(1:n.bars)-0.5, pos=4, cex=5/6, xpd=TRUE,
       labels=prettyNum(pdat, big.mark=","))
  axis(2, at=(1:n.bars)-0.5, tick=F, line=-0.8, las=2, cex=5/6,
       labels=nc.anno.family.names[names(pdat)])
  
  # Add cleanup line to left Y-axis
  axis(2, at=c(-10e10, 10e10), labels=NA, tck=0, col=blueblack)
}

# Scatterplot of mean element size vs. number of elements in class
track.scatter <- function(stats, pt.cex=0.2, blue.bg=TRUE, parmar=c(2.25, 2.5, 0.25, 1)){
  # Get plot data
  plot.df <- stats[which(stats$cnv=="DEL"), c("mean_size", "n_elements", "family")]
  plot.df <- plot.df[which(plot.df$n_elements>0), ]
  set.seed(2020)
  plot.df <- plot.df[sample(1:nrow(plot.df), nrow(plot.df), replace=F), ]
  plot.df[, 1:2] <- apply(plot.df[, 1:2], 2, log10)
  xlims <- range(plot.df[, 1], na.rm=T)
  ylims <- range(plot.df[, 2], na.rm=T)
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
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=ylims, xaxt="n", yaxt="n", xlab="", ylab="")
  ax.at <- -10:10
  
  # Add background shading
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty=plot.bty, border=plot.border, col=plot.bg)
  abline(h=ax.at, v=ax.at, col=grid.col)
  
  # Add points
  points(plot.df[, 1:2], pch=19, cex=pt.cex, col=nc.anno.family.colors[plot.df$family])

  # Add axes
  axis(1, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(1, at=log10(logscale.minor[which(logscale.minor >= 1)]), tck=-0.015, col=blueblack, labels=NA, lwd=0.7)
  axis(1, at=ax.at, tck=-0.025, col=blueblack, labels=NA)
  sapply(1:length(logscale.major.bp), function(i){
    axis(1, at=log10(logscale.major.bp[i]), tick=F, line=-0.8, cex.axis=5.5/6, 
         labels=logscale.major.bp.labels[i])
  })
  mtext(1, line=1.2, text="Mean Element Size")
  axis(2, at=log10(logscale.minor[which(logscale.minor >= 1)]), tck=-0.015, col=blueblack, labels=NA, lwd=0.7)
  axis(2, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(2, at=ax.at, tck=-0.025, col=blueblack, labels=NA)
  sapply(ax.at, function(i){
    axis(2, at=i, tick=F, line=-0.65, las=2, cex.axis=5.5/6,
         labels=bquote(10^.(i)))
  })
  mtext(2, line=1.6, text="Elements in Class")
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
args <- parse_args(OptionParser(usage=paste("%prog stats.tsv out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: stats.tsv and out.prefix\n"))
}

# Writes args & opts to vars
stats.in <- args$args[1]
out.prefix <- args$args[2]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# stats.in <- "~/scratch/rCNV.burden_stats.tsv.gz"
# out.prefix <- "~/scratch/track_stats_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/noncoding_association/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load track stats & subset to single CNV type (no double counting)
stats <- load.track.stats(stats.in)
stats <- stats[which(stats$cnv=="DEL"), ]

# Barplot of number of tracks per annotation family
pdf(paste(out.prefix, "track_stats.barplot.pdf", sep="."),
    height=2.5, width=4.1)
family.barplot(stats)
dev.off()

# Scatterplot of mean element size vs number of elements in track
pdf(paste(out.prefix, "track_stats.scatterplot.pdf", sep="."),
    height=2.75, width=2.9)
track.scatter(stats, blue.bg=FALSE)
dev.off()

