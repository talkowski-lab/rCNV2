#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Flexible script to generate UpSet plots for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load gene lists from tsv input
load.genes <- function(tsv.in){
  inputs <- read.table(tsv.in, header=F, sep="\t", comment.char="")
  genes <- lapply(inputs[, 2], function(path){
    unique(as.character(read.table(path, header=F)[, 1]))
  })
  names(genes) <- as.character(inputs[, 1])
  return(genes)
}

# Compute all possible pairwise intersections between lists of genes
intersect.sets <- function(gsets){
  setnames <- names(gsets)
  all.genes <- sort(unique(unlist(gsets)))
  gmat <- data.frame("gene"=all.genes)
  for(setname in setnames){
    gmat[setname] <- 0
    gmat[which(gmat$gene %in% gsets[[setname]]), setname] <- 1
  }
  counts <- sort(table(apply(gmat[, -1], 1, paste, collapse="")), decreasing=T)
  gmat <- as.data.frame(apply(cbind(counts, do.call("rbind", sapply(names(counts), strsplit, split=""))), 2, as.numeric))
  colnames(gmat) <- c("count", setnames)
  return(gmat)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Master UpSet plot
plot.upset <- function(gmat, min.highlight, light.color, dark.color,
                       top.bottom.ratio=3/2, button.cex=1.5, link.lwd=3,
                       parmar=c(0.2, 10, 1.5, 0.2)){
  # Get plot parameters
  n.rows <- ncol(gmat)-1
  n.cols <- nrow(gmat)
  top.height <- top.bottom.ratio * n.rows
  colors <- rep(light.color, n.cols)
  highlight.idxs <- which(apply(gmat[, -1], 1, sum) >= min.highlight)
  colors[highlight.idxs] <- dark.color
  
  # Transform bar heights
  count.max <- max(gmat$count, na.rm=T)
  count.scalar <- top.height / count.max
  bheights <- gmat$count * count.scalar
  
  # Get Y-axis scaling
  y.ax.log <- floor(log10(count.max))
  tick.space <- round(round(count.max, -y.ax.log) / 5, -y.ax.log+1)
  y.ax.lab <- seq(0, ceiling(count.max), tick.space)
  y.ax.at <- y.ax.lab * count.scalar
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, n.cols), ylim=c(-n.rows, top.height+1),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=0, col=bluewhite, bty="n", border=NA)
  # segments(x0=0:n.cols, x1=0:n.cols, y0=par("usr")[3], y1=0, col=blueblack, lend="butt")
  # abline(h=-1:-n.rows, col=blueblack, lend="butt")
  
  # Add buttons
  segments(x0=(1:n.cols)-0.5, x1=(1:n.cols)-0.5,
           y0=-n.rows+0.5, y1=-0.5, 
           lwd=link.lwd, col="white")
  sapply(1:n.cols, function(x){
    sapply(1:n.rows, function(y){
      points(x=x-0.5, y=-y+0.5, pch=19, cex=button.cex, col="white")
    })
  })
  
  # Add bars & labels
  rect(xleft=(1:n.cols)-0.9, xright=(1:n.cols)-0.1, ybottom=0, ytop=bheights,
       border=NA, bty="n", col=colors)
  text(x=(1:n.cols)-0.5, y=bheights, cex=5/6, xpd=T, srt=90, adj=c(0, 0.5),
       labels=paste(" ", prettyNum(gmat$count, big.mark=","), sep=""))
  
  # Add highlight arrows & count total number of genes
  lspace <- abs(0.0375*diff(par("usr")[3:4]))
  arrow.bufs <- sapply(strsplit(prettyNum(gmat$count[highlight.idxs], big.mark=","), split=""), length)
  arrow.starts <- bheights[highlight.idxs]+(lspace*(2+arrow.bufs))
  arrow.ends <- rep(max(bheights[highlight.idxs]+(lspace*(4+arrow.bufs))), length(highlight.idxs))
  Arrows(x0=highlight.idxs-0.5, x1=highlight.idxs-0.5,
         y0=arrow.starts, y1=arrow.ends, col=dark.color, lwd=2, 
         arr.type="triangle", code=1, arr.length=0.2, arr.width=0.15)
  segments(x0=min(highlight.idxs)-0.5, x1=max(highlight.idxs)-0.5,
           y0=max(arrow.ends), y1=max(arrow.ends), col=dark.color, lwd=2)
  text(x=mean(range(highlight.idxs) - 0.5), y=max(arrow.ends),
       labels=paste("N=", prettyNum(sum(gmat$count[highlight.idxs]), big.mark=","), "\ngenes", sep=""),
       col=dark.color, pos=3, cex=5.5/6, xpd=T)
  
  # Shade buttons
  sapply(1:n.cols, function(i){
    hits <- which(as.numeric(gmat[i, -1]) > 0)
    segments(x0=i-0.5, x1=i-0.5,
             y0=-min(hits)+0.5, y1=-max(hits)+0.5, 
             lwd=link.lwd, col=colors[i])
    points(x=rep(i-0.5, length(hits)), y=-hits+0.5, 
           pch=19, cex=button.cex, col=colors[i])
  })
  
  # Add axes
  abline(h=0, col=blueblack)
  axis(2, at=-(1:n.rows)+0.5, tick=F, line=-0.85, las=2, cex=5.5/6, labels=colnames(gmat)[-1])
  axis(2, at=c(0, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, tck=-0.025, labels=NA, col=blueblack)
  axis(2, at=y.ax.at[-1], tick=F, las=2, line=-0.65, labels=prettyNum(y.ax.lab[-1], big.mark=","), cex.axis=5.5/6)
  axis(2, at=par("usr")[4]/2, tick=F, line=2, label="Genes")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(shape, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
  make_option(c("--min-highlight"), default=2, metavar="integer",
              help="Minimum intersection size to highlight."),
  make_option(c("--dark-color"), default="black", metavar="character",
              help="Color to highlight --min-highlight categories."),
  make_option(c("--light-color"), default="gray70", metavar="character",
              help="Color for non-highlighted categories."),
  make_option(c("--cnv-coloring"), default=NULL, metavar="character", 
              help="Override color specifications based on rCNV color palettes"),
  make_option(c("--height"), default=2.25, help="Height of .pdf"),
  make_option(c("--width"), default=3.2, help="Width of .pdf")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog input.tsv output.pdf"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: input.tsv and output.pdf\n"))
}

# Writes args & opts to vars
tsv.in <- args$args[1]
outfile <- args$args[2]
rcnv.config <- opts$`rcnv-config`
min.highlight <- opts$`min-highlight`
dark.color <- opts$`dark-color`
light.color <- opts$`light-color`
cnv.coloring <- opts$`cnv-coloring`
pdf.height <- opts$`height`
pdf.width <- opts$`width`

# # DEV PARAMETERS
# tsv.in <- "~/scratch/hi.gs.upset.input.tsv"
# outfile <- "~/scratch/hi.gs.upset.pdf"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# min.highlight <- 2
# cnv.coloring <- "DEL"
# pdf.height <- 2.25
# pdf.width <- 3.2

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Set colors
if(!is.null(cnv.coloring)){
  dark.color <- cnv.colors[cnv.coloring]
  light.color <- control.cnv.colors[cnv.coloring]
}

# Load gene lists
gsets <- load.genes(tsv.in)

# Compute intersections
gmat <- intersect.sets(gsets)

# Plot
pdf(outfile, height=pdf.height, width=pdf.width)
plot.upset(gmat, min.highlight, light.color, dark.color)
dev.off()
