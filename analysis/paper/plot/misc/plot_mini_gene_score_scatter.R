#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot mini version of pHaplo & pTriplo scores for graphical abstract of rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Prepare color grid for gradient-based scatter
make.color.grid <- function(){
  yl.pal <- colorRampPalette(c(ns.color.light, control.cnv.colors[2], cnv.colors[2]))(101)
  yr.pal <- colorRampPalette(c(cnv.colors[1],  cnv.colors[3]))(101)
  # xb.pal <- colorRampPalette(c(ns.color, control.cnv.colors[1], cnv.colors[1]))(101)
  # xt.pal <- colorRampPalette(c(cnv.colors[2], cnv.colors[3]))(101)
  do.call("rbind", lapply(1:101, function(y){
    colorRampPalette(c(yl.pal[y], yr.pal[y]))(101)
  }))
}

# Get gene color based on (pHaplo, pTriplo) x,y pairs
query.color.grid <- function(gene, scores, color.grid){
  x <- round(100*as.numeric(scores$pHaplo[which(scores$gene==gene)]))
  y <- round(100*as.numeric(scores$pTriplo[which(scores$gene==gene)]))
  color.grid[y+1, x+1]
}

# Load scores and assign colors
load.scores <- function(scores.in, gradient.color=FALSE){
  scores <- read.table(scores.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(scores)[1] <- gsub("#", "", colnames(scores)[1], fixed=T)
  scores[, -1] <- apply(scores[, -1], 2, as.numeric)
  if(gradient.color==TRUE){
    color.grid <- make.color.grid()
    scores$color <- as.character(sapply(scores$gene, query.color.grid, scores, color.grid))
  }else{
    scores$color <- apply(scores[, -1], 1, function(vals){
      if(vals[1]>=0.9 & vals[2]<0.9){
        cnv.colors[which(names(cnv.colors) == "DEL")]
      }else if(vals[1]<0.9 & vals[2]>=0.9){
        cnv.colors[which(names(cnv.colors) == "DUP")]
      }else if(vals[1]>=0.9 & vals[2]>=0.9){
        cnv.colors[which(names(cnv.colors) == "CNV")]
      }else{
        "gray70"
      }
    })
  }
  return(scores)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
plot.mini.scatter <- function(scores, n.points=1000, shading=TRUE, pt.cex=0.25, seed=2020){
  # Downsample points
  set.seed(seed)
  scores <- scores[sample(1:nrow(scores), n.points, replace=F), ]
  # Prep plot area
  par(mar=c(1.1, 1.1, 0.2, 0.2), bty="n", bg=NA)
  plot(NA, xlim=c(0, 1), ylim=c(0, 1),
       xaxs="i", xlab="", xaxt="n",
       yaxs="i", ylab="", yaxt="n")
  # Add shading
  if(shading==TRUE){
    rect(xleft=par("usr")[1], xright=0.9,
         ybottom=par("usr")[3], ytop=0.9,
         bty="n", border=NA, col=bluewhite)
    rect(xleft=par("usr")[1], xright=0.9,
         ybottom=0.9, ytop=par("usr")[4],
         bty="n", border=NA, col=control.cnv.colors[2])
    rect(xleft=0.9, xright=par("usr")[2],
         ybottom=0.9, ytop=par("usr")[4],
         bty="n", border=NA, col=control.cnv.colors[3])
    rect(xleft=0.9, xright=par("usr")[2],
         ybottom=par("usr")[3], ytop=0.9,
         bty="n", border=NA, col=control.cnv.colors[1])
  }else{
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=par("usr")[3], ytop=par("usr")[4],
         border=NA, bty="n", col=bluewhite)
    abline(h=seq(0, 1, 0.2), v=seq(0, 1, 0.2), col="white")
  }
  
  # Add points
  points(x=scores$pHaplo, y=scores$pTriplo, pch=19, cex=pt.cex, 
         col=scores$color)
  # Axes
  # axis(2, labels=NA, tck=-0.03, col=blueblack)
  # axis(3, labels=NA, tck=-0.03, col=blueblack)
  mtext(2, line=0.1, text="Duplication", col=blueblack)
  mtext(1, line=0.1, text="Deletion", col=blueblack)
  
  # Cleanup
  if(shading==TRUE){
    abline(h=0.9, col=blueblack, lty=2)
    abline(v=0.9, col=redblack, lty=2)
  }
  box(bty="o", xpd=T, col=blueblack)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv out.png", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: scores.tsv and output.png\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
out.png <- args$args[2]

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# out.png <- "~/scratch/test_mini_scatter.png"

# Load scores
scores <- load.scores(scores.in, gradient.color=TRUE)

# Plot scores
png(out.png, height=1.15*300, width=1.15*300, res=300, family="sans")
plot.mini.scatter(scores, n.points=1250, pt.cex=0.25, shading=FALSE)
dev.off()
