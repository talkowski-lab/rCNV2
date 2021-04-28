#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Barplot of samples per HPO term per metacohort


# Set parameters
options(scipen=1000, stringsAsFactors=F)
meta.colors <- c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#984464", "#C0AFFB")


# Plotting functions
plot.labels <- function(labels){
  par(bty="n")
  plot(x=c(0, 1), y=c(0, length(labels)), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  text(x=par("usr")[2], y=(1:length(labels))-0.5, 
       pos=2, xpd=T, labels=rev(labels))
}
total.barplot <- function(totals){
  par(bty="n")
  plot(x=c(0, 1.25*max(totals)), y=c(0, length(totals)), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=0, xright=rev(totals), ybottom=(1:length(totals))-0.15, ytop=(1:length(totals))-0.85,
       col="gray15")
  text(x=rev(totals), y=(1:length(totals))-0.5, pos=4, xpd=T, 
       labels=prettyNum(rev(totals), big.mark=","))
  axis(3, labels=NA, line=-1.25)
  axis(3, at=axTicks(3), line=-1.25, labels=paste(axTicks(3)/1000, "k", sep=""))
}
fraction.barplot <- function(counts, colors){
  fracs <- t(sapply(1:nrow(counts), function(i){
    counts[i, grep("meta", colnames(counts))]/counts$Total[i]
  }))
  par(bty="n")
  plot(x=c(0, 1), y=c(0, nrow(fracs)), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  sapply(1:nrow(fracs), function(i){
    rect(xleft=c(0, cumsum(fracs[i, -ncol(fracs)])), xright=cumsum(fracs[i, ]),
    ybottom=nrow(fracs)-i+0.15, ytop=nrow(fracs)-i+0.85, col=colors)
  })
  axis(3, at=seq(0, 1, 0.2), labels=NA, line=-1.25)
  axis(3, at=seq(0, 1, 0.2), line=-1.25, 
       labels=paste(seq(0, 100, 20), "%", sep=""))
}


# Read arguments
args <- commandArgs(trailingOnly=T)
tsv.in <- as.character(args[1])
plot.out <- as.character(args[2])


# Read data
counts <- read.table(tsv.in, comment.char="", header=T, sep="\t")
counts$label <- paste(counts[, 1], " (", counts$description, ")", sep="")


# Plot data
jpeg(plot.out, res=300, height=7*300, width=9*300)
layout(matrix(1:3, nrow=1), widths=c(3, 2, 2))
par(bty="n", mar=c(0.25, 0.25, 2.5, 0.25))
plot.labels(counts$label)
mtext(3, line=1, text="HPO term", font=2)
fraction.barplot(counts, colors=meta.colors)
mtext(3, line=1, text="Metacohort composition", font=2)
total.barplot(counts$Total)
mtext(3, line=1, text="Total samples", font=2)
legend("right", legend=paste("meta", 1:6, sep=""),
       fill=meta.colors)
dev.off()
