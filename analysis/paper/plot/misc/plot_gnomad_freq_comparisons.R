#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot CNV frequency comparisons vs. gnomAD-SV for a single cohort


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load & deduplicate gnomAD comparison data
load.comp.dat <- function(data.in){
  dat <- read.table(data.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(dat)[1] <- gsub("#","", colnames(dat)[1], fixed=T)
  do.call("rbind", lapply(unique(dat$vid), function(vid){
    hits <- dat[which(dat$vid == vid), ]
    if(nrow(hits) > 1){
      hits$d <- abs(hits$freq - hits$gnomad)
      keeper <- with(hits, order(-overlap, d))[1]
      hits$d <- NULL
      hits <- hits[keeper, ]
    }
    return(hits)
  }))
}

# Group data for plotting
get.plot.dat <- function(dat){
  pdat <- lapply(c("DEL", "DUP"), function(cnv){
    cnv.idxs <- which(dat$cnv == cnv)
    res <- list(dat$freq[intersect(which(is.na(dat$gnomad)), cnv.idxs)],
                dat$freq[intersect(which(dat$gnomad < 1/10000), cnv.idxs)],
                dat$freq[intersect(which(dat$gnomad >= 1/10000 & dat$gnomad < 1/1000), cnv.idxs)],
                dat$freq[intersect(which(dat$gnomad >= 1/1000), cnv.idxs)])
  })
  names(pdat) <- c("DEL", "DUP")
  return(pdat)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot size distributions for a single CNV type
plot.size.dist <- function(dat, cnv, parmar=c(3, 3, 0.2, 0.5)){
  # Get plot data
  dat$size <- dat$end - dat$start
  h <- hist(log10(dat$size[which(dat$cnv==cnv)]), plot=F)
  h$density <- h$counts / sum(h$counts)

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(5, log10(20000000)), ylim=c(0, max(h$density)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")

  # Add histogram
  rect(xleft=h$breaks[-length(h$breaks)], xright=h$breaks[-1],
       ybottom=0, ytop=h$density, border=NA, col=cnv.colors[cnv])

  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=log10(logscale.major.bp), tck=-0.025, col=blueblack, labels=NA)
  x.lab.at <- log10(logscale.major.bp)
  x.lab.at.keep <- which(x.lab.at >= par("usr")[1] & x.lab.at <= par("usr")[2])
  text(x=x.lab.at[x.lab.at.keep]+(0.075*diff(par("usr")[1:2])),
       y=par("usr")[3]-(0.05*diff(par("usr")[3:4])),
       labels=logscale.major.bp.labels[x.lab.at.keep],
       srt=45, pos=2, xpd=T, cex=5/6)
  mtext(1, line=2, text="CNV Size")
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, tck=-0.025, labels=NA, col=blueblack)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=5/6, line=-0.7,
       labels=paste(round(100 * axTicks(2), 0), "%", sep=""))
  mtext(2, line=2, text=paste("Fraction of ", cnv, "s", sep=""))
}

# Plot frequency distributions for a single CNV type
plot.freq.dist <- function(dat, cnv, parmar=c(3, 3, 0.2, 0.5)){
  # Get plot data
  h <- hist(log10(dat$freq[which(dat$cnv==cnv)]), plot=F)
  h$density <- h$counts / sum(h$counts)

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(min(h$breaks), -2), ylim=c(0, max(h$density)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")

  # Add histogram
  rect(xleft=h$breaks[-length(h$breaks)], xright=h$breaks[-1],
       ybottom=0, ytop=h$density, border=NA, col=cnv.colors[cnv])

  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=log10(logscale.major), tck=-0.025, col=blueblack, labels=NA)
  x.lab.at <- log10(logscale.major)
  x.lab.at <- x.lab.at[which(x.lab.at >= par("usr")[1] & x.lab.at <= par("usr")[2])]
  text(x=x.lab.at+(0.075*diff(par("usr")[1:2])),
       y=par("usr")[3]-(0.05*diff(par("usr")[3:4])),
       labels=paste(100 * (10^x.lab.at), "%", sep=""),
       srt=45, pos=2, xpd=T, cex=5/6)
  mtext(1, line=2, text="Frequency in Cohort")
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, tck=-0.025, labels=NA, col=blueblack)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=5/6, line=-0.7,
       labels=paste(round(100 * axTicks(2), 0), "%", sep=""))
  mtext(2, line=2, text=paste("Fraction of ", cnv, "s", sep=""))
}

# Plot summary results from gnomAD comparison for a single CNV type
plot.res <- function(pdat, cnv, ymax=-2, swarm.max=50, bar.hex=0.3,
                     parmar=c(3.6, 3.5, 0.1, 0.1)){
  # Get plot dimensions & data
  subdat <- lapply(pdat[[cnv]], log10)
  n.x <- length(subdat)
  ymin <- min(unlist(subdat), na.rm=T)
  yrange <- abs(ymin - ymax)
  ymax.new <- ymax + (bar.hex * yrange)
  counts <- sapply(subdat, length)
  counts <- (counts / max(counts, na.rm=T)) * bar.hex * yrange

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, n.x), ylim=c(ymin, ymax.new),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  abline(h=log10(1/10847), lty=5, col=blueblack)

  # Add categorical X axis
  axis(1, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(1, at=(1:n.x)-0.5, tick=F, las=2, cex.axis=5/6, line=-0.9,
       labels=c("Not\nFound", "<0.01%", "0.01-\n0.1%", "0.1-\n1%"))
  mtext(1, line=2.6, text="Frequency in gnomAD")

  # Add lower Y axis
  axis(2, at=c(-10e10, ymax), col=blueblack, tck=0, labels=NA)
  lower.y.tick.at <- log10(logscale.major)[which(logscale.major <= 10^ymax)]
  axis(2, at=lower.y.tick.at, col=blueblack, tck=-0.025, labels=NA)
  sapply(lower.y.tick.at, function(y){
    axis(2, at=y, tick=F, las=2, cex.axis=5/6, line=-0.7,
         labels=paste(100 * (10^y), "%", sep=""))
  })
  axis(2, at=mean(c(ymin, ymax)), tick=F, line=1.6, labels="Frequency in Cohort")

  # Add vioplots (or beeswarms, if optioned)
  sapply(1:n.x, function(x){
    if(length(subdat[[x]]) <= swarm.max) {
      beeswarm(subdat[[x]], add=T, at=x-0.5, cex=1/3, col=cnv.colors[cnv], pch=19,
               corral="wrap", corralWidth=5/6)
    }else{
      vioplot(subdat[[x]], add=T, at=x-0.5, col=cnv.colors[cnv], drawRect=F)
    }
    points(x=x-0.5, y=mean(subdat[[x]]), pch=23, bg="white", col=cnv.blacks[cnv],
           lwd=2, cex=1.25)
  })

  # Add upper panel
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=ymax, ytop=ymax.new,
       bty="o", col="white", border="white")
  abline(h=ymax, col=blueblack)
  rect(xleft=(1:n.x)-(5/6), xright=(1:n.x)-(1/6), ybottom=ymax, ytop=ymax+counts,
       bty="n", border=NA, col=cnv.colors[cnv])
  axis(2, at=(ymax + (3*ymax.new))/4, las=2, line=-1.2,
       label=paste("# ", cnv, "s", sep=""), tick=F, col.axis=cnv.colors[cnv])
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(vioplot, quietly=T)
require(beeswarm, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog data.tsv out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: data.tsv and out.prefix\n"))
}

# Writes args & opts to vars
data.in <- args$args[1]
out.prefix <- args$args[2]

# # DEV PARAMETERS
# data.in <- "~/scratch/meta1.cnv_vs_gnomad.bed.gz"
# out.prefix <- "~/scratch/rCNV_gnomAD_comparison_test"

# Load & transform data
dat <- load.comp.dat(data.in)
pdat <- get.plot.dat(dat)

# Plot size histograms for deletions and duplications separately
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(out.prefix, "size_distrib", cnv, "pdf", sep="."),
      height=2, width=2.25)
  plot.size.dist(dat, cnv)
  dev.off()
})

# Plot frequency histograms for deletions and duplications separately
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(out.prefix, "freq_distrib", cnv, "pdf", sep="."),
      height=2, width=2.25)
  plot.freq.dist(dat, cnv)
  dev.off()
})

# Plot gnomAD comparisons for deletions and duplications separately
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(out.prefix, "vs_gnomAD_bins", cnv, "pdf", sep="."),
      height=2.75, width=2.4)
  plot.res(pdat, cnv)
  dev.off()
})

