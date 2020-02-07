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
par(family="sans")
cnv.colors <- c("DEL"="#D43925", "DUP"="#2376B2")
logscale.size <- log10(as.vector(sapply(2:6, function(e){(1:9)*10^e})))
# logscale.size.labels.at <- log10(as.vector(sapply(2:6, function(e){c(1, 5)*10^e})))
# logscale.size.labels <- c("100bp", "500bp", "1kb", "5kb", "10kb", 
#                           "50kb", "100kb", "500kb", "1Mb", "5Mb")
logscale.size.labels.at <- 2:6
logscale.size.labels <- c("100bp", "1kb", "10kb", "100kb", "1Mb")
logscale.rate <- log10(as.vector(sapply(-4:2, function(e){(1:9)*10^e})))
logscale.rate.labels.at <- log10(as.vector(sapply(-3:2, function(e){c(1, 5)*10^e})))


# Plotting functions
add.labels <- function(x, y, labels, lab.cex=0.8){
  # Assign points to quadrants
  x.med <- median(x, na.rm=T)
  y.med <- median(y, na.rm=T)
  quads <- sapply(1:length(x), function(i){
    if(!is.na(x[i]) & !is.na(y[i])){
      if(y[i]>=y.med){
        yq <- "t"
      }else{
        yq <- "b"
      }
      if(x[i]>=x.med){
        xq <- "r"
      }else{
        xq <- "l"
      }
      paste(yq, xq, sep="")
    }else{
      return(NA)
    }
  })
  # Add labels
  sapply(1:length(x), function(i){
    if(!is.na(quads[i])){
      if(quads[i]=="tl"){
        text(x=x[i], y=y[i], pos=2, cex=lab.cex, labels=labels[i])
      }else if(quads[i]=="tr"){
        text(x=x[i], y=y[i], pos=4, cex=lab.cex, labels=labels[i])
      }else if(quads[i]=="bl"){
        text(x=x[i], y=y[i], pos=2, cex=lab.cex, labels=labels[i])
      }else if(quads[i]=="br"){
        text(x=x[i], y=y[i], pos=4, cex=lab.cex, labels=labels[i])
      }
    }
  })
}
scatter_stats <- function(dat, cnv, pt.cex=1.5,
                          xlims=c(1000, 1000000), ylims=c(0.01, 100)){
  # Prep plot values
  xlims <- log10(xlims)
  ylims <- log10(ylims)
  plot.dat <- data.frame("cohort"=dat$cohort,
                         "case_size"=log10(dat[, which(colnames(dat)==paste("med_case", cnv, "size", sep="_"))]),
                         "case_rate"=log10(dat[, which(colnames(dat)==paste(cnv, "per_case", sep="_"))]),
                         "ctrl_size"=log10(dat[, which(colnames(dat)==paste("med_ctrl", cnv, "size", sep="_"))]),
                         "ctrl_rate"=log10(dat[, which(colnames(dat)==paste(cnv, "per_ctrl", sep="_"))]))
  plot.col <- cnv.colors[which(names(cnv.colors)==cnv)]
  
  # Prep plot area
  par(mar=c(2.75, 3.25, 1.5, 1.5))
  plot(x=xlims, y=ylims, type="n",
       xaxt="n", xlab="", yaxt="n", ylab="")
  
  # Add points & labels
  # segments(x0=plot.dat$case_size, x1=plot.dat$ctrl_size,
  #          y0=plot.dat$case_rate, y1=plot.dat$ctrl_rate,
  #          lwd=2, col="gray90")
  add.labels(x=c(plot.dat$ctrl_size, plot.dat$case_size),
             y=c(plot.dat$ctrl_rate, plot.dat$case_rate),
             labels=c(plot.dat$cohort, plot.dat$cohort))
  points(x=plot.dat$ctrl_size, y=plot.dat$ctrl_rate,
         pch=22, lwd=1.5, cex=pt.cex, bg="white",col=plot.col)
  points(x=plot.dat$case_size, y=plot.dat$case_rate,
         pch=16, cex=pt.cex, col=plot.col)
  
  # Add axes
  axis(1, at=logscale.size, labels=NA, tck=-0.015, lwd=0.8)
  axis(1, at=logscale.size.labels.at, labels=NA, lwd=1.2, tck=-0.025)
  sapply(1:length(logscale.size.labels), function(x){
    axis(1, at=logscale.size.labels.at[x], line=-0.5, tick=F,
         labels=logscale.size.labels[x], cex.axis=0.9)
  })
  mtext(1, text=paste("Median", cnv, "Size"), line=1.75)
  axis(2, at=logscale.rate, labels=NA, tck=-0.015, lwd=0.8)
  axis(2, at=logscale.rate.labels.at, labels=NA, lwd=1.2, tck=-0.025)
  sapply(2:length(logscale.rate.labels.at), function(x){
    axis(2, at=logscale.rate.labels.at[x], line=-0.25, tick=F,
         labels=round(10^logscale.rate.labels.at[x], 2), cex.axis=0.9, las=2)
  })
  mtext(2, text=paste(cnv, "per Sample"), line=2)
  mtext(3, text=paste(cnv, "Dataset Properties"), line=0.1, font=2)
  
  # Add legend
  # legend("topright", legend=c("Case", "Control", "Same cohort"),
  #        pch=c(16, 22, NA), col=c(rep(plot.col, 2), "gray90"),
  #        lwd=c(NA, NA, 2), pt.lwd=c(NA, 1.5, NA), pt.cex=c(1.5, 1.5, NA),
  #        cex=0.9)
  legend("topright", legend=c("Case", "Control"),
         pch=c(16, 22), col=rep(plot.col, 2),
         pt.lwd=c(NA, 1.5), pt.cex=1.5, cex=0.9)
}
ratio_bars <- function(dat){
  # Prep plot values
  del.fracs <- c(dat$case_DEL_DUP_ratio/(1+dat$case_DEL_DUP_ratio),
                 dat$ctrl_DEL_DUP_ratio/(1+dat$ctrl_DEL_DUP_ratio))
  names(del.fracs) <- c(paste(dat$cohort, "(case)"),
                       paste(dat$cohort, "(control)"))
  del.fracs <- del.fracs[which(!is.na(del.fracs))]
  del.fracs <- del.fracs[order(names(del.fracs))]
  
  # Prep plot area
  par(mar=c(2.75, 6, 1.5, 1.5), bty="n")
  plot(x=c(0, 1), y=c(0, -length(del.fracs)), type="n",
       xaxt="n", xlab="", yaxt="n", ylab="", yaxs="i")
  abline(v=0.5)
  
  # Add bars
  rect(xleft=0, xright=1, 
       ybottom=-(1:length(del.fracs))+0.15,
       ytop=-(1:length(del.fracs))+0.85,
       border=NA, bty="n", col=cnv.colors[2])
  rect(xleft=0, xright=del.fracs, 
       ybottom=-(1:length(del.fracs))+0.15,
       ytop=-(1:length(del.fracs))+0.85,
       border=NA, bty="n", col=cnv.colors[1])
  
  # Add axes
  axis(1, at=seq(0, 1, 0.25), tck=-0.025, labels=NA)
  axis(1, at=seq(0, 1, 0.25), tick=F, line=-0.5, cex.axis=0.9,
       labels=paste(seq(0, 100, 25), "%", sep=""))
  mtext(1, text="Fraction of CNVs", line=1.75)
  axis(2, at=-(1:length(del.fracs))+0.5, tick=F, line=-0.9,
       labels=names(del.fracs), cex.axis=0.9, las=2)
  mtext(3, text="DEL:DUP Ratio", line=0.1, font=2)
  
  # Add legend
  text(x=c(0, 1), y=-0.5, pos=c(4, 2),
       col="white", labels=c("DEL", "DUP"), cex=0.9)
}


# Read arguments
args <- commandArgs(trailingOnly=T)
tsv.in <- as.character(args[1])
plot.out <- as.character(args[2])


# Read data
dat <- read.table(tsv.in, header=T, sep="\t")


# Plot data
jpeg(plot.out, height=3.25*300, width=10*300, res=300)
par(mfrow=c(1, 3))
scatter_stats(dat, "DEL")
scatter_stats(dat, "DUP")
ratio_bars(dat)
dev.off()
