#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Common functions for gene dosage sensitivity score analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load dataframe of all gene scores
load.scores <- function(scores.in){
  scores <- read.table(scores.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(scores) <- c("gene", "pHI", "pTS")
  scores[, -1] <- apply(scores[, -1], 2, as.numeric)
  return(scores)
}

# Partition genes into subgroups based on scores
classify.genes <- function(scores, hc.cutoff=0.9, lc.cutoff=0.5){
  ds.hc <- scores$gene[which(scores$pHI>=hc.cutoff & scores$pTS>=hc.cutoff)]
  ds.lc <- scores$gene[which(scores$pHI>=lc.cutoff & scores$pTS>=lc.cutoff & !(scores$gene %in% ds.hc))]
  hi.hc <- scores$gene[which(scores$pHI>=hc.cutoff & scores$pTS<lc.cutoff)]
  hi.lc <- scores$gene[which(scores$pHI>=lc.cutoff & scores$pTS<lc.cutoff & !(scores$gene %in% hi.hc))]
  ts.hc <- scores$gene[which(scores$pHI<lc.cutoff & scores$pTS>=hc.cutoff)]
  ts.lc <- scores$gene[which(scores$pHI<lc.cutoff & scores$pTS>=lc.cutoff & !(scores$gene %in% ts.hc))]
  ns <- scores$gene[which(scores$pHI<lc.cutoff & scores$pTS<lc.cutoff)]
  list("ds.hc" = ds.hc, "ds.lc" = ds.lc,
       "hi.hc" = hi.hc, "hi.lc" = hi.lc,
       "ts.hc" = ts.hc, "ts.lc" = ts.lc,
       "ns" = ns)
}

# Get gene color based on membership in ds.groups
get.gene.color.byscore <- function(gene, ds.groups){
  if(gene %in% ds.groups$ds.hc){
    cnv.colors[3]
  }else if(gene %in% ds.groups$ds.lc){
    control.cnv.colors[3]
  }else if(gene %in% ds.groups$hi.hc){
    cnv.colors[1]
  }else if(gene %in% ds.groups$hi.lc){
    control.cnv.colors[1]
  }else if(gene %in% ds.groups$ts.hc){
    cnv.colors[2]
  }else if(gene %in% ds.groups$ts.lc){
    control.cnv.colors[2]
  }else{
    ns.color
  }
}

# Load gene features
load.features <- function(features.in, norm=F){
  feats <- read.table(features.in, header=T, sep="\t", comment.char="", check.names=F)[, -c(1:3)]
  if(norm==T){
    feats[, -1] <- apply(feats[, -1], 2, function(vals){
      scale(as.numeric(vals, scale=T, center=T))
    })
  }
  return(feats)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot distributions of values by gene score category
plot.feature.bydsgroup <- function(feats, ds.groups, feat.idx=2,
                                   parmar=c(2.25, 2.5, 1.25, 0.5)){
  # Get plot values
  vals <- lapply(ds.groups, function(genes){
    as.numeric(feats[which(feats$gene %in% genes), feat.idx])
  })
  names(vals) <- names(ds.groups)
  vals <- list("ns"=vals$ns,
               "hi.hc"=vals$hi.hc, "ds.hc"=vals$ds.hc, "ts.hc"=vals$ts.hc, 
               "hi.lc"=vals$hi.lc, "ds.lc"=vals$ds.lc, "ts.lc"=vals$ts.lc)
  ylims <- range(unlist(vals), na.rm=T)
  x.at <- c(1, 3:5, 7:9)-0.5
  plot.colors <- c(ns.color, cnv.colors[c(1, 3, 2)], control.cnv.colors[c(1, 3, 2)])
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 9), ylim=ylims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty="n", border=NA, col=bluewhite)
  abline(h=median(vals[[1]], na.rm=T), lty=2, col=blueblack)
  
  # Add axes
  x.labs <- c("NS", rep(c("HI", "DS", "TS"), 2))
  sapply(1:7, function(i){
    axis(1, at=x.at[i], tick=F, line=-0.9,
         labels=x.labs[i])
  })
  axis(1, at=c(2.1, 4.9), tck=0, col=blueblack, labels=NA, line=1.1)
  axis(1, at=c(6.1, 8.9), tck=0, col=blueblack, labels=NA, line=1.1)
  axis(1, at=3.5, line=0.1, labels="High Conf.", tick=F)
  axis(1, at=7.5, line=0.1, labels="Low Conf.", tick=F)
  axis(2, at=c(-10e10, 10e10), tck=0, col=blueblack)
  axis(2, tck=-0.03, labels=NA, col=blueblack)
  axis(2, tick=F, las=2, line=-0.65)
  mtext(2, text="Z-Score", line=1.25)
  mtext(3, text=colnames(feats)[feat.idx], font=2)
  
  # Add values
  sapply(1:7, function(i){
    ivals <- vals[[i]][which(!is.na(vals[[i]]))]
    if(length(unique(ivals)) > 2){
      vioplot(ivals, at=x.at[i], col=plot.colors[i], border=blueblack,
              add=T, drawRect=F)
      vmed <- median(ivals)
      segments(x0=x.at[i]+0.2, x1=x.at[i]-0.2, y0=vmed, y1=vmed, 
               lend="round", col=blueblack, lwd=2)
    }
  })
}

