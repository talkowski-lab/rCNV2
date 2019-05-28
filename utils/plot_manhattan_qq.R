#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Generate Manhattan and QQ plots for an input set of association statistics


options(scipen=1000, stringsAsFactors=F)


# Read an input file of association statistics
read.stats <- function(stats.in, p.col.name="p", p.is.phred=F){
  stats <- read.table(stats.in, header=T, comment.char="", sep="\t")
  
  # Get coordinates
  chr <- as.character(stats[,1])
  if("pos" %in% colnames(stats)){
    pos <- stats$pos
  }else if(all(c("start", "end") %in% colnames(stats))){
    pos <- (stats$start + stats$end) / 2
  }else{
    stop("Unable to identify locus coordinate info. Must supply columns for either 'pos' or 'start' and 'end'.")
  }
  pos <- as.numeric(pos)
  
  # Get untransformed p-values
  if(p.col.name %in% colnames(stats)){
    p <- stats[, which(colnames(stats) == p.col.name)]
    p <- as.numeric(p)
    if(p.is.phred == T){
      p <- 10^-p
    }
  }else{
    stop(paste("Unable to identify p-value column by header name, ",
               p.col.name, ".", sep=""))
  }
  
  data.frame("chr"=chr, "pos"=pos, "p"=p)
}


# Read intervals to highlight
read.highlight.bed <- function(highlight.in){
  bed <- read.table(highlight.in, sep="\t", header=F, comment.char="#")[, 1:3]
  colnames(bed) <- c("chr", "start", "end")
  return(bed)
}


# Manhattan plot
manhattan <- function (df, cutoff=1e-08, highlights=NULL, 
                       highlight.color="#4EE69A",
                       highlight.name="Positive Controls"){
  colors <- c("gray15", "gray65")
  
  contigs <- unique(df[, 1])
  contigs <- contigs[which(!(is.na(contigs)))]
  
  indexes <- as.data.frame(t(sapply(contigs, function(chr){
    return(c(chr, 0, max(df[which(df[, 1] == chr), 2])))
  })))
  indexes$sum <- cumsum(indexes[, 3])
  indexes$bg <- rep(colors[1:2], ceiling(nrow(indexes)/2))[1:nrow(indexes)]
  indexes[, 2:4] <- apply(indexes[, 2:4], 2, as.numeric)
  
  df.plot <- as.data.frame(t(apply(df, 1, function(row){
    return(c(row[1], 
             as.numeric(row[2]) + indexes[as.numeric(row[1]), 4] - indexes[as.numeric(row[1]), 3], 
             as.numeric(row[3]), 
             indexes[as.numeric(row[1]), 5]))
  })))
  df.plot[, 2] <- as.numeric(as.character(df.plot[, 2]))
  df.plot[, 3] <- as.numeric(as.character(df.plot[, 3]))
  colnames(df.plot) <- c("chr", "pos", "p", "color")
  
  if(!is.null(highlights)){
    hits <- sort(unique(unlist(apply(highlights, 1, function(coords){
      if(coords[1] %in% contigs){
        if(which(as.character(coords[1]) == contigs) > 1){
          coords[2] <- coords[2] + indexes[as.integer(coords[1]) - 1, 4]
          coords[3] <- coords[3] + indexes[as.integer(coords[1]) - 1, 4]
        }
        which(df.plot$chr == as.character(coords[1])
              & df.plot$pos <= as.numeric(coords[3])
              & df.plot$pos >= as.numeric(coords[2]))
      }
    }))))
    df.plot$color[hits] <- highlight.color
  }
  
  ymax <- -log10(df[which(!(is.na(df[, 3]))), 3])
  ymax <- max(max(ymax[which(!(is.infinite(ymax)))]), -log10(cutoff) + 2)
  
  par(mar=c(2.1, 3.1, 0.6, 1), bty="n")
  plot(x=c(0, max(indexes[, 4])), y=c(0, 1.1 * ymax), type="n", 
       yaxs="i", xaxt="n", yaxt="n", xlab="", 
       ylab="")
  abline(h=-log10(cutoff), col="red")
  points(df.plot[, 2], -log10(df.plot[, 3]), 
         cex=0.2, pch=19, 
         col=as.character(df.plot[, 4]))
  
  if(!is.null(highlights)){
    if(any(df.plot$color == highlight.color)){
      points(df.plot[which(df.plot$color == highlight.color), 2], 
             -log10(df.plot[which(df.plot$color == highlight.color), 3]), 
             cex=0.2, pch=19, 
             col=highlight.color)
    }
    legend("top", pch=19, col=highlight.color, cex=0.5, 
           legend=paste(highlight.name, " (n=",
                        prettyNum(nrow(highlights), big.mark=","),
                        ")", sep=""))
  }
  
  midpoints <- sapply(1:length(indexes[, 2]), function(i){
    return(mean(c(indexes[i, 4], indexes[i, 4] - indexes[i, 3])))
  })
  axis(1, at=c(0, max(indexes$sum)), tck=0, labels=NA)
  axis(1, at=midpoints, labels=indexes[, 1], 
       tick=F, line=-1.1, cex.axis=0.75)
  mtext(1, text="Chromosome", line=0.9)
  
  y.at <- seq(0, ceiling(par("usr")[4]), by=ceiling(par("usr")[4]/6))
  axis(2, at=y.at, labels=NA, tck=-0.02)
  axis(2, at=y.at, tick=F, line=-0.5, labels=y.at, 
       cex.axis=0.75, las=2)
  mtext(2, text=expression(-log[10](italic(p))), line=1.5)
}


# Quantile-quantile plot
qq <- function (stats, cutoff=NULL, highlights=NULL, 
                highlight.color="#4EE69A",
                highlight.name="Positive Controls",
                print.stats=T, legend=T){
  p <- as.numeric(stats$p)
  colors <- rep("grey15", times=length(p))
  if(!is.null(highlights)){
    hits <- sort(unique(unlist(apply(highlights, 1, function(coords){
        which(stats$chr == as.character(coords[1])
              & stats$pos <= as.numeric(coords[3])
              & stats$pos >= as.numeric(coords[2]))
    }))))
    colors[hits] <- highlight.color
  }
  
  if (!is.numeric(p)){
    stop("P values must be numeric.")
  }
  keep.idxs <- which(!is.na(p) & !is.nan(p) & !is.null(p) & 
                       is.finite(p) & p < 1 & p > 0)
  p <- p[keep.idxs]
  colors <- colors[keep.idxs]
  
  new.order <- order(p)
  p <- p[new.order]
  colors <- colors[new.order]
  
  expected <- ppoints(length(p))
  qqconf <- function (p.expected){
    n <- length(p.expected)
    mpts <- matrix(nrow=n * 2, ncol=2)
    for (i in seq(from=1, to=n)) {
      mpts[i, 1] <- -log10((i - 0.5)/n)
      mpts[i, 2] <- -log10(qbeta(0.975, i, n - i))
      mpts[n * 2 + 1 - i, 1] <- -log10((i - 0.5)/n)
      mpts[n * 2 + 1 - i, 2] <- -log10(qbeta(0.025, i, n - 
                                               i))
    }
    mpts <- as.data.frame(mpts)
    return(mpts)
  }
  conf.int <- qqconf(expected)
  
  ks.p <- suppressWarnings(ks.test(p, "punif")$p.value)
  lambda <- dchisq(median(p), df=1)/dchisq(median(expected), df=1)
  
  p <- -log10(p)
  expected <- -log10(expected)
  if (is.null(cutoff)){
    cutoff <- nominal/length(p)
  }
  
  ymax <- max(max(p[which(!(is.infinite(p)))]), -log10(cutoff) + 2)
  
  par(mar=c(2.1, 3.1, 0.6, 1), bty="n")
  plot(x=expected, y=p, type="n", xaxt="n", yaxt="n", 
       xlab="", ylab="", xaxs="i", yaxs="i",
       xlim=c(0, 1.1 * max(expected)), ylim=c(0, 1.1 * ymax))
  polygon(x=conf.int[, 1], y=conf.int[, 2], col="gray90", border=NA)
  abline(0, 1, col="gray50")
  abline(h=-log10(cutoff), col="red")
  
  points(x=expected, y=p, pch=19, col=colors, cex=0.2)
  
  if(!is.null(highlights)){
    if(any(colors == highlight.color)){
      points(x=expected[which(colors == highlight.color)],
             y=p[which(colors == highlight.color)], 
             col=colors[which(colors == highlight.color)], 
             pch=19, cex=0.2)
    }
    if(legend==T){
      legend("left", pch=19, col=highlight.color, cex=0.5, 
             legend=paste(highlight.name, " (n=",
                          prettyNum(nrow(highlights), big.mark=","),
                          ")", sep=""))
    }
  }
  
  axis(1, at=axTicks(1), labels=NA, tck=-0.02)
  axis(1, at=axTicks(1), cex.axis=0.75, tick=F, line=-0.9)
  mtext(1, text=expression(Expected ~ ~-log[10](italic(p))), line=1.2)
  
  axis(2, at=axTicks(2), labels=NA, tck=-0.02)
  axis(2, at=axTicks(2), cex.axis=0.75, tick=F, line=-0.6, las=2)
  mtext(2, text=expression(Observed ~ ~-log[10](italic(p))), line=1.5)
  
  if (print.stats == T){
    # text(par("usr")[1], 0.85 * par("usr")[4], pos=4, 
    #      labels=(paste("K-S p=", formatC(ks.p, format="E", digits=2), "\n", sep="")))
    text(par("usr")[1], 0.9 * par("usr")[4], pos=4, font=2, 
         labels=bquote(lambda ~ .(paste("=", sprintf("%.2f", lambda), sep=""))))
  }
}



# Dev parameters:
stats.in <- "~/scratch/NDD_burden_tmp/meta1.HP0012759.rCNV.DEL.sliding_window.stats.bed.gz"
p.col.name <- "fisher_phred_p"
p.is.phred <- T
cutoff <- 10^-5
out.prefix <- "~/scratch/meta1.NDD.test"
highlight.in <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz"
highlight.name <- "Positive Controls"


stats <- read.stats(stats.in, p.col.name, p.is.phred)

highlights <- read.highlight.bed(highlight.in)

png(paste(out.prefix, "manhattan.png", sep="."),
    height=1000, width=1800, res=400)
manhattan(stats, cutoff, highlights,
          highlight.name=highlight.name)
dev.off()

png(paste(out.prefix, "qq.png", sep="."),
    height=1000, width=1000, res=400)
qq(stats, cutoff, highlights,
   highlight.name=highlight.name,
   legend=F)
dev.off()

png(paste(out.prefix, "manhattan_with_qq.png", sep="."),
    height=1000, width=2800, res=400)
layout(matrix(c(1,2), nrow=1), widths=c(18, 10))
manhattan(stats, cutoff, highlights,
          highlight.name=highlight.name)
qq(stats, cutoff, highlights,
   highlight.name=highlight.name,
   legend=F)
dev.off()
