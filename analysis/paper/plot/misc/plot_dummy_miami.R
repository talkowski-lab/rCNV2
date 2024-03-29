#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot dummy Miami for graphical abstract of rCNV2 paper


options(stringsAsFactors=F, scipen=1000, family="sans")
require(rCNV2, quietly=T)


######################
### DATA FUNCTIONS ###
######################
# Load summary stats
load.sumstats <- function(stats.in, chrom.colors, sig.color, 
                          p.col.name="meta_neg_log10_p"){
  # Read data
  stats <- read.table(stats.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(stats)[1] <- gsub("#", "", colnames(stats)[1])
  
  # Get coordinates
  chr <- as.numeric(stats[, 1])
  if("pos" %in% colnames(stats)){
    pos <- stats$pos
  }else if(all(c("start", "end") %in% colnames(stats))){
    pos <- (stats$start + stats$end) / 2
  }else{
    stop("Unable to identify locus coordinate info. Must supply columns for either 'pos' or 'start' and 'end'.")
  }
  pos <- as.numeric(pos)
  
  # Get p-values
  if(p.col.name %in% colnames(stats)){
    p <- stats[, which(colnames(stats) == p.col.name)]
    p <- as.numeric(p)
  }else{
    stop(paste("Unable to identify p-value column by header name, ",
               p.col.name, ".", sep=""))
  }
  
  data.frame("chr"=chr, "pos"=pos, "p"=p)
}

# Get Manhattan plotting df
get.manhattan.plot.df <- function(df, chrom.colors, sig.color,
                                  max.p=12, gw.sig=6, sig.buffer=500000){
  contigs <- unique(df[, 1])
  contigs <- contigs[which(!(is.na(contigs)))]
  
  indexes <- as.data.frame(t(sapply(contigs, function(chr){
    return(c(chr, 0, max(df[which(df[, 1] == chr), 2])))
  })))
  indexes$sum <- cumsum(indexes[, 3])
  indexes$bg <- chrom.colors[(contigs %% 2) + 1]
  indexes[, 2:4] <- apply(indexes[, 2:4], 2, as.numeric)
  
  sig.idx <- which(df[, 3] >= gw.sig)
  sig.idx.all <- sort(unique(unlist(lapply(sig.idx, function(idx){
    which(df[, 1] == df[idx, 1] 
          & df[, 2] <= df[idx, 2] + sig.buffer 
          & df[, 2] >= df[idx, 2] - sig.buffer)
  }))))
  
  df.plot <- as.data.frame(t(sapply(1:nrow(df), function(idx){
    row <- df[idx, ]
    contig.idx <- which(indexes[, 1]==as.numeric(row[1]))
    pval <- as.numeric(row[3])
    if(is.na(pval)){
      pval <- 0
    }
    if(pval>max.p){
      pval <- max.p
    }
    if(idx %in% sig.idx.all){
      # pt.color <- sig.color
      # Disable significant peak highlighting
      pt.color <- indexes$bg[which(indexes[, 1] == contig.idx)]
      sig <- T
    }else{
      pt.color <- indexes$bg[which(indexes[, 1] == contig.idx)]
      sig <- T
    }
    return(c(row[1], 
             as.numeric(row[2]) + indexes[contig.idx, 4] - indexes[contig.idx, 3], 
             pval, 
             pt.color,
             sig))
  })))
  df.plot[, c(1, 4, 5)] <- apply(df.plot[, c(1, 4, 5)], 2, unlist)
  df.plot[, 2] <- as.numeric(as.character(df.plot[, 2]))
  df.plot[, 3] <- as.numeric(as.character(df.plot[, 3]))
  colnames(df.plot) <- c("chr", "pos", "p", "color", "sig")
  
  # Drop rows with NA p-values & return
  df.plot[which(!is.na(df.plot$p)), ]
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Custom Miami plot
mini.miami <- function(del, dup, max.p=10, 
                       gw.sig=6, sig.buffer=500000, sig.color=graphabs.green,
                       middle.axis.width=1,
                       parmar=c(0.5, 1.7, 0.5, 0.3)){
  # Get plotting data for dels & dups
  del.df <- get.manhattan.plot.df(del, 
                                  chrom.colors=c(cnv.colors[1], control.cnv.colors[1]),
                                  sig.color=sig.color,
                                  max.p=max.p, gw.sig=gw.sig, sig.buffer=sig.buffer)
  dup.df <- get.manhattan.plot.df(dup, 
                                  chrom.colors=c(cnv.colors[2], control.cnv.colors[2]),
                                  sig.color=sig.color,
                                  max.p=max.p, gw.sig=gw.sig, sig.buffer=sig.buffer)
  
  # Modify data to account for shared x-axis
  dup.df$p = dup.df$p + middle.axis.width/2
  del.df$p = -del.df$p - middle.axis.width/2
  
  # Get plot values
  x.range <- range(c(del.df$pos, dup.df$pos), na.rm=T)
  y.range <- range(c(del.df$p, dup.df$p), na.rm=T)
  
  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=x.range, ylim=y.range,
       xaxt="n", xlab="", yaxt="n", ylab="")
  abline(h=c(gw.sig + middle.axis.width/2, -gw.sig - middle.axis.width/2),
         lty=2, lwd=2, col=sig.color)
  
  # Add points
  points(x=del.df$pos, y=del.df$p, 
         col=del.df$color, pch=19, cex=0.275, xpd=T)
  points(x=dup.df$pos, y=dup.df$p, 
         col=dup.df$color, pch=19, cex=0.275, xpd=T)
  
  # Add axes
  y.at <- seq(0, ceiling(par("usr")[4]), by=ceiling(par("usr")[4]/6))
  axis(2, at=y.at+middle.axis.width/2, labels=NA, tck=0, col=blueblack, lwd=2)
  # axis(2, at=y.at+middle.axis.width/2, tick=F, line=-0.5, labels=abs(y.at), las=2)
  axis(2, at=-y.at-middle.axis.width/2, labels=NA, tck=0, col=blueblack, lwd=2)
  # axis(2, at=-y.at-middle.axis.width/2, tick=F, line=-0.5, labels=abs(y.at), las=2)
  axis(2, at=middle.axis.width/2 + (par("usr")[4]-par("usr")[3])/4,
       tick=F, line=-0.9, labels=bquote(-log[10](italic(P)) ~ "Duplication"),
       col.axis=blueblack)
  axis(2, at=-middle.axis.width/2 - (par("usr")[4]-par("usr")[3])/4,
       tick=F, line=-0.9, labels=bquote(-log[10](italic(P)) ~ "Deletion"),
       col.axis=blueblack)
  segments(x0=rep(par("usr")[1], 2), 
           x1=rep(par("usr")[2], 2),
           y0=c(0.5, -0.5) * middle.axis.width, 
           y1=c(0.5, -0.5) * middle.axis.width, 
           col=blueblack, lwd=2, xpd=T, lend="round")
  text(x=mean(par("usr")[1:2]), y=0, col=blueblack, labels="Chromosomes")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog del_stats.bed dup_stats.bed out.png", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: del_stats.bed, dup_stats.bed, and output.png\n", sep=" "))
}

# Writes args & opts to vars
del.in <- args$args[1]
dup.in <- args$args[2]
out.png <- args$args[3]

# Load sumstats
del <- load.sumstats(del.in, p.col.name="meta_neg_log10_p")
dup <- load.sumstats(dup.in, p.col.name="meta_neg_log10_p")

# Plot miami
png(out.png, height=4.1*300, width=3.75*300, res=300, family="sans", bg=NA)
mini.miami(del, dup, middle.axis.width=1.25, sig.col=graphabs.green, max.p=8)
dev.off()
