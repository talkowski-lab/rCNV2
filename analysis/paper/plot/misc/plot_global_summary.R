#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot summary of global rCNV effects across HPOs


options(stringsAsFactors=F, scipen=1000)
require(rCNV2, quietly=T)


######################
### DATA FUNCTIONS ###
######################
# Load table of global meta-analysis stats
load.global.stats <- function(stats.in){
  stats <- read.table(stats.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(stats)[1] <- gsub("#", "", colnames(stats)[1], fixed=T)
  return(stats)
}

# Collect effect sizes for a single panel
get.effect.sizes <- function(stats, category, hpos, return.log2=TRUE){
  res <- lapply(hpos, function(hpo){
    if(hpo %in% stats$HPO){
      df <- as.data.frame(t(sapply(c("DUP", "CNV", "DEL"), function(cnv){
        as.numeric(stats[which(stats$HPO==hpo & stats$category==category & stats$CNV==cnv),
                         c("meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper", "meta_neg_log10_p")])
      })))
    }else{
      df <- data.frame("meta_lnOR"=rep(NA, 3),
                       "meta_lnOR_lower"=rep(NA, 3),
                       "meta_lnOR_upper"=rep(NA, 3),
                       "meta_neg_log10_p"=rep(NA, 3),
                       row.names=c("DUP", "CNV", "DEL"))
    }
    df[, 1:3] <- apply(df[, 1:3], 1, exp)
    colnames(df) <- c("OR", "OR_lower", "OR_upper", "neg_log10_P")
    return(df)
  })
  if(return.log2){
    res <- lapply(res, function(df){
      df[, 1:3] <- apply(df[, 1:3], 1, log2)
      colnames(df)[1:3] <- paste("log2", colnames(df)[1:3], sep="_")
      return(df)
    })
  }
  names(res) <- hpos
  return(res)
}

# Get point formatting parameters for a single data frame generated by get.effect.sizes()
format.points <- function(df, bonf.p){
  bonf <- 10^-df[, "neg_log10_P"] <= bonf.p
  border <- control.cnv.colors[rownames(df)]
  border[which(bonf)] <- cnv.colors[rownames(df)][which(bonf)]
  fill <- cnv.whites[rownames(df)]
  fill[which(bonf)] <- cnv.colors[rownames(df)][which(bonf)]
  lwd <- rep(0.75, 3)
  lwd[which(bonf)] <- 1.25
  cex.adj <- rep(0.85, 3)
  cex.adj[which(bonf)] <- 1
  return(list("fill" = fill, "border" = border, "lwd" = lwd, "cex.adj" = cex.adj))
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot descriptive text columns
plot.text.columns <- function(hpos, meta, x.at=c(0, 2, 8), y.lab.buffer=0.05,
                              background=T, y.buffer=0.1, y.top=-0.6){
  # Get basic plotting info
  nhpos <- nrow(meta)
  y.at <- (1:nhpos)-0.5
  hpo.labels <- meta$HPO
  hpo.labels[which(hpo.labels == "UNKNOWN")] <- "N/A"

  # Prep plot area & add y-axis titles
  plot(x=NA, y=NA, xlim=c(0, 10), ylim=c(nhpos, y.top),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  sapply(1:length(x.at), function(i){
    axis(3, at=c(x.at[i]+0.125, (c(x.at[-1], par("usr")[2]) - 0.25)[i]),
         tck=0, labels=NA, xpd=T, col=blueblack, line=y.top+y.buffer)
  })
  height <- par("usr")[4]-par("usr")[3]
  text(x=x.at, y=par("usr")[4]+(y.lab.buffer*height)+(y.top*(sum(par("mar")[c(1,3)])/height)),
       labels=c("HPO", "Description", "Cases"),
       xpd=T, font=2, pos=4)

  # Plot text
  if(background==T){
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=y.at[which(!is.na(meta$HPO))]-0.5+y.buffer,
         ytop=y.at[which(!is.na(meta$HPO))]+0.5-y.buffer,
         border=NA, bty="n", col=bluewhite)
  }
  text(x=sapply(x.at, function(x){rep(x, nhpos)}), y=rep(y.at, 3),
       labels=c(hpo.labels, meta$description,
                sapply(meta$samples, function(x){if(is.na(x)){NA}else{prettyNum(x, big.mark=",")}})),
       pos=4, xpd=T)
}

# Plot a single vertical panel of effect sizes
plot.ors <- function(stats, category, hpos, xlims=NULL, title=NULL,
                     ytop=-0.6, y.buffer=0.1, y.lab.buffer=0.015,
                     pt.cex=0.85, pt.vex=0.15){
  # Get plot data
  pdat <- get.effect.sizes(stats, category, hpos)
  if(is.null(xlims)){
    xlims <- range(unlist(sapply(pdat, function(x){x[, 1:3]})), na.rm=T)
  }
  xleft <- min(xlims); xright <- max(xlims)
  na.idxs <- which(sapply(hpos, is.na))

  # Prep plot area
  par(bty="n")
  plot(NA, xlim=xlims, ylim=c(length(hpos), ytop),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  rect(xleft=xleft, xright=xright, ybottom=length(pdat), ytop=0,
       border=NA, bty="n", col=bluewhite)
  rect(xleft=xleft, xright=xright,
       ybottom=c(0:(length(pdat)-1))+y.buffer,
       ytop=c(1:length(pdat))-y.buffer,
       border=NA, bty="n", col="white")
  segments(x0=0, x1=0, y0=0, y1=length(hpos), col=ns.color)
  if(length(na.idxs) > 0){
    rect(xleft=xleft, xright=xright, ybottom=na.idxs-1, ytop=na.idxs,
         border="white", bty="o", col="white")
  }
  rect(xleft=xleft, xright=xright,
       ybottom=c(0, na.idxs), ytop=c(na.idxs-1, length(hpos)),
       border=blueblack, col=NA, xpd=T)
  height <- par("usr")[4]-par("usr")[3]
  text(x=mean(par("usr")[1:2]), y=par("usr")[4]+(y.lab.buffer*height),
       labels=title, xpd=T, pos=3)

  # Add points
  point.ymods <- -0.5+c(-pt.vex, 0, pt.vex)
  sapply(1:length(pdat), function(y){
    df <- pdat[[y]]
    fmt <- format.points(df, bonf.p=0.05/nrow(stats))
    segments(x0=df[, 2], x1=df[, 3], y0=y+point.ymods, y1=y+point.ymods,
             lwd=fmt$lwd, col=fmt$border, lend="round")
    points(x=df[, 1], y=y+point.ymods, pch=21, cex=pt.cex*fmt$cex.adj,
           col=fmt$border, lwd=fmt$lwd, bg=fmt$fill)
  })

  # Add top Y-axis
  axis(3, at=-10:10, col=blueblack, tck=-0.025, labels=NA, line=ytop)
  sapply(-10:10, function(x){
    axis(3, at=x, tick=F, labels=2^x, line=ytop-0.5)
  })
}

# Plot a single horizontal panel of effect sizes
plot.ors.horiz <- function(stats, category, hpos, ylims=NULL, title=NULL,
                           x.buffer=0.1, pt.cex=1, pt.wex=0.15,
                           channels.above=FALSE, invert.channels.above=FALSE,
                           channels.below=FALSE, invert.channels.below=FALSE,
                           channel.step.above=NULL, channel.step.below=NULL,
                           shade.channels=FALSE, channel.step=5,
                           hline=NA, hline.color=blueblack){
  # Get plot data
  pdat <- get.effect.sizes(stats, category, hpos)
  if(is.null(ylims)){
    ylims <- range(unlist(sapply(pdat, function(x){x[, 1:3]})), na.rm=T)
  }
  ybottom <- min(ylims); ytop <- max(ylims)
  na.idxs <- which(sapply(hpos, is.na))
  channel.idxs <- seq(0, length(hpos)+channel.step, by=channel.step)
  if(!is.null(channel.step.above)){
    top.channel.idxs <- seq(0, length(hpos)+channel.step.above, by=channel.step.above)
  }else{
    top.channel.idxs <- channel.idxs
  }
  if(!is.null(channel.step.below)){
    bottom.channel.idxs <- seq(0, length(hpos)+channel.step.below, by=channel.step.below)
  }else{
    bottom.channel.idxs <- channel.idxs
  }
  hpo.colors.by.severity <- adjustcolor(sapply(hpos, get.hpo.color, color.by="severity"),
                                        alpha=0.3)

  # Prep plot area
  par(bty="n")
  plot(NA, xlim=c(0, length(hpos)), ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  y.mid <- mean(par("usr")[3:4])
  if(channels.above){
    if(invert.channels.above){
      top.channel.col <- "white"
      top.channel.border <- hpo.colors.by.severity
      rect(xleft=0, xright=length(hpos), ybottom=y.mid, ytop=10e10, xpd=T,
           col=bluewhite, border=NA, bty="n")
    }else{
      top.channel.col <- hpo.colors.by.severity
      top.channel.border <- "white"
    }
    rect(xleft=top.channel.idxs[-length(top.channel.idxs)]+x.buffer,
         xright=top.channel.idxs[-1]-x.buffer,
         ybottom=y.mid, ytop=10e10, xpd=T,
         col=top.channel.col, border=top.channel.border, bty="o")
  }
  if(channels.below){
    if(invert.channels.below){
      bottom.channel.col <- "white"
      bottom.channel.border <- hpo.colors.by.severity
      rect(xleft=0, xright=length(hpos), ybottom=y.mid, ytop=-10e10, xpd=T,
           col=bluewhite, border=NA, bty="n")
    }else{
      bottom.channel.col <- hpo.colors.by.severity
      bottom.channel.border <- "white"
    }
    rect(xleft=bottom.channel.idxs[-length(bottom.channel.idxs)]+x.buffer,
         xright=bottom.channel.idxs[-1]-x.buffer,
         ybottom=y.mid, ytop=-10e10, xpd=T,
         col=bottom.channel.col, border=bottom.channel.border, bty="n")
  }
  rect(xleft=0, xright=length(pdat), ybottom=ybottom, ytop=ytop,
       border=NA, bty="n", col="white")
  if(shade.channels){
    rect(xleft=0, xright=length(pdat), ybottom=ybottom, ytop=ytop,
         border=NA, bty="n", col=bluewhite)
    rect(xleft=channel.idxs[-length(channel.idxs)]+x.buffer,
         xright=channel.idxs[-1]-x.buffer,
         ybottom=ybottom, ytop=ytop,
         border=NA, bty="n", col="white")
  }else{
    axis(1, at=0:length(pdat), tck=0.02, col=bluewhite, labels=NA)
    axis(3, at=0:length(pdat), tck=0.02, col=bluewhite, labels=NA)
  }
  segments(x0=0, x1=length(hpos), y0=0, y1=0, col=ns.color)
  if(!is.na(hline)){
    abline(h=hline, col=hline.color, lty=5, lwd=0.75)
  }
  if(length(na.idxs) > 0){
    rect(xleft=na.idxs-1, xright=na.idxs, ybottom=ybottom, ytop=ytop,
         border="white", bty="o", col="white")
  }
  rect(xleft=c(0, na.idxs), xright=c(na.idxs-1, length(hpos)),
       ybottom=ybottom, ytop=ytop,
       border=blueblack, col=NA, xpd=T)
  mtext(3, text=title, font=1)

  # Add points
  point.xmods <- -0.5+c(-pt.wex, 0, pt.wex)
  sapply(1:length(pdat), function(x){
    df <- pdat[[x]]
    fmt <- format.points(df, bonf.p=0.05/nrow(stats))
    segments(x0=x+point.xmods, x1=x+point.xmods, y0=df[, 2], y1=df[, 3],
             lwd=fmt$lwd, col=fmt$border, lend="round")
    points(x=x+point.xmods, y=df[, 1], pch=21, cex=pt.cex*fmt$cex.adj,
           col=fmt$border, lwd=pt.cex*fmt$lwd, bg=fmt$fill)
  })

  # Add left Y-axis
  axis(2, at=-10:10, col=blueblack, tck=-0.025, labels=NA, line=0)
  sapply(-10:10, function(y){
    axis(2, at=y, tick=F, labels=2^y, line=-0.5, las=2)
  })
  mtext(2, text="Odds Ratio", line=2.1)
}

# Plot diagonal text labels
plot.text.diag <- function(hpos, meta, y.at=4, background=T,
                           x.buffer=0.1, box.vex=0.5, angle=40,
                           angle.text.mod=12, text.xadj=-0.8){
  # Get basic plotting info
  nhpos <- length(hpos)
  x.at <- (1:nhpos)-0.5
  hpo.labels <- hpos
  hpo.labels[which(hpo.labels == "UNKNOWN")] <- "N/A"
  diag.xmod <- 1/tan(angle*pi/180)

  # Prep plot area & add y-axis titles
  par(bty="n")
  plot(x=NA, y=NA, xlim=c(0, nhpos), ylim=c(0, 11.2),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")

  # Get shading colors
  hpo.shading.colors <- sapply(hpos, get.hpo.color, color.by="severity")

  # First row of square indicators
  rect(xleft=(1:nhpos)-1+x.buffer, xright=(1:nhpos)-x.buffer,
       ybottom=10, ytop=par("usr")[4], border=NA, bty="n",
       col=adjustcolor(hpo.shading.colors, alpha=0.3))
  hpo.box.colors <- sapply(hpos, get.hpo.color, color.by="neuro")
  rect(xleft=(1:nhpos)-1+x.buffer, xright=(1:nhpos)-x.buffer,
       ybottom=10.5+c(-box.vex/2), ytop=10.5-c(-box.vex/2),
       xpd=T, border=blueblack, col=hpo.box.colors)

  # Draw diagonal shading
  sapply(1:nhpos, function(x){
    polygon(x=c(x-1+x.buffer, x-x.buffer, x-x.buffer-(10*diag.xmod),
                x-1+x.buffer-(10*diag.xmod)),
            y=c(10, 10, 0, 0), xpd=T, col=adjustcolor(hpo.shading.colors[x], alpha=0.3),
            border=NA, bty="n")
  })
  abline(h=y.at, col="white", xpd=T, lwd=1.5)
  col.header.y.at <- c(-0.15, y.at, 10)
  sapply(1:2, function(i){
    y0 <- col.header.y.at[i]+0.15
    y1 <- col.header.y.at[i+1]-0.15
    segments(x0=(diag.xmod*(-10+y0))-0.25, x1=(diag.xmod*(-10+y1))-0.25,
             y0=y0, y1=y1, xpd=T, col=blueblack)
  })

  # Add text labels
  text(x=(1:nhpos)-0.5-((10-y.at)*diag.xmod), y=rep(y.at-0.5, nhpos), pos=2,
       srt=angle + angle.text.mod, xpd=T,
       labels=prettyNum(meta$samples[match(hpos, meta$HPO)], big.mark=","))
  text(x=(1:nhpos)-(0.25*diag.xmod), y=rep(10-0.25, nhpos), pos=2,
       srt=angle + angle.text.mod, labels=hpo.abbrevs[hpos], xpd=T)

  # Add column titles
  text(x=c(-((10-y.at)*diag.xmod), 0)-0.5, y=c(y.at, 10),
       font=2, pos=c(2, 2), srt=angle + angle.text.mod, xpd=T, cex=1.25,
       labels=c("# Cases", "Phenotype"))
}


# Wrapper to plot full table of rCNV effect sizes organized by HPO hierarchy
plot.global.effects.byhpo <- function(stats, meta, panel.widths,
                                      y.top=-0.6, pt.cex=0.85,
                                      parmar=c(0.2, 0.1, 4, 0.2)){
  # Subset stats to HPOs in meta and to categories being plotted
  stats <- stats[which(stats$HPO %in% meta$HPO & stats$category %in% c(2, 3, 5, 8, 9)), ]
  hpos <- meta$HPO

  # Prep layout
  layout(matrix(c(1:length(panel.widths)), nrow=1, byrow=T), widths=panel.widths)
  par(bty="n", mar=parmar)

  # Plot dendrogram in left margin
  plot.dendro(meta, color=blueblack, lwd=2, box.wex=1.5, y.top=y.top)

  # Add text labels
  plot.text.columns(hpos, meta, x.at=c(0, 1.7, 8.75),
                    y.top=y.top, y.lab.buffer=0.0025)

  # Get effect size limits
  or.lims <- 1.1 * range(log2(exp(stats$meta_lnOR)))

  # Category 2: All Genic rCNVs
  plot.ors(stats, 2, meta$HPO, xlims=or.lims, ytop=y.top, pt.cex=pt.cex,
           title="Any Gene")

  # Category 3: Known GDs
  plot.ors(stats, 3, meta$HPO, xlims=or.lims, ytop=y.top, pt.cex=pt.cex,
           title="Established\nGenomic Disorders (GDs)")

  # Category 5: Known DS genes outside of GDs
  plot.ors(stats, 5, meta$HPO, xlims=or.lims, ytop=y.top, pt.cex=pt.cex,
           title="Dosage Sensitive (DS) Genes\n(Excluding GDs)")

  # Category 8: constrained genes genes outside of known GDs & DS genes
  plot.ors(stats, 8, meta$HPO, xlims=or.lims, ytop=y.top, pt.cex=pt.cex,
           title="PTV-Constrained Genes\n(Exclud. GDs + DS Genes)")

  # Category 9: all remaining genes
  plot.ors(stats, 9, meta$HPO, xlims=or.lims, ytop=y.top, pt.cex=pt.cex,
           title="All Other Genes\n(Excl. GDs, DS, and Constr.)")
}

# Wrapper to plot condensed figure of rCNV effect sizes sorted by effect size
plot.global.effects.byOR <- function(stats, meta, panel.heights,
                                     vertical.gridlines=FALSE, gridline.step=5,
                                     pt.cex=0.85, parmar=c(0, 3, 2, 0.2)){
  # Subset stats to HPOs in meta and only categories being plotted
  stats <- stats[which(stats$HPO %in% meta$HPO & stats$category %in% c(3, 7)), ]
  hpos <- meta$HPO[which(!is.na(meta$HPO))]

  # Prep layout
  layout(matrix(c(1:length(panel.heights)), ncol=1, byrow=T), heights=panel.heights)
  par(bty="n", mar=parmar)

  # Get effect size range of categories 3 and 7
  cat3.ors <- do.call("rbind", lapply(get.effect.sizes(stats, 3, hpos), function(df){df[, 1]}))
  cat7.ors <- do.call("rbind", lapply(get.effect.sizes(stats, 7, hpos), function(df){df[, 1]}))
  or.lims <- 1.1 * range(cbind(cat3.ors, cat7.ors), na.rm=T)

  # Sort HPOs by effect size of constrained gene deletions
  hpos.ordered <- hpos[order(-cat7.ors[, 3])]
  last.dev.idx <- max(which(hpos.ordered %in% developmental.hpos))

  # Category 3: Known GDs
  plot.ors.horiz(stats, 3, hpos.ordered, pt.cex=pt.cex, ylims=or.lims,
                 channels.below=FALSE, invert.channels.below=TRUE,
                 shade.channels=vertical.gridlines, channel.step=gridline.step)

  # Category 7: Constrained genes outside of known GDs
  plot.ors.horiz(stats, 7, hpos.ordered, pt.cex=pt.cex, ylims=or.lims, hline=log2(2),
                 channels.above=FALSE, channels.below=TRUE,
                 invert.channels.above=TRUE, channel.step.below=1,
                 shade.channels=vertical.gridlines, channel.step=gridline.step)
  abline(v=last.dev.idx+c(-0.035, 0.035), lwd=1, col=severity.colors[2:3])

  # Add HPO information in bottom half of plot
  parmar[3] <- 0
  par(mar=parmar)
  plot.text.diag(hpos.ordered, meta, y.at=3.25, box.vex=0.8,
                 angle=48, angle.text.mod=-3.5, text.xadj=-0.82)
  segments(x0=c(par("usr")[1], last.dev.idx),
           x1=c(last.dev.idx, par("usr")[2]),
           y0=par("usr")[4], y1=par("usr")[4],
           lwd=4, col=severity.colors[2:3], xpd=T, lend="butt")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog hpos.txt metadata.tsv counts.tsv stats.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: hpos.txt, metadata.tsv, counts.tsv, stats.tsv, out_prefix\n", sep=" "))
}

# Writes args & opts to vars
hpos.in <- args$args[1]
meta.in <- args$args[2]
sample.counts.in <- args$args[3]
stats.in <- args$args[4]
out.prefix <- args$args[5]

# # DEV PARAMETERS
# hpos.in <- "~/scratch/rCNV2_analysis_d2.reordered_hpos.txt"
# meta.in <- "~/scratch/phenotype_groups.HPO_metadata.txt"
# sample.counts.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# stats.in <- "~/scratch/rCNV2_analysis_d2.global_burden_stats.tsv.gz"
# out.prefix <- "~/scratch/test"

# Read data
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
meta <- load.hpo.metadata(meta.in, sample.counts.in, hpos)
stats <- load.global.stats(stats.in)

# Plot full figure ordered by HPOs
panel.widths <- c(1.1, 10, rep(3.5, 5))
pdf(paste(out.prefix, "ordered_by_hpo.pdf", sep="."),
    height=7.75, width=12)
plot.global.effects.byhpo(stats, meta, panel.widths,
                          parmar=c(0.2, 0.1, 3, 0.2))
dev.off()

# Plot slimmer figure ordered by average effect size
panel.heights <- c(4, 4, 3.5)
pdf(paste(out.prefix, "ordered_by_effect_size.pdf", sep="."),
    height=5, width=9.8)
plot.global.effects.byOR(stats, meta, panel.heights,
                         vertical.gridlines=T, gridline.step=5,
                         pt.cex=1.15, parmar=c(0.2, 12, 1.5, 0.2))
dev.off()

