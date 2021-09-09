#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plotting functions used for large segment analyses


#' Generic segment scatterplot
#'
#' Generic segment scatterplot function
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param x numeric vector for X axis
#' @param y numeric vector for Y axis
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#' @param xlims X axis limits
#' @param ylims Y axis limits
#' @param add.lm add linear trendline \[default: TRUE\]
#' @param pt.cex expansion factor for points \[default: 1\]
#' @param blue.bg add light blue background \[default: TRUE\]
#' @param horiz.lines.at numeric vector indicating where to draw horizontal lines \[default: no lines\]
#' @param horiz.lines.lty lty parameter for `horiz.lines.at` \[default: 1\]
#' @param horiz.lines.color colors for `horiz.lines.at`
#' @param abline.a `a` value for embedded [abline()] call
#' @param abline.b `b` value for embedded [abline()] call
#' @param abline.lty lty parameter for embedded [abline()] call
#' @param xtitle title of X axis
#' @param x.title.line line of X axis title
#' @param x.at tick positions for X axis
#' @param x.labs labels for X axis ticks
#' @param x.labs.at if different from `x.at`, specify locations of `x.labs`
#' @param parse.x.labs format X axis labels with [parse()] \[default: FALSE\]
#' @param x.title.cex expansion factor for X axis title \[default: 1\]
#' @param ytitle title of Y axis
#' @param y.title.line line of Y axis title
#' @param y.at tick positions for Y axis
#' @param y.labs labels for Y axis ticks
#' @param y.labs.at if different from `y.at`, specify locations of `y.labs`
#' @param parse.y.labs format Y axis labels with [parse()] \[default: FALSE\]
#' @param y.title.cex expansion factor for Y axis title \[default: 1\]
#' @param parmar numeric vector of margins passed to [par()]
#'
#' @seealso [load.segment.table()]
#'
#' @export segs.scatter
#' @export
segs.scatter <- function(segs, x, y, subset_to_regions=NULL,
                         xlims=NULL, ylims=NULL, add.lm=T, pt.cex=1, blue.bg=TRUE,
                         horiz.lines.at=NULL, horiz.lines.lty=1, horiz.lines.color=NULL,
                         abline.a=NULL, abline.b=NULL, abline.lty=1,
                         xtitle=NULL, x.title.line=1.75, x.at=NULL, x.labs=NULL, x.labs.at=NULL, parse.x.labs=FALSE, x.title.cex=1,
                         ytitle=NULL, y.title.line=1.75, y.at=NULL, y.labs=NULL, y.labs.at=NULL, parse.y.labs=FALSE, y.title.cex=1,
                         parmar=c(3, 3, 0.8, 0.8)){
  # Get plot values
  if(!is.null(subset_to_regions)){
    keep.idx <- which(segs$region_id %in% subset_to_regions)
    segs <- segs[keep.idx, ]
    x <- x[keep.idx]
    y <- y[keep.idx]
  }
  if(is.null(xlims)){
    xlims <- range(x[which(!is.infinite(x))], na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- range(y[which(!is.infinite(y))], na.rm=T)
  }
  del.idx <- which(segs$cnv=="DEL")
  dup.idx <- which(segs$cnv=="DUP")
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
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims, xlab="", ylab="", xaxt="n", yaxt="n")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg, xpd=TRUE)

  # Add gridlines
  if(is.null(x.at)){
    x.at <- axTicks(1)
  }
  if(is.null(y.at)){
    y.at <- axTicks(2)
  }
  abline(v=x.at, h=y.at, col=grid.col)
  if(!is.null(horiz.lines.at)){
    n.horiz.lines <- length(horiz.lines.at)
    if(is.null(horiz.lines.color) | length(horiz.lines.color) != n.horiz.lines){
      horiz.lines.color <- rep(blueblack, n.horiz.lines)
    }
    sapply(1:length(horiz.lines.at), function(i){
      abline(h=horiz.lines.at[i], lty=horiz.lines.lty[i], col=horiz.lines.color[i])
    })
  }
  if(!is.null(abline.a)){
    sapply(1:length(abline.a), function(i){
      abline(a=abline.a[i], b=abline.b[i], lty=abline.lty[i], col=blueblack)
    })
  }

  # Add linear fits, if optioned
  if(add.lm==T){
    del.fit <- robust.lm(x=x[del.idx], y=y[del.idx], conf=0.95)
    dup.fit <- robust.lm(x=x[dup.idx], y=y[dup.idx], conf=0.95)
    polygon(x=c(del.fit$ci$x, rev(del.fit$ci$x)),
            y=c(del.fit$ci$lower, rev(del.fit$ci$upper)),
            border=NA, col=adjustcolor(cnv.colors[1], alpha=0.2))
    polygon(x=c(dup.fit$ci$x, rev(dup.fit$ci$x)),
            y=c(dup.fit$ci$lower, rev(dup.fit$ci$upper)),
            border=NA, col=adjustcolor(cnv.colors[2], alpha=0.2))
    abline(del.fit$fit, lwd=2, col=cnv.colors[1])
    abline(dup.fit$fit, lwd=2, col=cnv.colors[2])
  }

  # Add points (always add gw-sig last)
  gw.idx <- which(segs$gw_sig)
  points(x[-gw.idx], y[-gw.idx], pch=segs$pt.pch[-gw.idx],
         bg=segs$pt.bg[-gw.idx], col=segs$pt.border[-gw.idx],
         cex=pt.cex)
  points(x[gw.idx], y[gw.idx], pch=segs$pt.pch[gw.idx],
         bg=segs$pt.bg[gw.idx], col=segs$pt.border[gw.idx],
         cex=pt.cex)

  # Add axis ticks
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=x.at, labels=NA, tck=-0.03, col=blueblack)
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.at, labels=NA, tck=-0.03, col=blueblack)

  # Add axis labels
  if(is.null(x.labs)){
    x.labs <- x.at
  }
  if(is.null(x.labs.at)){
    x.labs.at <- x.at
  }
  if(is.null(y.labs)){
    y.labs <- y.at
  }
  if(is.null(y.labs.at)){
    y.labs.at <- y.at
  }
  sapply(1:length(x.labs.at), function(i){
    if(parse.x.labs==TRUE){
      axis(1, at=x.labs.at[i], labels=parse(text=x.labs[i]), tick=F, line=-0.6)
    }else{
      axis(1, at=x.labs.at[i], labels=x.labs[i], tick=F, line=-0.6)
    }
  })
  sapply(1:length(y.labs.at), function(i){
    if(parse.y.labs==TRUE){
      axis(2, at=y.labs.at[i], labels=parse(text=y.labs[i]), tick=F, line=-0.6, las=2)
    }else{
      axis(2, at=y.labs.at[i], labels=y.labs[i], tick=F, line=-0.6, las=2)
    }
  })

  # Add axis titles
  mtext(1, text=xtitle, line=x.title.line, cex=x.title.cex)
  mtext(2, text=ytitle, line=y.title.line, cex=y.title.cex)
}


#' Segment swarm/boxplot
#'
#' Generic swarm/boxplot function for large segment analyses
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param x.bool boolean vector for subgrouping values on X axis
#' @param y numeric vector for Y axis
#' @param cnv.split split by CNV type \[default: TRUE\]
#' @param ylims Y axis limits
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#' @param add.pvalue print P-value on plot \[default: FALSE\]
#' @param stat.test statistical test to apply to groups \[default: `wilcoxon`\]
#' @param alternative alternative hypothesis for `stat.test` \[default: `two.sided`\]
#' @param print.fold print fold-enrichment to stdout \[default: FALSE\]
#' @param xtitle title of X axis
#' @param x.labs character vector for labels on X axis
#' @param add.y.axis add Y axis to plot \[default: TRUE\]
#' @param ytitle title of Y axis
#' @param y.title.line line of Y axis title
#' @param y.at tick positions for Y axis
#' @param y.labs labels for Y axis ticks
#' @param y.labs.at if different from `y.at`, specify locations of `y.labs`
#' @param parse.y.labs format Y axis labels with [parse()] \[default: FALSE\]
#' @param violin plot violins instead of boxes \[default: FALSE\]
#' @param pt.cex expansion factor for points\[default: 1\]
#' @param parmar numeric vector of margins passed to [par()]
#'
#' @export segs.swarm
#' @export
segs.swarm <- function(segs, x.bool, y, cnv.split=TRUE, ylims=NULL, subset_to_regions=NULL,
                       add.pvalue=FALSE, stat.test="wilcoxon", alternative="two.sided",
                       print.fold=FALSE, xtitle=NULL, x.labs=c("FALSE", "TRUE"),
                       add.y.axis=TRUE, ytitle=NULL, y.title.line=1.75,
                       y.at=NULL, y.labs=NULL, y.labs.at=NULL,
                       parse.y.labs=FALSE, violin=FALSE, pt.cex=1,
                       parmar=c(2.3, 3, 0.5, 0.5)){

  require(beeswarm, quietly=T)
  if(violin==T){
    require(vioplot, quietly=T)
  }

  # Subset to ID list, if optioned
  if(!is.null(subset_to_regions)){
    keep.idx <- which(segs$region_id %in% subset_to_regions)
    segs <- segs[keep.idx, ]
    x.bool <- x.bool[keep.idx]
    y <- y[keep.idx]
  }

  # Restrict to non-NA, finite values
  keep.idx <- which(!is.na(y) & !is.infinite(y))
  segs <- segs[keep.idx, ]
  x.bool <- x.bool[keep.idx]
  y <- y[keep.idx]

  # Get plot values
  if(is.null(ylims)){
    ylims <- range(y[which(!is.infinite(y))], na.rm=T)
  }
  if(cnv.split==TRUE){
    del.idx <- which(segs$cnv=="DEL")
    dup.idx <- which(segs$cnv=="DUP")
    x.at <- c(0.3, 0.7, 1.3, 1.7)
    width <- 0.2
    y.vals <- list(y[intersect(which(!x.bool), del.idx)],
                   y[intersect(which(!x.bool), dup.idx)],
                   y[intersect(which(x.bool), del.idx)],
                   y[intersect(which(x.bool), dup.idx)])
    if(print.fold==TRUE){
      del.fold <- mean(y.vals[[3]], na.rm=T)/mean(y.vals[[1]], na.rm=T)
      cat(paste("Deletions:", round(del.fold, 3), "fold\n"))
      dup.fold <- mean(y.vals[[4]], na.rm=T)/mean(y.vals[[2]], na.rm=T)
      cat(paste("Duplications:", round(dup.fold, 3), "fold\n"))
    }
    pt.color.list <- list(segs$pt.bg[intersect(which(!x.bool), del.idx)],
                          segs$pt.bg[intersect(which(!x.bool), dup.idx)],
                          segs$pt.bg[intersect(which(x.bool), del.idx)],
                          segs$pt.bg[intersect(which(x.bool), dup.idx)])
    pt.pch.list <- list(segs$pt.pch[intersect(which(!x.bool), del.idx)],
                        segs$pt.pch[intersect(which(!x.bool), dup.idx)],
                        segs$pt.pch[intersect(which(x.bool), del.idx)],
                        segs$pt.pch[intersect(which(x.bool), dup.idx)])
    pt.border.list <- list(segs$pt.border[intersect(which(!x.bool), del.idx)],
                           segs$pt.border[intersect(which(!x.bool), dup.idx)],
                           segs$pt.border[intersect(which(x.bool), del.idx)],
                           segs$pt.border[intersect(which(x.bool), dup.idx)])
    boxplot.colors <- rep(cnv.blacks[1:2], 2)
    boxplot.fill <- rep(cnv.whites[1:2], 2)
  }else{
    x.at <- c(0.5, 1.5)
    width <- 0.4
    y.vals <- list(y[which(!x.bool)],
                   y[which(x.bool)])
    if(print.fold==TRUE){
      cat(paste("True-vs-false:",
                round(mean(y.vals[[2]], na.rm=T)/mean(y.vals[[1]], na.rm=T), 3),
                "fold\n"))
    }
    pt.pch.list <- list(segs$pt.pch[which(!x.bool)],
                        segs$pt.pch[which(x.bool)])
    pt.color.list <- list(segs$pt.bg[which(!x.bool)],
                          segs$pt.bg[which(x.bool)])
    pt.border.list <- list(segs$pt.border[which(!x.bool)],
                           segs$pt.border[which(x.bool)])
    boxplot.colors <- rep(blueblack, 2)
    boxplot.fill <- rep(bluewhite, 2)
  }

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 2), ylim=ylims, type="n", xlab="", ylab="", xaxt="n", yaxt="n")

  # Add boxplots (or violins, if optioned)
  if(violin==T){
    lapply(1:length(y.vals), function(i){
      vioplot(y.vals[[i]], at=x.at[i], add=T, border=boxplot.colors[i], col=boxplot.fill[i],
              names=NA, drawRect=F, wex=2*width)
    })
    y.meds <- sapply(y.vals, median, na.rm=T)
    segments(x0=x.at-(2*width/3), x1=x.at+(2*width/3), y0=y.meds, y1=y.meds,
             lwd=2, col=boxplot.colors, lend="round")
  }else{
    boxplot(y.vals, at=x.at, outline=F, lty=1, add=T,
            outwex=width, staplewex=width, boxwex=width,
            border=boxplot.colors, col=boxplot.fill,
            xaxt="n", yaxt="n")
  }

  # Add swarms
  sapply(1:length(y.vals), function(i){
    beeswarm(y.vals[[i]], add=T, at=x.at[i],
             pwbg=pt.color.list[[i]], pwcol=pt.border.list[[i]],
             pwpch=pt.pch.list[[i]], cex=pt.cex,
             corral="random", corralWidth=width)
  })

  # Add x-axis
  sapply(1:2, function(x){
    axis(1, at=x-c(0.1, 0.9), tck=0, labels=NA, col=blueblack)
    axis(1, at=x-0.5, line=-0.9, labels=x.labs[x], tick=F)
  })
  if(!is.null(xtitle)){
    axis(1, at=c(0.1, 1.9), tck=0, labels=NA, line=1.2, col=blueblack)
    mtext(1, line=1.3, text=xtitle)
  }

  # Add y-axis, if optioned
  if(add.y.axis==T){
    if(is.null(y.at)){
      y.at <- axTicks(2)
    }
    axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
    axis(2, at=y.at, labels=NA, tck=-0.03, col=blueblack)
    if(is.null(y.labs)){
      y.labs <- y.at
    }
    if(is.null(y.labs.at)){
      y.labs.at <- y.at
    }
    sapply(1:length(y.labs.at), function(i){
      if(parse.y.labs==TRUE){
        axis(2, at=y.labs.at[i], labels=parse(text=y.labs[i]), tick=F, line=-0.6, las=2)
      }else{
        axis(2, at=y.labs.at[i], labels=y.labs[i], tick=F, line=-0.6, las=2)
      }
    })
    mtext(2, text=ytitle, line=y.title.line)
  }

  # Add P-values, if optioned
  if(add.pvalue==T){
    if(cnv.split==T){
      if(stat.test=="wilcoxon"){
        pvals <- sapply(c("DEL", "DUP"), function(cnv){
          stat.res <- wilcox.test(y[which(segs$cnv==cnv)] ~ x.bool[which(segs$cnv==cnv)], alternative=alternative)
          pval <- stat.res$p.value
          print(stat.res)
          return(pval)
        })
        axis(3, at=c(0.5, 1.5)-width, tck=0.03, col=blueblack, labels=NA, line=0.5)
        axis(3, at=mean(c(0.5, 1.5)-width), labels=format.pval(pvals[1]), line=-0.6, tick=F)
        axis(3, at=c(0.5, 1.5)+width, tck=0.03, col=blueblack, labels=NA, line=1.6)
        axis(3, at=mean(c(0.5, 1.5)+width), labels=format.pval(pvals[2]), line=0.5, tick=F)
      }
    }else{
      if(stat.test=="wilcoxon"){
        stat.res <- wilcox.test(y ~ x.bool, alternative=alternative)
        pval <- stat.res$p.value
        print(stat.res)
      }
      axis(3, at=c(0.5, 1.5), tck=0.03, col=blueblack, labels=NA, line=0.5)
      mtext(3, text=format.pval(pval), line=0.5)
    }
  }
}


#' Simple segment violin/swarmplot
#'
#' Simpler generic vioplot/swarmplot hybrid for plotting values split by DEL & DUP
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param y numeric vector for Y axis
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#' @param add.y.axis add Y axis to plot \[default: TRUE\]
#' @param ytitle title of Y axis
#' @param y.title.line line of Y axis title
#' @param add.pvalue print P-value on plot \[default: FALSE\]
#' @param stat.test statistical test to apply to groups \[default: `wilcoxon`\]
#' @param alternative alternative hypothesis for `stat.test` \[default: `two.sided`\]
#' @param pt.cex expansion factor for points\[default: 1\]
#' @param parmar numeric vector of margins passed to [par()]
#'
#' @seealso [segs.swarm()]
#'
#' @export segs.simple.vioswarm
#' @export
segs.simple.vioswarm <- function(segs, y, subset_to_regions=NULL,
                                 add.y.axis=T, ytitle=NULL, ytitle.line=1.75,
                                 add.pvalue=FALSE, stat.test="wilcoxon", alternative="two.sided",
                                 pt.cex=1, parmar=c(1.3, 3, 0.3, 0.3)){
  # Load necessary libraries
  require(vioplot, quietly=T)
  require(beeswarm, quietly=T)

  # Restrict to non-NA, finite values
  keep.idx <- which(!is.na(y) & !is.infinite(y))
  segs <- segs[keep.idx, ]
  y <- y[keep.idx]

  # Get plot data
  if(!is.null(subset_to_regions)){
    keepers <- which(segs$region_id %in% subset_to_regions)
    segs <- segs[keepers, ]
    y <- y[keepers]
  }
  ylims <- range(y, na.rm=T)
  plot.dat <- lapply(c("DEL", "DUP"), function(cnv){y[which(segs$cnv==cnv)]})
  pt.bg <- lapply(c("DEL", "DUP"), function(cnv){segs$pt.bg[which(segs$cnv==cnv)]})
  pt.color <- lapply(c("DEL", "DUP"), function(cnv){segs$pt.border[which(segs$cnv==cnv)]})
  pt.pch <- pt.colors <- lapply(c("DEL", "DUP"), function(cnv){segs$pt.pch[which(segs$cnv==cnv)]})

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 2), ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")

  # Add violins & swarms
  sapply(1:2, function(x){
    vioplot(plot.dat[[x]], at=x-0.5, add=T, wex=0.8, drawRect=F,
            col=cnv.whites[[x]], border=cnv.blacks[[x]])
    segments(x0=x-0.7, x1=x-0.3,
             y0=median(plot.dat[[x]], na.rm=T),
             y1=median(plot.dat[[x]], na.rm=T),
             lend="round", lwd=2, col=cnv.blacks[x])
    # boxplot(plot.dat[[x]], at=x-0.5, col=cnv.blacks[x], lty=1, outline=F,
    #         staplewex=0, add=T, boxwex=0.2)
    beeswarm(plot.dat[[x]], at=x-0.5, add=T, corralWidth=0.8, corral="wrap",
             pwcol=pt.color[[x]], pwbg=pt.bg[[x]], pwpch=pt.pch[[x]], cex=pt.cex)
  })

  # Add x-axis
  sapply(1:2, function(x){
    axis(1, at=x-c(0.1, 0.9), tck=0, labels=NA, col=blueblack)
    axis(1, at=x-0.5, line=-0.9, labels=c("DEL", "DUP")[x], tick=F)
  })

  # Add automatic Y-axis, if optioned
  if(add.y.axis==T){
    axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
    y.at <- axTicks(2)
    axis(2, at=y.at, labels=NA, tck=-0.03, col=blueblack)
    sapply(1:length(y.at), function(i){
      axis(2, at=y.at[i], labels=y.at[i], tick=F, line=-0.6, las=2)
    })
    mtext(2, text=ytitle, line=ytitle.line)
  }
  # Add P-values, if optioned
  if(add.pvalue==T){
    if(stat.test=="wilcoxon"){
      stat.res <- wilcox.test(plot.dat[[1]], plot.dat[[2]], alternative=alternative)
      pval <- stat.res$p.value
      print(stat.res)
    }
    axis(3, at=c(0.5, 1.5), tck=0.03, col=blueblack, labels=NA, line=0.5)
    mtext(3, text=format.pval(pval), line=0.5)
  }
}


#' Add a violin-histogram hybrid to a plot
#'
#' Helper function to add a single mirrored violin-histogram hybrid of values to an existing plot
#'
#' @param perm.dat.vals numeric vector of values
#' @param bins bins for histogram
#' @param y.at Y position for center of violin
#' @param width width of violin
#' @param obs.val observed value
#' @param obs.color color for observed value
#' @param obs.border border for observed value
#' @param obs.pch shape for observed value
#' @param color color for permuted distribution
#' @param border border for permuted distribution
#' @param diamond.cex expansion factor for observed value
#' @param y.title title for Y axis
#' @param left.ax.line add line to left axis \[default: FALSE\]
#'
#' @seealso [plot.seg.perms()]
#'
#' @export plot.viohist
#' @export
plot.viohist <- function(perm.dat.vals, bins, y.at, width=0.8,
                         obs.val=NA, obs.color=NA, obs.border=NA, obs.pch=23,
                         color=bluewhite, border=blueblack,
                         diamond.cex=4, y.title=NULL, left.ax.line=F){
  perm.hist <- hist(perm.dat.vals, breaks=bins, plot=F)
  perm.mean <- mean(as.numeric(perm.dat.vals[which(!is.infinite(perm.dat.vals))]), na.rm=T)
  values <- perm.hist$counts
  # Convert zero-bins to NAs for the outermost 5% of the distribution
  outer.lims <- quantile(perm.dat.vals, probs=c(0.025, 0.975), na.rm=T)
  outer.zero.bins <- intersect(which(perm.hist$mids < outer.lims[1] | perm.hist$mids > outer.lims[2]),
                               which(values==0))
  values[outer.zero.bins] <- NA
  values <- values / (max(values, na.rm=T) * 2/width)
  rect(xleft=bins[-length(bins)], xright=bins[-1],
       ybottom=-values + y.at, ytop=values + y.at,
       col=color, border="white")
  segments(x0=rep(bins[-length(bins)], 2),
           x1=rep(bins[-1], 2),
           y0=c(values + y.at, -values + y.at),
           y1=c(values + y.at, -values + y.at),
           col=border)
  segments(x0=rep(bins, 2), x1=rep(bins, 2),
           y0=c(y.at, values + y.at, y.at, -values + y.at),
           y1=c(values + y.at, y.at, -values + y.at, y.at),
           col=border)
  segments(x0=perm.mean, x1=perm.mean,
           y0=y.at - 0.2, y1=y.at + 0.2,
           col=blueblack, lwd=3, lend="round")
  segments(x0=obs.val, x1=obs.val,
           y0=y.at - 0.2, y1=y.at + 0.2,
           col=obs.color, lwd=3, lend="round")
  points(x=obs.val, y=y.at, pch=obs.pch, bg=obs.color,
         col=obs.border, cex=diamond.cex)
  if(!is.null(y.title)){
    axis(2, at=y.at, line=-0.8, tick=F, las=2, labels=y.title)
    if(left.ax.line==T){
      axis(2, at=c(y.at-(width/2), y.at+(width/2)), tck=0, labels=NA, col=blueblack)
    }
  }
}


#' Plot segment permutation results
#'
#' Function to plot segment permutation test results
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param perms list of permutation results
#' @param feature name of feature to evaluate
#' @param measure statistic to evaluate \[default: mean\]
#' @param norm normalize values before plotting \[default: FALSE\]
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#' @param n.bins ideal number of bins for histograms \[default: 100\]
#' @param min.bins minimum number of bins \[default: 10\]
#' @param xtitle title of X axis
#' @param x.title.line line of X axis title
#' @param xlims limits for X axis
#' @param xmin minimum for X axis
#' @param xmax maximum for X axis
#' @param diamond.pch shape for observed values
#' @param diamond.cex expansion for observed values
#' @param parmar numeric vector of margins passed to [par()]
#'
#' @seealso [plot.viohist()]
#'
#' @export plot.seg.perms
#' @export
plot.seg.perms <- function(segs, perms, feature, measure, norm=F,
                           subset_to_regions=NULL, n.bins=100, min.bins=10,
                           x.title=NULL, x.title.line=1.2, x.title.cex=1,
                           xlims=NULL, xmin=NULL, xmax=NULL,
                           diamond.pch=23, diamond.cex=1.25, parmar=c(2.25, 2, 0.5, 0.5)){
  # Get plot data
  if(!is.null(subset_to_regions)){
    segs <- segs[which(segs$region_id %in% subset_to_regions), ]
    perms <- lapply(perms, function(df){df[which(df$region_id %in% subset_to_regions), ]})
  }
  perm.dat <- perm.summary(perms, feature, measure)
  segs.dat <- calc.segs.dat(segs, feature, measure)

  # Normalize data, if optioned
  perm.dat.raw <- perm.dat
  segs.dat.raw <- segs.dat
  if(norm==T){
    perm.means <- apply(perm.dat, 2, mean, na.rm=T)
    perm.sds <- apply(perm.dat, 2, sd, na.rm=T)
    for(i in 1:3){
      perm.dat[, i] <- (perm.dat[, i] - perm.means[i]) / perm.sds[i]
      segs.dat[i] <- (segs.dat[i] - perm.means[i]) / perm.sds[i]
    }
  }

  # Determine value range and binning
  val.range <- range(perm.dat[which(!is.infinite(perm.dat) & !is.na(perm.dat) & !is.nan(perm.dat))], na.rm=T)
  val.range <- c(floor(val.range[1]), ceiling(val.range[2]))
  unique.vals <- length(unique(as.numeric(perm.dat)))
  if(unique.vals <= min.bins){
    bins <- seq(val.range[1], val.range[2], length.out=min.bins)
  }else if(unique.vals <= n.bins){
    bins <- seq(val.range[1], val.range[2], length.out=floor(unique.vals/2))
    # }else if(val.range[2] - val.range[1] <= n.bins){
    #   bins <- val.range[1]:val.range[2]
  }else{
    bins <- seq(val.range[1], val.range[2], length.out=n.bins)
  }
  bin.width <- bins[2]-bins[1]
  if(is.null(xlims)){
    xrange <- range(rbind(perm.dat, segs.dat), na.rm=T)
    xlims <- c(xrange[1], 1.25 * xrange[2])
  }
  if(!is.null(xmin)){
    xlims[1] <- xmin
  }
  if(!is.null(xmax)){
    xlims[2] <- xmax
  }

  # Gather more misc data for plotting
  perm.means <- apply(perm.dat, 2, mean, na.rm=T)
  perm.pvals <- sapply(1:3, function(i){calc.perm.p(perm.vals=perm.dat[, i], obs.val=segs.dat[i])})
  vio.colors <- rep(bluewhite, 3)
  vio.borders <- rep(blueblack, 3)
  # vio.colors <- c(purplewhite, redwhite, bluewhite)
  # vio.borders <- c(purpleblack, redblack, blueblack)
  row.colors <- c(cnv.colors[c(3, 1:2)])
  row.borders <- rep("black", 3)
  # row.borders <- c(purpleblack, redblack, blueblack)
  row.labels <- c("All\nCNV", "DEL", "DUP")
  stats.cex <- 0.85

  # Prep global plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=c(3, 0), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")

  # Add viohists
  sapply(1:3, function(i){
    plot.viohist(perm.dat[, i], bins, i-0.5,
                 color=vio.colors[i], border=vio.borders[i],
                 y.title=row.labels[i], diamond.cex=diamond.cex, obs.val=segs.dat[i],
                 obs.color=row.colors[i], obs.border=row.borders[i], obs.pch=diamond.pch)
    # segments(x0=perm.means[i], x1=perm.means[i],
    #          y0=i-0.7, y1=i-0.3, lwd=3,
    #          col=vio.borders[i], lend="round")
    if(segs.dat[i] >= perm.means[i]){
      text(x=segs.dat[i]-(0.03*(par("usr")[2]-par("usr")[1])), y=i-0.7, pos=4,
           labels=perm.pvals[2, ][[i]], xpd=T, cex=stats.cex)
    }else{
      text(x=perm.means[i]-(0.03*(par("usr")[2]-par("usr")[1])), y=i-0.7, pos=4,
           labels=perm.pvals[2, ][[i]], xpd=T, cex=stats.cex)
    }
    perm.mean <- mean(perm.dat.raw[, i], na.rm=T)
    fold <- round(segs.dat.raw[i]/perm.mean, 2)
    cat(paste(prettyNum(fold, small.interval=2), "-fold; obs=",
              prettyNum(segs.dat.raw[i], small.interval=4), "; exp=",
              prettyNum(perm.mean, small.interval=4), "\n", sep=""))
  })

  # Axes & cleanup
  axis(1, at=c(-10e10, 10e10), labels=NA, col=blueblack, tck=0)
  axis(1, at=unique(c(0, axTicks(1))), labels=NA, col=blueblack, tck=-0.03)
  sapply(1:length(axTicks(1)), function(i){
    axis(1, at=axTicks(1)[i], labels=prettyNum(axTicks(1)[i], big.mark=","), line=-0.65, tick=F)
  })
  mtext(1, text=x.title, line=x.title.line, cex=x.title.cex)
}


#' Multi-panel segment permutation plot
#'
#' Multi-panel plot of permutation results for segment subsets
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param gw.perms list of permutation results for empirical segments
#' @param lit.perms list of permutation results for literature segments
#' @param union.perms list of permutation results for all combined segments
#' @param feature name of feature to evaluate
#' @param measure statistic to evaluate \[default: mean\]
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#' @param norm normalize values before plotting \[default: FALSE\]
#' @param n.bins ideal number of bins for histograms \[default: 100\]
#' @param min.bins minimum number of bins \[default: 10\]
#' @param xtitle title of X axis
#' @param x.title.line line of X axis title
#' @param xlims limits for X axis
#' @param xmin minimum for X axis
#' @param xmax maximum for X axis
#' @param max.x.ticks maxmimum ticks for X axis
#' @param inner.axis.cex expansion of inner axis
#' @param diamond.cex expansion for observed values
#' @param parmar numeric vector of margins passed to [par()]
#'
#' @seealso [plot.seg.perms()]
#'
#' @export plot.seg.perms.multi
#' @export
plot.seg.perms.multi <- function(segs, gw.perms, lit.perms, union.perms,
                                 feature, measure, subset_to_regions=NULL,
                                 norm=F, n.bins=100, min.bins=10,
                                 x.title=NULL, xlims=NULL, xmin=NULL, xmax=NULL, max.x.ticks=5,
                                 outer.axis.cex=1, inner.axis.cex=0.85, diamond.cex=1.25,
                                 parmar=c(2.25, 6, 0.5, 0.5)){
  # Get ID subsets
  if(is.null(subset_to_regions)){
    subset_to_regions <- segs$region_id
  }
  all.sig.ids <- intersect(subset_to_regions, segs$region_id[which(segs$any_sig)])
  all.gw.ids <- intersect(subset_to_regions, segs$region_id[which(segs$gw_sig)])
  all.fdr.ids <- intersect(subset_to_regions, segs$region_id[which(segs$fdr_sig)])
  all.gd.ids <- intersect(subset_to_regions, segs$region_id[which(segs$any_gd)])
  gd.nonsig.ids <- intersect(subset_to_regions, setdiff(all.gd.ids, all.sig.ids))
  union.ids <- intersect(subset_to_regions, segs$region_id[which(segs$gw_sig | segs$any_gd)])

  # Get permutation plot data
  perm.dat <- list("union"=perm.summary(union.perms, feature=feature, measure=measure,
                                        subset_to_regions=union.ids),
                   "gw"=perm.summary(gw.perms, feature=feature, measure=measure,
                                     subset_to_regions=all.gw.ids),
                   "fdr"=perm.summary(gw.perms, feature=feature, measure=measure,
                                     subset_to_regions=all.fdr.ids),
                   "gd.nonsig"=perm.summary(lit.perms, feature=feature, measure=measure,
                                            subset_to_regions=gd.nonsig.ids))
  n.cats <- length(perm.dat)

  # Get observed plot data, and convert boolean valueas back to numeric, if needed
  segs.dat <- list("union"=calc.segs.dat(segs, feature, measure,
                                         subset_to_regions=union.ids),
                   "gw"=calc.segs.dat(segs, feature, measure,
                                      subset_to_regions=all.gw.ids),
                   "fdr"=calc.segs.dat(segs, feature, measure,
                                      subset_to_regions=all.fdr.ids),
                   "gd.nonsig"=calc.segs.dat(segs, feature, measure,
                                             subset_to_regions=gd.nonsig.ids))
  # Normalize data, if optioned
  perm.dat.raw <- perm.dat
  segs.dat.raw <- segs.dat
  if(norm==T){
    perm.means <- lapply(perm.dat, function(df){apply(df, 2, mean, na.rm=T)})
    perm.sds <- lapply(perm.dat, function(df){apply(df, 2, sd, na.rm=T)})
    for(i in 1:n.cats){
      for(j in 1:3){
        perm.dat[[i]][, j] <- (perm.dat[[i]][, j] - perm.means[[i]][j]) / perm.sds[[i]][j]
        segs.dat[[i]][j] <- (segs.dat[[i]][j] - perm.means[[i]][j]) / perm.sds[[i]][j]
      }
    }
  }

  # Determine value range and binning
  val.range <- range(unlist(perm.dat), na.rm=T)
  val.range <- c(floor(val.range[1]), ceiling(val.range[2]))
  unique.vals <- length(unique(as.numeric(unlist(perm.dat))))
  if(unique.vals <= min.bins){
    bins <- seq(val.range[1], val.range[2], length.out=min.bins)
  }else if(unique.vals <= n.bins){
    bins <- seq(val.range[1], val.range[2], length.out=floor(unique.vals/2))
    # }else if(val.range[2] - val.range[1] <= n.bins){
    #   bins <- val.range[1]:val.range[2]
  }else{
    bins <- seq(val.range[1], val.range[2], length.out=n.bins)
  }
  bin.width <- bins[2]-bins[1]
  if(is.null(xlims)){
    xrange <- range(rbind(do.call("rbind", perm.dat), do.call("rbind", segs.dat)), na.rm=T)
    xlims <- c(xrange[1], 1.25 * xrange[2])
  }
  if(!is.null(xmin)){
    xlims[1] <- xmin
  }
  if(!is.null(xmax)){
    xlims[2] <- xmax
  }

  # Gather more misc data for plotting
  perm.means <- lapply(perm.dat, function(df){apply(df, 2, mean, na.rm=T)})
  perm.pvals <- lapply(1:n.cats, function(i){
    sapply(1:3, function(j){calc.perm.p(perm.vals=perm.dat[[i]][, j], obs.val=segs.dat[[i]][j])})
  })
  vio.colors <- rep(bluewhite, n.cats)
  vio.borders <- rep(blueblack, n.cats)
  row.colors <- c(cnv.colors[c(3, 1:2)])
  row.borders <- rep("black", 3)
  # row.borders <- c(purpleblack, redblack, blueblack)
  outer.row.labels <- c("All\nCNV", "DEL", "DUP")
  inner.row.labels <- c("All Segs.", "GW Sig.", "FDR Sig.", "Lit. Only")
  stats.cex <- 0.85
  obs.pch <- c(23, 22, 23, 21)

  # Prep global plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=c(13, 0), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")

  # Add viohists
  sapply(1:n.cats, function(i){
    sapply(1:3, function(j){
      y.at <- (n.cats*j)-(n.cats-i)-0.5+(0.5*(j-1))
      plot.viohist(perm.dat[[i]][, j], bins, y.at,
                   color=vio.colors[j], border=vio.borders[j],
                   y.title=inner.row.labels[i], diamond.cex=diamond.cex, obs.val=segs.dat[[i]][j],
                   obs.pch=obs.pch[i], obs.color=row.colors[j], obs.border=row.borders[j],
                   left.ax.line=T)
      # segments(x0=perm.means[[i]][j], x1=perm.means[[i]][j],
      #          y0=y.at-0.2, y1=y.at+0.2, lwd=3,
      #          col=vio.borders[j], lend="round")
      if(segs.dat[[i]][j] >= perm.means[[i]][j]){
        text(x=segs.dat[[i]][j]-(0.015*(par("usr")[2]-par("usr")[1])), y=y.at-0.2, pos=4,
             labels=perm.pvals[[i]][2, ][[j]], xpd=T, cex=stats.cex)
      }else{
        text(x=perm.means[[i]][j]-(0.015*(par("usr")[2]-par("usr")[1])), y=y.at-0.2, pos=4,
             labels=perm.pvals[[i]][2, ][[j]], xpd=T, cex=stats.cex)
      }
    })
  })

  # Axes & cleanup
  x.at <- axTicks(1)
  if(length(x.at) > max.x.ticks){
    x.at <- x.at[seq(1, length(x.at), 2)]
  }
  axis(1, at=c(-10e10, 10e10), labels=NA, col=blueblack, tck=0)
  axis(1, at=x.at, labels=NA, col=blueblack, tck=-0.03)
  sapply(1:length(x.at), function(i){
    axis(1, at=x.at[i], labels=prettyNum(x.at[i], big.mark=","), line=-0.65, tick=F)
  })
  mtext(1, text=x.title, line=1.3)
  sapply(1:3, function(i){
    ax.at <- ((n.cats*i)+(0.5*(i-1)))-c(n.cats-0.1, 0.1)
    axis(2, at=ax.at, tck=0, labels=NA, col=blueblack, line=4)
    axis(2, at=mean(ax.at), tick=F, line=3.25, las=2, labels=outer.row.labels[i], cex.axis=outer.axis.cex)
  })
  segments(x0=rep(par("usr")[1], 2), x1=rep(par("usr")[2], 2),
           y0=c(4.25, 8.75), y1=c(4.25, 8.75),
           col=bluewhite, lend="round")
}


#' Plot all segment permutation results
#'
#' Wrapper function to plot four views on segment permutation results
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param gw.perms list of permutation results for empirical segments
#' @param lit.perms list of permutation results for literature segments
#' @param feature name of feature to evaluate
#' @param measure statistic to evaluate \[default: mean\]
#' @param outdir path to output directory
#' @param prefix filename prefix for output plots
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#' @param norm normalize values before plotting \[default: FALSE\]
#' @param norm.multi normalize values before plotting for multi-panel figure \[default: FALSE\]
#' @param n.bins.single ideal number of bins for single-panel histograms \[default: 100\]
#' @param n.bins.multi ideal number of bins for multi-panel histograms \[default: 100\]
#' @param min.bins minimum number of bins \[default: 10\]
#' @param xtitle title of X axis
#' @param x.title.line line of X axis title
#' @param xlims limits for X axis
#' @param xmin minimum for X axis
#' @param xmax maximum for X axis
#' @param diamond.cex expansion for observed values
#' @param pdf.dims.single dimensionf of pdf for single-panel plots
#' @param parmar.single numeric vector of margins passed to [par()] for single-panel plots
#' @param pdf.dims.multi dimensions for multi-panel plots
#' @param parmar.multi numeric vector of margins passed to [par()] for multi-panel plots

#'
#' @seealso [plot.seg.perms.multi()]
#'
#' @export plot.all.perm.res
#' @export
plot.all.perm.res <- function(segs, gw.perms, lit.perms,
                              feature, measure, outdir, prefix,
                              subset_to_regions=NULL, norm=F, norm.multi=F,
                              n.bins.single=100, n.bins.multi=100, min.bins=10,
                              x.title=NULL, xlims=NULL, xmin=NULL, xmax=NULL,
                              diamond.cex=1.25,
                              pdf.dims.single=c(2.2, 2.4),
                              parmar.single=c(2.25, 2, 0.5, 0.5),
                              pdf.dims.multi=c(2.2, 4.8),
                              parmar.multi=c(2.25, 6, 0.5, 0.5)){
  # Restrict segments to subset of regions, if optioned
  if(is.null(subset_to_regions)){
    subset_to_regions <- segs$region_id
  }
  segs <- segs[which(segs$region_id %in% subset_to_regions), ]

  # Get ID subsets
  all.sig.ids <- intersect(subset_to_regions, segs$region_id[which(segs$any_sig)])
  all.gw.ids <- intersect(subset_to_regions, segs$region_id[which(segs$gw_sig)])
  all.fdr.ids <- intersect(subset_to_regions, segs$region_id[which(segs$fdr_sig)])
  all.gd.ids <- intersect(subset_to_regions, segs$region_id[which(segs$any_gd)])
  gd.nonsig.ids <- intersect(subset_to_regions, setdiff(all.gd.ids, all.sig.ids))
  union.ids <- intersect(subset_to_regions, segs$region_id[which(segs$any_sig | segs$any_gd)])

  # Merge perm results
  union.perms <- lapply(1:length(gw.perms), function(i){
    shared.columns <- intersect(colnames(gw.perms[[i]]), colnames(lit.perms[[i]]))
    as.data.frame(rbind(gw.perms[[i]][shared.columns], lit.perms[[i]][shared.columns]))
  })

  # Prep output directory
  subdir <- paste(outdir, "/", prefix, "_", feature, "_", measure, sep="")
  if(!dir.exists(subdir)){
    dir.create(subdir)
  }

  # Plot union of all segments
  cat("Union of all segments:\n")
  pdf(paste(subdir, "/", prefix, ".", feature, ".", measure, ".union_all_segs.pdf", sep=""),
      height=pdf.dims.single[1], width=pdf.dims.single[2])
  plot.seg.perms(segs, union.perms, feature=feature, measure=measure,
                 subset_to_regions=union.ids,
                 n.bins=n.bins.single, min.bins=min.bins, norm=norm,
                 x.title=x.title, xlims=xlims, xmin=xmin, xmax=xmax,
                 parmar=parmar.single)
  dev.off()

  # Plot gw + FDR alone
  cat("Genome-wide + FDR significant:\n")
  pdf(paste(subdir, "/", prefix, ".", feature, ".", measure, ".gw_plus_fdr.pdf", sep=""),
      height=pdf.dims.single[1], width=pdf.dims.single[2])
  plot.seg.perms(segs, gw.perms, feature=feature, measure=measure,
                 subset_to_regions=all.sig.ids,
                 n.bins=n.bins.single, min.bins=min.bins, norm=norm,
                 x.title=x.title, xlims=xlims, xmin=xmin, xmax=xmax,
                 diamond.pch=23, parmar=parmar.single)
  dev.off()

  # Plot gw-sig alone
  cat("Genome-wide significant alone:\n")
  pdf(paste(subdir, "/", prefix, ".", feature, ".", measure, ".gw_sig.pdf", sep=""),
      height=pdf.dims.single[1], width=pdf.dims.single[2])
  plot.seg.perms(segs, gw.perms, feature=feature, measure=measure,
                 subset_to_regions=all.gw.ids,
                 n.bins=n.bins.single, min.bins=min.bins, norm=norm,
                 x.title=x.title, xlims=xlims, xmin=xmin, xmax=xmax,
                 diamond.pch=22, parmar=parmar.single)
  dev.off()

  # Plot FDR-sig alone
  cat("FDR significant alone:\n")
  pdf(paste(subdir, "/", prefix, ".", feature, ".", measure, ".fdr_sig.pdf", sep=""),
      height=pdf.dims.single[1], width=pdf.dims.single[2])
  plot.seg.perms(segs, gw.perms, feature=feature, measure=measure,
                 subset_to_regions=all.fdr.ids,
                 n.bins=n.bins.single, min.bins=min.bins, norm=norm,
                 x.title=x.title, xlims=xlims, xmin=xmin, xmax=xmax,
                 diamond.pch=23, parmar=parmar.single)
  dev.off()

  # Plot non-significant lit GDs alone
  cat("Literature GDs below genome-wide or FDR significance:\n")
  pdf(paste(subdir, "/", prefix, ".", feature, ".", measure, ".lit_gds_nogw_noFDR.pdf", sep=""),
      height=pdf.dims.single[1], width=pdf.dims.single[2])
  plot.seg.perms(segs, lit.perms, feature=feature, measure=measure,
                 subset_to_regions=gd.nonsig.ids,
                 n.bins=n.bins.single, min.bins=min.bins, norm=norm,
                 x.title=x.title, xlims=xlims, xmin=xmin, xmax=xmax,
                 diamond.pch=21, parmar=parmar.single)
  dev.off()

  # Plot combined analysis of all subsets
  pdf(paste(subdir, "/", prefix, ".", feature, ".", measure, ".multipanel.pdf", sep=""),
      height=pdf.dims.multi[1], width=pdf.dims.multi[2])
  plot.seg.perms.multi(segs, gw.perms, lit.perms, union.perms,
                       feature, measure, subset_to_regions=subset_to_regions,
                       n.bins=n.bins.multi, min.bins=min.bins, norm=norm.multi,
                       x.title=x.title, xlims=xlims, xmin=xmin, xmax=xmax,
                       diamond.cex=1, parmar=parmar.multi)
  dev.off()
}

