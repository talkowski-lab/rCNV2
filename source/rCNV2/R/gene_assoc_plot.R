#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plotting functions used for gene association & fine-mapping analyses


#' Credible set scatterplot
#'
#' Generic credsets scatterplot function
#'
#' @param credsets data.frame of credible sets as loaded by [load.credsets()]
#' @param x values to plot on X axis
#' @param y values to plot on Y axis
#' @param subset_to_regions vector of credset IDs to include \[default: include all credsets\]
#' @param xlims X axis limits
#' @param ylims Y axis limits
#' @param add.lm add linear trendline \[default: TRUE\]
#' @param pt.cex expansion factor for points \[default: 1\]
#' @param horiz.lines.at numeric vector indicating where to draw horizontal lines \[default: no lines\]
#' @param horiz.lines.lty lty parameter for `horiz.lines.at` \[default: 1\]
#' @param horiz.lines.color colors for `horiz.lines.at`
#' @param abline.a `a` value for embedded [abline()] call
#' @param abline.b `b` value for embedded [abline()] call
#' @param abline.lty lty parameter for embedded [abline()] call
#' @param blue.bg add light blue background \[default: TRUE\]
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
#' @seealso [load.credsets()]
#'
#' @export credsets.scatter
#' @export
credsets.scatter <- function(credsets, x, y, subset_to_regions=NULL,
                             xlims=NULL, ylims=NULL, add.lm=T, pt.cex=1,
                             horiz.lines.at=NULL, horiz.lines.lty=1, horiz.lines.color=NULL,
                             abline.a=NULL, abline.b=NULL, abline.lty=1, blue.bg=TRUE,
                             xtitle=NULL, x.title.line=1.5, x.at=NULL, x.labs=NULL,
                             x.labs.at=NULL, parse.x.labs=FALSE, x.title.cex=1,
                             ytitle=NULL, y.title.line=1.7, y.at=NULL, y.labs=NULL,
                             y.labs.at=NULL, parse.y.labs=FALSE, y.title.cex=1,
                             parmar=c(2.7, 2.7, 0.25, 0.25)){
  # Get plot values
  if(!is.null(subset_to_regions)){
    keep.idx <- which(credsets$region_id %in% subset_to_regions)
    credsets <- credsets[keep.idx, ]
    x <- x[keep.idx]
    y <- y[keep.idx]
  }
  if(is.null(xlims)){
    xlims <- range(x[which(!is.infinite(x))], na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- range(y[which(!is.infinite(y))], na.rm=T)
  }
  del.idx <- which(credsets$cnv=="DEL")
  dup.idx <- which(credsets$cnv=="DUP")
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
       border=plot.border, bty=plot.bty, col=plot.bg)

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
  points(x, y, pch=credsets$pt.pch, cex=pt.cex,
         bg=credsets$pt.bg, col=credsets$pt.border)

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
