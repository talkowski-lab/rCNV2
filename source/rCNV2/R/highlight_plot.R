#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plotting functions for locus highlight figures


#' Add idiogram line
#'
#' Add stick figure for idiogram in locus highlights
#'
#' @param genome.in BEDTools-style .genome file
#' @param chrom chromosome to be plotted
#' @param y.at where the idiogram should be drawn, in Y coordinates
#' @param start left-most coordinate of highlight region \[default: no highlight\]
#' @param end right-most coordinate of highlight region \[default: no highlight\]
#' @param cytobands .bed file of cytobands \[default: no cytoband annotation\]
#' @param tick.height relative height of idiogram ticks \[default: -0.03\]
#' @param wex relative width expansion factor \[default: 0.925\]
#' @param y.label boolean indicator to add chromosome label to the Y axis \[default: TRUE\]
#'
#' @export plot.idio.stick
#' @export
plot.idio.stick <- function(genome.in, chrom, y.at, start=NA, end=NA, cytobands=NULL,
                            tick.height=-0.03, wex=0.925, y.label=TRUE){
  # Load & scale genomic coordinates for chromosome of interest
  g <- read.table(genome.in, sep="\t", header=F)
  len <- as.numeric(g[which(g[, 1] == chrom), 2])
  pbuf <- ((1-wex)/2) * diff(par("usr")[1:2])
  pstart <- par("usr")[1] + pbuf
  pend <- par("usr")[2] - pbuf

  # Scale highlight coordinates
  iwidth <- pend - pstart
  highlight <- (!is.na(start) & !is.na(end))
  if(highlight){
    h.start <- pstart + ((start / len) * iwidth)
    h.end <- pstart + ((end / len) * iwidth)
  }

  # Plot idiogram stick + highlight box
  segments(x0=pstart, x1=pend, y0=y.at, y1=y.at, col="gray50")
  if(highlight){
    rect(xleft=h.start, xright=h.end, ybottom=y.at-tick.height,
         ytop=y.at+tick.height, col="red")
  }
  if(y.label){
    axis(2, at=y.at, tick=F, line=-0.85, labels=paste("chr", chrom, sep=""), las=2)
  }

  # Add cytoband annotation, if optioned

  if(!is.null(cytobands)){
    cyto.lab <- get.cytoband.range(chrom, start, end, cytobands)
    text(x=mean(c(h.start, h.end)), y=y.at+(3*abs(tick.height)), pos=3,
         labels=cyto.lab, cex=5/6, xpd=T)
  }
}


#' Plot coordinate line
#'
#' Add a rounded coordinate line to a locus highlight plot
#'
#' @param start left-most coordinate
#' @param end right-most coordinate
#' @param y0 position on Y axis for coordinate line
#' @param highlight.start left-most coordinate(s) to highlight
#' @param highlight.end right-most coordinate(s) to highlight
#' @param highlight.col color for highlights \[default: "red"\]
#' @param tick.height relative height for vertical ticks \[default: 0.1\]
#' @param max.tick maximum number of ticks across the entire coordinate line \[default 6\]
#' @param vlines boolean indicator to add vertical lines for each tick \[default: FALSE\]
#' @param vlines.bottom lowest Y coordinate to extend `vlines` \[default: bottom of user coordinates\]
#' @param lab.cex character expansion scalar for labels \[default: 4.5/6\]
#'
#' @export plot.coord.line
#' @export
plot.coord.line <- function(start, end, y0, highlight.start, highlight.end, highlight.col="red",
                            tick.height=0.1, max.tick=6, vlines=FALSE, vlines.bottom=NULL, lab.cex=4.5/6){
  require(shape, quietly=TRUE)

  # Get plot coordinates
  log.base <- floor(log10((end-start)/6))
  tick.space <- ceiling((end-start)/(10^log.base))*(10^(log.base-1))
  tick.at <- seq(0, 10e8, by=tick.space)
  tick.at <- tick.at[which(tick.at >= start & tick.at <= end)]
  if(length(tick.at) > max.tick){
    tick.at <- tick.at[seq(1, length(tick.at), 2)]
    tick.space <- tick.at[2] - tick.at[1]
  }
  n.ticks <- length(tick.at)

  # Set tick labels
  if(log.base < 3){
    suffix <- "kb"
    denom <- 1000
    nsmall <- (3 - log.base) + 1
  }else{
    suffix <- "Mb"
    denom <- 1000000
    nsmall <- (6 - log.base) + 1
  }
  tick.labels <- paste(format(round(tick.at / denom, nsmall), nsmall=nsmall), suffix)
  tick.labels.y.at <- rep(y0, n.ticks)

  # Draw vertical lines for ticks, if optioned
  if(vlines){
    if(is.null(vlines.bottom)){
      vlines.bottom <- par("usr")[3]
    }
    segments(x0=tick.at, x1=tick.at, y0=y0, y1=vlines.bottom, col=bluewhite)
  }

  # Draw coordinate line & highlight box
  abline(h=y0, col=blueblack, lwd=2)
  rect(xleft=highlight.start, xright=highlight.end,
       ybottom=y0 - (0.5 * tick.height), ytop=y0 + (0.5 * tick.height),
       bty="n", border=NA, col=highlight.col)
  shape::Arrowhead(x0=par("usr")[1:2], y0=rep(y0, 2), angle=c(180, 0),
                   arr.length=0.3, arr.width=0.2, arr.col=blueblack, arr.type="triangle",
                   xpd=T, lcol=NA)
  segments(x0=tick.at, x1=tick.at,
           y0=rep(y0-tick.height, n.ticks), y1=rep(y0+tick.height, n.ticks),
           col=blueblack)

  # Add tick labels
  text(x=tick.at, y=tick.labels.y.at, cex=lab.cex, labels=tick.labels, pos=3, xpd=T)
}


#' Plot gene bodies
#'
#' Add gene bodies to locus highlight
#'
#' @param genes gene features to be plotted, loaded by [load.genes.from.gtf()]
#' @param y0 center position on Y axis for gene panel
#' @param transcripts boolean indicator to plot all transcripts per gene
#' \[default: collapse all transcripts per gene\]
#' @param mark.tss boolean indicator to mark TSS of each gene \[default: FALSE\]
#' @param n.rows number of rows of genes \[default: 1\]
#' @param panel.height relative height of gene panel \[default: 0.2\]
#' @param col color for genes \[default: ns.color\]
#' @param y.axis.title title for Y axis \[default: no title\]
#' @param y.axis.title.col color for `y.axis.title` \[default: no title\]
#' @param label.genes vector of gene symbols to be labeled on plot
#'
#' @seealso [load.genes.from.gtf]
#'
#' @export plot.gene.bodies
#' @export
plot.gene.bodies <- function(genes, y0, transcripts=FALSE, mark.tss=FALSE,
                             n.rows=1, panel.height=0.2, col=ns.color,
                             y.axis.title=NULL, y.axis.title.col=NULL, label.genes=c()){
  require(shape, quietly=T)
  if(transcripts){
    genes$gene <- genes$transcript
  }
  gene.names <- unique(genes$gene)
  n.genes <- length(gene.names)

  # Get y scaling
  panel.bottom <- y0 - (0.5*panel.height)
  panel.top <- y0 + (0.5*panel.height)
  row.height <- panel.height / n.rows
  exon.height <- (2/3) * row.height
  utr.height <- (1/3) * row.height
  if(mark.tss){
    exon.height <- 0.75 * exon.height
    utr.height <- 0.75 * utr.height
  }
  row.mids <- seq(panel.bottom + (0.5*row.height),
                  panel.top - (0.5*row.height),
                  length.out=n.rows)

  # Iterate over genes and plot each one at a time
  sapply(1:n.genes, function(i){
    gene <- gene.names[i]
    row.idx <- (i %% n.rows) + 1
    tx.coords <- genes[which(genes$gene==gene & genes$feature=="transcript"), 2:3]
    tx.coords$end[which(tx.coords$end > par("usr")[2])] <- par("usr")[2]
    tx.coords$start[which(tx.coords$start < par("usr")[1])] <- par("usr")[1]
    tx.strand <- genes$strand[which(genes$gene==gene & genes$feature=="transcript")]
    ex.coords <- genes[which(genes$gene==gene & genes$feature=="exon"), 2:3]
    utr.coords <- genes[which(genes$gene==gene & genes$feature=="UTR"), 2:3]
    if(nrow(tx.coords) > 0){
      if(mark.tss){
        tss.arrow.buffer <- 0.015 * (end-start)
        tss.arrow.height <- 0.5 * 0.75 * row.height
        tss.arrow.y1 <- row.mids[row.idx] + tss.arrow.height
        if(tx.strand=="+"){
          tss.arrow.start <- as.numeric(tx.coords[1])
          tss.arrow.end <- tss.arrow.start + tss.arrow.buffer
        }else{
          tss.arrow.start <- as.numeric(tx.coords[2])
          tss.arrow.end <- tss.arrow.start - tss.arrow.buffer
        }
        segments(x0=tss.arrow.start, x1=tss.arrow.start,
                 y0=row.mids[row.idx], y1=tss.arrow.y1,
                 col=blueblack)
        Arrows(x0=tss.arrow.start, x1=tss.arrow.end,
               y0=tss.arrow.y1, y1=tss.arrow.y1,
               arr.type="triangle", arr.length=0.12, arr.width=0.08,
               col=blueblack)
      }
      segments(x0=tx.coords[, 1], x1=tx.coords[, 2],
               y0=row.mids[row.idx], y1=row.mids[row.idx],
               lend="round", col=col)
    }
    if(nrow(ex.coords) > 0){
      rect(xleft=ex.coords[, 1], xright=ex.coords[, 2],
           ybottom=row.mids[row.idx]-(0.5*exon.height),
           ytop=row.mids[row.idx]+(0.5*exon.height),
           border=col, lwd=1, col=col)
    }
    if(nrow(utr.coords) > 0){
      rect(xleft=utr.coords[, 1], xright=utr.coords[, 2],
           ybottom=row.mids[row.idx]-(0.5*utr.height),
           ytop=row.mids[row.idx]+(0.5*utr.height),
           border=col, lwd=1, col=col)
    }
    if(length(label.genes) > 0){
      if(gene %in% label.genes){
        if(row.mids[row.idx] >= y0){
          glabel.pos <- 3
          glabel.y <- row.mids[row.idx]-(0.25*row.height)
        }else{
          glabel.pos <- 1
          glabel.y <- row.mids[row.idx]+(0.25*row.height)
        }
        glab.coords <- as.numeric(tx.coords[1, 1:2])
        glab.coords[1] <- max(c(par("usr")[1], glab.coords[1]))
        glab.coords[2] <- min(c(par("usr")[2], glab.coords[2]))
        text(x=mean(glab.coords), y=glabel.y, xpd=T,
             cex=5/6, font=3, labels=gene, pos=glabel.pos)
      }
    }
  })

  # Add y-axis title (if optioned)
  if(is.null(y.axis.title.col)){
    y.axis.title.col <- col
  }
  axis(2, at=y0, tick=F, line=-0.9, las=2, labels=y.axis.title, col.axis=y.axis.title.col)
}


#' Plot square bracket
#'
#' Add square bracket to locus highlight
#'
#' @param xleft left-most coordinate for bracket
#' @param xright right-most coordinate for bracket
#' @param y0 vertical midpoint for bracket
#' @param height height of bracket
#' @param col color for bracket \[default: blueblack\]
#' @param staple.wex relative width of bracket staples \[default: 0.025\]
#' @param lwd line width \[default: 1\]
#' @param left.label label to be printed at the left end of the bracket
#' \[default: no label]
#'
#' @export plot.bracket
#' @export
plot.bracket <- function(xleft, xright, y0, height, col=blueblack, staple.wex=0.025,
                         lwd=1, left.label=NULL){
  # Get various dimensions for plotting
  staple.width <- staple.wex*diff(par("usr")[1:2])
  half.height <- height / 2
  ytop <- y0 - half.height
  ybottom <- y0 + half.height

  # Plot staples
  segments(x0=rep(c(xleft, xright), 2), x1=rep(c(xleft-staple.width, xright+staple.width), 2),
           y0=c(ytop, ytop, ybottom, ybottom), y1=c(ytop, ytop, ybottom, ybottom),
           col=col, lwd=lwd, lend="round")

  # Plot vertical connectors
  segments(x0=c(xleft-staple.width, xright+staple.width), x1=c(xleft-staple.width, xright+staple.width),
           y0=rep(ybottom, 2), y1=rep(ytop, 2),
           col=col, lwd=lwd, lend="round")

  # Add labels, if optioned
  text(x=xleft-staple.width, y=y0, pos=2, col=col, labels=left.label)
}


#' Plot P-values for locus highlight
#'
#' Add panel of neglog10-scaled P-values to locus highlight
#'
#' @param ss summary statistics to be plotted as loaded by [load.sumstats.for.region]
#' @param y0 vertical midpoint for panel
#' @param cnv.type CNV type to be plotted
#' @param panel.height relative height of panel \[default: 0.6\]
#' @param pt.cex expansion scalar for points \[default: 0.6\]
#' @param gw.sig value to annotate as genome-wide significance
#' \[default: do not annotate genome-wide significance\]
#' @param min.y minimum -log10(P) for the maximum Y axis \[default: 9\]
#' @param gw.sig.label text label for `gw.sig` \[default "Genome-wide significance"\]
#' @param gw.sig.label.side orientation of `gw.sig.label` relative to `gw.sig` line
#' \[default: "above"\]
#'
#' @details `gw.sig.label.side` accepts either "above" or "below"
#'
#' @seealso [load.sumstats.for.region]
#'
#' @export plot.pvalues.for.highlight
#' @export
plot.pvalues.for.highlight <- function(ss, y0, cnv.type, panel.height=0.2,
                                       pt.cex=0.6, gw.sig=NULL, min.y=9,
                                       gw.sig.label="Genome-wide significance",
                                       gw.sig.label.side="above"){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height

  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    pt.col <- cnv.blacks[1]
  }else if(cnv.type=="DUP"){
    pt.col <- cnv.blacks[2]
  }

  # Infer plotting mode
  if("start" %in% colnames(ss)){
    pmode <- "rects"
  }else{
    pmode <- "points"
  }

  # Scale p-values according to y0 and panel.height
  if(pmode == "rects"){
    starts <- as.numeric(ss$start)
    ends <- as.numeric(ss$end)
  }else{
    pos <- as.numeric(ss$pos)
  }
  pvals.orig <- as.numeric(ss$meta_neg_log10_p)
  max.pval.orig <- max(c(9, (ceiling(max(pvals.orig, na.rm=T)) + 1)))
  pval.scale.factor <- (panel.height / max.pval.orig)
  pvals.scaled <- pval.scale.factor * pvals.orig
  pvals <- pvals.scaled + y0 - half.height

  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")

  # Add marker for genome-wide significance, if optioned
  if(!is.null(gw.sig)){
    gw.sig.y <- (-log10(gw.sig) * pval.scale.factor) + y0 - half.height
    abline(h=gw.sig.y, col=graphabs.green, lty=5)
    if(gw.sig.label.side=="above"){
      gw.sig.label.y <- gw.sig.y+((2/3)*diff(y.ax.tick.spacing[1:2]))
    }else{
      gw.sig.label.y <- gw.sig.y-((2/3)*diff(y.ax.tick.spacing[1:2]))
    }
    text(x=par("usr")[1], y=gw.sig.label.y, cex=5/6, font=3, col=graphabs.green,
         labels=gw.sig.label, pos=4)
  }

  # Add points
  if(pmode == "rects"){
    segments(x0=starts, x1=ends, y0=pvals, y1=pvals, lwd=4, lend="square", col=pt.col)
  }else{
    points(x=pos, y=pvals, pch=19, col=pt.col, cex=pt.cex)
  }

  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+c(-half.height, half.height), tick=F, las=2, line=-0.65,
       labels=c(0, max.pval.orig), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, las=2,
       labels=bquote(-log[10](italic("P")[.(cnv.type)])))

  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop,
           col=blueblack, xpd=T)
}


#' Plot odds ratios for locus highlight
#'
#' Add panel of log-odds ratios to locus highlight
#'
#' @param ss summary statistics to be plotted as loaded by [load.sumstats.for.region]
#' @param y0 vertical midpoint for panel
#' @param cnv.type CNV type to be plotted
#' @param max.distinct maximum number of distinct points to plot before
#' switching to a continuous line \[default: 50\]
#' @param dx with of CI shading for distinct points supplied as a
#' fraction of the X axis \[default: 0.01\]
#' @param panel.height relative height of panel \[default: 0.2\]
#' @param pt.cex expansion scalar for points \[default: 0.7\]
#'
#' @seealso [load.sumstats.for.region]
#'
#' @export plot.ors.for.highlight
#' @export
plot.ors.for.highlight <- function(ss, y0, cnv.type, max.distinct=50, dx=0.01,
                                   panel.height=0.2, pt.cex=0.8){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height

  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    line.col <- cnv.colors[1]
    ci.col <- control.cnv.colors[1]
  }else if(cnv.type=="DUP"){
    line.col <- cnv.colors[2]
    ci.col <- control.cnv.colors[2]
  }

  # Infer plotting mode
  if("start" %in% colnames(ss)){
    pmode <- "rects"
  }else{
    pmode <- "points"
  }

  # Scale odds ratios according to y0 and panel.height
  or.col.idxs <- grep("meta_lnOR", colnames(ss), fixed=T)
  ss[, or.col.idxs] <- apply(ss[, or.col.idxs], 2, as.numeric)
  ss <- na.omit(ss)
  if(pmode == "rects"){
    ss <- ss[order(ss$start), ]
    starts <- as.numeric(ss$start)
    ends <- as.numeric(ss$end)
  }else{
    ss <- ss[order(ss$pos), ]
    pos <- as.numeric(ss$pos)
  }
  ors.orig <- ss[, or.col.idxs]
  ors.orig <- log2(exp(ors.orig))
  ors.scalar <- max(ors.orig[, 1], na.rm=T)
  ors.scaled <- (panel.height / (ceiling(ors.scalar) + 1)) * ors.orig
  ors.scaled <- as.data.frame(t(apply(ors.scaled, 1, function(vals){
    sapply(vals, function(x){max(c(min(c(x, panel.height)), 0))})
  })))
  ors <- ors.scaled + y0 - half.height

  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")

  # Add points & shading
  if(pmode == "rects"){
    rect(xleft=starts, xright=ends, ybottom=ors[, 2], ytop=ors[, 3],
         col=adjustcolor(ci.col, alpha=0.5), border=NA, bty="n")
    segments(x0=starts, x1=ends, y0=ors[, 1], y1=ors[, 1],
             lwd=4, lend="square", col=line.col)
  }else{
    if(length(pos) <= max.distinct){
      xbuf <- 0.5 * dx * diff(par("usr")[1:2])
      rect(xleft=pos-xbuf, xright=pos+xbuf, ybottom=ors[, 2], ytop=ors[, 3],
           col=adjustcolor(ci.col, alpha=0.5), border=NA, bty="n")
      points(x=pos, y=ors[, 1], col=line.col, pch=15, cex=pt.cex)
    }else{
      polygon(x=c(pos, rev(pos)), y=c(ors[, 2], rev(ors[, 3])),
              col=adjustcolor(ci.col, alpha=0.5), border=NA, bty="n")
      points(x=pos, y=ors[, 1], lwd=3, col=line.col, type="l")
    }
  }

  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+c(-half.height, half.height), tick=F, las=2, line=-0.65,
       labels=c(0, ceiling(ors.scalar) + 1), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, labels=bquote("log" [2] * ("OR"[.(cnv.type)])), las=2)

  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)
}


#' Plot hatched N/A rectangle
#'
#' Add hatched NA rectangle to locus highlight
#'
#' @param xleft left-most coordinate(s) for rectangles
#' @param xright right-most coordinate(s) for rectangles
#' @param y0 vertical midpoint for panel
#' @param panel.height relative height of panel \[default: 0.2\]
#' @param text text label to superimpose on top of rectangles \[default: no text\]
#'
#' @export plot.na.rect
#' @export
plot.na.rect <- function(xleft, xright, y0, panel.height=0.2, text=NULL){
  half.height <- 0.5*panel.height
  rect(xleft=xleft, xright=xright, ybottom=y0 - half.height, ytop=y0 + half.height,
       border=blueblack, col=control.cnv.colors[2], density=15)
  text(x=mean(c(xleft, xright)), y=y0, labels=text, cex=5/6, font=3)
}


#' Plot constraint track for locus highlight
#'
#' Add plot of LOEUF values for all genes for locus highlight
#'
#' @param constr.df data.frame of constraint features per gene as loaded by [load.features()]
#' @param genes data.frame of gene features to plot as loaded by [load.genes.from.gtf()]
#' @param y0 vertical midpoint for panel
#' @param panel.height relative height of panel \[default: 0.2\]
#' @param label.genes vector of gene symbols to label \[default: don't label genes\]
#'
#' @seealso [load.features], [load.genes.from.gtf]
#'
#' @export plot.constraint.track
#' @export
plot.constraint.track <- function(constr.df, genes, y0, panel.height=0.2,
                                  label.genes=NULL){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height

  # Merge gene coordinates & pLI values
  gcoords <- genes[which(genes$gene %in% constr.df$gene & genes$feature=="transcript"),
                   c("start", "end", "gene")]
  colnames(constr.df)[which(colnames(constr.df) == "gnomad_oe_lof_upper")] <- "LOEUF"
  pdat <- merge(gcoords, constr.df[, c("gene", "LOEUF")], all=F, sort=F, by="gene")
  pdat <- pdat[which(!duplicated(pdat)), ]

  # Assign colors to genes based on LOEUF range consistent with gnomAD browser
  gcols <- sapply(pdat$LOEUF, function(x){
    if(x<0.33){"#FF2600"}else if(x<0.66){"#FF9300"}else if(x<1){"#FFC000"}else{ns.color}
  })

  # Scale LOEUF according to y0 and panel.height
  pdat$LOEUF <- sapply(pdat$LOEUF, function(x){1-min(c(1, x))})
  max.loeuf <- max(c(1, ceiling(11*max(pdat$LOEUF, na.rm=T)) / 10))
  pdat$y <- (pdat$LOEUF * panel.height / min(c(1, max.loeuf))) + y0 - half.height

  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")

  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+half.height, tick=F, las=2, line=-0.65,
       labels=0, cex.axis=y.ax.label.cex)
  axis(2, at=y0-half.height, tick=F, las=2, line=-0.65,
       labels=bquote("" >= 1), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, labels="LOEUF", las=2)

  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)

  # Add segments for each gene
  segments(x0=pdat$start, x1=pdat$end, y0=pdat$y, y1=pdat$y,
           lwd=4, lend="square", col=gcols)

  # Label genes, if optioned
  if(!is.null(label.genes)){
    sapply(label.genes[which(label.genes %in% pdat$gene)], function(gene){
      gidx <- which(pdat$gene == gene)
      if(pdat$LOEUF[gidx] <= 0.5 * max.loeuf){
        pos <- 3
        y.buf <- -0.1*panel.height
      }else{
        pos <- 1
        y.buf <- 0.1*panel.height
      }
      text(x=mean(as.numeric(pdat[gidx, c("start", "end")])),
           y=pdat$y[gidx]+y.buf, labels=gene, font=3, cex=5/6, pos=pos)
    })
  }
}


#' Plot PIPs for locus highlight
#'
#' Add plot of PIPs for all genes for locus highlight
#'
#' @param pips data.frame of PIPs per gene as loaded by [load.pips.for.genelist()]
#' @param genes data.frame of gene features to plot as loaded by [load.genes.from.gtf()]
#' @param y0 vertical midpoint for panel
#' @param panel.height relative height of panel \[default: 0.2\]
#' @param label.genes vector of gene symbols to label \[default: don't label genes\]
#' @param col color for plotting genes \[default: blueblack\]
#' @param highlight.genes vector of gene symbols to plot in a different color
#' \[default: plot all genes using `col`]
#' @param highlight.col color to apply to `highlight.genes` \[default: "red"\]
#'
#' @seealso [load.pips.for.genelist], [load.genes.from.gtf]
#'
#' @export plot.pips.for.highlight
#' @export
plot.pips.for.highlight <- function(pips, genes, y0, panel.height=0.2,
                                    label.genes=NULL, col=blueblack,
                                    highlight.genes=NULL, highlight.col="red"){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height

  # Merge gene coordinates & PIPs
  gcoords <- genes[which(genes$gene %in% pips$gene & genes$feature=="transcript"),
                   c("start", "end", "gene")]
  pdat <- merge(gcoords, pips[, c("gene", "PIP")], all=F, sort=F, by="gene")
  pdat <- pdat[which(!duplicated(pdat)), ]

  # Scale PIPs according to y0 and panel.height
  max.pip <- min(c(1, ceiling(11*max(pdat$PIP, na.rm=T)) / 10))
  pdat$y <- (pdat$PIP * panel.height / min(c(1, max.pip))) + y0 - half.height

  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")

  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+c(-half.height, half.height), tick=F, las=2, line=-0.65,
       labels=c(0, max.pip), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, labels="PIP", las=2)

  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)

  # Add segments for each gene
  if(!is.null(highlight.genes)){
    highlight.idxs <- which(pdat$gene %in% highlight.genes)
    other.idxs <- which(!(pdat$gene %in% highlight.genes))
    segments(x0=pdat$start[highlight.idxs], x1=pdat$end[highlight.idxs],
             y0=pdat$y[highlight.idxs], y1=pdat$y[highlight.idxs],
             lwd=4, lend="square", col=highlight.col)
  }else{
    other.idxs <- 1:nrow(pdat)
  }
  segments(x0=pdat$start[other.idxs], x1=pdat$end[other.idxs],
           y0=pdat$y[other.idxs], y1=pdat$y[other.idxs],
           lwd=4, lend="square", col=col)

  # Label any genes, if optioned
  if(!is.na(label.genes)){
    sapply(label.genes[which(label.genes %in% pdat$gene)], function(gene){
      best.pip <- pdat$PIP[which(pdat$gene==gene)]
      best.y <- pdat$y[which(pdat$gene==gene)]
      if(best.pip > 0.6 * max.pip){
        best.pos <- 1
        best.y <- best.y + 0.1*panel.height
      }else{
        best.pos <- 3
        best.y <- best.y - 0.1*panel.height
      }
      text(x=mean(as.numeric(pdat[which(pdat$gene==gene), c("start", "end")])),
           y=best.y, labels=gene, font=3, cex=5/6, pos=best.pos)
    })
  }
}


#' Plot simple rectangles for locus highlight
#'
#' Add panel of simple colored rectangles to locus highlight
#'
#' @param xlefts left-most position(s) for rectangle(s)
#' @param xrights right-most position(s) for rectangle(s)
#' @param col color(s) to fill rectangle(s)
#' @param border color(s) for rectangle outline(s)
#' @param panel.height relative height of panel \[default: 0.2\]
#' @param rect.height.cex relative vertical expansion scalar for rectangles
#' \[default: 1\]
#' @param y.axis.title title for Y axis \[default: no title\]
#'
#' @export plot.rects.for.highlight
#' @export
plot.rects.for.highlight <- function(xlefts, xrights, y0, col=blueblack,
                                     border=blueblack, panel.height=0.2,
                                     rect.height.cex=1, y.axis.title=NULL){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- rect.height.cex * (y0 - half.height)
  ytop <- rect.height.cex * (y0 + half.height)

  # Add rectangles
  rect(xleft=xlefts, xright=xrights, ybottom=ybottom, ytop=ytop, col=col, border=border)

  # Add Y-axis title
  axis(2, at=y0, tick=F, line=-0.9, las=2, labels=y.axis.title)
}


#' Plot ChromHMM tracks
#'
#' Add panel of ChromHMM tracks to locus highlight
#'
#' @param chmm.tracks list of chromhmm tracks to be plotted as loaded by []
#' @param y0 vertical midpoint for panel
#' @param panel.height relative height of panel \[default: 0.2\]
#' @param y.axis.title title for Y axis \[default: no title\]
#'
#' @seealso [load.chromhmm.tracks]
#'
#' @export plot.chromhmm.tracks
#' @export
plot.chromhmm.tracks <- function(chmm.tracks, y0, panel.height=0.2, y.axis.title=NULL){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height
  n.tracks <- length(chmm.tracks)
  track.breaks <- rev(seq(ybottom, ytop, length.out=n.tracks+1))
  track.ybottoms <- track.breaks[1:n.tracks]
  track.ytops <- track.breaks[(1:n.tracks)+1]

  # Prep background
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=ybottom, ytop=ytop,
       col="white", border=NA, bty="n")

  # Add rectangles
  sapply(1:n.tracks, function(i){
    track <- chmm.tracks[[i]]
    rect(xleft=track$start, xright=track$end,
         ybottom=track.ybottoms[i], ytop=track.ytops[[i]],
         border=track$color, col=track$color)
  })

  # Add Y-axis title
  axis(2, at=y0, tick=F, line=-0.8, las=2, labels=y.axis.title)

  # Add cleanup gridlines & border
  segments(x0=par("usr")[1], x1=par("usr")[2], y0=track.breaks, y1=track.breaks, col=bluewhite)
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=ybottom, ytop=ytop, col=NA, border=blueblack, xpd=T)
}


#' Plot quantitative bedgraph-style feature
#'
#' Add track of BEDGraph-style quantitative feature as coordinate-based barplot to locus highlight
#'
#' @param bed BED-style data.frame of elements to plot as loaded by
#' [load.feature.bed.for.highlight()]
#' @param y0 vertical midpoint for panel
#' @param col color to fill bars
#' @param panel.height relative height of panel \[default: 0.2\]
#' @param y.axis.title title for Y axis \[default: no title\]
#'
#' @seealso [load.feature.bed.for.highlight]
#'
#' @export plot.features.from.bedgraph
#' @export
plot.features.from.bedgraph <- function(bed, y0, col=blueblack, panel.height=0.2,
                                        ytitle=NULL){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height

  # Scale feature values according to y0 and panel.height
  vals.orig <- as.numeric(bed$value)
  max.vals.orig <- max(vals.orig, na.rm=T)
  vals.scaled <- (panel.height / max.vals.orig) * vals.orig
  vals <- vals.scaled + y0 - half.height

  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")

  # Add rectangles
  rect(xleft=bed$start, xright=bed$end, ybottom=ybottom, ytop=vals,
       col=col, border=col)

  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+c(-half.height, half.height), tick=F, las=2, line=-0.65,
       labels=c(0, ceiling(max.vals.orig)+1), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, labels=ytitle, las=2)

  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)
}


#' Plot CNV panel for locus highlight
#'
#' Add mirrored case & control CNV pileups for a single cohort to a locus highlight
#'
#' @param cnvs cnvs to plot as loaded by [load.cnvs.from.region()]
#' @param n.case case sample size
#' @param n.ctrl control sample size
#' @param y0 vertical midpoint for panel
#' @param cnv.type type of CNV to plot
#' @param highlight.hpo highlight CNVs from this HPO in a different color
#' \[default: plot all CNVs in the same color\]
#' @param max.freq maximum frequency to display on Y-axis \[default: auto-set\]
#' @param start left-most coordinate to plot \[default: plot all CNVs\]
#' @param end right-most coordinate to plot \[default: plot all CNVs\]
#' @param y.axis.title title for Y axis \[default: "CNV\\nFreq."\]
#' @param expand.pheno.label boolean indicator whether to expand phenotype label \[default: TRUE\]
#' @param case.legend.side side to plot case legend \[default: "left"\]
#' @param case.legend.topbottom vertical alignment for case legend \[default "top"]
#' @param ctrl.legend.side side to plot control legend \[default: "left"\]
#' @param cc.legend.colors vector of colors to apply to case & control legends \[default: blueblack\]
#' @param cnv.lwd line width for each CNV \[default: 0.05\]
#' @param add.cohort.label boolean indicator to print `cohort.label` on plot \[default: FALSE\]
#' @param cohort.label text of cohort label \[default: no label\]
#' @param panel.height relative height of panel \[default: 2\]
#' @param dx number of steps across X axis \[default: 100\]
#'
#' @details `case.legend.side` and `control.legend.side` both accept either
#' "top" or "bottom".
#' Vertical alignment of control legend will automatically be set to opposite
#' of `case.legend.topbottom`.
#'
#' @seealso [load.cnvs.from.region], [pileup.cnvs.for.highlight]
#'
#' @export plot.cnv.panel.for.highlight
#' @export
plot.cnv.panel.for.highlight <- function(cnvs, n.case, n.ctrl, y0, cnv.type,
                                         highlight.hpo=NA, max.freq=NULL,
                                         start=NULL, end=NULL, y.axis.title="CNV\nFreq.",
                                         expand.pheno.label=TRUE, case.legend.side="left",
                                         case.legend.topbottom="top", ctrl.legend.side="left",
                                         cc.legend.colors=rep(blueblack, 2), cnv.lwd=0.05,
                                         add.cohort.label=FALSE, cohort.label=NULL,
                                         panel.height=2, dx=100){
  # Standardize inputs
  n.case <- as.numeric(n.case)
  n.ctrl <- as.numeric(n.ctrl)

  # Set horizontal dimensions
  if(is.null(start)){
    start <- par("usr")[1]
  }
  if(is.null(end)){
    end <- par("usr")[2]
  }

  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    col.case.highlight <- cnv.blacks[1]
    col.case.other <- cnv.colors[1]
    col.ctrl <- control.cnv.colors[1]
    col.midline <- cnv.blacks[1]
  }else if(cnv.type=="DUP"){
    col.case.highlight <- cnv.blacks[2]
    col.case.other <- cnv.colors[2]
    col.ctrl <- control.cnv.colors[2]
    col.midline <- cnv.blacks[2]
  }

  # Collect raw (unscaled) CNV pileups
  raw.pileups <- lapply(cnvs, pileup.cnvs.for.highlight, start=start, end=end, dx=dx)
  max.case.n <- max(unlist(lapply(raw.pileups$case$cnvs, function(l){max(l$y)})))
  max.ctrl.n <- max(unlist(lapply(raw.pileups$ctrl$cnvs, function(l){max(l$y)})))

  # Set scaling
  half.height <- panel.height / 2
  ytop <- y0 + half.height
  ybottom <- y0 - half.height
  if(is.null(max.freq)){
    max.case.freq <- (max.case.n + 1) / n.case
    max.ctrl.freq <- (max.ctrl.n + 1) / n.ctrl
    max.freq <- max(c(max.case.freq, max.ctrl.freq))
  }
  case.cnv.height <- half.height / ceiling(n.case * max.freq)
  ctrl.cnv.height <- half.height / ceiling(n.ctrl * max.freq)

  # Gather scaled CNV pileups
  case.pileup <- pileup.cnvs.for.highlight(cnvs$case, start=start, end=end, dx=dx,
                                           cnv.height=case.cnv.height,
                                           col=col.case.other,
                                           highlight.hpo=highlight.hpo,
                                           highlight.col=col.case.highlight)
  ctrl.pileup <- pileup.cnvs.for.highlight(cnvs$ctrl, start=start, end=end, dx=dx,
                                           cnv.height=ctrl.cnv.height, col=col.ctrl)

  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=7)
  abline(h=c(y0, y0 + y.ax.tick.spacing, y0 - y.ax.tick.spacing), col="white")

  # Plot midline, pileups, and outlines
  lapply(case.pileup$cnvs, function(l){polygon(l$x, y0 + l$y, border="white",
                                               col=l$color, lwd=cnv.lwd)})
  points(case.pileup$counts[, 1], case.pileup$counts[, 2] + y0, type="l", col=col.case.other)
  lapply(ctrl.pileup$cnvs, function(l){polygon(l$x, y0 - l$y, border="white",
                                               col=l$color, lwd=cnv.lwd)})
  points(ctrl.pileup$counts[, 1], -ctrl.pileup$counts[, 2] + y0, type="l", col=col.ctrl)
  abline(h=y0, col=col.midline)

  # Add Y axes
  y.ax.label.cex <- 5/6
  axis(2, at=y0 + y.ax.tick.spacing, labels=NA, col=blueblack, tck=-0.0075)
  max.freq.fmt <- format.scientific(max.freq, nsmall=0, max.decimal=0)
  axis(2, at=y0-half.height, tick=F, las=2, line=-0.65, labels=max.freq.fmt,
       cex.axis=y.ax.label.cex)
  axis(2, at=y0, tick=F, las=2, line=-0.65, labels=0,
       cex.axis=y.ax.label.cex)
  axis(2, at=y0+half.height, tick=F, las=2, line=-0.65, labels=max.freq.fmt,
       cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=0.25, tick=F, labels=y.axis.title)

  # Add case legend
  legend.y.cex <- 0.7
  legend.text.cex <- 5/6
  if(!is.na(case.legend.side)){
    case.legend.y <- y0 + (legend.y.cex * half.height)
    if(case.legend.side == "left"){
      case.legend.x <- start
      case.legend.pos <- 4
    }else{
      case.legend.x <- end
      case.legend.pos <- 2
    }
    if(expand.pheno.label == TRUE){
      case.legend.label <- paste(hpo.abbrevs[case.hpo], " cases (N=", prettyNum(n.case, big.mark=","), ")", sep="")
    }else{
      case.legend.label <- paste("Cases (N=", prettyNum(n.case, big.mark=","), ")", sep="")
    }

    text(x=case.legend.x, y=case.legend.y, labels=case.legend.label,
         pos=case.legend.pos, col=cc.legend.colors[1], cex=legend.text.cex)
  }

  # Add control legend
  if(!is.na(ctrl.legend.side)){
    ctrl.legend.y <- y0 - (legend.y.cex * half.height)
    if(ctrl.legend.side == "left"){
      ctrl.legend.x <- start
      ctrl.legend.pos <- 4
    }else{
      ctrl.legend.x <- end
      ctrl.legend.pos <- 2
    }
    ctrl.legend.label <- paste("Controls (N=", prettyNum(n.ctrl, big.mark=","), ")", sep="")
    text(x=ctrl.legend.x, y=ctrl.legend.y, labels=ctrl.legend.label,
         pos=ctrl.legend.pos, col=cc.legend.colors[2], cex=legend.text.cex)
  }

  # Add cohort label (upper panel, opposite case.legend.side)
  if(add.cohort.label){
    if(case.legend.topbottom=="top"){
      cohort.label.y.adj <- 1
    }else{
      cohort.label.y.adj <- -1
    }
    cohort.label.y <- y0 + (cohort.label.y.adj * legend.y.cex * half.height)
    if(case.legend.side=="left"){
      cohort.label.x <- end
      cohort.label.pos <- 2
    }else{
      cohort.label.x <- start
      cohort.label.pos <- 4
    }
    text(x=cohort.label.x, y=cohort.label.y, labels=cohort.label,
         pos=cohort.label.pos, col=blueblack, cex=legend.text.cex)
  }

  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)
}


#' Plot CNV key for locus highlight
#'
#' Add CNV key to locus highlight
#'
#' @param cnv.type type of CNV
#' @param y0 midpoint for key panel
#' @param total.n.ctrls total number of all controls
#' @param all.case.hpos all HPOs to be considered as cases
#' @param total.n.cases total number of all cases
#' @param highlight.case.hpo HPO to be highlighted \[default: no highlight\]
#' @param total.n.cases.highlight total number of cases matching `highlight.case.hpo`
#' \[default no highlight\]
#' @param panel.height relative height of key panel \[default: 0.2\]
#' @param text.cex character expansion scalar for text \[default: 5/6\]
#' @param pt.cex expansion scalar for points \[default: 1.3\]
#'
#' @export plot.cnv.key.for.highlight
#' @export
plot.cnv.key.for.highlight <- function(cnv.type, y0, total.n.ctrls, all.case.hpos,
                                       total.n.cases, highlight.case.hpo=NA,
                                       total.n.cases.highlight=NA, panel.height=0.2,
                                       text.cex=5/6, pt.cex=1.3){
  # Get plot data
  x.at <- par("usr")[1] + (c(0.075, 0.65) * diff(par("usr")[1:2]))
  y.buf <- panel.height/3

  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    col.case.highlight <- cnv.blacks[1]
    col.case.other <- cnv.colors[1]
    col.ctrl <- control.cnv.colors[1]
    cnv.label <- "Deletions:"
  }else if(cnv.type=="DUP"){
    col.case.highlight <- cnv.blacks[2]
    col.case.other <- cnv.colors[2]
    col.ctrl <- control.cnv.colors[2]
    cnv.label <- "Duplications:"
  }

  # Allow everything to plot beyond the boundaries
  par(xpd=T)

  # Add left-most text
  axis(2, at=y0, tick=F, line=-2, las=2, labels=cnv.label, cex=text.cex)

  # Add case labels
  all.case.label <- paste(hpo.abbrevs[all.case.hpos[1]], " cases (total N=",
                          prettyNum(total.n.cases, big.mark=","), ")", sep="")
  all.case.label <- gsub("cases cases", "cases", all.case.label, fixed=T)
  if(is.na(highlight.case.hpo)){
    points(x=x.at[1], y=y0, pch=15, col=col.case.other, cex=pt.cex)
    text(x.at[1], y=y0, pos=4, xpd=T, cex=text.cex, labels=all.case.label)
  }else{
    points(x=x.at[1], y=y0+y.buf, pch=15, col=col.case.other, cex=pt.cex)
    text(x.at[1], y=y0+y.buf, pos=4, xpd=T, cex=text.cex, labels=all.case.label)
    points(x=x.at[1], y=y0-y.buf, pch=15, col=col.case.highlight, cex=pt.cex)
    text(x.at[1], y=y0-y.buf, pos=4, xpd=T, cex=text.cex,
         labels=paste(hpo.abbrevs[highlight.case.hpo], " cases (total N=",
                      prettyNum(total.n.cases.highlight, big.mark=","), ")", sep=""))
  }

  # Add control label
  ctrl.label <- paste("Controls (total N=",
                      prettyNum(total.n.ctrls, big.mark=","),
                      ")", sep="")
  points(x=x.at[2], y=y0, pch=15, col=col.ctrl, cex=pt.cex)
  text(x.at[2], y=y0, pos=4, xpd=T, cex=text.cex, labels=ctrl.label)
}

