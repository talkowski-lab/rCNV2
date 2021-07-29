#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Functions for Manhattan & QQ analyses of association statistics


#' Load association statistics
#'
#' Loads a generic input file of locus-based association statistics
#'
#' @param stats.in path to input file of association statistics
#' @param p.col.name name of column in header for P-value column \[default: 'p'\]
#' @param p.is.phred boolean indicator if P-values as provided are -log10 scaled \[default: FALSE\]
#' @param min.p minimum P-value to retain as numerically meaningful \[default: 10^-100\]
#'
#' @return data frame with three columns: chromosome, position, and P-value
#'
#' @export load.manhattan.stats
#' @export
load.manhattan.stats <- function(stats.in, p.col.name="p", p.is.phred=F, min.p=10^-100){
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

  # Cap p-values at global min
  p[which(p<min.p)] <- min.p

  data.frame("chr"=chr, "pos"=pos, "p"=p)
}


#' Manhattan plot
#'
#' Generate a Manhattan plot from association statistics
#'
#' @param df data frame of association statistics as read by [read.manhattan.stats()]
#' @param cutoff smallest P-value to render \[default: 10^-8\]
#' @param highlights data frame of regions to highlight
#' @param highlight.color color to use for highlighted points
#' @param highlight.name label for highlighted points
#' @param mark.highlights boolean indicated to mark highlighted regions with small arrows \[default: FALSE\]
#' @param colors vector of colors to assign sequentially to all chromosomes
#' @param cutoff.color color to use for P-value cutoff line
#' @param lab.prefix prefix to prepend to various labels
#' @param ymax maximum Y-value to plot \[default: same as smallest P-value\]
#' @param reflection boolean indicator to invert plot (useful for Miami plots)
#' @param pt.cex scaling factor for points \[default: 0.2\]
#' @param label.cex scaling factor for label text
#' @param y.ax.cex scaling factor for labels on Y-axis \[default: 0.75\]
#'
#' @return None
#'
#' @export plot.manhattan
#' @export
plot.manhattan <- function(df, cutoff=1e-08, highlights=NULL,
                           highlight.color="#4EE69A",
                           highlight.name="Highlighted Loci",
                           mark.highlights=FALSE,
                           colors=c("gray15", "gray70"),
                           cutoff.color="#D01C8B",
                           lab.prefix=NULL,
                           ymax=NULL, reflection=F,
                           pt.cex=0.2, label.cex=1, y.ax.cex=0.75){
  contigs <- unique(df[, 1])
  contigs <- contigs[which(!(is.na(contigs)))]

  indexes <- as.data.frame(t(sapply(contigs, function(chr){
    return(c(chr, 0, max(df[which(df[, 1] == chr), 2])))
  })))
  indexes$sum <- cumsum(indexes[, 3])
  indexes$bg <- colors
  indexes[, 2:4] <- apply(indexes[, 2:4], 2, as.numeric)

  df.plot <- as.data.frame(t(apply(df, 1, function(row){
    contig.idx <- which(indexes[, 1]==row[1])
    return(c(as.character(row[1]),
             as.numeric(row[2]) + indexes[contig.idx, 4] - indexes[contig.idx, 3],
             as.numeric(row[3]),
             indexes[contig.idx, 5]))
  })))
  df.plot[, 2] <- as.numeric(as.character(df.plot[, 2]))
  df.plot[, 3] <- as.numeric(as.character(df.plot[, 3]))
  colnames(df.plot) <- c("chr", "pos", "p", "color")

  # Drop rows with NA p-values
  df.plot <- df.plot[which(!is.na(df.plot$p)), ]

  if(is.null(ymax)){
    ymax <- -log10(df[which(!(is.na(df[, 3]))), 3])
    ymax <- max(max(ymax[which(!(is.infinite(ymax)))], na.rm=T), -log10(cutoff) + 2, na.rm=T)
  }

  if(reflection == F){
    df.plot$phred <- -log10(df.plot[, 3])
    log.cutoff <- -log10(cutoff)
    legend.pos <- "top"
    x.ax <- 1
    mars <- c(2.1, 3.1, 0.5, 1)
    x.title <- T
  }else{
    df.plot$phred <- log10(df.plot[, 3])
    log.cutoff <- log10(cutoff)
    legend.pos <- "bottom"
    x.ax <- 3
    ymax <- -ymax
    mars <- c(0.5, 3.1, 0.6, 1)
    x.title <- F
  }

  par(mar=mars, bty="n")
  plot(x=c(0, max(indexes[, 4])), y=c(0, 1.1 * ymax), type="n",
       yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="")

  if(!is.null(highlights)){
    highlights <- t(apply(highlights, 1, function(coords){
      if(coords[1] %in% contigs){
        if(which(as.character(coords[1]) == contigs) > 1){
          coords[2] <- coords[2] + indexes[as.integer(coords[1]) - 1, 4]
          coords[3] <- coords[3] + indexes[as.integer(coords[1]) - 1, 4]
        }
      }
      return(coords)
    }))

    if(mark.highlights == TRUE){
      min.hbin.width <- 0.005 * (par("usr")[2] - par("usr")[1])
      highlight.bins <- t(apply(highlights, 1, function(coords){
        if(coords[3] - coords[2] < min.hbin.width){
          mid <- mean(coords[2:3])
          coords[2] <- mid - (0.5 * min.hbin.width)
          coords[3] <- mid + (0.5 * min.hbin.width)
        }
        return(coords)
      }))
      rect(xleft=highlight.bins[,2], xright=highlight.bins[,3],
           ybottom=par("usr")[3], ytop=par("usr")[4], bty="n", border=NA,
           col=adjustcolor(highlight.color, alpha=0.1))
      # rect(xleft=highlight.bins[,2], xright=highlight.bins[,3],
      #      ybottom=par("usr")[4] - (0.01 * (par("usr")[4] - par("usr")[3])),
      #      ytop=par("usr")[4], bty="n", border=NA,
      #      col=highlight.color)
    }

    hits <- sort(unique(unlist(apply(highlights, 1, function(coords){
      which(df.plot$chr == as.character(coords[1])
            & df.plot$pos <= as.numeric(coords[3])
            & df.plot$pos >= as.numeric(coords[2]))
    }))))
    df.plot$color[hits] <- highlight.color
  }

  abline(h=log.cutoff, col=cutoff.color, lty=2)

  points(df.plot[, 2], df.plot[, 5],
         cex=pt.cex, pch=19,
         col=as.character(df.plot[, 4]))

  if(!is.null(highlights)){
    if(any(df.plot$color == highlight.color)){
      points(df.plot[which(df.plot$color == highlight.color), 2],
             df.plot[which(df.plot$color == highlight.color), 5],
             cex=0.2, pch=19,
             col=highlight.color)
    }
    legend(legend.pos, pch=19, col=highlight.color, cex=0.5,
           legend=paste(highlight.name, " (n=",
                        prettyNum(nrow(highlights), big.mark=","),
                        ")", sep=""))
  }

  axis(x.ax, at=c(0, max(indexes$sum)), tck=0, labels=NA)
  midpoints <- sapply(1:length(indexes[, 2]), function(i){
    return(mean(c(indexes[i, 4], indexes[i, 4] - indexes[i, 3])))
  })
  xlab.frac <- 0.075
  xlab.bins <- seq(0, max(indexes$sum), xlab.frac * max(indexes$sum))
  xlab.contigs <- unique(sapply(xlab.bins, function(pos){
    # which(pos<=indexes$sum & pos>=cumsum(c(0,indexes[-1, 3])))[1]
    which(abs(midpoints-pos)==min(abs(midpoints-pos)))
  }))
  long.contigs <- which(indexes[, 3]/max(indexes$sum) > 0.0475)
  xlab.contigs <- unique(c(xlab.contigs, long.contigs))
  sapply(xlab.contigs, function(k){
    axis(x.ax, at=midpoints[k], labels=indexes[k, 1],
         tick=F, line=-1.1, cex.axis=0.75,
         col.axis=indexes$bg[k])
  })
  if(x.title == T){
    mtext(x.ax, text="Chromosome", line=0.95, cex=label.cex)
  }

  if(reflection == T){
    y.at <- seq(0, floor(par("usr")[3]), by=floor(par("usr")[3]/6))
  }else{
    y.at <- seq(0, ceiling(par("usr")[4]), by=ceiling(par("usr")[4]/6))
  }
  axis(2, at=y.at, labels=NA, tck=-0.02)
  axis(2, at=y.at, tick=F, line=-0.5, labels=abs(y.at),
       cex.axis=y.ax.cex, las=2)
  if(!is.null(lab.prefix)){
    mtext(2, text=bquote(-log[10](italic(P)) ~ .(lab.prefix)),
          line=1.5, cex=label.cex)
  }else{
    mtext(2, text=bquote(-log[10](italic(P))),
          line=1.5, cex=label.cex)
  }
}


#' Quantile-quantile plot
#'
#' Generate a Q-Q plot of association stats vs. uniform null
#'
#' @param stats data frame of association statistics as read by [read.manhattan.stats()]
#' @param cutoff smallest P-value to render
#' @param highlights data frame of regions to highlight
#' @param highlight.color color to use for highlighted points
#' @param highlight.name label for highlighted points
#' @param print.stats boolean indicator to print stats on plot \[default: TRUE\]
#' @param echo.lambdas boolean indicator to print genomic control statistic to
#' stdout (for logging purposes) \[default: FALSE\]
#' @param legend boolean indicator to add legend \[default: TRUE\]
#' @param ymax maximum Y-value to plot \[default: same as smallest P-value\]
#' @param reflection boolean indicator to invert plot (useful for Miami plots)
#' @param label.cex scaling factor for label text
#'
#' @return None
#'
#' @export plot.qq
#' @export
plot.qq <- function(stats, cutoff=NULL, highlights=NULL,
                    highlight.color="#4EE69A",
                    highlight.name="Positive Controls",
                    print.stats=T, echo.lambdas=F,
                    legend=T, ymax=NULL, reflection=F,
                    label.cex=1){

  # Dummy function to plot N/A QQ for analyses with no cases
  empty.qq <- function(){
    par(mar=c(2.1, 3.1, 0.6, 1), bty="n")
    plot(x=0:1, y=0:1, type="n", xaxt="n", yaxt="n",
         xlab="", ylab="", xaxs="i", yaxs="i")
    axis(1, at=axTicks(1), labels=NA, tck=-0.02)
    mtext(1, text=expression(Expected ~ ~-log[10](italic(p))),
          line=1.2, cex=label.cex)
    axis(2, at=axTicks(2), labels=NA, tck=-0.02)
    mtext(2, text=expression(Observed ~ ~-log[10](italic(p))),
          line=1.5, cex=label.cex)
    text(x=0.5,y=0.5,labels="N/A")
  }

  # Main QQ function
  p <- as.numeric(stats$p)
  if(all(p==1 | is.na(p) | is.nan(p) | is.infinite(p))){
    empty.qq()
  }else{
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
                         is.finite(p) & p <= 1 & p >= 0)
    p <- p[keep.idxs]
    colors <- colors[keep.idxs]

    new.order <- order(p)
    p <- p[new.order]
    colors <- colors[new.order]
    hits <- which(colors == highlight.color)

    expected <- ppoints(length(p))
    qqconf <- function (p.expected, reflection=F){
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
      if(reflection == T){
        mpts[, 2] <- -mpts[, 2]
      }
      return(mpts)
    }
    conf.int <- qqconf(expected, reflection)

    lambda <- dchisq(median(p), df=1)/dchisq(median(expected), df=1)
    if(!is.null(highlights)){
      expected.hits <- ppoints(length(p[hits]))
      lambda.hits <- dchisq(median(p[hits]), df=1)/dchisq(median(expected.hits), df=1)
      lambda.ratio <- lambda.hits / lambda
      if(echo.lambdas==T){
        cat(paste(round(lambda, 4), "\t",
                  round(lambda.hits, 4), "\t",
                  round(lambda.ratio, 4), "\n",
                  sep=""))
      }
    }

    if(is.null(ymax)){
      ymax <- max(max(-log10(p[which(!(is.infinite(p)))])), -log10(cutoff) + 2)
    }

    if(is.null(cutoff)){
      cutoff <- 0.05/length(p)
    }

    expected <- -log10(expected)

    if(reflection == F){
      p <- -log10(p)
      log.cutoff <- -log10(cutoff)
      ab.end <- 1
      x.ax <- 1
      mars <- c(2.1, 3.1, 0.5, 1)
      lambda.ypos <- c(0.9, 0.825, 0.75)
      x.title <- T
    }else{
      p <- log10(p)
      log.cutoff <- log10(cutoff)
      ab.end <- -1
      ymax <- -ymax
      x.ax <- 3
      mars <- c(0.5, 3.1, 0.6, 1)
      lambda.ypos <- rev(c(0.9, 0.825, 0.75))
      x.title <- F
    }

    par(mar=mars, bty="n")
    plot(x=expected, y=p, type="n", xaxt="n", yaxt="n",
         xlab="", ylab="", xaxs="i", yaxs="i",
         xlim=c(0, 1.1 * max(expected)), ylim=range(c(0, 1.1 * ymax)))
    polygon(x=conf.int[, 1], y=conf.int[, 2], col="gray90", border=NA)
    abline(0, ab.end, col="gray50")
    if(!is.null(cutoff)){
      abline(h=log.cutoff, col="#D01C8B", lty=2)
    }

    if (print.stats == T){
      xpos.adjust <- 0.025
      stat.label.color <- "gray20"
      if(reflection == T){
        ypos.ref <- par("usr")[3]
      }else{
        ypos.ref <- par("usr")[4]
      }
      if(is.null(highlights)){
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[1] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda ~ .(paste("=", sprintf("%.2f", lambda), sep=""))))
      }else{
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[1] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda[italic("All")] ~ .(paste("=", sprintf("%.2f", lambda), sep=""))))
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[2] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda[italic("Cyan")] ~ .(paste("=", sprintf("%.2f", lambda.hits), sep=""))))
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[3] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda[italic("C")] / lambda[italic("A")] ~ .(paste("=", sprintf("%.2f", lambda.ratio), sep=""))))
      }
    }

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

    axis(x.ax, at=axTicks(1), labels=NA, tck=0)
    axis(x.ax, at=axTicks(1), tick=F, line=-1.1, cex.axis=0.75, labels=abs(axTicks(1)))
    if(x.title == T){
      mtext(x.ax, text=expression(Expected ~ ~-log[10](italic(p))),
            line=1.15, cex=label.cex)
    }

    if(reflection == T){
      y.at <- seq(0, floor(par("usr")[3]), by=floor(par("usr")[3]/6))
    }else{
      y.at <- seq(0, ceiling(par("usr")[4]), by=ceiling(par("usr")[4]/6))
    }
    axis(2, at=y.at, labels=NA, tck=-0.02)
    axis(2, at=y.at, cex.axis=0.75, tick=F, line=-0.6, las=2, labels=abs(y.at))
    mtext(2, text=expression(Observed ~ ~-log[10](italic(p))),
          line=1.5, cex=label.cex)
  }
}

