#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot pHaplo/pTriplo score scatterplot for gene score analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000, family="sans")


##########################
### PLOTTING FUNCTIONS ###
##########################
# Prepare color grid for gradient-based scatter
make.color.grid <- function(){
  yl.pal <- colorRampPalette(c(ns.color.light, control.cnv.colors[2],
                               cnv.colors[2]))(101)
  yr.pal <- colorRampPalette(c(cnv.colors[1],  cnv.colors[3]))(101)
  do.call("rbind", lapply(1:101, function(y){
    colorRampPalette(c(yl.pal[y], yr.pal[y]))(101)
  }))
}

# Get gene color based on (pHaplo, pTriplo) x,y pairs
query.color.grid <- function(gene, scores, color.grid){
  x <- round(100*as.numeric(scores$pHaplo[which(scores$gene==gene)]))
  y <- round(100*as.numeric(scores$pTriplo[which(scores$gene==gene)]))
  color.grid[y+1, x+1]
}

# Helper function to add marginal density plot
marginal.density <- function(vals, colors, scale=0.25, bw=0.02, rotate=F,
                             hc.cutoff=0.9, lc.cutoff=0.5, gradient=F, pal=NA){
  # Gather plot values
  d <- density(vals, bw=bw)
  idx <- d$x
  idx[which(idx<0)] <- 0; idx[which(idx>1)] <- 1
  h <- d$y
  h <- scale*(1/max(h, na.rm=T))*h
  h <- h + 1

  # Split indexes into bins
  if(gradient==FALSE){
    ns.idx <- which(idx<lc.cutoff)
    lc.idx <- which(idx>=lc.cutoff & idx<hc.cutoff)
    hc.idx <- which(idx>=hc.cutoff)
  }

  # Compute polygon coordinate vectors based on rotation
  if(rotate==T){
    x <- h; y <- idx
    if(gradient==FALSE){
      ns.df <- data.frame("x"=c(x[ns.idx], rep(1, length(ns.idx))),
                          "y"=c(y[ns.idx], rev(y[ns.idx])))
      lc.df <- data.frame("x"=c(x[lc.idx], rep(1, length(lc.idx))),
                          "y"=c(y[lc.idx], rev(y[lc.idx])))
      hc.df <- data.frame("x"=c(x[hc.idx], rep(1, length(hc.idx))),
                          "y"=c(y[hc.idx], rev(y[hc.idx])))
    }
    all.df <- data.frame("x"=c(x, rep(1, length(x))), "y"=c(y, rev(y)))
  }else{
    x <- idx; y <- h
    if(gradient==FALSE){
      ns.df <- data.frame("x"=c(x[ns.idx], rev(x[ns.idx])),
                          "y"=c(y[ns.idx], rep(1, length(ns.idx))))
      lc.df <- data.frame("x"=c(x[lc.idx], rev(x[lc.idx])),
                          "y"=c(y[lc.idx], rep(1, length(lc.idx))))
      hc.df <- data.frame("x"=c(x[hc.idx], rev(x[hc.idx])),
                          "y"=c(y[hc.idx], rep(1, length(hc.idx))))
    }
    all.df <- data.frame("x"=c(x, rev(x)), "y"=c(y, rep(1, length(y))))
  }

  # Plot polygons in ascending order
  if(gradient==FALSE){
    polygon(x=ns.df$x, y=ns.df$y, border=colors[3], col=colors[3], xpd=T)
    polygon(x=lc.df$x, y=lc.df$y, border=colors[2], col=colors[2], xpd=T)
    polygon(x=hc.df$x, y=hc.df$y, border=colors[1], col=colors[1], xpd=T)
  }else{
    n.poly <- length(x)-1
    pal <- pal(n.poly)
    sapply(1:n.poly, function(i){
      if(rotate==TRUE){
        polygon(x=x[c(1, i, i+1, 1)], y=y[c(i, i, i+1, i+1)],
                border=pal[i], col=pal[i], xpd=T)
      }else{
        polygon(x=x[c(i, i, i+1, i+1)], y=y[c(1, i, i+1, 1)],
                border=pal[i], col=pal[i], xpd=T)
      }
    })
  }
  polygon(x=all.df$x, y=all.df$y, col=NA, border=colors[1], bty="n", xpd=T)
}

# Scatterplot of genes by scores with options for coloring and marginal densities
scores.scatterplot <- function(scores, pt.colors, add.cor=TRUE,
                               hi.hc=0.9, hi.lc=0.5, ts.hc=0.9, ts.lc=0.5,
                               margin.dens.height=0.175, margin.dens.gradient=FALSE,
                               pt.cex=0.1, pt.alpha=1, ax.tick=-0.03,
                               bg.col="white", gridlines.col=bluewhite,
                               category.lines.col=blueblack,
                               x.ax.title="Haploinsufficiency Score (pHaplo)",
                               y.ax.title="Triplosensitivity Score (pTriplo)",
                               ax.at=seq(0, 1, 0.2), parmar=c(2.7, 2.7, 1, 1)){
  # Prep plot area
  par(mar=parmar, bty="n", family="sans")
  plot(NA, xlim=c(0, 1 + margin.dens.height), ylim=c(0, 1 + margin.dens.height),
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  rect(xleft=par("usr")[1], xright=1,
       ybottom=par("usr")[3], ytop=1,
       border=NA, bty="n", col=bg.col)
  segments(x0=ax.at, x1=ax.at,
           y0=par("usr")[3], y1=1, col=gridlines.col)
  segments(x0=par("usr")[1], x1=1,
           y0=ax.at, y1=ax.at, col=gridlines.col)

  # Add points
  if(pt.alpha < 1){
    pt.colors <- sapply(pt.colors, adjustcolor, alpha=pt.alpha)
    points(scores$pHaplo, scores$pTriplo, cex=pt.cex, col=pt.colors, pch=19)
  }else{
    points(scores$pHaplo, scores$pTriplo, cex=pt.cex, col=pt.colors, pch=19)
  }

  # Add marginal densities
  if(margin.dens.height > 0){
    if(margin.dens.gradient==TRUE){
      marginal.density(scores$pHaplo, colors=cnv.colors[1], gradient=TRUE,
                       pal=colorRampPalette(c(ns.color.light, cnv.whites[1],
                                              control.cnv.colors[1], cnv.colors[1])),
                       bw=0.01, scale=margin.dens.height)
      marginal.density(scores$pTriplo, colors=cnv.colors[2], gradient=TRUE,
                       pal=colorRampPalette(c(ns.color.light, cnv.whites[2],
                                              control.cnv.colors[2], cnv.colors[2])),
                       bw=0.01, scale=margin.dens.height, rotate=T)
    }else{
      marginal.density(scores$pHaplo, hc.cutoff=hi.hc, lc.cutoff=hi.lc,
                       colors=c(cnv.colors[1], control.cnv.colors[1], redwhite),
                       bw=0.01, scale=margin.dens.height)
      marginal.density(scores$pTriplo, hc.cutoff=ts.hc, lc.cutoff=ts.lc,
                       colors=c(cnv.colors[2], control.cnv.colors[2], bluewhite),
                       bw=0.01, scale=margin.dens.height, rotate=T)
    }
  }

  # Add delimiting lines
  segments(x0=c(hi.hc, hi.lc, hi.hc, hi.lc),
           x1=c(hi.hc, hi.lc, hi.hc, hi.lc),
           y0=c(ts.hc, ts.lc, rep(par("usr")[3], 2)),
           y1=c(1, 1, rep(ts.lc, 2)),
           lwd=1, lend="round", col=category.lines.col, lty=2)
  segments(y0=c(ts.hc, ts.lc, ts.hc, ts.lc),
           y1=c(ts.hc, ts.lc, ts.hc, ts.lc),
           x0=c(hi.hc, hi.lc, rep(par("usr")[3], 2)),
           x1=c(1, 1, rep(hi.lc, 2)),
           lwd=1, lend="round", col=category.lines.col, lty=2)

  # Add axes
  axis(1, at=ax.at, col=blueblack, labels=NA, tck=ax.tick)
  axis(2, at=ax.at, col=blueblack, labels=NA, tck=ax.tick)
  sapply(ax.at, function(i){
    axis(1, at=i, col=blueblack, tick=F, line=-0.65)
    axis(2, at=i, col=blueblack, tick=F, line=-0.65, las=2)
  })
  axis(1, at=0.5, tick=F, line=0.3, labels=x.ax.title)
  axis(2, at=0.5, tick=F, line=0.75, labels=y.ax.title)

  # Add correlation coefficient
  if(add.cor==TRUE){
    r2 <- cor(scores$pHaplo, scores$pTriplo)^2
    r2.fmt <- formatC(round(r2, 2), small.interval=2)
    text(x=1.05+(margin.dens.height/2),
         y=1.065+(margin.dens.height/2),
         xpd=T, cex=0.9, srt=45,
         labels=bquote(italic(R)^2 * "=" * .(r2.fmt)))
  }

  # Cleanup
  rect(xleft=par("usr")[1], xright=1,
       ybottom=par("usr")[3], ytop=1,
       border=blueblack, bty="o", col=NA, xpd=T)
}

# Plot legend of genes per group
plot.scores.scatter.legend <- function(scores, ds.groups){
  # Get plot data
  colors <- c(cnv.colors[3], control.cnv.colors[3],
              cnv.colors[1], control.cnv.colors[1],
              cnv.colors[2], control.cnv.colors[2],
              ns.color)

  # Prep plot area
  par(mar=c(0, 0.1, 1.1, 0.1), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(7, 0),
       xaxt="n", xaxs="i", xlab="", yaxt="n", yaxs="i", ylab="")

  # Add points
  points(x=rep(0.05, 7), y=0.5:6.5, pch=18, col=colors, cex=2.25, xpd=T)

  # Add N genes
  gene.ax.at <- c(0.175, 0.35)
  axis(3, at=gene.ax.at, tck=0, labels=NA, col=blueblack)
  axis(3, at=gene.ax.at[2], tick=F, line=-1, labels="Genes", hadj=1)
  sapply(1:length(ds.groups), function(i){
    text(x=0.4, y=i-0.5, pos=2, labels=prettyNum(length(ds.groups[[i]]), big.mark=","))
  })

  # Add labels
  cat.ax.at <- c(0.4, 0.95)
  axis(3, at=cat.ax.at, tck=0, labels=NA, col=blueblack)
  axis(3, at=cat.ax.at[1], tick=F, line=-1, labels="Classification", hadj=0)
  text(x=rep(0.35, 7), y=0.5:6.5, pos=4, xpd=T,
       labels=c(as.vector(sapply(c("DS", "HI", "TS"),
                                 function(x){paste(x, " (", c("high", "low"),
                                                   " conf.)", sep="")})),
                "NS"))
}

# Plot histogram for color gradient
gradient.hist <- function(gradient, palette, x.ax.title="Gradient",
                          y.ax.title="\"# Genes\"", outline.lwd=1,
                          parmar=c(1, 1, 0.1, 0.1)){
  h <- hist(gradient, breaks=seq(0, 100, 1), plot=F)
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(0, 100), ylim=c(0, max(h$counts)),
       xaxt="n", xlab="", yaxt="n", ylab="", yaxs="i")
  rect(xleft=0:99, xright=1:100, ybottom=0, ytop=h$counts,
       col=palette, border=palette, lwd=0.3, xpd=T)
  segments(x0=h$breaks, x1=h$breaks,
           y0=c(0, h$counts), y1=c(h$counts, 0),
           col=blueblack, xpd=T, lwd=outline.lwd)
  segments(x0=h$breaks[-length(h$breaks)], x1=h$breaks[-1],
           y0=h$counts, y1=h$counts, col=blueblack, xpd=T, lwd=outline.lwd)
  axis(1, at=c(-100, 200), tck=0, labels=NA, col=blueblack)
  mtext(1, text=parse(text=x.ax.title), line=0.15)
  axis(2, at=c(-100, 2*par("usr")[4]), tck=0, labels=NA, col=blueblack)
  mtext(2, text=parse(text=y.ax.title))
}

# Plot legend for color gradient
gradient.legend <- function(palette){
  par(mar=rep(0.1, 4), bty="n")
  n.bins <- length(palette)
  plot(NA, xlim=c(0, n.bins), ylim=c(0, 1),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  rect(xleft=0:(n.bins-1), xright=1:n.bins,
       ybottom=0, ytop=1, border=palette, lwd=0.5, col=palette)
  rect(xleft=0, xright=n.bins, ybottom=0, ytop=1,
       col=NA, border=blueblack, xpd=T)
  box(bty="o", col=blueblack, xpd=T)
}

# Wrapper to generate all gradient plots
plot.gradients <- function(scores, gradient.norm, gradient.pal, null.x=NA,
                           null.y=NA, hist.x.ax.title, sub.out.prefix){
  # Scatterplot
  pt.colors.gradient <- gradient.pal[gradient.norm + 1]
  pdf(paste(out.prefix, "gene_scores_scatterplot",
            sub.out.prefix,"pdf", sep="."),
      height=2, width=2)
  scores.scatterplot(scores, pt.colors.gradient, add.cor=F,
                     category.lines.col=NA, margin.dens.height=0, pt.cex=0.125,
                     x.ax.title="pHaplo", y.ax.title="pTriplo",
                     ax.at=seq(0, 1, 0.25), parmar=c(2.7, 2.7, 0.4, 0.4))
  if(!is.na(null.x) & !is.na(null.y)){
    null.pt.idxs <- which(scores$pHaplo<null.x & scores$pTriplo<null.y)
    points(x=scores$pHaplo[null.pt.idxs], y=scores$pTriplo[null.pt.idxs],
           cex=0.125, col="gray95")
    segments(x0=c(0, 0), x1=c(null.x, null.x),
             y0=c(0, null.y), y1=c(null.y, 0),
             lend="butt", lwd=2, col="gray90")
    segments(x0=c(0, null.x), x1=rep(null.x, 2),
             y0=c(null.y, 0), y1=rep(null.y, 2),
             lty=2, col=blueblack, lend="round")
    box(bty="o", col=blueblack, xpd=T)
  }
  dev.off()

  # Histogram
  pdf(paste(out.prefix, "gene_scores_scatterplot", sub.out.prefix,
            "hist.pdf", sep="."),
      height=1.2, width=1.8)
  if(!is.na(null.x) & !is.na(null.y)){
    gradient.hist(gradient.norm[-null.pt.idxs], gradient.pal,
                  x.ax.title=hist.x.ax.title)
  }else{
    gradient.hist(gradient.norm, gradient.pal, x.ax.title=hist.x.ax.title)
  }

  dev.off()

  # Legend
  pdf(paste(out.prefix, "gene_scores_scatterplot", sub.out.prefix, "legend.pdf", sep="."),
      height=0.2, width=1.5)
  gradient.legend(gradient.pal)
  dev.off()
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option("--score-cutoffs", metavar="tsv", default=NULL,
              help=".tsv of gene score cutoffs generated by empirical_score_cutoffs.R")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: scores.tsv and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
out.prefix <- args$args[2]
score.cutoffs.in <- opts$`score-cutoffs`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# out.prefix <- "~/scratch/test_gene_score_feature_regressions"
# score.cutoffs.in <- "~/scratch/rCNV2.gene_score_cutoff.test.tsv"

# Load scores & classify genes into subgroups based on scores
scores <- load.scores(scores.in)
load.score.cutoffs(score.cutoffs.in)
ds.groups <- classify.genes.by.score(scores, hi.hc=hi.hc, hi.lc=hi.lc,
                                     ts.hc=ts.hc, ts.lc=ts.lc)
pt.colors.groupwise <- sapply(scores$gene, get.gene.color.byscore, ds.groups)

# Plot scores colored by category
pdf(paste(out.prefix, "gene_scores_scatterplot.pdf", sep="."),
    height=3.25, width=3.25)
scores.scatterplot(scores, pt.colors.groupwise, margin.dens.height=0.1,
                   hi.hc=hi.hc, hi.lc=hi.lc, ts.hc=ts.hc, ts.lc=ts.lc,
                   pt.cex=0.15, ax.tick=-0.02, parmar=c(2.7, 2.7, 1.5, 1.5))
dev.off()

# Plot legend for category-based coloring
pdf(paste(out.prefix, "gene_scores_scatterplot.legend.pdf", sep="."),
    height=2, width=2)
plot.scores.scatter.legend(scores, ds.groups)
dev.off()

# Plot scores without categories and colored by density
# Optional: code below can be used to color by mixture of x,y coordinates
# color.grid <- make.color.grid()
# pt.colors.xybasis <- sapply(scores$gene, query.color.grid, scores, color.grid)
# Color by density
pt.colors.density.df.unsorted <- color.points.by.density(scores$pHaplo, scores$pTriplo)
pt.colors.density <- pt.colors.density.df.unsorted$col[order(as.numeric(rownames(pt.colors.density.df.unsorted)))]
pdf(paste(out.prefix, "gene_scores_scatterplot.no_categories.pdf", sep="."),
    height=3.25, width=3.25)
scores.scatterplot(scores, pt.colors.density,
                   margin.dens.height=0.1, margin.dens.gradient=TRUE,
                   bg.col="white", gridlines.col=NA, category.lines.col=NA,
                   pt.cex=0.15, ax.tick=-0.02, parmar=c(2.7, 2.7, 1.5, 1.5))
dev.off()

# Plot scores colored by DS gradient
ds.gradient <- apply(scores[, -1], 1, function(vals){
  min(as.numeric(vals))
})
ds.gradient.norm <- round(100*ds.gradient)
plot.gradients(scores, ds.gradient.norm, ds.gradient.pal,
               hist.x.ax.title="italic(min) * \"(pHaplo, pTriplo)\"",
               sub.out.prefix="ds_gradient")

# Plot scores colored by HI/TS gradient
hits.gradient <- scores$pTriplo - scores$pHaplo
hits.gradient.norm <- round(((100*hits.gradient) + 100) / 2)
plot.gradients(scores, hits.gradient.norm, hits.gradient.pal,
               hist.x.ax.title="\"pTriplo\" - \"pHaplo\"",
               sub.out.prefix="hi_ts_gradient")
