#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot pHI/pTS score scatterplot for gene score analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000, family="sans")


##########################
### PLOTTING FUNCTIONS ###
##########################
# Helper function to add marginal density plot
marginal.density <- function(vals, colors, scale=0.25, bw=0.02, rotate=F,
                             hc.cutoff=0.9, lc.cutoff=0.5){
  # Gather plot values
  d <- density(vals, bw=bw)
  idx <- d$x
  idx[which(idx<0)] <- 0; idx[which(idx>1)] <- 1
  h <- d$y
  h <- scale*(1/max(h, na.rm=T))*h
  h <- h + 1
  
  # Split indexes into bins
  ns.idx <- which(idx<lc.cutoff)
  lc.idx <- which(idx>=lc.cutoff & idx<hc.cutoff)
  hc.idx <- which(idx>=hc.cutoff)
  
  # Compute polygon coordinate vectors based on rotation
  if(rotate==T){
    x <- h; y <- idx
    ns.df <- data.frame("x"=c(x[ns.idx], rep(1, length(ns.idx))),
                        "y"=c(y[ns.idx], rev(y[ns.idx])))
    lc.df <- data.frame("x"=c(x[lc.idx], rep(1, length(lc.idx))),
                        "y"=c(y[lc.idx], rev(y[lc.idx])))
    hc.df <- data.frame("x"=c(x[hc.idx], rep(1, length(hc.idx))),
                        "y"=c(y[hc.idx], rev(y[hc.idx])))
    all.df <- data.frame("x"=c(x, rep(1, length(x))), "y"=c(y, rev(y)))
  }else{
    x <- idx; y <- h
    ns.df <- data.frame("x"=c(x[ns.idx], rev(x[ns.idx])),
                        "y"=c(y[ns.idx], rep(1, length(ns.idx))))
    lc.df <- data.frame("x"=c(x[lc.idx], rev(x[lc.idx])),
                        "y"=c(y[lc.idx], rep(1, length(lc.idx))))
    hc.df <- data.frame("x"=c(x[hc.idx], rev(x[hc.idx])),
                        "y"=c(y[hc.idx], rep(1, length(hc.idx))))
    all.df <- data.frame("x"=c(x, rev(x)), "y"=c(y, rep(1, length(y))))
  }
  
  # Plot polygons in ascending order
  polygon(x=ns.df$x, y=ns.df$y, border=colors[3], col=colors[3], xpd=T)
  polygon(x=lc.df$x, y=lc.df$y, border=colors[2], col=colors[2], xpd=T)
  polygon(x=hc.df$x, y=hc.df$y, border=colors[1], col=colors[1], xpd=T)
  polygon(x=all.df$x, y=all.df$y, col=NA, border=colors[1], bty="n", xpd=T)
}

# Scatterplot of genes by scores, colored by group
scores.scatterplot <- function(scores, ds.groups,
                               hc.cutoff=0.9, lc.cutoff=0.5,
                               margin.dens.height=0.175,
                               parmar=c(2.7, 2.7, 1, 1)){
  # Get plot.data
  pt.colors <- sapply(scores$gene, get.gene.color.byscore, ds.groups)
  ax.at <- seq(0, 1, 0.2)

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 1 + margin.dens.height), ylim=c(0, 1 + margin.dens.height),
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  rect(xleft=par("usr")[1], xright=1,
       ybottom=par("usr")[3], ytop=1,
       border=NA, bty="n", col=bluewhite)
  segments(x0=ax.at, x1=ax.at,
           y0=par("usr")[3], y1=1, col="white")
  segments(x0=par("usr")[1], x1=1,
           y0=ax.at, y1=ax.at, col="white")
  
  # Add points
  points(scores$pHI, scores$pTS, cex=0.1, col=pt.colors)
  
  # Add marginal densities
  marginal.density(scores$pHI, colors=c(cnv.colors[1], control.cnv.colors[1], redwhite),
                   bw=0.01, scale=margin.dens.height)
  marginal.density(scores$pTS, colors=c(cnv.colors[2], control.cnv.colors[2], bluewhite), 
                   bw=0.01, scale=margin.dens.height, rotate=T)
  
  # Add delimiting lines
  segments(x0=c(hc.cutoff, lc.cutoff, hc.cutoff, lc.cutoff),
           x1=c(hc.cutoff, lc.cutoff, hc.cutoff, lc.cutoff),
           y0=c(hc.cutoff, lc.cutoff, rep(par("usr")[3], 2)),
           y1=c(1, 1, rep(lc.cutoff, 2)),
           lwd=1, lend="round", col=blueblack, lty=3)
  segments(y0=c(hc.cutoff, lc.cutoff, hc.cutoff, lc.cutoff),
           y1=c(hc.cutoff, lc.cutoff, hc.cutoff, lc.cutoff),
           x0=c(hc.cutoff, lc.cutoff, rep(par("usr")[3], 2)),
           x1=c(1, 1, rep(lc.cutoff, 2)),
           lwd=1, lend="round", col=blueblack, lty=3)
  
  # Add axes
  axis(1, at=ax.at, col=blueblack, labels=NA, tck=-0.03)
  axis(2, at=ax.at, col=blueblack, labels=NA, tck=-0.03)
  sapply(ax.at, function(i){
    axis(1, at=i, col=blueblack, tick=F, line=-0.65)
    axis(2, at=i, col=blueblack, tick=F, line=-0.65, las=2)
  })
  axis(1, at=0.5, tick=F, line=0.6, labels=bquote(italic(P) * "(Haploinsufficient [HI])"))
  axis(2, at=0.5, tick=F, line=0.6, labels=bquote(italic(P) * "(Triplosensitive [TS])"))
  
  # Add correlation coefficient
  r2 <- cor(scores$pHI, scores$pTS)^2
  r2.fmt <- formatC(round(r2, 2), small.interval=2)
  text(x=1.05+(margin.dens.height/2), 
       y=1.065+(margin.dens.height/2), 
       xpd=T, cex=0.9, srt=45,
       labels=bquote(italic(R)^2 * "=" * .(r2.fmt)))
  
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
  par(mar=c(0.1, 0.1, 1.25, 0.1), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(7, 0),
       xaxt="n", xaxs="i", xlab="", yaxt="n", yaxs="i", ylab="")
  
  # Add points
  points(x=rep(0.05, 7), y=0.5:6.5, pch=19, col=colors, cex=2.25)
  
  # Add N genes
  gene.ax.at <- c(0.125, 0.3)
  axis(3, at=gene.ax.at, tck=0, labels=NA, col=blueblack)
  axis(3, at=gene.ax.at[2], tick=F, line=-1, labels="Genes", hadj=1)
  sapply(1:length(ds.groups), function(i){
    text(x=0.35, y=i-0.5, pos=2, labels=prettyNum(length(ds.groups[[i]]), big.mark=","))
  })
  
  # Add labels
  cat.ax.at <- c(0.35, 0.95)
  axis(3, at=cat.ax.at, tck=0, labels=NA, col=blueblack)
  axis(3, at=cat.ax.at[1], tick=F, line=-1, labels="Classification", hadj=0)
  text(x=rep(0.3, 7), y=0.5:6.5, pos=4, 
       labels=c(as.vector(sapply(c("DS", "HI", "TS"), 
                                 function(x){paste(x, " (", c("high", "low"), " confidence)", sep="")})),
                "Not sensitive [NS]"))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Three positional arguments required: scores.tsv and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
out.prefix <- args$args[2]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# out.prefix <- "~/scratch/test_gene_score_feature_regressions"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load scores & classify genes into subgroups based on scores
scores <- load.scores(scores.in)
ds.groups <- classify.genes(scores, hc.cutoff=0.9, lc.cutoff=0.5)

# Plot scores
png(paste(out.prefix, "gene_scores_scatterplot.png", sep="."),
    height=3.25*300, width=3.25*300, res=300, bg=NA, family="sans")
scores.scatterplot(scores, ds.groups, margin.dens.height=0.1,
                   parmar=c(2.7, 2.7, 1.5, 1.5))
dev.off()

# Plot legend
pdf(paste(out.prefix, "gene_scores_scatterplot.legend.pdf", sep="."),
    height=2.6, width=2.4, family="sans")
plot.scores.scatter.legend(scores, ds.groups)
dev.off()
