#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot CDS optimization analyses for rCNV supplement


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load precomputed CDS optimization data and determine optimal cutoff
load.opt.stats <- function(stats.in, method="both"){
  # Read precomputed stats
  opt.df <- read.table(stats.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(opt.df)[1] <- gsub("#", "", colnames(opt.df)[1], fixed=T)

  # Find optimal cutoff
  best.idx <- which(opt.df$pareto.dist == min(opt.df$pareto.dist, na.rm=T))

  # Return data for plotting
  return(list("data" = opt.df, "best.idx" = best.idx))
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot optimization data
plot.opt <- function(opt.res, power.col="#00A4CCFF", effectsize.col="#F95700FF"){
  # Get plot data
  opt.df <- opt.res$data
  best.idx <- opt.res$best.idx
  best.cds <- opt.df$min.cds[best.idx]

  # Prep plot area
  par(bty="n", mar=c(2.5, 2.7, 1, 2.7))
  plot(NA, xlim=c(0, 1), ylim=c(0, 1),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  rect(xleft=best.cds-0.005, xright=best.cds+0.005,
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=adjustcolor(highlight.color, alpha=0.3),
       border=NA, bty="n", xpd=T)
  points(x=best.cds, y=par("usr")[4]+(0.05*diff(par("usr")[3:4])), xpd=T, pch=25,
         col=highlight.color, bg=highlight.color, cex=1.2)
  text(x=best.cds, y=par("usr")[4]+(0.065*diff(par("usr")[3:4])), xpd=T,
       labels=best.cds, cex=5/6)

  # Plot power & effect size
  points(opt.df$min.cds, opt.df$or.frac, lwd=4, col=effectsize.col, type="l", xpd=T)
  points(opt.df$min.cds, opt.df$power.frac, lwd=4, col=power.col, type="l", xpd=T)

  # Add axes
  axis(1, labels=NA, tck=-0.02)
  sapply(axTicks(1), function(x){
    axis(1, at=x, tick=F, line=-0.75, cex.axis=5/6)
  })
  mtext(1, text="Minimum CDS Cutoff", line=1.25)
  axis(2, labels=NA, tck=-0.02, col=power.col)
  axis(2, tick=F, line=-0.65, cex.axis=0.85, las=2, col.axis=power.col)
  mtext(2, text="Power Retained", line=1.5, col=power.col, cex.axis=5.5/6)
  axis(4, labels=NA, tck=-0.02, col=effectsize.col)
  axis(4, tick=F, line=-0.65, cex.axis=5/6, las=2, col.axis=effectsize.col)
  mtext(4, text="Effect Size Retained", line=1.5, col=effectsize.col, cex.axis=5.5/6)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog DEL.stats.tsv DUP.stats.tsv out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to vars
del.stats.in <- args$args[1]
dup.stats.in <- args$args[2]
out.prefix <- args$args[3]

# Load precomputed optimization stats
opt.stats <- list("DEL" = load.opt.stats(del.stats.in),
                  "DUP" = load.opt.stats(dup.stats.in))

# Plot optimization data
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(out.prefix, "cds_opt", cnv, "pdf", sep="."),
      height=2.3, width=2.6)
  plot.opt(opt.stats[[cnv]], control.cnv.colors[cnv], cnv.colors[cnv])
  dev.off()
})
