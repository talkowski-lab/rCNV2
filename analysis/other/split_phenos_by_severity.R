#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Bifurcate HPOs by effect size of whole-gene deletions of constrained genes


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=T)


# Plot lnORs with cutoff & classification
plot.dat <- function(dat, high.hpos, low.hpos, threshold){
  # Get plot data
  colors <- sapply(dat$hpo, get.hpo.color, color.by="severity")
  ylims <- range(c(dat$lnOR_lower, dat$lnOR_upper), na.rm=T)
  y.ax.labs <- c(0, 2^(0:10))
  y.ax.at <- log(y.ax.labs)
  x.at <- (1:nrow(dat)) - 0.5
  half.idx <- ceiling(nrow(dat)/2)
  y.buffer <- 0.035*(ylims[2]-ylims[1])

  # Prep plot area
  par(mar=c(3, 3.5, 1.5, 0.5), bty="n")
  plot(x=c(0, nrow(dat)), y=ylims, type="n",
       xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=c(log(threshold), par("usr")[3]),
       ytop=c(par("usr")[4], log(threshold)),
       bty="n", border=NA,
       col=sapply(severity.colors[c("developmental", "adult")],
                  adjustcolor, alpha=0.2))
  abline(h=c(0, log(threshold)), lty=c(1, 5))
  text(x=rep(par("usr")[2], 2), y=log(threshold)+(c(-1, 1)*y.buffer),
       labels=c("Adult-Onset", "Developmental"), font=2, pos=2,
       col=severity.colors[3:2])
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA)
  mtext(1, line=0, padj=1,
        text="HPOs ordered by odds ratio\nof constrained gene deletions")
  axis(2, at=c(-10e10, y.ax.at), labels=NA)
  axis(2, at=y.ax.at, labels=y.ax.labs, tick=F, las=2, line=-0.65)
  mtext(2, line=1.5, text="Odds ratio of constrained gene deletions")
  mtext(3, line=0.1, font=2, text="HPOs Split by Effect Size of Constrained Gene Deletions")

  # Add points, 95% CIs, and labels
  segments(x0=x.at, x1=x.at, y0=dat$lnOR_lower, y1=dat$lnOR_upper, col=colors, lwd=2)
  points(x=x.at, y=dat$lnOR, pch=21, bg=colors)
  text(x=x.at[1:half.idx], y=dat$lnOR[1:half.idx], pos=4, cex=0.75,
       labels=dat$description[1:half.idx])
  text(x=x.at[(half.idx+1):nrow(dat)], y=dat$lnOR[(half.idx+1):nrow(dat)],
       pos=2, cex=0.75, labels=dat$description[(half.idx+1):nrow(dat)])
}


# Get command-line arguments & options
option_list <- list(
  make_option(c("-l", "--use-lower"), action="store_true", default=FALSE,
              help="use lower bound of 95% confidence interval [default: use point estimate]"),
  make_option(c("-p", "--out-prefix"), type="character", default="rCNV2_pheno_split",
              help="prefix for output files [default %default]", metavar="string"),
  make_option(c("-t", "--threshold"), type="numeric", default=2,
              help="odds ratio threshold for splitting HPOs [default %default]",
              metavar="numeric")
)
args <- parse_args(OptionParser(usage="%prog data.tsv",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options


# Writes args & opts to vars
dat.in <- args$args[1]
use.lower <- opts$`use-lower`
out.prefix <- opts$`out-prefix`
threshold <- opts$`threshold`


# Read data
dat <- read.table(dat.in, header=T, sep="\t", comment.char="")
dat$description <- hpo.abbrevs[dat$hpo]


# Split HPOs by OR
if(use.lower==TRUE){
  high.hpos <- dat$hpo[which(!is.na(dat$lnOR) & dat$lnOR_lower>=log(threshold))]
  low.hpos <- dat$hpo[which(!is.na(dat$lnOR) & dat$lnOR_lower<log(threshold))]
}else{
  high.hpos <- dat$hpo[which(!is.na(dat$lnOR) & dat$lnOR>=log(threshold))]
  low.hpos <- dat$hpo[which(!is.na(dat$lnOR) & dat$lnOR<log(threshold))]
}


# Write classifications to file
write.table(high.hpos, paste(out.prefix, "developmental.list", sep="."),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(low.hpos, paste(out.prefix, "adult.list", sep="."),
            col.names=F, row.names=F, sep="\t", quote=F)
all.df <- data.frame("HPO"=c(high.hpos, low.hpos),
                     "group"=c(rep("developmental", length(high.hpos)),
                               rep("adult", length(low.hpos))))
colnames(all.df)[1] <- "#HPO"
write.table(all.df, paste(out.prefix, "all_hpos.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

# Plot classification
pdf(paste(out.prefix, "or_distribs.pdf", sep="."),
    height=4.5, width=8)
plot.dat(dat, high.hpos, low.hpos, threshold)
dev.off()
