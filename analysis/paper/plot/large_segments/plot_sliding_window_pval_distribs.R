#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot mini QQs per phenotype and primary vs secondary P-values from sliding window analysis for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000, check.names=F)


######################
### DATA FUNCTIONS ###
######################
# Load a matrix of p-values
load.pval.matrix <- function(matrix.in, has.coords=T, p.is.phred=T){
  x <- read.table(matrix.in, header=T, sep="\t", comment.char="")
  if(has.coords == T){
    colnames(x)[1:3] <- c("chrom", "start", "end")
    x[, 1] <- as.character(x[, 1])
    x[, -1] <- apply(x[, -1], 2, as.numeric)
    coords <- as.data.frame(x[, 1:3])
    pvals <- as.data.frame(x[, -c(1:3)])
  }else{
    coords <- NULL
    pvals <- as.data.frame(apply(x, 2, as.numeric))
  }
  if(p.is.phred == T){
    pvals <- as.data.frame(apply(pvals, 2, function(x){10^-x}))
  }
  expected <- ppoints(nrow(pvals))
  lambdas <- apply(pvals, 2, function(obs){dchisq(median(obs, na.rm=T), df=1)/dchisq(median(expected), df=1)})
  return(list("coords" = coords, "pvals" = pvals,
              "expected" = expected, "lambdas" = lambdas))
}

# Load HPO sample sizes
load.hpo.samplesize <- function(hpo.samplesize.in, hpos=NULL){
  hpo.n <- read.table(hpo.samplesize.in, header=T, sep="\t", comment.char="")
  out.v <- as.numeric(hpo.n$Total)
  names(out.v) <- hpo.n[, 1]
  if(!is.null(hpos)){
    out.v <- out.v[hpos]
  }
  return(out.v)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Mini quantile-quantile plot for deletions & duplications for a single phenotype
mini.qq <- function (del, dup, hpo, ymax=NULL, 
                     pt.cex=0.1, title.cex=1, axis.cex=1,
                     parmar=c(1.85, 1.85, 1, 1)){
  # Reformat HPO
  nocolon <- gsub(":", "", hpo, fixed=T)
  
  # Get plotting data
  p.del <- as.numeric(del$pvals[, which(colnames(del$pvals) == paste(nocolon, "DEL", sep="_"))])
  exp.del <- as.numeric(del$expected)
  p.dup <- as.numeric(dup$pvals[, which(colnames(dup$pvals) == paste(nocolon, "DUP", sep="_"))])
  exp.dup <- as.numeric(dup$expected)
  plot.dat <- list("DEL" = data.frame("exp" = -log10(sort(exp.del, decreasing=T)),
                                      "obs" = -log10(sort.int(p.del, decreasing=T, na.last=F))),
                   "DUP" = data.frame("exp" = -log10(sort(exp.dup, decreasing=T)),
                                      "obs" = -log10(sort.int(p.dup, decreasing=T, na.last=F))))
  xmax <- max(plot.dat$DEL$exp, plot.dat$DUP$exp, na.rm=T)
  if(is.null(ymax)){
    ymax <- max(max(plot.dat$DEL[which(!is.infinite(plot.dat$DEL$obs)), ], na.rm=T), 
                max(plot.dat$DUP[which(!is.infinite(plot.dat$DUP$obs)), ], na.rm=T))
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, xmax), ylim=c(0, ymax),
       xaxt="n", xaxs="i", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=bluewhite, border=NA, bty="n")
  x.ax.at <- seq(0, max(par("usr")[2]), by=ceiling(max(par("usr")[2]) / 4))
  y.ax.at <- seq(0, max(par("usr")[4]), by=ceiling(max(par("usr")[4]) / 4))
  # abline(v=x.ax.at, h=y.ax.at, col="white")
  abline(0, 1, col=blueblack)

  # Add points
  points(plot.dat$DUP, pch=19, cex=pt.cex, col=cnv.colors[2])
  points(plot.dat$DEL, pch=19, cex=pt.cex, col=cnv.colors[1])
  
  # Add title & axes
  mtext(3, line=0.1, text=hpo, cex=title.cex, font=2)
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=x.ax.at, tck=-0.05, labels=NA, col=blueblack)
  sapply(x.ax.at, function(x){
    axis(1, at=x, tick=F, line=-1.05, cex.axis=axis.cex)
  })
  mtext(1, line=0.7, text=bquote("-log"[10] * "(" * italic("P")["Exp"] * ")"), cex=title.cex)
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, tck=-0.05, labels=NA, col=blueblack)
  sapply(y.ax.at, function(y){
    axis(2, at=y, tick=F, line=-0.65, las=2, cex.axis=axis.cex)
  })
  mtext(2, line=0.9, text=bquote("-log"[10] * "(" * italic("P")["Obs"] * ")"), cex=title.cex)
}

# Horizontal dotplot of lambdas
plot.lambdas <- function(del, dup, hpos,
                         parmar=c(0.25, 10, 2.5, 0.25)){
  # Get plot data
  plot.dat <- as.data.frame(t(sapply(hpos, function(hpo){
    nocolon <- gsub(":", "", hpo, fixed=T)
    lambda.del <- del$lambdas[paste(nocolon, "DEL", sep="_")]
    lambda.dup <- dup$lambdas[paste(nocolon, "DUP", sep="_")]
    c(hpo, hpo.abbrevs[hpo], lambda.del, lambda.dup)
  })))
  plot.dat[, 3:4] <- apply(plot.dat[, 3:4], 2, as.numeric)
  colnames(plot.dat) <- c("hpo", "abbrev", "del", "dup")
  plot.dat <- plot.dat[order(apply(plot.dat[, 3:4], 1, min, na.rm=T)), ]
  xlims <- range(plot.dat[, 3:4], na.rm=T)
  y.at <- (1:nrow(plot.dat)) - 0.5
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=c(nrow(plot.dat), 0),
       xaxt="n", xlab="", yaxt="n", ylab="", yaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=y.at-0.35, ytop=y.at+0.35, border=NA, bty="n", col=bluewhite)
  abline(v=axTicks(3), col="white")
  abline(v=1, col=blueblack)
  
  # Add points
  points(x=plot.dat$del, y=y.at-0.1, pch=21, col=cnv.blacks[1], bg=cnv.colors[1])
  points(x=plot.dat$dup, y=y.at+0.1, pch=21, col=cnv.blacks[2], bg=cnv.colors[2])
  
  # Add axes
  y.labels <- sapply(1:nrow(plot.dat), function(i){
    paste(plot.dat$abbrev[i], " (", plot.dat$hpo[i], ")", sep="")
  })
  axis(2, at=y.at, line=-0.8, las=2, tick=F, labels=y.labels)
  axis(3, at=c(-10e10, 10e10), tck=0, labels=NA)
  axis(3, tck=-0.025, labels=NA)
  axis(3, tick=F, line=-0.65)
  mtext(3, line=1.2, text=bquote("Genomic Inflation Parameter," ~ lambda))
}

# Scatterplot of lambdas vs sample size
lambda.scatter <- function(del, dup, hpo.n,
                           parmar=c(2.5, 3, 0.25, 0.25)){
  # Get plot data
  plot.dat <- as.data.frame(t(sapply(hpos, function(hpo){
    nocolon <- gsub(":", "", hpo, fixed=T)
    lambda.del <- del$lambdas[paste(nocolon, "DEL", sep="_")]
    lambda.dup <- dup$lambdas[paste(nocolon, "DUP", sep="_")]
    c(hpo, hpo.n[hpo], lambda.del, lambda.dup)
  })))
  plot.dat[, 2:4] <- apply(plot.dat[, 2:4], 2, as.numeric)
  colnames(plot.dat) <- c("hpo", "n", "del", "dup")
  # xlims <- range(log10(plot.dat$n), na.rm=T)
  xlims <- log10(c(1000, 500000))
  ylims <- range(plot.dat[, 3:4], na.rm=T)
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       border=NA, bty="n", col=bluewhite)
  x.ax.at <- log10(logscale.major)
  y.ax.at <- seq(0, 2, 0.1)
  abline(h=y.ax.at, v=x.ax.at, col="white")
  abline(h=1, col=blueblack)
  
  # Add linear fits
  del.fit <- robust.lm(log10(plot.dat$n), plot.dat$del)
  del.fit.p <- format.pval(cor.test(log10(plot.dat$n), plot.dat$del)$p.value)
  polygon(x=c(del.fit$ci$x, rev(del.fit$ci$x)),
          y=c(del.fit$ci$lower, rev(del.fit$ci$upper)),
          border=NA, bty="n", col=adjustcolor(cnv.colors[1], alpha=0.15))
  abline(del.fit$fit, lwd=2, col=cnv.colors[1])
  dup.fit <- robust.lm(log10(plot.dat$n), plot.dat$dup)
  dup.fit.p <- format.pval(cor.test(log10(plot.dat$n), plot.dat$dup)$p.value)
  polygon(x=c(dup.fit$ci$x, rev(dup.fit$ci$x)),
          y=c(dup.fit$ci$lower, rev(dup.fit$ci$upper)),
          border=NA, bty="n", col=adjustcolor(cnv.colors[2], alpha=0.15))
  abline(dup.fit$fit, lwd=2, col=cnv.colors[2])
  
  # Add points
  sapply(1:nrow(plot.dat), function(i){
    points(x=log10(plot.dat$n)[i], y=plot.dat$dup[i], pch=21, 
           bg=cnv.colors[2], col=cnv.blacks[2])
    points(x=log10(plot.dat$n)[i], y=plot.dat$del[i], pch=21, 
           bg=cnv.colors[1], col=cnv.blacks[1])
  })
  
  # Add axes
  axis(1, at=log10(logscale.minor), tck=-0.0125, labels=NA, col=blueblack)
  axis(1, at=x.ax.at, tck=-0.025, labels=NA, col=blueblack)
  axis(1, at=x.ax.at, tick=F, line=-0.7, labels=logscale.major/1000)
  mtext(1, line=1.25, text="Cases (Thousands)")
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, tck=-0.025, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, tick=F, line=-0.65, las=2)
  mtext(2, line=1.75, text=bquote("Genomic Inflation," ~ lambda))
  
  # Add P-values
  x.length <- par("usr")[2] - par("usr")[1]
  y.length <- ylims[2] - ylims[1]
  text(x=par("usr")[2] + 0.05 * x.length, 
       y=par("usr")[3] + 0.175 * y.length,
       pos=2, col=cnv.colors[1], labels=bquote(.(del.fit.p)))
  text(x=par("usr")[2] + 0.05 * x.length, 
       y=par("usr")[3] + 0.05 * y.length,
       pos=2, col=cnv.colors[2], labels=bquote(.(dup.fit.p)))
}

# Primary vs. secondary P-value scatterplot
primary.vs.secondary.scatter <- function(pvals.1, pvals.2, keep.idx=NULL,
                                         cutoff.1=6, cutoff.2=-log10(0.05),
                                         sig.color="black", nonsig.color="gray75",
                                         pt.cex=0.2, parmar=c(2.25, 2.25, 0.25, 0.25)){
  # Subset p-values to indexes to keep (if optioned)
  if(!is.null(keep.idx)){
    pvals.1 <- pvals.1[keep.idx, ]
    pvals.2 <- pvals.2[keep.idx, ]
  }
  
  # Gather plot data
  pvals <- data.frame("primary"=-log10(as.vector(as.matrix(pvals.1))),
                      "secondary"=-log10(as.vector(as.matrix(pvals.2))))
  pvals <- pvals[which(!is.na(pvals$primary) & !is.na(pvals$secondary) &
                         !is.infinite(pvals$primary) & !is.infinite(pvals$secondary)), ]
  ax.lims <- range(pvals, na.rm=T)
  pt.colors <- apply(pvals, 1, function(vals){
    vals <- as.numeric(vals)
    if(vals[1] >= cutoff.1 & vals[2] >= cutoff.2){
      sig.color
    }else{
      nonsig.color
    }
  })
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=ax.lims, ylim=ax.lims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       border=NA, bty="n", col=bluewhite)
  abline(h=axTicks(2), v=axTicks(1), col="white")
  # abline(0, 1, col=blueblack)
  
  # Add points
  points(pvals, pch=19, cex=pt.cex, col=pt.colors)
  
  # Add annotation lines
  abline(v=cutoff.1, h=cutoff.2, lty=2, col=blueblack)
  
  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, tck=-0.03, labels=NA, col=blueblack)
  axis(1, tick=F, line=-0.75)
  mtext(1, line=1.4, text=bquote("-log"[10] * "(" * italic("P")["Primary"] * ")"))
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, tck=-0.03, labels=NA, col=blueblack)
  axis(2, tick=F, line=-0.65, las=2)
  mtext(2, line=1.25, text=bquote("-log"[10] * "(" * italic("P")["Secondary"] * ")"))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
  make_option(c("--del-cutoff"), default=10e-6, help="Significant P-value cutoff for deletions."),
  make_option(c("--dup-cutoff"), default=10e-6, help="Significant P-value cutoff for duplications.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog primary.del.pvals.bed primary.dup.pvals.bed",
                                            "secondary.del.pvals.bed secondary.dup.pvals.bed",
                                            "hpos.list hpo.samplesize.tsv output_prefix", 
                                            sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 7){
  stop(paste("Seven positional arguments required: primary.del.pvals.bed,",
             "primary.dup.pvals.bed, secondary.del.pvals.bed,",
             "secondary.dup.pvals.bed, hpos.list, hpo_samplesizes.tsv, and output_prefix\n"))
}

# Writes args & opts to vars
primary.del.pvals.in <- args$args[1]
primary.dup.pvals.in <- args$args[2]
secondary.del.pvals.in <- args$args[3]
secondary.dup.pvals.in <- args$args[4]
hpos.in <- args$args[5]
hpo.samplesize.in <- args$args[6]
out.prefix <- args$args[7]
rcnv.config <- opts$`rcnv-config`
del.cutoff <- -log10(as.numeric(opts$`del-cutoff`))
dup.cutoff <- -log10(as.numeric(opts$`dup-cutoff`))

# # DEV PARAMETERS
# primary.del.pvals.in <- "~/scratch/rCNV2_analysis_d1.DEL.meta_phred_p.all_hpos.bed.gz"
# primary.dup.pvals.in <- "~/scratch/rCNV2_analysis_d1.DUP.meta_phred_p.all_hpos.bed.gz"
# secondary.del.pvals.in <- "~/scratch/rCNV2_analysis_d1.DEL.meta_phred_p_secondary.all_hpos.bed.gz"
# secondary.dup.pvals.in <- "~/scratch/rCNV2_analysis_d1.DUP.meta_phred_p_secondary.all_hpos.bed.gz"
# hpos.in <- "~/scratch/rCNV2_analysis_d1.reordered_hpos.txt"
# hpo.samplesize.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/test_pval_qc"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"
# del.cutoff <- 6
# dup.cutoff <- 6

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load pvalue matrices & compute lambdas
del.1 <- load.pval.matrix(primary.del.pvals.in)
dup.1 <- load.pval.matrix(primary.dup.pvals.in)
del.2 <- load.pval.matrix(secondary.del.pvals.in)
dup.2 <- load.pval.matrix(secondary.dup.pvals.in)

# Read reordered HPOs
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
hpo.n <- load.hpo.samplesize(hpo.samplesize.in, hpos)

# Plot grid of QQs
n.plots.wide <- 10
n.plots.tall <- ceiling(length(hpos) / n.plots.wide)
dim.scalar <- 7.5 / 10
png(paste(out.prefix, "qq_grid_byHPO.png", sep="."),
    height=dim.scalar*300*n.plots.tall, width=dim.scalar*300*n.plots.wide,
    res=300)
par(mfrow=c(n.plots.tall, n.plots.wide))
for(hpo in hpos){
  mini.qq(del.1, dup.1, hpo, title.cex=0.5, axis.cex=0.7)
}
dev.off()

# Horizontal dotplot of primary lambdas
pdf(paste(out.prefix, "primary_lambdas_byHPO.pdf", sep="."),
    height=6, width=5)
plot.lambdas(del.1, dup.1, hpos, parmar=c(0.25, 14, 2.25, 0.75))
dev.off()

# Scatterplot of lambdas vs sample size
pdf(paste(out.prefix, "primary_lambdas_vs_sampleSize.pdf", sep="."),
    height=2.3, width=2.4)
lambda.scatter(del.1, dup.1, hpo.n)
dev.off()

# Scatterplots of primary vs. secondary P-values
png(paste(out.prefix, "primary_vs_secondary_pvalue.DEL.png", sep="."),
    height=2.25*300, width=2.25*300, res=300)
primary.vs.secondary.scatter(del.1$pvals, del.2$pvals, cutoff.1=del.cutoff,
                             sig.color=cnv.colors[1], nonsig.color=control.cnv.colors[1],
                             pt.cex=0.175, parmar=c(2.5, 2.5, 0.15, 0.15))
dev.off()
png(paste(out.prefix, "primary_vs_secondary_pvalue.DUP.png", sep="."),
    height=2.25*300, width=2.25*300, res=300)
primary.vs.secondary.scatter(dup.1$pvals, dup.2$pvals, cutoff.1=dup.cutoff,
                             sig.color=cnv.colors[2], nonsig.color=control.cnv.colors[2],
                             pt.cex=0.175, parmar=c(2.5, 2.5, 0.15, 0.15))
dev.off()

