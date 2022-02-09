#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot mini QQs per phenotype and primary vs secondary P-values from sliding window analysis for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000, check.names=F)


######################
### DATA FUNCTIONS ###
######################
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
mini.qq <- function (del, dup, hpo, ymax=NULL, blue.bg=TRUE,
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
  if(blue.bg==TRUE){
    plot.bg <- bluewhite
    plot.border <- NA
    plot.bty <- "n"
  }else{
    plot.bg <- "white"
    plot.border <- NA
    plot.bty <- "n"
  }

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, xmax), ylim=c(0, ymax),
       xaxt="n", xaxs="i", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=plot.bg, border=plot.border, bty=plot.bty)
  x.ax.at <- seq(0, max(par("usr")[2]), by=ceiling(max(par("usr")[2]) / 4))
  y.ax.at <- seq(0, max(par("usr")[4]), by=ceiling(max(par("usr")[4]) / 4))
  # abline(v=x.ax.at, h=y.ax.at, col="white")
  abline(0, 1, col=blueblack)

  # Add points
  points(plot.dat$DUP, pch=19, cex=pt.cex, col=cnv.colors[2])
  points(plot.dat$DEL, pch=19, cex=pt.cex, col=cnv.colors[1])

  # Add title & axes
  mtext(3, line=0.1, text=hpo.abbrevs[hpo], cex=title.cex, xpd=T)
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
plot.lambdas <- function(del, dup, hpos, hpo.cex=0.75,
                         parmar=c(0.25, 10, 2.5, 10)){
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
  max.x.diff <- max(abs(1-plot.dat[, 3:4]), na.rm=T)
  xlims <- c(1-max.x.diff, 1+max.x.diff)
  y.at <- (1:nrow(plot.dat)) - 0.5
  y.evens <- seq(1, nrow(plot.dat), 2)
  y.odds <- seq(2, nrow(plot.dat), 2)

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=c(nrow(plot.dat), 0),
       xaxt="n", xlab="", yaxt="n", ylab="", yaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=y.at-0.35, ytop=y.at+0.35, border=NA, bty="n", col=bluewhite)
  # abline(v=axTicks(3), col="white")
  abline(v=1, col=blueblack)

  # Add points
  points(x=plot.dat$del, y=y.at-0.1, pch=21, col=cnv.blacks[1], bg=cnv.colors[1])
  points(x=plot.dat$dup, y=y.at+0.1, pch=21, col=cnv.blacks[2], bg=cnv.colors[2])

  # Add axes
  sapply(1:nrow(plot.dat), function(i){
    y.label <- paste(plot.dat$abbrev[i], " (", plot.dat$hpo[i], ")", sep="")
    if(i %% 2 > 0){
      side <- 2
    }else{
      side <- 4
    }
    axis(side, at=y.at[i], line=-0.8, las=2, tick=F, labels=y.label, cex.axis=hpo.cex)
  })
  axis(3, at=c(-10e10, 10e10), tck=0, labels=NA)
  axis(3, tck=-0.015, labels=NA)
  axis(3, tick=F, line=-0.65)
  mtext(3, line=1.2, text=bquote("Genomic Inflation Statistic," ~ lambda["GC"]))
}

# Scatterplot of lambdas vs sample size
lambda.scatter <- function(del, dup, hpo.n, blue.bg=TRUE, pt.cex=0.85,
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
  plot(NA, xlim=xlims, ylim=ylims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  x.ax.at <- log10(logscale.major)
  y.ax.at <- seq(0, 2, 0.1)
  abline(h=y.ax.at, v=x.ax.at, col=grid.col)
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

  # Add P-values
  x.length <- par("usr")[2] - par("usr")[1]
  y.length <- ylims[2] - ylims[1]
  text(x=par("usr")[1] - 0.03 * x.length,
       y=par("usr")[4] - 0.15 * y.length,
       xpd=T, pos=4, col=cnv.colors[1], labels=bquote(.(del.fit.p)))
  text(x=par("usr")[1] - 0.03 * x.length,
       y=par("usr")[4] - 0.05 * y.length,
       xpd=T, pos=4, col=cnv.colors[2], labels=bquote(.(dup.fit.p)))

  # Add points
  sapply(1:nrow(plot.dat), function(i){
    points(x=log10(plot.dat$n)[i], y=plot.dat$dup[i], pch=21, cex=pt.cex,
           bg=cnv.colors[2], col=cnv.blacks[2])
    points(x=log10(plot.dat$n)[i], y=plot.dat$del[i], pch=21, cex=pt.cex,
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
  mtext(2, line=1.75, text=bquote("Genomic Inflation," ~ lambda["GC"]))
}

# Primary vs. secondary P-value scatterplot
primary.vs.secondary.scatter <- function(pvals.1, pvals.2, keep.idx=NULL,
                                         cutoff.1=6, cutoff.2=-log10(0.05),
                                         sig.color="black", nonsig.color="gray75",
                                         min.point.plot=-log10(0.5), pt.cex=0.2,
                                         blue.bg=TRUE,
                                         parmar=c(2.25, 2.25, 0.25, 0.25)){
  # Fill NA P-values with P=1
  pvals.1 <- as.data.frame(apply(pvals.1, 2, function(ps){ps[which(is.na(ps))] <- 1; return(ps)}))
  pvals.2 <- as.data.frame(apply(pvals.2, 2, function(ps){ps[which(is.na(ps))] <- 1; return(ps)}))

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
  cor.res <- cor.test(pvals$primary, pvals$secondary, use="complete.obs")
  pvals <- pvals[which(apply(pvals, 1, function(ps){max(ps, na.rm=T) >= min.point.plot})), ]
  ax.lims <- range(pvals, na.rm=T)
  pt.colors <- apply(pvals, 1, function(vals){
    vals <- as.numeric(vals)
    if(vals[1] >= cutoff.1 & vals[2] >= cutoff.2){
      sig.color
    }else{
      nonsig.color
    }
  })
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
  plot(NA, xlim=ax.lims, ylim=ax.lims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  abline(h=axTicks(2), v=axTicks(1), col=grid.col)
  rect(xleft=0, xright=min.point.plot,
       ybottom=0, ytop=min.point.plot,
       border=NA, bty="n", col=nonsig.color)
  points(x=seq(0, min.point.plot, length.out=25),
         y=rep(0, 25), pch=19, cex=pt.cex, col=pt.colors)
  points(y=seq(0, min.point.plot, length.out=25),
         x=rep(0, 25), pch=19, cex=pt.cex, col=pt.colors)

  # Add points
  points(pvals, pch=19, cex=pt.cex, col=pt.colors)

  # Add annotation lines
  abline(v=cutoff.1, h=cutoff.2, lty=5, col=graphabs.green)

  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, tck=-0.03, labels=NA, col=blueblack)
  axis(1, tick=F, line=-0.75)
  mtext(1, line=1.4, text=bquote("-log"[10] * "(" * italic("P")["Primary"] * ")"))
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, tck=-0.03, labels=NA, col=blueblack)
  axis(2, tick=F, line=-0.65, las=2)
  mtext(2, line=1.25, text=bquote("-log"[10] * "(" * italic("P")["Secondary"] * ")"))

  # Annotate with P-value and correlation coefficient
  y.length <- diff(par("usr")[3:4])
  text(x=mean(par("usr")[1:2]),
       y=par("usr")[4] - 0.025 * y.length,
       xpd=T, pos=1, cex=5/6,
       labels=bquote(R^2 == .(as.numeric(formatC(cor.res$estimate ^ 2, digits=3)))))
  text(x=mean(par("usr")[1:2]),
       y=par("usr")[4] - 0.125 * y.length,
       xpd=T, pos=1, cex=5/6,
       labels=format.pval(cor.res$p.value))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--del-cutoff"), default=10e-6, help="Significant P-value cutoff for deletions."),
  make_option(c("--dup-cutoff"), default=10e-6, help="Significant P-value cutoff for duplications."),
  make_option(c("--del-nomsig-bed"), default="./del_nomsig.bed",
              help="Path to BED file with coordinates of nominally significant deletion windows."),
  make_option(c("--dup-nomsig-bed"), default="./dup_nomsig.bed",
              help="Path to BED file with coordinates of nominally significant duplication windows."),
  make_option(c("--primary-vs-secondary-png-dim"), default=2.25,
              help="Height and width of primary vs. secondary P-value scatterplot, in inches.")
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
del.cutoff <- -log10(as.numeric(opts$`del-cutoff`))
dup.cutoff <- -log10(as.numeric(opts$`dup-cutoff`))
del.nomsig.bed <- opts$`del-nomsig-bed`
dup.nomsig.bed <- opts$`dup-nomsig-bed`
prim.vs.sec.png.dim <- opts$`primary-vs-secondary-png-dim`

# # DEV PARAMETERS
# primary.del.pvals.in <- "~/scratch/rCNV2_analysis_d2.DEL.meta_neg_log10_p.all_hpos.bed.gz"
# primary.dup.pvals.in <- "~/scratch/rCNV2_analysis_d2.DUP.meta_neg_log10_p.all_hpos.bed.gz"
# secondary.del.pvals.in <- "~/scratch/rCNV2_analysis_d2.DEL.meta_neg_log10_p_secondary.all_hpos.bed.gz"
# secondary.dup.pvals.in <- "~/scratch/rCNV2_analysis_d2.DUP.meta_neg_log10_p_secondary.all_hpos.bed.gz"
# hpos.in <- "~/scratch/rCNV2_analysis_d2.reordered_hpos.txt"
# hpo.samplesize.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/test_pval_qc"
# del.cutoff <- 6
# dup.cutoff <- 6
# del.nomsig.bed <- "~/scratch/del_nomsig.bed"
# dup.nomsig.bed <- "~/scratch/dup_nomsig.bed"
# prim.vs.sec.png.dim <- 2.25

# Load pvalue matrices & compute lambdas
del.1 <- load.pval.matrix(primary.del.pvals.in)
dup.1 <- load.pval.matrix(primary.dup.pvals.in)
del.2 <- load.pval.matrix(secondary.del.pvals.in)
dup.2 <- load.pval.matrix(secondary.dup.pvals.in)

# Write BED files of windows with at least one nominally significant association
del.nomsig.coords <- del.1$coords[which(apply(del.1$pvals, 1, function(vals){any(!is.na(vals) & vals<0.05)})), ]
write.table(del.nomsig.coords, del.nomsig.bed, col.names=F, row.names=F, sep="\t", quote=F)
system(paste("bgzip -f", del.nomsig.bed), wait=T)
dup.nomsig.coords <- dup.1$coords[which(apply(dup.1$pvals, 1, function(vals){any(!is.na(vals) & vals<0.05)})), ]
write.table(dup.nomsig.coords, dup.nomsig.bed, col.names=F, row.names=F, sep="\t", quote=F)
system(paste("bgzip -f", dup.nomsig.bed), wait=T)

# Read reordered HPOs
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
hpo.n <- load.hpo.samplesize(hpo.samplesize.in, hpos)

# Print median lambda
cat(paste("Median lambda, deletions only:", round(median(del.1$lambdas), 5), "\n"))
cat(paste("Median lambda, duplications only:", round(median(dup.1$lambdas), 5), "\n"))
cat(paste("Median lambda, deletions & duplications only:", round(median(c(del.1$lambdas, dup.1$lambdas)), 5), "\n"))

# Plot grid of QQs
n.plots.wide <- 9
n.plots.tall <- ceiling(length(hpos) / n.plots.wide)
dim.scalar <- 7.5 / 10
png(paste(out.prefix, "qq_grid_byHPO.png", sep="."),
    height=dim.scalar*300*n.plots.tall, width=dim.scalar*300*n.plots.wide,
    res=300)
par(mfrow=c(n.plots.tall, n.plots.wide))
for(hpo in hpos){
  mini.qq(del.1, dup.1, hpo, title.cex=0.42, axis.cex=0.7, blue.bg=FALSE)
}
dev.off()

# Horizontal dotplot of primary lambdas
pdf(paste(out.prefix, "primary_lambdas_byHPO.pdf", sep="."),
    height=6, width=8)
plot.lambdas(del.1, dup.1, hpos, parmar=c(0.25, 11, 2.25, 11))
dev.off()

# Scatterplot of lambdas vs sample size
pdf(paste(out.prefix, "primary_lambdas_vs_sampleSize.pdf", sep="."),
    height=2.4, width=2.6)
lambda.scatter(del.1, dup.1, hpo.n, blue.bg=FALSE)
dev.off()

# Scatterplots of primary vs. secondary P-values
png(paste(out.prefix, "primary_vs_secondary_pvalue.DEL.png", sep="."),
    height=prim.vs.sec.png.dim*300, width=prim.vs.sec.png.dim*300, res=300)
primary.vs.secondary.scatter(del.1$pvals, del.2$pvals, cutoff.1=del.cutoff,
                             sig.color=cnv.colors[1], nonsig.color=control.cnv.colors[1],
                             pt.cex=0.175, blue.bg=FALSE, parmar=c(2.5, 2.5, 0.15, 0.15))
dev.off()
png(paste(out.prefix, "primary_vs_secondary_pvalue.DUP.png", sep="."),
    height=prim.vs.sec.png.dim*300, width=prim.vs.sec.png.dim*300, res=300)
primary.vs.secondary.scatter(dup.1$pvals, dup.2$pvals, cutoff.1=dup.cutoff,
                             sig.color=cnv.colors[2], nonsig.color=control.cnv.colors[2],
                             pt.cex=0.175, blue.bg=FALSE, parmar=c(2.5, 2.5, 0.15, 0.15))
dev.off()
