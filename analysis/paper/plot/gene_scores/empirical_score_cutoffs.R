#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Determine optimal pHaplo & pTriplo score cutoffs based on empirical odds ratios


# Set parameters & define constants
options(stringsAsFactors=F, scipen=1000, family="sans")


######################
### DATA FUNCTIONS ###
######################
# Load genes impacted by CNVs from a single cohort
load.cohort.cnvs <- function(cnvs.in){
  x <- read.table(cnvs.in, header=T, sep="\t", comment.char="")
  x <- x[which(x$ngenes>0), c("phenos", "genes")]
  x$iscase <- strsplit(x$phenos, split=";") != c("HEALTHY_CONTROL")
  x$phenos <- NULL
  x$genes <- strsplit(x$genes, split=";")
  return(x)
}

# Load all data required for on-the-fly meta-analysis
load.meta.dat <- function(meta.inputs.in){
  inputs <- read.table(meta.inputs.in, sep="\t", header=T)
  cohorts <- as.character(inputs$cohort)
  n_case <- as.numeric(inputs$n_case)
  n_control <- as.numeric(inputs$n_control)
  names(n_control) <- names(n_case) <- cohorts
  cnvs <- apply(inputs[, c("DEL_path", "DUP_path")], 1, function(paths){
    return(list("DEL" = load.cohort.cnvs(paths[1]),
                "DUP" = load.cohort.cnvs(paths[2])))
  })
  names(cnvs) <- cohorts
  return(list("cnvs" = cnvs, "n_case" = n_case,
              "n_control" = n_control, "cohorts" = cohorts))
}

# Compute on-the-fly meta-analysis of all CNVs that hit a list of genes
gene.meta.otf <- function(meta.dat, genes, cnv){
  stats.merged <- do.call("cbind", lapply(1:length(meta.dat$cnvs), function(i){
    cnvs <- meta.dat$cnvs[[i]][[cnv]]
    cohort <- meta.dat$cohorts[i]
    ncase.all <- meta.dat$n_case[i]
    nctrl.all <- meta.dat$n_control[i]
    hits <- which(sapply(cnvs$genes, function(glist){length(intersect(glist, genes)) > 0}))
    hits.iscase <- table(cnvs$iscase[hits])
    for(lab in c("TRUE", "FALSE")){
      if(!(lab %in% names(hits.iscase))){
        names <- c(names(hits.iscase), lab)
        hits.iscase <- c(hits.iscase, 0)
        names(hits.iscase) <- names
      }
    }
    ncase.alt <- hits.iscase["TRUE"]
    ncase.ref <- ncase.all - ncase.alt
    nctrl.alt <- hits.iscase["FALSE"]
    nctrl.ref <- nctrl.all - nctrl.alt
    outvals <- data.frame(ncase.alt, ncase.ref, nctrl.alt, nctrl.ref)
    colnames(outvals) <- paste(cohort, c("case_alt", "case_ref", "control_alt",
                                         "control_ref"), sep=".")
    return(outvals)
  }))
  stats.merged$exclude_cohorts <- ""
  meta.single(stats.merged, meta.dat$cohorts, row.idx=1)
}

# Compute effect size per score bin
lnor.per.score.bin <- function(meta.dat, scores, score, cnv, xlist=c(),
                               bins=100, start=1, end=0){
  score.breaks <- seq(start, end, length.out=bins+1)
  res <- as.data.frame(t(sapply(1:bins, function(i){
    upper.genes <- scores$gene[which(scores[, score] >= score.breaks[i+1])]
    upper.genes <- setdiff(upper.genes, xlist)
    n.upper.genes <- length(upper.genes)
    upper.lnor <- gene.meta.otf(meta.dat, upper.genes, cnv)[1:3]
    bin.genes <- intersect(upper.genes, scores$gene[which(scores[, score] < score.breaks[i])])
    n.bin.genes <- length(bin.genes)
    bin.lnor <- gene.meta.otf(meta.dat, bin.genes, cnv)[1:3]
    return(c(score.breaks[i+1], n.upper.genes, upper.lnor, n.bin.genes, bin.lnor))
  })))
  colnames(res) <- c("cutoff",
                     "cumul_genes", "cumul_lnOR", "cumul_lnOR_lower", "cumul_lnOR_upper",
                     "margin_genes", "margin_lnOR", "margin_lnOR_lower", "margin_lnOR_upper")
  return(res)
}

# Derive comparable score cutoff vs constrained gene deletions
get.cutoff <- function(score.lnor.full, score.lnor.fine, min.lnor){
  full.idxs <- which(score.lnor.full$cumul_lnOR >= min.lnor)
  if(length(full.idxs) > 0){
    return(score.lnor.full$cutoff[max(full.idxs)])
  }else{
    fine.idxs <- which(score.lnor.fine$cumul_lnOR >= min.lnor)
    return(score.lnor.fine$cutoff[max(fine.idxs)])
  }
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot empirical effect sizes vs score bin
plot.score.vs.or <- function(bins, baseline, score.cutoff, score,
                             avg.pt.cex=1, avg.genes.per.bin=NULL,
                             x.ax.at=NULL, ylims=NULL, xtitle=NULL, ytitle=NULL,
                             blue.bg=TRUE, ax.tck=-0.025,
                             parmar=c(2.4, 2.4, 0.25, 0.25)){
  # Get plot values
  x <- bins$cutoff
  pt.y <- bins$margin_lnOR
  line.y <- bins$cumul_lnOR

  # Set plot parameters
  if(score == "pHaplo"){
    col.idx <- 1
    cnv.label <- "Deletion"
  }else{
    col.idx <- 2
    cnv.label <- "Duplication"
  }
  pt.col <- rep(control.cnv.colors[col.idx], nrow(bins))
  line.col <- cnv.colors[col.idx]
  baseline.col <- graphabs.green
  highlight.idxs <- which(bins$cutoff >= score.cutoff)
  pt.col[highlight.idxs] <- cnv.colors[col.idx]
  highlight.line.col <- cnv.blacks[col.idx]
  if(is.null(avg.genes.per.bin)){
    avg.genes.per.bin <- mean(bins$margin_genes, na.rm=T)
  }
  pt.cex <- avg.pt.cex*sqrt(bins$margin_genes/avg.genes.per.bin)
  xlims <- range(x, na.rm=T)
  if(is.null(ylims)){
    ylims <- range(c(pt.y, line.y), na.rm=T)
  }
  if(is.null(xtitle)){
    xtitle <- paste("Genes Binned by", score)
  }
  if(is.null(ytitle)){
    ytitle <- paste(cnv.label, "Odds Ratio")
  }
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
       xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=score.cutoff,
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=plot.bg, border=plot.border, bty=plot.bty)
  rect(xleft=score.cutoff, xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=adjustcolor(highlight.color, alpha=0.3), border=NA, bty="n")
  y.ax.at <- log(2^(-10:10))
  abline(h=y.ax.at, col=grid.col)
  abline(v=score.cutoff, lty=1, col=highlight.color)

  # Add points & lines
  abline(h=baseline, lty=5, col=baseline.col)
  points(x, pt.y, col=pt.col, cex=pt.cex, pch=19)
  points(x, line.y, type="l", col=line.col, lwd=2)
  points(x[highlight.idxs], line.y[highlight.idxs],
         type="l", col=highlight.line.col, lwd=2)
  points(x=score.cutoff, y=baseline, pch=23, cex=1.5*avg.pt.cex,
         bg=highlight.color, col=baseline.col, lwd=2)

  # Add axes
  if(is.null(x.ax.at)){
    x.ax.at <- axTicks(1)
  }
  axis(1, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(1, at=x.ax.at, tck=ax.tck, labels=NA, col=blueblack)
  axis(1, at=x.ax.at, tick=F, line=-0.7)
  mtext(1, line=1.3, text=xtitle, xpd=T)
  axis(2, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(2, at=y.ax.at, tck=ax.tck, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, tick=F, labels=exp(y.ax.at), las=2, line=-0.7)
  mtext(2, line=1.35, text=ytitle)
  axis(3, at=c(score.cutoff, par("usr")[2]), col=blueblack, tck=-ax.tck)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(metafor, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv meta_inputs.tsv constr.genes exclude.list out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: scores.tsv, meta_inputs.tsv, constr.genes.list, exclude.genes, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
meta.inputs.in <- args$args[2]
constr.genes.in <- args$args[3]
xlist.in <- args$args[4]
out.prefix <- args$args[5]

# # DEV PARAMETERS
# setwd("~/scratch/")
# scores.in <- "rCNV.gene_scores.tsv.gz"
# meta.inputs.in <- "empirical_score_cutoff.meta_inputs.tsv"
# constr.genes.in <- "gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
# xlist.in <- "rCNV.gene_scoring.excluded_training_genes.list"
# out.prefix <- "test_gene_score_empirical_cutoffs"

# Load scores
scores <- load.scores(scores.in)

# Load gene lists
constr.genes <- unique(as.character(read.table(constr.genes.in, header=F, sep="\t")[, 1]))
xlist <- unique(as.character(read.table(xlist.in, header=F, sep="\t")[, 1]))

# Load data necessary for on-the-fly meta-analyses
meta.dat <- load.meta.dat(meta.inputs.in)

# Compute lnOR for deletions of constrained genes
del.constr.lnor <- gene.meta.otf(meta.dat, setdiff(constr.genes, xlist), "DEL")[1]

# TODO: KEEP IMPLEMENTING THIS

# Compute lnOR estimates for genes by pHaplo & pTriplo bin
phi.lnor.full <- lnor.per.score.bin(meta.dat, scores, "pHaplo", "DEL", xlist,
                                    bins=100, start=1, end=0)
pts.lnor.full <- lnor.per.score.bin(meta.dat, scores, "pTriplo", "DUP", xlist,
                                    bins=100, start=1, end=0)
# pts.lnor.fine <- avg.lnor.per.score.bin(dup.lnors, scores, "pTriplo", bins=100, start=1, end=0.9)
# pts.lnor.fine.forplot <- avg.lnor.per.score.bin(dup.lnors, scores, "pTriplo", bins=50, start=1, end=0.95)

# Derive cutoffs
phi.cutoff <- get.cutoff(phi.lnor.full, phi.lnor.fine, del.constr.lnor[1])
pts.cutoff <- get.cutoff(pts.lnor.full, pts.lnor.fine, del.constr.lnor[1])

# Print derived cutoffs
cat(paste("\nAverage rare deletion of constrained gene confers odds ratio =",
          round(exp(del.constr.lnor[1]), 4), "\n"))
cat(paste("\nComparable pHaplo cutoff >=", phi.cutoff, "(includes",
          prettyNum(length(which(scores$pHaplo>=phi.cutoff)), big.mark=","),
          "genes)\n"))
cat(paste("\nComparable pTriplo cutoff >=", pts.cutoff, "(includes",
          prettyNum(length(which(scores$pTriplo>=pts.cutoff)), big.mark=","),
          "genes)\n"))

# Set standardized parameters for all plots
ylims <- quantile(c(phi.lnor.full$margin_lnOR, pts.lnor.full$margin_lnOR),
                  probs=c(0.005, 0.995), na.rm=T)
pt.cex <- 3/4
avg.genes.per.bin <- mean(c(phi.lnor.full$margin_genes, pts.lnor.full$margin_genes), na.rm=T)
pdf.height <- 2.25
pdf.width <- 2.75
pdf.parmar <- c(2.3, 2.3, 0.3, 0.3)

# Plot full pHaplo vs lnOR
pdf(paste(out.prefix, "phi_vs_effect_size.pdf", sep="."),
    height=pdf.height, width=pdf.width)
plot.score.vs.or(phi.lnor.full, del.constr.lnor[1], phi.cutoff, "pHaplo",
                 ylims=ylims, avg.pt.cex=pt.cex, avg.genes.per.bin=avg.genes.per.bin,
                 blue.bg=FALSE, parmar=pdf.parmar)
dev.off()

# Plot full pTriplo vs lnOR
pdf(paste(out.prefix, "pts_vs_effect_size.pdf", sep="."),
    height=pdf.height, width=pdf.width)
plot.score.vs.or(pts.lnor.full, del.constr.lnor[1], pts.cutoff, "pTriplo",
                 ylims=ylims, avg.pt.cex=pt.cex, avg.genes.per.bin=avg.genes.per.bin,
                 blue.bg=FALSE, parmar=pdf.parmar)
dev.off()

# # Plot fine pTriplo vs lnOR
# pdf(paste(out.prefix, "pts_vs_effect_size.fine.pdf", sep="."),
#     height=(5/6)*pdf.height, width=(1/2)*pdf.width)
# plot.score.vs.or(pts.lnor.fine.forplot, del.constr.lnor[1], pts.cutoff, "pTriplo",
#                  ylims=ylims, avg.pt.cex=pt.cex, xtitle="pTriplo", ytitle=NA,
#                  blue.bg=FALSE, parmar=c(pdf.parmar[1], 1.3, pdf.parmar[3:4]))
# axis(1, at=1, tick=F, line=-0.7)
# dev.off()

