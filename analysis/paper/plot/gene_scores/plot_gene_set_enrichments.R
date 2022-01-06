#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot gene set enrichments vs. dosage sensitivity scores for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load multiple gene sets from .tsv input
load.genesets <- function(genesets.in, elig.genes=NULL){
  glist <- read.table(genesets.in, header=F, sep="\t")
  gsets <- lapply(1:nrow(glist), function(i){
    path <- glist[i, 2]
    genes <- sort(unique(as.character(read.table(path, header=F)[, 1])))
    if(!is.null(elig.genes)){
      genes <- genes[which(genes %in% elig.genes)]
    }
    return(list("genes"=genes, "train"=as.logical(glist[i, 3])))
  })
  names(gsets) <- glist[, 1]
  return(gsets)
}

# Split genes into groups based on score bin
bin.genes <- function(scores, score, n.bins=10){
  breaks <- seq(0, 1, length.out=n.bins+1)
  lapply(2:length(breaks), function(i){
    scores$gene[which(scores[, score] >= breaks[i-1]
                      & scores[, score] <= breaks[i])]
  })
}

# Compute the enrichment of a single set of genes vs. score bins
calc.enrichment.single <- function(gset, score.groups){
  genes <- gset$genes
  n.genes <- length(genes)
  n.univ <- length(unlist(score.groups))
  baseline <- n.genes / n.univ
  enrich <- as.data.frame(t(sapply(score.groups, function(denom){
    n.denom <- length(denom)
    n.ovr <- length(which(genes %in% denom))
    frac <- n.ovr / n.genes
    fold <- (n.ovr / n.denom) / baseline
    binom.p <- binom.test(n.ovr, n.denom, p=baseline)$p.value
    return(c(n.ovr, n.denom, frac, fold, binom.p))
  })))
  colnames(enrich) <- c("hits", "denom", "frac", "fold", "binom.p")
  enrich$bonf <- (enrich$binom.p <= 0.05 / length(score.groups))
  return(list("enrich"=enrich, "train"=gset$train))
}

# Compute enrichments of all gene sets vs. score bins
calc.enrichments <- function(gsets, scores, score, n.bins=10){
  score.groups <- bin.genes(scores, score, n.bins)
  enrich.df <- lapply(gsets, calc.enrichment.single, score.groups)
  names(enrich.df) <- names(gsets)
  return(enrich.df)
}

# K-S test of scores for all gene sets
geneset.ks.tests <- function(gsets, scores, score){
  sapply(1:length(gsets), function(i){
    cat(paste("K-S test for ", names(gsets)[i], ":\n", sep=""))
    scores.hits <- as.numeric(scores[which(scores$gene %in% gsets[[i]]$genes), score])
    ks.p <- ks.test(scores[, score], scores.hits)$p.value
    cat(paste("P =", format(ks.p, scientific=T, digits=3, nsmall=2), "\n"))
  })
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot enrichments as heatmap
plot.enrichments <- function(enrich.df, score, max.pct=0.5, x.tck=-0.025,
                             parmar=c(0, 11.5, 2.15, 0.2)){
  # Set score-dependent plotting parameters
  if(score == "pHaplo"){
    score.pal <- colorRampPalette(c("white", cnv.colors[1]))(101)
    bonf.col <- cnv.blacks[1]
  }else if(score == "pTriplo"){
    score.pal <- colorRampPalette(c("white", cnv.colors[2]))(101)
    bonf.col <- cnv.blacks[2]
  }

  # Collect plot parameters
  n.rows <- length(enrich.df)
  n.cols <- nrow(enrich.df[[1]]$enrich)
  box.lefts <- 0:(n.cols-1)
  box.rights <- 1:n.cols
  x.ax.at <- seq(0, n.cols, length.out=6)
  x.ax.labels <- seq(0, 1, 0.2)

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(0, n.cols), ylim=c(n.rows, 0),
       xaxt="n", yaxt="n", xlab="", ylab="")

  # Add rectangles
  sapply(1:length(enrich.df), function(i){
    col.idxs <- (floor(100*enrich.df[[i]]$enrich$frac)+1)/max.pct
    col.idxs[which(col.idxs>=101)] <- 101
    rect(xleft=box.lefts, xright=box.rights,
         ybottom=i-1, ytop=i,
         border=bluewhite, col=score.pal[col.idxs])
  })

  # Add significance markers
  sapply(1:length(enrich.df), function(i){
    up.idx <- which(enrich.df[[i]]$enrich$bonf & enrich.df[[i]]$enrich$fold > 1)
    down.idx <- which(enrich.df[[i]]$enrich$bonf & enrich.df[[i]]$enrich$fold < 1)
    arrow.colors <- rep(bonf.col, n.cols)
    # maxed.idx <- which(enrich.df[[i]]$enrich$frac >= 0.8 * max.pct)
    # if(length(maxed.idx) > 0){
    #   arrow.colors[maxed.idx] <- "white"
    # }
    if(length(up.idx) > 0){
      Arrows(x0=up.idx-0.5, x1=up.idx-0.5, y0=i-0.25, y1=i-0.75, col=arrow.colors[up.idx],
             arr.type="triangle", arr.length=0.12, arr.width=0.1, arr.adj=1)
    }
    if(length(down.idx) > 0){
      Arrows(x0=down.idx-0.5, x1=down.idx-0.5, y0=i-0.75, y1=i-0.25, col=arrow.colors[down.idx],
             arr.type="triangle", arr.length=0.12, arr.width=0.1, arr.adj=1)
    }
  })

  # Add labels
  sapply(1:n.rows, function(i){
    y.label <- names(enrich.df)[i]
    if(enrich.df[[i]]$train==TRUE){
      y.label <- paste(y.label, "*", sep="")
    }
    n.genes <- sum(enrich.df[[i]]$enrich$hits)
    y.label <- paste(y.label, " (n=", prettyNum(n.genes, big.mark=","), ")", sep="")
    axis(2, at=i-0.5, tick=F, line=-1, las=2, cex=5.5/6, labels=y.label)
  })

  # Add x-axis label
  axis(3, at=x.ax.at, labels=NA, tck=x.tck, col=blueblack)
  sapply(1:length(x.ax.at), function(i){
    axis(3, at=x.ax.at[i], tick=F, line=-0.8, labels=x.ax.labels[i])
  })
  mtext(3, line=1.2, text=score)
}

# Plot color gradient legend
plot.gradient.legend <- function(score, max.pct=0.5){
  # Set score-dependent plotting parameters
  if(score == "pHaplo"){
    score.pal <- colorRampPalette(c("white", cnv.colors[1]))(101)
  }else if(score == "pTriplo"){
    score.pal <- colorRampPalette(c("white", cnv.colors[2]))(101)
  }else if(score == "both"){
    score.pal <- list(colorRampPalette(c("white", cnv.colors[1]))(101),
                       colorRampPalette(c("white", cnv.colors[2]))(101))
  }

  # Plot legend
  par(bty="n", mar=c(0.05, 1.2, 1, 2.1))
  plot(NA, xlim=c(0, 100), ylim=c(0, 1),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  if(score == "both"){
    rect(xleft=0:99, xright=1:100, ybottom=0.5, ytop=1, col=score.pal[[1]], border=score.pal[[1]], lwd=0.5)
    rect(xleft=0:99, xright=1:100, ybottom=0, ytop=0.5, col=score.pal[[2]], border=score.pal[[2]], lwd=0.5)
  }else{
    rect(xleft=0:99, xright=1:100, ybottom=0, ytop=1, col=score.pal, border=score.pal, lwd=0.5)
  }
  box(bty="o", col=blueblack, xpd=T)

  # Add labels
  # axis(3, at=c(0, 100), tck=-0.25, col=blueblack, labels=NA)
  axis(2, at=0.5, line=-0.9, tick=F, labels="0%", cex.axis=5/6, las=2)
  axis(4, at=0.5, line=-1.1, tick=F, labels=bquote("" >= .(round(100*max.pct, 0)) * "%"), cex.axis=5/6, las=2)
  mtext(3, text="Fraction of gene set", cex=5.5/6, line=0.1)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(shape, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--height"), help="Height of output pdf, in inches.", default=2)
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv genesets.tsv score",
                                            "out_prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, genesets.tsv,",
             "score, and out_prefix\n"))
}

# Writes args & opts to vars
scores.in <- args$args[1]
genesets.in <- args$args[2]
score <- args$args[3]
out.prefix <- args$args[4]
pdf.height <- opts$height

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# genesets.in <- "~/scratch/phi_vs_gene_sets.input.tsv"
# score <- "pHaplo"
# out.prefix <- "~/scratch/geneset_enrichments_test"
# pdf.height <- 1.8

# Load scores
scores <- load.scores(scores.in)

# Load gene sets
gsets <- load.genesets(genesets.in, elig.genes=scores$gene)

# Perform K-S tests for all gene sets
geneset.ks.tests(gsets, scores, score)

# Stratify gene sets per score bin
enrich.df <- calc.enrichments(gsets, scores, score, n.bins=10)

# Plot enrichments
pdf(paste(out.prefix, score, "geneset_enrichments.pdf", sep="."),
    height=pdf.height, width=4)
plot.enrichments(enrich.df, score, x.tck=-0.05/pdf.height)
dev.off()

# Plot legend
pdf(paste(out.prefix, score, "geneset_enrichments.legend.pdf", sep="."),
    height=0.35, width=1.8)
plot.gradient.legend(score)
dev.off()

# Plot split legend only for pHaplo
if(score=="pHaplo"){
  pdf(paste(out.prefix, score, "geneset_enrichments.double_legend.pdf", sep="."),
      height=0.35, width=1.8)
  plot.gradient.legend("both")
  dev.off()
}
