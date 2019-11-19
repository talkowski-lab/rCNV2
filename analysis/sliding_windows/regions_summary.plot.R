#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot summary schematics for final segments from sliding window analysis


# Set global parameters
options(scipen=100000, stringsAsFactors=F)
cnv.colors <- c("DEL"="#D43925",
                "DUP"="#2376B2")
logscale.size <- log10(as.vector(sapply(5:7, function(e){(1:9)*10^e})))
logscale.size.labels.at <- 5:7
logscale.size.labels <- c("100kb", "1Mb", "10Mb")


#################
### FUNCTIONS ###
#################
# Load a BED file of loci
load.loci <- function(path){
  dat <- read.table(path, header=T, sep="\t", comment.char="")
  dat$npheno <- unlist(lapply(strsplit(dat$hpos, split=";", fixed=T), length))
  dat$ngenes <- unlist(lapply(strsplit(dat$genes, split=";", fixed=T), length))
  colnames(dat)[1] <- "chr"
  return(dat)
}

# Plot swarms of region sizes
plot.size <- function(DEL, DUP, title.cex=1){
  require(beeswarm)
  par(mar=c(2, 4, 1.5, 0.5))
  beeswarm(list("DEL"=log10(DEL$size),
                "DUP"=log10(DUP$size)), 
           col=cnv.colors, pch=19, ylim=c(5, 7),
           xaxt="n", xlab="", yaxt="n", ylab="",
           corral="wrap", corralWidth=0.4)
  size.means <- sapply(list(DEL, DUP), function(x){mean(log10(x$size))})
  segments(x0=(1:2)-0.2, x1=(1:2)+0.2, y0=size.means, y1=size.means, lend="round", lwd=2)
  axis(2, at=logscale.size, labels=NA, tck=-0.02, col="gray50")
  axis(2, at=logscale.size.labels.at, labels=NA, tck=-0.04)
  axis(2, at=logscale.size.labels.at, labels=logscale.size.labels,
       las=2, line=-0.4, tick=F)
  mtext(2, text="Credible Region Size", line=2.25, cex=title.cex)
  axis(1, at=1, tick=F, line=-0.8, font=2, col.axis=cnv.colors[1], labels="DEL")
  axis(1, at=2, tick=F, line=-0.8, font=2, col.axis=cnv.colors[2], labels="DUP")
  mtext(3, line=0.2, font=2, text="Credible Region Size", cex=title.cex)
  box()
}

# Plot swarms of count of genes
plot.genes <- function(DEL, DUP, title.cex=1){
  require(beeswarm)
  par(mar=c(2, 4, 1.5, 0.5))
  beeswarm(list("DEL"=DEL$ngenes,
                "DUP"=DUP$ngenes), 
           col=cnv.colors, pch=19, 
           xaxt="n", xlab="", yaxt="n", ylab="",
           corral="wrap", corralWidth=0.4)
  gene.means <- c(mean(DEL$ngenes), mean(DUP$ngenes))
  segments(x0=(1:2)-0.2, x1=(1:2)+0.2, y0=gene.means, y1=gene.means, lend="round", lwd=2)
  axis(2, labels=NA, tck=-0.04)
  axis(2, las=2, line=-0.4, tick=F)
  mtext(2, text="Genes per Region", line=2.25, cex=title.cex)
  axis(1, at=1, tick=F, line=-0.8, font=2, col.axis=cnv.colors[1], labels="DEL")
  axis(1, at=2, tick=F, line=-0.8, font=2, col.axis=cnv.colors[2], labels="DUP")
  mtext(3, line=0.2, font=2, text="Genes per Region", cex=title.cex)
  box()
}

# Plot swarms of count of HPOs
plot.hpos <- function(DEL, DUP, title.cex=1){
  require(beeswarm)
  par(mar=c(2, 4, 1.5, 0.5))
  beeswarm(list("DEL"=DEL$npheno,
                "DUP"=DUP$npheno), 
           col=cnv.colors, pch=19, 
           xaxt="n", xlab="", yaxt="n", ylab="",
           corral="wrap", corralWidth=0.4)
  gene.means <- c(mean(DEL$npheno), mean(DUP$npheno))
  segments(x0=(1:2)-0.2, x1=(1:2)+0.2, y0=gene.means, y1=gene.means, lend="round", lwd=2)
  axis(2, labels=NA, tck=-0.04)
  axis(2, las=2, line=-0.4, tick=F)
  mtext(2, text="Associated HPOs", line=2.25, cex=title.cex)
  axis(1, at=1, tick=F, line=-0.8, font=2, col.axis=cnv.colors[1], labels="DEL")
  axis(1, at=2, tick=F, line=-0.8, font=2, col.axis=cnv.colors[2], labels="DUP")
  mtext(3, line=0.2, font=2, text="HPOs per Region", cex=title.cex)
  box()
}

# Plot scatter of size vs genes
plot.sizeVsGenes <- function(DEL, DUP, title.cex=1){
  par(mar=c(3, 4, 1.5, 0.5))
  sizes <- c(DEL$size, DUP$size)
  genes <- c(DEL$ngenes, DUP$ngenes)
  plot(x=sizes, y=genes, pch=19, col=cnv.colors,
       xaxt="n", xlab="", yaxt="n", ylab="",
       panel.first=c(abline(lm(genes ~ sizes), lty=2, col="gray50")),
       ylim=c(0, max(genes)))
  axis(1, labels=NA, tck=-0.04)
  axis(1, at=axTicks(1), labels=axTicks(1)/1000000, tick=F, line=-0.5)
  mtext(1, line=1.5, text="Region Size (Mb)", cex=title.cex)
  axis(2, labels=NA, tck=-0.04)
  axis(2, las=2, line=-0.4, tick=F)
  mtext(2, text="Genes per Region", line=2.25, cex=title.cex)
  text(x=par("usr")[2], y=0.025*par("usr")[4], pos=2,
       labels=bquote(italic(R)^2 == .(round(cor(sizes, genes), 3))))
  mtext(3, line=0.2, font=2, text="Size vs. Genes", cex=title.cex)
}

# Plot scatter of size vs HPOs
plot.sizeVsHpos <- function(DEL, DUP, title.cex=1){
  par(mar=c(3, 4, 1.5, 0.5))
  sizes <- c(DEL$size, DUP$size)
  hpos <- c(DEL$npheno, DUP$npheno)
  plot(x=sizes, y=hpos, pch=19, col=cnv.colors,
       xaxt="n", xlab="", yaxt="n", ylab="",
       panel.first=c(abline(lm(hpos ~ sizes), lty=2, col="gray50")),
       ylim=c(0, max(hpos)))
  axis(1, labels=NA, tck=-0.04)
  axis(1, at=axTicks(1), labels=axTicks(1)/1000000, tick=F, line=-0.5)
  mtext(1, line=1.5, text="Region Size (Mb)", cex=title.cex)
  axis(2, labels=NA, tck=-0.04)
  axis(2, las=2, line=-0.4, tick=F)
  mtext(2, text="Associated HPOs", line=2.25, cex=title.cex)
  text(x=par("usr")[2], y=0.025*par("usr")[4], pos=2,
       labels=bquote(italic(R)^2 == .(round(cor(sizes, hpos), 3))))
  mtext(3, line=0.2, font=2, text="Size vs. HPOs", cex=title.cex)
}

# Plot scatter of genes vs HPOs
plot.genesVsHpos <- function(DEL, DUP, title.cex=1){
  par(mar=c(3, 4, 1.5, 0.5))
  genes <- c(DEL$ngenes, DUP$ngenes)
  hpos <- c(DEL$npheno, DUP$npheno)
  plot(x=genes, y=hpos, pch=19, col=cnv.colors,
       xaxt="n", xlab="", yaxt="n", ylab="", 
       panel.first=c(abline(lm(hpos ~ genes), lty=2, col="gray50")),
       ylim=c(0, max(hpos)))
  axis(1, labels=NA, tck=-0.04)
  axis(1, tick=F, line=-0.5)
  mtext(1, line=1.5, text="Genes per Region", cex=title.cex)
  axis(2, labels=NA, tck=-0.04)
  axis(2, las=2, line=-0.4, tick=F)
  mtext(2, text="Associated HPOs", line=2.25, cex=title.cex)
  text(x=par("usr")[2], y=0.025*par("usr")[4], pos=2,
       labels=bquote(italic(R)^2 == .(round(cor(sizes, hpos), 3))))
  mtext(3, line=0.2, font=2, text="Genes vs. HPOs", cex=title.cex)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("-o", "--out-prefix"), type="character", default="./final_regions.",
              help="prefix for writing out all results. [default %default]", metavar="path")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog DEL_regions DUP_regions",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Must supply DEL and DUP regions as positional arguments.\n")
}

# Writes args & opts to vars
DEL.in <- args$args[1]
DUP.in <- args$args[2]
out.prefix <- opts$`out-prefix`

# # DEV PARAMETERS:
# DEL.in <- "~/scratch/rCNV.DEL.final_regions.loci.bed.gz"
# DUP.in <- "~/scratch/rCNV.DUP.final_regions.loci.bed.gz"
# out.prefix <- "~/scratch/sig_regions.test."

# Load data
DEL <- load.loci(DEL.in)
DUP <- load.loci(DUP.in)

# Plot size
pdf(paste(out.prefix, "region_sizes.pdf", sep=""),
    height=3, width=3)
plot.size(DEL, DUP)
dev.off()

# Plot count of genes
pdf(paste(out.prefix, "gene_count.pdf", sep=""),
    height=3, width=3)
plot.genes(DEL, DUP)
dev.off()

# Plot count of HPOs
pdf(paste(out.prefix, "hpo_count.pdf", sep=""),
    height=3, width=3)
plot.hpos(DEL, DUP)
dev.off()

# Plot size vs genes
pdf(paste(out.prefix, "size_vs_genes.pdf", sep=""),
    height=3, width=3)
plot.sizeVsGenes(DEL, DUP)
dev.off()

# Plot size vs HPOs
pdf(paste(out.prefix, "size_vs_hpos.pdf", sep=""),
    height=3, width=3)
plot.sizeVsHpos(DEL, DUP)
dev.off()

# Plot genes vs HPOs
pdf(paste(out.prefix, "genes_vs_hpos.pdf", sep=""),
    height=3, width=3)
plot.genesVsHpos(DEL, DUP)
dev.off()

# Plot combined six-panel figure
pdf(paste(out.prefix, "multipanel_summary.pdf", sep=""),
    height=4, width=6)
layout(matrix(1:6, nrow=2, byrow=T))
plot.size(DEL, DUP, title.cex=0.75)
plot.genes(DEL, DUP, title.cex=0.75)
plot.hpos(DEL, DUP, title.cex=0.75)
plot.sizeVsGenes(DEL, DUP, title.cex=0.75)
plot.sizeVsHpos(DEL, DUP, title.cex=0.75)
plot.genesVsHpos(DEL, DUP, title.cex=0.75)
dev.off()
