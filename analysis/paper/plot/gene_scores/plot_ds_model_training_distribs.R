#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distributions of training data for dosage sensitivity scoring for rCNV2 final paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load meta-analysis stats and extract effect size estimates
load.lnors <- function(meta.in, xlist){
  meta <- read.table(meta.in, sep="\t", header=T, comment.char="")
  meta[which(!(meta$gene %in% xlist)), c("gene", "meta_lnOR")]
}

# Load BFDPs and subset to genes passing xlist
load.bfdps <- function(bfdps.in, xlist){
  stats <- read.table(bfdps.in, sep="\t", header=T, comment.char="")
  stats[which(!(stats$gene %in% xlist)), c("gene", "bfdp")]
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Violin plots of effect sizes per gene per training subgroup
plot.lnors <- function(lnors, true.genes, false.genes, CNV,
                       blue.bg=TRUE, pt.cex=0.2, 
                       parmar=c(1.1, 2.7, 0.2, 0.2)){
  # Collect plot data
  plot.dat <- list(lnors$meta_lnOR[which(lnors$gene %in% false.genes)],
                   lnors$meta_lnOR[which(!(lnors$gene %in% c(false.genes, true.genes)))],
                   lnors$meta_lnOR[which(lnors$gene %in% true.genes)])
  plot.dat <- lapply(plot.dat, function(lnors){log2(exp(lnors))})
  ylim <- quantile(lnors$meta_lnOR, probs=c(0.005, 0.995), na.rm=T)
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

  # Set CNV type-specific parameters
  if(CNV == "DEL"){
    colors <- c(ns.color, control.cnv.colors[1], cnv.colors[1])
    labels <- c("HS", "Other", "HI")
  }else{
    colors <- c(ns.color, control.cnv.colors[2], cnv.colors[2])
    labels <- c("TI", "Other", "TS")
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 3), ylim=ylim, 
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  y.ax.at <- -3:10
  abline(h=y.ax.at, col=grid.col)
  
  # Add swarms
  sapply(1:3, function(i){
    beeswarm(plot.dat[[i]], add=T, at=i-0.5, col=colors[i], pch=19, cex=pt.cex,
             corral="wrap", corralWidth=0.9, xpd=T)
    boxplot(plot.dat[[i]], add=T, at=i-0.5, border=cnv.blacks[CNV], col=NA, lty=1,
            outline=F, width=0.175, staplewex=0, yaxt="n", xaxt="n", ylab="", xlab="")
    # segments(x0=i-0.75, x1=i-0.25,
    #          y0=mean(plot.dat[[i]], na.rm=T), y1=mean(plot.dat[[i]], na.rm=T),
    #          col=cnv.blacks[CNV], lend="round")
    axis(1, tick=F, line=-0.9, at=i-0.5, labels=labels[i], col.axis=colors[i])
  })
  
  # Add y-axis
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(2, at=y.ax.at, tick=F, las=2, line=-0.7, col=blueblack)
  mtext(2, line=1.25, text=bquote(log[2]("Odds Ratio")))
}

# Scatterplot of BFDP vs meta-analysis effect size
plot.bfdp.vs.lnor <- function(bfdps, lnors, true.genes, false.genes, CNV,
                              add.legend=TRUE, pt.cex=0.2, blue.bg=TRUE,
                              parmar=c(2.4, 2.7, 0.3, 0.2)){
  # Get plot data
  g.dat <- merge(lnors, bfdps, all=F, by="gene", sort=F)
  plot.dat <- list(g.dat[which(g.dat$gene %in% false.genes), -1],
                   g.dat[which(!(g.dat$gene %in% c(false.genes, true.genes))), -1],
                   g.dat[which(g.dat$gene %in% true.genes), -1])
  plot.dat <- lapply(plot.dat, function(df){
    df$meta_lnOR <- log2(exp(df$meta_lnOR))
    return(df)
  })
  # xlim <- quantile(lnors$meta_lnOR, probs=c(0.001, 0.999), na.rm=T)
  xlim <- range(log2(exp(lnors$meta_lnOR)), na.rm=T)
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
  
  # Set CNV type-specific parameters
  if(CNV == "DEL"){
    colors <- c(ns.color, control.cnv.colors[1], cnv.colors[1])
    labels <- c("HS", "Other", "HI")
  }else{
    colors <- c(ns.color, control.cnv.colors[2], cnv.colors[2])
    labels <- c("TI", "Other", "TS")
  }

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlim, ylim=c(0, 1), 
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  x.ax.at <- axTicks(1)
  y.ax.at <- axTicks(2)
  abline(v=x.ax.at, h=y.ax.at, col=grid.col)
  
  # Add points
  sapply(c(2, 1, 3), function(i){
    points(plot.dat[[i]], pch=19, cex=pt.cex, col=colors[i])
  })
  sapply(c(2, 1, 3), function(i){
    avgs <- apply(plot.dat[[i]], 2, mean, na.rm=T)
    points(avgs[1], avgs[2], pch=23, col=cnv.blacks[CNV], bg=colors[i])
  })
  
  # Add X-axis
  axis(1, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(1, at=x.ax.at, col=blueblack, tck=-0.025, labels=NA)
  sapply(x.ax.at, function(x){
    axis(1, at=x, tick=F, line=-0.75)
  })
  mtext(1, line=1.5, text=bquote(log[2]("Odds Ratio")))
  
  # Add Y-axis
  axis(2, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(2, at=y.ax.at, col=blueblack, tck=-0.025, labels=NA)
  axis(2, at=y.ax.at, tick=F, line=-0.65, las=2)
  mtext(2, line=1.75, text="Bayesian FDP")
  
  # Add legend (if optioned)
  if(add.legend==TRUE){
    points(x=par("usr")[1]+(0.04*diff(par("usr")[1:2])),
           y=0.05, pch=23, col=cnv.blacks[CNV], bg="white")
    text(x=par("usr")[1]+(0.02*diff(par("usr")[1:2])), y=0.04, 
         pos=4, labels="Group mean", cex=5/6)
  }
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(beeswarm, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog meta_stats.tsv bfdps.tsv exclude.list true.genes false.genes CNV out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 7){
  stop(paste("Seven positional arguments required: meta_stats.tsv, bfdps.tsv, true.genes, false.genes, exclude.genes, CNV, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
meta.in <- args$args[1]
bfdps.in <- args$args[2]
true.genes.in <- args$args[3]
false.genes.in <- args$args[4]
xlist.in <- args$args[5]
CNV <- args$args[6]
out.prefix <- args$args[7]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# meta.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
# bfdps.in <- "~/scratch/rCNV.DEL.gene_abfs.tsv"
# true.genes.in <- "~/scratch/gold_standard.haploinsufficient.genes.list"
# false.genes.in <- "~/scratch/gold_standard.haplosufficient.genes.list"
# xlist.in <- "~/scratch/rCNV.DEL.training_blacklist.genes.list"
# CNV <- "DEL"
# out.prefix <- "~/scratch/test_gene_score_train_distribs"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load gene lists
true.genes <- unique(as.character(read.table(true.genes.in, header=F, sep="\t")[, 1]))
false.genes <- unique(as.character(read.table(false.genes.in, header=F, sep="\t")[, 1]))
xlist <- unique(as.character(read.table(xlist.in, header=F, sep="\t")[, 1]))

# Extract meta-analysis odds ratios
lnors <- load.lnors(meta.in, xlist)

# Plot odds ratios
pdf(paste(out.prefix, "gene_scoring_training_distribs.effect_sizes.pdf", sep="."),
    height=2.25, width=1.85)
plot.lnors(lnors, true.genes, false.genes, CNV, pt.cex=0.15, blue.bg=FALSE)
dev.off()

# Load BFDPs
bfdps <- load.bfdps(bfdps.in, xlist)

# Plot odds ratios
pdf(paste(out.prefix, "gene_scoring_training_distribs.bfdps.pdf", sep="."),
    height=2.25, width=2.4)
plot.bfdp.vs.lnor(bfdps, lnors, true.genes, false.genes, CNV, pt.cex=0.15, blue.bg=FALSE)
dev.off()
