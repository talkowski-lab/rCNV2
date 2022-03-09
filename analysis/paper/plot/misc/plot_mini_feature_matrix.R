#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot mini version of gene feature matrix for graphical abstract of rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot mini heatmap of features
plot.mini.heat <- function(feats, n.genes=500, n.feats=20, seed=2020){
  # Drop gene column & ensure numerics
  feats <- apply(feats[, -1], 2, as.numeric)
  
  # Downsample genes
  if(nrow(feats) > n.genes){
    set.seed(seed)
    feats <- feats[sample(1:nrow(feats), n.genes, replace=F), ]
  }
  if(ncol(feats) > n.feats){
    set.seed(seed)
    feats <- feats[, sample(1:ncol(feats), n.feats, replace=F)]
  }
  
  # Prep plot area
  par(mar=c(1.2, 1.2, 0.2, 0.2), bty="n")
  plot(NA, xlim=c(0, ncol(feats)), ylim=c(0, nrow(feats)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  
  # Plot cells
  sapply(1:ncol(feats), function(i){
    vals <- (feats[, i] - min(feats[, i], na.rm=T))
    # pal <- viridis
    pal <- colorRampPalette(c(graphabs.green.darker, graphabs.green, viridis(5)[4:5]))
    # pal <- colorRampPalette(c(cnv.colors[1], cnv.whites[1:2], cnv.colors[2]))
    colors <- pal(101)[round(100 * vals / max(vals, na.rm=T)) + 1]
    rect(xleft=i-1, xright=i, ybottom=(1:nrow(feats))-1, ytop=1:nrow(feats),
         border=colors, lwd=0.25, col=colors)
  })
  
  # Add cleanup
  box(col=blueblack, lty="solid", bty="o")
  mtext(1, line=0, text="Features", cex=5/6, col=blueblack, font=3)
  mtext(2, line=0.1, text="Genes", cex=5/6, col=blueblack, font=3)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog features.bed.gz out.pdf", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: features.bed.gz and output.pdf\n", sep=" "))
}

# Writes args & opts to vars
features.in <- args$args[1]
out.pdf <- args$args[2]

# # DEV PARAMETERS
# features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz"
# out.pdf <- "~/scratch/test_mini_feat_matrix.pdf"

# Load features
feats <- load.features(features.in, norm=T)

# Plot feature matrix
pdf(out.pdf, height=1, width=0.85)
plot.mini.heat(feats, n.genes=60, n.feats=25)
dev.off()
