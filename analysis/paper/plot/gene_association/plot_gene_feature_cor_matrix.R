#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot correlation matrix of gene features vs. eigenfeatures


options(stringsAsFactors=F, scipen=1000, family="sans")


######################
### DATA FUNCTIONS ###
######################
# Compute correlation matrix between two sets of features
make.cor.mat <- function(f1, f2, method="spearman"){
  f1.members <- colnames(f1)[-1]
  f2.members <- colnames(f2)[-1]
  f12 <- merge(f1, f2, by="gene", all=F, sort=F)
  cor.mat <- cor(f12[, -1], method=method)
  cor.mat[which(rownames(cor.mat) %in% f1.members),
          which(colnames(cor.mat) %in% f2.members)]
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Add colored heatmap cells to a preexisting plot()
add.heat.cells <- function(mat, pal, x.start=0, y.start=0, lwd=1, outer.lwd=1){
  # Note: assumes cells of 1 unit x 1 unit
  # Note: value of mat must match dimensions of pal
  n.rows <- nrow(mat)
  n.cols <- ncol(mat)
  sapply(1:n.cols, function(i){
    rect(xleft=rep(x.start+i-1, n.rows), xright=rep(x.start+i, n.rows),
         ybottom=y.start+(1:n.rows)-1, ytop=y.start+(1:n.rows),
         col=pal[mat[, i]], border="white", lwd=lwd)
  })
  rect(xleft=x.start, xright=x.start+n.cols,
       ybottom=y.start, ytop=y.start+n.rows,
       col=NA, border=blueblack, lwd=outer.lwd, xpd=T)
}

# Heatmap of feature correlation coefficients
plot.feat.cor <- function(feat2feat.mat, feat2eigen.mat, meta, pal,
                          dist.method="euclidean", pdf.scalar=1, spacer.ex=0.015, 
                          feat.label.line=-0.7, parmar=c(0.5, 16, 16, 0.5)){
  # Get plot dimensions
  n.cols <- nrow(feat2feat.mat)
  n.rows.top <- ncol(feat2feat.mat)
  n.rows.bottom <- ncol(feat2eigen.mat)
  n.rows <- n.rows.top + n.rows.bottom
  spacer.rows <- spacer.ex*n.rows
  top.y0 <- 0
  bottom.y0 <- n.rows.top + spacer.rows
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(-2, n.cols), ylim=c(n.rows + spacer.rows, -2),
       xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")
  
  # Transpose matrixes and hierarchically cluster columns
  feat2feat.mat <- t(feat2feat.mat)
  feat2eigen.mat <- t(feat2eigen.mat)
  column.names <- sort(colnames(feat2feat.mat))
  merged.mat <- rbind(feat2feat.mat[, column.names], feat2eigen.mat[, column.names])
  col.order <- column.names[hclust(dist(t(merged.mat), method=dist.method))$order]
  feat.row.names <- unlist(sapply(sort(unique(meta$category)), 
                                  function(c){mixedsort(meta$name[which(meta$category==c & meta$feature %in% column.names)])}))
  feat.row.order <- as.character(sapply(feat.row.names, function(n){meta$feature[which(meta$name==n)]}))
  feat2feat.mat <- feat2feat.mat[feat.row.order, col.order]
  feat2eigen.mat <- feat2eigen.mat[, col.order]
  
  # Scale matrixes for plotting
  feat2feat.pmat <- round(100 * feat2feat.mat) + 101
  feat2eigen.pmat <- round(100 * feat2eigen.mat) + 101
  
  # Add heatmap cells
  add.heat.cells(feat2feat.pmat, pal, x.start=0, y.start=0, 
                 lwd=pdf.scalar/2, outer.lwd=pdf.scalar)
  add.heat.cells(feat2eigen.pmat, pal, x.start=0, y.start=n.rows.top+spacer.rows, 
                 lwd=pdf.scalar/2, outer.lwd=pdf.scalar)
  
  # Add column labels
  sapply(1:n.cols, function(i){
    axis(3, at=i-0.5, labels=meta$name[which(meta$feature==col.order[i])], 
         las=2, tick=F, line=feat.label.line, cex=5/6)
    rect(xleft=i-1, xright=i, ybottom=-1, ytop=-2,
         col=gene.feat.category.colors[meta$category[which(meta$feature==col.order[i])]],
         border=blueblack, lwd=pdf.scalar, xpd=T)
  })
  
  # Add row labels
  sapply(1:n.rows.top, function(i){
    axis(2, at=i-0.5, labels=meta$name[which(meta$feature==feat.row.order[i])],
         las=2, tick=F, line=feat.label.line, cex=5/6)
    rect(xleft=-2, xright=-1, ybottom=i-1, ytop=i,
         col=gene.feat.category.colors[meta$category[which(meta$feature==feat.row.order[i])]],
         border=blueblack, lwd=pdf.scalar, xpd=T)
  })
  sapply(1:n.rows.bottom, function(i){
    axis(2, at=n.rows.top+spacer.rows+i-0.5, 
         labels=meta$name[which(meta$feature==rownames(feat2eigen.mat)[i])],
         las=2, tick=F, line=feat.label.line, cex=5/6)
    rect(xleft=-2, xright=-1, ybottom=n.rows.top+spacer.rows+i-1, 
         ytop=n.rows.top+spacer.rows+i, col="white", border=blueblack, lwd=pdf.scalar, xpd=T)
  })
  
  # Add cleanup lines
  axis(3, at=c(0, n.cols), tck=0, col=blueblack, labels=NA, lwd=pdf.scalar)
  axis(2, at=c(0, n.rows.top), tck=0, col=blueblack, labels=NA, lwd=pdf.scalar)
  axis(2, at=c(n.rows.top+spacer.rows, n.rows+spacer.rows+1), tck=0, col=blueblack, labels=NA, lwd=pdf.scalar)
}

# Plot correlation coefficient scale for heatmap
plot.feat.cor.scale <- function(mpal){
  # Prep plot area
  par(mar=c(0.25, 0.5, 2.2, 0.5), bty="n")
  plot(NA, xlim=c(0, 201), ylim=c(0, 1),
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  
  # Add gradient
  rect(xleft=0:200, xright=1:201, ybottom=0, ytop=1,
       col=pal, border=pal)
  rect(xleft=0, xright=201, ybottom=0, ytop=1,
       col=NA, border=blueblack, xpd=T)
  
  # Add labels
  top.at <- c(0, 100.5, 201)
  axis(3, at=top.at, tck=-0.15, labels=NA, col=blueblack)
  axis(3, at=top.at, tick=F, line=-0.85, labels=c(-1, 0, 1))
  mtext(3, line=1, text=bquote("Spearman's" ~ rho))
}

# Plot legend for gene feature categories
plot.feat.cor.legend <- function(pt.cex=2){
  n.cats <- length(gene.feat.category.colors)
  
  # Prep plot area
  par(mar=c(0.25, 0.25, 0.25, 5), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, n.cats),
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  
  # Add points
  points(x=rep(0.5, n.cats), y=(1:n.cats)-0.5, pch=22, cex=pt.cex,
         col=blueblack, bg=gene.feat.category.colors)
  
  # Add labels
  sapply(1:n.cats, function(i){
    axis(4, at=i-0.5, las=2, tick=F, line=-1,
         labels=str_to_title(names(gene.feat.category.colors)[i]))
  })
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(gtools, quietly=T)
require(stringr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog features.bed eigenfeatures.bed metadata.tsv out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: features.bed, eigenfeatures.bed, feature_metadata.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
features.in <- args$args[1]
eigen.in <- args$args[2]
feature.metadata.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz"
# eigen.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz"
# feature.metadata.in <- "~/scratch/gene_feature_metadata.tsv"
# out.prefix <- "~/scratch/gene_feature_corplot"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_association/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load features
features <- load.features(features.in, fill="mean", norm=F)
n.features <- ncol(features) - 1
feat.meta <- load.gene.feature.metadata(feature.metadata.in)

# Remove features based on PCs or std dev
features <- features[, -grep("_component_", colnames(features), fixed=T)]
features <- features[, -grep("_sd$", colnames(features))]

# Load eigenfeatures
eigen <- load.features(eigen.in, fill="mean", norm=F)
n.eigen <- ncol(eigen) - 1
colnames(eigen)[-1] <- paste("eigenfeature", 1:n.eigen, sep="_")
eigen.meta <- as.data.frame(t(sapply(1:n.eigen, function(i){
  c(paste(c("eigenfeature_", "Eigenfeature "), i, sep=""), NA)
})))
colnames(eigen.meta) <- colnames(feat.meta)

# Subset to top 10 eigenfeatures
eigen <- eigen[, 1:11]

# Merge feature metadata
all.meta <- rbind(feat.meta, eigen.meta)

# Compute Spearman correlation matrixes
feat2feat.mat <- cor(features[-1], method="spearman")
feat2eigen.mat <- make.cor.mat(features, eigen)


# Set color palette
# pal <- colorRampPalette(c(cnv.colors[1], control.cnv.colors[1], redwhite, 
#                           bluewhite, control.cnv.colors[2], cnv.colors[2]))(201)
pal <- colorRampPalette(c("#D42BD4", "white", "#22DD22"))(201)

# Plot heatmap
pdf.scalar <- 2
pdf(paste(out.prefix, "gene_feature_cor.heatmap.pdf", sep="."),
    height=pdf.scalar*8, width=pdf.scalar*7)
plot.feat.cor(feat2feat.mat, feat2eigen.mat, all.meta, pal,
              dist.method="manhattan", pdf.scalar=pdf.scalar,
              parmar=c(0.5, 14, 14, 0.5))
dev.off()

# Plot legends
pdf(paste(out.prefix, "gene_feature_cor.scale.pdf", sep="."),
    height=0.65, width=2.1)
plot.feat.cor.scale(pal)
dev.off()
pdf(paste(out.prefix, "gene_feature_cor.legend.pdf", sep="."),
    height=1.25, width=1.25)
plot.feat.cor.legend()
dev.off()
