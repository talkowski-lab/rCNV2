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
# Load table of feature metadata (including plain English labels)
load.feature.metadata <- function(feature.metadata.in){
  meta <- read.table(feature.metadata.in, header=T, sep="\t", comment.char="")
  colnames(meta)[1] <- gsub("^X.", "", colnames(meta)[1])
  return(meta)
}

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
# Heatmap of feature correlation coefficients
plot.cor.mat <- function(cor.mat, row.meta, col.meta){
  n.rows <- nrow(cor.mat)
  n.cols <- ncol(cor.mat)
  pal <- colorRampPalette(c(cnv.colors[1], control.cnv.colors[1], redwhite, 
                            bluewhite, control.cnv.colors[2], cnv.colors[2]))(201)
  plot(NA, xlim=c(n.cols, 0), ylim=c(n.rows, 0))
  sapply(1:n.rows, function(r){
    colors <- pal[round(100 * cor.mat[r, ]) + 101]
    rect(xleft=0:(n.cols-1), xright=1:n.cols,
         ybottom=r-1, ytop=r,
         col=colors, border=colors, lwd=0.2)
  })
}



#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

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
feat.meta <- load.feature.metadata(feature.metadata.in)

# Load eigenfeatures
eigen <- load.features(eigen.in, fill="mean", norm=F)
n.eigen <- ncol(eigen) - 1
colnames(eigen)[-1] <- paste("eigenfeature", 1:n.eigen, sep="_")
eigen.meta <- as.data.frame(t(sapply(1:n.eigen, function(i){
  c(paste(c("eigenfeature_", "Eigenfeature "), i, sep=""), NA)
})))
colnames(eigen.meta) <- colnames(feat.meta)

# Compute Spearman correlation matrix
cor.mat <- make.cor.mat(features, eigen)
  
  