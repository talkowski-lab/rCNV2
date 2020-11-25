#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distributions of gene features for dosage sensitive gene categories for rCNV2 final paper


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(beeswarm, quietly=T)
require(vioplot, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv features.bed feature.meta.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, features.bed, feature.meta.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
features.in <- args$args[2]
feature.meta.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz"
# feature.meta.in <- "~/scratch/gene_feature_metadata.tsv"
# out.prefix <- "~/scratch/test_gene_score_feature_regressions"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load scores & classify genes into subgroups based on scores
scores <- load.scores(scores.in)
ds.groups <- classify.genes(scores, hc.cutoff=0.9, lc.cutoff=0.5)

# Load features (leave unnormalized, but make some adjustments for interpretability)
feats <- load.features(features.in, norm=T)
feats <- mod.features(feats)
feat.meta <- load.gene.feature.metadata(feature.meta.in)

# Iterate over features and plot one comparison for each
sapply(2:ncol(feats), function(i){
  out.pdf <- paste(out.prefix, gsub("/", "", colnames(feats)[i], fixed=T),
                   "feature_distribs.pdf", sep=".")
  pdf(out.pdf, height=2.4, width=3)
  plot.feature.bydsgroup(feats, ds.groups, feat.idx=i, 
                         title=feat.meta$name[which(feat.meta$feature==colnames(feats)[i])],
                         swarm.max=1200)
  dev.off()
})

