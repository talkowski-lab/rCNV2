#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Run logistic regression comparisons of dosage sensitive vs tolerant genes


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Perform logistic ridge regression of gene features between two subgroups
gene.logit <- function(feats, zero.genes, one.genes){
  all.genes <- unique(c(zero.genes, one.genes))
  logit.df <- feats[which(feats$gene %in% all.genes), ]
  Y <- sapply(logit.df$gene, function(g){if(g %in% one.genes){1}else{0}})
  logit.df <- apply(logit.df[, -1], 2, function(vals){as.vector(as.numeric(unlist(vals)))})
  fit <- cv.glmnet(logit.df, Y, family="binomial", alpha=0.5, nfolds=10, type.measure="class")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(glmnet, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv features.bed out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: scores.tsv, features.bed, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
features.in <- args$args[2]
out.prefix <- args$args[3]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# features.in <- "~/scratch/refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz"
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

# Load & normalize features
feats <- load.features(features.in, norm=T)



