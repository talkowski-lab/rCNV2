#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Determine gene features predictive of dosage sensitivity with regression


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# # Perform logistic ridge regression of gene features between two subgroups
# gene.logit <- function(feats, zero.genes, one.genes, glm.alpha=0.5){
#   all.genes <- unique(c(zero.genes, one.genes))
#   logit.df <- feats[which(feats$gene %in% all.genes), ]
#   Y <- sapply(logit.df$gene, function(g){if(g %in% one.genes){1}else{0}})
#   logit.df <- apply(logit.df[, -1], 2, function(vals){as.vector(as.numeric(unlist(vals)))})
#   fit <- cv.glmnet(logit.df, Y, family="binomial", alpha=glm.alpha, nfolds=10, type.measure="class")
# }

# Perform linear elastic net regression of gene features versus a continuous response
gene.glm <- function(scores, feats, response.name, glm.alpha=0.5, seed=2020){
  glm.df <- feats[which(feats$gene %in% scores$gene), ]
  if(!(response.name %in% colnames(scores))){
    stop(paste("Error: response.name \"", response.name, "\" not found in scores."))
  }
  Y <- as.numeric(sapply(as.vector(glm.df$gene), function(g){
    scores[which(scores$gene==g), which(colnames(scores)==response.name)]
  }))
  Y <- as.vector(as.numeric(scale(Y)))
  glm.df <- apply(glm.df[, -1], 2, function(vals){as.vector(as.numeric(unlist(vals)))})
  # glm.df <- cbind(Y, glm.df)
  set.seed(seed)
  # model <- caret::train(Y ~ ., data=glm.df, method="glmnet", 
  #                       trControl=trainControl("cv", number=10),
  #                       tuneLength=10)
  # fit <- model$bestTune
  # coefs <- coef(model$finalModel, model$bestTune$lambda)
  # coefs <- coefs[order(coefs[, 1]), ]
  cv <- cv.glmnet(glm.df, Y, alpha=glm.alpha, nfolds=10)
  fit <- glmnet(glm.df, Y, alpha=glm.alpha, lambda=cv$lambda.1se)
  coefs <- coef(fit)
  coefs <- coefs[order(coefs[, 1]), ]
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(glmnet, quietly=T)
require(caret, quietly=T)

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
# features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz"
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

# Compute regression gradients
scores$ds.gradient <- apply(scores[, -1], 1, function(vals){
  min(as.numeric(vals), na.rm=T)
})
scores$hits.gradient <- scores$pTS - scores$pHI

# Load & normalize features
feats <- load.features(features.in, fill="mean", norm=T)



