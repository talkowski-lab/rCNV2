#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Common functions for gene association and fine-mapping analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load credible sets from BED
load.credsets <- function(credsets.in){
  # Read data & clean columns
  credsets <- read.table(credsets.in, header=T, sep="\t", comment.char="")
  colnames(credsets)[1] <- gsub("X.", "", colnames(credsets)[1], fixed=T)
  
  # Ensure numerics
  numeric.cols <- c("start", "end", "mean_control_freq", "mean_case_freq",
                    "pooled_ln_or", "pooled_ln_or_ci_lower", "pooled_ln_or_ci_upper",
                    "best_pvalue", "n_genes")
  credsets[, numeric.cols] <- apply(credsets[, numeric.cols], 2, as.numeric)
  
  # Split list-style columns
  credsets$all_genes <- strsplit(credsets$all_genes, split=";")
  credsets$vconf_genes <- strsplit(credsets$vconf_genes, split=";")
  credsets$conf_genes <- strsplit(credsets$conf_genes, split=";")
  
  return(credsets)
}

# Load all individual gene associations from BED
load.associations <- function(assocs.in){
  # Read data & clean columns
  assocs <- read.table(assocs.in, header=T, sep="\t", comment.char="")
  colnames(assocs)[1] <- gsub("X.", "", colnames(assocs)[1], fixed=T)
  
  # Ensure numerics
  numeric.cols <- c("start", "end", "control_freq", "case_freq",
                    "ln_or", "ln_or_ci_lower", "ln_or_ci_upper",
                    "pvalue", "pip")
  assocs[, numeric.cols] <- apply(assocs[, numeric.cols], 2, as.numeric)
  
  return(assocs)
}


##########################
### PLOTTING FUNCTIONS ###
##########################

