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

# Categorize genes based on their rank in credible set & PIP
categorize.genes <- function(credsets){
  g.cats <- lapply(c("DEL", "DUP"), function(cnv){
    cnv.idxs <- which(credsets$cnv==cnv)
    g.top <- unique(sort(credsets$top_gene[cnv.idxs]))
    g.nottop <- setdiff(unique(sort(unlist(credsets$all_genes[cnv.idxs]))),
                        g.top)
    g.vconf <- unique(sort(unlist(credsets$vconf_genes[cnv.idxs])))
    g.conf <- setdiff(unique(sort(unlist(credsets$conf_genes[cnv.idxs]))),
                      g.vconf)
    g.notconf <- setdiff(setdiff(unique(sort(unlist(credsets$all_genes[cnv.idxs]))),
                                 g.vconf), g.conf)
    g.top.vconf <- intersect(g.top, g.vconf)
    g.nottop.vconf <- intersect(g.nottop, g.vconf)
    g.top.conf <- intersect(g.top, g.conf)
    g.nottop.conf <- intersect(g.nottop, g.conf)
    g.top.notconf <- intersect(g.top, g.notconf)
    g.nottop.notconf <- intersect(g.nottop, g.notconf)
    list("top.vconf"=g.top.vconf, "nottop.vconf"=g.nottop.vconf,
         "top.conf"=g.top.conf, "nottop.conf"=g.nottop.conf,
         "top.notconf"=g.top.notconf, "nottop.notconf"=g.nottop.notconf)
  })
  names(g.cats) <- c("DEL", "DUP")
  return(g.cats)
}

# Get all genes for a single quadrant of 2x2 grid
get.quadrant.genes <- function(gene.groups, top, conf){
  group.name <- paste(top, conf, sep=".")
  sort(unique(c(gene.groups$DEL[[group.name]],
                gene.groups$DUP[[group.name]])))
}


##########################
### PLOTTING FUNCTIONS ###
##########################

