#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Utility functions used for gene association & fine-mapping analyses


#' Load credible sets
#'
#' Load summary dataframe for credible set association stats
#'
#' @param credsets.in path to BED file of credible sets
#'
#' @return data.frame
#'
#' @export
load.credsets <- function(credsets.in){
  # Read data & clean columns
  credsets <- read.table(credsets.in, header=T, sep="\t", comment.char="")
  colnames(credsets)[1] <- gsub("X.", "", colnames(credsets)[1], fixed=T)

  # Ensure numerics
  numeric.cols <- c("start", "end", "mean_control_freq", "mean_case_freq",
                    "pooled_ln_or", "pooled_ln_or_ci_lower", "pooled_ln_or_ci_upper",
                    "best_pvalue", "n_genes", "n_sig_genes", "n_hpos")
  credsets[, numeric.cols] <- apply(credsets[, numeric.cols], 2, as.numeric)

  # Get pointwise parameters
  pw.params <- get.sig.gene.pw.params(credsets$best_sig_level, credsets$cnv)
  credsets$pt.bg <- pw.params$pt.bg
  credsets$pt.border <- pw.params$pt.border

  # Split list-style columns
  list.cols <- c("all_genes", "sig_genes", "top_gene",
                 "vconf_genes", "conf_genes", "hpos")
  credsets[, list.cols] <- apply(credsets[, list.cols], 2, strsplit, split=";")

  # Compute number of significant genes in credible set
  credsets$n_sig_in_credset <- apply(credsets[, c("all_genes", "sig_genes")],
                                     1, function(glists){
    length(glists$sig_genes %in% glists$all_genes)
  })

  return(credsets)
}


#' Load gene associations
#'
#' Load all gene association stats from a BED file
#'
#' @param assocs.in path to BED file of gene-level association stats
#'
#' @return data.frame
#'
#' @export
load.gene.associations <- function(assocs.in){
  # Read data & clean columns
  assocs <- read.table(assocs.in, header=T, sep="\t", comment.char="")
  colnames(assocs)[1] <- gsub("X.", "", colnames(assocs)[1], fixed=T)

  # Ensure numerics
  numeric.cols <- c("start", "end", "control_freq", "case_freq",
                    "ln_or", "ln_or_ci_lower", "ln_or_ci_upper",
                    "pvalue", "pip_final", "pip_HPO_specific")
  assocs[, numeric.cols] <- apply(assocs[, numeric.cols], 2, as.numeric)

  # Convert boolean columns
  assocs$top_gene <- as.logical(assocs$top_gene)

  # Get pointwise parameters
  pw.params <- get.sig.gene.pw.params(assocs$sig_level, assocs$cnv)
  assocs$pt.bg <- pw.params$pt.bg
  assocs$pt.border <- pw.params$pt.border

  return(assocs)
}


#' Label gene associations by significance
#'
#' Generate pointwise graphical annotations for gene associations by significance level
#'
#' @param sig character vector indicating significance level
#' @param cnv character vector of CNV type
#'
#' @return list of pt.bg and pt.border
#'
#' @export
get.sig.gene.pw.params <- function(sig, cnv){
  pt.bg <- rep(NA, length(sig))
  pt.border <- rep(NA, length(sig))

  gw.idxs <- which(sig == "exome_wide")
  pt.bg[gw.idxs] <- cnv.colors[cnv[gw.idxs]]
  pt.border[gw.idxs] <- cnv.blacks[cnv[gw.idxs]]

  fdr.idxs <- which(sig == "FDR")
  pt.bg[fdr.idxs] <- control.cnv.colors[cnv[fdr.idxs]]
  pt.border[fdr.idxs] <- cnv.colors[cnv[fdr.idxs]]

  return(list("pt.bg"=pt.bg, "pt.border"=pt.border))
}


#' Categorize fine-mapped genes
#'
#' Categorize genes based on their fine-mapped PIP and their rank in each credible set
#'
#' @param credsets data.frame of credible sets as read by [load.credsets()]
#'
#' @return list of gene categories split by CNV type
#'
#' @seealso [load.credsets()]
#'
#' @export
categorize.genes <- function(credsets){
  g.cats <- lapply(c("DEL", "DUP", "CNV"), function(cnv){
    if(cnv=="CNV"){
      cnv.idxs <- 1:nrow(credsets)
    }else{
      cnv.idxs <- which(credsets$cnv==cnv)
    }
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
  names(g.cats) <- c("DEL", "DUP", "CNV")
  return(g.cats)
}


#' Summarize categories of fine-mapped genes
#'
#' Summarize fine-mapping categories based on PIP & rank in each credible set
#'
#' @param g.cats categories of genes split by PIP & rank as produced by [categorize.genes()]
#' @param cnv CNV type to evaluate \[default: "CNV"\]
#'
#' @export summarize.finemapping.categories
#' @export
summarize.finemapping.categories <- function(g.cats, cnv="CNV"){
  # Extract groups for specified CNV type
  groups <- g.cats[[cnv]]
  cat(paste("Very confident genes & top in credible set =",
            prettyNum(length(groups$top.vconf)), "\n"))
  cat(paste("Very confident genes, not top in credible set =",
            prettyNum(length(groups$nottop.vconf)), "\n"))
  cat(paste("Confident genes & top in credible set =",
            prettyNum(length(groups$top.conf)), "\n"))
  cat(paste("Confident genes, not top in credible set =",
            prettyNum(length(groups$nottop.conf)), "\n"))
  cat(paste("Unlikely genes & top in credible set =",
            prettyNum(length(groups$top.notconf)), "\n"))
  cat(paste("Unlikely genes, not top in credible set =",
            prettyNum(length(groups$nottop.notconf)), "\n"))
}


#' Subset finemapped genes
#'
#' Get all genes for a single quadrant of 2x2 grid of fine-mapped genes
#'
#' @param gene.groups categories of genes split by PIP & rank as produced by [categorize.genes()]
#' @param top rank of genes in credible set
#' @param conf confidence level for subset
#'
#' @return vector of gene symbols
#'
#' @export get.quadrant.genes
#' @export
get.quadrant.genes <- function(gene.groups, top, conf){
  group.name <- paste(top, conf, sep=".")
  sort(unique(c(gene.groups$DEL[[group.name]],
                gene.groups$DUP[[group.name]])))
}

