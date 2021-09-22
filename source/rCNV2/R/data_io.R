#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Functions for handling I/O of various preformatted datasets


#' Load BED3
#'
#' Load a simple BED3 file
#'
#' @param bed.in path to input BED fil
#' @param header boolean indicator if BED has header \[default: TRUE\]
#'
#' @return data frame of coordinates
#'
#' @export load.bed3
#' @export
load.bed3 <- function(bed.in, header=T){
  bed <- read.table(bed.in, header=header, sep="\t", comment.char="", check.names=F)
  colnames(bed) <- c("chr", "start", "end")
  bed[, -1] <- apply(bed[, -1], 2, as.numeric)
  return(bed)
}


#' Load single-cohort association statistics
#'
#' Load association statistics for a single cohort from an input file
#'
#' @param stats.in path to input BED file with association statistics
#' @param prefix cohort name (to be appended to columns)
#' @param p.is.phred boolean indicator of the P-value being -log10-scaled in `stats.in`
#' @param keep.n.cols number of columns from original BED format to retain
#' @param keep.or.confint boolean indicator to retain confidence intervals for
#' odds ratio estimates \[default: only keep point estimate\]
#'
#' @return data frame of formatted association stats
#'
#' @export
read.assoc.stats.single <- function(stats.in, prefix, p.is.phred, keep.n.cols=3,
                                    keep.or.confint=FALSE){
  # Read data & subset to necessary columns
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")
  colnames(stats)[1] <- "chr"
  cols.to.keep <- c(colnames(stats)[1:keep.n.cols], "case_alt", "case_ref",
                    "control_alt", "control_ref", "fisher_phred_p",
                    "fisher_OR", "fisher_OR_lower", "fisher_OR_upper")
  stats <- stats[, which(colnames(stats) %in% cols.to.keep)]
  colnames(stats)[which(colnames(stats)=="fisher_phred_p")] <- "p_value"
  if(p.is.phred==T){
    stats$p_value <- 10^-stats$p_value
  }
  colnames(stats)[-(1:keep.n.cols)] <- paste(prefix, colnames(stats)[-(1:keep.n.cols)], sep=".")
  return(stats)
}


#' Load meta-analysis association statistics
#'
#' Load association statistics from a meta-analysis
#'
#' @param stats.in path to input BED file with association statistics
#' @param keep.n.cols number of columns from original BED format to retain
#'
#' @return data frame of formatted association stats
#'
#' @export
load.meta.stats <- function(stats.in, p.is.phred, keep.n.cols=3){
  # Read data & subset to necessary columns
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(stats)[1] <- gsub("#", "", colnames(stats)[1], fixed=T)
  mandatory.cols.idx <- which(colnames(stats) == "n_nominal_cohorts"):ncol(stats)
  if(keep.n.cols < 1){
    cols.to.keep <- colnames(stats)[mandatory.cols.idx]
  }else{
    cols.to.keep <- colnames(stats)[unique(sort(c(1:keep.n.cols, mandatory.cols.idx)))]
  }
  stats <- stats[, which(colnames(stats) %in% cols.to.keep)]
  return(stats)
}


#' Load P-value matrix
#'
#' Loads a precomputed matrix of P-values
#'
#' @param matrix.in path to matrix.tsv
#' @param has.coords boolean indicator if `matrix.in` is BED3+ formatted \[default: `TRUE`\]
#' @param ncols.coords number of leading columns to treat as coordinates (only
#' used if `has.coords = TRUE`) \[default: `3`\]
#' @param p.is.phred boolean indicator if P-values in matrix are `-log10`-scaled
#'
#' @details This function was designed to process the output from genome-wide
#' permutation analyses to parameterize false-discovery rates
#'
#' @return list of four objects:
#' 1. `$coords` : data frame of coordinates (`NULL` if `has.coords = FALSE`)
#' 2. `$pvals` : data frame of p-values
#' 3. `$expected` : expected P-value distribution under true null
#' 4. `$lambdas` : genomic inflation statistic for each column in `$pvals`
#'
#' @export
load.pval.matrix <- function(matrix.in, has.coords=T, ncols.coords=3, p.is.phred=T){
  x <- read.table(matrix.in, header=T, sep="\t", comment.char="")
  if(has.coords == T){
    colnames(x)[1:ncols.coords] <- c("chrom", "start", "end", "gene")[1:ncols.coords]
    x[, 1] <- as.character(x[, 1])
    x[, 2:3] <- apply(x[, 2:3], 2, as.numeric)
    x[, -c(1:ncols.coords)] <- apply(x[, -c(1:ncols.coords)], 2, as.numeric)
    coords <- as.data.frame(x[, 1:ncols.coords])
    pvals <- as.data.frame(x[, -c(1:ncols.coords)])
  }else{
    coords <- NULL
    pvals <- as.data.frame(apply(x, 2, as.numeric))
  }
  if(p.is.phred == T){
    pvals <- as.data.frame(apply(pvals, 2, function(x){10^-x}))
  }
  expected.full <- ppoints(nrow(pvals))
  lambdas <- apply(pvals, 2, function(obs){
    expected <- ppoints(length(obs[which(!is.na(obs))]))
    dchisq(median(obs, na.rm=T), df=1)/dchisq(median(expected), df=1)
  })
  return(list("coords" = coords, "pvals" = pvals,
              "expected" = expected.full, "lambdas" = lambdas))
}



#' Load gene features BED
#'
#' Load a BED4+ formatted .tsv of gene features
#'
#' @param features.in path to BED4+ file with feature columns (see `Details`)
#' @param fill specify fill behavior for missing values. See `Details` for
#' options. \[default: NA\]
#' @param norm boolean indicator whether standard normalization should be
#' applied to features upon loading. \[default: FALSE\]
#'
#' @details `features.in` must be BED4+, with column 4 corresponding to gene symbol
#'
#' `fill` accepts the following values:
#' * `mean` : fill missing values with column-wise mean
#' * `median` : fill missing values with column-wise median
#' * `{numeric}` : fill missing values with a user-specified numeric constant
#' * `NA` : do not fill missing values
#'
#' @return data frame of features
#' @export
load.features <- function(features.in, fill=NA, norm=F){
  feats <- read.table(features.in, header=T, sep="\t", comment.char="", check.names=F)[, -c(1:3)]
  feats[, -1] <- apply(feats[, -1], 2, as.numeric)
  if(!is.na(fill)){
    feats[, -1] <- apply(feats[, -1], 2, function(vals){
      na.idxs <- which(is.na(vals) | is.infinite(vals))
      if(length(na.idxs) > 0){
        if(fill=="mean"){
          vfill <- mean(vals[-na.idxs])
        }else if(fill=="median"){
          vfill <- median(vals[-na.idxs])
        }else if(is.numeric(fill)){
          vfill <- fill
        }
        vals[na.idxs] <- vfill
      }
      return(vals)
    })
  }
  if(norm==T){
    feats[, -1] <- apply(feats[, -1], 2, function(vals){
      scale(vals, scale=T, center=T)
    })
  }
  return(feats)
}


#' Load gene feature metadata
#'
#' Load table of gene feature metadata (including plain-text descriptions)
#'
#' @param feature.metadata.in path to gene feature metadata .tsv
#'
#' @return data.frame of contents from `feature.metadata.in`
#' @export
load.gene.feature.metadata <- function(feature.metadata.in){
  meta <- read.table(feature.metadata.in, header=T, sep="\t", comment.char="")
  colnames(meta)[1] <- gsub("^X.", "", colnames(meta)[1])
  return(meta)
}


#' Load gene list
#'
#' Loads a single list of genes from a flat text file
#'
#' @param path path to gene list
#'
#' @return character vector of gene symbols
#'
#' @export load.genelist
#' @export
load.genelist <- function(path){
  if(!is.null(path)){
    sort(unique(as.character(read.table(path, header=F)[, 1])))
  }else{
    NULL
  }
}

