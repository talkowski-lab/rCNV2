#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Utility functions used for noncoding association analyses


#' Load CRB sumstats
#'
#' Extract summary statistics from a single gene or CRB BED
#'
#' @param stats.in path to .bed file of CRB association stats
#'
#' @export
load.crb.sumstats <- function(stats.in){
  # Read data & drop coordinates
  stats <- read.table(stats.in, header=T, sep="\t", check.names=F, comment.char="")[, -c(1:3)]

  # Clean data types
  stats[, -c(1, 4, 6)] <- apply(stats[, -c(1, 4, 6)], 2, as.numeric)

  return(stats)
}


#' Assign noncoding annotation families
#'
#' Assign noncoding annotation family based on annotation track & source
#'
#' @param track track name
#' @param source name of track source
#'
#' @return string of track family
#'
#' @export
assign.track.family <- function(track, source){
  mappings <- c("boca" = "dhs",
                "dbSUPER" = "super.enhancers",
                "encode_dnaaccessibility" = "dhs",
                "encode_histone_mods" = "histone",
                "encode_tad_boundaries" = "tads",
                "encode_tfbs" = "tfbs",
                "encode_transcription" = "transcription",
                "EnhAtlas2" = "enhancers",
                "fantom_enh" = "enhancers",
                "hacer_enh" = "enhancers",
                "jaspar" = "tfbs",
                "roadmap_chromhmm" = "chromhmm",
                "UCNEbase" = "other")
  if(source %in% names(mappings)){
    return(mappings[source])
  }else if(any(sapply(c("SEdb_SE.", "SEA."), grepl, x=track, fixed=T))){
    return("super.enhancers")
  }else if(any(sapply(c("SEdb_TE.", "DENdb.", "plac_", "vista_", "enhancers", "_pREs"), grepl, x=track, fixed=T))){
    return("enhancers")
  }else if(any(sapply(c("rao_2014", "dixon_2012", "gao_2016", "tad_boundaries"), grepl, x=track, fixed=T))){
    return("tads")
  }else if(any(sapply(c("gencode", "all_TARs"), grepl, x=track, fixed=T))){
    return("transcription")
  }else if(grepl("_OCRs", track, fixed=T)){
    return("dhs")
  }else if(grepl("H3K27ac_peaks", track, fixed=T)){
    return("histone")
  }else{
    return(NA)
  }
}


#' Load track burden sumstats
#'
#' Load burden test statistics for all tracks
#'
#' @param stats.in path to track burden summary statistics file
#' @param spa update P-values with saddlepoint re-approximation of the null \[default: FALSE\]
#'
#' @return data.frame
#'
#' @export
load.track.stats <- function(stats.in, spa=F){
  # Read data
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")
  colnames(stats)[1] <- gsub("X.", "", colnames(stats)[1], fixed=T)

  # Update P-values with saddlepoint reapproximation, if optioned
  # Do this separately for DEL and DUP
  if(spa){
    for(cnv in c("DEL", "DUP")){
      require(EQL, quietly=T)
      cnv.idxs <- which(stats$cnv == cnv)
      spa.res <- saddlepoint.adj(stats$meta.zscore[cnv.idxs])
      stats$meta.zscore[cnv.idxs] <- spa.res$zscores
      # Must recompute P-values for Z<0 because saddlepoint.adj only does one-tailed tests
      stats$meta.neglog10_p[cnv.idxs] <- sapply(spa.res$zscores, function(z){
        if(z >= 0){
          -pnorm(z, lower.tail=F, log.p=T)/log(10)
        }else{
          -pnorm(z, lower.tail=T, log.p=T)/log(10)
        }
      })
      stats$meta.neglog10_fdr_q[cnv.idxs] <- -log10(p.adjust(spa.res$pvalues, method="fdr"))
    }
  }

  # Partition feature source name from track name, and assign family
  stats$source <- sapply(stats$trackname, function(str){unlist(strsplit(str, split=".", fixed=T))[1]})
  stats$trackname <- sapply(stats$trackname, function(str){paste(unlist(strsplit(str, split=".", fixed=T))[-1], collapse=".")})
  stats$family <- apply(stats[, c("trackname", "source")], 1, function(vals){assign.track.family(vals[1], vals[2])})

  # Drop unnecessary columns
  stats <- stats[, -c(grep("meta[0-9]+", colnames(stats)),
                      which(colnames(stats) %in% c("original_path")))]

  return(stats)
}


#' Load CRBs
#'
#' Load BED of CRBs
#'
#' @param crbs.in path to CRB .bed file
#'
#' @return data.frame of CRBs
#'
#' @export
load.crbs <- function(crbs.in){
  # Read data
  crbs <- read.table(crbs.in, header=T, sep="\t", comment.char="")
  colnames(crbs)[1] <- gsub("X.", "", colnames(crbs)[1], fixed=T)

  # Compute size
  crbs$size <- crbs$end - crbs$start

  return(crbs)
}


