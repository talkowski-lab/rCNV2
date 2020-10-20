#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Common functions for noncoding association analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Extract summary statistics from a single gene or CRB BED
load.sumstats <- function(stats.in){
  # Read data & drop coordinates
  stats <- read.table(stats.in, header=T, sep="\t", check.names=F, comment.char="")[, -c(1:3)]
  
  # Clean data types
  stats[, -c(1, 4, 6)] <- apply(stats[, -c(1, 4, 6)], 2, as.numeric)
  
  return(stats)
}

# Assign noncoding annotation family based on annotation track & source
assign.family <- function(track, source){
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

# Load burden test statistics for all tracks
load.track.stats <- function(stats.in){
  # Read data
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")
  colnames(stats)[1] <- gsub("X.", "", colnames(stats)[1], fixed=T)
  
  # Partition feature source name from track name, and assign family
  stats$source <- sapply(stats$trackname, function(str){unlist(strsplit(str, split=".", fixed=T))[1]})
  stats$trackname <- sapply(stats$trackname, function(str){paste(unlist(strsplit(str, split=".", fixed=T))[-1], collapse=".")})
  stats$family <- apply(stats[, c("trackname", "source")], 1, function(vals){assign.family(vals[1], vals[2])})
  
  # Drop unnecessary columns
  stats <- stats[, -c(grep("meta[1-4]", colnames(stats)), 
                      which(colnames(stats) %in% c("original_path")))]
  
  return(stats)
}

# Load BED of CRBs
load.crbs <- function(crbs.in){
  # Read data
  crbs <- read.table(crbs.in, header=T, sep="\t", comment.char="")
  colnames(crbs)[1] <- gsub("X.", "", colnames(crbs)[1], fixed=T)
  
  # Compute size
  crbs$size <- crbs$end - crbs$start
  
  return(crbs)
}


##########################
### PLOTTING FUNCTIONS ###
##########################

