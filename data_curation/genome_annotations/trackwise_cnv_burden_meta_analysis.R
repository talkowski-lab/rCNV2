#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Burden meta-analysis of noncoding rCNV counts across annotation tracks


options(scipen=1000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Load track stats 
load.stats <- function(stats.in, cnv.split=T){
  stats <- read.table(stats.in, sep="\t", comment.char="", 
                      check.names=F, header=T)
  colnames(stats)[1] <- gsub("#", "", colnames(stats)[1], fixed=T)
  if(cnv.split==T){
    list("DEL" = stats[which(stats$cnv=="DEL"), ],
         "DUP" = stats[which(stats$cnv=="DUP"), ])
  }else{
    return(stats)
  }
}

# Extract cohort names from an input stats file
extract.cohorts <- function(stats){
  as.vector(sapply(colnames(stats)[grep("_case_alt", colnames(stats), fixed=T)],
                   function(str){
                     unlist(strsplit(str, split="_"))[1]
                   }))
}

# Apply saddlepoint approximation to vector of Z-scores to generate adjusted P-values and Z-scores
saddlepoint.adj <- function(zscores, phred=T, alternative="two.sided"){
  mu.hat <- mean(zscores, na.rm=T)
  sd.hat <- sd(zscores, na.rm=T)
  cumuls <- gaussianCumulants(mu.hat, sd.hat)
  dx <- 0.01
  x <- seq(min(min(zscores, na.rm=T), -40), 
           max(max(zscores, na.rm=T), 40), dx)
  saddle.pdf <- saddlepoint(x, 1, cumuls)$approx
  saddle.cdf <- cumsum(saddle.pdf * 0.01)
  calc.saddle.p.z <- function(z){
    if(!is.na(z)){
      if(alternative=="greater"){
        p <- 1 - tail(saddle.cdf[which(x<z)], 1)
        new.z <- qnorm(1 - p)
      }else if(alternative=="less"){
        p <- tail(saddle.cdf[which(x<z)], 1)
        new.z <- qnorm(p)
      }else{
          p <- tail(saddle.cdf[which(x<z)], 1)
          new.z <- qnorm(p)
          if(p>0.5){
            p <- 1 - p
          }
      }
    }else{
      p <- NA
      new.z <- NA
    }
    return(c(new.z, p))
  }
  new.stats <- t(sapply(zscores, calc.saddle.p.z))
  colnames(new.stats) <- c("zscore", "p")
  if(phred==T){
    new.stats[, 2] <- -log10(new.stats[, 2])
    colnames(new.stats)[2] <- "phred_p"
  }
  return(as.data.frame(new.stats))
}

# Calculate odds ratios, variance, and Z-scores per cohort (with SPA, if optioned)
calc.ors <- function(stats, cohorts, spa=T){
  or.df <- as.data.frame(do.call("cbind", lapply(cohorts, function(cohort){
    ors <- as.data.frame(t(sapply(1:nrow(stats), function(i){
      as.numeric(metafor::escalc("OR", 
                      ai=stats[i, paste(cohort, "control", "ref", sep="_")],
                      bi=stats[i, paste(cohort, "case", "ref", sep="_")],
                      ci=stats[i, paste(cohort, "control", "alt", sep="_")],
                      di=stats[i, paste(cohort, "case", "alt", sep="_")]))
    })))
    colnames(ors) <- c("lnOR", "var")
    # Z-score computed against per-category variance and empirical mean across all categories
    # Assumes the average category should _not_ have an effect, so any systematic deviance from
    # mean = 0 is due to filtering biases
    ors$zscore <- (ors$lnOR - mean(ors$lnOR, na.rm=T)) / sqrt(ors$var)
    # Adjust Z-score with SPA, if optioned
    if(spa==T){
      adj.z.p <- saddlepoint.adj(ors$zscore, alternative="two.sided", phred=T)
      ors$zscore <- adj.z.p$zscore
      ors$phred_p <- adj.z.p$phred_p
    }
    colnames(ors) <- paste(cohort, colnames(ors), sep=".")
    return(ors)
  })))
  return(as.data.frame(cbind(stats, or.df)))
}

# Meta analysis of Z-scores weighted by inverse variance
weighted.z <- function(stats, cohorts, spa=T){
  meta.z <- sapply(1:nrow(stats), function(i){
    z <- stats[i, grep("zscore", colnames(stats), fixed=T)]
    var <- stats[i, grep("var", colnames(stats), fixed=T)]
    inv.var <- 1/var
    weighted.mean(z, inv.var, na.rm=T)
  })
  if(spa==T){
    meta.z.p <- saddlepoint.adj(meta.z, phred=T)
  }else{
    meta.z.p <- data.frame("zscore" = meta.z,
                           "phred_p" = -log10(pnorm(abs(meta.z), lower.tail=F)))
  }
  meta.z.p$lnOR <-  sapply(1:nrow(stats), function(i){
    lnor <- stats[i, grep("lnOR", colnames(stats), fixed=T)]
    var <- stats[i, grep("var", colnames(stats), fixed=T)]
    inv.var <- 1/var
    weighted.mean(lnor, inv.var, na.rm=T)
  })
  meta.z.p$phred_fdr_q <- p.adjust(10^-meta.z.p$phred_p, method="fdr")
  colnames(meta.z.p) <- paste("meta", colnames(meta.z.p), sep=".")
  as.data.frame(cbind(stats, meta.z.p))
}


################
###RSCRIPT BLOCK
################
# Load required libraries
require(optparse, quietly=T)
require(metafor, quietly=T)
require(EQL, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--p-cutoff"), type="numeric", default=0.05, 
              help="Meta-analysis P-value cutoff to consider a track significant [default '%default']",
              metavar="numeric"),
  make_option(c("--signif-tracks"), type="character",  
              help="output file for significant tracks", metavar="string")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog infile tracklist outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
stats.in <- args$args[1]
tracklist.in <- args$args[2]
outfile <- args$args[3]
p.cutoff <- -log10(as.numeric(opts$`p-cutoff`))
signif.outfile <- opts$`signif-tracks`

# # Dev parameters
# stats.in <- "~/scratch/rCNV.all.merged_stats.with_counts.tsv.gz"
# tracklist.in <- "~/scratch/rCNV.all_tracks.list"
# outfile <- "~/scratch/rCNV.chromhmm_plus_encode.test.tsv"
# p.cutoff <- -log10(0.05)
# signif.outfile <- "~/scratch/rCNV.chromhmm_plus_encode.test.signif_tracks.list"

# Read full tracklist
tracklist <- read.table(tracklist.in, header=F, sep="\t", comment.char="")[, 1]

# Read track stats and split by CNV type
stats <- load.stats(stats.in, cnv.split=T)
cohorts <- extract.cohorts(stats[[1]])

# Calculate case:control ORs, variance, and Z-scores for each track
stats <- lapply(stats, calc.ors, cohorts, spa=T)

# Weighted Z-score meta-analysis
meta.res.split <- lapply(stats, weighted.z, cohorts)

# Combine DEL & DUP, and write to outfile
meta.res <- as.data.frame(do.call("rbind", meta.res.split))
colnames(meta.res)[1] <- paste("#", colnames(meta.res)[1], sep="")
write.table(meta.res, outfile, sep="\t",
            row.names=F, col.names=T, quote=F)

# Extract significant track names
if(!is.null(signif.outfile)){
  sig.idx <- which(meta.res$meta.phred_p >= p.cutoff & meta.res$meta.zscore > 0)
  if(length(sig.idx) > 0){
    sig.names <- as.character(sort(unique(meta.res[sig.idx, 1])))
    sig.tracks <- as.data.frame(do.call("rbind", lapply(sig.names, function(name){
      c(tracklist[grep(name, tracklist, fixed=T)], name)})))
    write.table(sig.tracks, signif.outfile,
                col.names=F, row.names=F, quote=F, sep="\t")
  }
}
