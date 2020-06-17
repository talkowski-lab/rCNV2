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
load.stats <- function(stats.in){
  stats <- read.table(stats.in, sep="\t", comment.char="", 
                      check.names=F, header=T)
  colnames(stats)[1] <- gsub("#", "", colnames(stats)[1], fixed=T)
  return(stats)
}

# Extract cohort names from an input stats file
extract.cohorts <- function(stats){
  as.vector(sapply(colnames(stats)[grep("_case_alt", colnames(stats), fixed=T)],
                   function(str){
                     unlist(strsplit(str, split="_"))[1]
                   }))
}

# Apply empirical continuity correction to meta-analysis data frame
# Per Sweeting et al., Stat. Med., 2004 (section 3.3)
sweeting.correction <- function(meta.df, cc.sum=0.01){
  # Count number of carriers & non-carriers
  n.alt <- sum(meta.df[, grep("_alt", colnames(meta.df), fixed=T)])
  n.ref <- sum(meta.df[, grep("_ref", colnames(meta.df), fixed=T)])
  # Require at least one CNV to be observed
  if(n.alt>0){
    nt <- n.alt
    R <- n.ref/n.alt
    # Pooled odds ratio estimate of all non-zero studies with at least one case sample
    nonzero.studies <- intersect(which(apply(meta.df[, grep("_alt", colnames(meta.df), fixed=T)], 1, sum)>0),
                                 which(apply(meta.df[, grep("case_", colnames(meta.df), fixed=T)], 1, sum)>0))
    nonzero.case.odds <- sum(meta.df$case_alt[nonzero.studies])/sum(meta.df$case_ref[nonzero.studies])
    nonzero.control.odds <- sum(meta.df$control_alt[nonzero.studies])/sum(meta.df$control_ref[nonzero.studies])
    if(!is.nan(nonzero.case.odds) & !is.nan(nonzero.control.odds)){
      if(nonzero.control.odds>0){
        ohat <- nonzero.case.odds/nonzero.control.odds
        # Otherwise, apply standard continuity correction of 0.5 to pooled estimate if no CNVs observed in controls
      }else{
        nonzero.case.odds <- (sum(meta.df$case_alt[nonzero.studies])+0.5)/(sum(meta.df$case_ref[nonzero.studies])+0.5)
        nonzero.control.odds <- (sum(meta.df$control_alt[nonzero.studies])+0.5)/(sum(meta.df$control_ref[nonzero.studies])+0.5)
        ohat <- nonzero.case.odds/nonzero.control.odds
      }
      # Otherwise, apply standard continuity correction to pooled estimate of *all* studies
    }else{
      nonzero.case.odds <- (sum(meta.df$case_alt)+0.5)/(sum(meta.df$case_ref)+0.5)
      nonzero.control.odds <- (sum(meta.df$control_alt)+0.5)/(sum(meta.df$control_ref)+0.5)
      ohat <- nonzero.case.odds/nonzero.control.odds
    }
    
    # Solve for kc & kt
    kc <- R/(R+ohat)
    kt <- ohat/(R+ohat)
    # Compute continuity corrections
    cor.case_alt <- cc.sum * kt
    cor.case_ref <- cc.sum * kc
    cor.control_alt <- cc.sum * (nt + kt)
    cor.control_ref <- cc.sum * ((nt*R) + kc)
    # Apply continuity corrections
    meta.df$case_alt <- meta.df$case_alt + cor.case_alt
    meta.df$case_ref <- meta.df$case_ref + cor.case_ref
    meta.df$control_alt <- meta.df$control_alt + cor.control_alt
    meta.df$control_ref <- meta.df$control_ref + cor.control_ref
  }
  return(meta.df)
}

# Make meta-analysis data frame for a single track
make.meta.df <- function(stats, cohorts, row.idx, empirical.continuity=T){
  ncohorts <- length(cohorts)
  meta.df <- data.frame("cohort"=1:ncohorts,
                        "control_ref"=as.numeric(stats[row.idx, grep("control_ref", colnames(stats), fixed=T)]),
                        "case_ref"=as.numeric(stats[row.idx, grep("case_ref", colnames(stats), fixed=T)]),
                        "control_alt"=as.numeric(stats[row.idx, grep("control_alt", colnames(stats), fixed=T)]),
                        "case_alt"=as.numeric(stats[row.idx, grep("case_alt", colnames(stats), fixed=T)]),
                        "cohort_name"=cohorts)
  if(empirical.continuity==T){
    meta.df <- sweeting.correction(meta.df)
  }
  return(meta.df)
}

# Perform meta-analysis for a single track
meta.single <- function(stats, cohorts, row.idx, model="fe", empirical.continuity=T){
  # If all CNVs are ref or all are alt, return all NAs
  if(all(sum(stats[row.idx, grep("_alt", colnames(stats), fixed=T)])>0,
         sum(stats[row.idx, grep("_ref", colnames(stats), fixed=T)])>0)){
    meta.df <- make.meta.df(stats, cohorts, row.idx, empirical.continuity)
    # If strictly zero case CNVs are observed, unable to estimate effect size
    if(all(meta.df$case_alt==0)){
      out.v <- c(rep(NA, 4), 0)
    }else{
      # Meta-analysis
      if(model=="re"){
        meta.res <- tryCatch(rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                     measure="OR", data=meta.df, method="REML", random = ~ 1 | cohort, slab=cohort_name,
                                     add=0, drop00=F, correct=F, digits=5, control=list(maxiter=100, stepadj=0.5)),
                             error=function(e){
                               print(paste("row", row.idx, "failed to converge. Retrying with more iterations...", sep=" "))
                               rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                       measure="OR", data=meta.df, method="REML", random = ~ 1 | cohort, slab=cohort_name,
                                       add=0, drop00=F, correct=F, digits=5, control=list(maxiter=10000, stepadj=0.4))
                             })
        out.v <- as.numeric(c(meta.res$b[1,1], meta.res$ci.lb, meta.res$ci.ub,
                              meta.res$zval, -log10(meta.res$pval)))
      }else if(model=="mh"){
        meta.res <- rma.mh(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                           measure="OR", data=meta.df, slab=cohort_name,
                           add=0, drop00=F, correct=F)
        out.v <- as.numeric(c(meta.res$b, meta.res$ci.lb, meta.res$ci.ub,
                              meta.res$zval, -log10(meta.res$MHp)))
      }else if(model=="fe"){
        meta.res <- tryCatch(rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                     measure="OR", data=meta.df, method="FE", slab=cohort_name,
                                     add=0, drop00=F, correct=F, digits=5, control=list(maxiter=100, stepadj=0.5)),
                             error=function(e){
                               print(paste("row", row.idx, "failed to converge. Retrying with more iterations...", sep=" "))
                               rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                       measure="OR", data=meta.df, method="FE", slab=cohort_name,
                                       add=0, drop00=F, correct=F, digits=5, control=list(maxiter=10000, stepadj=0.4))
                             })
        out.v <- as.numeric(c(meta.res$b[1,1], meta.res$ci.lb, meta.res$ci.ub,
                              meta.res$zval, -log10(meta.res$pval)))
      }
    }
    return(out.v)
  }else{
    rep(NA, 5)
  }
}

# Make meta-analysis lookup table to shorten time required to run full meta-analysis
make.meta.lookup.table <- function(stats, cohorts, model, empirical.continuity=T){
  unique.counts.df <- unique(stats[, sort(unique(c(grep("_ref", colnames(stats), fixed=T),
                                                   grep("_alt", colnames(stats), fixed=T))))])
  
  unique.stats <- t(sapply(1:nrow(unique.counts.df), function(i){
    meta.single(unique.counts.df, cohorts, i, model, empirical.continuity)
  }))
  
  lookup.table <- cbind(unique.counts.df, unique.stats)
  stat.colnames <- c("meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper", "meta_z", "meta_phred_p")
  colnames(lookup.table)[(ncol(lookup.table)-4):ncol(lookup.table)] <- stat.colnames
  
  return(lookup.table)
}

# Apply saddlepoint approximation to vector of Z-scores to generate adjusted P-values
saddlepoint.adj <- function(zscores, phred=T){
  mu.hat <- mean(zscores, na.rm=T)
  sd.hat <- sd(zscores, na.rm=T)
  cumuls <- gaussianCumulants(mu.hat, sd.hat)
  dx <- 0.01
  x <- seq(-40, 40, dx)
  saddle.pdf <- saddlepoint(x, 1, cumuls)$approx
  saddle.cdf <- cumsum(saddle.pdf * 0.01)
  calc.saddle.p <- function(z){if(!is.na(z)){1 - tail(saddle.cdf[which(x<z)], 1)}else{NA}}
  new.pvals <- sapply(zscores, calc.saddle.p)
  if(phred==T){
    return(-log10(new.pvals))
  }else{
    return(new.pvals)
  }
}

# Wrapper function to perform a meta-analysis on all tracks
meta <- function(stats, cohorts, model="fe", saddle=T){
  # Make meta-analysis lookup table
  meta.lookup.table <- make.meta.lookup.table(stats, cohorts, model, 
                                              empirical.continuity=T)
  
  # Merge stats into full list
  meta.res <- merge(stats, meta.lookup.table, sort=F, all.x=T, all.y=F)
  meta.res <- meta.res[with(meta.res, order(trackname)), ]
  
  # Adjust P-values using saddlepoint approximation of null distribution, if optioned
  if(saddle==T){
    meta.res$meta_phred_p <- saddlepoint.adj(meta.res$meta_z)
  }
  
  # Add FDR-adjusted q-value
  meta.res$meta_phred_fdr_q <- -log10(p.adjust(10^-meta.res$meta_phred_p, method="fdr"))

  # Format output
  return(meta.res[, -unique(c(grep("_ref", colnames(meta.res), fixed=T),
                              grep("_alt", colnames(meta.res), fixed=T)))])
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
  make_option(c("--model"), type="character", default="fe", 
              help="specify meta-analysis model ('re': random effects, 'fe': fixed effects, 'mh': Mantel-Haenszel) [default '%default']",
              metavar="string"),
  make_option(c("--spa"), action="store_true", default=FALSE, 
              help="apply saddlepoint approximation of null distribution [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog infile outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
stats.in <- args$args[1]
outfile <- args$args[2]
model <- opts$model
spa <- opts$spa

# # Dev parameters
# stats.in <- "~/scratch/all_tracks.stats.with_counts.tsv.gz"
# outfile <- "~/scratch/all_tracks.burden_stats.tsv.gz"
# model <- "fe"
# spa <- T

# Read track stats
stats <- load.stats(stats.in)
cohorts <- extract.cohorts(stats)

# Run meta-analysis for all tracks
meta.res <- meta(stats, cohorts, model, saddle=spa)
colnames(meta.res)[1] <- paste("#", colnames(meta.res)[1], sep="")
write.table(meta.res, outfile, sep="\t",
            row.names=F, col.names=T, quote=F)
