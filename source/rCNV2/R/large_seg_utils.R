#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Utility functions used for large segment analyses


#' Load locus-level association stats
#'
#' Load summary dataframe for locus-level association stats
#'
#' @param loci.in path to loci BED file
#'
#' @return data.frame
#'
#' @export
load.loci <- function(loci.in){
  # Read data
  loci <- read.table(loci.in, sep="\t", header=T, comment.char="")
  colnames(loci)[1] <- "chr"

  # Split list-style columns
  loci$hpos <- strsplit(loci$hpos, split=";")
  loci$constituent_assocs <- strsplit(loci$constituent_assocs, split=";")
  loci$cred_interval_coords <- strsplit(loci$cred_interval_coords, split=";")
  loci$genes <- strsplit(loci$genes, split=";")

  # Convert numeric columns to numerics
  numeric.cols <- c("start_min", "end_max", "pooled_control_freq", "pooled_case_freq",
                    "pooled_ln_or", "pooled_ln_or_ci_lower", "pooled_ln_or_ci_upper",
                    "min_ln_or", "max_ln_or", "n_hpos", "n_constituent_assocs",
                    "n_cred_intervals", "cred_intervals_size", "n_genes")
  for(col in numeric.cols){
    cidx <- which(colnames(loci) == col)
    loci[, cidx] <- as.numeric(loci[, cidx])
  }

  # Add formatted locus names and sizes
  loci$formatted_name <- sapply(strsplit(loci$region_id, split="_"), function(parts){parts[4]})
  loci$formatted_size <- paste(prettyNum(round(loci$cred_intervals_size/1000, 0), big.mark=","), "kb", sep=" ")

  return(loci)
}


#' Load main segment table
#'
#' Load summary table for all segments (including lit. GDs & NAHR regions)
#'
#' @param segs.in path to segments BED file
#'
#' @return data.frame
#'
#' @export
load.segment.table <- function(segs.in){
  # Read data
  segs <- read.table(segs.in, header=T, sep="\t", comment.char="")
  colnames(segs)[1] <- "chr"

  # Split list-style columns
  listcol.idxs <- c(which(colnames(segs) %in% c("coords", "genes")),
                    intersect(grep("^n_", colnames(segs), invert=T, fixed=F),
                              grep("_genes$", colnames(segs), fixed=F)))
  segs[, listcol.idxs] <- apply(segs[, listcol.idxs], 2, strsplit, split=";")

  # Convert numeric columns to numerics
  numcol.idxs <- unique(c(which(colnames(segs) %in% c("start", "end", "size", "meta_best_p", "Redin_BCAs")),
                          grep("^n_", colnames(segs), fixed=F),
                          grep("_dnm_", colnames(segs), fixed=T),
                          grep("_express", colnames(segs), fixed=T)))
  segs[, numcol.idxs] <- apply(segs[, numcol.idxs], 2, as.numeric)

  # Convert boolean dummy columns to logicals
  boolcol.idxs <- which(colnames(segs) %in% c("any_sig", "gw_sig", "fdr_sig",
                                              "nom_sig", "nom_neuro", "nom_dev",
                                              "hc_gd", "mc_gd", "lc_gd", "any_gd",
                                              "pathogenic", "benign", "nahr",
                                              "pleiotropic"))
  segs[, boolcol.idxs] <- apply(segs[, boolcol.idxs], 2, function(vals){
    sapply(vals, function(val){if(val==1){TRUE}else{FALSE}})})

  # Add normalized columns
  segs$gnomAD_constrained_prop <- segs$n_gnomAD_constrained_genes / segs$n_genes
  segs$prop_ubiquitously_expressed <- segs$n_ubiquitously_expressed_genes / segs$n_genes
  segs <- normalize.dnms(segs)
  segs <- normalize.dnms(segs, dnm.cohorts=c("DDD_noSig", "ASC_noSig", "ASC_unaffected_noSig"))

  # Add formatted sizes
  segs$formatted_size <- paste(prettyNum(round(segs$size/1000, 0), big.mark=","), "kb", sep=" ")

  # Add graphical columns
  gw.sig.idx <- which(segs$gw_sig)
  fdr.sig.idx <- which(segs$fdr_sig)
  nomsig.gd.idx <- which(segs$any_gd & !segs$any_sig & segs$nom_sig)
  nonsig.gd.idx <- which(segs$any_gd & !segs$nom_sig)
  segs$pt.pch <- NA; segs$pt.border <- NA; segs$pt.bg <- NA
  segs$pt.pch[gw.sig.idx] <- 22
  segs$pt.border[gw.sig.idx] <- cnv.blacks[segs$cnv[gw.sig.idx]]
  segs$pt.bg[gw.sig.idx] <- cnv.colors[segs$cnv[gw.sig.idx]]
  segs$pt.pch[fdr.sig.idx] <- 23
  segs$pt.border[fdr.sig.idx] <- cnv.colors[segs$cnv[fdr.sig.idx]]
  segs$pt.bg[fdr.sig.idx] <- control.cnv.colors[segs$cnv[fdr.sig.idx]]
  segs$pt.pch[nomsig.gd.idx] <- 21
  segs$pt.border[nomsig.gd.idx] <- control.cnv.colors[segs$cnv[nomsig.gd.idx]]
  segs$pt.bg[nomsig.gd.idx] <- cnv.whites[segs$cnv[nomsig.gd.idx]]
  segs$pt.pch[nonsig.gd.idx] <- 21
  segs$pt.border[nonsig.gd.idx] <- control.cnv.colors[segs$cnv[nonsig.gd.idx]]
  segs$pt.bg[nonsig.gd.idx] <- NA

  # Return cleaned dataframe
  return(segs)
}


#' Get list of neuro-associated loci
#'
#' Compile master list of neuro-associated loci
#'
#' @param loci locus association dataframe (imported with [load.loci()])
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param sig.only only return genome-wide or FDR significant segments \[default: FALSE\]
#'
#' @return vector
#'
#' @seealso [load.loci()], [load.segment.table()]
#'
#' @export get.neuro.region_ids
#' @export
get.neuro.region_ids <- function(loci, segs, sig.only=FALSE){
  sig.hits <- loci$region_id[which(sapply(loci$hpos, function(hpos){any(hpos %in% neuro.hpos)}))]
  lit.hits <- segs$region_id[which(segs$any_gd & segs$nom_neuro & !segs$any_sig)]
  if(sig.only){
    sort(unique(sig.hits))
  }else{
    sort(unique(c(sig.hits, lit.hits)))
  }
}


#' Get list of developmental or adult loci
#'
#' Compile master list of loci associated with developmental or adult phenotypes
#'
#' @param loci locus association dataframe (imported with [load.loci()])
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param adult return loci only associated with adult phenotypes \[default: FALSE\]
#' @param sig.only only return genome-wide or FDR significant segments \[default: FALSE\]
#'
#' @return vector
#'
#' @seealso [load.loci()], [load.segment.table()]
#'
#' @export get.developmental.region_ids
#' @export
get.developmental.region_ids <- function(loci, segs, adult=FALSE, sig.only=FALSE){
  if(adult){
    sig.hits <- loci$region_id[which(sapply(loci$hpos, function(hpos){!any(hpos %in% developmental.hpos)}))]
    lit.hits <- segs$region_id[which(segs$any_gd & !segs$nom_dev & !segs$any_sig)]
  }else{
    sig.hits <- loci$region_id[which(sapply(loci$hpos, function(hpos){any(hpos %in% developmental.hpos)}))]
    lit.hits <- segs$region_id[which(segs$any_gd & segs$nom_dev & !segs$any_sig)]
  }
  if(sig.only){
    sort(unique(sig.hits))
  }else{
    sort(unique(c(sig.hits, lit.hits)))
  }
}


#' Get list of neurodevelopmental loci
#'
#' Compile master list of loci associated with the intersection of neurological and developmental phenotypes
#'
#' @param loci locus association dataframe (imported with [load.loci()])
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param sig.only only return genome-wide or FDR significant segments \[default: FALSE\]
#'
#' @return vector
#'
#' @seealso [load.loci()], [load.segment.table()]
#'
#' @export get.ndd.region_ids
#' @export
get.ndd.region_ids <- function(loci, segs, sig.only=FALSE){
  ndd.hpos <- intersect(developmental.hpos, neuro.hpos)
  sig.hits <- loci$region_id[which(sapply(loci$hpos, function(hpos){any(hpos %in% ndd.hpos)}))]
  lit.hits <- segs$region_id[which(segs$any_gd & segs$nom_dev & segs$nom_neuro & !segs$any_sig)]
  if(sig.only){
    sort(unique(sig.hits))
  }else{
    sort(unique(c(sig.hits, lit.hits)))
  }
}


#' Normalize DNM counts
#'
#' Normalize segment DNM counts vs. synonymous inflation
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param dnm.cohorts vector of cohorts to evaluate (default: DDD, ASC, ASC_unaffected)
#'
#' @return data.frame
#'
#' @export
normalize.dnms <- function(segs, dnm.cohorts=c("DDD", "ASC", "ASC_unaffected")){
  for(cohort in dnm.cohorts){
    syn.obs.colname <- paste(cohort, "dnm_syn_obs_wMu", sep="_")
    syn.exp.colname <- paste(cohort, "dnm_syn_exp_wMu", sep="_")
    if(syn.obs.colname %in% colnames(segs) & syn.exp.colname %in% colnames(segs)){
      syn.obs <- segs[, which(colnames(segs)==syn.obs.colname)]
      syn.exp <- segs[, which(colnames(segs)==syn.exp.colname)]
      syn.d <- (syn.obs - syn.exp) / segs$n_genes
      for(csq in c("lof", "mis")){
        dam.obs.colname <- paste(cohort, "dnm", csq, "obs_wMu", sep="_")
        dam.exp.colname <- paste(cohort, "dnm", csq, "exp_wMu", sep="_")
        if(dam.obs.colname %in% colnames(segs) & dam.exp.colname %in% colnames(segs)){
          dam.obs <- segs[, which(colnames(segs)==dam.obs.colname)]
          dam.exp <- segs[, which(colnames(segs)==dam.exp.colname)]
          dam.d <- (dam.obs - dam.exp) / segs$n_genes
          # Outlier-robust linear fit of obs-exp damaging ~ obs-exp synonymous
          fit <- robust.lm(syn.d, dam.d)$fit
          coeffs <- as.numeric(fit$coefficients)
          dam.d.adj <- dam.d - coeffs[1] - (coeffs[2] * syn.d)
          segs[paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_")] <- dam.d.adj
        }
      }
    }
  }
  return(segs)
}


#' Merge locus and segment files
#'
#' Merge locus association data with master segments table
#'
#' @param loci locus association dataframe (imported with [load.loci()])
#' @param segs segment dataframe (imported with [load.segment.table()])
#'
#' @return vector
#'
#' @seealso [load.loci()], [load.segment.table()]
#'
#' @export merge.loci.segs
#' @export
merge.loci.segs <- function(loci, segs){
  shared.columns <- intersect(colnames(loci), colnames(segs))
  shared.to.drop.from.segs <- shared.columns[which(shared.columns != "region_id")]
  segs.sub <- segs[, -which(colnames(segs) %in% shared.to.drop.from.segs)]
  gw <- merge(loci, segs.sub, by="region_id", all.x=T, all.y=F,
              sort=F, suffixes=c(".l", ".s"))
  gw$color <- cnv.colors[sapply(gw$cnv, function(cnv){which(names(cnv.colors)==cnv)})]
  gw$black <- cnv.blacks[sapply(gw$cnv, function(cnv){which(names(cnv.blacks)==cnv)})]
  return(gw)
}

#' Check for genomic disorder overlap
#'
#' Lookup overlap with any known GDs for a set of coordinates
#'
#' @param chrom chromosome
#' @param start start coordinate
#' @param end end coordinate
#' @param segs segment dataframe (imported with [load.segment.table()])
#'
#' @return vector of GD CNV types
#'
#' @export get.gd.overlap
#' @export
get.gd.overlap <- function(chrom, start, end, segs){
  sort(unique(segs$cnv[which(segs$any_gd
                             & segs$chr==chrom
                             & segs$end>=start
                             & segs$start<=end)]))
}

#' Summarize segment permutation results
#'
#' Summarize permutation results across all permutations for a single feature
#'
#' @param perms list of permutation results
#' @param feature name of feature to evaluate
#' @param measure statistic to evaluate \[default: mean\]
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#'
#' @details Valid options for `measure` are mean, median, sum, and frac.any
#'
#' @return vector of statistic split by all, deletions, and duplications
#'
#' @export perm.summary
#' @export
perm.summary <- function(perms, feature, measure="mean", subset_to_regions=NULL){
  if(!is.null(subset_to_regions)){
    perms <- lapply(perms, function(df){df[which(df$region_id %in% subset_to_regions), ]})
  }
  do.call("rbind", lapply(perms, function(df){
    if(measure == "mean"){
      c("ALL" = mean(df[, which(colnames(df)==feature)], na.rm=T),
        "DEL" = mean(df[which(df$cnv=="DEL"),
                        which(colnames(df)==feature)], na.rm=T),
        "DUP" = mean(df[which(df$cnv=="DUP"),
                        which(colnames(df)==feature)], na.rm=T))
    }else if(measure == "median"){
      c("ALL" = median(df[, which(colnames(df)==feature)], na.rm=T),
        "DEL" = median(df[which(df$cnv=="DEL"),
                          which(colnames(df)==feature)], na.rm=T),
        "DUP" = median(df[which(df$cnv=="DUP"),
                          which(colnames(df)==feature)], na.rm=T))
    }else if(measure == "sum"){
      c("ALL" = sum(df[, which(colnames(df)==feature)]),
        "DEL" = sum(df[which(df$cnv=="DEL"),
                       which(colnames(df)==feature)]),
        "DUP" = sum(df[which(df$cnv=="DUP"),
                       which(colnames(df)==feature)]))
    }else if(measure == "frac.any"){
      hits <- which(df[, which(colnames(df)==feature)] > 0)
      100 * c("ALL" = length(hits) / nrow(df),
              "DEL" = length(intersect(hits, which(df$cnv=="DEL"))) / length(which(df$cnv=="DEL")),
              "DUP" = length(intersect(hits, which(df$cnv=="DUP"))) / length(which(df$cnv=="DUP")))
    }
  }))
}

#' Summarize a segment metric
#'
#' Helper function to compute various statistics across a subset of segments
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param feature name of feature to evaluate
#' @param measure statistic to evaluate \[default: mean\]
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#'
#' @details Valid options for `measure` are mean, median, sum, and frac.any
#'
#' @return vector of statistic split by all, deletions, and duplications
#'
#' @export calc.segs.dat
#' @export
calc.segs.dat <- function(segs, feature, measure, subset_to_regions=NULL){
  if(!is.null(subset_to_regions)){
    segs <- segs[which(segs$region_id %in% subset_to_regions), ]
  }
  # Convert boolean to numeric, if needed
  segs.vals <- segs[, which(colnames(segs)==feature)]
  if(all(sapply(segs.vals, function(x){is.logical(x) | is.na(x)}))){
    segs.vals <- sapply(segs.vals, function(x){if(is.na(x)){NA}else if(x==T){1}else if(x==F){0}else{NA}})
  }
  if(measure == "mean"){
    c("ALL" = mean(segs.vals, na.rm=T),
      "DEL" = mean(segs.vals[which(segs$cnv=="DEL")], na.rm=T),
      "DUP" = mean(segs.vals[which(segs$cnv=="DUP")], na.rm=T))
  }else if(measure == "median"){
    c("ALL" = median(segs.vals, na.rm=T),
      "DEL" = median(segs.vals[which(segs$cnv=="DEL")], na.rm=T),
      "DUP" = median(segs.vals[which(segs$cnv=="DUP")], na.rm=T))
  }else if(measure == "sum"){
    c("ALL" = sum(segs.vals, na.rm=T),
      "DEL" = sum(segs.vals[which(segs$cnv=="DEL")], na.rm=T),
      "DUP" = sum(segs.vals[which(segs$cnv=="DUP")], na.rm=T))
  }else if(measure == "frac.any"){
    100 * c("ALL" = length(which(segs.vals > 0)) / length(segs.vals),
            "DEL" = length(which(segs.vals[which(segs$cnv=="DEL")] > 0)) / length(which(segs$cnv=="DEL")),
            "DUP" = length(which(segs.vals[which(segs$cnv=="DUP")] > 0)) / length(which(segs$cnv=="DUP")))
  }
}

