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

  # Rename columns with typos
  colnames(segs) <- sapply(colnames(segs), gsub, pattern="CLinGen",
                           replacement="ClinGen", fixed=T)

  # Add columns for union of ClinGen & DECIPHER gene sets
  segs$clinical_LoF_genes <- sapply(1:nrow(segs), function(i){
    u <- union(unlist(segs$ClinGen_HI_genes[i]), unlist(segs$DECIPHER_LoF_genes[i]))
    list(u[which(!is.na(u))])
  })
  segs$n_clinical_LoF_genes <- sapply(segs$clinical_LoF_genes, function(gstr){length(unlist(gstr))})
  segs$clinical_GoF_genes <- sapply(1:nrow(segs), function(i){
    u <- union(unlist(segs$ClinGen_TS_genes[i]), unlist(segs$DECIPHER_GoF_genes[i]))
    list(u[which(!is.na(u))])
  })
  segs$n_clinical_GoF_genes <- sapply(segs$clinical_GoF_genes, function(gstr){length(unlist(gstr))})

  # Convert numeric columns to numerics
  numcol.idxs <- unique(c(which(colnames(segs) %in% c("start", "end", "size", "meta_best_p",
                                                      "min_LOEUF", "min_MisOEUF",
                                                      "total_LoF_OE", "total_mis_OE",
                                                      "Redin_BCAs")),
                          grep("^n_", colnames(segs), fixed=F),
                          grep("_dnm_", colnames(segs), fixed=T),
                          grep("_express", colnames(segs), fixed=T)))
  segs[, numcol.idxs] <- apply(segs[, numcol.idxs], 2, as.numeric)

  # Convert boolean dummy columns to logicals
  boolcol.idxs <- which(colnames(segs) %in% c("any_sig", "gw_sig", "fdr_sig", "bonf_sig_gd",
                                              "nom_sig", "nom_neuro", "nom_dev",
                                              "hc_gd", "mc_gd", "lc_gd", "any_gd",
                                              "pathogenic", "benign", "nahr",
                                              "terminal", "pleiotropic"))
  segs[, boolcol.idxs] <- apply(segs[, boolcol.idxs], 2, function(vals){
    sapply(vals, function(val){if(val==1){TRUE}else{FALSE}})})

  # Extract cytoband column from region identifier
  segs[, "cytoband"] <- unlist(sapply(segs$region_id, function(region.id){
    lab.parts <- unlist(strsplit(region.id, split="_"))
    lab.parts <- lab.parts[grep("[A-Z]", lab.parts, invert=T)]
    lab.parts[length(lab.parts)]
  }))

  # Add normalized columns
  segs$gnomAD_constrained_prop <- segs$n_gnomAD_constrained_genes / segs$n_genes
  segs$prop_ubiquitously_expressed <- segs$n_ubiquitously_expressed_genes / segs$n_genes
  segs <- normalize.dnms(segs)
  segs <- normalize.dnms(segs, dnm.cohorts=c("DDD_noSig", "ASC_noSig", "ASC_unaffected_noSig"))

  # Add columns for combined DNM excess between DDD and ASC
  dnm.base.cols <- sapply(colnames(segs)[grep("^DDD_", colnames(segs))],
                          gsub, pattern="^DDD_", replacement="")
  for(col.base in dnm.base.cols){
    new.col <- paste("DDD_plus_ASC", col.base, sep="_")
    segs[, new.col] <- apply(segs[, paste(c("DDD", "ASC"), col.base, sep="_")], 1, sum)
  }

  # Add formatted sizes
  segs$formatted_size <- paste(prettyNum(round(segs$size/1000, 0), big.mark=","), "kb", sep=" ")

  # Add graphical columns
  gw.sig.idx <- which(segs$gw_sig)
  fdr.sig.idx <- which(segs$fdr_sig)
  bonfsig.gd.idx <- which(segs$any_gd & !segs$any_sig & segs$bonf_sig_gd)
  nonsig.gd.idx <- which(segs$any_gd & !segs$bonf_sig_gd)
  segs$pt.pch <- NA; segs$pt.border <- NA; segs$pt.bg <- NA
  segs$pt.pch[gw.sig.idx] <- 22
  segs$pt.border[gw.sig.idx] <- cnv.blacks[segs$cnv[gw.sig.idx]]
  segs$pt.bg[gw.sig.idx] <- cnv.colors[segs$cnv[gw.sig.idx]]
  segs$pt.pch[fdr.sig.idx] <- 23
  segs$pt.border[fdr.sig.idx] <- cnv.colors[segs$cnv[fdr.sig.idx]]
  segs$pt.bg[fdr.sig.idx] <- control.cnv.colors[segs$cnv[fdr.sig.idx]]
  segs$pt.pch[bonfsig.gd.idx] <- 21
  segs$pt.border[bonfsig.gd.idx] <- control.cnv.colors[segs$cnv[bonfsig.gd.idx]]
  segs$pt.bg[bonfsig.gd.idx] <- cnv.whites[segs$cnv[bonfsig.gd.idx]]
  segs$pt.pch[nonsig.gd.idx] <- 13
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


#' Split segments by effect size
#'
#' Split segments into groups based on effect size quantiles
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param quantiles number of quantiles to divide \[default: 2\]
#' @param effect.size.colname name of column to reference for effect sizes
#' \[default: `meta_best_lnor`\]
#'
#' @return list of length `quantiles` with region IDs for each quantile
#'
#' @details uses values from `effect.size.colname` \(see above\) while
#' matching on CNV type to compute quantile cutoffs
#'
#' @export split.regions.by.effect.size
#' @export
split.regions.by.effect.size <- function(segs, quantiles=2,
                                         effect.size.colname="meta_best_lnor"){
  # Compute quantiles per CNV type
  quants <- lapply(c("DEL", "DUP"), function(cnv){
    cutoffs <- as.numeric(quantile(segs[which(segs$cnv==cnv), effect.size.colname],
                                   probs=seq(0, 1, length.out=quantiles+1)))
    cutoffs[c(1, length(cutoffs))] <- c(-Inf, Inf)
    return(cutoffs)
  })
  names(quants) <- c("DEL", "DUP")

  lapply(1:quantiles, function(i){
    unique(as.character(unlist(sapply(c("DEL", "DUP"), function(cnv){
      segs$region_id[which(segs$cnv==cnv
                           & segs$meta_best_lnor >= quants[[cnv]][i]
                           & segs$meta_best_lnor < quants[[cnv]][i+1])]
    }))))
  })
}


#' Normalize DNM counts
#'
#' Normalize segment DNM counts vs. expected and (optional) synonymous inflation
#'
#' @param segs segment dataframe (imported with [load.segment.table()])
#' @param dnm.cohorts vector of cohorts to evaluate (default: DDD, ASC, ASC_unaffected)
#' @param synonymous.adjustment indicator to perform outlier-robust linear regression
#' for inflation in synonymous counts \[default: TRUE\]
#' @param is.data.table indicator that `segs` is a data.table \[default: FALSE\]
#'
#' @return data.frame
#'
#' @export
normalize.dnms <- function(segs, dnm.cohorts=c("DDD", "ASC", "ASC_unaffected"),
                           synonymous.adjustment=TRUE, is.data.table=FALSE){
  for(cohort in dnm.cohorts){
    syn.obs.colname <- paste(cohort, "dnm_syn_obs_wMu", sep="_")
    syn.exp.colname <- paste(cohort, "dnm_syn_exp_wMu", sep="_")
    if(syn.obs.colname %in% colnames(segs) & syn.exp.colname %in% colnames(segs)){
      new.colname <- paste(cohort, "dnm_syn_raw_excess_per_gene", sep="_")
      if(is.data.table){
        segs[, eval(new.colname) := (get(syn.obs.colname) - get(syn.exp.colname)) / n_genes]
      }else{
        syn.obs <- segs[, which(colnames(segs)==syn.obs.colname)]
        syn.exp <- segs[, which(colnames(segs)==syn.exp.colname)]
        syn.d <- (syn.obs - syn.exp) / segs$n_genes
        segs[, new.colname] <- syn.d
      }
      for(csq in c("lof", "mis")){
        dam.obs.colname <- paste(cohort, "dnm", csq, "obs_wMu", sep="_")
        dam.exp.colname <- paste(cohort, "dnm", csq, "exp_wMu", sep="_")
        if(dam.obs.colname %in% colnames(segs) & dam.exp.colname %in% colnames(segs)){
          new.colname <- paste(cohort, "dnm", csq, "raw_excess_per_gene", sep="_")
          if(is.data.table){
            segs[, eval(new.colname) := (get(dam.obs.colname) - get(dam.exp.colname)) / n_genes]
          }else{
            dam.obs <- segs[, which(colnames(segs)==dam.obs.colname)]
            dam.exp <- segs[, which(colnames(segs)==dam.exp.colname)]
            dam.d <- (dam.obs - dam.exp) / segs$n_genes
            segs[, new.colname] <- dam.d
          }

          # Outlier-robust linear fit of obs-exp damaging ~ obs-exp synonymous
          if(synonymous.adjustment){
            new.colname <- paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_")
            if(is.data.table){
              syn.colname <- paste(cohort, "dnm_syn_raw_excess_per_gene", sep="_")
              dam.colname <- paste(cohort, "dnm", csq, "raw_excess_per_gene", sep="_")
              adj.dam.d <- function(dt){
                fit <- robust.lm(dt[, get(syn.colname)], dt[, get(dam.colname)])$fit
                coeffs <- as.numeric(fit$coefficients)
                dt[, get(dam.colname) - coeffs[1] - (coeffs[2] * get(syn.colname))]
              }
              segs[, (new.colname) := numeric()]
              sapply(1:max(segs$perm_idx), function(i){
                segs[perm_idx == i, (new.colname) := adj.dam.d(segs[perm_idx == i])]
              })
            }else{
              fit <- robust.lm(syn.d, dam.d)$fit
              coeffs <- as.numeric(fit$coefficients)
              dam.d.adj <- dam.d - coeffs[1] - (coeffs[2] * syn.d)
              segs[new.colname] <- dam.d.adj
            }
          }
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

#' Load gene-based permutation data
#'
#' Load all permuted segments matched on the number of genes
#'
#' @param perm.res.in path to .tsv of all permutation results
#' @param subset_to_regions optional vector of region IDs to keep \[default: keep all regions\]
#'
#' @return data.table() of permutation results
#'
#' @export
load.perms.bygene <- function(perm.res.in, subset_to_regions=NULL){
  # Read full permutation table
  perms <- data.table::fread(perm.res.in, header=T, sep="\t", check.names=F)
  perms[, perm_idx := as.numeric(perms$perm_idx)]
  colnames(perms)[1] <- gsub("#", "", colnames(perms)[1], fixed=T)
  if(!is.null(subset_to_regions)){
    perms <- subset(perms, region_id %in% subset_to_regions)
  }

  # Add normalized columns
  perms[, gnomAD_constrained_prop := n_gnomAD_constrained_genes / n_genes]
  if("n_ubiquitously_expressed_genes" %in% colnames(perms)){
    perms[, prop_ubiquitously_expressed := n_ubiquitously_expressed_genes / n_genes]
  }

  # Normalize DNM excess per permutation
  perms <- normalize.dnms(perms, is.data.table=TRUE,
                          dnm.cohorts=c("DDD", "ASC", "ASC_unaffected",
                                        "DDD_noSig", "ASC_noSig",
                                        "ASC_unaffected_noSig"))

  # Add columns for combined DNM excess between DDD and ASC
  dnm.base.cols <- sapply(colnames(perms)[grep("^DDD_", colnames(perms))],
                          gsub, pattern="^DDD_", replacement="")
  for(col.base in dnm.base.cols){
    new.col <- paste("DDD_plus_ASC", col.base, sep="_")
    ddd.col <- paste("DDD", col.base, sep="_")
    asc.col <- paste("ASC", col.base, sep="_")
    perms[, eval(new.col) := get(ddd.col) + get(asc.col)]
  }

  # Drop unnecessary columns
  extra.dnm.col.idxs <- setdiff(grep("_dnm_", colnames(perms)),
                                grep("_norm_excess_per_gene", colnames(perms)))
  cols.to.drop <- c(colnames(perms)[extra.dnm.col.idxs])
  perms[, (cols.to.drop) := NULL]

  return(perms)
}


#' Summarize segment permutation results
#'
#' Summarize permutation results across all permutations for a single feature
#'
#' @import data.table
#'
#' @param perms.orig data.table of permutation results
#' @param feature name of feature to evaluate
#' @param measure statistic to evaluate \[default: mean\]
#' @param subset_to_regions vector of region IDs to include \[default: include all regions\]
#'
#' @details Valid options for `measure` are mean, median, sum, and frac.any
#'
#' @return matrix of statistic per permutation split by all, deletions, and duplications
#'
#' @export perm.summary
#' @export
perm.summary <- function(perms.orig, feature, measure="mean", subset_to_regions=NULL){
  if(!is.null(subset_to_regions)){
    # Make copy of permutation results to avoid overwriting by data.table
    perms <- perms.orig[region_id %in% subset_to_regions]
  }else{
    perms <- copy(perms.orig)
  }
  if(measure == "mean"){
    ALL <- perms[, mean(get(feature), na.rm=T), by=perm_idx]
    DEL <- perms[cnv == "DEL", mean(get(feature), na.rm=T), by=perm_idx]
    DUP <- perms[cnv == "DUP", mean(get(feature), na.rm=T), by=perm_idx]
  }else if(measure == "median"){
    ALL <- perms[, median(get(feature), na.rm=T), by=perm_idx]
    DEL <- perms[cnv == "DEL", median(get(feature), na.rm=T), by=perm_idx]
    DUP <- perms[cnv == "DUP", median(get(feature), na.rm=T), by=perm_idx]
  }else if(measure == "sum"){
    ALL <- perms[, sum(get(feature), na.rm=T), by=perm_idx]
    DEL <- perms[cnv == "DEL", sum(get(feature), na.rm=T), by=perm_idx]
    DUP <- perms[cnv == "DUP", sum(get(feature), na.rm=T), by=perm_idx]
  }else if(measure == "frac.any"){
    frac.any <- function(x){100 * length(which(x[!is.na(x)]>0)) / length(x[which(!is.na(x))])}
    ALL <- perms[, frac.any(get(feature)), by=perm_idx]
    DEL <- perms[cnv == "DEL", frac.any(get(feature)), by=perm_idx]
    DUP <- perms[cnv == "DUP", frac.any(get(feature)), by=perm_idx]
  }
  if(nrow(DEL) == 0){
    DEL <- copy(ALL)
    DEL[, V1 := 0]
  }
  if(nrow(DUP) == 0){
    DUP <- copy(ALL)
    DUP[, V1 := 0]
  }
  setnames(ALL, "V1", "ALL")
  setnames(DEL, "V1", "DEL")
  setnames(DUP, "V1", "DUP")
  as.matrix(as.data.frame(merge(merge(ALL, DEL), DUP))[, -1])
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

