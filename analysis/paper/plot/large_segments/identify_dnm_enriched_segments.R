#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Identify rCNV segments with significant excesses of DNMs vs. permuted expectation


options(stringsAsFactors=F, scipen=1000)
csqs <- c("PTVs" = "lof", "Mis." = "mis")


######################
### DATA FUNCTIONS ###
######################
# Compile table of significance of DNM enrichment for every segment
tabulate.dnm.enrichment.stats <- function(segs, union.perms, feature,
                                          measure="sum", subset_to_regions=NULL){
  # Get list of segments to evaluate
  if(!is.null(subset_to_regions)){
    region.ids <- subset_to_regions
  }else{
    region.ids <- segs$region_id
  }

  # Gather data from each segment separately
  res <- lapply(region.ids, function(rid){
    obs <- calc.segs.dat(segs, feature=feature, measure=measure,
                         subset_to_regions=rid)[1]
    if(is.na(obs)){
      list(c(rid, NA, NA, NA, NA),
           rep(NA, times=nrow(union.perms)))
    }else{
      perm.vals <- perm.summary(union.perms, feature=feature, measure=measure,
                                subset_to_regions=rid)[, 1]
      n.greater <- length(which(perm.vals >= obs))
      pval <- (n.greater + 1) / (length(perm.vals) + 1)
      list(c(rid, obs, mean(perm.vals), pval),
           as.numeric(perm.vals))
    }
  })

  # Make data.frame() of permuted values
  perm.vals <- as.data.frame(do.call("cbind", lapply(res, function(l){l[[2]]})))
  colnames(perm.vals) <- region.ids

  # Make data.frame() of summary stats
  dnm.stats <- as.data.frame(do.call("rbind", lapply(res, function(l){l[[1]]})))
  colnames(dnm.stats) <- c("region.id", "obs", "exp", "delta", "pval")
  dnm.stats[, -1] <- apply(dnm.stats[, -1], 2, as.numeric)

  # Return both data.frames
  list("stats" = dnm.stats, "perm.vals" = perm.vals)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(data.table, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv perm_res.bed lit_GD_perm_res.bed outdir prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Six positional arguments required: loci.bed, segs.tsv, perm_res.bed, lit_GD_perm_res.bed, outdir, and prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
perm.res.in <- args$args[3]
lit.perm.res.in <- args$args[4]
outdir <- args$args[5]
prefix <- args$args[6]

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# perm.res.in <- "~/scratch/rCNV2_analysis_d2.1000_permuted_segments_bygene.tsv.gz"
# lit.perm.res.in <- "~/scratch/rCNV2_analysis_d2.lit_GDs.1000_permuted_segments_bygene.tsv.gz"
# outdir <- "~/scratch"
# prefix <- "test_perm_bygene"

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
# Also must restrict to NDD-associated segs to match DNM data
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]
segs <- segs[which(segs$region_id %in% get.ndd.region_ids(loci, segs)), ]
analysis.ids <- segs$region_id

# Merge loci & segment data for genome-wide or FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs)

# Load permutation results
perms <- load.perms.bygene(perm.res.in, subset_to_regions=analysis.ids)
lit.perms <- load.perms.bygene(lit.perm.res.in, subset_to_regions=analysis.ids)

# Merge perm results
shared.columns <- intersect(colnames(perms), colnames(lit.perms))
union.perms <- rbind(perms[, ..shared.columns], lit.perms[, ..shared.columns])
rm(perms); rm(lit.perms)

# Collect DNM excess data
lof.res <- tabulate.dnm.enrichment.stats(segs, union.perms, "DDD_plus_ASC_dnm_lof_norm_excess_per_gene", "mean")
mis.res <- tabulate.dnm.enrichment.stats(segs, union.perms, "DDD_plus_ASC_dnm_mis_norm_excess_per_gene", "mean")
dnm.stats <- merge(lof.res$stats, mis.res$stats, by="region.id",
                   sort=F, suffixes=c(".lof", ".mis"))
dnm.perm.vals <- list("lof" = lof.res$perm.vals, "mis" = mis.res$perm.vals)

### TODO: FINISH THIS
