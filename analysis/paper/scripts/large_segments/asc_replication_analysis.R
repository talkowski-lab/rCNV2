#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Assess replication of rCNV segments from ASC/SPARK gCNV callset


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Run replication analysis for a single set of values
replicate.segment <- function(vals, n_control, n_case){
  if(any(sapply(vals, is.na))){
    c(rep(NA, 4), rep(0, 4))
  }else{
    control_alt <- as.numeric(vals[1])
    control_ref <- n_control - control_alt
    case_alt <- as.numeric(vals[2])
    case_ref <- n_case - case_alt
    fisher.mat <- matrix(c(control_ref, case_ref, control_alt, case_alt),
                         byrow=T, nrow=2)
    fisher.res <- fisher.test(fisher.mat, alternative="greater")
    fstats <- as.numeric(c(fisher.res$p.value, fisher.res$estimate,
                           fisher.res$conf.int))
    c(case_alt, control_alt, case_alt / n_case, control_alt / n_control, fstats)
  }
}

# Compute replication P-values and odds ratios for all segments
calc.replication.stats <- function(segs, n_case=13786-2155, n_control=5098){
  rep.res <- t(apply(segs[, c("asc_spark_replication_controls",
                            "asc_spark_replication_cases")],
                   1, replicate.segment, n_control=n_control, n_case=n_case))
  rep.res <- as.data.frame(rep.res)
  rep.res <- cbind(segs$region_id, rep.res)
  colnames(rep.res) <- c("region_id", "case_alt", "control_alt", "case_rate",
                         "control_rate", "p_value", "OR", "OR_lower", "OR_upper")
  return(rep.res)
}

# Summarize replication outcome for a subset of segments
replication.summary <- function(rep.res, seg.ids, subset.title){
  hits <- which(rep.res$region_id %in% seg.ids)
  n.all <- length(hits)
  testable <- intersect(hits, which(rep.res$case_alt + rep.res$control_alt > 0))
  n.testable <- length(testable)
  suggestive <- which(rep.res[testable, "OR"] >= 1)
  n.suggestive <- length(suggestive)
  nomsig <- which(rep.res[testable, "p_value"] <= 0.05)
  n.nomsig <- length(nomsig)

  # Print results to screen
  cat(paste("\n\n", subset.title, " replication vs. ASC/SPARK:\n", sep=""))
  cat(paste("All segments:", n.all, "\n"))
  cat(paste("Segments able to be evaluated: ", n.testable, "/", n.all, " (",
            round(100*n.testable/n.all, 2), "%)\n", sep=""))
  cat(paste("Segments with consistent direction of effects: ", n.suggestive,
            "/", n.testable, " (", round(100*n.suggestive/n.testable, 2), "%)\n", sep=""))
  cat(paste("Segments at nominal significance: ", n.nomsig, "/", n.testable,
            " (", round(100*n.nomsig/n.testable, 2), "%)\n", sep=""))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: loci.bed and segs.tsv\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]

# Collect replication P-values and odds ratios for all analysis segments
rep.res <- calc.replication.stats(segs)

# Make additional subset strictly of NDD and ASD segments (to match ASC phenotypes)
dev.seg.ids <- get.developmental.region_ids(loci, segs)
neuro.seg.ids <- get.neuro.region_ids(loci, segs)
ndd.seg.ids <- intersect(dev.seg.ids, neuro.seg.ids)
asd.seg.ids <- loci$region_id[which(sapply(loci$hpos, function(h){"HP:0100753" %in% h}))]

# Make additional subsets split by discovery significance
gw.seg.ids <- segs$region_id[which(segs$gw_sig)]
fdr.seg.ids <- segs$region_id[which(segs$fdr_sig)]

# Summarize replication stats for several subsets of segments
replication.summary(rep.res, segs.all$region_id, "All discovery segments")
replication.summary(rep.res, gw.seg.ids, "Genome-wide significant segments")
replication.summary(rep.res, fdr.seg.ids, "FDR-significant segments")
replication.summary(rep.res, ndd.seg.ids, "NDD-associated segments")
replication.summary(rep.res, intersect(ndd.seg.ids, gw.seg.ids),
                    "NDD-associated segments at G-W sig.")
replication.summary(rep.res, intersect(ndd.seg.ids, fdr.seg.ids),
                    "NDD-associated segments at G-W sig.")
replication.summary(rep.res, asd.seg.ids, "ASD-associated segments")
