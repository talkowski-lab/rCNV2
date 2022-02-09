#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Reformat noncoding track association stats for supplementary table in rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Reformat track stats table
reformat.table <- function(stats){
  # Sort rows & reorder columns
  row.order <- with(stats, order(family, source, trackname, cnv))
  col.order <- c("family", "source", "trackname", "cnv", "n_elements", "min_size",
                 "median_size", "mean_size", "max_size", "total_bp",
                 "meta.neglog10_p", "meta.lnOR")
  stats <- stats[row.order, col.order]
  stats$meta.neglog10_p <- 10^-stats$meta.neglog10_p
  stats$meta.lnOR <- exp(stats$meta.lnOR)

  # Rename columns
  colnames(stats)[which(colnames(stats) == "family")] <- "annotation_family"
  colnames(stats)[which(colnames(stats) == "trackname")] <- "annotation"
  colnames(stats)[which(colnames(stats) == "total_bp")] <- "total_basepairs"
  colnames(stats)[which(colnames(stats) == "meta.neglog10_p")] <- "cnv_p_value"
  colnames(stats)[which(colnames(stats) == "meta.lnOR")] <- "cnv_odds_ratio"

  return(stats)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--saddlepoint-adj"), action="store_true", default=FALSE,
              help="Update Z-scores and P-values with saddlepoint re-approximation of the null [default: %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog stats.tsv out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: stats.tsv and out.prefix\n"))
}

# Writes args & opts to vars
stats.in <- args$args[1]
out.prefix <- args$args[2]
spa <- opts$`saddlepoint-adj`

# # DEV PARAMETERS
# stats.in <- "~/scratch/rCNV.burden_stats.tsv.gz"
# out.prefix <- "~/scratch/track_stats_test"
# spa <- TRUE

# Load track stats
stats <- load.track.stats(stats.in, spa)

# Reformat stats table
stable <- reformat.table(stats)

# Write table to outfile
outfile <- paste(out.prefix, "annotation_burden_stats.tsv", sep=".")
write.table(stable, outfile, sep="\t", col.names=T, row.names=F, quote=F)
system(paste("gzip -f", outfile), wait=T)
