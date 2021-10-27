#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
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
                 "meta.neg.log10_p", "meta.lnOR")
  stats <- stats[row.order, col.order]
  stats$meta.neg.log10_p <- 10^-stats$meta.neg.log10_p
  stats$meta.lnOR <- exp(stats$meta.lnOR)
  
  # Rename columns
  colnames(stats)[which(colnames(stats) == "family")] <- "annotation_family"
  colnames(stats)[which(colnames(stats) == "trackname")] <- "annotation"
  colnames(stats)[which(colnames(stats) == "total_bp")] <- "total_basepairs"
  colnames(stats)[which(colnames(stats) == "meta.neg.log10_p")] <- "cnv_p_value"
  colnames(stats)[which(colnames(stats) == "meta.lnOR")] <- "cnv_odds_ratio"
  
  return(stats)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
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
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# stats.in <- "~/scratch/rCNV.burden_stats.tsv.gz"
# out.prefix <- "~/scratch/track_stats_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/noncoding_association/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(gsub("scripts/", "plot/", script.dir, fixed=T), "common_functions.R", sep="/"))

# Load track stats
stats <- load.track.stats(stats.in)

# Reformat stats table
stable <- reformat.table(stats)

# Write table to outfile
outfile <- paste(out.prefix, "annotation_burden_stats.tsv", sep=".")
write.table(stable, outfile, sep="\t", col.names=T, row.names=F, quote=F)
system(paste("gzip -f", outfile), wait=T)
