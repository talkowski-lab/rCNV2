#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot large segment permutation results (gene set-based test)


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load & parse permutation results
load.perms <- function(perm.res.in){
  # Read full permutation table
  perms <- read.table(perm.res.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(perms)[1] <- gsub("#", "", colnames(perms)[1], fixed=T)
  perm.range <- sort(unique(perms$perm_idx))
  
  # Split each permutation into a list of dfs (one per perm)
  lapply(perm.range, function(i){
    as.data.frame(perms[which(perms$perm_idx==i), -(which(colnames(perms)=="perm_idx"))])
  })
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(MASS, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv perm_res.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: loci.bed, segs.tsv, perm_res.tsv, output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
perm.res.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# perm.res.in <- "~/scratch/gene_perm_test.annotated.tsv.gz"
# out.prefix <- "~/scratch/test_perm_bygene"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Load permutation results
perms <- load.perms(perm.res.in)


plot.seg.perms(gw, perms, feature="n_HPOmatched_genes", measure="sum", n.bins=50,
               x.title="HPO-Matched Genes",
               diamond.cex=1.25, parmar=c(2.2, 2, 0, 1.6))


