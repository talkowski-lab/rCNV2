#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute & plot case:control effect sizes before and after loose noncoding filtering for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################


##########################
### PLOTTING FUNCTIONS ###
##########################


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
  make_option(c("--cnv"), help="CNV type.", default="DEL")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog sumstats.bed whitelist.genes out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: sumstats.bed, whitelist.genes,",
             "and out.prefix\n"))
}

# Writes args & opts to vars
sumstats.in <- args$args[1]
whitegenes.in <- args$args[2]
out.prefix <- args$args[3]
rcnv.config <- opts$`rcnv-config`
cnv <- opts$`cnv`

# # DEV PARAMETERS
# sumstats.in <- "~/scratch/HP0000118.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
# whitegenes.in <- "~/scratch/loose_noncoding_whitelist.genes.list"
# out.prefix <- "~/scratch/unconstrained_lnORs_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/noncoding_association/"
# cnv <- "DEL"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load likely unconstrained genes
whitegenes <- unique(as.character(read.table(whitegenes.in, header=F)[, 1]))

# Load sumstats & annotate with unconstrained status
ss <- load.sumstats(sumstats.in)
ss$unconstrained <- ss$gene %in% whitegenes
