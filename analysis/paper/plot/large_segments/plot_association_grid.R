#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot grid summarizing large segment association across phenotypes


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load summary dataframe for loci
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
  loci$name <- sapply(strsplit(loci$region_id, split="_"), function(parts){parts[4]})
  loci$size <- paste(prettyNum(round(loci$cred_intervals_size/1000, 0), big.mark=","), "kb", sep=" ")
  
  return(loci)
}

# Load flat table of effect sizes per region per phenotype per CNV
load.lnors <- function(lnors.in){
  lnors <- read.table(lnors.in, header=T, sep="\t")
  lnors[, 4:6] <- apply(lnors[, 4:6], 2, as.numeric)
  return(lnors)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot a single semicircle
plot.semicircle <- function(x, y, r, side="top", color="black"){
  # Adopted from https://stackoverflow.com/questions/31538534/plotting-half-circles-in-r
  rs <- seq(0, pi, len=100)
  xm <- r*cos(rs)
  ym <- r*sin(rs)
  if(side == "bottom"){
    sign <- -1
  }else{
    sign <- 1
  }
  xc <- x + (sign * xm)
  yc <- y + (sign * ym)
  polygon(xc, yc, col=color, border=NA)
}

# Plot data for a single locus-phenotype pair
plot.locus.single_hpo <- function(loci, lnors, region_id, hpo, x, y, max.lnor=5, ns.col="gray60"){
  # Scale radii
  radii <- as.numeric(sapply(c("DUP", "DEL"), function(cnv){
    x <- lnors$lnor[which(lnors$region_id==region_id & lnors$hpo==hpo & lnors$cnv==cnv)]
    min(c(max(c(0, x), na.rm=T), max.lnor)) / max.lnor
  }))
  # Determine significance
  colors <- sapply(c("DUP", "DEL"), function(cnv){
    if(hpo %in% unlist(loci$hpos[which(loci$region_id==region_id)])){
      as.character(cnv.colors[which(names(cnv.colors)==cnv)])
    }else{
      ns.col
    }
  })
  # Plot points
  plot.semicircle(x, y, r=radii[1], side="top", color=colors[1])
  plot.semicircle(x, y, r=radii[2], side="bottom", color=colors[2])
}

# Plot a row of all hpos per region
plot.locus.row <- function(loci, lnors, region_id, y, max.lnor=5, ns.col="gray60"){
  
}




#####################
### RSCRIPT BLOCK ###
#####################
require(optparse)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed lnors.tsv", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Seven positional arguments: loci.bed, lnors.tsv, hpos.tsv\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
lnors.in <- args$args[2]
hpos.in <- args$args[3]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# lnors.in <- "~/scratch/rCNV.final_segments.loci.all_effect_sizes.tsv.gz"
# hpos.in <- "~/scratch/hpos.txt"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Read data
loci <- load.loci(loci.in)
lnors <- load.lnors(lnors.in)
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
