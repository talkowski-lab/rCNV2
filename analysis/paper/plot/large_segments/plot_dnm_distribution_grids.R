#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Quantitative assessment of DNM distributions for large segments for rCNV paper

# DEV NOTE: THIS SCRIPT AND ANALYSIS ARE CURRENTLY UNDER DEVELOPMENT 

options(stringsAsFactors=F, scipen=1000)
csqs <- c("lof", "mis", "syn")
cohorts <- c("ddd", "asc", "asc_control")


######################
### DATA FUNCTIONS ###
######################
# Load table of mutation rates
load.mutrates <- function(mutrates.in){
  mutrates <- read.table(mutrates.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(mutrates)[1] <- "gene"
  return(mutrates)
}

# Load count of de novo mutations for a single study and residualize vs. expected
load.dnms <- function(dnms.in, mutrates){
  dnms <- read.table(dnms.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(dnms)[1] <- "gene"
  dnms <- merge(dnms, mutrates, all.x=T, all.y=F, sort=F, by="gene")
  for(csq in csqs){
    counts <- dnms[, which(colnames(dnms) == csq)]
    mus <- dnms[, which(colnames(dnms) == paste("mu", csq, sep="_"))]
    n.dnms <- sum(counts, na.rm=T)
    expected <- (mus / sum(mus, na.rm=T)) * n.dnms
    dnms[paste("expected", csq, sep="_")] <- expected
    residuals <- counts - expected
    dnms[paste("excess", csq, sep="_")] <- residuals
    dnms[paste("fold_enrich", csq, sep="_")] <- counts/expected
  }
  return(dnms)
}

# Compute grid of DNM excess segment count for a single cohort, consequence, and CNV type
calc.excess.dnm.grid.single <- function(segs, dnms, csq, cnv.type="CNV", fold.cutoffs=2^(1:4)){
  # Subset to regions with a specific CNV type, if specified
  if(cnv.type %in% c("DEL", "DUP")){
    segs <- segs[which(segs$cnv==cnv.type), ]
  }
  
  # Create matrix of segments with Y genes that are at least X-fold enriched for DNMs
  res <- do.call("cbind", lapply(fold.cutoffs, function(f){
    counts <- sapply(segs$genes, function(genes){
      folds <- dnms[which(dnms$gene %in% genes), paste("fold_enrich", csq, sep="_")]
      folds <- folds[which(!is.na(folds) & !is.infinite(folds))]
      length(which(folds >= f))
    })
    c(length(which(counts==0)), length(which(counts==1)), 
      length(which(counts==2)), length(which(counts>2)))
  }))
  colnames(res) <- paste(fold.cutoffs, "fold", sep=".")
  rownames(res) <- c(paste(0:2, "gene", sep="."), "more")
  
  return(res)
}

# Compute grid of expected (null) DNM excess segment count for a single cohort, consequence, and CNV type
calc.null.dnm.grid.single <- function(segs, dnms, csq, cnv.type="CNV", fold.cutoffs=2^(1:4)){
  # Subset to regions with a specific CNV type, if specified
  if(cnv.type %in% c("DEL", "DUP")){
    segs <- segs[which(segs$cnv==cnv.type), ]
  }
  
  # Create matrix of expected number of segments with Y genes that are at least X-fold enriched for DNMs
  res <- do.call("cbind", lapply(fold.cutoffs, function(f){
    # Compute baseline odds of a gene drawn at random being at least f-fold enriched for DNMs
    universe <- (dnms[, paste("fold_enrich", csq, sep="_")] >= f)
    m <- length(which(universe==TRUE))
    n <- length(which(universe==FALSE))
    baseline <- length(which(universe==TRUE)) / length(which(!is.na(universe)))
    
    # For each segment, compute the hypergeometric probability that the segment 
    # has exactly 0, 1, 2, or >2 genes at least f-fold enriched
    probs <- t(sapply(segs$n_genes, function(k){
      c(dhyper(0:2, m, n, k), phyper(2, m, n, k, lower.tail=FALSE))
    }))
    
    # Compute expected number of segments as the sum of probabilities
    apply(probs, 2, sum)
  }))
  colnames(res) <- paste(fold.cutoffs, "fold", sep=".")
  rownames(res) <- c(paste(0:2, "gene", sep="."), "more")
  
  return(res)
}


# Compute excess DNM grids for all consequences and cnv types for a single cohort
calc.excess.dnm.grids <- function(segs, dnms, fold.cutoffs=2^(1:4), null=FALSE){
  res <- lapply(c("CNV", "DEL", "DUP"), function(cnv.type){
    sublist <- lapply(csqs, function(csq){
      if(null==FALSE){
        calc.excess.dnm.grid.single(segs, dnms, csq, cnv.type, fold.cutoffs)
      }else{
        calc.null.dnm.grid.single(segs, dnms, csq, cnv.type, fold.cutoffs)
      }
    })
    names(sublist) <- csqs
    return(sublist)
  })
  names(res) <- c("CNV", "DEL", "DUP")
  return(res)
}

# Helper function to compute excess DNM grids for all cohorts
make.all.dnm.grids <- function(segs, dnms, null=FALSE){
  grids <- lapply(cohorts, function(cohort){
    calc.excess.dnm.grids(neuro.segs, dnms[[cohort]], null=null)
  })
  names(grids) <- cohorts
  return(grids)
}


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
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv ddd.tsv asc.tsv asc.control.tsv mutrates.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 7){
  stop(paste("Seven positional arguments required: loci.bed, segs.tsv, ddd.tsv, asc.tsv, asc.control.tsv, mutrates.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
ddd.dnms.in <- args$args[3]
asc.dnms.in <- args$args[4]
asc.control.dnms.in <- args$args[5]
gene.mutrates.in <- args$args[6]
out.prefix <- args$args[7]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# ddd.dnms.in <- "~/scratch/ddd_dnm_counts.tsv.gz"
# asc.dnms.in <- "~/scratch/asc_dnm_counts.tsv.gz"
# asc.control.dnms.in <- "~/scratch/asc_dnm_counts.unaffecteds.tsv.gz"
# gene.mutrates.in <- "~/scratch/gene_mutation_rates.tsv.gz"
# out.prefix <- "~/scratch/test_effect_sizes"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load dnm counts and residualize by mutation rates
mutrates <- load.mutrates(gene.mutrates.in)
ddd <- load.dnms(ddd.dnms.in, mutrates)
asc <- load.dnms(asc.dnms.in, mutrates)
asc.control <- load.dnms(asc.control.dnms.in, mutrates)
dnms <- list("ddd"=ddd, "asc"=asc, "asc_control"=asc.control)

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Restrict to segments nominally significant in at least one phenotype
segs <- segs[which(segs$nom_sig), ]

# Get list of neuro loci
neuro.region_ids <- get.neuro.region_ids(loci, segs)
neuro.segs <- segs[which(segs$region_id %in% neuro.region_ids), ]

# Compute grids of excess DNMs for neuro segs in every cohort
dnm.excess.grids <- make.all.dnm.grids(neuro.segs, dnms)
dnm.null.grids <- make.all.dnm.grids(neuro.segs, dnms, null=TRUE)
