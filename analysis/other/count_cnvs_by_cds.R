#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Summarize gene counts by CDS overlap
# (Subroutine of optimize_cds_overlap.sh)


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=T)


# Load annotated CNVs from BED and extract phenotypes + max CDS
load.annotated.cnvs <- function(cnvs.in){
  cnvs <- read.table(cnvs.in, header=T, sep="\t", comment.char="")[, c("phenos", "cds_per_gene")]
  cnvs$max_cds <- sapply(cnvs$cds_per_gene, function(vals){
    if(is.na(vals)){
      return(0)
    }else if(vals == ""){
      return(0)
    }else{
      return(max(as.numeric(unlist(strsplit(as.character(vals), split=";")), na.rm=T)))
    }
  })
  list("case" = cnvs$max_cds[which(cnvs$phenos != "HEALTHY_CONTROL")],
       "control" = cnvs$max_cds[which(cnvs$phenos == "HEALTHY_CONTROL")])
}


# Count CNV carriers and non-carriers at a minimum CDS threshold
count.cnvs.at.cds <- function(cds.steps, cds.vals, n.samp){
  res <- as.data.frame(t(sapply(cds.steps, function(min.cds){
    n.hits <- length(which(cds.vals >= min.cds))
    n.ref <- n.samp - n.hits
    c(min.cds, n.hits, n.ref)
  })))
  colnames(res) <- c("min.cds", "cnv", "ref")
  return(res)
}


# Merge & summarize CNV counts per CDS step
summarize.counts <- function(case.df, control.df){
  res <- merge(control.df, case.df, by="min.cds", sort=F, suffixes=c(".control", ".case"))
  res$"case_pct" <- res$cnv.case / (res$cnv.case + res$ref.case)
  res$"control_pct" <- res$cnv.control / (res$cnv.control + res$ref.control)
  for(pheno in c("control", "case")){
    colnames(res)[which(colnames(res) == paste("cnv", pheno, sep="."))] <- paste(pheno, "alt", sep="_")
    colnames(res)[which(colnames(res) == paste("ref", pheno, sep="."))] <- paste(pheno, "ref", sep="_")
  }
  fstats <- as.data.frame(do.call("rbind", (lapply(1:nrow(res), function(i){
    fisher.burden.test.single(as.numeric(as.vector(res[i, c("case_alt", "control_alt", "case_ref", "control_ref")])))
  })))[, c("p", "OR")])
  colnames(fstats) <- c("fisher_phred_p", "odds_ratio")
  fstats$fisher_phred_p <- -log10(as.numeric(fstats$fisher_phred_p))
  as.data.frame(cbind(res, fstats))
}


# List of command-line options
option_list <- list(
  make_option(c("--n-cases"), type="integer", default=NULL,
              help="number of case samples in cohort [required]",
              metavar="int"),
  make_option(c("--n-controls"), type="integer", default=NULL,
              help="number of control samples in cohort [required]",
              metavar="int"),
  make_option(c("-s", "--step-size"), type="numeric", default=0.01,
              help="step size for CDS increments [default: %default]",
              metavar="float")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog cnvs outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to vars
cnvs.in <- args$args[1]
outfile <- args$args[2]
n.cases <- opts$`n-cases`
n.controls <- opts$`n-controls`
step.size <- opts$`step-size`

# Load and process CNVs
cnvs <- load.annotated.cnvs(cnvs.in)

# Count case & control CNVs at each CDS step
cds.steps <- seq(step.size, 1, by=step.size)
case.df <- count.cnvs.at.cds(cds.steps, cnvs$case, n.cases)
control.df <- count.cnvs.at.cds(cds.steps, cnvs$control, n.controls)

# Merge & summarize properties
res <- summarize.counts(case.df, control.df)

# Reformat as dummy BED file for meta-analysis and write to outfile
coords <- data.frame("chr" = 100 * res$min.cds,
                     "start" = 100 * res$min.cds,
                     "end" = 100 * res$min.cds,
                     "min_cds" = paste("mincds", res$min.cds, sep="_"))
res <- cbind(coords, res[, -1])
colnames(res)[1] <- paste("#", colnames(res)[1], sep="")
write.table(res, outfile, col.names=T, row.names=F, sep="\t", quote=F)

