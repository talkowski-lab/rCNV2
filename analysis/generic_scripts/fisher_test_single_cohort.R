#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform Fisher's exact tests for rare CNVs from a single cohort


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=T)


# List of command-line options
option_list <- list(
  make_option(c("--pheno-table"), type="character", default=NULL,
              help="table with counts of samples per HPO term per cohort [default %default]",
              metavar="file"),
  make_option(c("--cohort-name"), type="character", default=NULL,
              help="name of cohort [default %default]",
              metavar="string"),
  make_option(c("--case-hpo"), type="character", default=NULL,
              help="HPO term to use for case samples [default %default]",
              metavar="string"),
  make_option(c("--control-hpo"), type="character", default='HEALTHY_CONTROL',
              help="HPO term to use for control samples [default %default]",
              metavar="string"),
  make_option(c("--case-column"), type="character", default='case_cnvs',
              help="name of column to use for case CNV counts [default %default]",
              metavar="string"),
  make_option(c("--control-column"), type="character", default='control_cnvs',
              help="name of column to use for control CNV counts [default %default]",
              metavar="string"),
  make_option(c("--keep-n-columns"), type="integer", default=3,
              help="retain first N columns from BED format [default %default]",
              metavar="integer"),
  make_option(c("--precision"), type="integer", default=6,
              help="level of precision for floats [default %default]",
              metavar="integer"),
  make_option(c("-z", "--bgzip"), action="store_true", default=FALSE,
              help="bgzip output BED file [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog bins outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

print(args)

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}
if(is.null(opts$`pheno-table`)){
  stop("Must provide --pheno-table\n")
}
if(is.null(opts$`cohort-name`)){
  stop("Must specify --cohort-name\n")
}
if(is.null(opts$`case-hpo`)){
  stop("Must specify --case-hpo\n")
}

# Writes args & opts to vars
bed.in <- args$args[1]
outfile <- args$args[2]
pheno.table.in <- opts$`pheno-table`
cohort.name <- opts$`cohort-name`
case.hpo <- opts$`case-hpo`
control.hpo <- opts$`control-hpo`
case.col.name <- opts$`case-column`
control.col.name <- opts$`control-column`
keep.n.cols <- opts$`keep-n-columns`
precision <- opts$precision

# Extract sample counts
sample.counts <- get.sample.counts(pheno.table.in, cohort.name, case.hpo, control.hpo)

# Process input BED
bed <- load.cc.cnv.counts(bed.in, case.col.name, sample.counts$case.n,
                          control.col.name, sample.counts$control.n,
                          keep.n.cols)

# Run burden tests
fisher.bed <- fisher.burden.test(bed, keep.n.cols, precision)

# Format output
colnames(fisher.bed)[1] <- "#chr"
if(length(grep(".gz", outfile, fixed=T)) > 0){
  outfile <- tools::file_path_sans_ext(outfile)
}
write.table(fisher.bed, outfile, sep="\t", col.names=T, row.names=F, quote=F)
if(opts$bgzip == TRUE){
  system(paste("bgzip -f", outfile),wait=T)
}
