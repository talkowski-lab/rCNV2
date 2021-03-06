#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform meta-analysis of CNV associations


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=T)
require(metafor, quietly=T)
require(EQL, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--or-corplot"), type="character", default=NULL,
              help="output .jpg file for pairwise odds ratio correlation plot [default %default]",
              metavar="path"),
  make_option(c("--model"), type="character", default="fe",
              help="specify meta-analysis model ('re': random effects, 'fe': fixed effects, 'mh': Mantel-Haenszel) [default '%default']",
              metavar="string"),
  make_option(c("--conditional-exclusion"), type="character", default=NULL,
              help="BED file annotated with cohorts to conditionally exclude on a locus-specific basis",
              metavar="path"),
  make_option(c("--p-is-phred"), action="store_true", default=FALSE,
              help="provided P-values are Phred-scaled (-log10(P)) [default %default]"),
  make_option(c("--spa"), action="store_true", default=FALSE,
              help="apply saddlepoint approximation of null distribution [default %default]"),
  make_option(c("--no-fdr"), action="store_true", default=FALSE,
              help="do not compute Benjamini-Hochberg FDR q-values [default: compute FDR]"),
  make_option(c("--no-secondary"), action="store_true", default=FALSE,
              help="do not compute secondary P-value [default %default]"),
  make_option(c("--keep-n-columns"), type="integer", default=3,
              help="retain first N columns from BED format [default %default]",
              metavar="integer")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog infile outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
infile <- args$args[1]
outfile <- args$args[2]
corplot.out <- opts$`or-corplot`
model <- opts$model
cond.excl.in <- opts$`conditional-exclusion`
p.is.phred <- opts$`p-is-phred`
spa <- opts$spa
calc.fdr <- !(opts$`no-fdr`)
secondary <- !(opts$`no-secondary`)
keep.n.cols <- opts$`keep-n-columns`

# # Dev parameters
# setwd("~/scratch/assoc_test_dev")
# infile <- "HP0001250.rCNV.loose_noncoding.DEL.crb_burden.meta_analysis.input.txt"
# outfile <- "HP0001250.rCNV.DEL.crb_meta_test_results.bed"
# corplot.out <- "corplot.test.jpg"
# model <- "fe"
# cond.excl.in <- "rCNV.crbs.cohort_exclusion.bed.gz"
# p.is.phred <- T
# spa <- T
# calc.fdr <- T
# secondary <- T
# keep.n.cols <- 4

# Read list of cohorts to meta-analyze
cohort.info <- read.table(infile, header=F, sep="\t")
ncohorts <- nrow(cohort.info)
stats.list <- lapply(1:ncohorts, function(i){
  read.assoc.stats.single(cohort.info[i, 2], cohort.info[i, 1], p.is.phred, keep.n.cols)
})
names(stats.list) <- cohort.info[, 1]

# Plot correlations of odds ratios between cohorts, if optioned
if(!is.null(corplot.out)){
  jpeg(corplot.out, res=300,
       height=300*(3.5+(ncohorts/2)),
       width=300*(4+(ncohorts/2)))
  or.corplot.grid(stats.list, pt.cex=0.25)
  dev.off()
}

# Conduct meta-analysis & write to file
stats.merged <- combine.single.cohort.assoc.stats(stats.list, cond.excl.in, keep.n.cols)
stats.meta <- meta(stats.merged, cohort.info[, 1], model=model,
                   saddle=spa, calc.fdr=calc.fdr, secondary=secondary,
                   keep.n.cols=keep.n.cols)
colnames(stats.meta)[1] <- "#chr"
write.table(stats.meta, outfile, sep="\t",
            row.names=F, col.names=T, quote=F)
