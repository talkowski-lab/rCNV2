#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform meta-analysis of CNV associations


# Load libraries
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(metafor, quietly=T)
require(EQL, quietly=T)

# Helper function to drop effect size confidence interval columns from per-cohort stats
drop.ci <- function(df){
  keep.idxs <- intersect(grep("OR_lower", colnames(df), fixed=T, invert=T),
                         grep("OR_upper", colnames(df), fixed=T, invert=T))
  df[, keep.idxs]
}

# List of command-line options
option_list <- list(
  make_option(c("--or-corplot"), type="character", default=NULL,
              help="output .jpg file for pairwise odds ratio correlation plot [default: %default]",
              metavar="path"),
  make_option(c("--model"), type="character", default="fe",
              help="specify meta-analysis model ('re': random effects, 'fe': fixed effects, 'mh': Mantel-Haenszel) [default: '%default']",
              metavar="string"),
  make_option(c("--conditional-exclusion"), type="character", default=NULL,
              help="BED file annotated with cohorts to conditionally exclude on a locus-specific basis",
              metavar="path"),
  make_option(c("--p-is-neg-log10"), action="store_true", default=FALSE,
              help="provided P-values are Phred-scaled (-log10(P)) [default: %default]"),
  make_option(c("--spa"), action="store_true", default=FALSE,
              help="apply saddlepoint approximation of null distribution [default: %default]"),
  make_option(c("--spa-exclude"), type="character", metavar="BED3+",
              help="BED file of regions to exclude when applying --spa option [default: include all regions"),
  make_option(c("--winsorize"), default=1, type="numeric", metavar="float",
              help="maximum quantile of meta-analysis statistics to include when applying --spa option [default: %default]"),
  make_option(c("--mirror-saddle"), action="store_true", default=FALSE,
              help="mirror bottom 50% of Z-scores when approximating null distribution via saddlepoint [default: %default]"),
  make_option(c("--adjust-biobanks"), action="store_true", default=FALSE,
              help="include biobank label as a covariate in meta-analysis [default: %default]"),
  make_option(c("--adjust-inflation"), action="store_true", default=FALSE,
              help="estimate baseline effect size inflation factor per cohort and include as a covariate in meta-analysis [default: %default]"),
  make_option(c("--min-cases"), default=1, type="numeric", metavar="integer",
              help="minimum number of cases required for a cohort to be included in meta-analysis [default: %default]"),
  make_option(c("--probe-counts"), type="character", metavar="path", default=NULL,
              help="BED file with one column of control probe counts per cohort (will be used as covariate) [default: %default]"),
  make_option(c("--no-fdr"), action="store_true", default=FALSE,
              help="do not compute Benjamini-Hochberg FDR q-values [default: compute FDR]"),
  make_option(c("--no-secondary"), action="store_true", default=FALSE,
              help="do not compute secondary P-value [default: %default]"),
  make_option(c("--keep-n-columns"), type="integer", default=3,
              help="retain first N columns from BED format [default: %default]",
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
p.is.neg.log10 <- opts$`p-is-neg-log10`
spa <- opts$spa
spa.xbed.in <- opts$`spa-exclude`
winsorize <- opts$winsorize
mirror.saddle <- opts$`mirror-saddle`
adjust.biobanks <- opts$`adjust-biobanks`
adjust.inflation <- opts$`adjust-inflation`
min.cases <- opts$`min-cases`
probe.counts.in <- opts$`probe-counts`
calc.fdr <- !(opts$`no-fdr`)
secondary <- !(opts$`no-secondary`)
keep.n.cols <- opts$`keep-n-columns`

# # Dev parameters
# setwd("~/scratch")
# infile <- "HP0012759.rCNV.DEL.gene_burden.meta_analysis.input.txt"
# outfile <- "HP0012759.rCNV.DEL.gene_burden.meta_analysis.bed"
# corplot.out <- "corplot.test.jpg"
# model <- "fe"
# cond.excl.in <- "gencode.v19.canonical.pext_filtered.cohort_exclusion.bed.gz"
# p.is.neg.log10 <- T
# spa <- T
# spa.xbed.in <- "DEL_GDs.bed.gz"
# winsorize <- 0.99
# mirror.saddle <- F
# adjust.biobanks <- F
# adjust.inflation <- F
# min.cases <- 300
# probe.counts.in <- NULL
# calc.fdr <- T
# secondary <- T
# keep.n.cols <- 4

# Read list of cohorts to meta-analyze
cohort.info <- read.table(infile, header=F, sep="\t")
ncohorts <- nrow(cohort.info)
stats.list <- lapply(1:ncohorts, function(i){
  read.assoc.stats.single(cohort.info[i, 2], cohort.info[i, 1], p.is.neg.log10,
                          keep.n.cols, keep.or.confint=TRUE)
})
names(stats.list) <- cohort.info[, 1]

# Calculate cohort-specific inflation factors, if optioned
if(adjust.inflation){
  cohort.inflation <- sapply(stats.list, estimate.cohort.inflation)
  names(cohort.inflation) <- names(stats.list)
}else{
  cohort.inflation <- NULL
}
stats.list <- lapply(stats.list, drop.ci)

# Read probe counts, if optioned
if(!is.null(probe.counts.in)){
  probe.counts <- read.table(probe.counts.in, header=T, sep="\t", comment.char="")
  colnames(probe.counts)[1] <- c("chr")
  # Sanity check to make sure exact same entries are represented in stats.list
  if(!identical(stats.list[[1]][1:keep.n.cols], probe.counts[1:keep.n.cols])){
    stop(paste("Issue with input files: First", keep.n.cols, "columns of input",
               "association stats files and --probe-counts file must be identical."))
  }
}else{
  probe.counts <- NULL
}

# Plot correlations of odds ratios between cohorts, if optioned
if(!is.null(corplot.out)){
  jpeg(corplot.out, res=300,
       height=300*(3.5+(ncohorts/1.5)),
       width=300*(4+(ncohorts/1.5)))
  or.corplot.grid(stats.list, pt.cex=0.4)
  dev.off()
}

# Merge data across cohorts
stats.merged <- combine.single.cohort.assoc.stats(stats.list, cond.excl.in,
                                                  min.cases, keep.n.cols)
# Load saddlepoint exclusion regions, if optioned
if(is.null(spa.xbed.in)){
  spa.xbed <- NULL
}else{
  spa.xbed <- load.bed3(spa.xbed.in)
}

# Conduct meta-analysis & write to file
stats.meta <- meta(stats.merged, cohort.info[, 1], model=model, saddle=spa,
                   adjust.biobanks=adjust.biobanks, cohort.inflation=cohort.inflation,
                   probe.counts=probe.counts, saddle.exclusion=spa.xbed,
                   winsorize=winsorize, mirror.saddle=mirror.saddle,
                   calc.fdr=calc.fdr, secondary=secondary, keep.n.cols=keep.n.cols)
colnames(stats.meta)[1] <- "#chr"
write.table(stats.meta, outfile, sep="\t",
            row.names=F, col.names=T, quote=F)
