#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Assess correlation of effect sizes across cohorts for a single phenotype


# Load libraries
require(rCNV2, quietly=T)
require(optparse, quietly=T)


# List of command-line options
option_list <- list(
  make_option(c("--conditional-exclusion"), type="character", default=NULL,
              help="BED file annotated with cohorts to conditionally exclude on a locus-specific basis",
              metavar="path"),
  make_option(c("--min-cases"), default=1, type="numeric", metavar="integer",
              help="minimum number of cases required for a cohort to be included in meta-analysis [default: %default]"),
  make_option(c("--keep-n-columns"), type="integer", default=3,
              help="retain first N columns from BED format [default: %default]",
              metavar="integer")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog del.infile dup.infile outfile.png",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
del.infile <- args$args[1]
dup.infile <- args$args[2]
outfile <- args$args[3]
cond.excl.in <- opts$`conditional-exclusion`
min.cases <- opts$`min-cases`
keep.n.cols <- opts$`keep-n-columns`
hpo <- gsub("HP", "HP:", unlist(strsplit(basename(del.infile), split=".", fixed=T))[1], fixed=T)

# # Dev parameters
# setwd("~/scratch")
# del.infile <- "HP0000118.DEL.cohort_effect_size_comparison.input.txt"
# dup.infile <- "HP0000118.DUP.cohort_effect_size_comparison.input.txt"
# outfile <- "test_cohort_or_comparisons.png"
# cond.excl.in <- "GRCh37.200kb_bins_10kb_steps.raw.cohort_exclusion.bed.gz"
# min.cases <- 300
# keep.n.cols <- 3
# hpo <- "HP:0000118"

# Load deletion & duplication stats
cnv.stats <- lapply(c(del.infile, dup.infile), function(infile){
  # Load stats from each cohort
  cohort.info <- read.table(infile, header=F, sep="\t")
  ncohorts <- nrow(cohort.info)
  stats.list <- lapply(1:ncohorts, function(i){
    read.assoc.stats.single(cohort.info[i, 2], cohort.info[i, 1], TRUE, keep.n.cols)
  })
  names(stats.list) <- cohort.info[, 1]
  n.cases <- sapply(stats.list, function(df){
    max(apply(df[, grep(".case_", colnames(df), fixed=T)], 1, sum), na.rm=T)
  })
  keep.cohorts <- which(n.cases >= min.cases)
  stats.list <- stats.list[keep.cohorts]

  # Merge across cohorts
  stats.merged <- combine.single.cohort.assoc.stats(stats.list, cond.excl.in,
                                                    min.cases, keep.n.cols)
  stats.merged <- stats.merged[, c(which(colnames(stats.merged) == "exclude_cohorts"),
                                   grep("fisher_OR", colnames(stats.merged)))]
  stats.merged <- stats.merged[, grep("OR_", colnames(stats.merged), invert=T)]

  # Apply conditional exclusion filter
  or.idxs <- grep("fisher_OR", colnames(stats.merged))
  for(i in or.idxs){
    cohort <- unlist(strsplit(colnames(stats.merged)[i], split=".", fixed=T))[1]
    na.rows <- grep(cohort, stats.merged$exclude_cohorts)
    stats.merged[na.rows, i] <- NA
  }
  stats.merged$exclude_cohorts <- NULL

  return(stats.merged)
})

# Merge deletion & duplication stats into a single data frame for plotting
plot.df <- as.data.frame(log2(do.call("rbind", cnv.stats)))
plot.df <- plot.df[which(!apply(plot.df, 1, function(vals){all(is.na(vals))})), ]

# Plot correlation grid
png(outfile, height=2100, width=2100, res=300)
par(mfrow=rep(ncol(plot.df), 2))
sapply(1:ncol(plot.df), function(x){
  sapply(1:ncol(plot.df), function(y){
    if(x < y){
      dens.scatter(plot.df[, x], plot.df[, y], pt.cex=0.1,
                   parmar=rep(0, 4))
    }else{
      par(bty="n", mar=rep(0.1, 4))
      plot(NA, ylim=0:1, xlim=0:1, xaxt="n", yaxt="n", xlab="", ylab="")
      if(x == y){
        if(x > 1){
          mtext(3, line=-1, cex=0.7, font=2,
                text=cohort.abbrevs[unlist(strsplit(colnames(plot.df)[x], split=".", fixed=T))[1]])
        }
        if(y < ncol(plot.df)){
          mtext(4, line=-1, cex=0.7, font=2,
                text=cohort.abbrevs[unlist(strsplit(colnames(plot.df)[y], split=".", fixed=T))[1]])
        }
      }
      par(bty="o")
    }
  })
})
dev.off()
