#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Bifurcate HPOs by effect size of whole-gene deletions of constrained genes


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=T)


# Get command-line arguments & options
option_list <- list(
  make_option(c("-p", "--out-prefix"), type="character", default="rCNV2_pheno_split",
              help="prefix for output files [default %default]", metavar="string"),
  make_option(c("-t", "--threshold"), type="numeric", default=2,
              help="odds ratio threshold for splitting HPOs [default %default]",
              metavar="numeric")
)
args <- parse_args(OptionParser(usage="%prog data hpo_dict.tsv",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options


# Writes args & opts to vars
dat.in <- args$args[1]
hpo_dict.in <- args$args[2]
out.prefix <- opts$`out-prefix`
threshold <- opts$`threshold`


# Read data
dat <- read.table(dat.in, header=T, sep="\t", comment.char="")
hpos <- read.table(hpo_dict.in, header=F, sep="\t", comment.char="")
colnames(hpos) <- c("hpo", "description")
x <- merge(dat, hpos, by="hpo", sort=F, all.x=T, all.y=F)


# Split HPOs by OR
high.hpos <- x$hpo[which(!is.na(x$lnOR) & x$lnOR>=log(threshold))]
low.hpos <- x$hpo[which(!is.na(x$lnOR) & x$lnOR<log(threshold))]


# Write classifications to file
write.table(high.hpos, paste(out.prefix, "developmental.list", sep="."),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(low.hpos, paste(out.prefix, "adult.list", sep="."),
            col.names=F, row.names=F, sep="\t", quote=F)
all.df <- data.frame("HPO"=c(high.hpos, low.hpos),
                     "group"=c(rep("developmental", length(high.hpos)),
                               rep("adult", length(low.hpos))))
colnames(all.df)[1] <- "#HPO"
write.table(all.df, paste(out.prefix, "all_hpos.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
