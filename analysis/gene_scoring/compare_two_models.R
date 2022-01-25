#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Helper script to compare performance of two different models for gene scoring


# Load arguments & options
options(stringsAsFactors=F, scipen=1000)
require(rCNV2, quietly=T)
args <- commandArgs(trailingOnly=T)

# # DEV parameters
# setwd("~/scratch/")
# args <- c("rCNV.DEL.gene_scores.neuralnet.baseline.tsv",
#           "rCNV.DEL.gene_scores.neuralnet.oneshot.tsv",
#           "gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz",
#           "gnomad_oe_lof_upper",
#           "test_score_comparison.png")

# Read data
s1 <- read.table(args[1], header=T, sep="\t", comment.char="")
s2 <- read.table(args[2], header=T, sep="\t", comment.char="")
x <- merge(s1, s2, by="gene", all=F, sort=F, suffixes=c(".s1", ".s2"))
x$d <- x$score.s2 - x$score.s1
features <- load.features(args[3])
x <- merge(x, features[, c("gene", args[4])], by="gene",
           all.x=T, all.y=F, sort=F)

# Plot difference in scores vs. optimization metric
png(args[5], height=900, width=900, res=300)
dens.scatter(x$d, x[, args[4]], parmar=c(3, 3, 1, 0.5), pt.cex=0.25)
axis(1)
mtext(side=1, text="Scores: File 2 - File 1", line=2)
axis(2)
mtext(side=2, text=args[4], line=2)
mtext(side=3, text=args[4], font=2)
dev.off()
