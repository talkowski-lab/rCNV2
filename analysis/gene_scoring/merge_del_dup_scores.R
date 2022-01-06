#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Merge del & dup scores


options(scipen=100000, stringsAsFactors=F)

# Read command line args
args <- commandArgs(trailingOnly=T)
del <- read.table(args[1], header=T, sep="\t")
dup <- read.table(args[2], header=T, sep="\t")

# Merge & reformat scores
x <- merge(del, dup, by="gene", suffix=c(".del", ".dup"), all=T, sort=F)
xo <- data.frame("gene"=x$gene, "pHaplo"=x$score.del, "pTriplo"=x$score.dup)
xo <- xo[rev(order(apply(xo[, 2:3], 1, max, na.rm=T))), ]

# Write out
colnames(xo)[1] <- "#gene"
write.table(xo, args[3], col.names=T, row.names=F, sep="\t", quote=F)
