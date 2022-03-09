#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Approximate P-value equivalent of FDR<1% from sumstats BED file

options(scipen=1000, stringsAsFactors=F)
require(rCNV2, quietly=T)

# Load sumstats from BED file
ss <- read.table(commandArgs(trailingOnly=T)[1], header=T, sep="\t", check.names=F, comment.char="")

# Find largest P value that has FDR Q-value <=1%
q1 <- max(10^-ss$meta_neg_log10_p[which(ss$meta_neg_log10_fdr_q >= -log10(0.01))], na.rm=T)

# Find smallest P-value that has FDR Q-value >1%
q2 <- min(10^-ss$meta_neg_log10_p[which(ss$meta_neg_log10_fdr_q < -log10(0.01))], na.rm=T)

# Return average of q1 and q2 as estimate of Q-value corresponding exactly to
cat(paste(mean(q1, q2), "\n", sep=""))
