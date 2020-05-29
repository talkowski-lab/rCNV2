#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Reformat gnomAD v2.1 mutation rates per gene


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
# Parse command-line arguments
args <- commandArgs(trailingOnly=T)
tsv.in <- as.character(args[1])
genes.in <- as.character(args[2])
tsv.out <- as.character(args[3])

# # Local dev parameters
# tsv.in <- "~/scratch/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
# genes.in <- "~/scratch/gencode.v19.canonical.pext_filtered.genes.list"
# tsv.out <- "~/scratch/cleaned_constraint_test.tsv"

# Read constraint metrics
df <- read.table(tsv.in, header=T, sep="\t", comment.char="", check.names=F)
df <- df[, which(colnames(df) %in% c("gene", "mu_lof", "mu_mis", "mu_syn"))]

# Read gene whitelist
genes <- read.table(genes.in, header=F, col.names="gene")

# Merge gene whitelist with mutation rates
x <- merge(genes, df, by="gene", all.x=T, all.y=F, sort=F)

# Normalize all mutation rates so that each column sums to P=1
x[, -1] <- apply(x[, -1], 2, function(vals){vals / sum(vals, na.rm=T)})

# Reorder columns according to mutational consequence
x <- x[c("gene", "mu_lof", "mu_mis", "mu_syn")]

# Write to tsv.out
colnames(x)[1] <- "#gene"
write.table(x, tsv.out, col.names=T, row.names=F, sep="\t", quote=F)
