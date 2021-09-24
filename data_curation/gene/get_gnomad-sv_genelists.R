#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Process gnomAD-SV v2.1 gene count data to generate gene lists


options(scipen=10000, stringsAsFactors=F)


# Read positional arguments
args <- commandArgs(trailingOnly=T)

# # DEV PARAMETERS
# args <- c("~/scratch/gencode.v19.canonical.pext_filtered.variation_features.bed.gz",
#           "~/scratch/gene_lists/gencode.v19.canonical.pext_filtered.genes.list",
#           "~/scratch/gnomad_sv.v2.1.nonneuro")

# Load data
g <- read.table(as.character(args[1]), sep="\t", comment.char="", header=T)
g <- g[order(g$gene), ]
elig <- read.table(as.character(args[2]))[, 1]
out.prefix <- as.character(args[3])

# Gather genes without any LoF deletions
no.lofdel <- g$gene[which(g$gnomad_sv_lof_del==0 & g$gene %in% elig)]
write.table(no.lofdel,
            paste(out.prefix, "no_lof_dels.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)

# Gather genes with at least one LoF deletion
yes.lofdel <- g$gene[which(g$gnomad_sv_lof_del>0 & g$gene %in% elig)]
write.table(yes.lofdel,
            paste(out.prefix, "has_lof_dels.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)

# Gather genes without any CG duplications
no.cg <- g$gene[which(g$gnomad_sv_cg==0 & g$gene %in% elig)]
write.table(no.cg,
            paste(out.prefix, "no_cg_dups.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)

# Gather genes with at least one CG duplication
yes.cg <- g$gene[which(g$gnomad_sv_cg>0 & g$gene %in% elig)]
write.table(yes.cg,
            paste(out.prefix, "has_cg_dups.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)
