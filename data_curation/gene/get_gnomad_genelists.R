#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Process gnomAD v2.1.1 constraint data to generate gene lists


options(scipen=10000, stringsAsFactors=F)


# Read positional arguments
args <- commandArgs(trailingOnly=T)

# # DEV PARAMETERS
# args <- c("~/scratch/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
#           "~/scratch/gene_lists/gencode.v19.canonical.pext_filtered.genes.list",
#           "~/scratch/gnomad.v2.1.1")

# Load data
g <- read.table(as.character(args[1]), sep="\t", comment.char="", header=T)
g <- g[order(g$gene), ]
g$oe_mis_upper_bin_6 <- floor((6*rank(g$oe_mis_upper)) / nrow(g))
elig <- read.table(as.character(args[2]))[, 1]
out.prefix <- as.character(args[3])

# Gather LoF constrained genes
lof.constrained <- g$gene[which((g$pLI >= 0.9 | g$oe_lof_upper_bin_6 == 0) & g$gene %in% elig)]
write.table(lof.constrained,
            paste(out.prefix, "lof_constrained.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)

# Gather missense constrained genes
mis.constrained <- g$gene[which((g$mis_z >= 3 | g$oe_mis_upper_bin_6 == 0) & g$gene %in% elig)]
write.table(mis.constrained,
            paste(out.prefix, "mis_constrained.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)

# Gather likely unconstrained genes
likely.unconstrained <- g$gene[which(g$oe_lof_upper >= 1 & g$oe_mis_upper >= 1 & 
                                       g$pLI <= 0.1 & g$syn_z >= -3 & g$syn_z <= 3 & 
                                       g$oe_lof >= median(g$oe_lof, na.rm=T) & 
                                       g$oe_mis >= median(g$oe_lof, na.rm=T) & 
                                       g$obs_lof > 0 & g$obs_mis > 0)]
write.table(likely.unconstrained,
            paste(out.prefix, "likely_unconstrained.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)

# Gather high-confidence mutationally tolerant genes
mut.tolerant <- g$gene[which(g$pLI <= 0.01 & g$oe_lof_upper_bin_6 >= 4 &
                               g$mis_z <= 0 & g$oe_mis_upper >= 1 &
                               g$syn_z >= -3 & g$syn_z <= 3)]
write.table(mut.tolerant,
            paste(out.prefix, "mutation_tolerant.genes.list", sep="."),
            col.names=F, row.names=F, quote=F)
