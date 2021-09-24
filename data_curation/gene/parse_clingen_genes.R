#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parse ClinGen gene curation list

options(scipen=1000, stringsAsFactors=F)

# Read positional arguments
args <- commandArgs(trailingOnly=T)
in.tsv <- as.character(args[1])
elig.in <- as.character(args[2])
outdir <- as.character(args[3])

# Load ClinGen data
dat <- read.table(in.tsv, header=T, skip=5, comment.char="", sep="\t")
colnames(dat)[1] <- "gene"
cols.to.keep <- c("gene", "Haploinsufficiency.Score", "Haploinsufficiency.Description",
                  "Triplosensitivity.Score", "Triplosensitivity.Description")
dat <- dat[, which(colnames(dat) %in% cols.to.keep)]

# Load list of eligible genes
elig <- read.table(elig.in, header=F, sep="\t")[, 1]

# Process haploinsufficient genes
write.table(sort(dat$gene[which(dat$Haploinsufficiency.Score == 3 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.hc_haploinsufficient.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Haploinsufficiency.Score == 2 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.mc_haploinsufficient.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Haploinsufficiency.Score %in% 2:3 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.hmc_haploinsufficient.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Haploinsufficiency.Score == 1 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.lc_haploinsufficient.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Haploinsufficiency.Score %in% 1:3 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.all_haploinsufficient.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)

# Process triplosensitive genes
write.table(sort(dat$gene[which(dat$Triplosensitivity.Score == 3 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.hc_triplosensitive.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Triplosensitivity.Score == 2 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.mc_triplosensitive.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Triplosensitivity.Score %in% 2:3 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.hmc_triplosensitive.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Triplosensitivity.Score == 1 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.lc_triplosensitive.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)
write.table(sort(dat$gene[which(dat$Triplosensitivity.Score %in% 1:3 & dat$gene %in% elig)]),
            paste(outdir, "ClinGen.all_triplosensitive.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)

