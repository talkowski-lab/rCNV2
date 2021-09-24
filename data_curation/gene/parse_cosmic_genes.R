#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parse COSMIC gene census list

# Set parameters
options(scipen=1000, stringsAsFactors=F)
gof.mutations <- c("A", "Mis", "O")

# Read positional arguments
args <- commandArgs(trailingOnly=T)
in.tsv <- as.character(args[1])
elig.in <- as.character(args[2])
outdir <- as.character(args[3])

# Load COSMIC data
dat <- read.table(in.tsv, header=T, sep=",")
colnames(dat)[1] <- "gene"
cols.to.keep <- c("gene", "Tier", "Molecular.Genetics", "Role.in.Cancer", 
                  "Mutation.Types", "Synonyms")
dat <- dat[, which(colnames(dat) %in% cols.to.keep)]

# Load list of eligible genes
elig <- read.table(elig.in, header=F, sep="\t")[, 1]

# Process oncogenes
hc.oncogene.idx <- intersect(intersect(intersect(which(dat$Tier == 1 & dat$gene %in% elig),
                          union(grep("Dom", dat$Molecular.Genetics, fixed=T), 
                               which(dat$Molecular.Genetics==""))),
                          grep("oncogene", dat$Role.in.Cancer, fixed=T)),
                          which(sapply(dat$Mutation.Types, function(mstr){length(intersect(gof.mutations, unlist(strsplit(mstr, split=", "))))}) > 0))
write.table(sort(dat$gene[hc.oncogene.idx]),
            paste(outdir, "COSMIC.hc_oncogenes.genes.list", sep="/"),
            col.names=F, row.names=F, sep="\t", quote=F)

