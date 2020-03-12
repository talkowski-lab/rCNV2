#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parse DECIPHER/DDG2P gene list

options(scipen=1000, stringsAsFactors=F)

# Read positional arguments
args <- commandArgs(trailingOnly=T)
in.csv <- as.character(args[1])
elig.in <- as.character(args[2])
outdir <- as.character(args[3])

# Load list of eligible genes
elig <- read.table(elig.in, header=F, sep="\t")[, 1]

# Load DDG2P genes and filter to eligible set
dat <- read.table(in.csv, sep=",", header=T)
colnames(dat)[1] <- "gene"
cols.to.keep <- c("gene", "DDD.category", "allelic.requirement", 
                  "mutation.consequence", "phenotypes")
dat <- dat[which(dat$gene %in% elig), which(colnames(dat) %in% cols.to.keep)]
colnames(dat)[ncol(dat)] <- c("hpos")

# Restrict to dominant genes
dat <- dat[grep("monoallelic|uncertain", dat$allelic.requirement, fixed=F), ]

# Function to extract a subset of genes based on mechanism & confidence
mech.map <- list("lof" = c("dominant negative", "loss of function"),
                 "gof" = c("activating", "gain of function", "increased gene dosage",
                           "part of contiguous gene duplication"),
                 "other" = c("all missense/in frame", "uncertain"))
conf.map <- list("hc"=c("confirmed"), "mc"=c("probable"), "hmc"=c("confirmed", "probable"), 
                 "lc"=c("possible"), "all"=c("confirmed", "probable", "possible"))
filter.genes <- function(dat, mech, conf){
  keep.idx <- intersect(which(dat$mutation.consequence %in% mech.map[[mech]]),
            which(dat$DDD.category %in% conf.map[[conf]]))
  return(dat[keep.idx, ])
}

# Write all gene lists to tables
for(mech in names(mech.map)){
  for(conf in names(conf.map)){
    gset <- filter.genes(dat, mech, conf)
    write.table(sort(gset$gene),
                paste(outdir, "/DDG2P.", conf, "_", mech, ".genes.list", sep=""),
                col.names=F, row.names=F, sep="\t", quote=F)
    colnames(gset)[1] <- "#gene"
    gset <- gset[which(gset$hpos != ""), ]
    write.table(gset[order(gset$`#gene`), which(colnames(gset) %in% c("#gene", "hpos"))],
                paste(outdir, "/DDG2P.", conf, "_", mech, ".genes_with_hpos.tsv", sep=""),
                col.names=T, row.names=F, sep="\t", quote=F)
  }
}
