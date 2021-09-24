#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parse ClinGen dosage sensitive regions

options(scipen=1000, stringsAsFactors=F)

# Read positional arguments
args <- commandArgs(trailingOnly=T)
in.tsv <- as.character(args[1])
out.prefix <- as.character(args[2])

# Load ClinGen data
dat <- read.table(in.tsv, header=T, skip=5, comment.char="", sep="\t")
cols.to.keep <- c("ISCA.Region.Name", "Genomic.Location", "Haploinsufficiency.Score", "Haploinsufficiency.Description",
                  "Triplosensitivity.Score", "Triplosensitivity.Description")
dat <- dat[, which(colnames(dat) %in% cols.to.keep)]

# Filter & reformat ClinGen data
format.coords <- function(cstr){
  cparts <- unlist(strsplit(cstr, split=":"))
  chrom <- gsub("chr", "", cparts[1])
  start <- as.numeric(unlist(strsplit(cparts[2], split="-"))[1])
  end <- as.numeric(unlist(strsplit(cparts[2], split="-"))[2])
  return(c(chrom, start, end))
}
extract.regions <- function(dat, cnv, elig.scores){
  if(cnv=="DEL"){
    score.idx <- which(colnames(dat)=="Haploinsufficiency.Score")
  }else if(cnv=="DUP"){
    score.idx <- which(colnames(dat)=="Triplosensitivity.Score")
  }
  dat <- dat[which(dat[, score.idx] %in% elig.scores), ]
  scores <- dat[, score.idx]
  coords <- as.data.frame(t(sapply(dat$Genomic.Location, format.coords)))
  res <- as.data.frame(cbind(coords, rep(cnv, nrow(dat)), 
                             dat$ISCA.Region.Name, scores))
  colnames(res) <- c("#chr", "start", "end", "cnv", "disorder", "score")
  rownames(res) <- NULL
  return(res)
}
write.table(extract.regions(dat, "DEL", 1:3),
            paste(out.prefix, "DEL_GD_all.bed", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
write.table(extract.regions(dat, "DEL", 2:3),
            paste(out.prefix, "DEL_GD_hmc.bed", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
write.table(extract.regions(dat, "DUP", 1:3),
            paste(out.prefix, "DUP_GD_all.bed", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
write.table(extract.regions(dat, "DUP", 2:3),
            paste(out.prefix, "DUP_GD_hmc.bed", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
