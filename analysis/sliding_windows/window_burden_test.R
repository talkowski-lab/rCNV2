#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform a Fisher's exact test for rare CNV burden for sliding windows


# Process an input bed file
import.bed <- function(bed.in, case.col.name, n.case, ctrl.col.name, n.ctrl){
  bed <- read.table(bed.in, sep="\t", header=T, comment.char="")
  
  case.col.idx <- grep(case.col.name, colnames(bed), fixed=T)
  if(length(case.col.idx) == 0){
    stop(paste("--case-column \"",case.col.name,"\" cannot be found ",
               "in BED header", sep=""))
  }
  
  ctrl.col.idx <- grep(ctrl.col.name, colnames(bed), fixed=T)
  if(length(ctrl.col.idx) == 0){
    stop(paste("--control-column \"",ctrl.col.idx,"\" cannot be found ",
               "in BED header", sep=""))
  }
  
  bed <- bed[,c(1:3, case.col.idx, ctrl.col.idx)]
  colnames(bed) <- c("chr", "start", "end", "case.CNV", "ctrl.CNV")
  bed$case.ref <- n.case - bed$case.CNV
  bed$ctrl.ref <- n.ctrl - bed$ctrl.CNV
  
  return(bed)
}


# Fisher's exact test for a single vector of CNV counts
burden.test.single <- function(counts){
  case.cnv <- as.integer(counts[1])
  ctrl.cnv <- as.integer(counts[2])
  case.ref <- as.integer(counts[3])
  ctrl.ref <- as.integer(counts[4])
  
  if(case.cnv == 0 & ctrl.cnv == 0){
    p <- 1
    or <- c(NA,NA,NA)
  }else{
    cnv.mat <- matrix(c(ctrl.ref, case.ref, ctrl.cnv, case.cnv),
                      byrow=T, nrow=2)
    
    p <- fisher.test(cnv.mat, alternative="greater")$p.value
    or <- fisher.test(cnv.mat)
    or <- c(or$estimate, or$conf.int)
  }
  
  f.res <- as.numeric(c(p, or))
  names(f.res) <- c("p", "OR", "OR.lower", "OR.upper")
  return(f.res)
}


# Fisher's exact test of case:control CNV burden per bin
burden.test <- function(bed){
  f.res <- t(apply(bed[ ,-c(1:3)], 1, burden.test.single))
  
  fisher.bed <- cbind(bed[ ,1:3], f.res)
  colnames(fisher.bed) <- c("chr", "start", "end",
                            "p", "OR", "OR.lower", "OR.upper")
  return(fisher.bed)
}


#################
### RSCRIPT BLOCK
#################
require(optparse,quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--case-column"), type="character", default=NULL,
              help="name of column to use for case CNV counts [default %default]",
              metavar="character"),
  make_option(c("--case-n"), type="integer", default=NULL,
              help="total case cohort sample size [default %default]",
              metavar="integer"),
  make_option(c("--control-column"), type="character", default=NULL,
              help="name of column to use for control CNV counts [default %default]",
              metavar="character"),
  make_option(c("--control-n"), type="integer", default=NULL,
              help="total control cohort sample size [default %default]",
              metavar="integer"),
  make_option(c("-z", "--bgzip"), action="store_true", default=FALSE,
              help="bgzip output BED file [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog bins outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}
if(is.null(opts$`case-column`)){
  stop("Must specify --case-column\n")
}
if(is.null(opts$`control-column`)){
  stop("Must specify --control-column\n")
}
if(is.null(opts$`case-n`)){
  stop("Must specify --case-n\n")
}
if(is.null(opts$`control-n`)){
  stop("Must specify --control-n\n")
}

# Writes args & opts to vars
bed.in <- args$args[1]
outfile <- args$args[2]
case.col.name <- opts$`case-column`
ctrl.col.name <- opts$`control-column`
n.case <- opts$`case-n`
n.ctrl <- opts$`control-n`

# Process input BED
bed <- import.bed(bed.in, case.col.name, n.case, ctrl.col.name, n.ctrl)

# Run burden tests
fisher.bed <- burden.test(bed)
  
# Format output
colnames(fisher.bed)[1] <- "#chr"
write.table(fisher.bed, outfile, sep="\t", col.names=T, row.names=F, quote=F)
if(opts$bgzip == TRUE){
  system(paste("bgzip -f", outfile),wait=T)
}
