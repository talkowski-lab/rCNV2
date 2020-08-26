#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parse ASC/SPARK pedigree file


options(stringsAsFactors=F, scipen=1000, family="sans")


# Constants
child.roles <- c("brother", "Child", "DD-ID_Proband", "Other_Proband",
                 "P", "proband", "Proband", "PROBAND", "Sibling")
p.map <- c("SVIP" = 0,
           "Control" = 1,
           "ASD" = 2,
           "DD" = 3,
           "ASD_DD" = 4,
           "Unknown" = 5)


######################
### DATA FUNCTIONS ###
######################
# Filter pedigree table
filter.ped <- function(p){
  keep.idx <- which(p$Sex_CN %in% c("XX", "XY")
                    & p$Affected_Status %in% c(1:4)
                    & !is.na(p$Person)
                    & p$Role %in% child.roles)
  p[keep.idx, ]
}

# Make table of child IDs with phenotypes
make.pheno.table <- function(p){
  phenos <- sapply(p$Affected_Status, function(k){names(p.map)[which(p.map==k)]})
  s.out <- data.frame("child_id" = p$Person,
             "phenotype" = phenos)
  colnames(s.out)[1] <- "#child_id"
  s.out[order(s.out), ]
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog infile.csv outfile", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: infile.csv and outfile\n", sep=" "))
}

# Writes args & opts to vars
ped.in <- args$args[1]
outfile <- args$args[2]

# # DEV PARAMETERS
# ped.in <- "~/scratch/asc_spark_denovo_cnvs/pedigree_cnv_07_27_2020.csv"
# outfile <- "~/scratch/asc_spark_child_phenos.list"

# Load pedigree .csv
p <- read.table(ped.in, header=T, sep=",", comment.char="")

# Filter pedigree
p <- filter.ped(p)

# Write out list of child IDs with phenotype
s.out <- make.pheno.table(p)
write.table(s.out, outfile, sep="\t", quote=F, col.names=T, row.names=F)
