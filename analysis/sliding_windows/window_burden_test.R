#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform a Fisher's exact test for rare CNV burden for sliding windows



#################
### RSCRIPT BLOCK
#################
require(optparse,quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--case-column"), type="character", default=NULL,
              help="name of column to use for case CNV counts [default %default]",
              metavar="character"),
  make_option(c("--control-column"), type="character", default=NULL,
              help="name of column to use for control CNV counts [default %default]",
              metavar="character")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog bins outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options
print(opts)

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}
if(is.null(opts$case_column)){
  stop("Must specify --case-column\n")
}
if(is.null(opts$control_column)){
  stop("Must specify --control-column\n")
}

# Writes args & opts to vars
bed.in <- args$args[1]
outfile <- args$args[2]
case_column <- opts$case_column
control_column <- opts$control_column

