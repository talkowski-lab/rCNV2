#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Identify windows with CNVs found at a much higher frequency in one cohort than any other


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--max-diff"), type="numeric", default=0.01, 
              help="Maxmimum difference in CNV frequency between highest and second-highest to flag. [default %default]"),
  make_option(c("--min-top"), type="numeric", default=0.01, 
              help="Minimum frequency in top cohort to flag. [default %default]"),
  make_option(c("--max-second"), type="numeric", default=0.01, 
              help="Maximum frequency in second cohort to flag. [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog freq_matrix.bed.gz outfile.bed",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Three positional arguments: freq_matrix.bed.gz and outfile.bedx\n")
}



