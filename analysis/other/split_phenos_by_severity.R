#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Bifurcate HPOs by effect size of whole-gene deletions of constrained genes


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=T)


# Get command-line arguments & options
option_list <- list(
  make_option(c("-p", "--out-prefix"), type="character", default="rCNV2_pheno_split",
              help="prefix for output files [default %default]", metavar="string")
)
args <- parse_args(OptionParser(usage="%prog data",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options
