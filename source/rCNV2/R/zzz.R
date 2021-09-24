#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Actions to perform upon loading rCNV2 package


.onLoad <- function(libname, pkgname){
  options(scipen=1000, stringsAsFactors=F, family="sans")
}

.onAttach <- function(libname, pkgname){
  load.rcnv.env()
}
