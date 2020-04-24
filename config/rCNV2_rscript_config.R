#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master parameters for rCNV2 manuscript formalized secondary analyses

options(scipen=1000, stringsAsFactors=F)


########
# DATA #
########
hpo.abbrevs <- c("HEALTHY_CONTROL" = "Control", 
                 "HP:0000118" = "All cases", 
                 "HP:0000707" = "Nervous system", 
                 "HP:0012638" = "Nervous sys. physiol.", 
                 "HP:0000708" = "Behavioral", 
                 "UNKNOWN" = "Other", 
                 "HP:0012639" = "Nervous sys. morph.", 
                 "HP:0002715" = "Immune", 
                 "HP:0012759" = "Neurodevelopmental", 
                 "HP:0002960" = "Autoimmune", 
                 "HP:0002011" = "CNS abnormality", 
                 "HP:0001626" = "Cardiovascular", 
                 "HP:0100753" = "Schizophrenia", 
                 "HP:0000729" = "Autistic behavior", 
                 "HP:0002597" = "Vascular", 
                 "HP:0100022" = "Movement", 
                 "HP:0001250" = "Seizures", 
                 "HP:0100545" = "Arterial", 
                 "HP:0000717" = "Autism", 
                 "HP:0000752" = "Hyperactivity", 
                 "HP:0001197" = "Birth defects", 
                 "HP:0000924" = "Skeletal", 
                 "HP:0031466" = "Personality", 
                 "HP:0000152" = "Head/Neck", 
                 "HP:0001627" = "Cardiac", 
                 "HP:0025031" = "Digestive", 
                 "HP:0001507" = "Growth abnormality", 
                 "HP:0100852" = "Anxiety", 
                 "HP:0012443" = "Brain morphology", 
                 "HP:0003011" = "Muscle", 
                 "HP:0011446" = "Cognition")

cohort.abbrevs <- c("meta1" = "Cohort 1",
                    "meta2" = "Cohort 2",
                    "meta3" = "Cohort 3",
                    "meta4" = "Cohort 4")

neuro.hpos <- c("HP:0000707", "HP:0012638", "HP:0000708", 
                "HP:0100753", "HP:0031466", "HP:0100852", 
                "HP:0000729", "HP:0000717", "HP:0100022", 
                "HP:0000752", "HP:0012759", "HP:0011446", 
                "HP:0001250", "HP:0012639", "HP:0002011", 
                "HP:0012443")

somatic.hpos <- c("HP:0000152", "HP:0003011", "HP:0001507", 
                  "HP:0000924", "HP:0002715", "HP:0002960", 
                  "HP:0025031", "HP:0001626", "HP:0001627", 
                  "HP:0002597", "HP:0100545", "HP:0001197")


##########
# COLORS #
##########
cnv.colors <- c("DEL" = "#D43925",
                "DUP" = "#2376B2")

control.cnv.colors <- c("DEL" = "#E69186",
                        "DUP" = "#79AACC")

require(viridisLite)
percentile.palette <- viridis(101)

blueblack <- "#003F6A"
bluewhite <- "#E8F3FB"

cohort.colors <- c("meta1" = "#0A5180",
                   "meta2" = "#1174B9",
                   "meta3" = "#51ACE8",
                   "meta4" = "#A3D3F2")


##########
# SCALES #
##########
logscale.major <- 10^(-10:10)
logscale.major.bp <- 10^(0:9)
logscale.major.bp.labels <- c(sapply(c("bp", "kb", "Mb"), 
                                     function(suf){paste(c(1, 10, 100), suf, sep="")}), 
                              "1 Gb")

logscale.demi <- as.numeric(sapply(logscale.major, function(e){c(1, 5)*e}))
logscale.demi.bp <- as.numeric(sapply(10^(0:9), function(e){c(1, 5)*e}))
logscale.demi.bp.labels <- c(paste(c(1, 5, 10, 50, 100, 500), "bp", sep=""),
                             paste(c(1, 5, 10, 50, 100, 500), "kb", sep=""),
                             paste(c(1, 5, 10, 50, 100, 500), "Mb", sep=""),
                             paste(c(1, 5), "Gb", sep=""))

logscale.minor <- as.numeric(sapply(logscale.major, function(e){(1:9)*e}))

  
  