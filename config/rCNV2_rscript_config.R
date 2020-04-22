#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master parameters for rCNV2 manuscript formalized secondary analyses

options(scipen=1000, stringsAsFactors=F)

cnv.colors <- c("DEL" = "#D43925",
                "DUP" = "#2376B2")

control.cnv.colors <- c("DEL" = "#E69186",
                        "DUP" = "#79AACC")

hpo.abbrevs <- c("HEALTHY_CONTROL" = "Control", 
                 "HP:0000118" = "All phenotypes", 
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
