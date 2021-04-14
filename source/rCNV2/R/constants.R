#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Load project-wide constants


#' Load rCNV2 environment variables
#'
#' Helper function to load all constants and variables used in rCNV2 analyses
#'
#' @param load.strings boolean indicator to load character vector constants
#' @param load.colors boolean indicator to load color encodings
#' @param load.sclaes boolean indicator to load scale-related constants
#'
#' @details `load.rcnv.env()` will export constants from an rCNV2-specific
#' environment into the global environment and can be referenced as any other
#' variable within your R session.
#'
#' @seealso .GlobalEnv
#'
#' @export
load.rcnv.env <- function(load.strings=TRUE, load.colors=TRUE, load.scales=TRUE){

  # Create environment to contain all variables
  rcnv.env <- new.env()

  # Character vectors
  if(load.strings == TRUE){
    # Define constants
    constants <- list(
      hpo.abbrevs = c("HEALTHY_CONTROL" = "Control",
                      "HP:0000118" = "All cases",
                      "HP:0000707" = "Nervous system",
                      "HP:0012638" = "Nervous sys. physiol.",
                      "HP:0000708" = "Behavioral",
                      "UNKNOWN" = "Unspecified",
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
                      "HP:0011446" = "Cognition"),
      cohort.abbrevs = c("meta1" = "Cohort 1",
                         "meta2" = "Cohort 2",
                         "meta3" = "Cohort 3",
                         "meta4" = "Cohort 4",
                         "meta5" = "Cohort 5"),
      neuro.hpos = c("HP:0000707", "HP:0012638", "HP:0000708",
                     "HP:0100753", "HP:0031466", "HP:0100852",
                     "HP:0000729", "HP:0000717", "HP:0100022",
                     "HP:0000752", "HP:0012759", "HP:0011446",
                     "HP:0001250", "HP:0012639", "HP:0002011",
                     "HP:0012443"),
      somatic.hpos = c("HP:0000152", "HP:0003011", "HP:0001507",
                       "HP:0000924", "HP:0002715", "HP:0002960",
                       "HP:0025031", "HP:0001626", "HP:0001627",
                       "HP:0002597", "HP:0100545", "HP:0001197"),
      pheno.abbrevs = c("Mixed", "Neuro.", "Non-Neuro."),
      ml.model.abbrevs = c("ensemble" = "Ensemble",
                           "lda" = "LDA",
                           "logit" = "Logistic Regression",
                           "naivebayes" = "Naive Bayes",
                           "neuralnet" = "Neural Net",
                           "randomforest" = "Random Forest",
                           "sgd" = "SGD",
                           "svm" = "SVM"),
      nc.anno.family.names = c("chromhmm" = "Inferred Chromatin States",
                               "dhs" = "Open Chromatin Regions",
                               "histone" = "Histone Modification Peaks",
                               "tfbs" = "Transcription Factor Sites",
                               "transcription" = "Transcribed Elements",
                               "tads" = "Chromatin Domain Boundaries",
                               "enhancers" = "Enhancers",
                               "super.enhancers" = "Super Enhancers",
                               "other" = "Other")
    )

    # Assign constants to rCNV2 environment
    for(variable in names(constants)){
      assign(variable, constants[[variable]], envir=rcnv.env)
    }
  }

  # Colors
  if(load.colors == TRUE){

    require(viridisLite, quietly=TRUE)

    # Define subset of constants recursively referenced in constants list
    ns.color <- "gray70"
    ns.color.light <- "#F1F1F1"
    cnv.colors <- c("DEL" = "#D43925",
                   "DUP" = "#2376B2",
                   "CNV" = "#7E4EB2")
    redblack <- "#4F1C14"
    blueblack <- "#003F6A"
    purpleblack <- "#3F2759"
    redwhite <- "#F0D6D3"
    bluewhite <-"#E8F3FB"
    purplewhite <- "#EADFF5"
    lof.color <- "#9D1309"
    mis.color <- "#FF6103"
    syn.color <- "#AAAAAA"

    # Define constants
    constants <- list(
      graphabs.green = "#027831",
      graphabs.green.darker = "#235020",
      gw.sig.color = "#FFB533",
      ns.color.dark = "gray50",
      ns.color = ns.color,
      ns.color.light = ns.color.light,
      highlight.color = "#FFCB00",
      cnv.colors = cnv.colors,
      control.cnv.colors = c("DEL" = "#E69186",
                             "DUP" = "#79AACC",
                             "CNV" = "#B488A1"),
      redblack = redblack,
      blueblack = blueblack,
      purpleblack = purpleblack,
      cnv.blacks = c("DEL" = redblack,
                     "DUP" = blueblack,
                     "CNV" = purpleblack),
      redwhite = redwhite,
      bluewhite = bluewhite,
      purplewhite = purplewhite,
      cnv.whites = c("DEL" = redwhite,
                     "DUP" = bluewhite,
                     "CNV" = purplewhite),
      cnv.color.palettes = list("DEL" = colorRampPalette(c("gray95", cnv.colors[1]))(101),
                                "DUP" = colorRampPalette(c("gray95", cnv.colors[2]))(101),
                                "CNV" = colorRampPalette(c("gray95", cnv.colors[3]))(101)),
      ds.gradient.pal = colorRampPalette(c(ns.color, ns.color.light, purplewhite, cnv.colors[3]))(101),
      hits.gradient.pal = colorRampPalette(c(cnv.colors[1], ns.color.light, cnv.colors[2]))(101),
      lof.color = lof.color,
      mis.color = mis.color,
      syn.color = syn.color,
      snv.colors = c("lof" = lof.color,
                     "mis" = mis.color,
                     "syn" = syn.color),
      percentile.palette = viridis(101),
      cohort.colors = c("meta1" = "#0A5180",
                        "meta2" = "#1174B9",
                        "meta3" = "#51ACE8",
                        "meta4" = "#A3D3F2"),
      pheno.colors = c("all" = "#808080",
                       "neuro" = "#F58F38",
                       "somatic" = "#854614"),
      gene.feat.category.colors = c("genomic" = "#027831",
                                    "expression" = "#FFA300",
                                    "chromatin" = "#6D3D84",
                                    "constraint" = "#F6313E"),
      h3k27ac.color = "#FFC34D",
      nc.anno.family.colors = c("chromhmm" = "#490C65",
                                "dhs" = "#D3441C",
                                "histone" = "#BA7FD0",
                                "tfbs" = "#46A040",
                                "transcription" = "#00441B",
                                "tads" = "#01AF99",
                                "enhancers" = "#FFA300",
                                "super.enhancers" = "#F6313E",
                                "other" = ns.color)
    )

    # Assign constants to rCNV2 environment
    for(variable in names(constants)){
      assign(variable, constants[[variable]], envir=rcnv.env)
    }
  }

  # Scales
  if(load.scales == TRUE){

    # Define constants
    logscale.major <- 10^(-10:10)
    constants <- list(
      logscale.major = logscale.major,
      logscale.major.bp = 10^(0:9),
      logscale.major.bp.labels = c(sapply(c("bp", "kb", "Mb"),
                                          function(suf){paste(c(1, 10, 100), suf, sep="")}),
                                   "1 Gb"),
      logscale.demi = as.numeric(sapply(logscale.major, function(e){c(1, 5)*e})),
      logscale.demi.bp = as.numeric(sapply(10^(0:9), function(e){c(1, 5)*e})),
      logscale.demi.bp.labels = c(paste(c(1, 5, 10, 50, 100, 500), "bp", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "kb", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "Mb", sep=""),
                                  paste(c(1, 5), "Gb", sep="")),
      logscale.minor = as.numeric(sapply(logscale.major, function(e){(1:9)*e}))
    )

    # Assign constants to rCNV2 environment
    for(variable in names(constants)){
      assign(variable, constants[[variable]], envir=rcnv.env)
    }
  }

  # Assign all constants from rcnv.env to global environment for use by the user
  for(variable in names(rcnv.env)){
    assign(variable, rcnv.env[[variable]], .GlobalEnv)
  }
}

