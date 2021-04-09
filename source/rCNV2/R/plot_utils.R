#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Project-wide utility functions for plotting


#' Get HPO-based color code
#'
#' Return hex color code for phenotype by HPO
#'
#' @param hpo HPO code
#'
#' @return color code
#' @export
get.hpo.color <- function(hpo){
  # TODO: automatically load pheno.colors from constants.R within this function
  if(hpo %in% neuro.hpos){
    pheno.colors[which(names(pheno.colors) == "neuro")]
  }else if(hpo %in% somatic.hpos){
    pheno.colors[which(names(pheno.colors) == "somatic")]
  }else{
    pheno.colors[which(names(pheno.colors) == "all")]
  }
}


#' Format P-value
#'
#' Format P-value for plotting
#'
#' @param p P-value
#' @param nsmall number of digits after the decimal to retain for scientific
#' notification \[default: 2\]
#' @param max.decimal convert all P-values requiring more digits after the decimal
#' to be converted to scientific notation \[default: 3\]
#' @param equality equality symbol to print after `P` \[default: '='\]
#' @param min.phred.p minimum order of magnitude to process before considering
#' P-value to be arbitrarily/meaninglessly small \[default: 100\]
#'
#' @return formatted P-value as character
#' @export
format.pval <- function(p, nsmall=2, max.decimal=3, equality="=", min.phred.p=100){
  if(-log10(p)>min.phred.p){
    bquote(italic(P) %~~% 0)
  }else if(ceiling(-log10(p)) > max.decimal){
    parts <- unlist(strsplit(format(p, scientific=T), split="e"))
    base <- gsub(" ", "", formatC(round(as.numeric(parts[1]), nsmall), digits=1+nsmall), fixed=T)
    exp <- gsub(" ", "", as.numeric(parts[2]), fixed=T)
    bquote(italic(P) ~ .(equality) ~ .(base) ~ "x" ~ 10 ^ .(exp))
  }else{
    bquote(italic(P) ~ .(equality) ~ .(formatC(round(p, max.decimal), digits=max.decimal)))
  }
}


#' Format scientific value
#'
#' Format a scientific value for plotting
#'
#' @param x numeric value to be formatted
#' @param nsmall number of digits after the decimal to retain for scientific
#' notification \[default: 2\]
#' @param max.decimal convert all values requiring more digits after the decimal
#' to be converted to scientific notation \[default: 3\]
#'
#' @return formatted value as character
#' @export
format.scientific <- function(x, nsmall=2, max.decimal=3){
  parts <- unlist(strsplit(format(x, scientific=T), split="e", fixed=T))
  base <- format(round(as.numeric(parts[1]), max.decimal), nsmall=nsmall)
  exp <- as.character(as.numeric(parts[2]))
  bquote(.(base) ~ "x" ~ 10 ^ .(exp))
}


#' Color points by density
#'
#' Generate colors for XY scatterplot based on point density
#'
#' @param x independent variable vector
#' @param y dependent variable vector
#' @param palette 256-color palette to be applied based on density \[default: `viridis()`\]
#'
#' @details Inspired by heatscatter.R from Colby Chiang:
#'  https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R
#'
#' @return dataframe of values to be plotted with density and colors
#'
#' @seealso `viridis()`
#' @export
color.points.by.density <- function(x, y, palette=NULL){
  # Based on heatscatter.R from Colby Chiang
  # (https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R)
  plot.df <- data.frame("x"=x, "y"=y)
  plot.df <- plot.df[which(!is.infinite(plot.df$x) & !is.infinite(plot.df$y)
                           & !is.na(plot.df$x) & !is.na(plot.df$y)), ]
  dens <- densCols(plot.df$x, plot.df$y, colramp=colorRampPalette(c("black", "white")))
  plot.df$dens <- col2rgb(dens)[1, ] + 1L
  if(is.null(palette)){
    require(viridis, quietly=TRUE)
    palette <- viridis(256)
  }
  plot.df$col <- palette[plot.df$dens]
  plot.df[order(plot.df$dens), ]
}
