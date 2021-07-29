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
#' @param color.by indicator to color by 'neuro' or 'severity'
#'
#' @return color code
#'
#' @export
get.hpo.color <- function(hpo, color.by="neuro"){
  # Ensure rCNV2 constants are loaded
  load.rcnv.env()

  # Color differently based on "color.by"
  if(color.by == "neuro"){
    if(hpo %in% neuro.hpos){
      pheno.colors[which(names(pheno.colors) == "neuro")]
    }else if(hpo %in% somatic.hpos){
      pheno.colors[which(names(pheno.colors) == "somatic")]
    }else{
      pheno.colors[which(names(pheno.colors) == "all")]
    }
  }else if(color.by == "severity"){
    if(hpo %in% developmental.hpos){
      severity.colors[which(names(severity.colors) == "developmental")]
    }else if(hpo %in% adult.hpos){
      severity.colors[which(names(severity.colors) == "adult")]
    }else{
      severity.colors[which(names(severity.colors) == "all")]
    }
  }else{
    stop(paste("get.hpo.color does not recognize color.by =", color.by))
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
#'
#' @export format.pval
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
#'
#' @export format.scientific
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
#' @seealso [viridis()]
#'
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

#' Density-colored scatterplot
#'
#' Generate a 2D scatterplot colored by local density
#'
#' @param x independent variable vector
#' @param y dependent variable vector
#' @param pt.cex scaling factor for points
#' @param parmar argument for option(parmar=...)
#' @param add.cor boolean indicator to add correlation coefficient to plot \[default: TRUE\]
#'
#' @return None
#'
#' @seealso [color.points.by.density()]
#'
#' @export dens.scatter
#' @export
dens.scatter <- function(x, y, pt.cex=1, parmar=rep(0.5, 4), add.cor=T){
  # Based on heatscatter.R from Colby Chiang
  # (https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R)
  par(mar=parmar, bty="o")
  plot.df <- data.frame("x"=x, "y"=y)
  plot.df <- plot.df[which(!is.infinite(plot.df$x) & !is.infinite(plot.df$y)
                           & !is.na(plot.df$x) & !is.na(plot.df$y)), ]
  dens <- densCols(plot.df$x, plot.df$y, colramp=colorRampPalette(c("black", "white")))
  plot.df$dens <- col2rgb(dens)[1, ] + 1L
  palette <- colorRampPalette(c("#440154", "#33638d", "#218f8d",
                                "#56c667", "#fde725"))(256)
  plot.df$col <- palette[plot.df$dens]
  plot.df <- plot.df[order(plot.df$dens), ]
  plot(plot.df$x, plot.df$y, type="n",
       xlab="", xaxt="n", ylab="", yaxt="n")
  abline(lm(y ~ x, data=plot.df), lwd=2)
  points(plot.df$x, plot.df$y, col=plot.df$col, pch=19, cex=pt.cex)
  if(add.cor==T){
    r <- cor(plot.df$x, plot.df$y)
    text(x=par("usr")[2], y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])),
         labels=bquote(italic(R) == .(formatC(round(r, 3), digits=3))), pos=2)
  }
}

