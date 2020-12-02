#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master parameters and common functions for rCNV2 manuscript formalized secondary analyses

options(scipen=1000, stringsAsFactors=F, family="sans")


########
# DATA #
########
hpo.abbrevs <- c("HEALTHY_CONTROL" = "Control", 
                 "HP:0000118" = "All cases", 
                 "HP:0000707" = "Nervous system", 
                 "HP:0012638" = "Nervous sys. physiol.", 
                 "HP:0000708" = "Behavioral", 
                 "UNKNOWN" = "Unpecified", 
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

pheno.abbrevs <- c("Mixed", "Neuro.", "Non-Neuro.")

ml.model.abbrevs <- c("ensemble" = "Ensemble",
                      "lda" = "LDA",
                      "logit" = "Logistic Regression",
                      "naivebayes" = "Naive Bayes",
                      "neuralnet" = "Neural Net",
                      "randomforest" = "Random Forest",
                      "sgd" = "SGD",
                      "svm" = "SVM")

nc.anno.family.names <- c("chromhmm" = "Inferred Chromatin States",
                          "dhs" = "Open Chromatin Regions",
                          "histone" = "Histone Modification Peaks",
                          "tfbs" = "Transcription Factor Sites",
                          "transcription" = "Transcribed Elements",
                          "tads" = "Chromatin Domain Boundaries",
                          "enhancers" = "Enhancers",
                          "super.enhancers" = "Super Enhancers",
                          "other" = "Other")


##########
# COLORS #
##########
graphabs.green <- "#027831"
graphabs.green.darker <- "#235020"
gw.sig.color <- "#FFB533"
ns.color.dark <- "gray50"
ns.color <- "gray70"
ns.color.light <- "#F1F1F1"
highlight.color <- "#FFCB00"

cnv.colors <- c("DEL" = "#D43925",
                "DUP" = "#2376B2",
                "CNV" = "#7E4EB2")

control.cnv.colors <- c("DEL" = "#E69186",
                        "DUP" = "#79AACC",
                        "CNV" = "#B488A1")

blueblack <- "#003F6A"
redblack <- "#4F1C14"
purpleblack <- "#3F2759"
cnv.blacks <- c("DEL" = redblack,
                "DUP" = blueblack,
                "CNV" = purpleblack)

bluewhite <- "#E8F3FB"
redwhite <- "#F0D6D3"
purplewhite <- "#EADFF5"
cnv.whites <- c("DEL" = redwhite,
                "DUP" = bluewhite,
                "CNV" = purplewhite)

cnv.color.palettes <- list("DEL" = colorRampPalette(c("gray95", cnv.colors[1]))(101),
                           "DUP" = colorRampPalette(c("gray95", cnv.colors[2]))(101),
                           "CNV" = colorRampPalette(c("gray95", cnv.colors[3]))(101))

ds.gradient.pal <- colorRampPalette(c(ns.color, ns.color.light, purplewhite, cnv.colors[3]))(101)

hits.gradient.pal <- colorRampPalette(c(cnv.colors[1], ns.color.light, cnv.colors[2]))(101)

lof.color = "#9D1309"
mis.color = "#FF6103"
syn.color = "#AAAAAA"
snv.colors <- c("lof" = lof.color,
                "mis" = mis.color,
                "syn" = syn.color)

require(viridisLite)
percentile.palette <- viridis(101)

cohort.colors <- c("meta1" = "#0A5180",
                   "meta2" = "#1174B9",
                   "meta3" = "#51ACE8",
                   "meta4" = "#A3D3F2")

pheno.colors <- c("all" = "#808080",
                  "neuro" = "#F58F38",
                  "somatic" = "#854614")

gene.feat.category.colors <- c("genomic" = "#027831",
                               "expression" = "#FFA300",
                               "chromatin" = "#6D3D84",
                               "constraint" = "#F6313E")

h3k27ac.color <- "#FFC34D"

nc.anno.family.colors <- c("chromhmm" = "#490C65",
                           "dhs" = "#D3441C",
                           "histone" = "#BA7FD0",
                           "tfbs" = "#46A040",
                           "transcription" = "#00441B",
                           "tads" = "#01AF99",
                           "enhancers" = "#FFA300",
                           "super.enhancers" = "#F6313E",
                           "other" = ns.color)


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


#############
# FUNCTIONS #
#############
# Return hex color code for phenotype by HPO
get.hpo.color <- function(hpo){
  if(hpo %in% neuro.hpos){
    pheno.colors[which(names(pheno.colors) == "neuro")]
  }else if(hpo %in% somatic.hpos){
    pheno.colors[which(names(pheno.colors) == "somatic")]
  }else{
    pheno.colors[which(names(pheno.colors) == "all")]
  }
}

# Format p-value for printing to plots
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

# Calculate & format permuted P-value
calc.perm.p <- function(perm.vals, obs.val, alternative="greater"){
  if(alternative=="greater"){
    p <- length(which(perm.vals >= obs.val)) / length(perm.vals) 
  }else{
    p <- length(which(perm.vals <= obs.val)) / length(perm.vals)
  }
  if(p==0){
    p.fmt <- format.pval(1/length(perm.vals), equality="<")
  }else{
    p.fmt <- format.pval(p)
  }
  return(list("p" = p, "formatted" = p.fmt))
}

# Fit a robust linear regression with confidence interval
robust.lm <- function(x, y, conf=0.95){
  require(MASS, quietly=T)
  fit.df <- data.frame("Y"=y, "X"=x)
  xrange <- range(x[which(!is.infinite(x))], na.rm=T)
  xspan <- xrange[2] - xrange[1]
  fit <- MASS::rlm(Y ~ X, data=fit.df)
  pred.df <- data.frame("X"=seq(xrange[1] - 2*xspan, xrange[2] + 2*xspan, length.out=1000))
  pred <- predict(fit, pred.df, interval="confidence", level=conf)
  pred.out <- data.frame("x"=pred.df$X, "lower"=pred[, 2], "upper"=pred[, 3])
  return(list("fit" = fit, "ci" = pred.out))
}

# Fit exponential decay to a list of points
fit.exp.decay <- function(x, y){
  # Following example on https://rpubs.com/mengxu/exponential-model
  
  train.df <- data.frame("x"=as.numeric(x), 
                         "y"=as.numeric(y))
  train.df <- train.df[which(train.df$y > 0), ]
  
  # Estimate initial parameters
  theta.0 <- 0.5 * min(train.df$y)
  model.0 <- lm(log(y - theta.0) ~ x, data=train.df)
  alpha.0 <- exp(coef(model.0)[1])
  beta.0 <- coef(model.0)[2]
  start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
  
  # Re-fit the model with estimated starting parameters
  return(nls(y ~ alpha * exp(beta * x) + theta, start = start, data=train.df,
             control=nls.control(maxiter=1000, warnOnly=T)))
}

# Load a matrix of p-values
load.pval.matrix <- function(matrix.in, has.coords=T, ncols.coords=3, p.is.phred=T){
  x <- read.table(matrix.in, header=T, sep="\t", comment.char="")
  if(has.coords == T){
    colnames(x)[1:ncols.coords] <- c("chrom", "start", "end", "gene")[1:ncols.coords]
    x[, 1] <- as.character(x[, 1])
    x[, 2:3] <- apply(x[, 2:3], 2, as.numeric)
    x[, -c(1:ncols.coords)] <- apply(x[, -c(1:ncols.coords)], 2, as.numeric)
    coords <- as.data.frame(x[, 1:ncols.coords])
    pvals <- as.data.frame(x[, -c(1:ncols.coords)])
  }else{
    coords <- NULL
    pvals <- as.data.frame(apply(x, 2, as.numeric))
  }
  if(p.is.phred == T){
    pvals <- as.data.frame(apply(pvals, 2, function(x){10^-x}))
  }
  expected <- ppoints(nrow(pvals))
  lambdas <- apply(pvals, 2, function(obs){dchisq(median(obs, na.rm=T), df=1)/dchisq(median(expected), df=1)})
  return(list("coords" = coords, "pvals" = pvals,
              "expected" = expected, "lambdas" = lambdas))
}

# Load a BED-like tsv of gene features
load.features <- function(features.in, fill=NA, norm=F){
  feats <- read.table(features.in, header=T, sep="\t", comment.char="", check.names=F)[, -c(1:3)]
  feats[, -1] <- apply(feats[, -1], 2, as.numeric)
  if(!is.na(fill)){
    feats[, -1] <- apply(feats[, -1], 2, function(vals){
      na.idxs <- which(is.na(vals) | is.infinite(vals))
      if(length(na.idxs) > 0){
        if(fill=="mean"){
          vfill <- mean(vals[-na.idxs])
        }else if(fill=="median"){
          vfill <- median(vals[-na.idxs])
        }else if(is.numeric(fill)){
          vfill <- fill
        }
        vals[na.idxs] <- vfill
      }
      return(vals)
    })
  }
  if(norm==T){
    feats[, -1] <- apply(feats[, -1], 2, function(vals){
      scale(vals, scale=T, center=T)
    })
  }
  return(feats)
}

# Load table of gene feature metadata (including plain English labels)
load.gene.feature.metadata <- function(feature.metadata.in){
  meta <- read.table(feature.metadata.in, header=T, sep="\t", comment.char="")
  colnames(meta)[1] <- gsub("^X.", "", colnames(meta)[1])
  return(meta)
}

# Format a scientific value for printing to plots
format.scientific <- function(x, nsmall=2, max.decimal=3){
  parts <- unlist(strsplit(format(x, scientific=T), split="e", fixed=T))
  base <- format(round(as.numeric(parts[1]), max.decimal), nsmall=nsmall)
  exp <- as.character(as.numeric(parts[2]))
  bquote(.(base) ~ "x" ~ 10 ^ .(exp))
}

# Generate colors for XY scatterplot based on density
color.points.by.density <- function(x, y, palette=NULL){
  # Based on heatscatter.R from Colby Chiang
  # (https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R)
  plot.df <- data.frame("x"=x, "y"=y)
  plot.df <- plot.df[which(!is.infinite(plot.df$x) & !is.infinite(plot.df$y)
                           & !is.na(plot.df$x) & !is.na(plot.df$y)), ]
  dens <- densCols(plot.df$x, plot.df$y, colramp=colorRampPalette(c("black", "white")))
  plot.df$dens <- col2rgb(dens)[1, ] + 1L
  if(is.null(palette)){
    palette <- viridis(256)
  }
  plot.df$col <- palette[plot.df$dens]
  plot.df[order(plot.df$dens), ]
}
