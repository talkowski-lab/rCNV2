#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot raw feature coefficients from gene scoring


options(stringsAsFactors=F, scipen=1000)
colors <- c("#D43925", "#2376B2", "#7459B2")


#################
### FUNCTIONS ###
#################
# Load gene score stats .tsv
load.stats <- function(path){
  x <- read.table(path, header=T, sep="\t", comment.char="")[, 1:3]
  colnames(x)[1] <- "gene"
  x[, 2:3] <- apply(x[, 2:3], 2, as.numeric)
  return(x)
}

# Load raw gene features & subset to genes with scores
load.features <- function(stats, features.in){
  f <- read.table(features.in, sep="\t", comment.char="", header=T)[, -c(1:3)]
  # Standard normalize all features
  f[, 2:ncol(f)] <- apply(f[, 2:ncol(f)], 2, function(vals){
    newvals <- as.numeric(scale(as.numeric(vals), center=T, scale=T))
    nas <- which(is.na(newvals))
    if(length(nas > 0)){
      newvals[nas] <- mean(newvals, na.rm=T)
    }
    return(newvals)
  })
  f <- f[which(f$gene %in% stats$gene), ]
  return(f)
}

# Get raw correlation of features vs PIPs
calc.feature.cors <- function(stats, features, score, logit=F){
  x.all <- merge(stats[, which(colnames(stats) %in% c("gene", score))],
                 features, by="gene", sort=F, all=F)
  rhos <- sapply(2:ncol(features), function(i){
    x <- data.frame("score"=x.all[, which(colnames(x.all) == score)],
                    "feature"=x.all[, which(colnames(x.all) == colnames(features)[i])])
    if(logit==T){
      glm(score ~ feature, family="binomial", data=x)$coefficients[-1]
    }else{
      as.numeric(cor.test(x$score, x$feature, method="spearman")$estimate)
    }
  })
  rhos.df <- data.frame("feature"=colnames(features)[2:ncol(features)],
                        "value"=as.numeric(rhos))
  return(rhos.df)
}


# Barplot of top and bottom N feature correlations
plot.feat.cors <- function(stats, features, score, logit=F, top.N=10, center.buffer=6){
  # Get plot data
  data <- calc.feature.cors(stats, features, score, logit)
  vals <- data$value
  labs <- data$feature
  labs <- labs[which(!is.na(vals))]
  vals <- vals[which(!is.na(vals))]
  new.order <- order(vals)
  vals <- vals[new.order]
  labs <- labs[new.order]
  if(score == "pHI"){
    color <- colors[1]
  }else{
    color <- colors[2]
  }
  
  # Set parameters
  n.pos.feats <- min(c(length(which(vals>0)), top.N))
  pos.idx <- tail(1:length(vals), n.pos.feats)
  n.neg.feats <- min(c(length(which(vals<0)), top.N))
  neg.idx <- head(1:length(vals), n.neg.feats)
  n.feats <- n.pos.feats + n.neg.feats
  valmax <- max(abs(vals), na.rm=T)
  
  par(mfrow=c(1, 2), bty="n")
  
  # Left panel — negative features
  par(mar=c(0.5, 1, 4, center.buffer))
  plot(x=c(-valmax, 0), y=c(0, n.feats), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  abline(v=axTicks(3), col="gray90")
  if(n.neg.feats > 0){
    rect(xleft=vals[neg.idx], xright=0,
         ybottom=(1:n.neg.feats)-0.85, ytop=(1:n.neg.feats)-0.15,
         col=color)
    axis(4, at=(1:n.neg.feats)-0.5, las=2, line=-0.9, tick=F,
         labels=labs[neg.idx])
  }
  axis(3)
  if(logit==T){
    mtext(3, text=expression("Standardized logit" ~ beta), line=2.5)
  }else{
    mtext(3, text="Spearman's Rho", line=2.5)
  }
  
  # Right panel — positive features
  par(mar=c(0.5, center.buffer, 4, 1))
  plot(x=c(0, valmax), y=c(0, n.feats), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  abline(v=axTicks(3), col="gray90")
  if(n.pos.feats > 0){
    rect(xleft=vals[pos.idx], xright=0,
         ybottom=(1:n.pos.feats)-0.85+n.neg.feats,
         ytop=(1:n.pos.feats)-0.15+n.neg.feats,
         col=color)
    axis(2, at=(1:n.pos.feats)-0.5+n.neg.feats, las=2, line=-0.9, tick=F,
         labels=labs[pos.idx])
  }
  axis(3)
  if(logit==T){
    mtext(3, text=expression("Standardized logit" ~ beta), line=2.5)
  }else{
    mtext(3, text="Spearman's Rho", line=2.5)
  }
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog scores.tsv features.bed.gz out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Three positional arguments: scores.tsv, features.bed.gz and out_prefix\n")
}

# Writes args & opts to vars
stats.in <- args$args[1]
features.in <- args$args[2]
out.prefix <- args$args[3]

# # DEV PARAMTERS
# setwd("~/scratch/")
# stats.in <- "rCNV.gene_scores.tsv.gz"
# features.in <- "gencode.v19.canonical.pext_filtered.all_features.bed.gz"
# out.prefix <- "gene_scores_vs_features"

# Load gene scores & features
stats <- load.stats(stats.in)
features <- load.features(stats, features.in)

# Plot correlations with raw features
for(score in c("pHI", "pTS")){
  pdf(paste(out.prefix, score, "rhos.pdf", sep="."),
      height=6, width=6)
  plot.feat.cors(stats, features, score=score, logit=F)
  dev.off()
}

# Plot logit betas
for(score in c("pHI", "pTS")){
  pdf(paste(out.prefix, score, "betas.pdf", sep="."),
      height=6, width=6)
  plot.feat.cors(stats, features, score=score, logit=T)
  dev.off()
}
