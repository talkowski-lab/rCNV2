#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot raw feature coefficients from fine-mapping


options(scipen=100000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Load a single finemap stats .tsv
load.stats.single <- function(path){
  x <- read.table(path, header=T, sep="\t", comment.char="")
  colnames(x)[1] <- "HPO"
  x[, 3:4] <- apply(x[, 3:4], 2, as.numeric)
  return(x)
}

# Load raw gene features & subset to genes with PIPs
load.features.single <- function(stats, features.in){
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
calc.feature.cors <- function(stats, features, logit=F){
  rhos <- sapply(2:ncol(features), function(i){
    x <- merge(features[, c(1, i)],
               stats[, which(colnames(stats) %in% c("gene", "PIP"))],
               all.x=F, all.y=T)
    x <- as.data.frame(t(sapply(sort(unique(x$gene)), function(g){
      apply(x[which(x$gene==g), -c(1)], 2, mean, na.rm=T)
    })))
    colnames(x) <- c("feature", "PIP")
    if(logit==T){
      glm(PIP ~ feature, family="binomial", data=x)$coefficients[-1]
    }else{
      as.numeric(cor.test(x$PIP, x$feature, method="spearman")$estimate)
    }
  })
  rhos.df <- data.frame("feature"=colnames(features)[2:ncol(features)],
                        "value"=as.numeric(rhos))
  return(rhos.df)
}

# Load a list of finemap stats and calculate feature coefficients
load.stats <- function(statslist.in){
  statslist <- read.table(statslist.in, header=F, sep="\t")
  colnames(statslist) <- c("name", "color", "path", "features")
  stats <- lapply(1:nrow(statslist), function(i){
    stats <- load.stats.single(statslist$path[i])
    features <- load.features.single(stats, statslist$features[i])
    rhos <- calc.feature.cors(stats, features, logit=F)
    betas <- calc.feature.cors(stats, features, logit=T)
    list("stats"=stats,
         "features"=features,
         "rhos"=rhos,
         "betas"=betas,
         "color"=statslist$color[i],
         "name"=statslist$name[i])
  })
  names(stats) <- statslist$name
  return(stats)
}

# Barplot of top and bottom N feature correlations
plot.feat.cors <- function(data, top.N=10, logit=F, center.buffer=6){
  # Get plot data
  if(logit==T){
    vals <- data$betas$value
    labs <- data$betas$feature
  }else{
    vals <- data$rhos$value
    labs <- data$rhos$feature
  }
  labs <- labs[which(!is.na(vals))]
  vals <- vals[which(!is.na(vals))]
  new.order <- order(vals)
  vals <- vals[new.order]
  labs <- labs[new.order]
  
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
         col=data$color)
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
         col=data$color)
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
args <- parse_args(OptionParser(usage="%prog finemap_res.tsv features.tsv out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Three positional arguments: inputs.tsv and out_prefix\n")
}

# Writes args & opts to vars
statslist.in <- args$args[1]
out.prefix <- args$args[2]

# # DEV PARAMTERS
# setwd("~/scratch/")
# statslist.in <- "finemap_feature_cor_input.tsv"
# out.prefix <- "finemap_features"

# Load stats & features, and calculate rhos/betas
stats <- load.stats(statslist.in)

# Plot correlations with raw features
sapply(1:length(stats), function(i){
  sanitized.name <- gsub(" ", "_", tolower(names(stats)[i]), fixed=T)
  pdf(paste(out.prefix, sanitized.name, "rhos.pdf", sep="."),
      height=6, width=6)
  plot.feat.cors(stats[[i]], logit=F)
  dev.off()
})

# Plot logit betas
sapply(1:length(stats), function(i){
  sanitized.name <- gsub(" ", "_", tolower(names(stats)[i]), fixed=T)
  pdf(paste(out.prefix, sanitized.name, "betas.pdf", sep="."),
      height=6, width=6)
  plot.feat.cors(stats[[i]], logit=T)
  dev.off()
})
