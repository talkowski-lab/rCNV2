#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Burden meta-analysis of noncoding rCNV counts across annotation tracks


options(scipen=10000, stringsAsFactors=F)


#################
### FUNCTIONS ###
#################
# Load track stats 
load.stats <- function(stats.in, cnv.split=T){
  stats <- read.table(stats.in, sep="\t", comment.char="", 
                      check.names=F, header=T)
  colnames(stats)[1] <- gsub("#", "", colnames(stats)[1], fixed=T)
  if(cnv.split==T){
    list("DEL" = stats[which(stats$cnv=="DEL"), ],
         "DUP" = stats[which(stats$cnv=="DUP"), ])
  }else{
    return(stats)
  }
}

# Extract cohort names from an input stats file
extract.cohorts <- function(stats){
  as.vector(sapply(colnames(stats)[grep("_case_alt", colnames(stats), fixed=T)],
                   function(str){
                     unlist(strsplit(str, split="_"))[1]
                   }))
}

# Apply saddlepoint approximation to vector of Z-scores to generate adjusted P-values and Z-scores
saddlepoint.adj <- function(zscores, phred=T, min.p=10e-300){
  mu.hat <- mean(zscores, na.rm=T)
  sd.hat <- sd(zscores, na.rm=T)
  zscores <- zscores - mu.hat
  cumuls <- gaussianCumulants(0, sd.hat)
  dx <- 0.01
  x <- seq(min(min(zscores, na.rm=T), -40), 
           max(max(zscores, na.rm=T), 40), dx)
  saddle.pdf <- saddlepoint(x, 1, cumuls)$approx
  saddle.cdf <- caTools::cumsumexact(saddle.pdf * dx)
  calc.saddle.p.z <- function(z){
    if(!is.na(z)){
      # negative absolute value trick to account for R floating point precision producing p=1 too early
      if(z>=0){
        s <- -1
      }else{
        s <- 1
      }
      mirror.z <- -abs(z)
      p <- max(c(tail(saddle.cdf[which(x<mirror.z)], 1), min.p), na.rm=T)
      new.z <- s * qnorm(p)
    }else{
      p <- NA
      new.z <- NA
    }
    return(c(new.z, p))
  }
  new.stats <- t(sapply(zscores, calc.saddle.p.z))
  colnames(new.stats) <- c("zscore", "p")
  if(phred==T){
    new.stats[, 2] <- -log10(new.stats[, 2])
    colnames(new.stats)[2] <- "phred_p"
  }
  return(as.data.frame(new.stats))
}

# Calculate odds ratios, variance, and Z-scores per cohort (with SPA, if optioned)
calc.ors <- function(stats, cohorts, spa=T){
  or.df <- as.data.frame(do.call("cbind", lapply(cohorts, function(cohort){
    ors <- as.data.frame(t(sapply(1:nrow(stats), function(i){
      as.numeric(metafor::escalc("OR", 
                      ai=stats[i, paste(cohort, "control", "ref", sep="_")],
                      bi=stats[i, paste(cohort, "case", "ref", sep="_")],
                      ci=stats[i, paste(cohort, "control", "alt", sep="_")],
                      di=stats[i, paste(cohort, "case", "alt", sep="_")]))
    })))
    colnames(ors) <- c("lnOR", "var")
    # Z-score computed against per-category variance and empirical mean across all categories
    # Assumes the average category should _not_ have an effect, so any systematic deviance from
    # mean = 0 is due to filtering biases
    ors$zscore <- (ors$lnOR - mean(ors$lnOR, na.rm=T)) / sqrt(ors$var)
    # Adjust Z-score with SPA, if optioned
    if(spa==T){
      adj.z.p <- saddlepoint.adj(ors$zscore, phred=T)
      ors$zscore <- adj.z.p$zscore
      ors$phred_p <- adj.z.p$phred_p
    }
    colnames(ors) <- paste(cohort, colnames(ors), sep=".")
    return(ors)
  })))
  return(as.data.frame(cbind(stats, or.df)))
}

# Meta analysis of Z-scores weighted by inverse variance
weighted.z <- function(stats, cohorts, spa=T){
  meta.z <- sapply(1:nrow(stats), function(i){
    z <- stats[i, grep("zscore", colnames(stats), fixed=T)]
    var <- stats[i, grep("var", colnames(stats), fixed=T)]
    inv.var <- 1/var
    weighted.mean(z, inv.var, na.rm=T)
  })
  if(spa==T){
    meta.z.p <- saddlepoint.adj(meta.z, phred=T)
  }else{
    meta.z.p <- data.frame("zscore" = meta.z,
                           "phred_p" = -log10(pnorm(abs(meta.z), lower.tail=F)))
  }
  meta.z.p$lnOR <-  sapply(1:nrow(stats), function(i){
    lnor <- stats[i, grep("lnOR", colnames(stats), fixed=T)]
    var <- stats[i, grep("var", colnames(stats), fixed=T)]
    inv.var <- 1/var
    weighted.mean(lnor, inv.var, na.rm=T)
  })
  meta.z.p$phred_fdr_q <- -log10(p.adjust(10^-meta.z.p$phred_p, method="fdr"))
  colnames(meta.z.p) <- paste("meta", colnames(meta.z.p), sep=".")
  as.data.frame(cbind(stats, meta.z.p))
}

# Volcano plot
volcano <- function(meta.res, p.cutoff=-log10(0.05), color="red", ymax=NULL){
  if(is.null(ymax)){
    ymax <- max(p.cutoff, max(meta.res$meta.phred_p, na.rm=T))
  }
  ylims <- c(0, ymax)
  xlims <- range(meta.res$meta.lnOR)
  sig <- which(meta.res$meta.zscore>0 & meta.res$meta.phred_p>=p.cutoff)
  par(mar=c(3.5, 3.5, 1.3, 1))
  plot(NA, xlim=xlims, ylim=ylims, xlab="", ylab="")
  abline(v=0)
  segments(x0=0, x1=par("usr")[2], y0=p.cutoff, y1=p.cutoff, lty=2, col=color)
  points(meta.res$meta.lnOR[-sig], meta.res$meta.phred_p[-sig], cex=0.2, col="gray30")
  points(meta.res$meta.lnOR[sig], meta.res$meta.phred_p[sig], cex=0.2, col=color)
  mtext(1, line=2.25, text=bquote(italic("ln") * (OR)))
  mtext(2, line=2, text=bquote(-log[10](italic(P))))
}


################
###RSCRIPT BLOCK
################
# Load required libraries
require(optparse, quietly=T)
require(metafor, quietly=T)
require(EQL, quietly=T)
require(caTools, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--cutoff"), type="numeric", default=0.05, 
              help="meta-analysis P-value (or FDR q-value) cutoff to consider a track significant [default '%default']",
              metavar="numeric"),
  make_option(c("--use-fdr"), action="store_true", default=FALSE,
              help="use FDR q-value for determining significance (rather than P-value) [default %default]"),
  make_option(c("--signif-tracks"), type="character",  
              help="output file for significant tracks", metavar="string"),
  make_option(c("--volcano"), type="character",  
              help="path to (optional) volcano plot .png", metavar="string")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog infile outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
stats.in <- args$args[1]
outfile <- args$args[2]
cutoff <- -log10(as.numeric(opts$cutoff))
use.fdr <- opts$`use-fdr`
signif.outfile <- opts$`signif-tracks`
volcano.out <- opts$`volcano`

# # Dev parameters
# stats.in <- "~/scratch/rCNV.all.merged_stats.with_counts.tsv.gz"
# outfile <- "~/scratch/rCNV.chromhmm_plus_encode_plus_enhancers.test.tsv"
# cutoff <- -log10(0.05)
# use.fdr <- F
# signif.outfile <- "~/scratch/rCNV.rCNV.chromhmm_plus_encode_plus_enhancers.test.signif_tracks.list"
# volcano.out <- "~/scratch/test.volcano.png"

# Read track stats and split by CNV type
stats <- load.stats(stats.in, cnv.split=T)
cohorts <- extract.cohorts(stats[[1]])

# Calculate case:control ORs, variance, and Z-scores for each track
stats <- lapply(stats, calc.ors, cohorts, spa=T)

# Weighted Z-score meta-analysis
meta.res.split <- lapply(stats, weighted.z, cohorts)
meta.res <- as.data.frame(do.call("rbind", meta.res.split))

# Extract significant track names
if(!is.null(signif.outfile)){
  if(use.fdr==T){
    sig.idx <- which(meta.res$meta.phred_fdr_q >= cutoff & meta.res$meta.zscore > 0)
  }else{
    sig.idx <- which(meta.res$meta.phred_p >= cutoff & meta.res$meta.zscore > 0)
  }
  if(length(sig.idx) > 0){
    sig.tracks <- unique(meta.res[sig.idx, which(colnames(meta.res) %in% c("trackname", "original_path"))])
    write.table(sig.tracks, signif.outfile,
                col.names=F, row.names=F, quote=F, sep="\t")
  }
}

# Write stats to outfile
colnames(meta.res)[1] <- paste("#", colnames(meta.res)[1], sep="")
write.table(meta.res, outfile, sep="\t",
            row.names=F, col.names=T, quote=F)

# Volcano plots, if optioned
if(!is.null(volcano.out)){
  png(volcano.out, height=4*300, width=8*300, res=300)
  par(mfrow=c(1, 2))
  volcano(meta.res.split$DEL, p.cutoff=cutoff, color="red", ymax=max(meta.res$meta.phred_p, na.rm=T))
  mtext(3, line=0.1, text="Deletions", col="red", font=2)
  volcano(meta.res.split$DUP, p.cutoff=cutoff, color="blue", ymax=max(meta.res$meta.phred_p, na.rm=T))
  mtext(3, line=0.1, text="Duplications", col="blue", font=2)
  dev.off()
}
