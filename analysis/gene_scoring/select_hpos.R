#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Select phenotypes for gene scoring


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############

# Load a single tsv of counts per HPO
load.data.single <- function(path){
  x <- read.table(path, sep="\t", comment.char="", header=T)
  colnames(x)[1] <- "hpo"
  res <- as.data.frame(t(sapply(unique(x$hpo), function(hpo){
    c(hpo,
      as.numeric(c(x$hit[which(x$hpo==hpo & x$gset=="pos")],
                   x$hit[which(x$hpo==hpo & x$gset=="neg")])))
  })))
  colnames(res) <- c("hpo", "pos", "neg")
  return(res)
}


# Load input tsvs for all cohorts
load.data <- function(infile){
  x <- read.table(infile, sep="\t", header=F)
  dlist <- lapply(1:nrow(x), function(i){
    load.data.single(x[i, 2])
  })
  names(dlist) <- as.character(x[, 1])
  return(dlist)
}


# Make meta-analysis data frame for a single gene
make.meta.df <- function(dlist, cohorts, hpo){
  ncohorts <- length(cohorts)
  meta.df <- data.frame("cohort"=1:ncohorts,
                        "control_neg"=as.numeric(sapply(dlist, function(df){df$neg[which(df$hpo=="HEALTHY_CONTROL")]})),
                        "case_neg"=as.numeric(sapply(dlist, function(df){df$neg[which(df$hpo==hpo)]})),
                        "control_pos"=as.numeric(sapply(dlist, function(df){df$pos[which(df$hpo=="HEALTHY_CONTROL")]})),
                        "case_pos"=as.numeric(sapply(dlist, function(df){df$pos[which(df$hpo==hpo)]})),
                        "cohort_name"=cohorts)
  return(meta.df)
}


# Perform meta-analysis for a single gene
meta.single <- function(dlist, cohorts, hpo){
  meta.df <- make.meta.df(dlist, cohorts, hpo)
  # If strictly zero case CNVs are observed, unable to estimate effect size
  if(all(meta.df$case_pos==0)){
    out.v <- c(rep(NA, 4), 0)
  }else{
    # Meta-analysis
    meta.res <- tryCatch(rma.uni(ai=control_neg, bi=case_neg, ci=control_pos, di=case_pos,
                                 measure="OR", data=meta.df, method="FE", slab=cohort_name,
                                 add=0, drop00=F, correct=F, digits=5, control=list(maxiter=100, stepadj=0.5)),
                         error=function(e){
                           rma.uni(ai=control_neg, bi=case_neg, ci=control_pos, di=case_pos,
                                   measure="OR", data=meta.df, method="FE", slab=cohort_name,
                                   add=0, drop00=F, correct=F, digits=5, control=list(maxiter=10000, stepadj=0.4))
                         })
    out.v <- as.numeric(c(meta.res$b[1,1], meta.res$ci.lb, meta.res$ci.ub,
                          meta.res$zval, -log10(meta.res$pval)))
    # Force to p-values neglecting Ha : OR > 1
    if(!is.na(out.v[1]) & !is.na(out.v[5])){
      if(out.v[1] < 0){
        out.v[5] <- 0
      }
    }
  }
  return(out.v)
}


# Wrapper to run meta-analysis for all hpos
meta <- function(dlist, cohorts){
  case.hpos <- setdiff(unique(dlist[[1]]$hpo), "HEALTHY_CONTROL")
  meta.res <- as.data.frame(cbind(case.hpos, t(sapply(case.hpos, meta.single, dlist=dlist, cohorts=cohorts))))
  meta.res[, -1] <- apply(meta.res[, -1], 2, as.numeric)
  meta.res <- as.data.frame(meta.res)
  colnames(meta.res) <- c("hpo", "lnOR", "lnOR.lower", "lnOR.upper", "zscore", "neglog10.pval")
  return(meta.res)
}


# Plot effect sizes for all phenotypes
plot.effects <- function(del.stats, dup.stats){
  # Get plot values
  hpos.all <- del.stats$hpo
  hpos.ordered <- hpos.all[order(sapply(hpos.all, function(hpo){min(c(del.stats$lnOR[which(del.stats$hpo==hpo)],
                                                                      dup.stats$lnOR[which(dup.stats$hpo==hpo)]))}))]
  xlims <- range(rbind(del.stats[, 2:4], dup.stats[, 2:4]))
  
  # Prep plot area
  par(mar=c(1, 6, 3, 1), bty="n")
  plot(NA, xlim=xlims, ylim=c(0, length(hpos.all)),
       xaxt="n", yaxt="n", xlab="", ylab="")
  
  # Add axis
  x.at <- log(c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16))
  x.labs <- c("1/16", "1/8", "1/4", "1/2", 1, 2, 4, 8, 16)
  axis(3, at=x.at, labels=NA)
  axis(3, at=x.at, line=-0.6, tick=F, labels=x.labs)
  abline(v=0)
  mtext(3, line=1.5, text="Odds Ratio")
  
  # Add data
  sapply(1:length(hpos.ordered), function(i){
    hpo <- hpos.ordered[i]
    del.vals <- as.numeric(del.stats[which(del.stats$hpo==hpo), -1])
    del.sig <- del.vals[length(del.vals)] > -log10(0.05)
    if(del.sig==T){
      del.fill <- "red"
    }else{
      del.fill <- "white"
    }
    dup.vals <- as.numeric(dup.stats[which(dup.stats$hpo==hpo), -1])
    dup.sig <- dup.vals[length(dup.vals)] >= -log10(0.05)
    if(dup.sig==T){
      dup.fill <- "blue"
    }else{
      dup.fill <- "white"
    }
    vals <- rbind(del.vals, dup.vals)
    y.at <- i-0.5+c(0.15, -0.15)
    # if(del.sig & dup.sig){
    #   rect(xleft=par("usr")[1], xright=par("usr")[2],
    #        ybottom=i-0.9, ytop=i-0.1,
    #        border=NA, bty="n", col="yellow")
    # }
    segments(x0=vals[, 2], x1=vals[, 3], y0=y.at, y1=y.at,
             col=c("red", "blue"))
    points(x=vals[, 1], y=y.at, pch=22, col=c("red", "blue"), bg=c(del.fill, dup.fill))
    if(del.sig & dup.sig & hpo != "HP:0000118"){
      font <- 2
      text.color <- "black"
    }else{
      font <- 1
      text.color <- "gray70"
    }
    axis(2, at=i-0.5, line=-1, tick=F, labels=hpo, las=2, font=font, col.axis=text.color)
  })
  
  legend("bottomright", legend=c("DEL", "DUP", "Sig.", "NS"),
         col=c("red", "blue", "black", "black"), 
         pt.bg=c("red", "blue", "black", "white"),
         pch=22)
}


################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)
require(metafor, quietly=T)

# List of command-line options
option_list <- list(
  # make_option(c("--max-eval"), default=80,
  #             help="Maximum number of genes per CNVs to evaluate")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog del.in dup.in hpos.keep summary.pdf",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
del.in <- args$args[1]
dup.in <- args$args[2]
hpos.keep.out <- args$args[3]
out.pdf <- args$args[4]
max.eval <- opts$`max-eval`

# # Dev parameters
# del.in <- "~/scratch/select_hpos.DEL.input.tsv"
# dup.in <- "~/scratch/select_hpos.DUP.input.tsv"
# hpos.keep.out <- "~/scratch/keep_hpos.list"
# out.pdf <- "~/scratch/test.pdf"
# setwd("~/scratch/")

# Load data
del <- load.data(del.in)
dup <- load.data(dup.in)
cohorts <- names(del)

# Run meta-analyses
del.stats <- meta(del, cohorts)
dup.stats <- meta(dup, cohorts)

# Get HPOs to keep and write to list
keep.hpos <- setdiff(intersect(del.stats$hpo[which(del.stats$neglog10.pval >= -log10(0.05))],
                               dup.stats$hpo[which(dup.stats$neglog10.pval >= -log10(0.05))]),
                     c("HP:0000118"))
write.table(keep.hpos, hpos.keep.out, col.names=F, row.names=F, quote=F)

# Plot results
pdf(out.pdf, height=6, width=8)
plot.effects(del.stats, dup.stats)
dev.off()
