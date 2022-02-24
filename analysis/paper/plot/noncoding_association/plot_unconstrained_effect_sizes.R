#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute & plot case:control effect sizes before and after loose noncoding filtering for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load table of sample sizes
load.samples.table <- function(pheno.table.in){
  samples <- read.table(pheno.table.in, header=T, sep="\t", comment.char="")
  colnames(samples)[1] <- "hpo"
  samples[, c(which(colnames(samples) == "hpo"),
            grep("meta", colnames(samples), fixed=T))]
}

# Load summary table of CNV counts per hpo by gene content
load.counts <- function(counts.in, samples){
  # Load data
  counts <- read.table(counts.in, header=T, sep="\t", comment.char="")
  colnames(counts) <- c("cohort", "hpo", "gset", "alt")

  # Add reference counts for each row
  counts$ref <- apply(counts, 1, function(vals){
    cohort <- vals[1]; hpo <- vals[2]; alt <- as.numeric(vals[4])
    max(c(0, samples[which(samples$hpo==hpo), cohort] - alt))
  })

  # Add fourth gene set reflecting unconstrained_only after removing
  # samples from "not_unconstrained" from denominator
  fourth.df <- as.data.frame(t(apply(counts[which(counts$gset == "unconstrained_only"), ],
                                     1, function(vals){
    cohort <- vals[1]; hpo <- vals[2]; alt <- vals[4]; ref <- as.numeric(vals[5])
    gset <- "unconstrained_only_denom_subbed"
    n_sub <- as.numeric(counts$alt[which(counts$cohort == cohort
                                         & counts$hpo == hpo
                                         & counts$gset == "not_unconstrained")])
    c(cohort, hpo, gset, alt, max(c(0, ref - n_sub)))
  })))
  colnames(fourth.df) <- colnames(counts)
  counts <- rbind(counts, fourth.df)
  counts[, c("alt", "ref")] <- apply(counts[, c("alt", "ref")], 2, as.numeric)


  # Reformat data to match expected inputs for meta-analysis functions
  hpos <- unique(counts$hpo[which(counts$hpo != "HEALTHY_CONTROL")])
  gsets <- unique(counts$gset)
  cohorts <- unique(counts$cohort)
  as.data.frame(do.call("rbind", lapply(hpos, function(hpo){
    sub.df <- as.data.frame(do.call("rbind", lapply(gsets, function(gset){
      vals <- unlist(lapply(cohorts, function(cohort){
        base.hits <- which(counts$gset == gset & counts$cohort == cohort)
        vals <- c(counts$ref[intersect(base.hits, which(counts$hpo == "HEALTHY_CONTROL"))],
                  counts$ref[intersect(base.hits, which(counts$hpo == hpo))],
                  counts$alt[intersect(base.hits, which(counts$hpo == "HEALTHY_CONTROL"))],
                  counts$alt[intersect(base.hits, which(counts$hpo == hpo))])
        names(vals) <- paste(c("control_ref", "case_ref", "control_alt", "case_alt"),
                             cohort, sep=".")
        return(vals)
      }))
    })))
    sub.df$gset <- gsets
    sub.df$hpo <- hpo
    sub.df$exclude_cohorts <- ""
    return(sub.df)
  })))
}

# Wrapper for meta-analysis of all hpos and gsets
meta.all <- function(counts, cohorts){
  hpos <- sort(unique(counts$hpo))
  hpos <- hpos[which(hpos != "HEALTHY_CONTROL")]
  gsets <- unique(counts$gset)
  meta.res <- as.data.frame(do.call("rbind", lapply(1:nrow(counts), function(i){
    c(meta.single(counts, cohorts, i, empirical.continuity=F),
      counts$hpo[i], counts$gset[i])
  })))
  meta.res[, 1:5] <- apply(meta.res[, 1:5], 2, as.numeric)
  colnames(meta.res) <- c("meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper",
                          "meta_z", "meta_neg_log10_p", "hpo", "gset")
  return(meta.res)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Dotplot of effect sizes based on genic context
lnor.dotplot <- function(lnors, hpo="HP:0000118",
                         parmar=c(0.25, 10, 2.25, 0.25)){
  # Get plot data
  plot.dat <- lapply(lnors, function(df){exp(df[which(df$hpo==hpo), 1:3])})
  xlims <- range(do.call("rbind", plot.dat), na.rm=T)

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=c(3, 0),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")

  # Add background shading and gridlines
  x.ax.at <- axTicks(1)
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, bty="n", col=bluewhite)
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=0:3-0.1, ytop=(0:3)+0.1,
       col="white", border=NA, bty="n")
  abline(v=x.ax.at, col="white")
  abline(v=1, col=blueblack)

  # Add points
  y.mods <- 0.1 * c(-1, 1)
  sapply(1:2, function(i){
    segments(x0=plot.dat[[i]][, 2], x1=plot.dat[[i]][, 3],
             y0=(1:3)-0.5+y.mods[i], y1=(1:3)-0.5+y.mods[i],
             lwd=1.5, lend="round", col=cnv.colors[i])
    points(x=plot.dat[[i]][, 1], y=(1:3)-0.5+y.mods[i],
           pch=19, col=cnv.colors[i])
  })

  # Add y-axis labels
  sapply(1:3, function(i){axis(2, at=c(i-0.9, i-0.1), tck=0, labels=NA, col=blueblack)})
  axis(2, at=(1:3)-0.5, tick=F, las=2, cex.axis=5/6, line=-0.8,
       labels=c("All rCNVs",
                "rCNVs overlapping exons\nfrom at least one\nnon-unconstrained gene",
                "Noncoding rCNVs +\nrCNVs overlapping exons\nfrom unconstrained genes"))

  # Add x-axis
  axis(3, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(3, at=x.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(3, at=x.ax.at, tick=F, line=-0.65, cex.axis=5/6)
  mtext(3, line=1.25, text=paste("Odds Ratio (", hpo.abbrevs[hpo], ")", sep=""))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(metafor, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog del.counts.tsv dup.counts.tsv",
                                            "pheno.table.tsv out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: del.counts.tsv, dup.counts.tsv,",
             "pheno.table.tsv, and out.prefix\n"))
}

# Writes args & opts to vars
del.counts.in <- args$args[1]
dup.counts.in <- args$args[2]
pheno.table.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# del.counts.in <- "~/scratch/unconstrained_cnv_counts.DEL.tsv.gz"
# dup.counts.in <- "~/scratch/unconstrained_cnv_counts.DUP.tsv.gz"
# pheno.table.in <- "~/scratch/HPOs_by_metacohort.w_DEV.table.tsv"
# out.prefix <- "~/scratch/unconstrained_lnORs_test"

# Load sample size table
samples <- load.samples.table(pheno.table.in)
cohorts <- colnames(samples)[which(colnames(samples) != "hpo")]

# Load CNV counts and infer number of reference samples
del.counts <- load.counts(del.counts.in, samples)
dup.counts <- load.counts(dup.counts.in, samples)

# Compute effect sizes per hpo and gset
lnors <- list("DEL"=meta.all(del.counts, cohorts),
              "DUP"=meta.all(dup.counts, cohorts))

# Dotplot of effect sizes by genic context
pdf(paste(out.prefix, "cnv_lnORs_by_genic_context.pdf", sep="."),
    height=2.5, width=4.5)
lnor.dotplot(lnors, hpo="HP:0012759")
dev.off()
