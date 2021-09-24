#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
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
  samples[, c("hpo", paste("meta", 1:4, sep=""))]
}

# Load summary table of CNV counts per hpo by gene content
load.counts <- function(counts.in, samples){
  counts <- read.table(counts.in, header=T, sep="\t", comment.char="")
  colnames(counts) <- c("cohort", "hpo", "gset", "alt")
  counts$ref <- apply(counts, 1, function(vals){
    cohort <- vals[1]; hpo <- vals[2]; alt <- as.numeric(vals[4])
    max(c(0, samples[which(samples$hpo==hpo), cohort] - alt))
  })
  counts[with(counts, order(cohort, hpo)), 
         c("cohort", "hpo", "gset", "ref", "alt")]
}

# Apply empirical continuity correction to meta-analysis data frame
# Per Sweeting et al., Stat. Med., 2004 (section 3.3)
sweeting.correction <- function(meta.df, cc.sum=0.01){
  # Count number of carriers & non-carriers
  n.alt <- sum(meta.df[, grep("_alt", colnames(meta.df), fixed=T)])
  n.ref <- sum(meta.df[, grep("_ref", colnames(meta.df), fixed=T)])
  # Require at least one CNV to be observed
  if(n.alt>0){
    nt <- n.alt
    R <- n.ref/n.alt
    # Pooled odds ratio estimate of all non-zero studies with at least one case sample
    nonzero.studies <- intersect(which(apply(meta.df[, grep("_alt", colnames(meta.df), fixed=T)], 1, sum)>0),
                                 which(apply(meta.df[, grep("case_", colnames(meta.df), fixed=T)], 1, sum)>0))
    nonzero.case.odds <- sum(meta.df$case_alt[nonzero.studies])/sum(meta.df$case_ref[nonzero.studies])
    nonzero.control.odds <- sum(meta.df$control_alt[nonzero.studies])/sum(meta.df$control_ref[nonzero.studies])
    if(!is.nan(nonzero.case.odds) & !is.nan(nonzero.control.odds)){
      if(nonzero.control.odds>0){
        ohat <- nonzero.case.odds/nonzero.control.odds
        # Otherwise, apply standard continuity correction of 0.5 to pooled estimate if no CNVs observed in controls
      }else{
        nonzero.case.odds <- (sum(meta.df$case_alt[nonzero.studies])+0.5)/(sum(meta.df$case_ref[nonzero.studies])+0.5)
        nonzero.control.odds <- (sum(meta.df$control_alt[nonzero.studies])+0.5)/(sum(meta.df$control_ref[nonzero.studies])+0.5)
        ohat <- nonzero.case.odds/nonzero.control.odds
      }
      # Otherwise, apply standard continuity correction to pooled estimate of *all* studies
    }else{
      nonzero.case.odds <- (sum(meta.df$case_alt)+0.5)/(sum(meta.df$case_ref)+0.5)
      nonzero.control.odds <- (sum(meta.df$control_alt)+0.5)/(sum(meta.df$control_ref)+0.5)
      ohat <- nonzero.case.odds/nonzero.control.odds
    }
    
    # Solve for kc & kt
    kc <- R/(R+ohat)
    kt <- ohat/(R+ohat)
    # Compute continuity corrections
    cor.case_alt <- cc.sum * kt
    cor.case_ref <- cc.sum * kc
    cor.control_alt <- cc.sum * (nt + kt)
    cor.control_ref <- cc.sum * ((nt*R) + kc)
    # Apply continuity corrections
    meta.df$case_alt <- meta.df$case_alt + cor.case_alt
    meta.df$case_ref <- meta.df$case_ref + cor.case_ref
    meta.df$control_alt <- meta.df$control_alt + cor.control_alt
    meta.df$control_ref <- meta.df$control_ref + cor.control_ref
  }
  return(meta.df)
}

# Make meta-analysis data frame for a single hpo & gset
make.meta.df <- function(counts, hpo, gset, empirical.continuity=T){
  meta.df <- data.frame("cohort"=1:4,
                        "control_ref"=as.numeric(counts[which(counts$hpo=="HEALTHY_CONTROL" & counts$gset==gset), "ref"]),
                        "case_ref"=as.numeric(counts[which(counts$hpo==hpo & counts$gset==gset), "ref"]),
                        "control_alt"=as.numeric(counts[which(counts$hpo=="HEALTHY_CONTROL" & counts$gset==gset), "alt"]),
                        "case_alt"=as.numeric(counts[which(counts$hpo==hpo & counts$gset==gset), "alt"]),
                        "cohort_name"=paste("meta", 1:4, sep=""))
  if(empirical.continuity==T){
    meta.df <- sweeting.correction(meta.df)
  }
  return(meta.df)
}

# Fixed effects meta-analysis for a single hpo & gset
meta.single <- function(counts, hpo, gset, empirical.continuity=T){
  meta.df <- make.meta.df(counts, hpo, gset, empirical.continuity)
  # Meta-analysis
  meta.res <- tryCatch(rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                               measure="OR", data=meta.df, method="FE", slab=cohort_name,
                               add=0, drop00=F, correct=F, digits=5, control=list(maxiter=100, stepadj=0.5)),
                       error=function(e){
                         print(paste("row", row.idx, "failed to converge. Retrying with more iterations...", sep=" "))
                         rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                 measure="OR", data=meta.df, method="FE", slab=cohort_name,
                                 add=0, drop00=F, correct=F, digits=5, control=list(maxiter=10000, stepadj=0.4))
                       })
  out.v <- as.numeric(c(meta.res$b[1,1], meta.res$ci.lb, meta.res$ci.ub,
                        meta.res$zval, -log10(meta.res$pval)))
  # Force to p-values reflecting Ha : OR > 1
  if(!is.na(out.v[1]) & !is.na(out.v[5])){
    if(out.v[1] < 0){
      out.v[5] <- 0
    }
  }
  return(out.v)
}

# Wrapper for meta-analysis of all hpos and gsets
meta.all <- function(counts){
  hpos <- sort(unique(counts$hpo))
  hpos <- hpos[which(hpos != "HEALTHY_CONTROL")]
  gsets <- unique(counts$gset)
  meta.res <- do.call("rbind", lapply(hpos, function(hpo){
    meta.res <- as.data.frame(do.call("rbind", lapply(gsets, meta.single, counts=counts, hpo=hpo)))
    meta.res$hpo <- hpo
    cbind(meta.res, gsets)
  }))
  colnames(meta.res) <- c("meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper", "meta_z", "meta_neg_log10_p", "hpo", "gset")
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
                "Noncoding rCNVs +\nrCNVs overlapping exons\nonly from unconstrained genes"))
  
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
require(funr, quietly=T)
require(metafor, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

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
# pheno.table.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/unconstrained_lnORs_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/noncoding_association/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load sample size table
samples <- load.samples.table(pheno.table.in)

# Load CNV counts and infer number of reference samples
del.counts <- load.counts(del.counts.in, samples)
dup.counts <- load.counts(dup.counts.in, samples)

# Compute effect sizes per hpo and gset
lnors <- list("DEL"=meta.all(del.counts), "DUP"=meta.all(dup.counts))

# Dotplot of effect sizes by genic context
pdf(paste(out.prefix, "cnv_lnORs_by_genic_context.pdf", sep="."),
    height=2.5, width=4.25)
lnor.dotplot(lnors)
dev.off()
