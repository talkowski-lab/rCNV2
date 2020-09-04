#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compare de novo CNVs from ASC/SPARK exomes vs. gene-level intolerance scores


options(stringsAsFactors=F, scipen=1000)


# Constants
case.phenos <- c("ASD", "ID", "ASD_ID")
control.pheno <- "Control"
all.score.names <- c("pHI" = "pHI",
                     "pTS" = "pTS",
                     "gnomad_pLI" = "pLI (Karczewski 2020)",
                     "gnomad_oe_mis_upper" = "Mis. OEUF (Karczewski 2020)",
                     "gnomad_oe_lof_upper" = "LOEUF (Karczewski 2020)",
                     "exac_cnv_z" = "CNV Z-score (Ruderfer 2016)",
                     "rvis_pct" = "RVIS (Petrovski 2013)",
                     "eds" = "EDS (Wang 2020)",
                     "hurles_hi" = "HI index (Huang 2010)")


######################
### DATA FUNCTIONS ###
######################
# Load ASC/SPARK phenotype metadata
load.phenos <- function(asd_phenos.in){
  phenos <- read.table(asd_phenos.in, header=T, sep="\t", comment.char="")
  colnames(phenos)[1] <- gsub("X.", "", colnames(phenos[1]))
  return(phenos)
}

# Load ASC/SPARK de novo CNVs and subset to samples present in phenotype metadata
load.cnvs <- function(asd_cnvs.in, phenos){
  cnvs <- read.table(asd_cnvs.in, header=T, sep="\t", comment.char="")[, -c(1:3)]
  keep.idxs <- which(cnvs$child_id %in% phenos$child_id)
  cnvs[keep.idxs, ]
}

# Annotate de novo CNVs with a single score
anno.cnv.single.score <- function(cnvs, score.data, score, operation="max"){
  as.numeric(as.vector(sapply(cnvs$genes, function(gstr){
    genes <- unlist(strsplit(gstr, split=";", fixed=T))
    keep.idxs <- which(score.data$gene %in% genes)
    if(length(keep.idxs) > 0){
      if(operation=="max"){
        max(score.data[keep.idxs, score], na.rm=T)
      }
    }else{
      NA
    }
  })))
}

# Annotate de novo CNVs with all scores
anno.cnv.scores <- function(cnvs, scores.rCNV, scores.other){
  cnvs$pHI <- anno.cnv.single.score(cnvs, scores.rCNV, "pHI")
  cnvs$pTS <- anno.cnv.single.score(cnvs, scores.rCNV, "pTS")
  cnvs.others <- sapply(colnames(scores.other)[-1], function(score){
    anno.cnv.single.score(cnvs, scores.other, score)
  })
  cnvs <- cbind(cnvs, cnvs.others)
  return(cnvs)
}

# Compute odds ratio for a given subset of CNVs
calc.or <- function(cnvs, phenos){
  n.case.all <- length(unique(phenos$child_id[which(phenos$phenotype %in% case.phenos)]))
  n.ctrl.all <- length(unique(phenos$child_id[which(phenos$phenotype == control.pheno)]))
  n.case.cnv <- length(unique(cnvs$child_id[which(cnvs$pheno %in% case.phenos)]))
  n.ctrl.cnv <- length(unique(cnvs$child_id[which(cnvs$pheno == control.pheno)]))
  n.case.ref <- n.case.all - n.case.cnv
  n.ctrl.ref <- n.ctrl.all - n.ctrl.cnv
  or.table <- matrix(c(n.ctrl.ref, n.case.ref, n.ctrl.cnv, n.case.cnv), nrow=2, byrow=T) + 0.5
  as.numeric(oddsratio.wald(or.table, correction=F)$measure[2, ])
}

# Compute odds ratios for CNVs binned by score (either percentile or absolute bin)
calc.or.by.scorebin <- function(cnvs, phenos, score, n.bins=3, cnv.pct=TRUE){
  bin.breaks <- seq(0, 1, length.out=n.bins+1)
  cnvs <- cnvs[which(!is.na(cnvs[, score])), ]
  if(cnv.pct==T){
    cnv.pct <- rank(cnvs[, score]) / nrow(cnvs)
  }else{
    cnv.pct <- cnvs[, score]
  }
  ors <- as.data.frame(t(sapply(1:n.bins, function(i){
    bin.cnvs <- cnvs[which(cnv.pct >= bin.breaks[i] & cnv.pct <= bin.breaks[i+1]), ]
    calc.or(bin.cnvs, phenos)
  })))
  colnames(ors) <- c("estimate", "lower", "upper")
  return(ors)
}

# Calculate odds ratios for dnCNVs stratified by high/low pHI & pTS
calc.or.stratified <- function(cnvs, phenos, high.cutoff=0.8, low.cutoff=0.5, 
                               log.trans=TRUE, norm.vs.baseline=TRUE){
  strat.ors <- lapply(c("DEL", "DUP"), function(cnvtype){
    ors <- t(data.frame("low.low"=calc.or(cnvs[which(cnvs$pHI<=low.cutoff & cnvs$pTS<=low.cutoff & cnvs$cnv==cnvtype), ], phenos),
    "high.low"=calc.or(cnvs[which(cnvs$pHI>=high.cutoff & cnvs$pTS<=low.cutoff & cnvs$cnv==cnvtype), ], phenos),
    "low.high"=calc.or(cnvs[which(cnvs$pHI<=low.cutoff & cnvs$pTS>=high.cutoff & cnvs$cnv==cnvtype), ], phenos),
    "high.high"=calc.or(cnvs[which(cnvs$pHI>=high.cutoff & cnvs$pTS>=high.cutoff & cnvs$cnv==cnvtype), ], phenos)))
    colnames(ors) <- c("estimate", "lower", "upper")
    if(norm.vs.baseline==TRUE){
      baseline <- calc.or(cnvs[which(cnvs$cnv==cnvtype), ], phenos)[1]
      ors <- ors / baseline
    }
    if(log.trans==TRUE){
      ors <- log2(ors)
    }
    return(ors)
  })
  names(strat.ors) <- c("DEL", "DUP")
  return(strat.ors)
}

# Compute ROC of a single score vs. proband/sibling labels
roc.asd_cnv <- function(cnvs, phenos, score, steps=seq(1, 0, -0.001)){
  # Make table of all samples with CNVs with highest score per sample
  carriers <- unique(cnvs$child_id)
  x <- data.frame("score"=c(), "true"=c(), "false"=c())
  for(id in carriers){
    child_pheno <- phenos$phenotype[which(phenos$child_id==id)]
    best.score <- max(cnvs[which(cnvs$child_id==id), score], na.rm=T)
    if(length(best.score) > 0){
      if(child_pheno == control.pheno){
        x <- rbind(x, data.frame("score"=best.score, "true"=FALSE, "false"=TRUE))
      }else if(child_pheno %in% case.phenos){
        x <- rbind(x, data.frame("score"=best.score, "true"=TRUE, "false"=FALSE))
      }
    }
  }
  roc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score >= k)
    ftrue <- length(which(x$true[idxs])) / length(which(x$true))
    ffalse <- length(which(x$false[idxs])) / length(which(x$false))
    fother <- length(which(!x$true[idxs])) / length(which(!x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, fother, ftrue, ffalse))
  })))
  roc_res <- rbind(c(-Inf, 0, 0, 0, 0),
                   roc_res,
                   c(Inf, 1, 1, 1, 1))
  colnames(roc_res) <- c("min_score", "frac_all", "frac_other", "frac_true", "frac_false")
  return(roc_res)
}

# Compute PRC of a single score vs. proband/sibling labels
prc.asd_cnv <- function(cnvs, phenos, score, steps=seq(1, 0, -0.001)){
  # Make table of all samples with CNVs with highest score per sample
  carriers <- unique(cnvs$child_id)
  x <- data.frame("score"=c(), "true"=c(), "false"=c())
  for(id in carriers){
    child_pheno <- phenos$phenotype[which(phenos$child_id==id)]
    best.score <- max(cnvs[which(cnvs$child_id==id), score], na.rm=T)
    if(length(best.score) > 0){
      if(child_pheno == control.pheno){
        x <- rbind(x, data.frame("score"=best.score, "true"=FALSE, "false"=TRUE))
      }else if(child_pheno %in% case.phenos){
        x <- rbind(x, data.frame("score"=best.score, "true"=TRUE, "false"=FALSE))
      }
    }
  }
  prc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score >= k)
    prec <- length(which(x$true[idxs])) / (length(which(x$true[idxs])) + length(which(x$false[idxs])))
    recall <- length(which(x$true[idxs])) / length(which(x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, prec, recall))
  })))
  prc_res <- rbind(c(-Inf, 0, 1, 0),
                   prc_res)
  colnames(prc_res) <- c("min_score", "frac_all", "precision", "recall")
  return(prc_res)
}

# Wrapper to calculate all performance stats a single score vs. proband/sibling labels
evaluate.score.asd_cnv <- function(cnvs, phenos, score){
  roc.res <- roc.asd_cnv(cnvs, phenos, score)
  roc.auc <- flux::auc(roc.res$frac_false, roc.res$frac_true)
  prc.res <- prc.asd_cnv(cnvs, phenos, score)
  prc.auc <- flux::auc(prc.res$recall, prc.res$precision)
  return(list("roc"=roc.res,
              "roc.auc"=roc.auc,
              "prc"=prc.res,
              "prc.auc"=prc.auc))
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot odds ratios for de novo CNVs by CNV score (either percentile or absolute bin)
plot.or.by.scorebin <- function(cnvs, phenos, score, n.bins=3, cnv.pct=TRUE,
                                x.label=NULL, parse.x.label=T, x.label.line=2,
                                x.ax.labels=NULL, parse.x.ax.labels=T,
                                cex.x.ax.labels=1, pt.color=blueblack, 
                                baseline.color=blueblack, null.color=bluewhite, 
                                parmar=c(3.5, 3.5, 0.5, 0.5)){
  # Gather plot data
  baseline <- log2(calc.or(cnvs, phenos)[1])
  ors <- log2(calc.or.by.scorebin(cnvs, phenos, score, n.bins, cnv.pct))
  or.vals <- as.numeric(unlist(ors))
  ylims <- range(or.vals[which(!is.infinite(or.vals) & !is.nan(or.vals) & !is.na(or.vals))])
  if(ylims[1] > 0){
    ylims[1] <- 0
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, n.bins), ylim=ylims, xaxt="n", yaxt="n", xlab="", ylab="")
  abline(h=c(0, baseline), lty=c(1, 2), col=c(null.color, baseline.color))
  text(x=par("usr")[1]-(0.05*(par("usr")[2]-par("usr")[1])),
       y=baseline+(0.035*(par("usr")[4]-par("usr")[3])), 
       labels="Baseline", cex=0.85, font=3, col=baseline.color, pos=4)
  
  # Add points
  pt.x.at <- (1:n.bins)-0.5
  segments(x0=pt.x.at, x1=pt.x.at, y0=ors$lower, y1=ors$upper, lend="round", 
           lwd=2.5, col=pt.color)
  points(x=pt.x.at, y=ors$estimate, pch=19, col=pt.color, cex=1.5)
  
  # Add X axis
  if(!is.null(x.ax.labels)){
    sapply(1:length(pt.x.at), function(i){
      if(parse.x.ax.labels==T){
        axis(1, at=pt.x.at[i], tick=F, line=-1, labels=parse(text=x.ax.labels[i]), 
             cex.axis=cex.x.ax.labels)
      }else{
        axis(1, at=pt.x.at[i], tick=F, line=-1, labels=x.ax.labels[i], 
             cex.axis=cex.x.ax.labels)
      }
    })
  }
  if(parse.x.label==T){
    mtext(1, text=parse(text=x.label), line=x.label.line)
  }else{
    mtext(1, text=x.label, line=x.label.line)
  }
  
  # Add Y axis
  y.ax.at <- axTicks(2)
  axis(2, at=c(-100, 100), tck=0, col=blueblack, labels=NA)
  axis(2, at=y.ax.at, col=blueblack, labels=NA, tck=-0.025)
  axis(2, at=y.ax.at, tick=F, labels=2^y.ax.at, las=2, line=-0.6)
  mtext(2, text=bquote("log"[2] * "(Odds Ratio for Autism)"), line=1.5)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(epitools, quietly=T)
require(flux, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog rCNV_scores.tsv constraint.meta.bed asd_cnvs.tsv asd_phenos.tsv out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: scores.tsv, constraint.meta.bed, asd_cnvs.tsv, asd_phenos.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
constraint.meta.in <- args$args[2]
asd_cnvs.in <- args$args[3]
asd_phenos.in <- args$args[4]
out.prefix <- args$args[5]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# constraint.meta.in <- "~/scratch/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz"
# asd_cnvs.in <- "~/scratch/asc_spark_denovo_cnvs.cleaned.b37.annotated.bed.gz"
# asd_phenos.in <- "~/scratch/asc_spark_child_phenos.list"
# out.prefix <- "~/scratch/test_gene_score_dnCNV_analyses"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load rCNV scores and other (non-rCNV) scores
scores.rCNV <- load.scores(scores.in)
scores.other <- load.other.scores(constraint.meta.in)
all.scores <- c(colnames(scores.rCNV)[-1], colnames(scores.other)[-1])

# Load ASD phenotypes & CNVs
phenos <- load.phenos(asd_phenos.in)
cnvs <- load.cnvs(asd_cnvs.in, phenos)

# Annotate CNVs with all scores
cnvs <- anno.cnv.scores(cnvs, scores.rCNV, scores.other)

# Plot odds ratios binned by pTS/pHI
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.odds_ratios.del.pdf", sep="."),
    height=2.75, width=2.1)
plot.or.by.scorebin(cnvs[which(cnvs$cnv=="DEL"), ], phenos, score="pHI", cnv.pct=T, n.bins=4, 
                    x.label="italic(\"De Novo\") ~ Deletions", 
                    parse.x.label=T, x.label.line=1.1,
                    x.ax.labels=paste("\"Q\"[", 1:4, "]", sep=""),
                    parse.x.ax.labels=T, cex.x.ax.labels=0.9, 
                    pt.color=cnv.colors[1], baseline.color=control.cnv.colors[1], 
                    null.color=blueblack, parmar=c(3.25, 2.75, 0.5, 0.5))
mtext(1, line=2.1, text="in Quartiles by pHI")
dev.off()
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.odds_ratios.dup.pdf", sep="."),
    height=2.75, width=2.1)
plot.or.by.scorebin(cnvs[which(cnvs$cnv=="DUP"), ], phenos, score="pTS", cnv.pct=T, n.bins=4, 
                    x.label="italic(\"De Novo\") ~ Duplications", 
                    parse.x.label=T, x.label.line=1.25,
                    x.ax.labels=paste("\"Q\"[", 1:4, "]", sep=""),
                    parse.x.ax.labels=T, cex.x.ax.labels=0.9, 
                    pt.color=cnv.colors[2], baseline.color=control.cnv.colors[2], 
                    null.color=blueblack, parmar=c(3.25, 2.75, 0.5, 0.5))
mtext(1, line=2.1, text="in Quartiles by pTS")

# Plot odds ratios stratified by high/low pHI & pTS
strat.ors <- calc.or.stratified(cnvs, phenos, high.cutoff=0.8, low.cutoff=0.5, norm.vs.baseline=F)
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.odds_ratios.stratified.pdf", sep="."),
    height=2.25, width=3)
plot.stratified.metric(strat.ors, y.title="\"log\"[2](\"ASD Odds Ratio\")")
dev.off()

# Evaluate performance of various scores vs. proband/sibling labels
del.evals <- lapply(all.scores, function(score){
  evaluate.score.asd_cnv(cnvs[which(cnvs$cnv=="DEL"), ], phenos, score=score)
})
dup.evals <- lapply(all.scores, function(score){
  evaluate.score.asd_cnv(cnvs[which(cnvs$cnv=="DUP"), ], phenos, score=score)
})

# Assign colors for score comparison
score.colors <- c(cnv.colors[1:2], rev(viridis(ncol(scores.other)-1)))
score.text.colors <- rev(c(rep("white", 2), 
                       rep("black", floor((ncol(scores.other)-1) / 2)), 
                       rep("white", ceiling((ncol(scores.other)-1) / 2))))

# Plot ROC
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.roc.del.pdf", sep="."),
    height=2.75, width=2.75)
plot.roc(del.evals, colors=score.colors, auc.text.colors=score.text.colors)
dev.off()
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.roc.dup.pdf", sep="."),
    height=2.75, width=2.75)
plot.roc(dup.evals, colors=score.colors, auc.text.colors=score.text.colors)
dev.off()

# Plot PRC
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.prc.del.pdf", sep="."),
    height=2.75, width=2.75)
plot.prc(del.evals, colors=score.colors, auc.text.colors=score.text.colors)
dev.off()
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.prc.dup.pdf", sep="."),
    height=2.75, width=2.75)
plot.prc(dup.evals, colors=score.colors, auc.text.colors=score.text.colors)
dev.off()

# Plot legend
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.roc_prc.legend.pdf", sep="."),
    height=2.2, width=2.7)
simple.legend(labels=rev(all.score.names[all.scores]), colors=rev(score.colors))
dev.off()

