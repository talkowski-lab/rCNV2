#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compare de novo CNVs from ASC/SPARK exomes vs. gene-level intolerance scores


options(stringsAsFactors=F, scipen=1000)
asc_spark_samplesizes = c("ASD" = 13694, "Control" = 5007)


######################
### DATA FUNCTIONS ###
######################
# Load ASC/SPARK de novo CNVs
load.cnvs <- function(asd_cnvs.in){
  read.table(asd_cnvs.in, header=T, sep="\t", comment.char="")[, -c(1:3)]
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
  cnvs$pHaplo <- anno.cnv.single.score(cnvs, scores.rCNV, "pHaplo")
  cnvs$pTriplo <- anno.cnv.single.score(cnvs, scores.rCNV, "pTriplo")
  cnvs.others <- sapply(colnames(scores.other)[-1], function(score){
    anno.cnv.single.score(cnvs, scores.other, score)
  })
  cnvs <- cbind(cnvs, cnvs.others)
  return(cnvs)
}

# Compute odds ratio for a given subset of CNVs
calc.or <- function(cnvs, quiet=TRUE, continuity.correction=0.5){
  n.case.all <- asc_spark_samplesizes["ASD"]
  n.ctrl.all <- asc_spark_samplesizes["Control"]
  n.case.cnv <- length(unique(cnvs$sample[which(cnvs$pheno == "ASD")]))
  n.ctrl.cnv <- length(unique(cnvs$sample[which(cnvs$pheno == "Control")]))
  n.case.ref <- n.case.all - n.case.cnv
  n.ctrl.ref <- n.ctrl.all - n.ctrl.cnv
  or.table <- matrix(c(n.ctrl.ref, n.case.ref, n.ctrl.cnv, n.case.cnv), nrow=2, byrow=T) + continuity.correction
  colnames(or.table) <- c("Control", "ASD")
  rownames(or.table) <- c("Ref", "CNV")
  if(!quiet){
    print(or.table - continuity.correction)
  }
  as.numeric(oddsratio.wald(or.table, correction=F)$measure[2, ])
}

# Compute odds ratios for CNVs binned by score (either percentile or absolute bin)
calc.or.by.scorebin <- function(cnvs, score, n.bins=3, cnv.pct=TRUE, marginal.ors=TRUE){
  bin.breaks <- seq(0, 1, length.out=n.bins+1)
  cnvs <- cnvs[which(!is.na(cnvs[, score])), ]
  if(cnv.pct==T){
    cnv.pct <- rank(cnvs[, score]) / nrow(cnvs)
  }else{
    cnv.pct <- cnvs[, score]
  }
  ors <- as.data.frame(t(sapply(1:n.bins, function(i){
    if(marginal.ors){
      bin.cnvs <- cnvs[which(cnv.pct >= bin.breaks[i] & cnv.pct <= bin.breaks[i+1]), ]
    }else{
      bin.cnvs <- cnvs[which(cnv.pct >= bin.breaks[i]), ]
    }
    calc.or(bin.cnvs)
  })))
  colnames(ors) <- c("estimate", "lower", "upper")
  return(ors)
}

# Calculate odds ratios for dnCNVs stratified by high/low pHaplo & pTriplo
calc.or.stratified <- function(cnvs, high.cutoff=0.8, low.cutoff=0.5,
                               log.trans=TRUE, norm.vs.baseline=TRUE){
  strat.ors <- lapply(c("DEL", "DUP"), function(cnvtype){
    ors <- t(data.frame("low.low"=calc.or(cnvs[which(cnvs$pHaplo<=low.cutoff & cnvs$pTriplo<=low.cutoff & cnvs$cnv==cnvtype), ]),
    "high.low"=calc.or(cnvs[which(cnvs$pHaplo>=high.cutoff & cnvs$pTriplo<=low.cutoff & cnvs$cnv==cnvtype), ]),
    "low.high"=calc.or(cnvs[which(cnvs$pHaplo<=low.cutoff & cnvs$pTriplo>=high.cutoff & cnvs$cnv==cnvtype), ]),
    "high.high"=calc.or(cnvs[which(cnvs$pHaplo>=high.cutoff & cnvs$pTriplo>=high.cutoff & cnvs$cnv==cnvtype), ])))
    colnames(ors) <- c("estimate", "lower", "upper")
    if(norm.vs.baseline==TRUE){
      baseline <- calc.or(cnvs[which(cnvs$cnv==cnvtype), ])[1]
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
evaluate.score.asd_cnv <- function(cnvs, score){
  # Add dummy column to CNVs to match rCNV2::roc() and rCNV2::prc()
  cnvs$gene <- cnvs$pheno
  roc.res <- rCNV2::roc(cnvs, score, "ASD", "Control")
  roc.auc <- flux::auc(roc.res$frac_false, roc.res$frac_true)
  prc.res <- rCNV2::prc(cnvs, score, "ASD", "Control")
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
plot.or.by.scorebin <- function(cnvs, score, n.bins=3, cnv.pct=TRUE, marginal.ors=TRUE, print.ors=TRUE,
                                pt.type="p", pt.lwd=1, pt.cex=1.5, pt.color=blueblack,
                                ci.type="bars", ci.color=pt.color,
                                add.baseline=T, baseline.color=blueblack, blue.bg=TRUE, null.color=bluewhite,
                                x.label=NULL, parse.x.label=T, x.label.line=2,
                                x.ax.labels=NULL, parse.x.ax.labels=T,
                                cex.x.ax.labels=1, x.ax.labels.line=-1,
                                fit.ylims.to="ci", ylim.buff.cex=1/3,
                                parmar=c(3.5, 3.5, 0.5, 0.5)){
  # Gather plot data
  baseline <- log2(calc.or(cnvs)[1])
  ors <- calc.or.by.scorebin(cnvs, score, n.bins, cnv.pct, marginal.ors)
  if(print.ors==TRUE){
    cat("ASD risk per score bin:\n")
    cat(paste(paste(round(ors[, 1], 3), collapse=", "), "\n"))
  }
  ors <- log2(ors)
  or.vals <- as.numeric(unlist(ors))
  if(fit.ylims.to == "ci"){
    ylims <- range(or.vals[which(!is.infinite(or.vals) & !is.nan(or.vals) & !is.na(or.vals))])
  }else if(fit.ylims.to == "mean"){
    ylims <- range(ors[which(!is.infinite(ors[, 1]) & !is.nan(ors[, 1]) & !is.na(ors[, 1])), 1])
    ylim.buff <- diff(ylims) * ylim.buff.cex
    ylims <- ylims + c(-ylim.buff, ylim.buff)
  }
  if(blue.bg==TRUE){
    plot.bg <- bluewhite
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- "white"
  }else{
    plot.bg <- "white"
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- NA
  }

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, n.bins), ylim=ylims, xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  y.ax.at <- -10:10
  abline(h=y.ax.at, col=grid.col)
  if(add.baseline){
    abline(h=c(0, baseline), lty=c(1, 2), col=c(null.color, baseline.color))
    text(x=par("usr")[1]-(0.05*(par("usr")[2]-par("usr")[1])),
         y=baseline+(0.035*(par("usr")[4]-par("usr")[3])),
         labels="Baseline", cex=0.85, font=3, col=baseline.color, pos=4)
  }

  # Add points
  pt.x.at <- (1:n.bins)-0.5
  if(ci.type == "bars"){
    segments(x0=pt.x.at, x1=pt.x.at, y0=ors$lower, y1=ors$upper, lend="round",
             lwd=2.5, col=ci.color)
  }else if(ci.type == "polygon"){
    polygon(x=c(pt.x.at, rev(pt.x.at)), y=c(ors$lower, rev(ors$upper)),
            bty="n", border=NA, col=ci.color)
  }
  points(x=pt.x.at, y=ors$estimate, type=pt.type,lwd=pt.lwd,
         pch=19, col=pt.color, cex=pt.cex, xpd=T)

  # Add X axis
  if(!is.null(x.ax.labels)){
    sapply(1:length(pt.x.at), function(i){
      if(parse.x.ax.labels==T){
        axis(1, at=pt.x.at[i], tick=F, line=x.ax.labels.line, labels=parse(text=x.ax.labels[i]),
             cex.axis=cex.x.ax.labels)
      }else{
        axis(1, at=pt.x.at[i], tick=F, line=x.ax.labels.line, labels=x.ax.labels[i],
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
  axis(2, at=c(-100, 100), tck=0, col=blueblack, labels=NA)
  axis(2, at=y.ax.at, col=blueblack, labels=NA, tck=-0.025)
  axis(2, at=y.ax.at, tick=F, labels=2^y.ax.at, las=2, line=-0.6)
  mtext(2, text="Odds Ratio for Autism", line=1.5)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(epitools, quietly=T)
require(flux, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog rCNV_scores.tsv constraint.meta.bed asd_cnvs.tsv out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, constraint.meta.bed, asd_cnvs.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
constraint.meta.in <- args$args[2]
asd_cnvs.in <- args$args[3]
out.prefix <- args$args[4]

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# constraint.meta.in <- "~/scratch/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz"
# asd_cnvs.in <- "~/scratch/asc_spark_2021_denovo_cnvs.cleaned.b37.annotated.nonNAHR.bed.gz"
# out.prefix <- "~/scratch/test_gene_score_dnCNV_analyses"

# Load rCNV scores and other (non-rCNV) scores
scores.rCNV <- load.scores(scores.in)
scores.other <- load.other.scores(constraint.meta.in)
all.scores <- c(colnames(scores.rCNV)[-1], colnames(scores.other)[-1])

# Load ASD CNVs
cnvs <- load.cnvs(asd_cnvs.in)

# Annotate CNVs with all scores
cnvs <- anno.cnv.scores(cnvs, scores.rCNV, scores.other)

# Plot odds ratios of CNVs ranked by pTriplo/pHaplo
lapply(list(c("DEL", "pHaplo", "Deletions", 2),
            c("DUP", "pTriplo", "Duplications", 2.15)),
       function(vals){
  pdf(paste(out.prefix, "asc_spark_denovo_cnvs.odds_ratios",
            tolower(vals[1]), "pdf", sep="."),
      height=3, width=2.25)
  plot.or.by.scorebin(cnvs[which(cnvs$cnv==vals[1]), ], score=vals[2], cnv.pct=T, n.bins=3,
                      x.label=paste("italic(\"De Novo\") ~", vals[3]),
                      parse.x.label=T, x.label.line=vals[4],
                      x.ax.labels=paste(c("Bottom", "Mid.", "Top"), "33%", sep="\n"),
                      parse.x.ax.labels=F, cex.x.ax.labels=0.9, x.ax.labels.line=-0.25,
                      pt.color=cnv.colors[vals[1]], baseline.color=control.cnv.colors[vals[1]],
                      null.color=blueblack, blue.bg=FALSE,
                      parmar=c(4, 2.5, 0.5, 0.25))
  mtext(1, line=2.9, text=paste("Ranked by", vals[2]))
  dev.off()
})

# Plot odds ratios of CNVs vs. sliding pTriplo/pHaplo thresholds
lapply(list(c("pHaplo", "DEL", 0.2),
            c("pTriplo", "DUP", 0.5)),
       function(vals){
pdf(paste(out.prefix, ".asc_spark_denovo_cnvs.odds_ratios.sliding_",
          vals[1], "_cutoff.", tolower(vals[2]), ".pdf", sep=""),
    height=2, width=2.5)
plot.or.by.scorebin(cnvs[which(cnvs$cnv==vals[2]), ], score=vals[1], cnv.pct=F,
                    marginal.ors=F, n.bins=100, print.ors=F, pt.type="l",
                    pt.lwd=3, pt.color=cnv.colors[vals[2]], ci.type="polygon",
                    ci.color=cnv.whites[vals[2]], x.label=paste(vals[1], "Cutoff"),
                    parse.x.label=F, x.label.line=1.1, add.baseline=F,
                    fit.ylims.to="mean", ylim.buff.cex=as.numeric(vals[3]),
                    null.color=blueblack, blue.bg=FALSE,
                    parmar=c(2.25, 2.5, 0.5, 0.25))
axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
axis(1, at=seq(0, 100, 25), tick=F, line=-0.9, labels=seq(0, 1, 0.25), cex.axis=0.9)
dev.off()
})

# Plot odds ratios stratified by high/low pHaplo & pTriplo
strat.ors <- calc.or.stratified(cnvs, high.cutoff=0.8, low.cutoff=0.5, norm.vs.baseline=F)
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.odds_ratios.stratified.pdf", sep="."),
    height=2.25, width=3)
plot.stratified.metric(strat.ors, y.title="\"log\"[2](\"ASD Odds Ratio\")")
dev.off()

# Evaluate performance of various scores vs. proband/sibling labels
del.evals <- lapply(all.scores, function(score){
  evaluate.score.asd_cnv(cnvs[which(cnvs$cnv=="DEL"), ], score=score)
})
dup.evals <- lapply(all.scores, function(score){
  evaluate.score.asd_cnv(cnvs[which(cnvs$cnv=="DUP"), ], score=score)
})

# Assign colors for score comparison
score.colors <- c(cnv.colors[1:2], rev(viridis(ncol(scores.other)-1)))
score.text.colors <- rev(c(rep("white", 2),
                       rep("black", floor((ncol(scores.other)-1) / 2)),
                       rep("white", ceiling((ncol(scores.other)-1) / 2))))

# Plot ROC
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.roc.del.pdf", sep="."),
    height=2.75, width=2.75)
plot.roc(del.evals, colors=score.colors, auc.text.colors=score.text.colors, grid.col=NA)
dev.off()
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.roc.dup.pdf", sep="."),
    height=2.75, width=2.75)
plot.roc(dup.evals, colors=score.colors, auc.text.colors=score.text.colors, grid.col=NA)
dev.off()

# Plot PRC
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.prc.del.pdf", sep="."),
    height=2.75, width=2.75)
plot.prc(del.evals, colors=score.colors, auc.text.colors=score.text.colors, grid.col=NA)
dev.off()
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.prc.dup.pdf", sep="."),
    height=2.75, width=2.75)
plot.prc(dup.evals, colors=score.colors, auc.text.colors=score.text.colors, grid.col=NA)
dev.off()

# Plot legend
pdf(paste(out.prefix, "asc_spark_denovo_cnvs.roc_prc.legend.pdf", sep="."),
    height=2.2, width=2.7)
simple.legend(labels=rev(score.names[all.scores]), colors=rev(score.colors))
dev.off()

