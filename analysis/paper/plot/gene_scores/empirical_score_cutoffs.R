#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Determine optimal pHI & pTS score cutoffs based on empirical odds ratios


# Set parameters & define constants
options(stringsAsFactors=F, scipen=1000, family="sans")


######################
### DATA FUNCTIONS ###
######################
# Load meta-analysis stats and extract effect size estimates w/relative variance
load.lnors <- function(meta.in, xlist){
  meta <- read.table(meta.in, sep="\t", header=T, comment.char="")
  meta[which(!(meta$gene %in% xlist)), c("gene", "meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper")]
  meta$meta_lnOR_se <- (meta$meta_lnOR_upper - meta$meta_lnOR_lower) / (2 * qnorm(0.975))
  meta$meta_lnOR_var <- meta$meta_lnOR_se^2
  return(meta[, c("gene", "meta_lnOR", "meta_lnOR_var")])
}

# Compute inverse-variance weighted mean and 95% CI
inv.var.avg <- function(lnors, vars, conf=0.95){
  keep.idx <- which(complete.cases(data.frame(lnors, vars)))
  numerator <- sum((lnors/vars)[keep.idx])
  denominator <- sum(1/vars[keep.idx])
  avg <- numerator/denominator
  pooled.se <- sqrt(1/denominator)
  lower.ci <- avg + qnorm((1-conf)/2)*pooled.se
  upper.ci <- avg + qnorm(conf+(1-conf)/2)*pooled.se
  return(c(avg, lower.ci, upper.ci))
}

# Compute average effect size per score bin
avg.lnor.per.score.bin <- function(lnor.df, scores, score, bins=100, 
                                   start=1, end=0){
  score.breaks <- seq(start, end, length.out=bins+1)
  res <- as.data.frame(t(sapply(1:bins, function(i){
    upper.genes <- scores$gene[which(scores[, score] >= score.breaks[i+1])]
    n.upper.genes <- length(upper.genes)
    upper.lnor <- inv.var.avg(lnor.df$meta_lnOR[which(lnor.df$gene %in% upper.genes)],
                             lnor.df$meta_lnOR_var[which(lnor.df$gene %in% upper.genes)])
    bin.genes <- intersect(upper.genes, scores$gene[which(scores[, score] < score.breaks[i])])
    n.bin.genes <- length(bin.genes)
    bin.lnor <- inv.var.avg(lnor.df$meta_lnOR[which(lnor.df$gene %in% bin.genes)],
                             lnor.df$meta_lnOR_var[which(lnor.df$gene %in% bin.genes)])
    return(c(score.breaks[i+1], n.upper.genes, upper.lnor, n.bin.genes, bin.lnor))
  })))
  colnames(res) <- c("cutoff", 
                     "cumul_genes", "cumul_lnOR", "cumul_lnOR_lower", "cumul_lnOR_upper",
                     "margin_genes", "margin_lnOR", "margin_lnOR_lower", "margin_lnOR_upper")
  return(res)
}

# Derive comparable score cutoff vs constrained gene deletions
get.cutoff <- function(score.lnor.full, score.lnor.fine, min.lnor){
  full.idxs <- which(score.lnor.full$cumul_lnOR >= min.lnor)
  if(length(full.idxs) > 0){
    return(score.lnor.full$cutoff[max(full.idxs)])
  }else{
    fine.idxs <- which(score.lnor.fine$cumul_lnOR >= min.lnor)
    return(score.lnor.fine$cutoff[max(fine.idxs)])
  }
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot empirical effect sizes vs score bin
plot.score.vs.or <- function(bins, baseline, score.cutoff, score,
                             avg.pt.cex=1, avg.genes.per.bin=NULL, 
                             x.ax.at=NULL, ylims=NULL, xtitle=NULL, ytitle=NULL,
                             ax.tck=-0.025, parmar=c(2.4, 2.4, 0.25, 0.25)){
  # Get plot values
  x <- bins$cutoff
  pt.y <- bins$margin_lnOR
  line.y <- bins$cumul_lnOR

  # Set plot parameters
  if(score == "pHI"){
    col.idx <- 1
    cnv.label <- "Deletion"
  }else{
    col.idx <- 2
    cnv.label <- "Duplication"
  }
  pt.col <- rep(control.cnv.colors[col.idx], nrow(bins))
  line.col <- cnv.colors[col.idx]
  baseline.col <- graphabs.green
  highlight.idxs <- which(bins$cutoff >= score.cutoff)
  pt.col[highlight.idxs] <- cnv.colors[col.idx]
  highlight.line.col <- cnv.blacks[col.idx]
  if(is.null(avg.genes.per.bin)){
    avg.genes.per.bin <- mean(bins$margin_genes, na.rm=T)
  }
  pt.cex <- avg.pt.cex*sqrt(bins$margin_genes/avg.genes.per.bin)
  xlims <- range(x, na.rm=T)
  if(is.null(ylims)){
    ylims <- range(c(pt.y, line.y), na.rm=T)
  }
  if(is.null(xtitle)){
    xtitle <- paste("Genes Binned by", score)
  }
  if(is.null(ytitle)){
    ytitle <- paste(cnv.label, "Odds Ratio")
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=score.cutoff,
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=bluewhite, border=NA, bty="n")
  rect(xleft=score.cutoff, xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=adjustcolor(highlight.color, alpha=0.6), border=NA, bty="n")
  y.ax.at <- log(2^(-10:10))
  abline(h=y.ax.at, col="white")
  abline(v=score.cutoff, lty=1, col=highlight.color)
  
  # Add points & lines
  abline(h=baseline, lty=5, col=baseline.col)
  points(x, pt.y, col=pt.col, cex=pt.cex, pch=19)
  points(x, line.y, type="l", col=line.col, lwd=2)
  points(x[highlight.idxs], line.y[highlight.idxs], 
         type="l", col=highlight.line.col, lwd=2)
  points(x=score.cutoff, y=baseline, pch=23, cex=1.5*avg.pt.cex, 
         bg=highlight.color, col=baseline.col, lwd=2)
  
  # Add axes
  if(is.null(x.ax.at)){
    x.ax.at <- axTicks(1)
  }
  axis(1, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(1, at=x.ax.at, tck=ax.tck, labels=NA, col=blueblack)
  axis(1, at=x.ax.at, tick=F, line=-0.7)
  mtext(1, line=1.3, text=xtitle, xpd=T)
  axis(2, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(2, at=y.ax.at, tck=ax.tck, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, tick=F, labels=exp(y.ax.at), las=2, line=-0.6)
  mtext(2, line=1.2, text=ytitle)
  axis(3, at=c(score.cutoff, par("usr")[2]), col=blueblack, tck=-ax.tck)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv meta_stats.del.tsv meta_stats.dup.tsv constr.genes exclude.list out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Six positional arguments required: scores.tsv, meta_stats.del.tsv, meta_stats.dup.tsv, constr.genes.list, exclude.genes, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
del.meta.in <- args$args[2]
dup.meta.in <- args$args[3]
constr.genes.in <- args$args[4]
xlist.in <- args$args[5]
out.prefix <- args$args[6]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# del.meta.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
# dup.meta.in <- "~/scratch/rCNV2_analysis_d1.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz"
# constr.genes.in <- "~/scratch/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
# xlist.in <- "~/scratch/rCNV.gene_scoring.training_gene_blacklist.bed.gz"
# out.prefix <- "~/scratch/test_gene_score_empirical_cutoffs"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load scores
scores <- load.scores(scores.in)

# Load gene lists
constr.genes <- unique(as.character(read.table(constr.genes.in, header=F, sep="\t")[, 1]))
xlist <- unique(as.character(read.table(xlist.in, header=F, sep="\t")[, 1]))

# Extract meta-analysis odds ratios
del.lnors <- load.lnors(del.meta.in, xlist)
dup.lnors <- load.lnors(dup.meta.in, xlist)

# Compute average lnOR for deletions of constrained genes
del.constr.idxs <- which(del.lnors$gene %in% constr.genes & !(del.lnors$gene %in% xlist))
del.constr.lnor <- inv.var.avg(del.lnors$meta_lnOR[del.constr.idxs], del.lnors$meta_lnOR_var[del.constr.idxs])

# Compute lnOR estimates for genes by pHI & pTS bin
phi.lnor.full <- avg.lnor.per.score.bin(del.lnors, scores, "pHI", bins=100, start=1, end=0)
phi.lnor.fine <- avg.lnor.per.score.bin(del.lnors, scores, "pHI", bins=100, start=1, end=0.9)
pts.lnor.full <- avg.lnor.per.score.bin(dup.lnors, scores, "pTS", bins=100, start=1, end=0)
pts.lnor.fine <- avg.lnor.per.score.bin(dup.lnors, scores, "pTS", bins=100, start=1, end=0.9)
pts.lnor.fine.forplot <- avg.lnor.per.score.bin(dup.lnors, scores, "pTS", bins=50, start=1, end=0.95)

# Derive cutoffs
phi.cutoff <- get.cutoff(phi.lnor.full, phi.lnor.fine, del.constr.lnor[1])
pts.cutoff <- get.cutoff(pts.lnor.full, pts.lnor.fine, del.constr.lnor[1])

# Print derived cutoffs
cat(paste("\nAverage rare deletion of constrained gene confers odds ratio =",
            round(exp(del.constr.lnor[1]), 4), "\n"))
cat(paste("\nComparable pHI cutoff >=", phi.cutoff, "(includes",
          prettyNum(length(which(scores$pHI>=phi.cutoff)), big.mark=","), 
          "genes)\n"))
cat(paste("\nComparable pTS cutoff >=", pts.cutoff, "(includes",
          prettyNum(length(which(scores$pTS>=pts.cutoff)), big.mark=","), 
          "genes)\n"))

# Set standardized parameters for all plots
ylims <- quantile(c(phi.lnor.full$margin_lnOR, pts.lnor.full$margin_lnOR, pts.lnor.fine$margin_lnOR), 
                  probs=c(0.005, 0.995), na.rm=T)
pt.cex <- 3/4
avg.genes.per.bin <- mean(c(phi.lnor.full$margin_genes, pts.lnor.full$margin_genes), na.rm=T)
pdf.height <- 2.25
pdf.width <- 2.75
pdf.parmar <- c(2.3, 2.3, 0.3, 0.3)

# Plot full pHI vs lnOR
pdf(paste(out.prefix, "phi_vs_effect_size.full.pdf", sep="."),
    height=pdf.height, width=pdf.width)
plot.score.vs.or(phi.lnor.full, del.constr.lnor[1], phi.cutoff, "pHI", 
                 ylims=ylims, avg.pt.cex=pt.cex, avg.genes.per.bin=avg.genes.per.bin,
                 parmar=pdf.parmar)
dev.off()

# Plot full pTS vs lnOR
pdf(paste(out.prefix, "pts_vs_effect_size.full.pdf", sep="."),
    height=pdf.height, width=pdf.width)
plot.score.vs.or(pts.lnor.full, del.constr.lnor[1], pts.cutoff, "pTS", 
                 ylims=ylims, avg.pt.cex=pt.cex, avg.genes.per.bin=avg.genes.per.bin,
                 parmar=pdf.parmar)
# text(x=par("usr")[1], y=del.constr.lnor[1]+(0.08*diff(par("usr")[3:4])), 
#      labels="Rare deletions of\nall constrained genes",
#      font=3, cex=5/6, col=graphabs.green, pos=4)
dev.off()

# Plot fine pTS vs lnOR
pdf(paste(out.prefix, "pts_vs_effect_size.fine.pdf", sep="."),
    height=(5/6)*pdf.height, width=(1/2)*pdf.width)
plot.score.vs.or(pts.lnor.fine.forplot, del.constr.lnor[1], pts.cutoff, "pTS", 
                 ylims=ylims, avg.pt.cex=pt.cex, xtitle="pTS", ytitle=NA, 
                 parmar=c(pdf.parmar[1], 1.3, pdf.parmar[3:4]))
axis(1, at=1, tick=F, line=-0.7)
dev.off()

