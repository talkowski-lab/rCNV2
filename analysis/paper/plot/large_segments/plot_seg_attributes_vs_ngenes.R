#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot relationships vs. number of genes per segment for rCNV2 formal analyses


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Split values by CNV & bins of number of genes
quantile.split <- function(segs, feature, breaks=NULL, break.names=NULL, probs=seq(0, 1, 0.2)){
  vals <- segs[, which(colnames(segs)==feature)]
  ngenes <- segs$n_genes
  if(is.null(breaks)){
    breaks <- floor(quantile(ngenes, probs=probs, na.rm=T))
    breaks[length(breaks)] <- breaks[length(breaks)] + 1
  }
  if(is.null(break.names)){
    break.names <- sapply(2:length(breaks), function(i){
      paste(breaks[i-1], breaks[i], sep="-")
    })
  }

  res <- lapply(c("DEL", "DUP"), function(cnv){
    cnv.vals <- vals[which(segs$cnv == cnv)]
    cnv.ngenes <- ngenes[which(segs$cnv == cnv)]
    cnv.res <- lapply(2:length(breaks), function(i){
      cnv.vals[which(cnv.ngenes>=breaks[i-1] & cnv.ngenes<breaks[i])]
    })
    names(cnv.res) <- break.names
    return(cnv.res)
  })
  names(res) <- c("DEL", "DUP")

  return(res)
}

# Compute Chi-square P-values for obs vs. expected constrained genes
constrained.obs_exp.test <- function(segs){
  segs$constrained_obs_exp_pvalue <- sapply(1:nrow(segs), function(i){
    n_genes <- segs$n_genes[i]
    exp <- segs$constrained_expected[i]
    obs <- segs$n_gnomAD_constrained_genes[i]
    chisq.mat <- matrix(c(n_genes - exp, n_genes - obs, exp, obs),
                        byrow=T, nrow=2)
    if(n_genes > 0){
      chisq.test(chisq.mat)$p.value
    }else{
      return(1)
    }
  })
  return(segs)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot values vs # of genes
scatter.vsGenes <- function(segs, feature, pt.cex=0.6, fit=NULL, rolling.span=1,
                            xlims=NULL, ylims=NULL, y.title=NULL,
                            y.title.line=2, y.pct=FALSE,
                            horiz.line=NULL, legend.pos=NULL, blue.bg=TRUE,
                            parmar=c(2.5, 3, 0.5, 0.5)){
  # Get plot values
  segs <- segs[which(segs$n_genes>0), ]
  ngenes <- segs$n_genes
  vals <- segs[, which(colnames(segs) == feature)]
  if(is.null(xlims)){
    xlims <- c(0, max(ngenes, na.rm=T))
  }
  if(is.null(ylims)){
    ylims <- range(vals, na.rm=T)
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
  plot(NA, xlim=xlims, ylim=ylims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg, xpd=T)
  abline(h=axTicks(2), v=axTicks(1), col=grid.col)
  if(!is.null(horiz.line)){
    abline(h=horiz.line, lty=2, col=blueblack)
  }

  # Add points
  points(ngenes, vals,
         pch=segs$pt.pch, cex=pt.cex,
         bg=segs$pt.bg, col=segs$pt.border)

  # Fit trendlines for DEL and DUP, if optioned
  if(!is.null(fit)){
    if(fit == "exp"){
      lapply(c("DEL", "DUP"), function(cnv){
        fit.model <- fit.exp.decay(x=ngenes[which(segs$cnv==cnv)],
                                 y=vals[which(segs$cnv==cnv)])
        fit.df <- data.frame("x"=seq(par("usr")[1], par("usr")[2], length.out=500))
        fit.df$y <- predict(fit.model, newdata=fit.df)
        points(fit.df, type="l", lwd=2, col=cnv.colors[which(names(cnv.colors)==cnv)])
      })
    }else if(fit %in% c("loess", "rollmean")){
      lapply(c("DEL", "DUP"), function(cnv){
        fit.df <- data.frame("x"=ngenes[which(segs$cnv==cnv)],
                             "y"=vals[which(segs$cnv==cnv)])
        fit.df <- fit.df[order(fit.df$x), ]
        if(fit == "loess"){
          fit <- loess(y ~ x, data=fit.df, span=rolling.span)
          points(x=fit.df$x, y=predict(fit), type="l", lwd=2, col=cnv.colors[which(names(cnv.colors)==cnv)])
        }else{
          fit <- rollapply(fit.df$y, rolling.span, mean, fill=NA, partial=FALSE)
          points(x=fit.df$x, y=fit, type="l", lwd=2, col=cnv.colors[which(names(cnv.colors)==cnv)])
        }
      })
    }
  }

  # Add axes
  x.ax.at <- axTicks(1)
  if(length(x.ax.at) > 6){
    x.ax.at <- x.ax.at[seq(1, length(x.ax.at), by=2)]
  }
  axis(1, at=c(-10e10, 10e10), col=blueblack, labels=NA)
  axis(1, at=x.ax.at, labels=NA, tck=-0.03, col=blueblack)
  sapply(x.ax.at, function(x){
    axis(1, at=x, tick=F, labels=prettyNum(x, big.mark=","), line=-0.65)
  })
  mtext(1, line=1.3, text="Genes in Segment")
  axis(2, at=c(-10e10, 10e10), col=blueblack, labels=NA)
  axis(2, labels=NA, tck=-0.03, col=blueblack)
  if(y.pct==T){
    y.labs <- paste(round(100*axTicks(2), 0), "%", sep="")
  }else{
    y.labs <- axTicks(2)
  }
  axis(2, at=axTicks(2), labels=y.labs, tick=F, las=2, line=-0.65)
  mtext(2, line=y.title.line, text=y.title)

  # Add legend, if optioned
  if(!is.null(legend.pos)){
    legend(legend.pos, legend=c("DEL", "DUP", "GW-Sig.", "Literature"),
           pch=c(19, 19, 22, 21), col=c(cnv.colors[1:2], "black", "gray50"),
           pt.bg=c(NA, NA, "gray30", "gray70"),
           box.col=blueblack, bg="white", cex=0.8, xpd=T)
  }
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(MASS, quietly=T)
require(zoo, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: loci.bed, segs.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
out.prefix <- args$args[3]

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# out.prefix <- "~/scratch/final_segs"

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs[which(segs$any_gd | segs$any_sig), ]

# Set expected values of constrained genes
prop_constrained.genome_avg <- 3036/18641
segs.all$constrained_expected <- segs.all$n_genes * prop_constrained.genome_avg

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]

# Subset to loci either in top or bottom third of effect size distribution
lnor.groups <- split.regions.by.effect.size(segs, quantiles=3)
segs.bylnor <- segs[which(segs$region_id %in% c(lnor.groups[[1]], lnor.groups[[3]])), ]

# Merge loci & segment data for genome-wide/FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs)

# Get list of neuro loci (sig + lit GDs)
neuro.plus.lit.ids <- get.neuro.region_ids(loci, segs.all)
neuro.segs <- segs[which(segs$region_id %in% neuro.plus.lit.ids), ]

# Get list of developmental loci (sig + lit GDs)
dev.plus.lit.ids <- get.developmental.region_ids(loci, segs.all)
dev.segs <- segs[which(segs$region_id %in% dev.plus.lit.ids), ]

# Get list of NDD loci (sig + lit GDs)
NDD.plus.lit.ids <- get.ndd.region_ids(loci, segs.all)
NDD.segs <- segs[which(segs$region_id %in% NDD.plus.lit.ids), ]

# Prepare output directory
outdir <- paste(dirname(out.prefix), "segs_vs_n_genes", sep="/")
if(!dir.exists(outdir)){
  dir.create(outdir)
}
out.prefix <- paste(outdir, "/", basename(out.prefix), sep="/")

# Plot proportion of constrained genes vs. # of genes
pdf(paste(out.prefix, "prop_constrained_vs_ngenes.all_phenos.pdf", sep="."),
    height=2.4, width=2.7)
scatter.vsGenes(segs, feature="gnomAD_constrained_prop", y.title="Constrained Genes",
                y.pct=T, pt.cex=0.6, blue.bg=FALSE,
                horiz.line=prop_constrained.genome_avg,
                y.title.line=2.3, parmar=c(2.45, 3.25, 0.5, 0.5))
dev.off()
pdf(paste(out.prefix, "prop_constrained_vs_ngenes.developmental_only.pdf", sep="."),
    height=2.4, width=2.7)
scatter.vsGenes(dev.segs, feature="gnomAD_constrained_prop", y.title="Constrained Genes",
                y.pct=T, pt.cex=0.6, blue.bg=FALSE,
                horiz.line=prop_constrained.genome_avg,
                y.title.line=2.3, parmar=c(2.45, 3.25, 0.5, 0.5))
dev.off()

# Loop over exome cohorts and generate enrichment plots
for(cohort in c("ASC", "DDD")){
  cohort.prefix <- paste(out.prefix, cohort, sep=".")
  dnm.plot.dims <- c(2.2, 2.3)
  if(cohort == "ASC"){
    seg.df <- NDD.segs
  }else{
    seg.df <- dev.segs
  }

  # Loop over consequences
  for(csq in c("lof", "mis")){

    # Plot enrichment while including all genes
    pdf(paste(cohort.prefix, csq, "vs_ngenes.pdf", sep="_"),
        height=dnm.plot.dims[1], width=dnm.plot.dims[2])
    scatter.vsGenes(seg.df,
                    feature=paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_"),
                    y.title=bquote("Excess" ~ italic("dn") * "PTVs / Gene"),
                    y.pct=F, pt.cex=0.6, horiz.line=0, blue.bg=FALSE,
                    y.title.line=1.75, parmar=c(2.5, 2.75, 0.5, 0.5))
    dev.off()

    # Plot enrichment while holding-out genes significant from each study
    pdf(paste(cohort.prefix, csq, "vs_ngenes.no_exomeSig_genes.pdf", sep="_"),
        height=dnm.plot.dims[1], width=dnm.plot.dims[2])
    scatter.vsGenes(seg.df,
                    feature=paste(cohort, "noSig_dnm", csq, "norm_excess_per_gene", sep="_"),
                    y.title=bquote("Excess" ~ italic("dn") * "PTVs / Gene"),
                    y.pct=F, pt.cex=0.6, horiz.line=0, blue.bg=FALSE,
                    y.title.line=1.75, parmar=c(2.5, 2.75, 0.5, 0.5))
    dev.off()
  }
}

# Number of constrained genes vs expected
constrained.max <- max(segs$constrained_expected, segs$n_gnomAD_constrained_genes)
pdf(paste(out.prefix, "constrained_genes_obs_vs_exp.pdf", sep="."),
    height=2.4, width=2.4)
segs.scatter(segs, segs$constrained_expected, segs$n_gnomAD_constrained_genes,
             xlims=c(0, constrained.max), ylims=c(0, constrained.max),
             xtitle="Constrained (Expected)", ytitle="Constrained (Observed)",
             pt.cex=0.6, blue.bg=FALSE)
dev.off()

