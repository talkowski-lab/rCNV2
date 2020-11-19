#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
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
scatter.vsGenes <- function(segs, feature, pt.cex=0.75, fit=NULL, rollmean.span=1,
                            xlims=NULL, ylims=NULL, y.title=NULL, 
                            y.title.line=2, y.pct=FALSE,
                            horiz.line=NULL, legend.pos=NULL,
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
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, bty="n", col=bluewhite)
  abline(h=axTicks(2), v=axTicks(1), col="white")
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
    }else if(fit == "rollmean"){
      lapply(c("DEL", "DUP"), function(cnv){
        fit.df <- data.frame("x"=ngenes[which(segs$cnv==cnv)],
                             "y"=vals[which(segs$cnv==cnv)])
        fit.df <- fit.df[order(fit.df$x), ]
        fit <- loess(y ~ x, data=fit.df, span=rollmean.span)
        points(x=fit.df$x, y=predict(fit), type="l", lwd=2, col=cnv.colors[which(names(cnv.colors)==cnv)])
      })
    }
  }
  
  # Add axes
  axis(1, at=c(-10e10, 10e10), col=blueblack, labels=NA)
  axis(1, labels=NA, tck=-0.03, col=blueblack)
  sapply(axTicks(1), function(x){
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
require(funr, quietly=T)
require(MASS, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

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
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# out.prefix <- "~/scratch/test_effect_sizes"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load loci & segment tables
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Subset to nominally significant segments
segs <- segs[which(segs$nom_sig), ]

# Set expected values of constrained genes
prop_constrained.genome_avg <- 3036/18641
segs$constrained_expected <- segs$n_genes * prop_constrained.genome_avg

# Get list of neuro loci plus lit GDs
neuro.plus.lit.ids <- get.neuro.region_ids(loci, segs)
neuro.segs <- segs[which(segs$region_id %in% neuro.plus.lit.ids), ]

# Plot proportion of constrained genes vs. # of genes
pdf(paste(out.prefix, "prop_constrained_vs_ngenes.pdf", sep="."),
    height=2.4, width=2.8)
scatter.vsGenes(segs, feature="gnomAD_constrained_prop", y.title="Constrained Genes",
                y.pct=T, pt.cex=0.85,
                horiz.line=prop_constrained.genome_avg,
                y.title.line=2.3, parmar=c(2.45, 3.25, 0.5, 3.5))
axis(4, at=prop_constrained.genome_avg, tck=-0.03, col=blueblack, labels=NA)
axis(4, at=prop_constrained.genome_avg, tick=F, line=-0.9, las=2,
     labels=c("Genome\naverage"), col.axis=blueblack, font=3, cex=0.8)
dev.off()

# NOTE: all DNM enrichment plots restricted to neuro-associated regions

# Plot average enrichment of ASC DNMs vs. # of genes
pdf(paste(out.prefix, "ASC_dnPTVs_vs_ngenes.pdf", sep="."),
    height=2.2, width=2.3)
scatter.vsGenes(neuro.segs, feature="ASC_dnm_lof_norm_excess_per_gene", 
                y.title=bquote("Excess" ~ italic("dn") * "PTVs / Gene"),
                y.pct=F, pt.cex=0.85, horiz.line=0,
                y.title.line=1.75, parmar=c(2.5, 2.75, 0.5, 0.5))
dev.off()
pdf(paste(out.prefix, "ASC_dnMis_vs_ngenes.pdf", sep="."),
    height=2.2, width=2.3)
scatter.vsGenes(neuro.segs, feature="ASC_dnm_mis_norm_excess_per_gene", 
                y.title=bquote("Excess" ~ italic("dn") * "Mis. / Gene"),
                y.pct=F, pt.cex=0.85, horiz.line=0,
                y.title.line=1.75, parmar=c(2.5, 2.75, 0.5, 0.5))
dev.off()

# Plot average enrichment of DDD DNMs vs. # of genes
pdf(paste(out.prefix, "DDD_dnPTVs_vs_ngenes.pdf", sep="."),
    height=2.2, width=2.3)
scatter.vsGenes(neuro.segs, feature="DDD_dnm_lof_norm_excess_per_gene", 
                y.title=bquote("Excess" ~ italic("dn") * "PTVs / Gene"),
                y.pct=F, pt.cex=0.85, horiz.line=0,
                y.title.line=1.75, parmar=c(2.5, 2.75, 0.5, 0.5))
dev.off()
pdf(paste(out.prefix, "DDD_dnMis_vs_ngenes.pdf", sep="."),
    height=2.2, width=2.3)
scatter.vsGenes(neuro.segs, feature="DDD_dnm_mis_norm_excess_per_gene", 
                y.title=bquote("Excess" ~ italic("dn") * "Mis. / Gene"),
                y.pct=F, pt.cex=0.85, horiz.line=0,
                y.title.line=1.75, parmar=c(2.5, 2.75, 0.5, 0.5))
dev.off()

# Number of constrained genes vs expected
constrained.max <- max(segs$constrained_expected, segs$n_gnomAD_constrained_genes)
pdf(paste(out.prefix, "constrained_genes_obs_vs_exp.pdf", sep="."),
    height=2.4, width=2.4)
segs.scatter(segs, segs$constrained_expected, segs$n_gnomAD_constrained_genes, 
             xlims=c(0, constrained.max), ylims=c(0, constrained.max),
             xtitle="Constrained (Expected)", ytitle="Constrained (Observed)",
             pt.cex=0.9)
dev.off()
