#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot basic distributions for large segments for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--gw-sig"), default=10e-6, help="P-value cutoff for genome-wide signficance.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: loci.bed, segs.tsv, output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
out.prefix <- args$args[3]
gw.sig <- -log10(opts$`gw-sig`)

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# out.prefix <- "~/scratch/final_segs"
# gw.sig <- -log10(0.000003742655)

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]

# Subset to loci either in top or bottom third of effect size distribution
lnor.groups <- split.regions.by.effect.size(segs, quantiles=3)
segs.bylnor <- segs[which(segs$region_id %in% c(lnor.groups[[1]], lnor.groups[[3]])), ]

# Merge loci & segment data for genome-wide/FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs)


# Loop over all pairwise comparisons for vioswarm plots
# Write each comparison to its own subdirectory
for(comp in c("DEL_vs_DUP", "gw_vs_FDR", "NAHR_vs_nonrecurrent",
              "terminal_vs_interstitial", "dev_vs_adult",
              "strong_vs_weak", "sig_vs_litGDs", "pleiotropy")){

  # Set plotting variables
  boxplot.colors <- NULL
  boxplot.fill <- NULL
  if(comp == "DEL_vs_DUP"){
    seg.df <- segs
    x.bool <- seg.df$cnv == "DUP"
    x.labs <- c("DEL", "DUP")
    boxplot.colors <- cnv.blacks[1:2]
    boxplot.fill <- cnv.whites[1:2]
  }else if(comp == "gw_vs_FDR"){
    seg.df <- segs.sig
    x.bool <- !(seg.df$gw_sig)
    x.labs <- c("GW Sig.", "FDR Sig.")
  }else if(comp == "NAHR_vs_nonrecurrent"){
    seg.df <- segs
    x.bool <- seg.df$nahr
    x.labs <- c("Nonrecurrent", "NAHR")
  }else if(comp == "terminal_vs_interstitial"){
    seg.df <- segs
    x.bool <- seg.df$terminal
    x.labs <- c("Interstitial", "Terminal")
  }else if(comp == "dev_vs_adult"){
    seg.df <- segs.sig
    x.bool <- !(seg.df$region_id %in% get.developmental.region_ids(loci, seg.df, sig.only=TRUE))
    x.labs <- c("Develop.", "Adult")
  }else if(comp == "strong_vs_weak"){
    seg.df <- segs.bylnor
    x.bool <- !(segs.bylnor$region_id %in% lnor.groups[[1]])
    x.labs <- c("Weak", "Strong")
  }else if(comp == "sig_vs_litGDs"){
    seg.df <- segs.all
    x.bool <- !(seg.df$any_sig)
    x.labs <- c("GW+FDR", "Lit. GDs")
  }else if(comp == "pleiotropy"){
    seg.df <- segs.sig
    x.bool <- seg.df$n_hpos > 1
    x.labs <- c("Single HPO", ">1 HPO")
  }

  # Create output subdirectory, if it doesn't exist
  out.subdir <- paste(dirname(out.prefix), comp, sep="/")
  if(!dir.exists(out.subdir)){
    dir.create(out.subdir)
  }
  out.prefix.sub.fname <- paste(basename(out.prefix), comp, sep=".")
  out.prefix.sub <- paste(out.subdir, out.prefix.sub.fname, sep="/")

  # Generate two plots per comparison, one for all CNV types and one split by CNV type
  for(cnv.split in c(TRUE, FALSE)){
    if(cnv.split){
      if(comp == "DEL_vs_DUP"){
        next
      }
      plot.suffix <- "cnv_split.pdf"
      plot.dims <- c(2.25, 2.45)
      parmar.y.mod <- 1
    }else{
      plot.suffix <- "pdf"
      plot.dims <- c(2.25, 2.3)
      parmar.y.mod <- 0
    }

    # Swarmplot of segment size
    pdf(paste(out.prefix.sub, "segment_size_distribs", plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2] + 0.3)
    segs.swarm(seg.df, x.bool=x.bool, y=log10(seg.df$size),
               x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
               boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
               add.y.axis=F, pt.cex=0.5, parmar=c(1.2, 4, 1.65 + parmar.y.mod, 0.2))
    axis(2, at=log10(logscale.minor), tck=-0.015, col=blueblack, labels=NA, lwd=0.7)
    axis(2, at=log10(logscale.demi), tck=-0.03, col=blueblack, labels=NA)
    axis(2, at=log10(logscale.demi.bp), tick=F, las=2, line=-0.65, labels=logscale.demi.bp.labels)
    mtext(2, line=2.75, text=bquote("log"[10] * "(Segment Size)"))
    dev.off()

    # Swarmplot of genes per segment
    pdf(paste(out.prefix.sub, "n_genes_per_segment", plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(seg.df, x.bool=x.bool, y=seg.df$n_genes,
               x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
               boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
               add.y.axis=T, ytitle="Genes in Segment", pt.cex=0.5,
               parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
    dev.off()

    # Swarmplot of gene density
    pdf(paste(out.prefix.sub, "basic_gene_density", plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(seg.df, x.bool=x.bool, y=100000*(seg.df$n_genes/seg.df$size),
               x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
               boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
               add.y.axis=T, ytitle="Genes per 100kb", pt.cex=0.5,
               parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
    dev.off()

    # Swarmplot of effect size
    pdf(paste(out.prefix.sub, "effect_size", plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(seg.df, x.bool=x.bool, y=log2(exp(seg.df$meta_best_lnor)),
               x.labs=x.labs, violin=T, add.pvalue=T,
               add.y.axis=T, pt.cex=0.5, cnv.split=cnv.split,
               boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
               ytitle=bquote("Max." ~  log[2]("Odds Ratio")), y.title.line=1.25,
               parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
    dev.off()

    # Swarmplot of peak P-value
    pdf(paste(out.prefix.sub, "peak_p_value", plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(seg.df, x.bool=x.bool, y=seg.df$meta_best_p,
               x.labs=x.labs, violin=T, add.pvalue=T,
               add.y.axis=T, pt.cex=0.5, cnv.split=cnv.split,
               boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
               ytitle=bquote("Max." ~ -log10(italic("P"))), y.title.line=1.25,
               parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
    dev.off()

    # Swarmplot of control CNV frequency
    if("pooled_control_freq" %in% colnames(seg.df)){
      pdf(paste(out.prefix.sub, "control_freq", plot.suffix, sep="."),
          height=plot.dims[1], width=plot.dims[2] + 0.15)
      segs.swarm(seg.df, x.bool=x.bool, y=log10(seg.df$pooled_control_freq),
                 x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
                 boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
                 add.y.axis=F, pt.cex=0.5,
                 parmar=c(1.2, 3, 1.65 + parmar.y.mod, 0.2))
      axis(2, at=log10(logscale.major), col=blueblack, tck=-0.03, labels=NA)
      sapply(1:length(logscale.major), function(x){
        axis(2, at=-log10(logscale.major[x]), tick=F, las=2, line=-0.65,
             labels=bquote(10^.(-log10(logscale.major[x]))))
      })
      mtext(2, text="CNV Freq. in Controls", line=2.1)
      dev.off()
    }

    # Swarmplot of average gene expression
    pdf(paste(out.prefix.sub, "avg_expression", plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2] + 0.15)
    segs.swarm(seg.df, x.bool, seg.df$gene_expression_harmonic_mean, add.y.axis=T,
               x.labs=x.labs, violin=T, add.pvalue=T,
               boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
               subset_to_regions=seg.df$region_id[which(seg.df$n_genes>0)],
               ytitle=bquote("Avg. Gene Expression"), y.title.line=1.75,
               pt.cex=0.5, cnv.split=cnv.split,
               parmar=c(1.2, 3, 1.65 + parmar.y.mod, 0.2))
    dev.off()

    # Swarmplot of proportion of ubiuqitously expressed genes
    pdf(paste(out.prefix.sub, "prop_ubi_expressed", plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2] + 0.15)
    segs.swarm(seg.df, x.bool, seg.df$prop_ubiquitously_expressed, add.y.axis=T,
               x.labs=x.labs, violin=T, add.pvalue=T,
               boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
               subset_to_regions=seg.df$region_id[which(seg.df$n_genes>0)],
               ytitle=bquote("Prop. Ubiquitously Expr."), y.title.line=1.75,
               pt.cex=0.5, cnv.split=cnv.split,
               parmar=c(1.2, 3, 1.65 + parmar.y.mod, 0.2))
    dev.off()

    # Gene set-based swarmplots
    for(gset in c("gnomAD_constrained", "clinical_LoF",
                  "clinical_GoF", "OMIM", "HPOmatched")){

      column <- paste("n", gset, "genes", sep="_")
      if(!(column %in% colnames(seg.df))){
        next
      }
      if(!any(seg.df[, column] > 0)){
        next
      }
      if(comp == "sig_vs_litGDs" & gset == "HPOmatched"){
        next
      }

      # Set gene set-specific parameters
      if(gset == "gnomAD_constrained"){
        glabel <- "Constrained"
      }else if(length(grep("clinical_", gset, fixed=T)) > 0){
        glabel <- paste(unlist(strsplit(gsub("clinical", "Known", gset), split="_")), collapse=" ")
      }else if(gset == "HPOmatched"){
        glabel <- "HPO-Matched"
      }else{
        glabel <- gset
      }

      # Swarmplot of number of genes
      pdf(paste(out.prefix.sub, column, plot.suffix, sep="."),
          height=plot.dims[1], width=plot.dims[2])
      segs.swarm(seg.df, x.bool=x.bool, y=seg.df[, column],
                 x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
                 boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
                 add.y.axis=T, ytitle=paste(glabel, "Genes"), pt.cex=0.5,
                 parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
      dev.off()

      # Swarmplot of gene density
      pdf(paste(out.prefix.sub, gsub("_genes$", "_gene_density", gsub("^n_", "", column)),
                plot.suffix, sep="."),
          height=plot.dims[1], width=plot.dims[2])
      segs.swarm(seg.df, x.bool=x.bool, y=1000000*(seg.df[, column]/seg.df$size),
                 subset_to_regions=seg.df$region_id[which(seg.df$n_genes>0)],
                 x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
                 boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
                 add.y.axis=T, ytitle=paste(glabel, "Genes / Mb"), pt.cex=0.5,
                 parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
      dev.off()
    }

    # Swarmplot of strongest constraint score
    for(constr in c("LOEUF", "MisOEUF")){
      column <- paste("min", constr, sep="_")
      pdf(paste(out.prefix.sub, column, plot.suffix, sep="."),
          height=plot.dims[1], width=plot.dims[2])
      segs.swarm(seg.df, x.bool=x.bool, y=seg.df[, column],
                 subset_to_regions=seg.df$region_id[which(seg.df$n_genes>0)],
                 x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
                 boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
                 add.y.axis=T, ytitle=bquote(min(.(constr))), pt.cex=0.5,
                 parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
      dev.off()
    }

    # Swarmplot of aggregated obs/exp in gnomAD for LoF and missense
    for(csq in c("LoF", "mis")){
      column <- paste("total", csq, "OE", sep="_")
      csq.name <- csq
      if(csq == "mis"){
        csq.name <- "Missense"
      }
      pdf(paste(out.prefix.sub, column, plot.suffix, sep="."),
          height=plot.dims[1], width=plot.dims[2])
      segs.swarm(seg.df, x.bool=x.bool, y=seg.df[, column],
                 subset_to_regions=seg.df$region_id[which(seg.df$n_genes>0)],
                 x.labs=x.labs, violin=T, add.pvalue=T, cnv.split=cnv.split,
                 boxplot.colors=boxplot.colors, boxplot.fill=boxplot.fill,
                 add.y.axis=T, ytitle=paste("Total", csq.name, "Obs/Exp"), pt.cex=0.5,
                 parmar=c(1.2, 2.75, 1.65 + parmar.y.mod, 0.2))
      dev.off()
    }
  }
}


# Scatterplots of size vs genes (log-log)
pdf(paste(out.prefix, "sig_plus_litGDs.genes_vs_size.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.all, x=log10(segs.all$size), y=log2(segs.all$n_genes), blue.bg=FALSE,
             subset_to_regions=segs.all$region_id[which(segs.all$n_genes>0)],
             x.at=log10(logscale.minor), x.labs=rep(NA, length(logscale.minor)),
             y.labs=NA, ytitle=bquote("log"[2] * "(Genes)"), y.title.line=1.5,
             pt.cex=0.5, parmar=c(2.75, 2.75, 0.3, 0.3))
axis(1, at=log10(c(200000, 1000000, 5000000)), tick=F,
     line=-0.7, labels=c("200kb", "1Mb", "5Mb"))
mtext(1, line=1.5, text=bquote("log"[10] * "(Segment Size)"))
axis(2, at=axTicks(2), tick=F, line=-0.65, labels=2^axTicks(2), las=2)
dev.off()
pdf(paste(out.prefix, "sig_only.genes_vs_size.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig, x=log10(segs.sig$size), y=log2(segs.sig$n_genes), blue.bg=FALSE,
             subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
             x.at=log10(logscale.minor), x.labs=rep(NA, length(logscale.minor)),
             y.labs=NA, ytitle=bquote("log"[2] * "(Genes)"), y.title.line=1.5,
             pt.cex=0.5, parmar=c(2.75, 2.75, 0.3, 0.3))
axis(1, at=log10(c(200000, 1000000, 5000000)), tick=F,
     line=-0.7, labels=c("200kb", "1Mb", "5Mb"))
mtext(1, line=1.5, text=bquote("log"[10] * "(Segment Size)"))
axis(2, at=axTicks(2), tick=F, line=-0.65, labels=2^axTicks(2), las=2)
dev.off()


# Scatterplot of peak P-value and corresponding lnOR for all sites
pdf(paste(out.prefix, "all_segs.best_p_vs_or.pdf", sep="."),
    height=2.9, width=3.3)
segs.best.l2or <- log2(exp(segs.all$meta_best_lnor))
x.max <- (7/6) * max(segs.best.l2or, na.rm=T)
segs.best.p <- segs.all$meta_best_p
segs.best.p[which(segs.best.p > 20)] <- 20
segs.scatter(segs.all, x=segs.best.l2or, y=segs.best.p,
             horiz.lines.at=c(gw.sig, -log10(0.05/95)), horiz.lines.lty=c(5, 2),
             horiz.lines.color=c(graphabs.green, blueblack), blue.bg=FALSE,
             xtitle="Max odds ratio, any phenotype",
             x.at=seq(0, x.max, 2), x.labs=c(2^seq(0, 6, 2), paste("2 ^", seq(8, x.max, 2))), parse.x.labs = T,
             ytitle=bquote("Max -log"[10] * (italic(P)) * ", any phenotype"),
             y.at=seq(0, 20, 4), y.labs=c(seq(0, 16, 4), expression(phantom(x) >= 20)), parse.y.labs=TRUE,
             x.title.line=1.4, y.title.line=1.9, xlims=c(0, x.max), ylims=c(0, 20),
             add.lm=F, pt.cex=0.75, parmar=c(2.4, 3, 0.4, 0.4))
x.bump <- 0.04 * (par("usr")[2] - par("usr")[1])
y.bump <- 0.03 * (par("usr")[4] - par("usr")[3])
text(x=par("usr")[2] + x.bump, y=-log10(0.05/95) + y.bump, labels=format.pval(0.05/95),
     cex=0.85, pos=2, col=blueblack, xpd=T)
text(x=par("usr")[2] + x.bump, y=gw.sig + y.bump, labels=format.pval(10^-gw.sig, nsmall=1),
     cex=0.85, pos=2, col=graphabs.green, xpd=T)
dev.off()
