#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot panels involving pleiotropy for large segments
# Note: does not use HPO clustering-based definition of pleiotropy
# Instead, uses simpler definition of ≥2 associated HPOs


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(MASS, quietly=T)

# List of command-line options
option_list <- list()

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

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# out.prefix <- "~/scratch/test_pleiotropy"

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]

# Merge loci & segment data for genome-wide/FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs)

# Get list of gw-sig loci
gw.region_ids <- segs$region_id[which(segs$gw_sig)]

# Set global swarmplot values
parmar <- c(2.2, 2.7, 1.65, 0.2)
x.bool <- segs.sig$n_hpos > 1
x.labs <- c("One", "Multiple")
x.title <- "Associated Phenos."

# Generate two swarmplots plots per comparison
# one for all CNV types and one split by CNV type
for(cnv.split in c(TRUE, FALSE)){
  if(cnv.split){
    plot.suffix <- "cnv_split.pdf"
    plot.dims <- c(2.25, 2.45)
    parmar.y.mod <- 1
  }else{
    plot.suffix <- "pdf"
    plot.dims <- c(2.25, 2.3)
    parmar.y.mod <- 0
  }

  # Swarmplot of segment size
  pdf(paste(out.prefix, "segment_size_distribs", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2] + 0.3)
  segs.swarm(segs.sig, x.bool=x.bool, y=log10(segs.sig$size),
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T, cnv.split=cnv.split,
             add.y.axis=F, pt.cex=0.5, parmar=parmar + c(0, 1.1, parmar.y.mod, 0))
  axis(2, at=log10(logscale.minor), tck=-0.015, col=blueblack, labels=NA, lwd=0.7)
  axis(2, at=log10(logscale.demi), tck=-0.03, col=blueblack, labels=NA)
  axis(2, at=log10(logscale.demi.bp), tick=F, las=2, line=-0.65, labels=logscale.demi.bp.labels)
  mtext(2, line=2.75, text=bquote("log"[10] * "(Segment Size)"))
  dev.off()

  # Swarmplot of genes per segment
  pdf(paste(out.prefix, "n_genes_per_segment", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2])
  segs.swarm(segs.sig, x.bool=x.bool, y=segs.sig$n_genes,
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T, cnv.split=cnv.split,
             add.y.axis=T, ytitle="Genes in Segment", pt.cex=0.5,
             parmar=parmar + c(0, 0, parmar.y.mod, 0))
  dev.off()

  # Swarmplot of gene density
  pdf(paste(out.prefix, "basic_gene_density", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2])
  segs.swarm(segs.sig, x.bool=x.bool, y=100000*(segs.sig$n_genes/segs.sig$size),
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T, cnv.split=cnv.split,
             add.y.axis=T, ytitle="Genes per 100kb", pt.cex=0.5,
             parmar=parmar + c(0, 0, parmar.y.mod, 0))
  dev.off()

  # Swarmplot of effect size
  pdf(paste(out.prefix, "effect_size", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2])
  segs.swarm(segs.sig, x.bool=x.bool, y=log2(exp(segs.sig$meta_best_lnor)),
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T,
             add.y.axis=T, pt.cex=0.5, cnv.split=cnv.split,
             ytitle=bquote("Max." ~  log[2]("Odds Ratio")), y.title.line=1.25,
             parmar=parmar + c(0, 0, parmar.y.mod, 0))
  dev.off()

  # Swarmplot of average gene expression
  pdf(paste(out.prefix, "avg_expression", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2] + 0.15)
  segs.swarm(segs.sig, x.bool, segs.sig$gene_expression_harmonic_mean, add.y.axis=T,
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T,
             subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
             ytitle=bquote("Avg. Gene Expression"), y.title.line=1.75,
             pt.cex=0.5, cnv.split=cnv.split,
             parmar=parmar + c(0, 0, parmar.y.mod, 0))
  dev.off()

  # Swarmplot of proportion of ubiuqitously expressed genes
  pdf(paste(out.prefix, "prop_ubi_expressed", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2] + 0.15)
  segs.swarm(segs.sig, x.bool, segs.sig$prop_ubiquitously_expressed, add.y.axis=T,
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T,
             subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
             ytitle=bquote("Prop. Ubiquitously Expr."), y.title.line=1.75,
             pt.cex=0.5, cnv.split=cnv.split,
             parmar=parmar + c(0, 0, parmar.y.mod, 0))
  dev.off()

  # Swarmplot of ubiuqitously expressed genes per segment
  pdf(paste(out.prefix, "n_ubi_expressed", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2] + 0.15)
  segs.swarm(segs.sig, x.bool, segs.sig$n_ubiquitously_expressed, add.y.axis=T,
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T,
             subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
             ytitle=bquote("Ubiquitously Expr. Genes"), y.title.line=1.75,
             pt.cex=0.5, cnv.split=cnv.split,
             parmar=parmar + c(0, 0, parmar.y.mod, 0))
  dev.off()

  # Swarmplot of density of ubiuqitously expressed genes
  pdf(paste(out.prefix, "ubi_expressed_gene_density", plot.suffix, sep="."),
      height=plot.dims[1], width=plot.dims[2] + 0.15)
  segs.swarm(segs.sig, x.bool, 1000000 * segs.sig$n_ubiquitously_expressed / segs.sig$size, add.y.axis=T,
             x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T,
             subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
             ytitle=bquote("Ubi. Expr. Genes / Mb"), y.title.line=1.75,
             pt.cex=0.5, cnv.split=cnv.split,
             parmar=parmar + c(0, 0, parmar.y.mod, 0))
  dev.off()

  # Gene set-based swarmplots
  for(gset in c("gnomAD_constrained", "clinical_LoF",
                "clinical_GoF", "OMIM")){

    column <- paste("n", gset, "genes", sep="_")
    if(!(column %in% colnames(segs.sig))){
      next
    }
    if(!any(segs.sig[, column] > 0)){
      next
    }

    # Set gene set-specific parameters
    if(gset == "gnomAD_constrained"){
      glabel <- "Constrained"
    }else if(length(grep("clinical_", gset, fixed=T)) > 0){
      glabel <- paste(unlist(strsplit(gsub("clinical", "Known", gset), split="_")), collapse=" ")
    }else{
      glabel <- gset
    }

    # Swarmplot of number of genes
    pdf(paste(out.prefix, column, plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(segs.sig, x.bool=x.bool, y=segs.sig[, column],
               x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T, cnv.split=cnv.split,
               add.y.axis=T, ytitle=paste(glabel, "Genes"), pt.cex=0.5,
               parmar=parmar + c(0, 0, parmar.y.mod, 0))
    dev.off()

    # Swarmplot of gene density
    pdf(paste(out.prefix, gsub("_genes$", "_gene_density", gsub("^n_", "", column)),
              plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(segs.sig, x.bool=x.bool, y=1000000*(segs.sig[, column]/segs.sig$size),
               subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
               x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T, cnv.split=cnv.split,
               add.y.axis=T, ytitle=paste(glabel, "Genes / Mb"), pt.cex=0.5,
               parmar=parmar + c(0, 0, parmar.y.mod, 0))
    dev.off()
  }

  # Swarmplot of strongest constraint score
  for(constr in c("LOEUF", "MisOEUF")){
    column <- paste("min", constr, sep="_")
    pdf(paste(out.prefix, column, plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(segs.sig, x.bool=x.bool, y=segs.sig[, column],
               subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
               x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T, cnv.split=cnv.split,
               add.y.axis=T, ytitle=paste("min(", constr, ")", sep=""), pt.cex=0.5,
               parmar=parmar + c(0, 0, parmar.y.mod, 0))
    dev.off()
  }

  # Swarmplot of aggregated obs/exp in gnomAD for LoF and missense
  for(csq in c("LoF", "mis")){
    column <- paste("total", csq, "OE", sep="_")
    csq.name <- csq
    if(csq == "mis"){
      csq.name <- "Missense"
    }
    pdf(paste(out.prefix, column, plot.suffix, sep="."),
        height=plot.dims[1], width=plot.dims[2])
    segs.swarm(segs.sig, x.bool=x.bool, y=segs.sig[, column],
               subset_to_regions=segs.sig$region_id[which(segs.sig$n_genes>0)],
               x.labs=x.labs, xtitle=x.title, violin=T, add.pvalue=T, cnv.split=cnv.split,
               add.y.axis=T, ytitle=paste("Total", csq.name, "Obs/Exp"), pt.cex=0.5,
               parmar=parmar + c(0, 0, parmar.y.mod, 0))
    dev.off()
  }
}


# Set global scatterplot values
scatter.dims <- c(2, 3)
scatter.parmar <- c(2.15, 2.3, 0.4, 5.1)

# Scatterplots of # HPOs vs. various features
sapply(c("size", "n_genes", "n_gnomAD_constrained_genes", "n_OMIM_genes",
         "min_LOEUF", "min_MisOEUF", "total_LoF_OE", "total_mis_OE",
         "gene_expression_harmonic_mean", "n_ubiquitously_expressed_genes",
         "meta_best_lnor", "n_clinical_LoF_genes", "n_clinical_GoF_genes",
         "gene_density", "constrained_gene_density", "OMIM_gene_density",
         "ubiquitously_expressed_gene_density", "clinical_LoF_gene_density",
         "clinical_GoF_gene_density"),
       function(feature){

  # Set global default parameters
  if(feature %in% colnames(segs.sig)){
    xvals <- segs.sig[, feature]
  }
  x.at <- NULL
  x.labs.at <- NULL
  x.labs  <- NULL
  x.title.line <- 1.1
  x.lims <- NULL
  parmar.mod <- rep(0, 4)
  dims.mod <- rep(0, 2)

  # Set feature-specific parameters
  if(feature == "size"){
    xvals <- log10(segs.sig[, feature])
    x.at <- log10(logscale.major.bp)
    x.labs.at <- log10(logscale.major.bp)
    x.labs <- logscale.major.bp.labels
    x.title <- bquote(log[10]("Segment Size"))
    x.title.line <- 1.3
    x.lims <- c(5, 7)
    parmar.mod <- c(0.1, 0, 0, 0)
  }else if(feature == "n_genes"){
    x.title <- "Genes per Segment"
  }else if(feature == "n_gnomAD_constrained_genes"){
    x.title <- "Constrained Genes"
  }else if(feature == "n_OMIM_genes"){
    x.title <- "OMIM Genes"
  }else if(feature == "min_LOEUF"){
    x.title <- "min(LOEUF)"
  }else if(feature == "min_MisOEUF"){
    x.title <- "min(MisOEUF)"
  }else if(feature == "total_LoF_OE"){
    x.title <- "Obs:Exp PTVs per Seg."
  }else if(feature == "total_mis_OE"){
    x.title <- "Obs:Exp Mis. per Seg."
  }else if(feature == "gene_expression_harmonic_mean"){
    x.title <- "Avg. Gene Expression"
  }else if(feature == "n_ubiquitously_expressed_genes"){
    x.title <- "Ubi. Expressed Genes"
  }else if(feature == "meta_best_lnor"){
    xvals <- log2(exp(segs.sig$meta_best_lnor))
    x.title <- bquote(log[2]("Odds Ratio"))
    x.title.line  <- 1.3
  }else if(feature == "n_clinical_LoF_genes"){
    x.title <- "Known HI Genes"
  }else if(feature == "n_clinical_GoF_genes"){
    x.title <- "Known GoF Genes"
  }else if(feature == "gene_density"){
    xvals <- 100000 * segs.sig$n_genes / segs.sig$size
    x.title <- "Genes per 100kb"
  }else if(feature == "constrained_gene_density"){
    xvals <- 1000000 * segs.sig$n_gnomAD_constrained_genes / segs.sig$size
    x.title <- "Constrained per Mb"
  }else if(feature == "OMIM_gene_density"){
    xvals <- 1000000 * segs.sig$n_OMIM_genes / segs.sig$size
    x.title <- "OMIM Genes per Mb"
  }else if(feature == "ubiquitously_expressed_gene_density"){
    xvals <- 1000000 * segs.sig$n_ubiquitously_expressed_genes / segs.sig$size
    x.title <- "Ubi. Expr. per Mb"
  }else if(feature == "clinical_LoF_gene_density"){
    xvals <- 1000000 * segs.sig$n_clinical_LoF_genes / segs.sig$size
    x.title <- "Known HI Genes per Mb"
  }else if(feature == "clinical_GoF_gene_density"){
    xvals <- 1000000 * segs.sig$n_clinical_GoF_genes / segs.sig$size
    x.title <- "GoF Genes per Mb"
  }

  # Generate scatterplot
  pdf(paste(out.prefix, ".nHPOs_vs_", feature, ".pdf", sep=""),
      height=scatter.dims[1] + dims.mod[1],
      width=scatter.dims[2] + dims.mod[2])
  segs.scatter(segs.sig,
               x=xvals,
               y=segs.sig$n_hpos,
               pt.cex=0.6, blue.bg=FALSE,
               xlims=x.lims,
               x.at=x.at,
               x.labs.at=x.labs.at,
               x.labs=x.labs,
               x.labs.line=-0.8,
               xtitle=x.title,
               ytitle="Associated HPOs",
               x.title.line=x.title.line, y.title.line=1.4,
               parmar=scatter.parmar + parmar.mod)
  dev.off()
})

