#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot panels involving pleiotropy for large segments
# Note: does not use HPO clustering-based definition of pleiotropy
# Instead, uses simpler definition of ≥2 associated HPOs


options(stringsAsFactors=F, scipen=1000)


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
  stop(paste("Three positional arguments required: loci.bed, segs.tsv, output_prefix\n", sep=" "))
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

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Restrict to segments nominally significant in at least one phenotype
segs <- segs[which(segs$nom_sig), ]

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Set global plotting values
parmar <- c(2.3, 3.0, 1.5, 0.5)

# Plot pleiotropy vs region size
pdf(paste(out.prefix, "pleiotropy_vs_size.pdf", sep="."),
    height=2.25, width=2)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=log10(gw$size),
         cnv.split=F,
         violin=T,
         add.pvalue=T,
         alternative="less",
         xtitle="Associated Phenos.",
         x.labs=c("One", "Multiple"), 
         y.at=log10(logscale.minor), 
         y.labs=logscale.demi.bp.labels, 
         y.labs.at=log10(logscale.demi.bp),
         ytitle=expression(italic("log")[10] * "(Size)"),
         parmar=parmar)
dev.off()

# Plot pleiotropy vs. number of genes & gene density
pdf(paste(out.prefix, "pleiotropy_vs_genes.pdf", sep="."),
    height=2.25, width=2)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=gw$n_genes,
         cnv.split=F,
         violin=T,
         add.pvalue=T,
         alternative="less",
         xtitle="Associated Phenos.",
         x.labs=c("One", "Multiple"), 
         ytitle="Genes Overlapped",
         parmar=parmar)
dev.off()
# Note: this plot is formatted differently because of placement in main figure
pdf(paste(out.prefix, "pleiotropy_vs_gene_density.pdf", sep="."),
    height=2.3, width=1.8)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=100000 * gw$n_genes / gw$size,
         cnv.split=F,
         violin=T,
         pt.cex=0.85,
         add.pvalue=T,
         alternative="less",
         xtitle="",
         x.labs=c("One", "Multiple"), 
         ytitle="Genes per 100kb",
         y.at=seq(0, 8, 2),
         y.title.line=1.25,
         parmar=c(2.3, 2.3, 1.5, 0.5))
dev.off()

# Plot pleiotropy vs. constrained genes & gene density
pdf(paste(out.prefix, "pleiotropy_vs_constrained_genes.pdf", sep="."),
    height=2.25, width=2)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=gw$n_gnomAD_constrained_genes,
         cnv.split=F,
         violin=T,
         add.pvalue=T,
         alternative="less",
         xtitle="Associated Phenos.",
         x.labs=c("One", "Multiple"), 
         ytitle="Constrained Genes",
         parmar=parmar)
dev.off()
# Note: this plot is formatted differently because of placement in main figure
pdf(paste(out.prefix, "pleiotropy_vs_constrained_gene_density.pdf", sep="."),
    height=2.3, width=1.8)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=1000000 * gw$n_gnomAD_constrained_genes / gw$size,
         cnv.split=F,
         violin=T,
         pt.cex=0.85,
         add.pvalue=T,
         alternative="less",
         xtitle="",
         x.labs=c("One", "Multiple"), 
         ytitle="Genes per 1Mb",
         y.title.line=1.35,
         parmar=c(2.3, 2.3, 1.5, 0.5))
dev.off()

# Plot pleiotropy vs. OMIM genes & gene density
pdf(paste(out.prefix, "pleiotropy_vs_OMIM_genes.pdf", sep="."),
    height=2.25, width=2)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=gw$n_OMIM_genes,
         cnv.split=F,
         violin=T,
         add.pvalue=T,
         alternative="less",
         xtitle="Associated Phenos.",
         x.labs=c("One", "Multiple"), 
         ytitle="OMIM Genes",
         parmar=parmar)
dev.off()
pdf(paste(out.prefix, "pleiotropy_vs_OMIM_gene_density.pdf", sep="."),
    height=2.25, width=2)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=1000000 * gw$n_OMIM_genes / gw$size,
         cnv.split=F,
         violin=T,
         add.pvalue=T,
         alternative="less",
         xtitle="Associated Phenos.",
         x.labs=c("One", "Multiple"), 
         ytitle="OMIM Genes per 1Mb",
         parmar=parmar)
dev.off()

# Plot pleiotropy vs. mechanism
pdf(paste(out.prefix, "pleiotropy_vs_mechanism.pdf", sep="."),
    height=2.25, width=2.5)
gw.swarm(gw, x.bool=gw$nahr, y=gw$n_hpos, add.pvalue=T, cnv.split=F,
         x.labs=c("Nonrecurrent", "NAHR"), violin=T, add.y.axis=T,
         ytitle="Associated HPOs", pt.cex=0.75,
         parmar=c(1.2, 3, 1.5, 0))
dev.off()

# Plot pleiotropy vs. average gene expression
pdf(paste(out.prefix, "pleiotropy_vs_gene_expression.pdf", sep="."),
    height=2.25, width=2)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=gw$gene_expression_harmonic_mean,
         cnv.split=F,
         violin=T,
         add.pvalue=T,
         alternative="two.sided",
         xtitle="Associated Phenos.",
         x.labs=c("One", "Multiple"), 
         ytitle="Avg. Gene Expression",
         parmar=parmar)
dev.off()

# Plot pleiotropy vs. proportion of ubiquitously expressed genes
pdf(paste(out.prefix, "pleiotropy_vs_prop_ubi_expressed.pdf", sep="."),
    height=2.25, width=2)
gw.swarm(gw,
         x.bool=gw$n_hpos > 1,
         y=gw$prop_ubiquitously_expressed,
         cnv.split=F,
         violin=T,
         add.pvalue=T,
         alternative="two.sided",
         xtitle="Associated Phenos.",
         x.labs=c("One", "Multiple"), 
         ytitle="Prop. Ubiquitously Expr.",
         parmar=parmar)
dev.off()

