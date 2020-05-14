#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot panels involving effect size for large segments


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

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Plot effect size vs control frequency
pdf(paste(out.prefix, "lnOR_vs_control_freq.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=-log10(gw$pooled_control_freq), 
           y=gw$pooled_ln_or,
           xlims=c(2, max(-log10(gw$pooled_control_freq))),
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle=expression(-italic("log")[10] * "(Control Freq.)"),
           x.at=log10(logscale.major),
           x.labs=paste("10 ^", -log10(logscale.major)),
           parse.x.labs=T,
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs case frequency
pdf(paste(out.prefix, "lnOR_vs_case_freq.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=-log10(gw$pooled_case_freq), 
           y=gw$pooled_ln_or,
           xlims=c(2, max(-log10(gw$pooled_case_freq))),
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle=expression(-italic("log")[10] * "(Case Freq.)"),
           x.at=log10(logscale.major),
           x.labs=paste("10 ^", -log10(logscale.major)),
           parse.x.labs=T,
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs segment size
pdf(paste(out.prefix, "lnOR_vs_size.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=log10(gw$size), 
           y=gw$pooled_ln_or,
           xlims=c(5, 7),
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle=expression(italic("log")[10] * "(Size)"),
           x.at=log10(logscale.major.bp),
           x.labs.at=log10(logscale.major.bp),
           x.labs=logscale.major.bp.labels,
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs number of genes & gene density
pdf(paste(out.prefix, "lnOR_vs_genes.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=gw$n_genes, 
           y=gw$pooled_ln_or,
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle="Genes Overlapped",
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()
pdf(paste(out.prefix, "lnOR_vs_gene_density.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=100000 * gw$n_genes / gw$size, 
           y=gw$pooled_ln_or,
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle="Genes per 100kb",
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs number of genes & gene density
pdf(paste(out.prefix, "lnOR_vs_constrained_genes.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=gw$n_gnomAD_constrained_genes, 
           y=gw$pooled_ln_or,
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle="LoF Constrained Genes",
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()
pdf(paste(out.prefix, "lnOR_vs_constrained_gene_density.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=1000000 * gw$n_gnomAD_constrained_genes / gw$size, 
           y=gw$pooled_ln_or,
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle="Constr. Genes per 1Mb",
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs number of genes & gene density
pdf(paste(out.prefix, "lnOR_vs_OMIM_genes.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=gw$n_OMIM_genes, 
           y=gw$pooled_ln_or,
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle="OMIM Genes",
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()
pdf(paste(out.prefix, "lnOR_vs_OMIM_gene_density.pdf", sep="."),
    height=2.25, width=2.25)
gw.scatter(gw, 
           x=1000000 * gw$n_OMIM_genes / gw$size, 
           y=gw$pooled_ln_or,
           ylims=c(0, max(gw$pooled_ln_or)),
           xtitle="OMIM Genes per 1Mb",
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs mechanism
pdf(paste(out.prefix, "lnOR_vs_NAHR.pdf", sep="."),
    height=2.25, width=2.5)
gw.swarm(gw, 
         x.bool=gw$nahr, 
         y=gw$pooled_ln_or,
         ylims=c(0, max(gw$pooled_ln_or)),
         x.labs=c("Nonrecurrent", "NAHR"),
         xtitle="Mechanism",
         ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()
