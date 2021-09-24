#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot panels involving effect size for large segments


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
# out.prefix <- "~/scratch/test_effect_sizes"

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Restrict to segments nominally significant in at least one phenotype
segs.all <- segs[which(segs$any_gd | segs$any_sig), ]
segs <- segs.all[which(segs.all$nom_sig), ]

# Merge loci & segment data for genome-wide/FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs[which(segs$any_sig), ])

# Get list of developmental and gw-sig loci
dev.region_ids <- get.developmental.region_ids(loci, segs)
gw.region_ids <- segs$region_id[which(segs$gw_sig)]

# Plot effect size vs control frequency
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_control_freq.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
           x=-log10(segs.sig$pooled_control_freq),
           y=segs.sig$pooled_ln_or,
           pt.cex=0.6,
           xlims=c(2, max(-log10(segs.sig$pooled_control_freq))),
           ylims=c(0, max(segs.sig$pooled_ln_or)),
           xtitle=expression(-italic("log")[10] * "(Control Freq.)"),
           x.at=log10(logscale.major),
           x.labs=paste("10 ^", -log10(logscale.major)),
           parse.x.labs=T,
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs case frequency
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_case_freq.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
           x=-log10(segs.sig$pooled_case_freq),
           y=segs.sig$pooled_ln_or,
           pt.cex=0.6,
           xlims=c(2, max(-log10(segs.sig$pooled_case_freq))),
           ylims=c(0, max(segs.sig$pooled_ln_or)),
           xtitle=expression(-italic("log")[10] * "(Case Freq.)"),
           x.at=log10(logscale.major),
           x.labs=paste("10 ^", -log10(logscale.major)),
           parse.x.labs=T,
           ytitle=expression("Effect Size" ~ (italic("ln") * " OR")))
dev.off()

# Plot effect size vs segment size
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_size.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
           x=log10(segs.sig$size),
           y=segs.sig$max_ln_or,
           pt.cex=0.6,
           xlims=c(5, 7),
           ylims=c(0, max(segs.sig$max_ln_or)),
           xtitle=expression(italic("log")[10] * "(Size)"),
           x.at=log10(logscale.major.bp),
           x.labs.at=log10(logscale.major.bp),
           x.labs=logscale.major.bp.labels,
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           x.title.line=1.5,
           y.title.line=1.25,
           parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_size.dev_only.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
             x=log10(segs.sig$size),
             y=segs.sig$max_ln_or,
             subset_to_regions=dev.region_ids,
             pt.cex=0.6,
             xlims=c(5, 7),
             ylims=c(0, max(segs.sig$max_ln_or)),
             xtitle=expression(italic("log")[10] * "(Size)"),
             x.at=log10(logscale.major.bp),
             x.labs.at=log10(logscale.major.bp),
             x.labs=logscale.major.bp.labels,
             ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
             x.title.line=1.5,
             y.title.line=1.25,
             parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()

# Plot effect size vs number of genes & gene density
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_genes.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
           x=segs.sig$n_genes,
           y=segs.sig$max_ln_or,
           pt.cex=0.6,
           ylims=c(0, max(segs.sig$max_ln_or)),
           xtitle="Genes Overlapped",
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           x.title.line=1.5,
           y.title.line=1.25,
           parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_gene_density.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
           x=100000 * segs.sig$n_genes / segs.sig$size,
           y=segs.sig$max_ln_or,
           pt.cex=0.6,
           ylims=c(0, max(segs.sig$max_ln_or)),
           xtitle="Genes per 100kb",
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           x.title.line=1.5,
           y.title.line=1.25,
           parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_gene_density.dev_only.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
             x=100000 * segs.sig$n_genes / segs.sig$size,
             y=segs.sig$max_ln_or,
             pt.cex=0.6,
             subset_to_regions=dev.region_ids,
             ylims=c(0, max(segs.sig$max_ln_or)),
             xtitle="Genes per 100kb",
             ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
             x.title.line=1.5,
             y.title.line=1.25,
             parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()

# Plot effect size vs number of constrained genes & gene density
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_constrained_genes.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
             x=segs.sig$n_gnomAD_constrained_genes,
             y=segs.sig$max_ln_or,
             pt.cex=0.6,
             ylims=c(0, max(segs.sig$max_ln_or)),
             xtitle="LoF Constrained Genes",
             ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
             x.title.line=1.5,
             y.title.line=1.25,
             parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_constrained_genes.dev_only.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
             x=segs.sig$n_gnomAD_constrained_genes,
             y=segs.sig$max_ln_or,
             subset_to_regions=dev.region_ids,
             pt.cex=0.6,
             ylims=c(0, max(segs.sig$max_ln_or)),
             xtitle="LoF Constrained Genes",
             ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
             x.title.line=1.5,
             y.title.line=1.25,
             parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()
pdf(paste(out.prefix, "gw_plus_fdr.lnOR_vs_any_constrained.pdf", sep="."),
    height=2.25, width=2)
segs.swarm(segs.sig,
           x.bool=segs.sig$n_gnomAD_constrained_genes>0,
           y=segs.sig$max_ln_or,
           pt.cex=0.4,
           add.pvalue=TRUE,
           violin=TRUE,
           xtitle="Constrained Genes",
           x.labs=c("0", "1+"),
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           y.title.line=1.15,
           parmar=c(2.3, 2.5, 2.5, 0.5))
dev.off()
pdf(paste(out.prefix, "gw_plus_fdr.lnOR_vs_any_constrained.dev_only.pdf", sep="."),
    height=2.25, width=2)
segs.swarm(segs.sig,
           x.bool=segs.sig$n_gnomAD_constrained_genes>0,
           y=segs.sig$max_ln_or,
           subset_to_regions=dev.region_ids,
           pt.cex=0.4,
           add.pvalue=TRUE,
           violin=TRUE,
           xtitle="Constrained Genes",
           x.labs=c("0", "1+"),
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           y.title.line=1.15,
           parmar=c(2.3, 2.5, 2.5, 0.5))
dev.off()
pdf(paste(out.prefix, "gw_only.lnOR_vs_any_constrained.pdf", sep="."),
    height=2.25, width=2)
segs.swarm(segs.sig,
           x.bool=segs.sig$n_gnomAD_constrained_genes>0,
           y=segs.sig$max_ln_or,
           subset_to_regions=gw.region_ids,
           pt.cex=0.5,
           add.pvalue=TRUE,
           violin=TRUE,
           xtitle="Constrained Genes",
           x.labs=c("0", "1+"),
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           y.title.line=1.15,
           parmar=c(2.3, 2.5, 2.5, 0.5))
dev.off()
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_constrained_gene_density.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
             x=100000 * segs.sig$n_gnomAD_constrained_genes / segs.sig$size,
             y=segs.sig$max_ln_or,
             pt.cex=0.6,
             ylims=c(0, max(segs.sig$max_ln_or)),
             xtitle="Constrained per 100kb",
             ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
             x.title.line=1.5,
             y.title.line=1.25,
             parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_constrained_gene_density.dev_only.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
             x=100000 * segs.sig$n_gnomAD_constrained_genes / segs.sig$size,
             y=segs.sig$max_ln_or,
             pt.cex=0.6,
             subset_to_regions=dev.region_ids,
             ylims=c(0, max(segs.sig$max_ln_or)),
             xtitle="Constrained per 100kb",
             ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
             x.title.line=1.5,
             y.title.line=1.25,
             parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()

# Plot effect size vs number of OMIM genes & gene density
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_OMIM_genes.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
           x=segs.sig$n_OMIM_genes,
           y=segs.sig$max_ln_or,
           pt.cex=0.6,
           ylims=c(0, max(segs.sig$max_ln_or)),
           xtitle="OMIM Genes",
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           x.title.line=1.5,
           y.title.line=1.25,
           parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_OMIM_gene_density.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs.sig,
           x=1000000 * segs.sig$n_OMIM_genes / segs.sig$size,
           y=segs.sig$max_ln_or,
           pt.cex=0.6,
           ylims=c(0, max(segs.sig$max_ln_or)),
           xtitle="OMIM Genes per 1Mb",
           ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
           x.title.line=1.5,
           y.title.line=1.25,
           parmar=c(2.5, 2.5, 0.8, 0.8))
dev.off()

# Plot effect size vs mechanism
pdf(paste(out.prefix, "gw_plus_FDR.lnOR_vs_NAHR.pdf", sep="."),
    height=2.25, width=2.35)
segs.swarm(segs.sig,
         x.bool=segs.sig$nahr,
         y=segs.sig$max_ln_or,
         pt.cex=0.4,
         add.pvalue=TRUE,
         violin=TRUE,
         ylims=c(0, max(segs.sig$max_ln_or)),
         x.labs=c("Nonrecurrent", "NAHR"),
         xtitle="Mechanism",
         ytitle=expression("Max Effect Size" ~ (italic("ln") * " OR")),
         y.title.line=1.15,
         parmar=c(2.3, 2.5, 2.5, 0.1))
dev.off()
