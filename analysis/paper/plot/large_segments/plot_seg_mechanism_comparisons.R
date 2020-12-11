#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot comparisons of NAHR vs. nonrecurrent large segments for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


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
# out.prefix <- "~/scratch/seg_mechanism_test"
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

# Get list of neuro loci
neuro.region_ids <- get.neuro.region_ids(loci, segs)

# Swarmplot of segment size vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.size.pdf", sep="."),
    height=2.25, width=2.7)
segs.swarm(segs, x.bool=segs$nahr, y=log10(segs$size), 
           x.labs=c("Non-NAHR", "NAHR"), violin=T, alternative="two.sided",
           add.pvalue=T, add.y.axis=F, pt.cex=0.75, parmar=c(1.2, 4, 2.5, 0))
axis(2, at=log10(logscale.minor), tck=-0.015, col=blueblack, labels=NA, lwd=0.7)
axis(2, at=log10(logscale.demi), tck=-0.03, col=blueblack, labels=NA)
axis(2, at=log10(logscale.demi.bp), tick=F, las=2, line=-0.65, labels=logscale.demi.bp.labels)
mtext(2, line=2.75, text=bquote("log"[10] * "(Segment Size)"))
dev.off()
cat(paste("Two-sided Wilcox test of NAHR vs non-NAHR size:",
          format(wilcox.test(segs$size ~ segs$nahr)$p.value, scientific=T),
          "\n"))

# Swarmplot of genes vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.n_genes.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs, x.bool=segs$nahr, y=segs$n_genes, 
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Genes in Segment", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of gene density vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.basic_gene_density.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs, x.bool=segs$nahr, y=100000*(segs$n_genes/segs$size), 
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Genes per 100kb", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of constrained gene density vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.constrained_gene_density.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs, x.bool=segs$nahr, y=1000000*(segs$n_gnomAD_constrained_genes/segs$size), 
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Constr. Genes per 1Mb", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of OMIM gene density vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.omim_gene_density.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs, x.bool=segs$nahr, y=1000000*(segs$n_OMIM_genes/segs$size), 
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="OMIM Genes per 1Mb", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of DDD dnLoFs vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.DDD_dnLoF_per_gene.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs[which(segs$n_genes>0), ], 
           x.bool=segs$nahr[which(segs$n_genes>0)], 
           y=segs$DDD_dnm_lof_norm_excess_per_gene,
           subset_to_regions=neuro.region_ids,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.y.axis=T, add.pvalue=T,
           ytitle=bquote("Excess" ~ italic("dn") * "PTV" ~ "/ Gene"), 
           pt.cex=0.75, parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of DDD dnMis vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.DDD_dnMis_per_gene.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs[which(segs$n_genes>0), ], 
           x.bool=segs$nahr[which(segs$n_genes>0)], 
           y=segs$DDD_dnm_mis_norm_excess_per_gene,
           subset_to_regions=neuro.region_ids,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.y.axis=T, add.pvalue=T,
           ytitle=bquote("Excess" ~ italic("dn") * "Mis" ~ "/ Gene"), 
           pt.cex=0.75, parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of ASC dnLoFs vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.ASC_dnLoF_per_gene.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs[which(segs$n_genes>0), ], 
           x.bool=segs$nahr[which(segs$n_genes>0)], 
           y=segs$ASC_dnm_lof_norm_excess_per_gene,
           subset_to_regions=neuro.region_ids,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.y.axis=T, add.pvalue=T,
           ytitle=bquote("Excess" ~ italic("dn") * "PTV" ~ "/ Gene"), 
           pt.cex=0.75, parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of ASC dnMis vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.ASC_dnMis_per_gene.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs[which(segs$n_genes>0), ], 
           x.bool=segs$nahr[which(segs$n_genes>0)], 
           y=segs$ASC_dnm_mis_norm_excess_per_gene,
           subset_to_regions=neuro.region_ids,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.y.axis=T, add.pvalue=T,
           ytitle=bquote("Excess" ~ italic("dn") * "Mis" ~ "/ Gene"), 
           pt.cex=0.75, parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of average gene expression vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.avg_gene_expression.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs, 
           x.bool=segs$nahr, 
           y=segs$gene_expression_harmonic_mean, 
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.y.axis=T, add.pvalue=T,
           ytitle=bquote("Avg. Gene Expression"), 
           pt.cex=0.75, parmar=c(1.2, 3, 2.5, 0))
dev.off()
