#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot large segment permutation results (gene set-based test)


options(stringsAsFactors=F, scipen=1000)
csqs <- c("PTVs" = "lof", "Mis." = "mis")


######################
### DATA FUNCTIONS ###
######################
# Load & parse permutation results
load.perms <- function(perm.res.in, subset_to_regions=NULL){
  # Read full permutation table
  perms <- read.table(perm.res.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(perms)[1] <- gsub("#", "", colnames(perms)[1], fixed=T)
  perm.range <- sort(unique(perms$perm_idx))
  if(!is.null(subset_to_regions)){
    perms <- perms[which(perms$region_id %in% subset_to_regions), ]
  }
  
  # Add normalized columns
  perms$gnomAD_constrained_prop <- perms$n_gnomAD_constrained_genes / perms$n_genes
  if("n_ubiquitously_expressed_genes" %in% colnames(perms)){
    perms$prop_ubiquitously_expressed <- perms$n_ubiquitously_expressed_genes / perms$n_genes
  }
  
  # Split each permutation into a list of dfs (one per perm)
  # Also performs per-permutation normalization of DNM excess
  lapply(perm.range, function(i){
    normalize.dnms(as.data.frame(perms[which(perms$perm_idx==i), -(which(colnames(perms)=="perm_idx"))]))
  })
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
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv perm_res.bed lit_GD_perm_res.bed outdir prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Six positional arguments required: loci.bed, segs.tsv, perm_res.bed, lit_GD_perm_res.bed, outdir, and prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
perm.res.in <- args$args[3]
lit.perm.res.in <- args$args[4]
outdir <- args$args[5]
prefix <- args$args[6]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# perm.res.in <- "~/scratch/rCNV2_analysis_d1.10000_permuted_segments_bygene.tsv.gz"
# lit.perm.res.in <- "~/scratch/rCNV2_analysis_d1.lit_GDs.10000_permuted_segments_bygene.tsv.gz"
# outdir <- "~/scratch"
# prefix <- "test_perm_bygene"
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

# Subset to segments nominally significant for at least one phenotype
segs <- segs[which(segs$nom_sig), ]
nomsig.ids <- segs$region_id

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Get list of neuro loci
neuro.region_ids <- get.neuro.region_ids(loci, segs)

# Load permutation results
perms <- load.perms(perm.res.in, subset_to_regions=nomsig.ids)
lit.perms <- load.perms(lit.perm.res.in, subset_to_regions=nomsig.ids)

# Fraction of gw-sig segments with at least one HPO-matched gene
print("Fraction of segments with at least one HPO-matched gene:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.frac_any.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(gw, perms, feature="n_HPOmatched_genes", measure="frac.any", 
               n.bins=15, x.title="Pct. w/HPO-Matched Gene", 
               diamond.pch=22, parmar=c(2.2, 2, 0, 1.8))
dev.off()

# Mean # of HPO-matched genes per gw-sig segment
print("Mean HPO-matched genes per segment:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.mean.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(gw, perms, feature="n_HPOmatched_genes", measure="mean", 
               n.bins=30, x.title="Mean HPO-Matched Genes",
               diamond.pch=22, parmar=c(2.2, 2, 0, 2))
dev.off()

# Fraction of segments with at least one constrained gene
print("Fraction of segments with at least one constrained gene:")
plot.all.perm.res(segs, perms, lit.perms, 
                  feature="n_gnomAD_constrained_genes", measure="frac.any",
                  outdir, prefix, norm=F, norm.multi=F,
                  n.bins.single=15, n.bins.multi=30,
                  xmax=100, x.title="Pct. w/Constrained Gene",
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.25, 0, 3.5))

# Mean number of constrained genes per segment
print("Mean constrained genes per segment:")
plot.all.perm.res(segs, perms, lit.perms,
                  feature="n_gnomAD_constrained_genes", measure="mean",
                  outdir, prefix, norm=F, norm.multi=F,
                  n.bins.single=30, n.bins.multi=50,
                  x.title="Mean Constrained Genes",
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 1.2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.25, 0, 1.75))

# Average gene expression per segment
print("Harmonic mean of gene expression per segment:")
plot.all.perm.res(segs, perms, lit.perms, 
                  feature="gene_expression_harmonic_mean", measure="mean",
                  outdir, prefix, norm=F, norm.multi=F,
                  n.bins.single=30, n.bins.multi=50,
                  x.title="Avg. Gene Expression", 
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 1.2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.25, 0, 0.5))

# Mean excess number of de novo PTVs & missense per gene in ASC & DDD
# When restricting gw-sig to neuro-associated loci
sapply(c("ASC", "DDD", "ASC_unaffected"), function(cohort){
  sapply(1:length(csqs), function(ci){
    csq <- csqs[ci]
    csq.abbrev <- names(csqs)[ci]
    print(paste(cohort, csq.abbrev, "DNMs per gene vs. expected:"))
    plot.all.perm.res(segs, perms, lit.perms, 
                      feature=paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_"),
                      measure="mean",
                      outdir, paste(prefix, "gw_neuro_plus_lit", sep="."), 
                      subset_to_regions=neuro.region_ids,
                      norm=F, norm.multi=F,
                      n.bins.single=100, n.bins.multi=100, min.bins=100,
                      x.title=bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev) ~ "/ Gene"), 
                      pdf.dims.single=c(2.2, 2.4),
                      parmar.single=c(2.25, 2, 0, 1.2),
                      pdf.dims.multi=c(4, 3.5),
                      parmar.multi=c(2.25, 6.05, 0, 2.3))
  })
})

