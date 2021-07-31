#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
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
require(MASS, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list()

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

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# perm.res.in <- "~/scratch/rCNV2_analysis_d2.10000_permuted_segments_bygene.tsv.gz"
# lit.perm.res.in <- "~/scratch/rCNV2_analysis_d2.lit_GDs.10000_permuted_segments_bygene.tsv.gz"
# outdir <- "~/scratch"
# prefix <- "test_perm_bygene"

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Subset to segments nominally significant for at least one phenotype
segs <- segs[which(segs$nom_sig), ]
nomsig.ids <- segs$region_id

# Merge loci & segment data for genome-wide or FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs)

# Split seg IDs by gw-sig vs FDR vs. lit GD
sig.ids <- segs$region_id[which(segs$any_sig)]
gw.ids <- segs$region_id[which(segs$gw_sig)]
fdr.ids <- segs$region_id[which(segs$fdr_sig)]
lit.ids <- segs$region_id[which(segs$any_gd & !segs$any_sig)]
all.lit.ids <- segs.all$region_id[which(segs.all$any_gd & !segs.all$any_sig)]

# Get list of developmental and NDD loci
dev.region_ids <- get.developmental.region_ids(loci, segs)
NDD.region_ids <- get.ndd.region_ids(loci, segs)

# Load permutation results
perms <- load.perms(perm.res.in, subset_to_regions=nomsig.ids)
lit.perms <- load.perms(lit.perm.res.in, subset_to_regions=nomsig.ids)


# Fraction of discovered segments with at least one HPO-matched gene
print("Fraction of all significant segments with at least one HPO-matched gene:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.frac_any.all_sig.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(segs.sig, perms, feature="n_HPOmatched_genes", measure="frac.any",
               n.bins=15, x.title="Pct. w/HPO-Matched Gene",
               diamond.pch=23, parmar=c(2.2, 2, 0, 1.8))
dev.off()
print("Fraction of GW significant segments with at least one HPO-matched gene:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.frac_any.gw_sig.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(segs.sig, perms, feature="n_HPOmatched_genes", measure="frac.any",
               subset_to_regions=gw.ids, n.bins=15, x.title="Pct. w/HPO-Matched Gene",
               diamond.pch=22, parmar=c(2.2, 2, 0, 1.8))
dev.off()
print("Fraction of FDR significant segments with at least one HPO-matched gene:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.frac_any.FDR_sig.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(segs.sig, perms, feature="n_HPOmatched_genes", measure="frac.any",
               subset_to_regions=fdr.ids, n.bins=15, x.title="Pct. w/HPO-Matched Gene",
               diamond.pch=23, parmar=c(2.2, 2, 0, 1.8))
dev.off()


# Mean # of HPO-matched genes per significant segment
print("Mean HPO-matched genes per significant segment:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.mean.all_sig.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(segs.sig, perms, feature="n_HPOmatched_genes", measure="mean",
               n.bins=30, x.title="Mean HPO-Matched Genes",
               diamond.pch=23, parmar=c(2.2, 2, 0, 2))
dev.off()
print("Mean HPO-matched genes per GW sig segment:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.mean.gw_sig.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(segs.sig, perms, feature="n_HPOmatched_genes", measure="mean",
               subset_to_regions=gw.ids,
               n.bins=30, x.title="Mean HPO-Matched Genes",
               diamond.pch=22, parmar=c(2.2, 2, 0, 2))
dev.off()
print("Mean HPO-matched genes per FDR sig segment:")
pdf(paste(outdir, "/", prefix, ".HPOmatched_genes.mean.fdr_sig.pdf", sep=""),
    height=2.2, width=2.4)
plot.seg.perms(segs.sig, perms, feature="n_HPOmatched_genes", measure="mean",
               subset_to_regions=fdr.ids,
               n.bins=30, x.title="Mean HPO-Matched Genes",
               diamond.pch=23, parmar=c(2.2, 2, 0, 2))
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
print("Fraction of DEVELOPMENTAL segments with at least one constrained gene:")
plot.all.perm.res(segs, perms, lit.perms, subset_to_regions=dev.region_ids,
                  feature="n_gnomAD_constrained_genes", measure="frac.any",
                  outdir, paste(prefix, "dev_only", sep="."), norm=F, norm.multi=F,
                  n.bins.single=15, n.bins.multi=30,
                  xmax=100, x.title="Pct. w/Constrained Gene",
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.25, 0, 3.5))
print("Fraction of NDD segments with at least one constrained gene:")
plot.all.perm.res(segs, perms, lit.perms, subset_to_regions=NDD.region_ids,
                  feature="n_gnomAD_constrained_genes", measure="frac.any",
                  outdir, paste(prefix, "NDD_only", sep="."), norm=F, norm.multi=F,
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

