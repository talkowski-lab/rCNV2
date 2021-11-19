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
  perms <- data.table::fread(perm.res.in, header=T, sep="\t", check.names=F)
  perms[, perm_idx := as.numeric(perms$perm_idx)]
  colnames(perms)[1] <- gsub("#", "", colnames(perms)[1], fixed=T)
  if(!is.null(subset_to_regions)){
    perms <- subset(perms, region_id %in% subset_to_regions)
  }

  # Add normalized columns
  perms[, gnomAD_constrained_prop := n_gnomAD_constrained_genes / n_genes]
  if("n_ubiquitously_expressed_genes" %in% colnames(perms)){
    perms[, prop_ubiquitously_expressed := n_ubiquitously_expressed_genes / n_genes]
  }

  # Normalize DNM excess per permutation
  normalize.dnms(perms, is.data.table=TRUE)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(MASS, quietly=T)
require(data.table, quietly=T)
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
# perm.res.in <- "~/scratch/rCNV2_analysis_d2.1000_permuted_segments_bygene.tsv.gz"
# lit.perm.res.in <- "~/scratch/rCNV2_analysis_d2.lit_GDs.1000_permuted_segments_bygene.tsv.gz"
# outdir <- "~/scratch"
# prefix <- "test_perm_bygene"

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]
analysis.ids <- segs$region_id

# Merge loci & segment data for genome-wide or FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs)

# Split seg IDs by gw-sig vs FDR vs. Bonf GD vs. nonsig GD
sig.ids <- segs$region_id[which(segs$any_sig)]
gw.ids <- segs$region_id[which(segs$gw_sig)]
fdr.ids <- segs$region_id[which(segs$fdr_sig)]
lit.ids <- segs$region_id[which(segs$any_gd & !segs$any_sig)]
all.lit.ids <- segs.all$region_id[which(segs.all$any_gd & !segs.all$any_sig)]
nonsig.lit.ids <- segs.all$region_id[which(segs.all$region_id %in% all.lit.ids & !segs.all$any_sig)]

# Get list of developmental loci (analysis set of segs only)
dev.seg.ids <- get.developmental.region_ids(loci, segs)
adult.seg.ids <- get.developmental.region_ids(loci, segs, adult=TRUE)
dev.segs <- segs[which(segs$region_id %in% dev.seg.ids), ]

# Get list of NDD loci (analysis set of segs only)
NDD.seg.ids <- get.ndd.region_ids(loci, segs)
NDD.segs <- segs[which(segs$region_id %in% NDD.seg.ids), ]

# Subset to loci either in top or bottom third of effect size distribution
lnor.groups <- split.regions.by.effect.size(segs, quantiles=3)

# Merge loci & segment data for significant sites only
segs.sig <- merge.loci.segs(loci, segs[which(segs$any_sig), ])

# Load permutation results
perms <- load.perms(perm.res.in, subset_to_regions=analysis.ids)
lit.perms <- load.perms(lit.perm.res.in, subset_to_regions=analysis.ids)

# Prepare logfile for permutation tests
stats.df <- data.frame()

# Set global plot properties
pdf.dims.single <- c(2.2, 2.4)
parmar.single <- c(2.25, 2, 0, 2)
pdf.dims.multi <- c(4, 3.5)
parmar.multi <- c(2.3, 6.05, 0, 1.5)
parmar.mod.frac.any <- c(0, 0, 0, 0)


# Generate permutation plots for HPO-matched OMIM genes among discovery segments
for(subset in c("all_sig", "gw_sig", "fdr_sig")){
  cat(paste("\nAnalyzing HPO-matched gene enrichment for", subset))
  if(subset == "all_sig"){
    region.ids <- sig.ids
    sig.label <- "GW_plus_FDR"
    diamond.pch <- 23
  }else if(subset == "gw_sig"){
    region.ids <- gw.ids
    sig.label <- "GW"
    diamond.pch <- 22
  }else if(subset == "fdr_sig"){
    region.ids <- fdr.ids
    sig.label <- "FDR"
    diamond.pch <- 23
  }

  # Prep output directory
  outdir.sub <- paste(outdir, "HPOmatched_genes", sep="/")
  if(!dir.exists(outdir.sub)){
    dir.create(outdir.sub)
  }

  # Single panel of permutation results
  for(measure in c("mean", "frac.any")){
    if(measure == "frac.any"){
      x.title <- "Pct. w/HPO-Matched Gene"
    }else{
      x.title <- "HPO-Matched Genes per Seg."
    }
    pdf(paste(paste(outdir.sub, prefix, sep="/"), "HPOmatched_genes", measure, subset, "pdf", sep="."),
        height=pdf.dims.single[1], width=pdf.dims.single[2])
    new.stats.df <- plot.seg.perms(segs, perms,
                                   feature="n_HPOmatched_genes", measure=measure,
                                   subset_to_regions=region.ids,
                                   n.bins=15, x.title=x.title,
                                   diamond.pch=diamond.pch, parmar=parmar.single)
    dev.off()
    new.stats.df$feature <- "n_HPOmatched_genes"
    new.stats.df$measure <- measure
    new.stats.df$sig <- sig.label
    new.stats.df$segs <- subset
    stats.df <- rbind(stats.df, new.stats.df)
  }
}



# Loop over all segment partitions and generate plots for each
for(subset in c("all_segs", "strong", "weak", "developmental", "adult",
                "terminal", "interstitial", "NAHR", "nonrecurrent")){

  stats.df.tmp <- data.frame()

  # Set subset-specific properties
  cat(paste("\nEvaluating other permutation tests for", subset))
  if(subset == "all_segs"){
    region.ids <- segs$region_id
  }else if(subset == "strong"){
    region.ids <- lnor.groups[[3]]
  }else if(subset == "weak"){
    region.ids <- lnor.groups[[1]]
  }else if(subset == "developmental"){
    region.ids <- dev.seg.ids
  }else if(subset == "adult"){
    region.ids <- adult.seg.ids
  }else if(subset == "terminal"){
    region.ids <- segs.all$region_id[which(segs.all$terminal)]
  }else if(subset == "interstitial"){
    region.ids <- segs.all$region_id[which(!segs.all$terminal)]
  }else if(subset == "NAHR"){
    region.ids <- segs.all$region_id[which(segs.all$nahr)]
  }else if(subset == "nonrecurrent"){
    region.ids <- segs.all$region_id[which(!segs.all$nahr)]
  }
  segs.sub <- segs[which(segs$region_id %in% region.ids), ]
  outdir.sub <- paste(outdir, paste(subset, "permByGene", sep="_"), sep="/")
  if(!dir.exists(outdir.sub)){
    dir.create(outdir.sub)
  }

  # Mean number of constrained genes per segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="n_gnomAD_constrained_genes", measure="mean",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Constrained Genes per Seg.",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  segs.df <- rbind(stats.df.tmp, new.stats.df)

  # Fraction of segments with at least one constrained gene
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="n_gnomAD_constrained_genes", measure="frac.any",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Pct. w/Constrained Gene",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  segs.df <- rbind(stats.df.tmp, new.stats.df)

  # Mean number of OMIM genes per segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="n_OMIM_genes", measure="mean",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="OMIM Genes per Seg.",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  segs.df <- rbind(stats.df.tmp, new.stats.df)

  # Fraction of segments with at least one OMIM gene
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="n_OMIM_genes", measure="frac.any",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Pct. w/OMIM Gene",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  segs.df <- rbind(stats.df.tmp, new.stats.df)

  ### TODO: FIX THIS TO ALLOW INVERSION OF X-AXIS
  ### THIS COULD BE ACCOMPLISHED BY JUST TAKING THE NEGATIVE OF THE FEATURE
  ### BUT FORCING THE X-AXIS FUNCTION TO TAKE THE ABSOLUTE VALUE
  # Average min(LOEUF) per Segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="min_LOEUF", measure="mean",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="min(LOEUF) per Segment",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  segs.df <- rbind(stats.df.tmp, new.stats.df)

  # Cleanup summary stats and append to existing log
  stats.df.tmp$segs <- subset
  stats.df <- rbind(stats.df, stats.df.tmp)
}

### TODO: REWORK ALL OF THE BELOW



### ANALYSES
# mean & sum of min_LOEUF
# mean & sum of min_MisOEUF
# mean & sum of total_LoF_OE
# mean & sum of total_mis_OE
# mean of gene_expression_harmonic_mean
# mean & frac.any of n_ubiquitously_expressed_genes
# All DNM enrichments with and without significant genes (NDD segments only)



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
                  x.title="Pct. w/Constrained Gene",
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.25, 0, 3.5))
print("Fraction of DEVELOPMENTAL segments with at least one constrained gene:")
plot.all.perm.res(segs, perms, lit.perms, subset_to_regions=dev.region_ids,
                  feature="n_gnomAD_constrained_genes", measure="frac.any",
                  outdir, paste(prefix, "dev_only", sep="."), norm=F, norm.multi=F,
                  n.bins.single=15, n.bins.multi=30,
                  x.title="Pct. w/Constrained Gene",
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.25, 0, 3.5))
print("Fraction of NDD segments with at least one constrained gene:")
plot.all.perm.res(segs, perms, lit.perms, subset_to_regions=NDD.region_ids,
                  feature="n_gnomAD_constrained_genes", measure="frac.any",
                  outdir, paste(prefix, "NDD_only", sep="."), norm=F, norm.multi=F,
                  n.bins.single=15, n.bins.multi=30,
                  x.title="Pct. w/Constrained Gene",
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
print("Mean constrained genes per DEVELOPMENTAL segment:")
plot.all.perm.res(segs, perms, lit.perms, subset_to_regions=dev.region_ids,
                  feature="n_gnomAD_constrained_genes", measure="mean",
                  outdir, paste(prefix, "dev_only", sep="."), norm=F, norm.multi=F,
                  n.bins.single=30, n.bins.multi=50,
                  x.title="Mean Constrained Genes",
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 1.2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.25, 0, 1.75))
print("Mean constrained genes per NDD segment:")
plot.all.perm.res(segs, perms, lit.perms, subset_to_regions=NDD.region_ids,
                  feature="n_gnomAD_constrained_genes", measure="mean",
                  outdir, paste(prefix, "NDD_only", sep="."), norm=F, norm.multi=F,
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
# Note: these analyses restricted to developmental & NDD loci only
sapply(c("ASC", "DDD", "ASC_unaffected"), function(cohort){
  sapply(1:length(csqs), function(ci){
    csq <- csqs[ci]
    csq.abbrev <- names(csqs)[ci]
    print(paste(cohort, csq.abbrev, "DNMs per gene vs. expected for DEVELOPMENTAL segments:"))
    plot.all.perm.res(segs, perms, lit.perms,
                      feature=paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_"),
                      measure="mean",
                      outdir, paste(prefix, "dev_only", sep="."),
                      subset_to_regions=dev.region_ids,
                      norm=F, norm.multi=F,
                      n.bins.single=100, n.bins.multi=100, min.bins=100,
                      x.title=bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev) ~ "/ Gene"),
                      pdf.dims.single=c(2.2, 2.4),
                      parmar.single=c(2.25, 2, 0, 1.2),
                      pdf.dims.multi=c(4, 3.5),
                      parmar.multi=c(2.25, 6.05, 0, 2.3))
    print(paste(cohort, csq.abbrev, "DNMs per gene vs. expected for NDD segments:"))
    plot.all.perm.res(segs, perms, lit.perms,
                      feature=paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_"),
                      measure="mean",
                      outdir, paste(prefix, "NDD_only", sep="."),
                      subset_to_regions=NDD.region_ids,
                      norm=F, norm.multi=F,
                      n.bins.single=100, n.bins.multi=100, min.bins=100,
                      x.title=bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev) ~ "/ Gene"),
                      pdf.dims.single=c(2.2, 2.4),
                      parmar.single=c(2.25, 2, 0, 1.2),
                      pdf.dims.multi=c(4, 3.5),
                      parmar.multi=c(2.25, 6.05, 0, 2.3))
  })
})


# Mean excess number of de novo PTVs & missense per gene in ASC & DDD
# AFTER removing exome-wide significant genes from their original publications
# Note: these analyses restricted to developmental & NDD loci only
sapply(c("ASC_noSig", "DDD_noSig", "ASC_unaffected_noSig"), function(cohort){
  sapply(1:length(csqs), function(ci){
    csq <- csqs[ci]
    csq.abbrev <- names(csqs)[ci]
    print(paste(cohort, csq.abbrev, "DNMs per gene vs. expected for DEVELOPMENTAL segments (significant genes REMOVED):"))
    plot.all.perm.res(segs, perms, lit.perms,
                      feature=paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_"),
                      measure="mean",
                      outdir, paste(prefix, "dev_only", sep="."),
                      subset_to_regions=dev.region_ids,
                      norm=F, norm.multi=F,
                      n.bins.single=100, n.bins.multi=100, min.bins=100,
                      x.title=bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev) ~ "/ Gene"),
                      pdf.dims.single=c(2.2, 2.4),
                      parmar.single=c(2.25, 2, 0, 1.2),
                      pdf.dims.multi=c(4, 3.5),
                      parmar.multi=c(2.25, 6.05, 0, 2.3))
    print(paste(cohort, csq.abbrev, "DNMs per gene vs. expected for NDD segments  (significant genes REMOVED):"))
    plot.all.perm.res(segs, perms, lit.perms,
                      feature=paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_"),
                      measure="mean",
                      outdir, paste(prefix, "NDD_only", sep="."),
                      subset_to_regions=NDD.region_ids,
                      norm=F, norm.multi=F,
                      n.bins.single=100, n.bins.multi=100, min.bins=100,
                      x.title=bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev) ~ "/ Gene"),
                      pdf.dims.single=c(2.2, 2.4),
                      parmar.single=c(2.25, 2, 0, 1.2),
                      pdf.dims.multi=c(4, 3.5),
                      parmar.multi=c(2.25, 6.05, 0, 2.3))
  })
})

