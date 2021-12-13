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
perms <- load.perms.bygene(perm.res.in, subset_to_regions=analysis.ids)
lit.perms <- load.perms.bygene(lit.perm.res.in, subset_to_regions=analysis.ids)

# Prepare logfile for permutation tests
stats.df <- data.frame()

# Set global plot properties
pdf.dims.single <- c(2.2, 2.4)
parmar.single <- c(2.25, 2, 0, 2.25)
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
    pdf(paste(paste(outdir.sub, prefix, sep="/"), "HPOmatched_genes_permByGene",
              measure, subset, "pdf", sep="."),
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
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

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
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

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
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

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
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Average min(LOEUF) per Segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="min_LOEUF", measure="mean", invert=T,
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="min(LOEUF) per Segment",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Average min(MisOEUF) per Segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="min_MisOEUF", measure="mean", invert=T,
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="min(MisOEUF) per Segment",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Average total LoF O/E per Segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="total_LoF_OE", measure="mean", invert=T,
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Obs:Exp PTVs per Seg.",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Average total missense O/E per Segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="total_mis_OE", measure="mean", invert=T,
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Obs:Exp Missense per Seg.",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Average gene expression
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="gene_expression_harmonic_mean", measure="mean",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Average Gene Expression",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Proportion of ubiquitously expressed genes
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="prop_ubiquitously_expressed", measure="mean",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Prop. Ubiquitously Expressed",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Fraction of segments with at least one ubiquitously expressed gene
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="n_ubiquitously_expressed_genes", measure="frac.any",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="% Segs w/Ubiq. Expr. Gene",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Cleanup summary stats and append to existing log
  stats.df.tmp$segs <- subset
  stats.df <- rbind(stats.df, stats.df.tmp)
}


# Loop over a subset of segment partitions and generate DNM enrichment plots for each
# NOTE: only NDD segments used for these enrichments to match the phenotypes in DDD & ASC
for(subset in c("all_segs", "strong", "weak", "terminal",
                "interstitial", "NAHR", "nonrecurrent")){

  stats.df.tmp <- data.frame()

  # Set subset-specific properties
  cat(paste("\nEvaluating DNM enrichments for", subset))
  if(subset == "all_segs"){
    region.ids <- intersect(NDD.seg.ids, segs$region_id)
  }else if(subset == "strong"){
    region.ids <- intersect(NDD.seg.ids, lnor.groups[[3]])
  }else if(subset == "weak"){
    region.ids <- intersect(NDD.seg.ids, lnor.groups[[1]])
  }else if(subset == "terminal"){
    region.ids <- intersect(NDD.seg.ids,
                            segs.all$region_id[which(segs.all$terminal)])
  }else if(subset == "interstitial"){
    region.ids <- intersect(NDD.seg.ids,
                            segs.all$region_id[which(!segs.all$terminal)])
  }else if(subset == "NAHR"){
    region.ids <- intersect(NDD.seg.ids,
                            segs.all$region_id[which(segs.all$nahr)])
  }else if(subset == "nonrecurrent"){
    region.ids <- intersect(NDD.seg.ids,
                            segs.all$region_id[which(!segs.all$nahr)])
  }
  segs.sub <- segs[which(segs$region_id %in% region.ids), ]
  outdir.sub <- paste(outdir, paste(subset, "permByGene", sep="_"), sep="/")
  if(!dir.exists(outdir.sub)){
    dir.create(outdir.sub)
  }

  # Loop over each cohort & consequence
  for(cohort in c("ASC", "DDD", "ASC_unaffected", "DDD_plus_ASC")){
    for(ci in 1:length(csqs)){
      csq <- csqs[ci]
      csq.abbrev <- names(csqs)[ci]

      # Mean excess DNMs per gene
      new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                        subset_to_regions=region.ids,
                                        feature=paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_"),
                                        measure="mean",
                                        outdir.sub, paste(prefix, "NDD_only", sep="."),
                                        norm=F, norm.multi=F,
                                        n.bins.single=100, n.bins.multi=100, min.bins=100,
                                        x.title=bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev) ~ "/ Gene"),
                                        pdf.dims.single=pdf.dims.single,
                                        parmar.single=parmar.single,
                                        pdf.dims.multi=pdf.dims.multi,
                                        parmar.multi=parmar.multi)
      stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

      # Mean excess DNMs per gene after excluding significant genes from each study
      new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                        subset_to_regions=region.ids,
                                        feature=paste(cohort, "noSig_dnm", csq, "norm_excess_per_gene", sep="_"),
                                        measure="mean",
                                        outdir.sub, paste(prefix, "NDD_only", sep="."),
                                        norm=F, norm.multi=F,
                                        n.bins.single=100, n.bins.multi=100, min.bins=100,
                                        x.title=bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev) ~ "/ Gene"),
                                        pdf.dims.single=pdf.dims.single,
                                        parmar.single=parmar.single,
                                        pdf.dims.multi=pdf.dims.multi,
                                        parmar.multi=parmar.multi)
      stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)
    }
  }

  # Cleanup summary stats and append to existing log
  stats.df.tmp$segs <- subset
  stats.df <- rbind(stats.df, stats.df.tmp)
}


# Reorganize all main plots
sapply(c("all_segs", "strong", "weak", "developmental", "adult",
         "terminal", "interstitial", "NAHR", "nonrecurrent"),
       function(subset){
         outdir.sub <- paste(outdir, paste(subset, "permByGene", sep="_"), sep="/")
         reorganize.perm.plots(outdir.sub)
       })


# Write permutation stat results to logfile
write.table(stats.df[, c("feature", "measure", "sig", "segs", "CNV", "obs", "exp", "fold", "pval")],
            paste(outdir, "/", prefix, ".permByGene.stats.tsv", sep=""),
            col.names=T, row.names=F, sep="\t", quote=F)
