#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot large segment permutation results (bedtools shuffle)


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load & parse permutation results
load.perms <- function(perm.res.in, subset_to_regions=NULL, constrained.genes=NULL){
  # Read full permutation table
  perms <- data.table::fread(perm.res.in, header=T, sep="\t", check.names=F)
  perms[, perm_idx := as.numeric(perms$perm_idx)]

  # Post hoc annotation of constrained genes, if optioned
  posthoc.constr <- function(gstr){
    if(is.na(gstr)){as.integer(0)}else{
      length(which(unlist(strsplit(gstr, split=";", fixed=T)) %in% constrained.genes))
    }
  }
  perms[, n_gnomAD_constrained_genes := posthoc.constr(genes), by=seq_len(nrow(perms))]

  # Drop unnecessary columns
  cols.to.drop <- c("#chr", "start", "end", "coords", "size", "genes")
  perms[, (cols.to.drop) := NULL]
  if(!is.null(subset_to_regions)){
    perms <- perms[region_id %in% subset_to_regions]
  }

  return(perms)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot scaled Venn diagram of GW-sig and literature-curated segs
cnv.venn <- function(segs, cnv, sig.labels, perms=NULL, margin=0.05){
  # Get plot values
  gw.idx <- which(sig.labels & segs$cnv==cnv)
  n.gw <- length(gw.idx)
  lit.idx <- which(segs$any_gd & segs$cnv==cnv)
  n.lit <- length(lit.idx)
  both.idx <- intersect(gw.idx, lit.idx)
  n.both <- length(both.idx)

  # Set colors
  colors <- c(control.cnv.colors[which(names(cnv.colors)==cnv)],
              cnv.colors[which(names(cnv.colors)==cnv)])

  # Calculate permuted P-value of overlap if permutation results are provided
  if(!is.null(perms)){
    perm.dat <- perm.summary(perms, feature="any_gd", measure="sum",
                             subset_to_regions=segs$region_id[gw.idx])
    perm.p <- length(which(perm.dat[, 1] >= n.both)) / nrow(perm.dat)
    if(perm.p == 0){
      p.fmt <- format.pval(1/nrow(perm.dat), equality="<")
    }else{
      p.fmt <- format.pval(perm.p)
    }
  }else{
    p.fmt <- NULL
  }

  # Set axis label dimensions
  sig.line <- -0.75
  if(cnv == "DEL"){
    lab.side <- 1
    sig.side <- 3
  }else{
    lab.side <- 3
    sig.side <- 1
  }

  # Prep plot area
  par(mar=rep(0.1, 4), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, 1), asp=1,
       xaxs="i", xaxt="n", xlab="",
       yaxs="i", yaxt="n", ylab="")
  draw.pairwise.venn(n.lit, n.gw, n.both,
                     col=colors, rotation.degree=180,
                     fill=sapply(colors, adjustcolor, alpha=0.5),
                     add=T, fontfamily="sans", margin=margin)
  mtext(lab.side, line=-1-margin, text=cnv, font=2, col=colors[2], xpd=T)
  mtext(sig.side, line=sig.line-margin, text=p.fmt)
  axis(sig.side, at=c(0.25, 0.75), tck=0.025, col=blueblack, labels=NA, line=sig.line-margin)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(MASS, quietly=T)
require(VennDiagram, quietly=T)
require(data.table, quietly=T)
require(rCNV2, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--constrained-genes"), help="List of constrained genes to be annotated.")
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
constrained.in <- opts$`constrained-genes`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# perm.res.in <- "~/scratch/rCNV2_analysis_d2.1000_permuted_segments.bed.gz"
# lit.perm.res.in <- "~/scratch/rCNV2_analysis_d2.lit_GDs.1000_permuted_segments.bed.gz"
# outdir <- "~/scratch/"
# prefix <- "test_seg_perm_res"
# constrained.in <- "~/scratch/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]
analysis.ids <- segs$region_id

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
constrained.genes <- load.genelist(constrained.in)
perms <- load.perms(perm.res.in, constrained.genes=constrained.genes)
lit.perms <- load.perms(lit.perm.res.in, constrained.genes=constrained.genes)

# Prepare logfile for permutation tests
stats.df <- data.frame()

# Set global plot properties
pdf.dims.single <- c(2.2, 2.4)
parmar.single <- c(2.25, 2, 0, 2)
pdf.dims.multi <- c(4, 3.5)
parmar.multi <- c(2.3, 6.05, 0, 1.5)
parmar.mod.frac.any <- c(0, 0, 0, 0)


# Generate permutation plots & Venns for overlap with genomic disorders
for(subset in c("all_sig", "gw_sig", "fdr_sig")){
  cat(paste("\nAnalyzing GD overlap for", subset))
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
  outdir.sub <- paste(outdir, "GD_overlap", sep="/")
  if(!dir.exists(outdir.sub)){
    dir.create(outdir.sub)
  }

  # Single panel of permutation results
  pdf(paste(paste(outdir.sub, prefix, sep="/"), "gd_overlap", subset, "pdf", sep="."),
      height=pdf.dims.single[1], width=pdf.dims.single[2])
  new.stats.df <- plot.seg.perms(segs, perms, feature="any_gd",
                 subset_to_regions=region.ids,
                 measure="sum", n.bins=100,
                 x.title="Known Genomic Disorders",
                 diamond.pch=diamond.pch,
                 parmar=parmar.single)
  dev.off()
  new.stats.df$feature <- "any_gd"
  new.stats.df$measure <- "sum"
  new.stats.df$sig <- sig.label
  new.stats.df$segs <- subset
  stats.df <- rbind(stats.df, new.stats.df)

  # Venn diagrams of overlap between discovery and lit GDs
  sapply(c("DEL", "DUP"), function(cnv){
    pdf(paste(paste(outdir.sub, prefix, sep="/"),
              "gd_overlap", subset, cnv, "venn.pdf", sep="."),
        height=1.5, width=1.5)
    cnv.venn(segs.all, cnv, segs.all$region_id %in% region.ids,
             perms, margin=0.05)
    dev.off()
  })
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
  outdir.sub <- paste(outdir, paste(subset, "permBySize", sep="_"), sep="/")
  if(!dir.exists(outdir.sub)){
    dir.create(outdir.sub)
  }

  # Mean number of genes per segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="n_genes", measure="mean",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Genes per Segment",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  segs.df <- rbind(stats.df.tmp, new.stats.df)

  # Fraction of segments with at least one gene
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="n_genes", measure="frac.any",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="% Overlapping Any Gene",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single + parmar.mod.frac.any,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi + parmar.mod.frac.any)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  if(!is.null(constrained.in)){
    # Average number of constrained genes
    new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                      subset_to_regions=region.ids,
                                      feature="n_gnomAD_constrained_genes", measure="mean",
                                      outdir.sub, prefix, norm=F, norm.multi=F,
                                      n.bins.single=15, n.bins.multi=30,
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
                                      n.bins.single=15, n.bins.multi=30,
                                      x.title="Pct. w/Constrained Gene",
                                      pdf.dims.single=pdf.dims.single,
                                      parmar.single=parmar.single + parmar.mod.frac.any,
                                      pdf.dims.multi=pdf.dims.multi,
                                      parmar.multi=parmar.multi + parmar.mod.frac.any)
    stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)
  }

  # Mean number of BCA breakpoints per segment
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="Redin_BCAs", measure="mean",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="BCA Breakpoints per Seg.",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Fraction of segments with at least one BCA
  new.stats.df <- plot.all.perm.res(segs.sub, perms, lit.perms,
                                    subset_to_regions=region.ids,
                                    feature="Redin_BCAs", measure="frac.any",
                                    outdir.sub, prefix, norm=F, norm.multi=F,
                                    n.bins.single=30, n.bins.multi=50,
                                    x.title="Pct. w/BCA Breakpoint",
                                    pdf.dims.single=pdf.dims.single,
                                    parmar.single=parmar.single + parmar.mod.frac.any,
                                    pdf.dims.multi=pdf.dims.multi,
                                    parmar.multi=parmar.multi + parmar.mod.frac.any)
  stats.df.tmp <- rbind(stats.df.tmp, new.stats.df)

  # Cleanup summary stats and append to existing log
  stats.df.tmp$segs <- subset
  stats.df <- rbind(stats.df, stats.df.tmp)
}


# Enrichment for top P-values for non-significant GDs
pdf(paste(outdir, "/", prefix, ".meta_best_p.nonsig_lit_gds.pdf", sep=""),
    height=pdf.dims.single[1], width=pdf.dims.single[2])
new.stats.df <- plot.seg.perms(segs.all[which(segs.all$region_id %in% nonsig.lit.ids), ],
               lit.perms, feature="meta_best_p",
               subset_to_regions=nonsig.lit.ids,
               measure="mean", n.bins=30, min.bins=10,
               x.title=bquote("Best" ~ -log[10](italic(P))),
               diamond.pch=21, x.title.line=1.4,
               parmar=parmar.single)
dev.off()
new.stats.df$feature <- "meta_best_P"
new.stats.df$measure <- "mean"
new.stats.df$sig <- "nonsig_litGDs"
new.stats.df$segs <- "nonsig_litGDs"
stats.df <- rbind(stats.df, stats.df.tmp)

# Write permutation stat results to logfile
write.table(stats.df[, c("feature", "measure", "sig", "segs", "CNV", "obs", "exp", "fold", "pval")],
            paste(outdir, "/", prefix, ".permBySize.stats.tsv", sep=""),
            col.names=T, row.names=F, sep="\t", quote=F)
