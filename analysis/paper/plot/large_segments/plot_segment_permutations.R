#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot large segment permutation results (bedtools shuffle)


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load gene list
load.genelist <- function(path){
  if(!is.null(path)){
    as.character(read.table(path, header=F)[, 1])
  }else{
    NULL
  }
}

# Load & parse permutation results
load.perms <- function(perm.res.in, subset_to_regions=NULL, constrained.genes=NULL){
  # Read full permutation table
  perms <- read.table(perm.res.in, header=T, sep="\t", check.names=F, comment.char="")
  perms$perm_idx <- as.numeric(perms$perm_idx)
  perm.range <- sort(unique(perms$perm_idx))
  
  # Post hoc annotation of constrained genes, if optioned
  if(!is.null(constrained.genes)){
    perms$n_gnomAD_constrained_genes <- as.numeric(sapply(perms$genes, function(gstr){
      if(is.na(gstr)){
        0
      }else{
        length(which(unlist(strsplit(gstr, split=";", fixed=T)) %in% constrained.genes))
      }
    }))
  }

  # Drop unnecessary columns
  cols.to.drop <- c("#chr", "start", "end", "coords", "size", "genes")
  perms <- perms[, -which(colnames(perms) %in% cols.to.drop)]
  if(!is.null(subset_to_regions)){
    perms <- perms[which(perms$region_id %in% subset_to_regions), ]
  }
  
  # Split each permutation into a list of dfs (one per perm)
  lapply(perm.range, function(i){
    as.data.frame(perms[which(perms$perm_idx==i), -(which(colnames(perms)=="perm_idx"))])
  })
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot scaled Venn diagram of GW-sig and literature-curated segs
cnv.venn <- function(segs, cnv, perms=NULL, margin=0.05){
  # Get plot values
  gw.idx <- which(segs$gw_sig & segs$cnv==cnv)
  n.gw <- length(gw.idx)
  lit.idx <- which(segs$any_gd & segs$cnv==cnv)
  n.lit <- length(lit.idx)
  both.idx <- intersect(gw.idx, lit.idx)
  n.both <- length(both.idx)
  
  # Set colors
  colors <- c(control.cnv.colors[which(names(cnv.colors)==cnv)],
              cnv.colors[which(names(cnv.colors)==cnv)])

  # # Calculate permuted P-value of overlap if permutation results are provided
  # if(!is.null(perms)){
  #   perm.dat <- perm.summary(perms, feature="any_gd", measure="sum", 
  #                            subset_to_regions=segs$region_id[gw.idx])
  #   perm.p <- length(which(perm.dat[, 1] >= n.both)) / nrow(perm.dat)
  #   if(perm.p == 0){
  #     p.fmt <- format.pval(1/nrow(perm.dat), equality="<")
  #   }else{
  #     p.fmt <- format.pval(perm.p)
  #   }
  # }else{
  #   p.fmt <- NULL
  # }
  # 
  # # Set axis label dimensions
  # if(cnv == "DEL"){
  #   lab.side <- 1
  #   lab.line <- -0.1
  # }else{
  #   lab.side <- 3
  #   lab.line <- -0.3
  # }
  
  # Prep plot area
  par(mar=rep(0.1, 4), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, 1), asp=1,
       xaxs="i", xaxt="n", xlab="", 
       yaxs="i", yaxt="n", ylab="")
  draw.pairwise.venn(n.lit, n.gw, n.both, 
                     col=colors, rotation.degree=180,
                     fill=sapply(colors, adjustcolor, alpha=0.5),
                     add=T, fontfamily="sans", margin=margin)
  # mtext(lab.side, line=lab.line, text=cnv, font=2, col=colors[2])
  # mtext(lab.side, line=lab.line, text=p.fmt)
  # axis(lab.side, at=c(0.05, 0.95), tck=0, col=blueblack, labels=NA)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(MASS, quietly=T)
require(VennDiagram, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
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
rcnv.config <- opts$`rcnv-config`
constrained.in <- opts$`constrained-genes`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# perm.res.in <- "~/scratch/rCNV2_analysis_d1.10000_permuted_segments.bed.gz"
# lit.perm.res.in <- "~/scratch/rCNV2_analysis_d1.lit_GDs.10000_permuted_segments.bed.gz"
# outdir <- "~/scratch/"
# prefix <- "test_seg_perm_res"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# constrained.in <- "~/scratch/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
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
segs.all <- load.segment.table(segs.in)

# Restrict to segments nominally significant in at least one phenotype
segs <- segs.all[which(segs.all$nom_sig), ]
nomsig.ids <- segs$region_id

# Split seg IDs by gw-sig vs GD (but not gw-sig)
gw.ids <- segs$region_id[which(segs$gw_sig)]
lit.ids <- segs$region_id[which(segs$any_gd & !segs$gw_sig)]
all.lit.ids <- segs.all$region_id[which(segs.all$any_gd & !segs.all$gw_sig)]

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Load permutation results
constrained.genes <- load.genelist(constrained.in)
perms <- load.perms(perm.res.in, subset_to_regions=nomsig.ids, 
                    constrained.genes=constrained.genes)
lit.perms <- load.perms(lit.perm.res.in, subset_to_regions=nomsig.ids,
                        constrained.genes=constrained.genes)

# Plot single panel of overlap with known genomic disorders
print("Overlap with known genomic disorders:")
pdf(paste(outdir, "/", prefix, ".gd_overlap.pdf", sep=""),
    height=2.2, width=2.5)
plot.seg.perms(segs, perms, feature="any_gd", 
               subset_to_regions=gw.ids,
               measure="sum", n.bins=100,
               x.title="Known Genomic Disorders",
               diamond.pch=22,
               parmar=c(2.2, 2, 0, 2.0))
dev.off()

# Venn diagrams of overlap between gw-sig and lit GDs
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(outdir, "/", prefix, ".gw_vs_lit.", cnv, ".venn.pdf", sep=""),
      height=1, width=1)
  cnv.venn(segs[which(segs$gw_sig | segs$any_gd), ], cnv, perms, margin=0.025)
  dev.off()
})

# Plot number of genes
print("Mean number of genes per segment:")
plot.all.perm.res(segs, perms, lit.perms, 
                  feature="n_genes", measure="mean",
                  outdir, prefix, norm=F, norm.multi=F,
                  n.bins.single=30, n.bins.multi=50,
                  x.title="Mean Genes per Segment", 
                  pdf.dims.single=c(2.2, 2.4),
                  parmar.single=c(2.25, 2, 0, 2.2),
                  pdf.dims.multi=c(4, 3.5),
                  parmar.multi=c(2.25, 6.05, 0, 2.2))

# Plot number of constrained genes
if(!is.null(constrained.in)){
  print("Fraction of segments with at least one constrained gene:")
  plot.all.perm.res(segs, perms, lit.perms,
                    feature="n_gnomAD_constrained_genes", measure="frac.any",
                    outdir, prefix, norm=F, norm.multi=F,
                    n.bins.single=15, n.bins.multi=30,
                    xmax=100, x.title="Pct. w/Constrained Gene",
                    pdf.dims.single=c(2.2, 2.4),
                    parmar.single=c(2.25, 2, 0, 3),
                    pdf.dims.multi=c(4, 3.5),
                    parmar.multi=c(2.25, 6.05, 0, 3.25))
  
}

# Plot single panel of enrichment for top P-values for non-significant GDs
print("Enrichment for top P-values among non-sig lit GDs:")
pdf(paste(outdir, "/", prefix, ".meta_best_p.nonsig_lit_gds.pdf", sep=""),
    height=2.3, width=2.3)
plot.seg.perms(segs.all[which(segs$region_id %in% all.lit.ids), ], 
               lit.perms, feature="meta_best_p", 
               subset_to_regions=all.lit.ids,
               measure="mean", n.bins=30, min.bins=10,
               x.title=bquote("Best" ~ -log[10](italic(P))),
               diamond.pch=21, x.title.line=1.4,
               parmar=c(2.4, 2, 0, 0.25))
dev.off()
