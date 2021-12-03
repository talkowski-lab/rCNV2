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
segs.all <- load.segment.table(segs.in)

# Add column of log2-scaled odds ratios for plotting convenience
segs.all$meta_best_log2_or <- log2(exp(segs.all$meta_best_lnor))

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]

# Merge loci & segment data for genome-wide/FDR significant sites only
segs.sig <- merge.loci.segs(loci, segs)

# Get list of developmental loci (analysis set of segs only)
dev.seg.ids <- get.developmental.region_ids(loci, segs)
adult.seg.ids <- get.developmental.region_ids(loci, segs, adult=TRUE)

# Set global parameters for plotting
scatter.dims <- c(2, 3)
scatter.parmar <- c(2.2, 2.5, 0.4, 5.1)

# Generate one set of plots for each subset of segments
sapply(c("all_segs", "developmental", "adult", "terminal",
         "interstitial", "NAHR", "nonrecurrent"),
       function(subset){

         # Set subset-specific properties
         if(subset == "all_segs"){
           region.ids <- segs$region_id
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

         # Create output subdirectory, if it doesn't exist
         out.subdir <- paste(dirname(out.prefix), subset, sep="/")
         if(!dir.exists(out.subdir)){
           dir.create(out.subdir)
         }
         out.prefix.sub.fname <- paste(basename(out.prefix), subset, sep=".")
         out.prefix.sub <- paste(out.subdir, out.prefix.sub.fname, sep="/")

         # Effect size vs. segment size
         pdf(paste(out.prefix.sub, "OR_vs_size.pdf", sep="."),
             height=scatter.dims[1], width=scatter.dims[2])
         segs.scatter(segs.sub, subset_to_regions=region.ids,
                      x=log10(segs.sub$size),
                      y=segs.sub$meta_best_log2_or,
                      pt.cex=0.6, xlims=c(5, 7), blue.bg=FALSE,
                      xtitle=expression(italic("log")[10] * "(Size)"),
                      x.at=log10(logscale.major.bp),
                      x.labs.at=log10(logscale.major.bp),
                      x.labs=logscale.major.bp.labels,
                      x.labs.line=-0.8,
                      ytitle=expression("Max" ~ log[2]("Odds Ratio")),
                      x.title.line=1.25, y.title.line=1.25,
                      parmar=scatter.parmar)
         dev.off()

         # Iterate over gene lists
         sapply(c("", "gnomAD_constrained", "OMIM"), function(gset){
           # Get generic values
           gfeat <- gsub("__", "_", paste("n", gset, "genes", sep="_"), fixed=T)
           ngene.xvals <- segs.sub[, gfeat]
           gdens.xvals <- 1000000 * ngene.xvals / segs.sub[, "size"]
           gdens.prefix <- gsub("^_", "", paste(gset, "gene_density", sep="_"))


           # Get gene set-specific values
           if(gset == ""){
             ngene.xlab <- "Genes per Segment"
             gdens.xlab <- "Genes per 1Mb"
           }else if(gset == "gnomAD_constrained"){
             ngene.xlab <- "Constrained Genes"
             gdens.xlab <- "Constrained per 1Mb"
           }else if(gset == "OMIM"){
             ngene.xlab <- "OMIM Genes"
             gdens.xlab <- "OMIM Genes per 1Mb"
           }

           # Effect size vs. number of genes in set
           pdf(paste(out.prefix.sub, ".OR_vs_", gfeat, ".pdf", sep=""),
               height=scatter.dims[1], width=scatter.dims[2])
           segs.scatter(segs.sub, subset_to_regions=region.ids,
                        x=ngene.xvals,
                        y=segs.sub$meta_best_log2_or,
                        pt.cex=0.6, blue.bg=FALSE,
                        xtitle=ngene.xlab,
                        x.labs.line=-0.8,
                        ytitle=expression("Max" ~ log[2]("Odds Ratio")),
                        x.title.line=1.25, y.title.line=1.25,
                        parmar=scatter.parmar)
           dev.off()

           # Effect size vs. gene density
           pdf(paste(out.prefix.sub, ".OR_vs_", gdens.prefix, ".pdf", sep=""),
               height=scatter.dims[1], width=scatter.dims[2])
           segs.scatter(segs.sub, subset_to_regions=region.ids,
                        x=gdens.xvals,
                        y=segs.sub$meta_best_log2_or,
                        pt.cex=0.6, blue.bg=FALSE,
                        xtitle=gdens.xlab,
                        x.labs.line=-0.8,
                        ytitle=expression("Max" ~ log[2]("Odds Ratio")),
                        x.title.line=1.25, y.title.line=1.25,
                        parmar=scatter.parmar)
           dev.off()
         })

         # Effect size vs. min(LOEUF) amd min(MisOEUF)
         sapply(c("min_LOEUF", "min_MisOEUF"), function(feature){
           pdf(paste(out.prefix.sub, ".OR_vs_", feature, ".pdf", sep=""),
               height=scatter.dims[1], width=scatter.dims[2])
           segs.scatter(segs.sub, subset_to_regions=region.ids,
                        x=segs.sub[, feature],
                        y=segs.sub$meta_best_log2_or,
                        pt.cex=0.6, blue.bg=FALSE,
                        xtitle=paste("min(", unlist(strsplit(feature, split="_"))[2], ")", sep=""),
                        x.labs.line=-0.8,
                        ytitle=expression("Max" ~ log[2]("Odds Ratio")),
                        x.title.line=1.25, y.title.line=1.25,
                        parmar=scatter.parmar)
           dev.off()
         })

         # Effect size vs. total LoF O/E
         pdf(paste(out.prefix.sub, ".OR_vs_total_LoF_OE.pdf", sep=""),
             height=scatter.dims[1], width=scatter.dims[2])
         segs.scatter(segs.sub, subset_to_regions=region.ids,
                      x=segs.sub$total_LoF_OE,
                      y=segs.sub$meta_best_log2_or,
                      pt.cex=0.6, blue.bg=FALSE,
                      xtitle="Obs:Exp PTVs per Seg.",
                      x.labs.line=-0.8, max.x.ticks=5,
                      ytitle=expression("Max" ~ log[2]("Odds Ratio")),
                      x.title.line=1.25, y.title.line=1.25,
                      parmar=scatter.parmar)
         dev.off()

         # Effect size vs. total LoF O/E
         pdf(paste(out.prefix.sub, ".OR_vs_total_mis_OE.pdf", sep=""),
             height=scatter.dims[1], width=scatter.dims[2])
         segs.scatter(segs.sub, subset_to_regions=region.ids,
                      x=segs.sub$total_mis_OE,
                      y=segs.sub$meta_best_log2_or,
                      pt.cex=0.6, blue.bg=FALSE,
                      xtitle="Obs:Exp Mis. per Seg.",
                      x.labs.line=-0.8, max.x.ticks=5,
                      ytitle=expression("Max" ~ log[2]("Odds Ratio")),
                      x.title.line=1.25, y.title.line=1.25,
                      parmar=scatter.parmar)
         dev.off()
       })

