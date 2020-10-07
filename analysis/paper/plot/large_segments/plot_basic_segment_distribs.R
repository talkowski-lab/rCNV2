#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot basic distributions for large segments for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
  make_option(c("--gw-sig"), default=10e-6, help="P-value cutoff for genome-wide signficance.")
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
gw.sig <- -log10(opts$`gw-sig`)

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# out.prefix <- "~/scratch/final_segs"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# gw.sig <- -log10(3.715428E-6)
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

# Subset to all GD segs & those with nominal significance
segs.all <- segs[which(segs$any_gd | segs$gw_sig), ]
segs <- segs.all[which(segs.all$nom_sig), ]

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Swarmplot of segment size
pdf(paste(out.prefix, "gw_plus_litGDGs.segment_size_distribs.pdf", sep="."),
    height=2.25, width=2.75)
segs.swarm(segs, x.bool=!(segs$gw_sig), y=log10(segs$size), 
           x.labs=c("GW-Sig.", "Lit. (Not sig.)"), violin=T, add.pvalue=T,
           add.y.axis=F, pt.cex=0.75, parmar=c(1.2, 4, 2.65, 0.2))
axis(2, at=log10(logscale.minor), tck=-0.015, col=blueblack, labels=NA, lwd=0.7)
axis(2, at=log10(logscale.demi), tck=-0.03, col=blueblack, labels=NA)
axis(2, at=log10(logscale.demi.bp), tick=F, las=2, line=-0.65, labels=logscale.demi.bp.labels)
mtext(2, line=2.75, text=bquote("log"[10] * "(Segment Size)"))
dev.off()

# Swarmplot of genes per segment
pdf(paste(out.prefix, "gw_plus_litGDGs.n_genes_per_segment.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs, x.bool=!(segs$gw_sig), y=segs$n_genes, 
           x.labs=c("GW-Sig.", "Lit. (Not sig.)"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Genes in Segment", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0.2))
dev.off()

# Swarmplot of gene density
pdf(paste(out.prefix, "gw_plus_litGDGs.basic_gene_density.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(segs, x.bool=!(segs$gw_sig), y=100000*(segs$n_genes/segs$size), 
           x.labs=c("GW-Sig.", "Lit. (Not sig.)"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Genes per 100kb", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0.2))
dev.off()

# Scatterplot of size vs genes (log-log)
pdf(paste(out.prefix, "gw_plus_litGDs.genes_vs_size.pdf", sep="."),
    height=2.25, width=2.25)
segs.scatter(segs, x=log10(segs$size), y=log2(segs$n_genes), 
             subset_to_regions=segs$region_id[which(segs$n_genes>0)],
             x.at=log10(logscale.minor), x.labs=rep(NA, length(logscale.minor)),
             y.labs=NA, ytitle=bquote("log"[2] * "(Genes)"), y.title.line=1.5,
             pt.cex=0.85, parmar=c(2.75, 2.75, 0.3, 0.3))
axis(1, at=log10(c(500000, 5000000)), tick=F, 
     line=-0.7, labels=c("500kb", "5Mb"))
mtext(1, line=1.5, text=bquote("log"[10] * "(Segment Size)"))
axis(2, at=axTicks(2), tick=F, line=-0.65, labels=2^axTicks(2), las=2)
dev.off()

# Swarmplot of del vs dup effect size
pdf(paste(out.prefix, "gw_del_vs_dup.lnOR.pdf", sep="."),
    height=2, width=2)
segs.simple.vioswarm(gw, gw$max_ln_or, add.y.axis=T,
                     ytitle=bquote("Max." ~ italic("ln") * ("Odds Ratio")), ytitle.line=1,
           add.pvalue=T, parmar=c(1.2, 2.25, 1.5, 0.2))
dev.off()

# Swarmplot of del vs dup control frequency
pdf(paste(out.prefix, "gw_del_vs_dup.control_freq.pdf", sep="."),
    height=2, width=2)
segs.simple.vioswarm(gw, -log10(gw$pooled_control_freq), add.y.axis=F,
                     add.pvalue=T, parmar=c(1.2, 3.25, 1.5, 0.2))
axis(2, at=-log10(logscale.major), col=blueblack, tck=-0.03, labels=NA)
sapply(1:length(logscale.major), function(x){
  axis(2, at=-log10(logscale.major[x]), tick=F, las=2, line=-0.65,
       labels=bquote(10^.(log10(logscale.major[x]))))
})
mtext(2, text="CNV Freq. in Controls", line=2.25)
dev.off()

# Swarmplot of del vs dup avg. expression levels
pdf(paste(out.prefix, "gw_del_vs_dup.avg_expression.pdf", sep="."),
    height=2, width=2.2)
segs.simple.vioswarm(segs, segs$gene_expression_harmonic_mean, add.y.axis=T,
                     ytitle=bquote("Avg. Gene Expression"), ytitle.line=1.75,
                     pt.cex=0.75, add.pvalue=T, parmar=c(1.2, 2.8, 1.5, 0.2))
dev.off()

# Swarmplot of prop ubiquitously expressed for del vs dup
pdf(paste(out.prefix, "gw_del_vs_dup.prop_ubi_expressed.pdf", sep="."),
    height=2, width=2.2)
segs.simple.vioswarm(segs, segs$prop_ubiquitously_expressed, add.y.axis=T,
                     ytitle=bquote("Prop. Ubiquitously Expr."), ytitle.line=1.75,
                     pt.cex=0.75, add.pvalue=T, parmar=c(1.2, 2.8, 1.5, 0.2))
dev.off()

# Scatterplot of peak P-value and corresponding lnOR for all sites
pdf(paste(out.prefix, "all_segs.best_p_vs_or.pdf", sep="."),
    height=2.3, width=2.4)
segs.scatter(segs.all, x=log2(exp(segs.all$meta_best_lnor)), y=segs.all$meta_best_p, 
             subset_to_regions=segs.all$region_id[which(!is.infinite(segs.all$meta_best_p))],
             horiz.lines.at=c(gw.sig, -log10(0.05)), horiz.lines.lty=c(2, 1),
             xtitle=bquote(ln("Odds Ratio")), ytitle=bquote("Best -log"[10] * (italic(P))),
             x.title.line=1.6, y.title.line=1.5,
             add.lm=F, pt.cex=0.75, parmar=c(2.75, 2.75, 0.2, 0.2))
x.bump <- 0.04 * (par("usr")[2] - par("usr")[1])
y.bump <- 0.04 * (par("usr")[4] - par("usr")[3])
text(x=par("usr")[2] + x.bump, y=-log10(0.05) + y.bump, labels=format.pval(0.05), cex=0.75, pos=2, col=blueblack)
text(x=par("usr")[2] + x.bump, y=gw.sig + 1.5*y.bump, labels=format.pval(10^-gw.sig, nsmall=1), cex=0.75, pos=2, col=blueblack)
text(x=par("usr")[2] + x.bump, y=8.3, labels="Genome-wide\nsignificance", cex=0.75, pos=2, font=3, col=blueblack)
# text(x=par("usr")[2] + x.bump, y=2.4, labels="Nominal\nsignificance", cex=0.65, pos=2, font=3, col=blueblack)
# text(x=par("usr")[2] + x.bump, y=-log10(0.05) + y.bump, labels=format.pval(0.05), cex=0.75, pos=2)
dev.off()
