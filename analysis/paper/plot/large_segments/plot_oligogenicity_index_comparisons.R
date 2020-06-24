#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Comparison of monogenicity vs. oligogenicity for large segments for rCNV paper


options(stringsAsFactors=F, scipen=1000)
csqs <- c("lof", "mis", "syn")


######################
### DATA FUNCTIONS ###
######################
# Load table of mutation rates
load.mutrates <- function(mutrates.in){
  mutrates <- read.table(mutrates.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(mutrates)[1] <- "gene"
  return(mutrates)
}

# Load count of de novo mutations for a single study and residualize vs. expected
load.dnms <- function(dnms.in, mutrates){
  dnms <- read.table(dnms.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(dnms)[1] <- "gene"
  dnms <- merge(dnms, mutrates, all.x=T, all.y=F, sort=F, by="gene")
  for(csq in csqs){
    counts <- dnms[, which(colnames(dnms) == csq)]
    mus <- dnms[, which(colnames(dnms) == paste("mu", csq, sep="_"))]
    n.dnms <- sum(counts, na.rm=T)
    scaled.mus <- (mus / sum(mus, na.rm=T)) * n.dnms
    residuals <- counts - scaled.mus
    dnms[paste("excess", csq, sep="_")] <- residuals
  }
  return(dnms)
}

# Calculate oligogenicity index (O_i)
# Defined as the minimum number of genes that captures X% of the total DNM excess in a segment
calc.oligo.index <- function(segs, dnms, csq, cutoff=0.90){
  sapply(segs$genes, function(genes){
    x <- dnms[which(dnms$gene %in% genes), which(colnames(dnms) == paste("excess", csq, sep="_"))]
    x <- x[which(!is.na(x) & x>0)]
    if(length(x) > 0){
      min(which(cumsum(x[order(-x)])/sum(x) >= cutoff))
    }else{
      NA
    }
  })
}


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
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv ddd.tsv asc.tsv mutrates.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Six positional arguments required: loci.bed, segs.tsv, ddd.tsv, asc.tsv, mutrates.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
ddd.dnms.in <- args$args[3]
asc.dnms.in <- args$args[4]
gene.mutrates.in <- args$args[5]
out.prefix <- args$args[6]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# ddd.dnms.in <- "~/scratch/ddd_dnm_counts.tsv.gz"
# asc.dnms.in <- "~/scratch/asc_dnm_counts.tsv.gz"
# gene.mutrates.in <- "~/scratch/gene_mutation_rates.tsv.gz"
# out.prefix <- "~/scratch/test_effect_sizes"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load dnm counts and residualize by mutation rates
mutrates <- load.mutrates(gene.mutrates.in)
ddd <- load.dnms(ddd.dnms.in, mutrates)
asc <- load.dnms(asc.dnms.in, mutrates)

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Restrict to segments nominally significant in at least one phenotype
segs <- segs[which(segs$nom_sig), ]

# Annotate segs with oligogenicity indexes
for(csq in csqs){
  segs[paste("ddd", csq, "oligo", sep=".")] <- calc.oligo.index(segs, ddd, csq=csq)
  segs[paste("asc", csq, "oligo", sep=".")] <- calc.oligo.index(segs, asc, csq=csq)
}

# Get list of neuro loci
neuro.region_ids <- get.neuro.region_ids(loci, segs)
neuro.segs <- segs[which(segs$region_id %in% neuro.region_ids), ]

# Swarmplots of DDD oligogenicity index for all deletions & duplications
pdf(paste(out.prefix, "neuro_segs.ddd_lof_oligogenicity_index.pdf", sep="."),
    height=2.25, width=2)
segs.simple.vioswarm(neuro.segs, y=neuro.segs$ddd.lof.oligo, add.pvalue=T, ytitle="",
                     pt.cex=0.65, parmar=c(1.2, 3.5, 1.5, 0))
mtext(2, line=2.35, text=bquote("Genes Capturing" >= "90%"))
mtext(2, line=1.45, text=bquote("of Excess" ~ italic("dn") * "PTVs"))
dev.off()
pdf(paste(out.prefix, "neuro_segs.ddd_mis_oligogenicity_index.pdf", sep="."),
    height=2.25, width=1.6)
segs.simple.vioswarm(neuro.segs, y=neuro.segs$ddd.mis.oligo, add.pvalue=T, 
                     ytitle="Oligogenicity Index",
                     pt.cex=0.65, parmar=c(1.2, 2.65, 1.5, 0))
dev.off()

# Scatterplots of DDD oligogenicity index vs genes in segment for all deletions & duplications
pdf(paste(out.prefix, "neuro_segs.ddd_lof_oligogenicity_index_vs_ngenes.pdf", sep="."),
    height=2.4, width=2.4)
segs.scatter(neuro.segs, x=neuro.segs$n_genes, y=neuro.segs$ddd.lof.oligo,
             xlims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             ylims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             xtitle="Genes in Segment", x.title.line=1.4, ytitle="",
             add.lm=T, abline.a=0, abline.b=0.9, abline.lty=2,
             pt.cex=0.85, parmar=c(3.5, 3.5, 1, 1))
text(x=mean(par("usr")[1:2]), y=mean(par("usr")[3:4])+(0.04*(par("usr")[4]-par("usr")[3])),
     col=control.cnv.colors[2], labels=bquote("Uniform" ~ italic("dn") * "PTV" ~ "distribution"), 
     srt=45*0.925, font=3, cex=0.9)
mtext(2, line=2.35, text=bquote("Genes Capturing" >= "90%"))
mtext(2, line=1.55, text=bquote("of Excess" ~ italic("dn") * "PTVs"))
dev.off()
pdf(paste(out.prefix, "neuro_segs.ddd_mis_oligogenicity_index_vs_ngenes.pdf", sep="."),
    height=2, width=2)
segs.scatter(neuro.segs, x=neuro.segs$n_genes, y=neuro.segs$ddd.mis.oligo,
             xlims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             ylims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             xtitle="Genes in Segment", x.title.line=1.4,
             ytitle="Oligogenicity Index", y.title.line=1.5,
             add.lm=T, abline.a=0, abline.b=0.9, abline.lty=2,
             pt.cex=0.85, parmar=c(2.6, 2.6, 0.2, 0.2))
text(x=mean(par("usr")[1:2]), y=mean(par("usr")[3:4])+(0.04*(par("usr")[4]-par("usr")[3])),
     col=control.cnv.colors[2], labels=bquote("Uniform" ~ italic("dn") * "Mis." ~ "distribution"), 
     srt=45*0.925, font=3, cex=0.9)
dev.off()

# Swarmplot of ASC oligogenicity index for all deletions & duplications
pdf(paste(out.prefix, "neuro_segs.asc_lof_oligogenicity_index.pdf", sep="."),
    height=2.25, width=1.6)
segs.simple.vioswarm(neuro.segs, y=neuro.segs$asc.lof.oligo, add.pvalue=T, 
                     ytitle="Oligogenicity Index",
                     pt.cex=0.65, parmar=c(1.2, 2.65, 1.5, 0))
dev.off()
pdf(paste(out.prefix, "neuro_segs.asc_mis_oligogenicity_index.pdf", sep="."),
    height=2.25, width=1.6)
segs.simple.vioswarm(neuro.segs, y=neuro.segs$asc.mis.oligo, add.pvalue=T, 
                     ytitle="Oligogenicity Index",
                     pt.cex=0.65, parmar=c(1.2, 2.65, 1.5, 0))
dev.off()

# Scatterplots of ASC oligogenicity index vs genes in segment for all deletions & duplications
pdf(paste(out.prefix, "neuro_segs.asc_lof_oligogenicity_index_vs_ngenes.pdf", sep="."),
    height=2, width=2)
segs.scatter(neuro.segs, x=neuro.segs$n_genes, y=neuro.segs$asc.lof.oligo,
             xlims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             ylims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             xtitle="Genes in Segment", x.title.line=1.4,
             ytitle="Oligogenicity Index", y.title.line=1.5,
             add.lm=T, abline.a=0, abline.b=0.9, abline.lty=2,
             pt.cex=0.85, parmar=c(2.6, 2.6, 0.2, 0.2))
text(x=mean(par("usr")[1:2]), y=mean(par("usr")[3:4])+(0.04*(par("usr")[4]-par("usr")[3])),
     col=control.cnv.colors[2], labels=bquote("Uniform" ~ italic("dn") * "PTV" ~ "distribution"), 
     srt=45*0.925, font=3, cex=0.9)
dev.off()
pdf(paste(out.prefix, "neuro_segs.asc_mis_oligogenicity_index_vs_ngenes.pdf", sep="."),
    height=2, width=2)
segs.scatter(neuro.segs, x=neuro.segs$n_genes, y=neuro.segs$asc.mis.oligo,
             xlims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             ylims=c(0, max(neuro.segs$n_genes, na.rm=T)), 
             xtitle="Genes in Segment", x.title.line=1.4,
             ytitle="Oligogenicity Index", y.title.line=1.5,
             add.lm=T, abline.a=0, abline.b=0.9, abline.lty=2,
             pt.cex=0.85, parmar=c(2.6, 2.6, 0.2, 0.2))
text(x=mean(par("usr")[1:2]), y=mean(par("usr")[3:4])+(0.04*(par("usr")[4]-par("usr")[3])),
     col=control.cnv.colors[2], labels=bquote("Uniform" ~ italic("dn") * "Mis." ~ "distribution"), 
     srt=45*0.925, font=3, cex=0.9)
dev.off()


# Swarmplot of DDD oligogenicity indexes vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.ddd_lof_oligogenicity_index.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(neuro.segs, x.bool=neuro.segs$nahr, y=neuro.segs$ddd.lof.oligo, 
           x.labs=c("Nonrecurrent", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()
pdf(paste(out.prefix, "segs_by_mechanism.ddd_mis_oligogenicity_index.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(neuro.segs, x.bool=neuro.segs$nahr, y=neuro.segs$ddd.mis.oligo, 
           x.labs=c("Nonrecurrent", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of ASC oligogenicity indexes vs. mechanism
pdf(paste(out.prefix, "segs_by_mechanism.asc_lof_oligogenicity_index.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(neuro.segs, x.bool=neuro.segs$nahr, y=neuro.segs$asc.lof.oligo, 
           x.labs=c("Nonrecurrent", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()
pdf(paste(out.prefix, "segs_by_mechanism.asc_mis_oligogenicity_index.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(neuro.segs, x.bool=neuro.segs$nahr, y=neuro.segs$asc.mis.oligo, 
           x.labs=c("Nonrecurrent", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.75, 
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

