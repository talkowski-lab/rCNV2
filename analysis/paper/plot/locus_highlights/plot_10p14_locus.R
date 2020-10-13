#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot 10p14 intergenic duplication rCNV locus highlight for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)
require(shape, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
  make_option(c("--gtf"), help="GTF for plotting gene bodies. Must be tabix indexed."),
  make_option(c("--ncrna-gtf"), help="GTF for plotting gene bodies. Must be tabix indexed."),
  make_option(c("--pips"), help="BED file of PIPs for all genes."),
  make_option(c("--rnaseq"), help="RNAseq BED file."),
  make_option(c("--chipseq"), help="ChIP-seq BED file."),
  make_option(c("--esc-chromhmm-tracks"), help=".tsv of ESC ChromHMM BED files to plot."),
  make_option(c("--adult-chromhmm-tracks"), help=".tsv of adult cortex ChromHMM BED files to plot."),
  make_option(c("--chromhmm-manifest"), help=".tsv of ChromHMM state manifest.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog cnvlist sumstats samples_per_hpo.tsv", 
                                            "genome.tsv out_prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: cnvlist.tsv, sumstats.bed,", 
             "samples_per_hpo.tsv, genome.tsv, and output_prefix\n"))
}

# Writes args & opts to vars
cnvlist.in <- args$args[1]
sumstats.in <- args$args[2]
sample.size.table <- args$args[3]
genome.in <- args$args[4]
out.prefix <- args$args[5]
rcnv.config <- opts$`rcnv-config`
gtf.in <- opts$gtf
ncrna.gtf.in <- opts$`ncrna-gtf`
pips.in <- opts$pips
rnaseq.in <- opts$rnaseq
chipseq.in <- opts$chipseq
esc.chromhmm.tracks.in <- opts$`esc-chromhmm-tracks`
adult.chromhmm.tracks.in <- opts$`adult-chromhmm-tracks`
chromhmm.manifest.in <- opts$`chromhmm-manifest`

# # DEV PARAMETERS
# cnvlist.in <- "~/scratch/cnvlist.tsv"
# sumstats.in <- "~/scratch/HP0012759.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz"
# sample.size.table <- "~/scratch/HPOs_by_metacohort.table.tsv"
# genome.in <- "~/scratch/GRCh37.genome"
# out.prefix <- "~/scratch/noncoding_locus_highlights_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# gtf.in <- "~/scratch/gencode.v19.canonical.pext_filtered.gtf.gz"
# ncrna.gtf.in <- "~/scratch/gencode.v19.annotation.chr10_ncRNAs.gtf.gz"
# chipseq.in <- "~/scratch/10p14.ESC_H3k27ac.1kb_bins.bed.gz"
# esc.chromhmm.tracks.in <- "~/scratch/10p14.chromhmm_paths.tsv"
# adult.chromhmm.tracks.in <- "~/scratch/CADM2.chromhmm_paths.tsv"
# chromhmm.manifest.in <- "~/scratch/REP_state_manifest.tsv"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/locus_highlights/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Set locus parameters
chrom <- 10
highlight.start <- 6620000
highlight.end <- 6920000
highlight.width <- highlight.end - highlight.start
crb.starts <- c(6756200, 6816600)
crb.ends <- c(6793711, 6838000)
viewpoint.buffer.left <- (1/3) * highlight.width
viewpoint.buffer.right <- (6/5) * highlight.width
start <- floor(max(c(0, highlight.start - viewpoint.buffer.left)))
end <- ceiling(highlight.end + viewpoint.buffer.right)
region <- paste(chrom, ":", start, "-", end, sep="")
cnv.type <- "DUP"
all.case.hpos <- c("HP:0000707")
highlight.case.hpo <- "HP:0012759"

# # Required for bedr::tabix in local Rstudio only:
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),
#                         "/Users/collins/anaconda3/envs/py3/bin",
#                         sep = ":"))

# Load genes
genes <- load.genes(gtf.in, region)
ncrnas <- load.genes(ncrna.gtf.in, region)

# Load meta-analysis summary stats
ss <- load.sumstats(sumstats.in, region)

# Load sample sizes
n.samples <- get.sample.sizes(sample.size.table, all.case.hpos)
n.samples.highlight <- get.sample.sizes(sample.size.table, highlight.case.hpo)

# Load H3K27ac ChIP-seq
h3k27ac <- load.feature.bed(chipseq.in, region, keep.col=6)

# Load ChromHMM BEDs & annotate with state colors
esc.chmm.tracks <- load.chromhmm.tracks(esc.chromhmm.tracks.in, chromhmm.manifest.in, region)
adult.chmm.tracks <- load.chromhmm.tracks(adult.chromhmm.tracks.in, chromhmm.manifest.in, region)

# Load all CNVs
cnvs <- load.cnvs.multi(cnvlist.in, region, cnv.type, all.case.hpos)

# Notes: panel ordering
# idiogram (above y=0)
# coordinate line (at y=0)
# CRBs
# protein-coding genes
# ncRNAs
# ESC H3K27ac chip-seq
# ESC ChromHMM
# Adult cortex ChromHMM
# cnv pileups
# Set plot values
coord.panel.height <- 0.2
coord.panel.space <- 0.05
coord.panel.y0 <- 0
idio.panel.height <- 0.1
idio.panel.space <- 0.3
idio.panel.y0 <- (0.5 * coord.panel.height) + idio.panel.space + (0.5*idio.panel.height)
# Top panels = those at or above y=0 (idiogram & coordinate line)
top.panels.height <- sum(c(idio.panel.height, idio.panel.space, 0.5*coord.panel.height))
# Upper panels = all panels below y=0 but above CNV pileups
crb.panel.space <- 0.05
crb.panel.height <- 0.1
crb.panel.y0 <- -((0.5*coord.panel.height) + coord.panel.space + (0.5*crb.panel.height))
pcgene.panel.space <- 0.15
pcgene.panel.height <- 0.1
pcgene.panel.y0 <- crb.panel.y0 - (0.5*crb.panel.height) - pcgene.panel.space - (0.5*pcgene.panel.height)
nc.panel.space <- 0.075
nc.panel.height <- 0.1
nc.panel.y0 <- pcgene.panel.y0 - (0.5*pcgene.panel.height) - nc.panel.space - (0.5*nc.panel.height)
h3k27ac.panel.space <- 0.1
h3k27ac.panel.height <- 0.3
h3k27ac.panel.y0 <- nc.panel.y0 - (0.5*nc.panel.height) - h3k27ac.panel.space - (0.5*h3k27ac.panel.height)
esc.chmm.panel.space <- 0.075
esc.chmm.panel.height <- 0.3
esc.chmm.panel.y0 <- h3k27ac.panel.y0 - (0.5*h3k27ac.panel.height) - esc.chmm.panel.space - (0.5*esc.chmm.panel.height)
adult.chmm.panel.space <- 0.1
adult.chmm.panel.height <- 0.3
adult.chmm.panel.y0 <- esc.chmm.panel.y0 - (0.5*esc.chmm.panel.height) - adult.chmm.panel.space - (0.5*adult.chmm.panel.height)
upper.panels.height <- -sum(c((0.5 * coord.panel.height), coord.panel.space,
                              crb.panel.height, crb.panel.space,
                              pcgene.panel.height, pcgene.panel.space, 
                              nc.panel.height, nc.panel.space, 
                              h3k27ac.panel.height, h3k27ac.panel.space,
                              esc.chmm.panel.height, esc.chmm.panel.space,
                              adult.chmm.panel.height, adult.chmm.panel.space))
# Lower panels = CNV pileup
cnv.panel.height <- 0.65
cnv.panel.spacing <- 0.15
cnv.panel.y0s <- upper.panels.height - sapply(1:4, function(i){((i-1)*(cnv.panel.height+cnv.panel.spacing)) + (0.5*cnv.panel.height)}) - (1.1*cnv.panel.spacing)
total.height <- min(cnv.panel.y0s - (0.5*cnv.panel.height))
cnv.key.height <- 0.25
cnv.key.spacing <- 0.15
cnv.key.y0 <- total.height - cnv.key.spacing - (0.5*cnv.key.height)
total.height.plus.key <- cnv.key.y0 + (0.5*cnv.key.height)

# Prep locus plot
pdf(paste(out.prefix, "locus_highlight.10p14.pdf", sep="."), height=3.5*2, width=3.25*2)
par(mar=c(0.7, 6.5, 0, 0.5), bty="n")
plot(NA, xlim=c(start, end), ylim=c(total.height.plus.key, top.panels.height), 
     yaxt="n", xaxt="n", ylab="", xlab="")

# Add stick idiogram
add.idio.stick(genome.in, chrom, highlight.start, highlight.end, idio.panel.y0,
               tick.height=-0.03)

# Add background shading
rect(xleft=highlight.start, xright=crb.starts[1], 
     ybottom=total.height, ytop=0, 
     col=adjustcolor(highlight.color, alpha=0.3), border=NA)
rect(xleft=crb.starts[1], xright=crb.ends[2], 
     ybottom=total.height, ytop=0, 
     col=adjustcolor(highlight.color, alpha=0.6), border=NA)
rect(xleft=crb.ends[2], xright=highlight.end, 
     ybottom=total.height, ytop=0, 
     col=adjustcolor(highlight.color, alpha=0.3), border=NA)
# rect(xleft=rss.starts, xright=rss.ends, 
#      ybottom=total.height, ytop=0, 
#      col=control.cnv.colors[3], border=control.cnv.colors[3])
blueshade.ybottoms <- c(h3k27ac.panel.y0+(0.5*h3k27ac.panel.height),
                        cnv.panel.y0s+(0.5*cnv.panel.height))
blueshade.ytops <- c(h3k27ac.panel.y0-(0.5*h3k27ac.panel.height),
                     cnv.panel.y0s-(0.5*cnv.panel.height))
rect(xleft=par("usr")[1], xright=par("usr")[2], 
     ybottom=blueshade.ybottoms,
     ytop=blueshade.ytops,
     border=NA, bty="n", col=bluewhite)

# Add coordinate line
add.coord.line(start, end, coord.panel.y0, highlight.start, highlight.end, 
               highlight.col=highlight.color, tick.height=0.04, vlines=TRUE, 
               vlines.bottom=total.height)

# Add CRBs
add.rects(crb.starts, crb.ends, crb.panel.y0, panel.height=crb.panel.height, 
          col=blueblack, border=blueblack, y.axis.title="")
text(x=min(crb.starts), y=crb.panel.y0-(0.1*crb.panel.height), pos=2, cex=5.5/6, 
     labels="Significant CRBs", col=blueblack)

# Add all genes
add.genes(genes, n.rows=1, transcripts=FALSE,
          y0=pcgene.panel.y0, panel.height=pcgene.panel.height, col=graphabs.green,
          y.axis.title="Coding Genes", label.genes=c("PRKCQ", "SFMBT2"))

# Add ncRNAs
add.genes(ncrnas, n.rows=1, transcripts=FALSE,
          y0=nc.panel.y0, panel.height=nc.panel.height, col=cnv.colors[3],
          y.axis.title="Noncoding RNAs")

# Add cortex H3k27ac ChIP-seq
add.feature.barplot(h3k27ac, y0=h3k27ac.panel.y0, col=h3k27ac.color,
                    panel.height=h3k27ac.panel.height, ytitle="ESC H3K27ac")

# Add ESC ChromHMM tracks
add.chromhmm(esc.chmm.tracks, y0=esc.chmm.panel.y0, panel.height=esc.chmm.panel.height,
             y.axis.title="ESC\nChromatin State")

# Add adult cortex ChromHMM tracks
add.chromhmm(adult.chmm.tracks, y0=adult.chmm.panel.y0, panel.height=adult.chmm.panel.height,
             y.axis.title="Adult Cortex\nChromatin State")

# Add CNV panels
if(cnv.type=="DEL"){
  cnv.title <- "Deletion Evidence per Cohort"
}else if(cnv.type=="DUP"){
  cnv.title <- "Duplication Evidence per Cohort"
}
text(x=mean(c(start, end)), y=upper.panels.height-(1.4*cnv.panel.spacing), labels=cnv.title, pos=3, cex=5.5/6)
legend.sides <- rep("right", 4)
topbottoms <- rep("top", 4)
sapply(1:4, function(i){
  add.cnv.panel(cnvs[[i]], n.case=n.samples$case[i], n.ctrl=n.samples$ctrl[i], 
                highlight.hpo=highlight.case.hpo, y0=cnv.panel.y0s[i], 
                cnv.type=cnv.type, panel.height=cnv.panel.height, 
                y.axis.title=paste(cnv.type, "Freq.", sep="\n"), expand.pheno.label=F,
                case.legend.side=legend.sides[i], case.legend.topbottom=topbottoms[i],
                ctrl.legend.side=legend.sides[i], add.cohort.label=TRUE, cohort.label=paste("Cohort", i))
})

# Add key for CNV panels
add.cnv.key(cnv.type, cnv.key.y0, total.n.ctrls=sum(n.samples$ctrl), 
            all.case.hpos, total.n.cases=sum(n.samples$case),
            highlight.case.hpo, total.n.cases.highlight=sum(n.samples.highlight$case),
            panel.height=cnv.key.height)

# Close pdf
dev.off()

