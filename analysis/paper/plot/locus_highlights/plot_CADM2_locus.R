#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot CADM2 intronic deletion rCNV locus highlight for rCNV2 paper


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
  make_option(c("--pips"), help="BED file of PIPs for all genes."),
  make_option(c("--rnaseq"), help="RNAseq BED file.")
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
pips.in <- opts$pips
rnaseq.in <- opts$rnaseq

# # DEV PARAMETERS
# cnvlist.in <- "~/scratch/cnvlist.tsv"
# sumstats.in <- "~/scratch/HP0000752.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz"
# sample.size.table <- "~/scratch/HPOs_by_metacohort.table.tsv"
# genome.in <- "~/scratch/GRCh37.genome"
# out.prefix <- "~/scratch/noncoding_locus_highlights_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# gtf.in <- "~/scratch/gencode.v19.annotation.CADM2.gtf.gz"
# rnaseq.in <- "~/scratch/CADM2.fetal_cortex_RNAseq.3kb_bins.bed.gz"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/locus_highlights/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Set locus parameters
chrom <- 3
highlight.start <- 85160000
highlight.end <- 85560000
highlight.width <- highlight.end - highlight.start
crb.starts <- c(85245894, 85300195)
crb.ends <- c(85278400, 85314000)
rss.starts <- c(85288006, 85560906)
rss.ends <- c(85288206, 85561148)
viewpoint.buffer.left <- 0.75 * highlight.width
viewpoint.buffer.right <- 1.5 * highlight.width
start <- floor(max(c(0, highlight.start - viewpoint.buffer.left)))
end <- ceiling(highlight.end + viewpoint.buffer.right)
region <- paste(chrom, ":", start, "-", end, sep="")
cnv.type <- "DEL"
all.case.hpos <- c("HP:0000707")
highlight.case.hpo <- "HP:0000752"
credset.genes <- c("CADM2")

# # Required for bedr::tabix in local Rstudio only:
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),
#                         "/Users/collins/anaconda3/envs/py3/bin",
#                         sep = ":"))

# Load genes
genes <- load.genes(gtf.in, region)
elig.tx <- unique(genes$transcript[which(genes$end - genes$start > 10000 & genes$gene=="CADM2" & genes$feature=="transcript")])
genes <- genes[which(genes$transcript %in% elig.tx), ]

# Load meta-analysis summary stats
ss <- load.sumstats(sumstats.in, region)
ss$meta_phred_p[which(is.na(as.numeric(ss$meta_phred_p)))] <- 0
ss$meta_lnOR[which(is.na(as.numeric(ss$meta_lnOR)))] <- 0
ss$meta_lnOR_lower[which(is.na(as.numeric(ss$meta_lnOR_lower)))] <- -100
ss$meta_lnOR_upper[which(is.na(as.numeric(ss$meta_lnOR_upper)))] <- 100

# Load sample sizes
n.samples <- get.sample.sizes(sample.size.table, all.case.hpos)
n.samples.highlight <- get.sample.sizes(sample.size.table, highlight.case.hpo)

# Load RNAseq BED
rnaseq <- load.feature.bed(rnaseq.in, keep.chrom=chrom, keep.col=6)
max.rnaseq.value <- quantile(rnaseq$value, 0.995)
rnaseq$value[which(rnaseq$value >= max.rnaseq.value)] <- max.rnaseq.value

# Load all CNVs
cnvs <- load.cnvs.multi(cnvlist.in, region, cnv.type, all.case.hpos)

# Notes: panel ordering
# idiogram (above y=0)
# coordinate line (at y=0)
# CRBs
# all transcripts
# RNAseq
# RS sites
# cnv pileups
# Set plot values
coord.panel.height <- 0.2
coord.panel.space <- 0.05
coord.panel.y0 <- 0
idio.panel.height <- 0.1
idio.panel.space <- 0.4
idio.panel.y0 <- (0.5 * coord.panel.height) + idio.panel.space + (0.5*idio.panel.height)
# Top panels = those at or above y=0 (idiogram & coordinate line)
top.panels.height <- sum(c(idio.panel.height, idio.panel.space, 0.5*coord.panel.height))
# Upper panels = all panels below y=0 but above CNV pileups
crb.panel.space <- 0.05
crb.panel.height <- 0.1
crb.panel.y0 <- -((0.5*coord.panel.height) + coord.panel.space + (0.5*crb.panel.height))
gene.panel.space <- 0.05
gene.panel.height <- 0.4
gene.panel.y0 <- crb.panel.y0 - (0.5*crb.panel.height) - gene.panel.space - (0.5*gene.panel.height)
rna.panel.space <- 0.075
rna.panel.height <- 0.35
rna.panel.y0 <- gene.panel.y0 - (0.5*gene.panel.height) - rna.panel.space - (0.5*rna.panel.height)
rss.panel.space <- 0.05
rss.panel.height <- 0.1
rss.panel.y0 <- rna.panel.y0 - (0.5*rna.panel.height) - rss.panel.space - (0.5*rss.panel.height)
upper.panels.height <- -sum(c((0.5 * coord.panel.height), coord.panel.space,
                              crb.panel.height, crb.panel.space,
                              gene.panel.height, gene.panel.space, 
                              rna.panel.height, rna.panel.space,
                              rss.panel.height, rss.panel.space))
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
pdf(paste(out.prefix, "locus_highlight.CADM2.pdf", sep="."), height=3.5*2, width=3.5*2)
par(mar=c(0.7, 7, 0, 0.5), bty="n")
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
blueshade.ybottoms <- c(rna.panel.y0+(0.5*rna.panel.height),
                        cnv.panel.y0s+(0.5*cnv.panel.height))
blueshade.ytops <- c(rna.panel.y0-(0.5*rna.panel.height),
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
add.genes(genes[which(genes$gene == "CADM2"), ], n.rows=3, transcripts=TRUE,
          y0=gene.panel.y0, panel.height=gene.panel.height, col=graphabs.green,
          y.axis.title="", y.axis.title.col="black", label.genes=c(), mark.tss=TRUE)
axis(2, at=gene.panel.y0+(0.05*gene.panel.height), tick=F, line=-0.9, las=2, labels=bquote(italic("CADM2")), padj=0)
axis(2, at=gene.panel.y0-(0.05*gene.panel.height), tick=F, line=-0.9, las=2, labels="Transcripts", padj=1)

# Add RNAseq
add.feature.barplot(rnaseq, y0=rna.panel.y0, col=graphabs.green, panel.height=rna.panel.height,
                    ytitle="Fetal Cortex\nRNA-seq")
Arrows(x0=(rss.starts[1]+rss.ends[1])/2, x1=(rss.starts[1]+rss.ends[1])/2,
       y0=rna.panel.y0+(0.4*rna.panel.height), y1=rna.panel.y0+(0.05*rna.panel.height),
       arr.type="triangle", col=blueblack, lwd=2, arr.adj=1, arr.length=0.175, arr.width=0.15)

# Add RSSs
add.rects(rss.starts, rss.ends, y0=rss.panel.y0, panel.height=rss.panel.height,
          col=graphabs.green, y.axis.title="")
text(x=min(rss.starts), y=rss.panel.y0-(0.1*rss.panel.height), pos=2, cex=5.5/6, 
     labels="Recursive Splice Sites", col=graphabs.green)


# Add CNV panels
if(cnv.type=="DEL"){
  cnv.title <- "Deletion Evidence per Cohort"
}else if(cnv.type=="DUP"){
  cnv.title <- "Duplication Evidence per Cohort"
}
text(x=mean(c(start, end)), y=upper.panels.height-(1.25*cnv.panel.spacing), labels=cnv.title, pos=3, cex=5.5/6)
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

