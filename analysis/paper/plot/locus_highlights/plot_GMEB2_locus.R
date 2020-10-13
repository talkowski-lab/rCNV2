#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot GMEB2 rCNV locus highlight for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced."),
  make_option(c("--gtf"), help="GTF for plotting gene bodies. Must be tabix indexed."),
  make_option(c("--pips"), help="BED file of PIPs for all genes."),
  make_option(c("--gw-sig"), type="numeric", default=3.715428E-6,
              help="P-value cutoff to mark as genome-wide significant.")
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
gw.sig <- opts$`gw-sig`

# # DEV PARAMETERS
# cnvlist.in <- "~/scratch/cnvlist.tsv"
# sumstats.in <- "~/scratch/HP0012639.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz"
# sample.size.table <- "~/scratch/HPOs_by_metacohort.table.tsv"
# genome.in <- "~/scratch/GRCh37.genome"
# out.prefix <- "~/scratch/noncoding_locus_highlights_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# gtf.in <- "~/scratch/gencode.v19.canonical.pext_filtered.gtf.gz"
# pips.in <- "~/scratch/rCNV.DUP.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/locus_highlights/"
# gw.sig <- 3.715428E-6

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Set locus parameters
chrom <- 20
highlight.start <- 62152076
highlight.end <- 62251229
highlight.width <- highlight.end - highlight.start
viewpoint.buffer <- 2 * highlight.width
start <- floor(max(c(0, highlight.start - viewpoint.buffer)))
end <- ceiling(highlight.end + viewpoint.buffer)
region <- paste(chrom, ":", start, "-", end, sep="")
cnv.type <- "DUP"
all.case.hpos <- c("HP:0000707")
highlight.case.hpo <- "HP:0012639"
credset.genes <- c("GMEB2", "HELZ2", "PPDPF", "PTK6")

# # Required for bedr::tabix in local Rstudio only:
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),
#                         "/Users/collins/anaconda3/envs/py3/bin",
#                         sep = ":"))

# Load genes
genes <- load.genes(gtf.in, region)

# Load meta-analysis summary stats
ss <- load.sumstats(sumstats.in, region)

# Load PIPs for genes in region for highlight HPO
pips <- load.pips(pips.in, unique(genes$gene), highlight.case.hpo)

# Load sample sizes
n.samples <- get.sample.sizes(sample.size.table, all.case.hpos)
n.samples.highlight <- get.sample.sizes(sample.size.table, highlight.case.hpo)

# Load all CNVs
cnvs <- load.cnvs.multi(cnvlist.in, region, cnv.type, all.case.hpos)


# Notes: panel ordering
# idiogram (above y=0)
# coordinate line (at y=0)
# p-values
# odds ratios
# all genes
# PIPs
# credible set genes
# other genes
# cnv pileups
# Set plot values
coord.panel.height <- 0.2
coord.panel.space <- 0.2
coord.panel.y0 <- 0
idio.panel.height <- 0.1
idio.panel.space <- 0.2
idio.panel.y0 <- (0.5 * coord.panel.height) + idio.panel.space + (0.5*idio.panel.height)
# Top panels = those at or above y=0 (idiogram & coordinate line)
top.panels.height <- sum(c(idio.panel.height, idio.panel.space, 0.5*coord.panel.height))
# Upper panels = all panels below y=0 but above CNV pileups
pval.panel.space <- 0.3
pval.panel.height <- 0.4
pval.panel.y0 <- -((0.5*coord.panel.height) + coord.panel.space + (0.5*coord.panel.height))
or.panel.space <- 0.15
or.panel.height <- 0.4
or.panel.y0 <- pval.panel.y0 - (0.5*pval.panel.height) - or.panel.space - (0.5*or.panel.height)
gene.panel.space <- 0.1
allgene.panel.height <- 0.1
allgene.panel.y0 <- or.panel.y0 - (0.5*or.panel.height) - gene.panel.space - (0.5*allgene.panel.height)
pip.panel.space <- 0.1
pip.panel.height <- 0.4
pip.panel.y0 <- allgene.panel.y0 - (0.5*allgene.panel.height) - pip.panel.space - (0.5*pip.panel.height)
credset.panel.space <- gene.panel.space + 0.1
credset.panel.height <- 0.1
credset.panel.y0 <- pip.panel.y0 - (0.5*pip.panel.height) - credset.panel.space - (0.5*credset.panel.height)
othergene.panel.height <- 0.1
othergene.panel.y0 <- credset.panel.y0 - (0.5*credset.panel.height) - gene.panel.space - (0.5*othergene.panel.height)
upper.panels.height <- -sum(c((0.5 * coord.panel.height), coord.panel.space,
                              pval.panel.space, pval.panel.height, or.panel.space, or.panel.height,
                              allgene.panel.height, gene.panel.space, pip.panel.space, pip.panel.height,
                              credset.panel.height, credset.panel.space, othergene.panel.height, gene.panel.space))
# Lower panels = CNV pileup
cnv.panel.height <- 0.65
cnv.panel.spacing <- 0.15
cnv.panel.y0s <- upper.panels.height - sapply(1:4, function(i){((i-1)*(cnv.panel.height+cnv.panel.spacing)) + (0.5*cnv.panel.height)}) + (1.1*cnv.panel.spacing)
total.height <- min(cnv.panel.y0s - (0.5*cnv.panel.height))
cnv.key.height <- 0.25
cnv.key.spacing <- 0.15
cnv.key.y0 <- total.height - cnv.key.spacing - (0.5*cnv.key.height)
total.height.plus.key <- cnv.key.y0 + (0.5*cnv.key.height)

# Prep locus plot
pdf(paste(out.prefix, "locus_highlight.GMEB2.pdf", sep="."), height=3.75*2, width=3.25*2)
par(mar=c(0.2, 6, 0, 0.5), bty="n")
plot(NA, xlim=c(start, end), ylim=c(total.height.plus.key, top.panels.height), 
     yaxt="n", xaxt="n", ylab="", xlab="")

# Add stick idiogram
add.idio.stick(genome.in, chrom, highlight.start, highlight.end, idio.panel.y0,
               tick.height=-0.03)

# Add background shading
causal.start <- min(genes$start[which(genes$gene=="GMEB2")], na.rm=T)
causal.end <- max(genes$end[which(genes$gene=="GMEB2")], na.rm=T)
rect(xleft=highlight.start, xright=causal.start, 
     ybottom=total.height, ytop=0, 
     col=adjustcolor(highlight.color, alpha=0.3), border=NA)
rect(xleft=causal.start, xright=causal.end, 
     ybottom=total.height, ytop=0, 
     col=adjustcolor(highlight.color, alpha=0.6), border=NA)
rect(xleft=causal.end, xright=highlight.end, 
     ybottom=total.height, ytop=0, 
     col=adjustcolor(highlight.color, alpha=0.3), border=NA)
blueshade.ybottoms <- c(cnv.panel.y0s+(0.5*cnv.panel.height),
                        pval.panel.y0+(0.5*pval.panel.height),
                        or.panel.y0+(0.5*or.panel.height),
                        pip.panel.y0+(0.5*pip.panel.height))
blueshade.ytops <- c(cnv.panel.y0s-(0.5*cnv.panel.height),
                     pval.panel.y0-(0.5*pval.panel.height),
                     or.panel.y0-(0.5*or.panel.height),
                     pip.panel.y0-(0.5*pip.panel.height))
rect(xleft=par("usr")[1], xright=par("usr")[2], 
     ybottom=blueshade.ybottoms,
     ytop=blueshade.ytops,
     border=NA, bty="n", col=bluewhite)

# Add coordinate line
add.coord.line(start, end, coord.panel.y0, highlight.start, highlight.end, 
               highlight.col=highlight.color, tick.height=0.04, vlines=TRUE, vlines.bottom=total.height)

# Add association summary stats
add.pvalues(ss, y0=pval.panel.y0, panel.height=pval.panel.height, cnv.type=cnv.type, gw.sig=gw.sig)
add.ors(ss, y0=or.panel.y0, panel.height=or.panel.height, cnv.type=cnv.type)

# Add all genes
add.genes(genes, n.rows=1, 
          y0=allgene.panel.y0, panel.height=allgene.panel.height, col=blueblack,
          y.axis.title="All Genes", y.axis.title.col="black")

# Add PIPs
add.pips(pips, genes, y0=pip.panel.y0, panel.height=pip.panel.height, label.genes=c("GMEB2"),
         col=ns.color, highlight.col=graphabs.green, highlight.genes=credset.genes)

# Add finemapped genes
add.bracket(xleft=min(as.numeric(genes[which(genes$gene %in% credset.genes), "start"])),
            xright=max(as.numeric(genes[which(genes$gene %in% credset.genes), "end"])),
            y0=credset.panel.y0, height=2*credset.panel.height, col=graphabs.green,
            lwd=2, staple.wex=0.015)
add.genes(genes[which(genes$gene %in% credset.genes), ], n.rows=1, 
          y0=credset.panel.y0, panel.height=credset.panel.height, col=graphabs.green,
          y.axis.title="Credible Set", label.genes=c("GMEB2"))
add.genes(genes[which(!(genes$gene %in% credset.genes)), ], n.rows=1, 
          y0=othergene.panel.y0, panel.height=othergene.panel.height, col=ns.color,
          y.axis.title="Other Genes")

# Add CNV panels
if(cnv.type=="DEL"){
  cnv.title <- "Deletion Evidence per Cohort"
}else if(cnv.type=="DUP"){
  cnv.title <- "Duplication Evidence per Cohort"
}
text(x=mean(c(start, end)), y=upper.panels.height+(0.85*cnv.panel.spacing), labels=cnv.title, pos=3, cex=5.5/6)
legend.sides <- c("right", "left", "right", "right")
topbottoms <- c("bottom", "top", "top", "bottom")
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

