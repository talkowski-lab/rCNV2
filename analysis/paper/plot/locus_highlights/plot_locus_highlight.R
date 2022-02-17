#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Flexible locus highlight plotting script


#############
### SETUP ###
#############
# Load librarires
require(optparse, quietly=T)
require(rCNV2, quietly=T)

# Set defaults not exposed to user
options(stringsAsFactors=F, scipen=1000)
coord.panel.height <- 0.2
coord.panel.space <- 0.2
coord.panel.y0 <- 0
idio.panel.height <- 0.1
idio.panel.space <- 0.2
idio.panel.y0 <- (0.5 * coord.panel.height) + idio.panel.space + (0.5*idio.panel.height)
# Top panels = those at or above y=0 (idiogram & coordinate line)
top.panels.height <- sum(c(idio.panel.height, idio.panel.space, 0.5*coord.panel.height))


#####################
### RSCRIPT BLOCK ###
#####################
# List of command-line options
option_list <- list(
  make_option(c("--case-hpos"), default="HP:0000118",
              help=paste("semicolon-delimited list of case HPOs to plot.",
                         "[default: use all cases]")),
  make_option(c("--highlight-hpo"), default=NA, type="character",
              help="HPO to highlight in CNV panels. [default: no highlighted HPO]"),
  make_option(c("--highlights"), help=paste("semicolon-delimited list of regions",
                                            "to highlight. [default: no highlights]")),
  make_option(c("--highlight-color"), default=highlight.color,
              help="color to use for region highlights. [default: %default]"),
  make_option(c("--sumstats"), help="BED of summary statistics."),
  make_option(c("--cytobands", help=paste("BED of cytobands. If provided, will",
                                          "annotate idiogram with cytoband.",
                                          "[default: no annotation]"))),
  make_option(c("--gtf"), help="GTF for plotting gene bodies. Must be tabix indexed."),
  make_option(c("--label-genes"), help=paste("semicolon-delimited list of gene",
                                             "symbols to label. [default: no labels]")),
  make_option(c("--autolabel-n-genes"), default=4, type="numeric",
              help=paste("automatically label all genes for loci involving no more",
                         "than this many genes [default: %default]")),
  make_option(c("--pips"), help="BED file of PIPs for all genes."),
  make_option(c("--min-cases"), default=300, type="numeric",
              help=paste("minimum number of cases required for a cohort to be",
                         "plotted. [default: %default]")),
  make_option(c("--gw-sig"), type="numeric", default=10^-6,
              help="P-value cutoff to mark as genome-wide significant. [default: %default]"),
  make_option(c("--gw-sig-label"), type="character", default="Genome-wide significance",
              help="Text label to annotate above --gw-sig [default: %default]"),
  make_option(c("--gw-sig-label-side"), default="above",
              help=paste("orientation of genome-wide significance label relative",
                         "to horizontal line. [default: %default]")),
  make_option(c("--standardize-frequencies"), action="store_true", default=FALSE,
              help="standardize Y axis limits for all CNV panels. [default: %default]"),
  make_option(c("--collapse-cohorts"), action="store_true", default=FALSE,
              help="collapse all cohorts into a single CNV panel. [default: %default]"),
  make_option(c("--pval-panel-space"), default=0.3, type="numeric",
              help=paste("spacing between P-value panel and next panel. Only used",
                         "if --sumstats are provided. [default: %default]")),
  make_option(c("--pval-panel-height"), default=0.4, type="numeric",
              help=paste("hight of P-value panel. Only used if --sumstats are",
                         "provided. [default: %default]")),
  make_option(c("--or-panel-space"), default=0.15, type="numeric",
              help=paste("spacing between odds ratio panel and next panel. Only",
                         "used if --sumstats are provided. [default: %default]")),
  make_option(c("--or-panel-height"), default=0.4, type="numeric",
              help=paste("hight of odds ratio panel. Only used if --sumstats are",
                         "provided. [default: %default]")),
  make_option(c("--gene-panel-space"), default=0.15, type="numeric",
              help=paste("spacing between gene strip panels and next panels. Only",
                         "used if --gtf is provided. [default: %default]")),
  make_option(c("--gene-panel-height"), default=0.1, type="numeric",
              help=paste("hight of gene strip panels. Only used if --gtf is",
                         "provided. [default: %default]")),
  make_option(c("--pip-panel-space"), default=0.1, type="numeric",
              help=paste("spacing between PIP panel and next panel. Only used",
                         "if --pips are provided. [default: %default]")),
  make_option(c("--pip-panel-height"), default=0.4, type="numeric",
              help=paste("hight of PIP panel. Only used if --pips are",
                         "provided. [default: %default]")),
  make_option(c("--cnv-panel-space"), default=0.15, type="numeric",
              help="spacing between each CNV panel and next panel. [default: %default]"),
  make_option(c("--cnv-panel-height"), default=0.65, type="numeric",
              help="hight of each CNV panel. [default: %default]"),
  make_option(c("--pdf-height"), help=paste("height of output .pdf (at 6pt font",
                                            "resolution). [default %default]"),
              default=4.75, type="numeric"),
  make_option(c("--pdf-width"), help=paste("width of output .pdf (at 6pt font",
                                           "resolution). [default %default]"),
              default=3.5, type="numeric")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog region cnvlist cnvtype",
                                            "samples_per_hpo.tsv genome.tsv",
                                            "out_prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Six positional arguments required: region cnvlist.tsv, cnv.type, ",
             "samples_per_hpo.tsv, genome.tsv, and output_prefix\n"))
}

# Writes args & opts to vars
region <- args$args[1]
cnvlist.in <- args$args[2]
cnv.type <- args$args[3]
sample.size.table <- args$args[4]
genome.in <- args$args[5]
out.prefix <- args$args[6]
case.hpo.str <- opts$`case-hpos`
highlight.hpo <- opts$`highlight-hpo`
highlights.str <- opts$`highlights`
highlight.color <- opts$`highlight-color`
sumstats.in <- opts$sumstats
cytobands.in <- opts$cytobands
gtf.in <- opts$gtf
label.genes.str <- opts$`label-genes`
autolabel.n.genes <- opts$`autolabel-n-genes`
pips.in <- opts$pips
min.cases <- opts$`min-cases`
gw.sig <- opts$`gw-sig`
gw.sig.label <- opts$`gw-sig-label`
gw.sig.label.side <- opts$`gw-sig-label-side`
standardize.freqs <- opts$`standardize-frequencies`
collapse.cohorts <- opts$`collapse-cohorts`
pval.panel.space <- if(!is.null(sumstats.in)){opts$`pval-panel-space`}else{0}
pval.panel.height <- if(!is.null(sumstats.in)){opts$`pval-panel-height`}else{0}
or.panel.space <- if(!is.null(sumstats.in)){opts$`or-panel-space`}else{0}
or.panel.height <- if(!is.null(sumstats.in)){opts$`or-panel-height`}else{0}
gene.panel.space <- if(!is.null(gtf.in)){opts$`gene-panel-space`}else{0}
gene.panel.height <- if(!is.null(gtf.in)){opts$`gene-panel-height`}else{0}
pip.panel.space <- if(!is.null(pips.in)){opts$`pip-panel-space`}else{0}
pip.panel.height <- if(!is.null(pips.in)){opts$`pip-panel-height`}else{0}
credset.panel.space <- if(!is.null(pips.in)){gene.panel.space}else{0}
credset.panel.height <- if(!is.null(pips.in)){gene.panel.height}else{0}
othergene.panel.space <- if(!is.null(pips.in)){gene.panel.space}else{0}
othergene.panel.height <- if(!is.null(pips.in)){gene.panel.height}else{0}
cnv.panel.space <- opts$`cnv-panel-space`
cnv.panel.height <- opts$`cnv-panel-height`

pdf.dims <- c(opts$`pdf-height`, opts$`pdf-width`)

# Parse coordinates
parse.region.for.highlight(region, genome.in)

# Parse highlight regions
if(!is.null(highlights.str)){
  highlight.df <- data.frame(do.call("rbind",
                                     lapply(unlist(strsplit(highlights.str, split=";")),
                                            parse.region.for.highlight,
                                            genome.in=genome.in, return.values=TRUE,
                                            export.values=FALSE)))
  colnames(highlight.df) <- c("chr", "start", "end")
  highlight.df[, 2:3] <- apply(highlight.df[, 2:3], 2, as.numeric)
}else{
  highlight.df <- data.frame("chr"=NA, "start"=NA, "end"=NA)
}
highlight.start <- min(highlight.df$start)
highlight.end <- max(highlight.df$end)

# Parse case HPOs
all.case.hpos <- sort(unique(unlist(strsplit(case.hpo.str, split=";"))))

# Load sample sizes
n.samples <- get.sample.sizes.for.highlight(sample.size.table, all.case.hpos)
if(!is.na(highlight.hpo)){
  n.samples.highlight <- get.sample.sizes.for.highlight(sample.size.table,
                                                        highlight.hpo)
}else{
  n.samples.highlight <- NA
}

# Load all CNVs
cnvs <- load.cnvs.from.region.multi(cnvlist.in, region, cnv.type, all.case.hpos)

# Drop cohorts below min.cases
keep.idxs <- which(n.samples$case >= min.cases)
cnvs <- cnvs[keep.idxs]
n.samples$case <- n.samples$case[keep.idxs]
n.samples$ctrl <- n.samples$ctrl[keep.idxs]

# Collapse all cohorts into a single panel, if optioned
if(collapse.cohorts){
  cnvlist <- read.table(cnvlist.in, header=F, sep="\t")[keep.idxs, ]
  cnvs <- list(load.cnvs.from.region(cnvlist[, 2], region, cnv.type, all.case.hpos))
  n.samples$case <- sum(n.samples$case)
  n.samples$ctrl <- sum(n.samples$ctrl)
}

# Get maximum CNV frequency to be plotted, if optioned
if(standardize.freqs & !collapse.cohorts){
  max.freq <- 1.025 * max(sapply(1:length(cnvs), function(i){
    sapply(1:2, function(k){
      max.n.cnvs <- max(pileup.cnvs.for.highlight(cnvs[[i]][[k]], start=start,
                                                  end=end)$counts$count, na.rm=T)
      max.n.cnvs / n.samples[[k]][i]
    })
  }), na.rm=T)
}else{
  max.freq <- NULL
}

# Load meta-analysis summary stats, if optioned
if(!is.null(sumstats.in)){
  ss <- load.sumstats.for.region(sumstats.in, region)
}

# Load genes, if optioned
if(!is.null(gtf.in)){
  genes <- load.genes.from.gtf(gtf.in, region)
}

# Load PIPs, if optioned
if(!is.null(pips.in)){
  pips <- load.pips.for.genelist(pips.in, unique(genes$gene))
  credset.genes <- unique(pips$gene[which(!is.na(pips$credible_set))])
}

# Determine gene labeling, if optioned
if(!is.null(label.genes.str)){
  label.genes <- unique(strsplit(label.genes.str, split=";"))
}else{
  if(!is.null(pips.in)){
    label.genes <- unique(pips$gene[which(pips$PIP == max(pips$PIP))])
  }else if(length(unique(genes$gene)) <= autolabel.n.genes){
    label.genes <- sort(unique(genes$gene))
  }else{
    label.genes <- c()
  }
}

# Configure panel layout
# Upper panels = all panels below y=0 but above CNV pileups
pval.panel.y0 <- -((0.5*coord.panel.height) + coord.panel.space + (0.5*coord.panel.height))
or.panel.y0 <- pval.panel.y0 - (0.5*pval.panel.height) - or.panel.space - (0.5*or.panel.height)
gene.panel.y0 <- or.panel.y0 - (0.5*or.panel.height) - gene.panel.space - (0.5*gene.panel.height)
pip.panel.y0 <- gene.panel.y0 - (0.5*gene.panel.height) - pip.panel.space - (0.5*pip.panel.height)
credset.panel.y0 <- pip.panel.y0 - (0.5*pip.panel.height) - gene.panel.space - (0.5*gene.panel.height)
othergene.panel.y0 <- credset.panel.y0 - (0.5*gene.panel.height) - gene.panel.space - (0.5*gene.panel.height)
upper.panels.height <- -sum(c((0.5 * coord.panel.height), coord.panel.space,
                              pval.panel.space, pval.panel.height,
                              or.panel.space, or.panel.height,
                              gene.panel.space, gene.panel.height,
                              pip.panel.space, pip.panel.height,
                              credset.panel.space, credset.panel.height,
                              othergene.panel.space, othergene.panel.height))
# Lower panels = CNV pileup
cnv.panel.y0s <- upper.panels.height - sapply(1:length(cnvs), function(i){
  ((i-1)*(cnv.panel.height+cnv.panel.space)) + (0.5*cnv.panel.height)
}) + (1.1*cnv.panel.space)
total.height <- min(cnv.panel.y0s - (0.5*cnv.panel.height))
cnv.key.height <- 0.25
cnv.key.spacing <- 0.15
cnv.key.y0 <- total.height - cnv.key.spacing - (0.5*cnv.key.height)
total.height.plus.key <- cnv.key.y0 + (0.5*cnv.key.height)

# Prepare output pdf
pdf(paste(out.prefix, "locus_highlight.pdf", sep="."),
    height=pdf.dims[1]*2, width=pdf.dims[2]*2)
par(mar=c(0.5, 6, 0, 0.5), bty="n")
plot(NA, xlim=c(start, end), ylim=c(total.height.plus.key, top.panels.height),
     yaxt="n", xaxt="n", ylab="", xlab="")

# Add stick idiogram
plot.idio.stick(genome.in, chrom, idio.panel.y0, start=highlight.start,
                end=highlight.end, cytobands=cytobands.in)

# Add background shading
rect(xleft=highlight.df$start, xright=highlight.df$end,
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
plot.coord.line(start, end, coord.panel.y0, highlight.df$start, highlight.df$end,
                highlight.col=highlight.color, tick.height=0.04, vlines=TRUE,
                vlines.bottom=total.height)

# Add association summary stats, if optioned
if(!is.null(sumstats.in)){
  plot.pvalues.for.highlight(ss, y0=pval.panel.y0, panel.height=pval.panel.height,
                             cnv.type=cnv.type, gw.sig=gw.sig,
                             gw.sig.label=gw.sig.label,
                             gw.sig.label.side=gw.sig.label.side)
  plot.ors.for.highlight(ss, y0=or.panel.y0, panel.height=or.panel.height,
                         cnv.type=cnv.type)
}


# Add all genes, if optioned
if(!is.null(gtf.in)){
  plot.gene.bodies(genes, n.rows=1, y0=gene.panel.y0,
                   panel.height=gene.panel.height, col=blueblack,
                   y.axis.title="All Genes", y.axis.title.col="black",
                   label.genes=label.genes)
}


# Add PIPs and finemapped genes, if optioned
if(!is.null(pips.in)){
  plot.pips.for.highlight(pips, genes, y0=pip.panel.y0, panel.height=pip.panel.height,
                          label.genes=label.genes, col=ns.color,
                          highlight.col=graphabs.green, highlight.genes=credset.genes)
  plot.bracket(xleft=min(as.numeric(genes[which(genes$gene %in% credset.genes), "start"])),
               xright=max(as.numeric(genes[which(genes$gene %in% credset.genes), "end"])),
               y0=credset.panel.y0, height=2*credset.panel.height, col=graphabs.green,
               lwd=2, staple.wex=0.015)
  plot.gene.bodies(genes[which(genes$gene %in% credset.genes), ], n.rows=1,
                   y0=credset.panel.y0, panel.height=credset.panel.height,
                   col=graphabs.green, y.axis.title="Credible Set",
                   label.genes=label.genes)
  plot.gene.bodies(genes[which(!(genes$gene %in% credset.genes)), ], n.rows=1,
                   y0=othergene.panel.y0, panel.height=othergene.panel.height,
                   col=ns.color, y.axis.title="Other Genes")
}

# Add labels for CNV panels
cnv.names <- c("DEL"="Deletion", "DUP"="Duplication")
if(collapse.cohorts){
  cnv.title <- paste(cnv.names[cnv.type], "Evidence Across All Cohorts")
}else{
  cnv.title <- paste(cnv.names[cnv.type], "Evidence per Cohort")
}
text(x=mean(c(start, end)), y=upper.panels.height+(0.85*cnv.panel.space),
     labels=cnv.title, pos=3, cex=5.5/6)
mtext(2, at=mean(cnv.panel.y0s), line=3.5,
      text=paste(cnv.names[cnv.type], "Frequency"))

# Plot CNV panels
sapply(1:length(cnvs), function(i){
  cohort.label <- if(collapse.cohorts){NA}else{cohort.abbrevs[names(cnvs)[i]]}
  plot.cnv.panel.for.highlight(cnvs[[i]], n.case=n.samples$case[i],
                               n.ctrl=n.samples$ctrl[i], y0=cnv.panel.y0s[i],
                               cnv.type=cnv.type, highlight.hpo=highlight.hpo,
                               max.freq=max.freq, panel.height=cnv.panel.height,
                               y.axis.title=NA, expand.pheno.label=F,
                               add.cohort.label=TRUE,
                               cohort.label=cohort.label)
})

# Add key for CNV panels
plot.cnv.key.for.highlight(cnv.type, cnv.key.y0, total.n.ctrls=sum(n.samples$ctrl),
                           all.case.hpos, total.n.cases=sum(n.samples$case),
                           highlight.case.hpo=highlight.hpo,
                           total.n.cases.highlight=sum(n.samples.highlight$case),
                           panel.height=cnv.key.height)

# Close output pdf
dev.off()

