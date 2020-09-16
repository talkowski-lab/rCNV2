#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot CADM2 rCNV locus highlight for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


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
ss.in <- args$args[2]
sample.size.table <- args$args[3]
genome.in <- args$args[4]
out.prefix <- args$args[5]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# cnvlist.in <- "~/scratch/cnvlist.tsv"
# sumstats.in <- "~/scratch/HP0012638.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz"
# sample.size.table <- "~/scratch/HPOs_by_metacohort.table.tsv"
# genome.in <- "~/scratch/GRCh37.genome"
# out.prefix <- "~/scratch/noncoding_locus_highlights_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
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
start <- 84832000
end <- 86148000
region <- paste(chrom, ":", start, "-", end, sep="")
highlight.start <- 85245894
highlight.end <- 85314000
cnv.type <- "DEL"
all.case.hpos <- c("HP:0012638")
highlight.case.hpo <- "HP:0000752"

# Load sample sizes
n.samples <- get.sample.sizes(sample.size.table, all.case.hpos)
n.samples.highlight <- get.sample.sizes(sample.size.table, highlight.case.hpo)

# Load all CNVs
# Required for bedr::tabix in local Rstudio only:
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),
#                         "/Users/collins/anaconda3/envs/py3/bin",
#                         sep = ":"))
cnvs <- load.cnvs.multi(cnvlist.in, region, cnv.type, all.case.hpos)

# Load meta-analysis summary stats
ss <- load.sumstats(sumstats.in, region)

# Notes: panel ordering
# idiogram (above y=0)
# coordinate line (at y=0)
# genes (TODO)
# odds ratios (TODO)
# p-values (TODO)
# cnv pileups

# Set plot values
coord.panel.height <- 0.2
coord.panel.space <- 0.2
coord.panel.y0 <- 0
idio.panel.height <- 0.15
idio.panel.space <- 0.3
idio.panel.y0 <- (0.5 * coord.panel.height) + idio.panel.space + (0.5*idio.panel.height)
# Top panels = those at or above y=0 (idiogram & coordinate line)
top.panels.height <- sum(c(idio.panel.height, idio.panel.space, 0.5*coord.panel.height))
# Upper panels = all panels below y=0 but above CNV pileups
upper.panels.height <- -sum(c((0.5 * coord.panel.height), coord.panel.space))
# Lower panels = CNV pileup
cnv.panel.height <- 0.85
cnv.panel.spacing <- 0.25
cnv.panel.y0s <- upper.panels.height - sapply(1:4, function(i){((i-1)*(cnv.panel.height+cnv.panel.spacing)) + (0.5*cnv.panel.height)})
total.height <- min(cnv.panel.y0s - (0.5*cnv.panel.height))
cnv.key.height <- 0.25
cnv.key.spacing <- 0.15
cnv.key.y0 <- total.height - cnv.key.spacing - (0.5*cnv.key.height)
total.height.plus.key <- cnv.key.y0 + (0.5*cnv.key.height)

# Prep locus plot
pdf(paste(out.prefix, "locus_highlights.CADM2.pdf", sep="."), height=3*2, width=3*2)
par(mar=c(0.5, 4, 0, 0.5), bty="n")
plot(NA, xlim=c(start, end), ylim=c(total.height.plus.key, top.panels.height), 
     yaxt="n", xaxt="n", ylab="", xlab="")

# Add background shading
rect(xleft=par("usr")[1], xright=par("usr")[2], 
     ybottom=cnv.panel.y0s+(0.5*cnv.panel.height),
     ytop=cnv.panel.y0s-(0.5*cnv.panel.height),
     border=NA, bty="n", col=bluewhite)
rect(xleft=highlight.start, xright=highlight.end, 
     ybottom=total.height, ytop=0, 
     col=adjustcolor(highlight.color, alpha=0.3), border=NA)

# Add stick idiogram
add.idio.stick(genome.in, chrom, highlight.start, highlight.end, idio.panel.y0,
               tick.height=-0.03)

# Add coordinate line
add.coord.line(start, end, coord.panel.y0, highlight.start, highlight.end, 
               highlight.col=highlight.color, tick.height=0.04)

# Add CNV panels
sapply(1:4, function(i){
  add.cnv.panel(cnvs[[i]], n.case=n.samples$case[i], n.ctrl=n.samples$ctrl[i], 
                highlight.hpo=highlight.case.hpo, y0=cnv.panel.y0s[i], 
                cnv.type=cnv.type, panel.height=cnv.panel.height, 
                y.axis.title="Carrier\nFreq.", expand.pheno.label=F,
                case.legend.side="left", ctrl.legend.side="left",
                add.cohort.label=TRUE, cohort.label=paste("Cohort", i))
})

# Add key for CNV panels
add.cnv.key(cnv.type, cnv.key.y0, total.n.ctrls=sum(n.samples$ctrl), 
            all.case.hpos, total.n.cases=sum(n.samples$case),
            highlight.case.hpo, total.n.cases.highlight=sum(n.samples.highlight$case),
            panel.height=cnv.key.height)

# Close pdf
dev.off()

