#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform a Fisher's exact test for rare CNV burden for genes


options(scipen=1000, stringsAsFactors=F)


# Extract sample counts from table
get.sample.counts <- function(pheno.table.in, cohort.name, case.hpo, control.hpo){
  ptab <- read.table(pheno.table.in, header=T, sep="\t", comment.char="")
  if(!any(colnames(ptab)==cohort.name)){
    stop(paste("Cohort \"", cohort.name, "\" cannot be found in header of ",
               pheno.table.in, sep=""))
  }
  case.n <- ptab[which(ptab[,1] == case.hpo),
                 which(colnames(ptab) == cohort.name)]
  control.n <- ptab[which(ptab[,1] == control.hpo),
                    which(colnames(ptab) == cohort.name)]
  return(list("case.n"=as.numeric(case.n),
              "control.n"=as.numeric(control.n)))
}


# Process an input bed file of CNVs
import.bed <- function(bed.in, case.col.name, case.n, control.col.name, control.n){
  bed <- read.table(bed.in, sep="\t", header=T, comment.char="")

  case.col.idx <- which(colnames(bed)==case.col.name)
  if(length(case.col.idx) == 0){
    stop(paste("--case-column \"", case.col.name, "\" cannot be found ",
               "in BED header", sep=""))
  }
  wcase.col.idx <- which(colnames(bed)==paste(case.col.name, weighted.suffix, sep=""))
  if(length(case.col.idx) == 0){
    stop(paste("Weighted case CNV column \"", case.col.name, weighted.suffix,
               "\" cannot be found in BED header", sep=""))
  }
  
  control.col.idx <- which(colnames(bed)==control.col.name)
  if(length(control.col.idx) == 0){
    stop(paste("--control-column \"", control.col.name, "\" cannot be found ",
               "in BED header", sep=""))
  }
  wcontrol.col.idx <- which(colnames(bed)==paste(control.col.name, weighted.suffix, sep=""))
  if(length(control.col.idx) == 0){
    stop(paste("Weighted control CNV column \"", control.col.name, weighted.suffix,
               "\" cannot be found in BED header", sep=""))
  }
  
  bed <- bed[,c(1:4, case.col.idx, control.col.idx, wcase.col.idx, wcontrol.col.idx)]
  colnames(bed) <- c("chr", "start", "end", "gene",
                     "case.CNV", "control.CNV", 
                     "case.CNV.w", "control.CNV.w")
  bed$case.ref <- case.n - bed$case.CNV
  bed$control.ref <- control.n - bed$control.CNV
  bed$case.CNV.freq <- bed$case.CNV / case.n
  bed$control.CNV.freq <- bed$control.CNV / control.n
  bed$case.CNV.w.norm <- bed$case.CNV.w / case.n
  bed$control.CNV.w.norm <- bed$control.CNV.w / control.n

  return(bed)
}


# # Binned Z-test case:control CNV burden per gene
# bin.burden.test <- function(bed, precision, bins=50){
#   # Step 1: calculate case proportions and quantiles for all genes
#   # Note: this also drops all double-zero genes
#   bed$total.CNV.w <- bed$case.CNV.w + bed$control.CNV.w
#   bed$delta.total <- bed$case.CNV.w - bed$control.CNV.w
#   bed$delta.norm <- 100*(bed$case.CNV.w.norm - bed$control.CNV.w.norm)
#   nonzero.idx <- which(bed$case.CNV + bed$control.CNV > 0)
#   tbed <- bed[nonzero.idx, ]
#   tbed$case.prop <- tbed$case.CNV.w / tbed$total.CNV.w
#   tbed$quantile <- floor(bins*rank(tbed$total.CNV.w)/nrow(tbed)) + 1
#   tbed$quantile[which(tbed$quantile==max(tbed$quantile))] <- max(tbed$quantile) - 1
#   
#   # Step 2: fit null distribution for total difference in normalized, weighted case:control counts
#   # Fit separately for each quantile
#   nonzero.idx <- which(tbed$case.CNV + tbed$control.CNV > 0)
#   more.controls.idx <- intersect(which(tbed$control.CNV.freq >= tbed$case.CNV.freq),
#                                  nonzero.idx)
#   sds <- sapply(1:bins, function(q){
#     null.vals.oneside <- abs(tbed$delta.norm[intersect(which(tbed$quantile==q), more.controls.idx)])
#     null.vals <- c(-null.vals.oneside, null.vals.oneside)
#     sd(null.vals, na.rm=T)
#   })
#   
#   #Working code
#   q <- 50
#   hist(tbed$delta.norm[which(tbed$quantile==q)], breaks=50, freq=F)
#   curve(dnorm(x, sd=sds[q]), add=T, col="red", lwd=2)
#   
#   # Step 3: calculate Z-score and p-value for each gene depending on quantile
#   sapply(1:nrow(tbed), function(i){
#     
#   })
# }


# Z-test for case:control CNV burden per gene
z.burden.test <- function(bed, use.unweighted.controls, min.total.CNV.w=1, precision){
  # Step 1: fit normal distribution to mirrored distribution of genes with more CNVs in controls than cases
  # Approach adopted from pLoF constraint calculations per Lek et al., Nature, 2016
  bed$total.CNV.w <- bed$case.CNV.w + bed$control.CNV.w
  nonzero.idx <- which(bed$case.CNV + bed$control.CNV > 0)
  more.controls.idx <- intersect(which(bed$control.CNV.freq >= bed$case.CNV.freq),
                                 nonzero.idx)
  if(use.unweighted.controls==T){
    bed$delta.total <- bed$case.CNV.w - bed$control.CNV
    bed$delta.norm <- 100*(bed$case.CNV.w.norm - bed$control.CNV.freq)
  }else{
    bed$delta.total <- bed$case.CNV.w - bed$control.CNV.w
    bed$delta.norm <- 100*(bed$case.CNV.w.norm - bed$control.CNV.w.norm)
  }
  null.vals.oneside <- abs(bed$delta.norm[more.controls.idx])
  # null.vals.oneside.capped <- null.vals.oneside[which(null.vals.oneside < quantile(null.vals.oneside, 0.995))]
  null.vals <- c(-null.vals.oneside, null.vals.oneside)
  null.mean <- mean(null.vals)
  null.sd <- sd(null.vals)

  # Step 2: compute Z-score and p-value for each gene according to null distribution
  bed$Zscore <- bed$delta.norm / null.sd
  bed$pvalue <- pnorm(bed$Zscore, lower.tail=F)
  
  # Step 3: overwrite all p-values below minimum case CNV count
  bed$pvalue[which(bed$total.CNV.w < min.total.CNV.w)] <- NA
  
  burden.bed <- data.frame("chr" = bed$chr,
                           "start" = bed$start,
                           "end" = bed$end,
                           "gene" = bed$gene,
                           "case_cnvs" = bed$case.CNV,
                           "case_freq" = round(bed$case.CNV.freq, precision),
                           "case_cnvs_weighted" = round(bed$case.CNV.w, precision),
                           "case_freq_weighted" = round(bed$case.CNV.w.norm, precision),
                           "control_cnvs" = bed$control.CNV,
                           "control_freq" = round(bed$control.CNV.freq, precision),
                           "control_cnvs_weighted" = round(bed$control.CNV.w, precision),
                           "control_freq_weighted" = round(bed$control.CNV.w.norm, precision),
                           "z_score" = round(bed$Zscore, precision),
                           "phred_p" = round(-log10(bed$pvalue), precision))
  return(burden.bed)
}


#################
### RSCRIPT BLOCK
#################
require(optparse,quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--pheno-table"), type="character", default=NULL,
              help="table with counts of samples per HPO term per cohort [default %default]",
              metavar="file"),
  make_option(c("--cohort-name"), type="character", default=NULL,
              help="name of cohort [default %default]",
              metavar="string"),
  make_option(c("--case-hpo"), type="character", default=NULL,
              help="HPO term to use for case samples [default %default]",
              metavar="string"),
  make_option(c("--control-hpo"), type="character", default='HEALTHY_CONTROL',
              help="HPO term to use for control samples [default %default]",
              metavar="string"),
  make_option(c("--case-column"), type="character", default='case_cnvs',
              help="name of column to use for raw (unweighted) case CNV counts [default %default]",
              metavar="string"),
  make_option(c("--control-column"), type="character", default='control_cnvs',
              help="name of column to use for raw (unweighted) control CNV counts [default %default]",
              metavar="string"),
  make_option(c("--weighted-suffix"), type="character", default='_weighted',
              help="suffix appended to the end of CNV count columns for weighted counts [default %default]",
              metavar="string"),
  make_option(c("--unweighted-controls"), action="store_true", default=FALSE,
              help="use unweighted control CNV counts for burden testing [default %default]"),
  make_option(c("--min-weighted-cnvs"), type="numeric", default=1,
              help="do not compute p-value for genes with fewer than N weighted CNVs [default %default]",
              metavar="numeric"),
  make_option(c("--precision"), type="integer", default=6,
              help="level of precision for floats [default %default]",
              metavar="integer"),
  make_option(c("-z", "--bgzip"), action="store_true", default=FALSE,
              help="bgzip output BED file [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog genecounts outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}
if(is.null(opts$`pheno-table`)){
  stop("Must provide --pheno-table\n")
}
if(is.null(opts$`cohort-name`)){
  stop("Must specify --cohort-name\n")
}
if(is.null(opts$`case-hpo`)){
  stop("Must specify --case-hpo\n")
}

# Writes args & opts to vars
bed.in <- args$args[1]
outfile <- args$args[2]
pheno.table.in <- opts$`pheno-table`
cohort.name <- opts$`cohort-name`
case.hpo <- opts$`case-hpo`
control.hpo <- opts$`control-hpo`
case.col.name <- opts$`case-column`
control.col.name <- opts$`control-column`
weighted.suffix <- opts$`weighted-suffix`
min.total.CNV.w <- opts$`min-weighted-cnvs`
use.unweighted.controls <- opts$`unweighted-controls`
precision <- opts$precision

# # DEV PARAMETERS:
# bed.in <- "~/scratch/meta1.HP0012759.uCNV.DEL.gene_burden.counts.bed.gz"
# outfile <- "~/scratch/test_gene_burden.stats.txt"
# pheno.table.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# cohort.name <- "meta1"
# case.hpo <- "HP:0012759"
# control.hpo <- "HEALTHY_CONTROL"
# case.col.name <- "case_cnvs"
# control.col.name <- "control_cnvs"
# weighted.suffix <- "_weighted"
# use.unweighted.controls <- F
# min.total.CNV.w <- 1
# precision <- 6

# Extract sample counts
sample.counts <- get.sample.counts(pheno.table.in, cohort.name, case.hpo, control.hpo)

# Process input BED
bed <- import.bed(bed.in, case.col.name, sample.counts$case.n,
                  control.col.name, sample.counts$control.n)

# Run burden tests
burden.bed <- z.burden.test(bed, use.unweighted.controls, min.total.CNV.w, precision)

# Format output
colnames(burden.bed)[1] <- "#chr"
if(length(grep(".gz", outfile, fixed=T)) > 0){
  outfile <- tools::file_path_sans_ext(outfile)
}
write.table(burden.bed, outfile, sep="\t", col.names=T, row.names=F, quote=F)
if(opts$bgzip == TRUE){
  system(paste("bgzip -f", outfile),wait=T)
}
