#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform burden test for rare CNVs in genes adjusted for gene-level covariates


options(scipen=1000, stringsAsFactors=F)


# Extract sample counts from table
get.sample.counts <- function(pheno.table.in, cohort.name, 
                              case.hpo, control.hpo){
  ptab <- read.table(pheno.table.in, header=T, 
                     sep="\t", comment.char="")
  if(!any(colnames(ptab)==cohort.name)){
    stop(paste("Cohort \"", cohort.name, 
               "\" cannot be found in header of ",
               pheno.table.in, sep=""))
  }
  case.n <- ptab[which(ptab[,1] == case.hpo),
                 which(colnames(ptab) == cohort.name)]
  control.n <- ptab[which(ptab[,1] == control.hpo),
                    which(colnames(ptab) == cohort.name)]
  return(list("case.n"=as.numeric(case.n),
              "control.n"=as.numeric(control.n)))
}


# Load an input bed file of gene features
import.features <- function(features.in){
  # Load raw features
  features <- read.table(features.in, sep="\t", header=T, comment.char="")
  colnames(features)[1] <- "chr"
  features[, 5:ncol(features)] <- apply(features[, 5:ncol(features)], 2, as.numeric)
  return(features)
}


# Process an input bed file of genes and counts
import.bed <- function(bed.in, features, 
                       case.col.name, case.n, 
                       control.col.name, control.n,
                       neighbor.dist=1000000){
  bed <- read.table(bed.in, sep="\t", header=T, comment.char="")
  
  case.col.idx <- which(colnames(bed)==case.col.name)
  if(length(case.col.idx) == 0){
    stop(paste("--case-column \"", case.col.name, 
               "\" cannot be found in BED header", sep=""))
  }
  wcase.col.idx <- which(colnames(bed)==paste(case.col.name, 
                                              weighted.suffix, sep=""))
  if(length(case.col.idx) == 0){
    stop(paste("Weighted case CNV column \"", case.col.name, 
               weighted.suffix, "\" cannot be found in BED header", 
               sep=""))
  }
  
  control.col.idx <- which(colnames(bed)==control.col.name)
  if(length(control.col.idx) == 0){
    stop(paste("--control-column \"", control.col.name, 
               "\" cannot be found in BED header", sep=""))
  }
  wcontrol.col.idx <- which(colnames(bed)==paste(control.col.name, 
                                                 weighted.suffix, sep=""))
  if(length(control.col.idx) == 0){
    stop(paste("Weighted control CNV column \"", control.col.name, 
               weighted.suffix, "\" cannot be found in BED header", 
               sep=""))
  }
  
  bed <- bed[,c(1:4, case.col.idx, control.col.idx, 
                wcase.col.idx, wcontrol.col.idx)]
  colnames(bed) <- c("chr", "start", "end", "gene",
                     "case.CNV", "control.CNV", 
                     "case.CNV.w", "control.CNV.w")
  bed$case.ref <- case.n - bed$case.CNV
  bed$control.ref <- control.n - bed$control.CNV
  bed$case.CNV.freq <- bed$case.CNV / case.n
  bed$control.CNV.freq <- bed$control.CNV / control.n
  bed$case.CNV.w.norm <- bed$case.CNV.w / case.n
  bed$control.CNV.w.norm <- bed$control.CNV.w / control.n
  bed$all.CNV <- bed$case.CNV + bed$control.CNV
  bed$all.CNV.freq <- bed$all.CNV / (case.n + control.n)
  bed$all.CNV.w <- bed$case.CNV.w + bed$control.CNV.w
  bed$all.CNV.w.norm <- bed$all.CNV.w / (case.n + control.n)
  
  # Add regional annotation of average # of CNVs per 
  # neighboring genes within a prespecified distance
  bed$avg_neighbor.all.CNV <- sapply(1:nrow(bed), function(i){
    counts <- bed$all.CNV[which(bed$chr == bed$chr[i]
                                & bed$start <= bed$end[i] + neighbor.dist
                                & bed$end >= bed$start[i] - neighbor.dist
                                & bed$gene != bed$gene[i])]
    if(length(counts) > 0){
      mean(counts)
    }else{
      return(0)
    }
  })
  
  # Append gene features
  bed <- merge(bed, features, by=c("chr", "start", "end", "gene"),
               all=F, sort=F, suffixes=c(".cnvs", ".features"))
  
  return(bed)
}


# Predict number of case CNVs per gene
predict.case.cnvs.singleChrom <- function(bed, test.chrom, pred.colName="case.CNV.w.norm"){
  # Clean data
  cols_to_keep <- c("control.CNV.w.norm")
  df.sub <- bed[, c(which(colnames(bed)==pred.colName),
                    which(colnames(bed) %in% cols_to_keep),
                    grep("eigenfeature", colnames(bed)))]
  colnames(df.sub)[1] <- "response"
  df.train <- df.sub[which(bed$chr!=test.chrom), ]
  df.test <- df.sub[which(bed$chr==test.chrom), ]
  
  # Fit model
  fit <- glm(response ~ ., data=df.train)
  
  # Predict on unseen data
  pred.vals <- predict.glm(fit, newdata=df.test)
  names(pred.vals) <- bed$gene[which(bed$chr==test.chrom)]
  return(pred.vals)
}
predict.case.cnvs <- function(bed, pred.colName="case.CNV.w.norm"){
  pred.vals <- unlist(sapply(unique(bed$chr), function(test.chrom){
    predict.case.cnvs.singleChrom(bed, test.chrom, pred.colName)
  }))
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)
require(MASS, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--pheno-table"), type="character", default=NULL,
              help="table with counts of samples per HPO term per cohort [default %default]",
              metavar="file"),
  make_option(c("--cohort-name"), type="character", default=NULL,
              help="name of cohort [default %default]",
              metavar="string"),
  make_option(c("--cnv"), type="character", default=NULL,
              help="CNV class to test (options: 'CNV', 'DEL', 'DUP') [default %default]",
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
  make_option(c("--precision"), type="integer", default=6,
              help="level of precision for floats [default %default]",
              metavar="integer"),
  make_option(c("-z", "--bgzip"), action="store_true", default=FALSE,
              help="bgzip output BED file [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog cnvcounts genefeatures outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
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
features.in <- args$args[2]
outfile <- args$args[3]
pheno.table.in <- opts$`pheno-table`
cohort.name <- opts$`cohort-name`
cnv <- opts$cnv
null.table.in <- opts$`null-table-in`
null.model <- opts$`null-model`
case.hpo <- opts$`case-hpo`
control.hpo <- opts$`control-hpo`
case.col.name <- opts$`case-column`
control.col.name <- opts$`control-column`
weighted.suffix <- opts$`weighted-suffix`
precision <- opts$precision

# DEV PARAMETERS:
bed.in <- "~/scratch/gene_burden_test/counts/meta2.HP0012638.uCNV.DEL.gene_burden.counts.bed.gz"
features.in <- "~/scratch/gene_burden_test/gencode.v19.canonical.genomic_features.eigenfeatures.bed.gz"
outfile <- "~/scratch/gene_burden_test/meta2.HP0012638.uCNV.DEL.gene_burden_stats.bed.gz"
pheno.table.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
cohort.name <- "meta2"
case.hpo <- "HP:0012638"
control.hpo <- "HEALTHY_CONTROL"
case.col.name <- "case_cnvs"
control.col.name <- "control_cnvs"
weighted.suffix <- "_weighted"
precision <- 6

# Extract sample counts
sample.counts <- get.sample.counts(pheno.table.in, cohort.name, 
                                   case.hpo, control.hpo)

# Read features
features <- import.features(features.in)

# Process input BED
bed <- import.bed(bed.in, features, case.col.name, sample.counts$case.n,
                  control.col.name, sample.counts$control.n)

## DEV EXPLORATORY ANALYSES
for(i in 19:ncol(bed)){
  if(length(unique(bed[, i])) > 1){
    plot(bed[, i], bed$all.CNV, main=colnames(bed)[i], cex=0.5,
         xlab=colnames(bed)[i], ylab="All CNVs")
    abline(lm(bed$all.CNV ~ bed[, i]), col="red")
    text(x=mean(par("usr")[1:2]), y=par("usr")[4], pos=1, col="red",
         labels=cor(bed[, i], bed$all.CNV))
  }
}

# Proxy code:
bed$other_genes_within_1Mb <- sapply(1:nrow(bed), function(i){
  length(which(bed$chr==bed$chr[i] & 
                 (abs(bed$end[i] - bed$start) <= 1000000 | abs(bed$end - bed$start[i]) <= 1000000))) - 1
})
bed$mean_control_CNVs_nearby <- sapply(1:nrow(bed), function(i){
  m <- mean(bed$control.CNV[setdiff(which(bed$chr==bed$chr[i] & 
                                            (abs(bed$end[i] - bed$start) <= 1000000 | abs(bed$end - bed$start[i]) <= 1000000)),
                                    i)])
  if(is.na(m)){
    m <- 0
  }
  return(m)
})
bed$mean_case_CNVs_nearby <- sapply(1:nrow(bed), function(i){
  m <- mean(bed$case.CNV[setdiff(which(bed$chr==bed$chr[i] & 
                                         (abs(bed$end[i] - bed$start) <= 1000000 | abs(bed$end - bed$start[i]) <= 1000000)),
                                 i)])
  if(is.na(m)){
    m <- 0
  }
  return(m)
})

fit <- glm(case.CNV.w ~ control.CNV.w + eigenfeature_1 + eigenfeature_2 + eigenfeature_3 + eigenfeature_4 + eigenfeature_5 + eigenfeature_6 + eigenfeature_7 + eigenfeature_8 + eigenfeature_9 + other_genes_within_1Mb + mean_case_CNVs_nearby, data=bed)
pred.case <- predict.glm(fit, data=bed, type="response")

# IDEA:
# predict: n_case_cnvs ~ n_control_cnvs_in_gene + n_case_cnvs_in_1Mb_region + gene_pcs1-N
# then:
# chisq test based on expected # case CNVs

# Maybe add features for:
# number of genes within 1Mb
# distance to nearest gene
# average CNV count for genes within 1Mb

# May not have to use weighted case CNV counts if factoring in features above (hopefully will capture segmental associations this way)

# For each chromosome:
# Train on all other chroms
# Predict on test chrom

# Maybe add features for:
# number of genes within 1Mb
# distance to nearest gene
# average CNV count for genes within 1Mb

# # Format output
# colnames(burden.bed)[1] <- "#chr"
# if(length(grep(".gz", outfile, fixed=T)) > 0){
#   outfile <- tools::file_path_sans_ext(outfile)
# }
# write.table(burden.bed, outfile, sep="\t", col.names=T, row.names=F, quote=F)
# if(opts$bgzip == TRUE){
#   system(paste("bgzip -f", outfile),wait=T)
# }
