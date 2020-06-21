#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform a Fisher's exact test for rare CNV burden for cis-regulatory blocks (CRBs)


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


# Apply empirical continuity correction to 2x2 contingency table
# Per Sweeting et al., Stat. Med., 2004 (section 3.3)
sweeting.correction <- function(control.ref, case.ref, control.alt, case.alt, cc.sum=0.01){
  # Count number of carriers & non-carriers
  n.alt <- control.alt + case.alt
  n.ref <- control.ref + case.ref
  # Require at least one CNV to be observed
  if(n.alt>0){
    nt <- n.alt
    R <- n.ref/n.alt
    # Preliminary odds ratio estimate
    nonzero.case.odds <- case.alt/case.ref
    nonzero.control.odds <- control.alt/control.ref
    if(nonzero.control.odds>0){
      ohat <- nonzero.case.odds/nonzero.control.odds
      # Apply standard continuity correction of 0.5 to pooled estimate if no CNVs observed in controls  
    }else{
      nonzero.case.odds <- (case.alt+(cc.sum/2))/(case.ref+(cc.sum/2))
      nonzero.control.odds <- (control.alt+(cc.sum/2))/(control.ref+(cc.sum/2))
      ohat <- nonzero.case.odds/nonzero.control.odds
    }
    # Solve for kc & kt
    kc <- R/(R+ohat)
    kt <- ohat/(R+ohat)
    # Compute continuity corrections
    cor.case_alt <- cc.sum * kt
    cor.case_ref <- cc.sum * kc
    cor.control_alt <- cc.sum * (nt + kt)
    cor.control_ref <- cc.sum * ((nt*R) + kc)
    # Apply continuity corrections
    case.alt <- case.alt + cor.case_alt
    case.ref <- case.ref + cor.case_ref
    control.alt <- control.alt + cor.control_alt
    control.ref <- control.ref + cor.control_ref
    return(list("control.ref"=control.ref, 
                "case.ref"=case.ref, 
                "control.alt"=control.alt, 
                "case.alt"=case.alt))
  }else{
    return(list("control.ref"=NA, 
                "case.ref"=NA, 
                "control.alt"=NA, 
                "case.alt"=NA))
  }
}


# Estimate odds ratio & CI with Sweeting correction (as necessary)
calc.or <- function(control.ref, case.ref, control.alt, case.alt, 
                    conf=0.95, cc.sum=0.01){
  corrected <- sweeting.correction(control.ref, case.ref, control.alt, 
                                   case.alt, cc.sum)
  case.odds <- corrected$case.alt/corrected$case.ref
  control.odds <- corrected$control.alt/corrected$control.ref
  or <- case.odds/control.odds
  log.or <- log(or)
  se <- sqrt(sum(unlist(lapply(corrected, function(x){1/x}))))
  log.or.lower <- log.or - qnorm((1-conf)/2, lower.tail=F)*se
  log.or.upper <- log.or + qnorm(0.95)*se
  return(c(log.or, log.or.lower, log.or.upper))
}


# Process an input bed file of CNVs
import.bed <- function(bed.in, case.col.name, case.n, 
                       control.col.name, control.n){
  bed <- read.table(bed.in, sep="\t", header=T, comment.char="")
  
  case.col.idx <- which(colnames(bed)==case.col.name)
  if(length(case.col.idx) == 0){
    stop(paste("--case-column \"", case.col.name, 
               "\" cannot be found in BED header", sep=""))
  }

  control.col.idx <- which(colnames(bed)==control.col.name)
  if(length(control.col.idx) == 0){
    stop(paste("--control-column \"", control.col.name, 
               "\" cannot be found in BED header", sep=""))
  }

  bed <- bed[,c(1:4, case.col.idx, control.col.idx)]
  colnames(bed) <- c("chr", "start", "end", "crb_id",
                     "case.CNV", "control.CNV")
  bed$case.ref <- case.n - bed$case.CNV
  bed$control.ref <- control.n - bed$control.CNV
  bed$case.CNV.freq <- bed$case.CNV / case.n
  bed$control.CNV.freq <- bed$control.CNV / control.n

  # Calculate odds ratio
  or.df <- as.data.frame(t(sapply(1:nrow(bed), function(i){
    calc.or(control.ref=bed$control.ref[i],
            case.ref=bed$case.ref[i],
            control.alt=bed$control.CNV[i],
            case.alt=bed$case.CNV[i],
            conf=0.95, cc.sum=0.01)
  })))
  colnames(or.df) <- c("ln.OR", "ln.OR.lower", "ln.OR.upper")
  bed <- as.data.frame(cbind(bed, or.df))
  
  return(bed)
}


# Fisher's exact test for a single vector of CNV counts
burden.test.single <- function(counts){
  case.cnv <- as.integer(counts[1])
  control.cnv <- as.integer(counts[2])
  case.ref <- as.integer(counts[3])
  control.ref <- as.integer(counts[4])
  
  if(case.cnv == 0 & control.cnv == 0){
    p <- 1
    or <- c(NA,NA,NA)
  }else{
    cnv.mat <- matrix(c(control.ref, case.ref, control.cnv, case.cnv),
                      byrow=T, nrow=2)
    
    p <- fisher.test(cnv.mat, alternative="greater")$p.value
    or <- fisher.test(cnv.mat)
    or <- c(or$estimate, or$conf.int)
  }
  
  f.res <- as.numeric(c(p, or))
  names(f.res) <- c("p", "OR", "OR.lower", "OR.upper")
  return(f.res)
}


# Compute Fisher's exact test lookup table
build.fisher.lookup.table <- function(bed){
  counts.df <- unique(bed[, which(colnames(bed) %in% c("case.CNV", "control.CNV",
                                                        "case.ref", "control.ref"))])
  counts.df <- counts.df[with(counts.df, order(control.CNV, case.CNV)),]
  
  f.stats <- t(apply(counts.df, 1, burden.test.single))
  
  f.table <- cbind(counts.df, f.stats)
  return(f.table)
}


# Fisher's exact test of case:control CNV burden per bin
burden.test <- function(bed, precision){
  f.table <- build.fisher.lookup.table(bed)
  
  f.res <- merge(bed, f.table, sort=F, all.x=T, all.y=F,
                 c("case.CNV", "control.CNV",
                   "case.ref", "control.ref"))
  f.res <- f.res[with(f.res, order(chr, start)), ]
  
  fisher.bed <- data.frame("chr" = f.res$chr,
                           "start" = f.res$start, 
                           "end" = f.res$end,
                           "crb_id" = f.res$crb_id,
                           "case_alt" = f.res$case.CNV,
                           "case_ref" = f.res$case.ref,
                           "case_freq" = round(f.res$case.CNV.freq, precision),
                           "control_alt" = f.res$control.CNV, 
                           "control_ref" = f.res$control.ref,
                           "control_freq" = round(f.res$control.CNV.freq, precision),
                           "fisher_phred_p" = round(-log10(f.res$p), precision),
                           "fisher_OR" = round(f.res$OR, precision), 
                           "fisher_OR_lower" = round(f.res$OR.lower, precision), 
                           "fisher_OR_upper" = round(f.res$OR.upper, precision))
  return(fisher.bed)
}

#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)

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
  make_option(c("--precision"), type="integer", default=6,
              help="level of precision for floats [default %default]",
              metavar="integer"),
  make_option(c("-z", "--bgzip"), action="store_true", default=FALSE,
              help="bgzip output BED file [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog crb_counts outfile",
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
cnv <- opts$cnv
case.hpo <- opts$`case-hpo`
control.hpo <- opts$`control-hpo`
case.col.name <- opts$`case-column`
control.col.name <- opts$`control-column`
precision <- opts$precision

# # Dev parameters:
# bed.in <- "~/scratch/meta2.UNKNOWN.rCNV.loose_noncoding.DEL.crb_burden.counts.bed.gz"
# outfile <- "~/scratch/meta2.UNKNOWN.rCNV.loose_noncoding.DEL.crb_burden.stats.bed"
# pheno.table.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# cohort.name <- "meta2"
# cnv <- "DEL"
# case.hpo <- "UNKNOWN"
# control.hpo <- "HEALTHY_CONTROL"
# case.col.name <- "case_cnvs"
# control.col.name <- "control_cnvs"
# precision <- 6

# Extract sample counts
sample.counts <- get.sample.counts(pheno.table.in, cohort.name, 
                                   case.hpo, control.hpo)

# Process input BED
bed <- import.bed(bed.in, case.col.name, sample.counts$case.n,
                  control.col.name, sample.counts$control.n)

# Run burden tests
fisher.bed <- burden.test(bed, precision)

# Format output
colnames(fisher.bed)[1] <- "#chr"
if(length(grep(".gz", outfile, fixed=T)) > 0){
  outfile <- tools::file_path_sans_ext(outfile)
}
write.table(fisher.bed, outfile, sep="\t", col.names=T, row.names=F, quote=F)
if(opts$bgzip == TRUE){
  system(paste("bgzip -f", outfile),wait=T)
}
