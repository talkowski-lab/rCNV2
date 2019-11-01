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


# Conservatively adjust weighted CNV counts to lower bound of confidence interval
adjust.counts <- function(counts, samples, ci=0.9){
  # nz.idx <- which(counts>0)
  # nz.counts <- counts[nz.idx]
  sigma <- sd(counts)
  ci.z <- qnorm((1-ci)/2, lower.tail=F)
  ci.mod <- ci.z * (sigma / sqrt(samples))
  sapply(counts - ci.mod, function(x){max(c(0, x))})
}


# Process an input bed file of CNVs
import.bed <- function(bed.in, case.col.name, case.n, 
                       control.col.name, control.n,
                       adjust.weighted.counts=F){
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
  if(adjust.weighted.counts==T){
    bed$case.CNV.w.adj <- adjust.counts(bed$case.CNV.w, case.n)
    bed$case.CNV.w.norm <- bed$case.CNV.w.adj / case.n
    bed$control.CNV.w.adj <- adjust.counts(bed$control.CNV.w, control.n)
    bed$control.CNV.w.norm <- bed$control.CNV.w.adj / control.n
  }else{
    bed$case.CNV.w.norm <- bed$case.CNV.w / case.n
    bed$control.CNV.w.norm <- bed$control.CNV.w / control.n
  }
  
  # Calculate odds ratio
  or.df <- as.data.frame(t(sapply(1:nrow(bed), function(i){
    calc.or(control.ref=control.n-bed$control.CNV.w[i],
            case.ref=case.n-bed$case.CNV.w[i],
            control.alt=bed$control.CNV.w[i],
            case.alt=bed$case.CNV.w[i],
            conf=0.95, cc.sum=0.01)
  })))
  colnames(or.df) <- c("ln.OR", "ln.OR.lower", "ln.OR.upper")
  bed <- cbind(bed, or.df)
  
  return(bed)
}


# Plot null distribution
plot.null.normal <- function(null.vals, null.mean, null.sd, cohort.name){
  plot.lims <- quantile(null.vals, c(0.001, 0.999))
  par(mar=c(4, 4, 2.5, 1), bty="n")
  hist(null.vals, col="gray90", freq=F, main="", xlim=plot.lims, breaks=50,
       xlab=bquote(Delta * "(Case, Control; %)"))
  curve(dnorm(x, null.mean, null.sd), add=T, lwd=2)
  mtext(3, line=0.2, font=2,
        text=paste("Gaussian null fit for\n",
                   prettyNum(length(null.vals)/2, big.mark=","),
                   "genes from", cohort.name))
}
plot.null.exp <- function(null.vals.oneside, null.exp.rate, cohort.name){
  par(mar=c(4, 4, 2.5, 1), bty="n")
  hist(null.vals.oneside, col="gray90", freq=F, main="", breaks=50,
       xlim=c(0, quantile(null.vals.oneside, 0.999)),
       xlab=bquote(Delta * "(Case, Control; %)"))
  curve(dexp(x, null.exp.rate), add=T, lwd=2)
  mtext(3, line=0.2, font=2,
        text=paste("Exponential null fit for\n",
                   prettyNum(length(null.vals.oneside), big.mark=","),
                   "genes from", cohort.name))
}
plot.null <- function(null.vals, null.mean, null.sd, null.exp.rate, cohort.name){
  par(mfrow=c(1, 2))
  plot.null.normal(null.vals, null.mean, null.sd, cohort.name)
  plot.null.exp(null.vals[which(null.vals>=0)], null.exp.rate, cohort.name)
}


# Build null distribution of CNVs in cases vs controls
build.null <- function(bed, use.unweighted.controls, min.total.CNV=NULL, 
                       precision, null.dist.plot.out, cohort.name){
  # Fit normal distribution to mirrored distribution of delta.norm 
  # for genes with at least as many CNVs in controls than cases
  # Approach adopted from pLoF constraint calculations per Lek et al., Nature, 2016
  bed$total.CNV <- bed$case.CNV + bed$control.CNV
  bed$total.CNV.w <- bed$case.CNV.w + bed$control.CNV.w
  nonzero.idx <- which(bed$total.CNV > 0)
  if(is.null(min.total.CNV)){
    min.total.CNV <- floor(quantile(bed$total.CNV[which(bed$total.CNV < quantile(bed$total.CNV, 0.95, na.rm=T))], 0.25))
    min.total.CNV <- max(c(1, min.total.CNV))
  }
  training.idx <- which(bed$control.CNV.freq >= bed$case.CNV.freq
                        & bed$total.CNV >= min.total.CNV)
  cat(paste("gene_burden_test.R::build.null identified",
            prettyNum(length(training.idx), big.mark=","),
            "genes for null distribution fitting with at least", 
            prettyNum(min.total.CNV, big.mark=","),
            "CNVs in total between cases & controls, and where",
            "control frequency is at least equal to case frequency\n"))
  if(use.unweighted.controls==T){
    bed$delta.total <- bed$case.CNV.w - bed$control.CNV
    bed$delta.norm <- 100*(bed$case.CNV.w.norm - bed$control.CNV.freq)
  }else{
    bed$delta.total <- bed$case.CNV.w - bed$control.CNV.w
    bed$delta.norm <- 100*(bed$case.CNV.w.norm - bed$control.CNV.w.norm)
  }
  null.vals.oneside <- abs(bed$delta.norm[training.idx])
  null.vals <- c(-null.vals.oneside, null.vals.oneside)
  # Fit gaussian null
  gaus.fit <- fitdistr(null.vals, "normal")
  null.mean <- 0
  # null.sd <- sd(null.vals)
  null.sd <- as.numeric(gaus.fit$estimate[2])
  gaus.loglik <- gaus.fit$loglik
  # Fit exponential null
  exp.fit <- fitdistr(null.vals.oneside, "exponential")
  null.exp.rate <- as.numeric(exp.fit$estimate)
  exp.loglik <- exp.fit$loglik  
  
  # Plot distributions & fits
  jpeg(null.dist.plot.out, height=300*3, width=300*7, res=300)
  plot.null(null.vals, null.mean, null.sd, null.exp.rate, cohort.name)
  dev.off()
  
  # Return fit values
  data.frame("cohort"=cohort.name,
             "genes_trained"=length(training.idx),
             "min_CNV_per_gene"=min.total.CNV,
             "gaussian_mean"=round(null.mean, precision),
             "gaussian_sd"=round(null.sd, precision),
             "gaussian_loglik"=round(gaus.loglik, precision),
             "exponential_rate"=round(null.exp.rate, precision),
             "exponential_loglik"=round(exp.loglik, precision))
}


# Read null model parameters from prespecified table (for burden testing)
parse.null.table <- function(null.table.in, cohort.name, cnv, null.model){
  nt <- read.table(null.table.in, sep="\t", header=T, comment.char="")
  if(null.model=="gaussian"){
    null.params <- c("mean"=nt$gaussian_mean[which(nt$CNV==cnv & nt$cohort==cohort.name)],
                     "sd"=nt$gaussian_sd[which(nt$CNV==cnv & nt$cohort==cohort.name)])
  }else if(null.model=="exponential"){
    null.params <- c("rate"=nt$exponential_rate[which(nt$CNV==cnv & nt$cohort==cohort.name)])
  }else{
    stop(paste("gene_burden_test.R::parse.null.table does not ",
               "recognize null.model '", null.model, "'", sep=""))
  }
  return(null.params)
}


# Z-test for case:control CNV burden per gene
z.burden.test <- function(bed, use.unweighted.controls, null.model, null.params, 
                          min.total.CNV, precision){
  # Prepare data for Z-test
  bed$total.CNV <- bed$case.CNV + bed$control.CNV
  bed$total.CNV.w <- bed$case.CNV.w + bed$control.CNV.w
  if(use.unweighted.controls==T){
    bed$delta.total <- bed$case.CNV.w - bed$control.CNV
    bed$delta.norm <- 100*(bed$case.CNV.w.norm - bed$control.CNV.freq)
  }else{
    bed$delta.total <- bed$case.CNV.w - bed$control.CNV.w
    bed$delta.norm <- 100*(bed$case.CNV.w.norm - bed$control.CNV.w.norm)
  }
  
  # Compute test statistic & p-value
  if(null.model=="gaussian"){
    bed$Zscore <- bed$delta.norm / null.params[2]
    bed$pvalue <- pnorm(bed$Zscore, lower.tail=F)
  }else if(null.model=="exponential"){
    bed$Zscore <- NA
    bed$pvalue <- pexp(bed$delta.norm, rate=null.params[1], lower.tail=F)
  }
  
  # Overwrite p-values for genes with fewer than min.total.CNV
  bed$pvalue[which(bed$total.CNV < min.total.CNV)] <- NA
  
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
                           "ln_odds_ratio" = round(bed$ln.OR, precision),
                           "ln_odds_ratio_lower" = round(bed$ln.OR.lower, precision),
                           "ln_odds_ratio_upper" = round(bed$ln.OR.upper, precision),
                           "z_score" = round(bed$Zscore, precision),
                           "phred_p" = round(-log10(bed$pvalue), precision))
  
  return(burden.bed)
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
  make_option(c("--null-table-in"), type="character", default=NULL,
              help="table with null distribution parameters [default %default]",
              metavar="file"),
  make_option(c("--null-model"), type="character", default="gaussian",
              help="specify null distribution to model (options: 'gaussian', 'exponential') [default %default]",
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
  make_option(c("--min-cnvs"), type="numeric", default=NULL,
              help="do not include genes with fewer than N total CNVs in model [default: infer cutoff from data]",
              metavar="numeric"),
  make_option(c("--build-null"), action="store_true", default=FALSE,
              help="build null distribution *instead* of performing burden testing [default %default]"),
  make_option(c("--null-dist-plot"), type="character", default="null_dist.jpg",
              help="path to plot null distribution (only used with --build-null) [default %default]",
              metavar="path"),
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
if(is.null(opts$`null-table-in`)
   & opts$`build-null` == FALSE){
  stop("Must provide --null-table-in or specify --build-null")
}
if(!(opts$`null-model` %in% c("gaussian", "exponential"))){
  stop("--null-model must be either 'gaussian' or 'exponential'")
}

# Writes args & opts to vars
bed.in <- args$args[1]
outfile <- args$args[2]
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
min.total.CNV <- opts$`min-cnvs`
use.unweighted.controls <- opts$`unweighted-controls`
do.build.null <- opts$`build-null`
null.dist.plot.out <- opts$`null-dist-plot`
precision <- opts$precision

# # DEV PARAMETERS (FOR NULL BUILDING):
# bed.in <- "~/scratch/meta1.HP0000118.uCNV.DEL.gene_burden.counts.bed.gz"
# outfile <- "~/scratch/test_gene_burden.null.txt"
# pheno.table.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# cohort.name <- "meta1"
# case.hpo <- "HP:0000118"
# control.hpo <- "HEALTHY_CONTROL"
# case.col.name <- "case_cnvs"
# control.col.name <- "control_cnvs"
# weighted.suffix <- "_weighted"
# use.unweighted.controls <- F
# min.total.CNV <- NULL
# do.build.null <- T
# null.dist.plot.out <- "~/scratch/null_dist.jpg"
# precision <- 6

# # DEV PARAMETERS (FOR BURDEN TESTING):
# bed.in <- "~/scratch/meta1.HP0001250.uCNV.DEL.gene_burden.counts.bed.gz"
# outfile <- "~/scratch/test_gene_burden.stats.txt"
# pheno.table.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# null.table.in <- "~/scratch/uCNV.gene_burden.all_null.fits.txt"
# null.model <- "gaussian"
# cohort.name <- "meta1"
# cnv <- "DEL"
# case.hpo <- "HP:0001250"
# control.hpo <- "HEALTHY_CONTROL"
# case.col.name <- "case_cnvs"
# control.col.name <- "control_cnvs"
# weighted.suffix <- "_weighted"
# use.unweighted.controls <- F
# min.total.CNV <- 1
# do.build.null <- F
# precision <- 6

# Extract sample counts
sample.counts <- get.sample.counts(pheno.table.in, cohort.name, 
                                   case.hpo, control.hpo)

# Process input BED
bed <- import.bed(bed.in, case.col.name, sample.counts$case.n,
                  control.col.name, sample.counts$control.n)

# Build null, if optioned
if(do.build.null == T){
  # Build null
  null.dist <- build.null(bed, use.unweighted.controls, min.total.CNV, 
                          precision, null.dist.plot.out, cohort.name)
  write.table(null.dist, outfile, sep="\t",
              col.names=T, row.names=F, quote=F)
  # Otherwise, run burden tests
}else{
  # Reset parameters
  if(is.null(min.total.CNV)){
    min.total.CNV <- 0
  }
  
  # Read null distribution
  null.params <- parse.null.table(null.table.in, cohort.name, 
                                  cnv, null.model)
  
  # Run burden tests
  burden.bed <- z.burden.test(bed, use.unweighted.controls, 
                              null.model, null.params, 
                              min.total.CNV, precision)
  
  # Format output
  colnames(burden.bed)[1] <- "#chr"
  if(length(grep(".gz", outfile, fixed=T)) > 0){
    outfile <- tools::file_path_sans_ext(outfile)
  }
  write.table(burden.bed, outfile, sep="\t", col.names=T, row.names=F, quote=F)
  if(opts$bgzip == TRUE){
    system(paste("bgzip -f", outfile),wait=T)
  }
}
