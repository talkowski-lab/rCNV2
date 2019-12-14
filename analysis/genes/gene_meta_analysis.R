#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform meta-analysis of gene-based CNV associations


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############

# Import all necessary cohort information
load.cohort.info <- function(infile, pheno.table.in, case.hpo, control.hpo){
  cohort.info <- read.table(infile, header=F, sep="\t")
  colnames(cohort.info) <- c("cohort", "path")
  pt <- read.table(pheno.table.in, header=T, sep="\t", comment.char="")  
  sample.counts <- t(sapply(cohort.info$cohort, function(cohort){
    as.numeric(c(pt[which(pt[, 1]==case.hpo),
       which(colnames(pt)==cohort)],
      pt[which(pt[, 1]==control.hpo),
         which(colnames(pt)==cohort)]))
  }))
  colnames(sample.counts) <- c("case.n", "control.n")
  cbind(cohort.info, sample.counts)
}


# Extract sample counts from table
get.sample.counts <- function(pheno.table.in, cohort.info, 
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


# Read an input file of association statistics
read.stats <- function(stats.in, prefix, case.n, control.n){
  # Read data & subset to necessary columns
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")
  colnames(stats)[1] <- "chr"
  cols.to.keep <- c("chr", "start", "end", "gene", 
                    "case_alt", "control_alt", 
                    "fisher_phred_p", "fisher_OR")
  stats <- stats[, which(colnames(stats) %in% cols.to.keep)]
  stats$case_ref <- case.n - stats$case_alt
  stats$control_ref <- control.n - stats$control_alt
  colnames(stats)[-(1:4)] <- paste(prefix, colnames(stats)[-(1:4)], sep=".")
  return(stats)
}


# Single density-colored scatterplot
dens.scatter <- function(x, y, pt.cex=1, parmar=rep(0.5, 4), add.cor=T){
  # Based on heatscatter.R from Colby Chiang
  # (https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R)
  par(mar=parmar, bty="o")
  plot.df <- data.frame("x"=x, "y"=y)
  plot.df <- plot.df[which(!is.infinite(plot.df$x) & !is.infinite(plot.df$y)
                           & !is.na(plot.df$x) & !is.na(plot.df$y)), ]
  dens <- densCols(plot.df$x, plot.df$y, colramp=colorRampPalette(c("black", "white")))
  plot.df$dens <- col2rgb(dens)[1, ] + 1L
  palette <- colorRampPalette(c("#440154", "#33638d", "#218f8d",
                                "#56c667", "#fde725"))(256)
  plot.df$col <- palette[plot.df$dens]
  plot.df <- plot.df[order(plot.df$dens), ]
  plot(plot.df$x, plot.df$y, type="n",
       xlab="", xaxt="n", ylab="", yaxt="n")
  abline(lm(y ~ x, data=plot.df), lwd=2)
  points(plot.df$x, plot.df$y, col=plot.df$col, pch=19, cex=pt.cex)
  if(add.cor==T){
    r <- cor(plot.df$x, plot.df$y)
    text(x=par("usr")[2], y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])),
         labels=bquote(italic(R) == .(formatC(round(r, 3), digits=3))), pos=2)
  }
}


# Scatterplot grid of log odds ratios between cohorts
or.corplot.grid <- function(stats.list, pt.cex=1){
  ncohorts <- length(stats.list)
  cohorts <- names(stats.list)
  par(mfrow=c(ncohorts, ncohorts))
  ymar <- 3
  xmar <- 3.7
  sapply(1:ncohorts, function(r){
    sapply(1:ncohorts, function(c){
      parmar <- c(ymar/2, xmar/2, ymar/2, xmar/2)
      # Set margins
      if(c==1){
        parmar[2] <- xmar-0.2; parmar[4] <- 0.2
      }
      if(c==ncohorts){
        parmar[2] <- 0.2; parmar[4] <- xmar-0.2
      }
      if(r==1){
        parmar[1] <- 0.2; parmar[3] <- ymar-0.2
      }
      if(r==ncohorts){
        parmar[1] <- ymar-0.2; parmar[3] <- 0.2
      }
      # Don't plot self-self correlations
      if(c==r){
        par(mar=parmar)
        plot(x=c(0, 1), y=c(0, 1), type="n",
             xaxt="n", xlab="", yaxt="n", ylab="")
        rect(xleft=par("usr")[1], xright=par("usr")[2],
             ybottom=par("usr")[3], ytop=par("usr")[4],
             bty="n", col="gray95")
        text(x=0.5, y=0.5, labels=bquote(italic(R)==1))
        # Otherwise, plot as normal
      }else{
        dens.scatter(x=log10(stats.list[[c]][, grep("fisher_OR", colnames(stats.list[[c]]), fixed=T)]), 
                     y=log10(stats.list[[r]][, grep("fisher_OR", colnames(stats.list[[r]]), fixed=T)]),
                     parmar=parmar, pt.cex=pt.cex)
      }
      # Add headers & axes
      if(c==1){
        mtext(2, line=0.2, text=cohorts[r], font=2)
      }
      if(r==1){
        mtext(3, line=0.2, text=cohorts[c], font=2)
      }
      if(c==ncohorts){
        axis(4, labels=NA)
        axis(4, tick=F, line=-0.5, las=2)
        mtext(4, line=2.25, text=bquote(log[10](OR)), cex=0.8)
      }
      if(r==ncohorts){
        axis(1, labels=NA)
        axis(1, tick=F, line=-0.5)
        mtext(1, line=2, text=bquote(log[10](OR)) ,cex=0.8)
      }
    })
  })
}


# Merge a list of association statistics
combine.stats <- function(stats.list){
  merged <- stats.list[[1]]
  for(i in 2:length(stats.list)){
    merged <- merge(merged, stats.list[[i]], 
                    by=c("chr", "start", "end", "gene"),
                    all=F, sort=F)
  }
  merged[, -c(1:4)] <- apply(merged[, -c(1:4)], 2, as.numeric)
  return(merged)
}


# Apply empirical continuity correction to meta-analysis data frame
# Per Sweeting et al., Stat. Med., 2004 (section 3.3)
sweeting.correction <- function(meta.df, cc.sum=0.01){
  # Count number of carriers & non-carriers
  n.alt <- sum(meta.df[, grep("_alt", colnames(meta.df), fixed=T)])
  n.ref <- sum(meta.df[, grep("_ref", colnames(meta.df), fixed=T)])
  # Require at least one CNV to be observed
  if(n.alt>0){
    nt <- n.alt
    R <- n.ref/n.alt
    # Pooled odds ratio estimate of all non-zero studies
    nonzero.studies <- which(apply(meta.df[, grep("_alt", colnames(meta.df), fixed=T)], 1, sum)>0)
    nonzero.case.odds <- sum(meta.df$case_alt[nonzero.studies])/sum(meta.df$case_ref[nonzero.studies])
    nonzero.control.odds <- sum(meta.df$control_alt[nonzero.studies])/sum(meta.df$control_ref[nonzero.studies])
    if(nonzero.control.odds>0){
      ohat <- nonzero.case.odds/nonzero.control.odds
      # Apply standard continuity correction of 0.5 to pooled estimate if no CNVs observed in controls  
    }else{
      nonzero.case.odds <- (sum(meta.df$case_alt[nonzero.studies])+0.5)/(sum(meta.df$case_ref[nonzero.studies])+0.5)
      nonzero.control.odds <- (sum(meta.df$control_alt[nonzero.studies])+0.5)/(sum(meta.df$control_ref[nonzero.studies])+0.5)
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
    meta.df$case_alt <- meta.df$case_alt + cor.case_alt
    meta.df$case_ref <- meta.df$case_ref + cor.case_ref
    meta.df$control_alt <- meta.df$control_alt + cor.control_alt
    meta.df$control_ref <- meta.df$control_ref + cor.control_ref
  }
  return(meta.df)
}


# Make meta-analysis data frame for a single window
make.meta.df <- function(stats.merged, cohorts, row.idx, empirical.continuity=T){
  ncohorts <- length(cohorts)
  meta.df <- data.frame("cohort"=1:ncohorts,
                        "control_ref"=as.numeric(stats.merged[row.idx, setdiff(grep("control_ref", colnames(stats.merged), fixed=T),
                                                                               grep("_weighted", colnames(stats.merged), fixed=T))]),
                        "case_ref"=as.numeric(stats.merged[row.idx, setdiff(grep("case_ref", colnames(stats.merged), fixed=T),
                                                                            grep("_weighted", colnames(stats.merged), fixed=T))]),
                        "control_alt"=as.numeric(stats.merged[row.idx, setdiff(grep("control_cnvs", colnames(stats.merged), fixed=T),
                                                                               grep("_weighted", colnames(stats.merged), fixed=T))]),
                        "case_alt"=as.numeric(stats.merged[row.idx, setdiff(grep("case_cnvs", colnames(stats.merged), fixed=T),
                                                                            grep("_weighted", colnames(stats.merged), fixed=T))]),
                        "cohort_name"=cohorts)
  if(empirical.continuity==T){
    meta.df <- sweeting.correction(meta.df)
  }
  return(meta.df)
}


# Perform meta-analysis for a single window
meta.single <- function(stats.merged, cohorts, row.idx, empirical.continuity=T){
  # If no CNVs are observed, return all NAs
  if(sum(stats.merged[row.idx, grep("_cnvs", colnames(stats.merged), fixed=T)])>0){
    meta.df <- make.meta.df(stats.merged, cohorts, row.idx, empirical.continuity)
    # Meta-analysis
    if(model=="re"){
      meta.res <- rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                          measure="OR", data=meta.df, random = ~ 1 | cohort, slab=cohort_name,
                          add=0, drop00=F, correct=F,
                          digits=5, control=list(maxiter=1000, stepadj=0.5))
      as.numeric(c(meta.res$b[1,1], meta.res$ci.lb, meta.res$ci.ub,
                   meta.res$zval, -log10(meta.res$pval)))
    }else if(model=="mh"){
      meta.res <- rma.mh(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                         measure="OR", data=meta.df, slab=cohort_name,
                         add=0, drop00=F, correct=F)
      as.numeric(c(meta.res$b, meta.res$ci.lb, meta.res$ci.ub,
                   meta.res$zval, -log10(meta.res$pval)))
      
    }
  }else{
    rep(NA, 5)
  }
}


# Wrapper function to perform a meta-analysis on all windows
meta <- function(stats.merged, cohorts, model="re"){
  meta.stats <- t(sapply(1:nrow(stats.merged), function(i){
    meta.single(stats.merged, cohorts, i, model)
  }))
  meta.res <- cbind(stats.merged[, 1:4], meta.stats)
  colnames(meta.res) <- c("chr", "start", "end", "gene",
                          "meta_OR", "meta_OR_lower", "meta_OR_upper",
                          "meta_z", "meta_phred_p")
  return(meta.res)
}


################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)
require(metafor, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--pheno-table"), type="character", default=NULL,
              help="table with counts of samples per HPO term per cohort [default %default]",
              metavar="file"),
  make_option(c("--case-hpo"), type="character", default=NULL,
              help="HPO term to use for case samples [default %default]",
              metavar="string"),
  make_option(c("--control-hpo"), type="character", default='HEALTHY_CONTROL',
              help="HPO term to use for control samples [default %default]",
              metavar="string"),
  make_option(c("--or-corplot"), type="character", default=NULL, 
              help="output .jpg file for pairwise odds ratio correlation plot [default %default]",
              metavar="path"),
  make_option(c("--model"), type="character", default="wz", 
              help="specify meta-analysis model ('wz': weighted Z, 'mh': Mantel-Haenzel) [default %default]",
              metavar="string")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog infile outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options
if(is.null(opts$`pheno-table`)){
  stop("Must provide --pheno-table\n")
}
if(is.null(opts$`case-hpo`)){
  stop("Must specify --case-hpo\n")
}

# Writes args & opts to variable
infile <- args$args[1]
outfile <- args$args[2]
pheno.table.in <- opts$`pheno-table`
case.hpo <- opts$`case-hpo`
control.hpo <- opts$`control-hpo`
corplot.out <- opts$`or-corplot`
model <- opts$model

# # Dev parameters
# infile <- "~/scratch/gene_burden_test/gene.meta_test.input.txt"
# outfile <- "~/scratch/gene_burden_test/gene.meta_test.results.bed"
# pheno.table.in <- "~/scratch/gene_burden_test/HPOs_by_metacohort.table.tsv"
# case.hpo <- "HP:0012638"
# control.hpo <- "HEALTHY_CONTROL"
# corplot.out <- "~/scratch/gene_burden_test/gene_corplot.test.jpg"
# model <- "mh"

# Read list of cohorts to meta-analyze
cohort.info <- load.cohort.info(infile, pheno.table.in, case.hpo, control.hpo)
ncohorts <- nrow(cohort.info)
stats.list <- lapply(1:ncohorts, function(i){read.stats(cohort.info[i, 2], cohort.info[i, 1],
                                                        cohort.info[i, 3], cohort.info[i, 4])})
names(stats.list) <- cohort.info[, 1]

# Plot correlations of odds ratios between cohorts, if optioned
if(!is.null(corplot.out)){
  jpeg(corplot.out, res=300,
       height=300*(3.5+(ncohorts/2)),
       width=300*(4+(ncohorts/2)))
  or.corplot.grid(stats.list, pt.cex=0.25)
  dev.off()
}

# Conduct meta-analysis & write to file
stats.merged <- combine.stats(stats.list)
stats.meta <- meta(stats.merged, cohort.info[, 1], model=model)
colnames(stats.meta)[1] <- "#chr"
write.table(stats.meta, outfile, sep="\t",
            row.names=F, col.names=T, quote=F)
