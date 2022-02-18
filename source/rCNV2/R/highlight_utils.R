#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Utility functions used for locus highlight plots


#' Parse region for locus highlight
#'
#' Parse a tabix-style coordinate string for locus highlight
#'
#' @param region tabix-style string of region to be plotted
#' @param genome.in path to BEDTools-style .genome file
#' @param export.values boolean indicator to export values to .GlobalEnv \[default: TRUE\]
#' @param return.values boolean indicator to return values as vector \[default: FALSE\]
#'
#' @details if necessary, start and end coordinates will be bounded to the global
#' start and end coordinates of that chromosome as per `genome.in`
#'
#' @return vector (if `return.values` is `TRUE`)
#'
#' @export parse.region.for.highlight
#' @export
parse.region.for.highlight <- function(region, genome.in, export.values=TRUE,
                                       return.values=FALSE){
  genome <- read.table(genome.in, sep="\t", header=F)
  region.parts <- unlist(strsplit(region, split=":"))
  chrom <- region.parts[1]
  chrom.end <- genome[which(genome[, 1] == chrom), 2]
  start <- max(c(1, as.numeric(unlist(strsplit(region.parts[2], split="-"))[1])))
  end <- min(c(as.numeric(unlist(strsplit(region.parts[2], split="-"))[2]), chrom.end))
  if(export.values){
    assign("chrom", chrom, envir=.GlobalEnv)
    assign("start", start, envir=.GlobalEnv)
    assign("end", end, envir=.GlobalEnv)
  }
  if(return.values){
    return(c(chrom, start, end))
  }
}


#' Load genes from GTF for plotting
#'
#' Extract gene features from a GTF for plotting
#'
#' @param gtf.in path to .gtf
#' @param region coordinates of region to be extracted
#' @param rstudio.local boolean to indicate local Rstudio environment \[default: FALSE\]
#'
#' @details `region` coordinates must be a `tabix`-compatbile string
#'
#' @return data.frame
#'
#' @export
load.genes.from.gtf <- function(gtf.in, region, rstudio.local=FALSE){
  # Required for bedr in local Rstudio only:
  if(rstudio.local){
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/anaconda3/envs/py3/bin", sep = ":"))
  }

  # Tabix region of interest
  if(!file.exists(paste(gtf.in, "tbi", sep="."))){
    stop(paste("tabix index not found for input file", gtf.in))
  }
  require(bedr, quietly=T)
  gtf <- bedr::tabix(region, gtf.in, check.chr=FALSE, verbose=FALSE)

  # Reformat entries
  if(!is.null(gtf)){
    colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score",
                       "strand", "frame", "attribute")
    gtf$gene <- sapply(gtf$attribute, function(atrs.str){
      atrs <- unlist(strsplit(atrs.str, split=";"))
      parts <- unlist(strsplit(atrs[grep("gene_name", atrs)], split=" "))
      gsub("\"", "", parts[length(parts)])
    })
    gtf$transcript <- sapply(gtf$attribute, function(atrs.str){
      atrs <- unlist(strsplit(atrs.str, split=";"))
      parts <- unlist(strsplit(atrs[grep("transcript_id", atrs)], split=" "))
      gsub("\"", "", parts[length(parts)])
    })
    gtf <- gtf[, c("chr", "start", "end", "gene", "strand", "feature", "transcript")]
    gtf[, c("start", "end")] <- apply(gtf[, c("start", "end")], 2, as.numeric)
  }else{
    gtf <- data.frame("chr"=character(), "start"=numeric(), "end"=numeric(),
                      "gene"=character(), "strand"=character(), "feature"=character(),
                      "transcript"=character())
  }

  return(gtf)
}


#' Load sumstats from BED
#'
#' Extract meta-analysis summary statistics from a single BED
#'
#' @param bedpath path to .bed file with summary statistics
#' @param region coordinates of region to be extracted
#' @param rstudio.local boolean to indicate local Rstudio environment \[default: FALSE\]
#'
#' @details `region` coordinates must be a `tabix`-compatbile string
#'
#' @return data.frame
#'
#' @export
load.sumstats.for.region <- function(bedpath, region, rstudio.local=FALSE){
  # Required for bedr in local Rstudio only:
  if(rstudio.local){
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/anaconda3/envs/py3/bin", sep = ":"))
  }

  # Tabix region of interest
  if(!file.exists(paste(bedpath, "tbi", sep="."))){
    stop(paste("tabix index not found for input file", bedpath))
  }
  require(bedr, quietly=T)
  ss <- bedr::tabix(region, bedpath, check.chr=FALSE, verbose=FALSE)

  # Add midpoint
  ss$pos <- (ss$start + ss$stop)/2

  # Return columns of interest
  return(ss[, c("chr", "pos", "meta_neg_log10_p", "meta_lnOR",
                "meta_lnOR_lower", "meta_lnOR_upper")])
}


#' Load PIPs for a list of genes
#'
#' Extract PIPs from BED for specified genes for locus highlight plotting
#'
#' @param pips.in path to .tsv of PIPs per gene
#' @param genes vector of gene symbols to retain
#' @param hpos HPO(s) of interest \[default: do not filter on HPO\]
#'
#' @return data.frame
#'
#' @export
load.pips.for.genelist <- function(pips.in, genes, hpos=NULL){
  pips <- read.table(pips.in, header=T, sep="\t", comment.char="", check.names=F)
  pips <- pips[which(pips$gene %in% genes),
               c("gene", "PIP_final", "credible_set", "#HPO")]
  colnames(pips) <- c("gene", "PIP", "credible_set", "HPO")
  if(!is.null(hpos)){
    pips <- pips[which(pips$HPO == hpo), ]
  }
  pips <- pips[, 1:3]
  return(pips[which(!duplicated(pips)), ])
}


#' Load sample sizes for locus highlights
#'
#' Load sample sizes for case/control contrast for locus highlight plotting
#'
#' @param table.in path to .tsv with sample sizes per HPO
#' @param case.hpos vector of case HPOs to evaluate \[default: c("HP:0000118")\]
#' @param ctrl.hpo control HPO \[default: "HEALTHY_CONTROL"\]
#'
#' @return list with the following two elements:
#' 1. `$case` vector of case sample sizes
#' 2. `$ctrl` vector of control sample sizes
#'
#' @export
get.sample.sizes.for.highlight <- function(table.in, case.hpos=c("HP:0000118"),
                                           ctrl.hpo="HEALTHY_CONTROL"){
  n <- read.table(table.in, header=T, sep="\t", comment.char="")
  n.case <- apply(n[which(n[, 1] %in% case.hpos),
                    grep("meta", colnames(n), fixed=T)],
                  2, max, na.rm=T)
  n.ctrl <- n[which(n[, 1] == ctrl.hpo), grep("meta", colnames(n), fixed=T)]
  return(list("case"=as.numeric(as.vector(n.case)),
              "ctrl"=as.numeric(as.vector(n.ctrl))))
}


#' Load quantitative feature
#'
#' Load a quantitative feature track from a BED file for locus highlights
#'
#' @param bedpath path to .bed file with feature information
#' @param region coordinates of region to be extracted
#' @param keep.col column number to use as feature values \[default: 4\]
#' @param rstudio.local boolean to indicate local Rstudio environment \[default: FALSE\]
#'
#' @details `region` coordinates must be a `tabix`-compatbile string
#'
#' @return data.frame
#'
#' @export
load.feature.bed.for.highlight <- function(bedpath, region, keep.col=4,
                                           rstudio.local=FALSE){
  # Required for bedr in local Rstudio only:
  if(rstudio.local){
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/anaconda3/envs/py3/bin", sep = ":"))
  }

  # Tabix region of interest
  if(!file.exists(paste(bedpath, "tbi", sep="."))){
    stop(paste("tabix index not found for input file", bedpath))
  }
  require(bedr, quietly=T)
  bed <- bedr::tabix(region, bedpath, check.chr=FALSE, verbose=FALSE)
  if(!is.null(bed)){
    bed <- bed[, c(1:3, keep.col)]
    colnames(bed) <- c("chr", "start", "end", "value")
    bed$value <- as.numeric(bed$value)
  }else{
    bed <- data.frame("chr"=character(), "start"=numeric(),
                      "end"=numeric(), value=numeric())
  }

  return(bed)
}


#' Load ChromHMM color code
#'
#' Load ChromHMM state color code as HEX from manifest .tsv
#'
#' @param chromhmm.manifest.in path to ChromHMM manifest .tsv
#'
#' @return named vector of colors per ChromHMM state
#'
#' @export
load.chromhmm.colors <- function(chromhmm.manifest.in){
  mfst <- read.table(chromhmm.manifest.in, sep="\t", comment.char="", header=T)
  colors <- sapply(mfst$COLOR.CODE, function(str){
    cvals <- as.numeric(unlist(strsplit(str, split=",")))
    rgb(cvals[1], cvals[2], cvals[3], maxColorValue=255)
  })
  names(colors) <- mfst[, 1]
  return(colors)
}


#' Load ChromHMM tracks for locus highlight
#'
#' Load a set of ChromHMM tracks from an input .tsv and apply color scheme
#'
#' @param chromhmm.tracks.in .tsv of ChromHMM tracks to be loaded
#' @param chromhmm.manifest.in path to ChromHMM manifest .tsv
#' @param region coordinates of region to be extracted
#'
#' @return data.frame
#'
#' @seealso [load.chromhmm.colors]
#'
#' @export
load.chromhmm.tracks <- function(chromhmm.tracks.in, chromhmm.manifest.in, region){
  tlist <- read.table(chromhmm.tracks.in)[, 1]
  chmm.colors <- load.chromhmm.colors(chromhmm.manifest.in)
  lapply(tlist, function(tpath){
    track <- load.feature.bed(tpath, region)
    track$color <- chmm.colors[track$value]
    return(track)
  })
}


#' Load CNVs for a region of interest
#'
#' Load CNVs from a single BED for a single region, and split by case/control
#'
#' @param bedpaths paths to one or more .bed files with CNVs
#' @param region coordinates of region to be extracted
#' @param cnv filter to this CNV type \[default: do not filter by CNV type\]
#' @param case.hpos vector of case HPOs to evaluate \[default: c("HP:0000118")\]
#' @param ctrl.hpo control HPO \[default: "HEALTHY_CONTROL"\]
#' @param rstudio.local boolean to indicate local Rstudio environment \[default: FALSE\]
#'
#' @details `region` coordinates must be a `tabix`-compatbile string. If multiple
#' .bed files are provided as a vector to `bedpaths`, the contents of these files
#' will be concatenated after loading.
#'
#' @return list with the following two elements:
#' 1. `$case` data.frame of case CNVs
#' 2. `$ctrl` data.frame of control CNVs
#'
#' @export
load.cnvs.from.region <- function(bedpaths, region, cnv=NULL,
                                  case.hpos=c("HP:0000118"),
                                  ctrl.hpo="HEALTHY_CONTROL",
                                  rstudio.local=FALSE){
  # Required for bedr in local Rstudio only:
  if(rstudio.local){
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/anaconda3/envs/py3/bin", sep = ":"))
  }

  # Tabix region of interest
  require(bedr, quietly=T)
  cnvs <- as.data.frame(do.call("rbind", lapply(bedpaths, function(bedpath){
    if(!file.exists(paste(bedpath, "tbi", sep="."))){
      stop(paste("tabix index not found for input file", bedpath))
    }
    bedr::tabix(region, bedpath, check.chr=FALSE, verbose=FALSE)
  })))

  if(nrow(cnvs) > 0){
    # Ensure consistent column names
    colnames(cnvs) <- c("chr", "start", "end", "cnv_id", "cnv", "pheno")

    # Sort & filter CNVs
    cnvs <- cnvs[with(cnvs, order(start, end)), ]
    if(!is.null(cnv)){
      cnvs <- cnvs[which(cnvs$cnv==cnv), ]
    }
    case.cnv.idxs <- which(sapply(cnvs$pheno, function(pstr){
      any(case.hpos %in% unlist(strsplit(pstr, split=";", fixed=T)))
    }))
    case.cnvs <- cnvs[case.cnv.idxs, ]
    ctrl.cnvs <- cnvs[grep(ctrl.hpo, cnvs$pheno, fixed=T), ]
  }else{
    empty.df <- data.frame("chr"=character(), "start"=numeric(),
                           "end"=numeric(), "cnv_id"=character(),
                           "cnv"=character(), "pheno"=character())
    case.cnvs <- empty.df
    ctrl.cnvs <- empty.df
  }

  return(list("case"=case.cnvs, "ctrl"=ctrl.cnvs))
}


#' Load multiple CNV BEDs
#'
#' Load all CNVs for a specific region from a list of BEDs
#'
#' @param cnvlist .tsv of CNV .bed files to load
#' @param region coordinates of region to be extracted
#' @param cnv filter to this CNV type \[default: do not filter by CNV type\]
#' @param case.hpos vector of case HPOs to evaluate \[default: c("HP:0000118")\]
#' @param ctrl.hpo control HPO \[default: "HEALTHY_CONTROL"\]
#' @param rstudio.local boolean to indicate local Rstudio environment \[default: FALSE\]
#'
#' @details `region` coordinates must be a `tabix`-compatbile string
#'
#' @return list of outputs from [load.cnvs.from.region()] for each row in `cnvlist`
#'
#' @seealso [load.cnvs.from.region]
#'
#' @export
load.cnvs.from.region.multi <- function(cnvlist, region, cnv=NULL,
                                        case.hpo="HP:0000118",
                                        ctrl.hpo="HEALTHY_CONTROL",
                                        rstudio.local=FALSE){
  # cnvlist should be a tsv of (cohort, path) pairs
  cnv.list <- read.table(cnvlist, header=F, sep="\t")
  cnvs <- lapply(1:nrow(cnv.list), function(i){
    load.cnvs.from.region(cnv.list[i, 2], region, cnv, case.hpo, ctrl.hpo, rstudio.local)
  })
  names(cnvs) <- cnv.list[, 1]
  return(cnvs)
}


#' Format CNVs for locus highlight
#'
#' Transform CNV coordinates to plotting values (with colors) for locus highlights
#'
#' @param cnvs data.frame of CNVs to be plotted
#' @param start left-most plotting coordinate
#' @param end right-most plotting coordinate
#' @param dx plotting resolution along X axis, specified as total number of steps \[default: 100\]
#' @param cnv.height relative height for each CNV \[default: 1\]
#' @param cnv.buffer relative buffer between adjacent CNVs \[default: 0\]
#' @param bevel.switch.pct fraction of start/end of each CNV to be beveled \[default: 0.025\]
#' @param col CNV color \[default: blueblack\]
#' @param highlight.hpo highlight CNVs from this HPO in a different color
#' \[default: plot all CNVs in the same color\]
#' @param highlight.color color to be used for `highlight.hpo`
#' \[default: plot all CNVs in the same color\]
#'
#' @return list ofwith the following two elements:
#' `$cnvs`: plotting values for each individual CNV
#' `$counts`: total number of CNVs overlapping each bin on the X axis
#'
#' @seealso [load.cnvs.from.region]
#'
#' @export
pileup.cnvs.for.highlight <- function(cnvs, start=NULL, end=NULL, dx=100,
                                      cnv.height=1, cnv.buffer=0,
                                      bevel.switch.pct=0.025, col=blueblack,
                                      highlight.hpo=NA, highlight.col=NULL){
  # Set range of values to evaluate
  if(is.null(start)){
    start <- min(cnvs$start)
  }
  if(is.null(end)){
    end <- max(cnvs$end)
  }
  x <- seq(start, end, length.out=dx+1)

  # Set other scaling parameters
  cnv.y.buf <- cnv.height * cnv.buffer

  # Build empty dataframe of CNV pileup
  counts <- data.frame("pos"=x, "count"=0, "idx"=1:length(x))

  # Create plotting values for each CNV
  # Note: must use for loop to increment counts after each CNV is added
  cnv.plot.values <- list()
  if(nrow(cnvs) > 0){
    for(i in 1:nrow(cnvs)){
      # Get CNV info and increment counts
      cnv.id <- cnvs$cnv_id[i]
      cnv.start <- cnvs$start[i]
      cnv.end <- cnvs$end[i]
      cnv.x.idxs <- which(x >= cnv.start & x <= cnv.end)
      counts$count[cnv.x.idxs] <- counts$count[cnv.x.idxs] + cnv.height
      cnv.hpos <- unlist(strsplit(cnvs$pheno[i], split=";", fixed=T))

      # Create plotting vectors for each CNV
      cnv.x <- c(x[cnv.x.idxs], rev(x[cnv.x.idxs]))
      cnv.y <- c(counts$count[cnv.x.idxs] - cnv.y.buf,
                 rev(counts$count[cnv.x.idxs] - cnv.height + cnv.y.buf))

      # Bevel edges by single dx
      # Bevel left edge according to the CNV that came before
      max.change.dist <- ceiling(bevel.switch.pct * dx)
      if(length(cnv.x.idxs) > 0){
        if(i == 1 | cnv.x.idxs[1] == 1){
          dist.to.nearest.change <- 10e10
        }else{
          same.height.idxs <- which(prev.counts$count == max(counts$count[cnv.x.idxs]))
          if(length(same.height.idxs) > 0){
            dist.to.nearest.change <- min(cnv.x.idxs) - max(same.height.idxs)
          }else{
            dist.to.nearest.change <- 10e10
          }
        }
        if(dist.to.nearest.change > max.change.dist){
          cnv.y[1] <- counts$count[cnv.x.idxs[1]] - cnv.height + cnv.y.buf
        }else{
          cnv.y[length(cnv.y)] <- counts$count[cnv.x.idxs[1]] - cnv.y.buf
        }
      }

      # Always bevel right edge sloping outward
      cnv.y[length(cnv.x.idxs)] <- counts$count[cnv.x.idxs[length(cnv.x.idxs)]] - cnv.height + cnv.y.buf

      # Assign color
      if(!is.na(highlight.hpo)){
        if(highlight.hpo %in% cnv.hpos){
          cnv.color <- highlight.col
        }else{
          cnv.color <- col
        }
      }else{
        cnv.color <- col
      }

      # Add CNV plotting values to output list
      cnv.plot.values[[cnv.id]] <- list("x"=cnv.x, "y"=cnv.y, "color"=cnv.color)

      # Save previous CNV's x indexes and counts for comparisons
      prev.x.idxs <- cnv.x.idxs
      prev.counts <- counts
    }
  }

  # Increment final counts by cnv.y.buf (for plotting)
  counts$counts  <- counts$count + cnv.y.buf

  return(list("cnvs"=cnv.plot.values, "counts"=counts))
}

