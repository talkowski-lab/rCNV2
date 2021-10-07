#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plotting functions used for overview global/cohort-level analyses


#' Load HPO metadata
#'
#' Load a table of HPO metadata
#'
#' @param meta.in path to metadata .tsv
#' @param counts.in path to .tsv of sample counts
#' @param hpos character vector specifying the desired order of output
#'
#' @return data frame of metadata
#'
#' @export load.hpo.metadata
#' @export
load.hpo.metadata <- function(meta.in, counts.in, hpos){
  # Read and merge metadata and counts
  meta <- read.table(meta.in, sep="\t", comment.char="", header=T)
  colnames(meta)[1] <- "HPO"
  counts <- read.table(counts.in, sep="\t", comment.char="", header=T)
  colnames(counts)[1] <- "HPO"
  meta <- merge(meta, counts[, c(1, grep("meta", colnames(counts), fixed=T))],
                by="HPO", all.x=T, all.y=F, sort=F)

  # Clean columns
  meta$parent_terms <- strsplit(meta$parent_terms, split=";")
  meta$child_terms <- strsplit(meta$child_terms, split=";")
  meta$description[which(meta$HPO=="UNKNOWN")] <- "Other / not otherwise specified"

  # Recompute tier based on number of parents, assuming unambiguous mappings
  meta$HPO_tier <- sapply(meta$parent_terms, function(vals){length(vals[which(!is.na(vals))])}) + 1

  # Correct any situations where a child has no parent <=1 tier away
  for(i in 1:nrow(meta)){
    term.tier <- meta$HPO_tier[i]
    if(term.tier > 1){
      lowest.parent.tier <- max(meta$HPO_tier[which(meta$HPO %in% unlist(meta$parent_terms[i]))])
      if(term.tier - lowest.parent.tier > 1){
        meta$HPO_tier[i] <- lowest.parent.tier + 1
      }
    }
  }

  #Sort according to hpos input vector
  meta <- meta[match(hpos, meta$HPO), ]

  # Insert null rows separating neuro & somatic phenotypes
  neuro.idx <- which(hpos %in% neuro.hpos)
  soma.idx <- which(hpos %in% somatic.hpos)
  meta.wnulls <- as.data.frame(rbind(meta[1, ],
                                     NA,
                                     meta[neuro.idx, ],
                                     NA,
                                     meta[soma.idx, ],
                                     NA,
                                     meta[(max(soma.idx) + 1):nrow(meta), ]))
  colnames(meta.wnulls) <- colnames(meta)

  return(meta.wnulls)
}


#' Plot HPO dendrogram
#'
#' Plot a vertical dendrogram of HPO hierarchy
#'
#' @param meta HPO metadata \(see [load.hpo.metadata()]\)
#' @param lwd line width
#' @param color line color
#' @param y.top top Y coordinate
#'
#' @export plot.dendro
#' @export
plot.dendro <- function(meta, lwd=1, color="black", y.top=-0.6, box.wex=1){
  # Get basic plotting info
  ntiers <- max(meta$HPO_tier, na.rm=T)
  nhpos <- nrow(meta)
  y.at <- (1:nhpos)-0.5

  # Prep plot area
  plot(x=NA, y=NA, xlim=c(0, ntiers+box.wex), ylim=c(nhpos, y.top),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")

  # Draw horizontal segments
  segments(x0=meta$HPO_tier-1, x1=ntiers,
           y0=y.at, y1=y.at,
           lwd=lwd, col=color, lend="round")

  # Draw vertical segments & colored boxes for leaves
  sapply(1:nhpos, function(i){
    direct.children <- which(sapply(meta$parent_terms,
                                    function(parents){meta$HPO[i] %in% parents})
                             & meta$HPO_tier == meta$HPO_tier[i] + 1)
    if(length(direct.children) > 0){
      segments(x0=meta$HPO_tier[i], x1=meta$HPO_tier[i],
               y0=y.at[i], y1=y.at[max(direct.children)],
               lwd=lwd, col=color)
    }
    # Above diagonal colored by neuro
    # Below diagonal colored by severity
    hpo.color.top <- get.hpo.color(meta$HPO[i], color.by="neuro")
    hpo.color.bottom <- get.hpo.color(meta$HPO[i], color.by="severity")
    if(!is.na(meta$HPO_tier[i])){
      polygon(x=ntiers+c(0, 0, box.wex, box.wex),
              y=y.at[i]+c(-0.4, 0.4, -0.4, -0.4),
              xpd=T, border=color, col=hpo.color.top)
      polygon(x=ntiers+c(0, box.wex, box.wex, 0),
              y=y.at[i]+c(0.4, 0.4, -0.4, 0.4),
              xpd=T, border=color, col=hpo.color.bottom)
    }
  })
}

