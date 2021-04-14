#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Generate Manhattan and QQ plots for an input set of association statistics


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=TRUE)


# Helper function to read intervals to highlight
read.highlight.bed <- function(highlight.in){
  bed <- read.table(highlight.in, sep="\t", header=F, comment.char="#")[, 1:3]
  colnames(bed) <- c("chr", "start", "end")
  return(bed)
}


# Helper function to generate title panel for composite plots
title.panel <- function(title){
  par(mar=c(0.1, 3.1, 0.1, 1), bty="n")
  plot(0:1, 0:1, xlab="", xaxt="n", ylab="", yaxt="n", type="n")
  subs <- unlist(strsplit(title, split="\n", fixed=T))
  main <- subs[1]
  text(x=0.5, y=1, pos=1, labels=main, font=2)
  if(length(subs) > 1){
    sub.cats <- paste(subs[-1], collapse="\n")
    text(x=0.5, y=0.5, pos=1, labels=sub.cats, cex=0.85, xpd=T)
  }
}


################
###RSCRIPT BLOCK
################
# List of command-line options
option_list <- list(
  make_option(c("--miami"), type="logical", default=F, action="store_true",
              help="Miami plot (instead of Manhattan). Requires [stats2] input. [default %default]"),
  make_option(c("--p-col-name"), type="character", default="P",
              help="column name corresponding to P-values [default %default]",
              metavar="string"),
  make_option(c("--p-is-phred"), type="logical", default=F, action="store_true",
              help="supplied P-values are Phred-scaled (-log10[P]) [default %default]"),
  make_option(c("--max-phred-p"), type="numeric", default=100,
              help="maximum P-value to report; more sigificant values will be rounded down [default %default]"),
  make_option(c("--cutoff"), type="numeric", default=10^-8,
              help="P-value of significance threshold [default %default]",
              metavar="numeric"),
  make_option(c("--highlight-bed"), type="character", default=NA,
              help="BED file of coordinates to highlight [default %default]",
              metavar="BED"),
  make_option(c("--highlight-name"), type="character", default="Highlighted regions",
              help="name for highlighted regions in legend [default %default]",
              metavar="string"),
  make_option(c("--highlight-color"), type="character", default="#4EE69A",
              help="color for highlighted regions [default %default]",
              metavar="string"),
  make_option(c("--cutoff-2"), type="numeric", default=NULL,
              help="P-value of significance threshold in second plot (--miami only) [default %default]",
              metavar="numeric"),
  make_option(c("--highlight-bed-2"), type="character", default=NA,
              help="BED file of coordinates to highlight in second plot (--miami only) [default: uses value passed as --cutoff]",
              metavar="BED"),
  make_option(c("--highlight-name-2"), type="character", default="Highlighted regions",
              help="name for highlighted regions in legend of second plot (--miami only) [default %default]",
              metavar="string"),
  make_option(c("--highlight-color-2"), type="character", default="#4EE69A",
              help="color for highlighted regions of second plot [default %default]",
              metavar="string"),
  make_option(c("--label-prefix"), type="character", default=NULL,
              help="prefix to append to some labels [default %default]",
              metavar="string"),
  make_option(c("--label-prefix-2"), type="character", default=NULL,
              help="prefix to append to labels for second plot (--miami only) [default %default]",
              metavar="string"),
  make_option(c("--title"), type="character", default=NULL,
              help="title for composite Manhattan/Miami & QQ plot [default %default]",
              metavar="string"),
  make_option(c("--echo-lambdas"), type="logical", default=F, action="store_true",
              help="Print QQ lambdas to stdout [default %default]"),
  make_option(c("--qq-only"), type="logical", default=F, action="store_true",
              help="Only plot QQ [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog stats [stats2] out_prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to vars
stats.in <- args$args[1]
miami <- opts$miami
if(miami == F){
  stats2.in <- NULL
  out.prefix <- args$args[2]
}else{
  stats2.in <- args$args[2]
  out.prefix <- args$args[3]
}
p.col.name <- opts$`p-col-name`
p.is.phred <- opts$`p-is-phred`
min.p <- 10^-(opts$`max-phred-p`)
cutoff <- opts$cutoff
highlight.in <- opts$`highlight-bed`
highlight.name <- opts$`highlight-name`
highlight.color <- opts$`highlight-color`
cutoff2 <- opts$`cutoff-2`
if(is.null(cutoff2)){
  cutoff2 <- cutoff
}
highlight2.in <- opts$`highlight-bed-2`
highlight2.name <- opts$`highlight-name-2`
highlight2.color <- opts$`highlight-color-2`
lab.prefix <- opts$`label-prefix`
lab2.prefix <- opts$`label-prefix-2`
title <- gsub("\\n", "\n", opts$title, fixed=T)
echo.lambdas <- opts$`echo-lambdas`
qq.only <- opts$`qq-only`

# Checks for appropriate positional arguments, depending on mode
if(miami == F){
  if(length(args$args) != 2){
    stop("Two positional arguments required for Manhattan (default) mode\n")
  }
}else{
  if(length(args$args) != 3){
    stop("Three positional arguments required for Miami mode\n")
  }
}

# Load data
stats <- load.manhattan.stats(stats.in, p.col.name, p.is.phred, min.p)
if(miami == T){
  stats2 <- load.manhattan.stats(stats2.in, p.col.name, p.is.phred, min.p)
}

if(!is.na(highlight.in)){
  highlights <- read.highlight.bed(highlight.in)
}else{
  highlights <- NULL
}
if(miami == T){
  if(!is.na(highlight2.in)){
    highlights2 <- read.highlight.bed(highlight2.in)
  }else{
    highlights2 <- NULL
  }
}

# Plotting protocol for Manhattan mode
if(miami == F){

  global.p.min <- min(c(stats$p, 10^-(-log10(cutoff) + 1)), na.rm=T)

  # Generate Manhattan plot
  if(qq.only==F){
    manhattan.png.out <- paste(out.prefix, "manhattan.png", sep=".")
    cat(paste("Printing Manhattan plot to", manhattan.png.out, "\n"))
    png(manhattan.png.out,
        height=1000, width=1800, res=400)
    plot.manhattan(stats, cutoff, highlights=highlights,
              highlight.name=highlight.name,
              highlight.color=highlight.color,
              lab.prefix=lab.prefix,
              ymax=-log10(global.p.min))
    dev.off()
  }

  # Generate QQ plot
  qq.png.out <- paste(out.prefix, "qq.png", sep=".")
  cat(paste("Printing QQ plot to", qq.png.out, "\n"))
  png(qq.png.out,
      height=1000, width=1000, res=400)
  plot.qq(stats, cutoff, highlights=highlights,
     highlight.name=highlight.name,
     highlight.color=highlight.color,
     echo.lambdas=echo.lambdas, legend=F,
     ymax=-log10(global.p.min))
  dev.off()

  # Generate combo Manhattan & QQ plot
  if(qq.only==F){
    combo.png.out <- paste(out.prefix, "manhattan_with_qq.png", sep=".")
    cat(paste("Printing combined Manhattan & QQ plots to",
              combo.png.out, "\n"))
    if(is.null(title)){
      png(combo.png.out,
          height=1000, width=2800, res=400)
      layout(matrix(c(1,2), nrow=1), widths=c(18, 10))
    }else{
      png(combo.png.out,
          height=1200, width=2800, res=400)
      layout(matrix(c(1,1,2,3), nrow=2, byrow=T),
             heights=c(2, 10), widths=c(18, 10))
      title.panel(title)
    }
    plot.manhattan(stats, cutoff, highlights=highlights,
              highlight.name=highlight.name,
              highlight.color=highlight.color,
              lab.prefix=lab.prefix,
              ymax=-log10(global.p.min))
    plot.qq(stats, cutoff, highlights=highlights,
       highlight.name=highlight.name,
       highlight.color=highlight.color,
       echo.lambdas=echo.lambdas, legend=F,
       ymax=-log10(global.p.min))
    dev.off()
  }

# Plotting protocol for Miami mode
}else{

  global.p.min <- min(c(stats$p, stats2$p, 10^-(-log10(cutoff) + 1)), na.rm=T)

  # Generate Miami plot
  if(qq.only==F){
    miami.png.out <- paste(out.prefix, "miami.png", sep=".")
    cat(paste("Printing Miami plot to", miami.png.out, "\n"))
    png(miami.png.out,
        height=1600, width=1800, res=400)
    layout(matrix(1:2, nrow=2), heights=c(1, 0.88))
    plot.manhattan(stats, cutoff, highlights=highlights,
              highlight.name=highlight.name,
              highlight.color=highlight.color,
              lab.prefix=lab.prefix,
              ymax=-log10(global.p.min))
    plot.manhattan(stats2, cutoff2, highlights=highlights2,
              highlight.name=highlight2.name,
              highlight.color=highlight2.color,
              lab.prefix=lab2.prefix,
              ymax=-log10(global.p.min),
              reflection=T)
    dev.off()
  }

  # Generate joint Miami plot with QQs
  combo.png.out <- paste(out.prefix, "miami_with_qq.png", sep=".")
  cat(paste("Printing combined Miami & QQ plots to",
            combo.png.out, "\n"))
  if(is.null(title)){
    png(combo.png.out,
        height=1600, width=2600, res=400)
    layout(matrix(1:4, nrow=2, byrow=T),
           heights=c(1, 0.88), widths=c(18, 8))
  }else{
    png(combo.png.out,
        height=0.8*1800, width=0.8*2600, res=400)
    layout(matrix(c(1, 1, 2:5), nrow=3, byrow=T),
           heights=c(0.235, 1, 0.94), widths=c(18, 8))
    title.panel(title)
  }
  label.cex <- 0.75
  plot.manhattan(stats, cutoff, highlights=highlights,
            highlight.name=highlight.name,
            highlight.color=highlight.color,
            lab.prefix=lab.prefix,
            ymax=-log10(global.p.min),
            label.cex=label.cex)
  plot.qq(stats, cutoff, highlights=highlights,
     highlight.name=highlight.name,
     highlight.color=highlight.color,
     echo.lambdas=echo.lambdas, legend=F,
     ymax=-log10(global.p.min),
     label.cex=0.9*label.cex)
  plot.manhattan(stats2, cutoff2, highlights=highlights2,
            highlight.name=highlight2.name,
            highlight.color=highlight2.color,
            lab.prefix=lab2.prefix,
            ymax=-log10(global.p.min),
            reflection=T,
            label.cex=label.cex)
  plot.qq(stats2, cutoff2, highlights=highlights2,
     highlight.name=highlight2.name,
     highlight.color=highlight2.color,
     echo.lambdas=echo.lambdas, legend=F,
     ymax=-log10(global.p.min),
     reflection=T,
     label.cex=0.9*label.cex)
  dev.off()
}

