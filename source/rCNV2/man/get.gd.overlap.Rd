% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/large_seg_utils.R
\name{get.gd.overlap}
\alias{get.gd.overlap}
\title{Check for genomic disorder overlap}
\usage{
get.gd.overlap(chrom, start, end, segs)
}
\arguments{
\item{chrom}{chromosome}

\item{start}{start coordinate}

\item{end}{end coordinate}

\item{segs}{segment dataframe (imported with \code{\link[=load.segment.table]{load.segment.table()}})}
}
\value{
vector of GD CNV types
}
\description{
Lookup overlap with any known GDs for a set of coordinates
}
