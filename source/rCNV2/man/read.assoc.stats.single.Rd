% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_io.R
\name{read.assoc.stats.single}
\alias{read.assoc.stats.single}
\title{Load single-cohort association statistics}
\usage{
read.assoc.stats.single(stats.in, prefix, p.is.phred, keep.n.cols = 3)
}
\arguments{
\item{stats.in}{path to input BED file with association statistics}

\item{prefix}{cohort name (to be appended to columns)}

\item{p.is.phred}{boolean indicator of the P-value being -log10-scaled in \code{stats.in}}

\item{keep.n.cols}{number of columns from original BED format to retain}
}
\value{
data frame of formatted associtation stats
}
\description{
Load association statistics for a single cohort from an input file
}