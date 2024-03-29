% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_io.R
\name{load.pval.matrix}
\alias{load.pval.matrix}
\title{Load P-value matrix}
\usage{
load.pval.matrix(
  matrix.in,
  has.coords = T,
  ncols.coords = 3,
  p.is.neg.log10 = T
)
}
\arguments{
\item{matrix.in}{path to matrix.tsv}

\item{has.coords}{boolean indicator if \code{matrix.in} is BED3+ formatted [default: \code{TRUE}]}

\item{ncols.coords}{number of leading columns to treat as coordinates (only
used if \code{has.coords = TRUE}) [default: \code{3}]}

\item{p.is.neg.log10}{boolean indicator if P-values in matrix are \code{-log10}-scaled}
}
\value{
list of four objects:
\enumerate{
\item \verb{$coords} : data frame of coordinates (\code{NULL} if \code{has.coords = FALSE})
\item \verb{$pvals} : data frame of p-values
\item \verb{$expected} : expected P-value distribution under true null
\item \verb{$lambdas} : genomic inflation statistic for each column in \verb{$pvals}
}
}
\description{
Loads a precomputed matrix of P-values
}
\details{
This function was designed to process the output from genome-wide
permutation analyses to parameterize false-discovery rates
}
