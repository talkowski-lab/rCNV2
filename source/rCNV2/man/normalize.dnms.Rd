% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/large_seg_utils.R
\name{normalize.dnms}
\alias{normalize.dnms}
\title{Normalize DNM counts}
\usage{
normalize.dnms(
  segs,
  dnm.cohorts = c("DDD", "ASC", "ASC_unaffected"),
  synonymous.adjustment = TRUE,
  is.data.table = FALSE
)
}
\arguments{
\item{segs}{segment dataframe (imported with \code{\link[=load.segment.table]{load.segment.table()}})}

\item{dnm.cohorts}{vector of cohorts to evaluate (default: DDD, ASC, ASC_unaffected)}

\item{synonymous.adjustment}{indicator to perform outlier-robust linear regression
for inflation in synonymous counts [default: TRUE]}

\item{is.data.table}{indicator that \code{segs} is a data.table [default: FALSE]}
}
\value{
data.frame
}
\description{
Normalize segment DNM counts vs. expected and (optional) synonymous inflation
}
