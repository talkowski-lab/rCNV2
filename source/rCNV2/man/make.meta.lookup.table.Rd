% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assoc.R
\name{make.meta.lookup.table}
\alias{make.meta.lookup.table}
\title{Create meta-analysis lookup table}
\usage{
make.meta.lookup.table(
  stats.merged,
  cohorts,
  model,
  adjust.biobanks = F,
  cohort.inflation = NULL,
  probe.counts = NULL,
  empirical.continuity = T
)
}
\arguments{
\item{stats.merged}{association statistics for all cohorts (see \code{\link[=combine.single.cohort.assoc.stats]{combine.single.cohort.assoc.stats()}})}

\item{cohorts}{vector of all cohort names}

\item{model}{specify meta-analysis model to use (see \code{Details})}

\item{adjust.biobanks}{boolean indicator to include biobank status as a categorical
covariate in meta-analysis model [default: FALSE]}

\item{cohort.inflation}{numeric vector of inflation terms for each cohort as
estimated by \code{\link[=estimate.cohort.inflation]{estimate.cohort.inflation()}} [default: NULL]}

\item{probe.counts}{data frame of control probe counts (see \code{Details}) [default: NULL]}

\item{empirical.continuity}{boolean indicator to apply Sweeting empirical
continuity correction [default: TRUE]}
}
\value{
data frame of all unique case & control CNV counts
}
\description{
Make meta-analysis lookup table to shorten time required to run full meta-analysis
}
\details{
\code{model} accepts several character inputs, including:
\itemize{
\item \code{fe} : fixed-effects model implemented by \code{\link[metafor:rma.uni]{metafor::rma.uni()}} [default]
\item \code{re} : random-effects model implemented by \code{\link[metafor:rma.uni]{metafor::rma.uni()}}
\item \code{mh} : Mantel-Haenszel model implemented by \code{\link[metafor:rma.mh]{metafor::rma.mh()}}
}

\code{probe.counts} must be a BED-style data frame matching the exact entries in
stats.merged with one numeric column of probe counts for each cohort.
}
