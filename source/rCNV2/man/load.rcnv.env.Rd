% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constants.R
\name{load.rcnv.env}
\alias{load.rcnv.env}
\title{Load rCNV2 environment variables}
\usage{
load.rcnv.env(load.strings = TRUE, load.colors = TRUE, load.scales = TRUE)
}
\arguments{
\item{load.strings}{boolean indicator to load character vector constants}

\item{load.colors}{boolean indicator to load color encodings}

\item{load.sclaes}{boolean indicator to load scale-related constants}
}
\description{
Helper function to load all constants and variables used in rCNV2 analyses
}
\details{
\code{load.rcnv.env()} will export constants from an rCNV2-specific
environment into the global environment and can be referenced as any other
variable within your R session.
}
\seealso{
\link{.GlobalEnv}
}
