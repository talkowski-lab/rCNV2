% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/large_seg_plot.R
\name{plot.viohist}
\alias{plot.viohist}
\title{Add a violin-histogram hybrid to a plot}
\usage{
\method{plot}{viohist}(
  perm.dat.vals,
  bins,
  y.at,
  width = 0.8,
  obs.val = NA,
  obs.color = NA,
  obs.border = NA,
  obs.pch = 23,
  color = bluewhite,
  border = blueblack,
  diamond.cex = 4,
  y.title = NULL,
  y.title.color = blueblack,
  left.ax.line = F
)
}
\arguments{
\item{perm.dat.vals}{numeric vector of values}

\item{bins}{bins for histogram}

\item{y.at}{Y position for center of violin}

\item{width}{width of violin}

\item{obs.val}{observed value}

\item{obs.color}{color for observed value}

\item{obs.border}{border for observed value}

\item{obs.pch}{shape for observed value}

\item{color}{color for permuted distribution}

\item{border}{border for permuted distribution}

\item{diamond.cex}{expansion factor for observed value}

\item{y.title}{title for Y axis}

\item{y.title.color}{color for Y-axis title [default: blueblack]}

\item{left.ax.line}{add line to left axis [default: FALSE]}
}
\description{
Helper function to add a single mirrored violin-histogram hybrid of values to an existing plot
}
\seealso{
\code{\link[=plot.seg.perms]{plot.seg.perms()}}
}
