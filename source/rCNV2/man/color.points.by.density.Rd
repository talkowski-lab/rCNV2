% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_utils.R
\name{color.points.by.density}
\alias{color.points.by.density}
\title{Color points by density}
\usage{
color.points.by.density(x, y, palette = NULL)
}
\arguments{
\item{x}{independent variable vector}

\item{y}{dependent variable vector}

\item{palette}{256-color palette to be applied based on density [default: \code{viridis()}]}
}
\value{
dataframe of values to be plotted with density and colors
}
\description{
Generate colors for XY scatterplot based on point density
}
\details{
Inspired by heatscatter.R from Colby Chiang:
https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R
}
\seealso{
\code{\link[=viridis]{viridis()}}
}
