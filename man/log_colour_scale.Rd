% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_single_cell.R
\name{log_colour_scale}
\alias{log_colour_scale}
\title{Calculate log colour scale}
\usage{
log_colour_scale(values = vvv3, scale = 0.05, midpoint = "white",
  interpolate = "spline", asymetric = T, print = F)
}
\arguments{
\item{values}{numeric vector, could be matrix or vector of values you will plot
used to calculate range for asymetrical}

\item{scale}{use this to rescale the range, to change what part of the log scale you're on}

\item{midpoint}{colour of central value (0)}

\item{interpolate}{passed to colorRampPalette}

\item{asymetric}{logical; should 0 be at the center of offset to account for range of values}

\item{print}{logical; should colour scale be printed or should vector of hex values be returned}
}
\description{
Instead of taking log of the values and plotting that,
log the colour scale instead
}
