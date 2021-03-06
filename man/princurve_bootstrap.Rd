% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_ordering.R
\name{princurve_bootstrap}
\alias{princurve_bootstrap}
\title{Bootstrap Principal Curve Analysis}
\usage{
princurve_bootstrap(input, df_start = 4, df_end = 10, plot = TRUE)
}
\arguments{
\item{input}{numeric matrix; a matrix of points in arbitrary dimension, e.g. column for PC1/Tsne1 and column for PC2/Tsne2}

\item{df_start}{integer; degrees of freedom to start at}

\item{df_end}{integer; degrees of freedom to finish at}

\item{plot}{logical; if TRUE (default) curve is plotted at end of each df iteration}
}
\value{
list object containing results of each sucessive princurve
}
\description{
Bootstrap Principal Curve Analysis
}
\details{
A principal curve is fitted first with df = df_start, then again for df_start+1 iteratively until reaching df_end
}
