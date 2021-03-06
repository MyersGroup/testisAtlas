% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_single_cell.R
\name{plot_cell_scores}
\alias{plot_cell_scores}
\title{Plot scores}
\usage{
plot_cell_scores(component = "V25", point_size = 0.6,
  cell_metadata = cell_data)
}
\arguments{
\item{component}{string ; component name as in cell_metadata object e.g. "V12"}

\item{point_size}{numeric}

\item{cell_metadata}{data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.}
}
\value{
ggplot object
}
\description{
Plot scores
}
\details{

}
