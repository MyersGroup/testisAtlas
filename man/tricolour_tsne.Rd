% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_single_cell.R
\name{tricolour_tsne}
\alias{tricolour_tsne}
\title{Plot tricolor tSNE}
\usage{
tricolour_tsne(genes, returndf = F, predict = "SDA",
  cell_metadata = cell_data, ptsize = 1.5, jitter = 0.25)
}
\arguments{
\item{genes;}{string vector of two or three gene names}

\item{returndf;}{logical, should the data be returned, or a plot (default)}

\item{predict;}{which factorisation to use for predicted expression}

\item{ptsize;}{numeric, size of the points in the plot}

\item{jitter;}{numeric, how jittered should be points be (to avoid overplotting)}
}
\value{
The Y matrix from Rtsne output, rotated by angle
}
\description{
Plot tricolor tSNE
}
\details{
print tsne with ternary colour scheme for expression
}
