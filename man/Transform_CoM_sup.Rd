% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Transform_CoM_sup.R
\name{Transform_CoM_sup}
\alias{Transform_CoM_sup}
\title{Sup function for matrix transformation}
\usage{
Transform_CoM_sup(RandCoM_directory, Ngenes, NTransformedMatrix = 1)
}
\arguments{
\item{RandCoM_directory}{The directory where simulated Co-regulation matrix saved}

\item{Ngenes}{Number of genes in the binary DEG matrix}

\item{NTransformedMatrix}{Number of GenePair vs. SimTimes matrix required}
}
\value{
A summary data frame helps to decide whether the transformation will be good
}
\description{
This function generates a summary of transforming gene-gene co-regulation matrix to GenePairs-Simulations matrix
}
\examples{
trans_summary <- Transform_CoM_sup(RandCoM_directory = "~/RandBM/", Ngenes = 10, NTransformedMatrix = 1)
}
