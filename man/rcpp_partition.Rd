% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_partition}
\alias{rcpp_partition}
\title{Graph partition}
\usage{
rcpp_partition(neighbors, centers)
}
\arguments{
\item{neighbors}{A matrix defining a the links between the graph nodes.}

\item{centers}{The desired centers to start the algorithm}
}
\value{
a vector indicating the cluster center for each graph node.
}
\description{
Partition of a graph given its centers
}
\details{
Each center in order takes up all its free neighbors. The process is repeated 
until there are no more free regions.
}
\author{
Marcelo Costa
}
