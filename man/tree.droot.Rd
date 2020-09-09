\name{tree.droot}
\alias{tree.droot}
\title{
Distance from root to tip(s) and node(s) on phylogenetic tree
}
\description{
To calculate the distance from root to tip(s) and node(s) on phylogenetic tree
}
\usage{
tree.droot(tree, range = NA, nworker = 4, output.path = FALSE)
}
\arguments{
  \item{tree}{Phylogenetic tree, an object of class "phylo".}
  \item{range}{NA or a vector of integer, to specify the numbering of the tips/nodes of which the distances to root will be calculated. The numbering corresponds to those in the element "edge" of the tree. Default is NA, means to calculate all tips and nodes.}
  \item{nworker}{integer, for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{output.path}{logic, this function will call the function tree.path, if output.path is TRUE, the result of tree.path will be included in the output. Default is FALSE.}
}
\details{
A tool to get distances to root, used in phylogenetic binning.
}
\value{
If output.path is FALSE, output is a matrix where the first column indicates the numbering of nodes/tips and the second column has the distance to root.
If output.path is TRUE, output is a list with two elements.
\item{droot}{matrix, the first column indicates the numbering of nodes/tips and the second column has the distance to root.}
\item{path}{result of tree.path, list of nodes and edge lengthes from root to each tip and/or node.}
}

\author{
Daliang Ning
}
\note{
Version 2: 2020.8.19, add example.
Version 1: 2015.8.19
}
\seealso{
\code{\link{tree.path}}
}
\examples{
tree=ape::rtree(4)
nworker=2 # parallel computing thread number
droot=tree.droot(tree = tree, nworker = nworker)
}

\keyword{phylogenetic}
