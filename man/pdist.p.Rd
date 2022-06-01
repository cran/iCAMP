\name{pdist.p}
\alias{pdist.p}
\title{
Pairwise phylogenetic distance matrix from small tree
}
\description{
Calculates between-species phylogenetic distance matrix from a tree. only deal with relatively small dataset.
}
\usage{
pdist.p(tree, nworker = 4, memory.G = 50, silent = FALSE, time.count = FALSE)
}
\arguments{
  \item{tree}{phylogenetic tree, an object of class "phylo".}
  \item{nworker}{for parallel computing the tree paths. a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memory.G}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb}
  \item{silent}{logic, whether to show messages. Default is FALSE, thus all messages will be showed.}
  \item{time.count}{logic, whether to count calculation time, default is FALSE.}
}
\details{
The cophenetic distance between each pair of taxa is calculated (Sokal and Rohlf 1962). Modified from the function "cophenetic" in package "ape" (Paradis & Schliep 2018), this function can calculate pairwise distance from phylogenetic tree quickly by parallel computing. If the tree has too many tips (taxa), please use another function pdist.big designed for large datasets.
}
\value{
Output is a data.frame object, a square matrix of pairwise phylogenetic distances. Row names are the same as column names, indicating taxa IDs. 
}
\references{
Sokal, R. R. & Rohlf, F. J.. (1962). The comparison of dendrograms by objective methods. Taxon, 11:33-40

Paradis, E. & Schliep, K. (2018). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526-528.
}
\author{
Daliang Ning
}
\note{
Version 1: 2021.9.24
}
\seealso{
\code{\link{pdist.big}}
}
\examples{
data("example.data")
tree=example.data$tree
nworker=2 # parallel computing thread number
pd=pdist.p(tree = tree, nworker = nworker)
}
\keyword{phylogenetic}
