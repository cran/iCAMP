\name{dist.bin.3col}
\alias{dist.bin.3col}
\title{
Convert a list of dist (or matrixes) to a matrix
}
\description{
Convert a list of distance matrixes (or square matrixes) with the same sample IDs into a matrix.
}
\usage{
dist.bin.3col(dist.bin, obj.name = NULL)
}
\arguments{
  \item{dist.bin}{a list, each element is a distance matrix or square matrix. all elements have exactly the same sample IDs (rownames and colnames) which are in the same order.}
  \item{obj.name}{a character, as a prefix of the bin names.}
}
\details{
A tool to facilitate format transformation in iCAMP analysis.
}
\value{
output is a matrix. The first two columns are sample IDs, and each of the following columns represent an element in the original list which usually is a bin in iCAMP analysis.
}
\author{
Daliang Ning
}
\note{
Version 2: 2020.8.18, add example
Version 1: 2015.8.30
}

\examples{
# let's see a very simple example
bin.dist=as.matrix(dist(1:10))
rownames(bin.dist)<-colnames(bin.dist)<-paste0("Sample",1:10)
dist.bins=list(bin1=bin.dist,bin2=bin.dist+1,bin3=bin.dist*2)
dis.3c=dist.bin.3col(dist.bins,obj.name="test")
}

\keyword{tool}
