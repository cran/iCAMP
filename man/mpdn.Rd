\name{mpdn}
\alias{mpdn}
\title{
Mean pairwise distance (MPD)
}
\description{
Calculate mean pairwise distance (MPD) in each community in a given community matrix.
}
\usage{
mpdn(comm, pd, abundance.weighted = TRUE, time.output = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are OTU ids.}
  \item{pd}{matrix, pairwise phylogenetic distance matrix.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
  \item{time.output}{logic, whether to count calculation time, default is FALSE.}
}
\details{
mean pairwise distance (MPD) in each community, which is the same index as 'mpd' in package 'picante', but calculated by matrix multiplication. 
}
\value{
result is a numeric vector with sample names
}
\references{
Webb C, Ackerly D, McPeek M, and Donoghue M. (2002). Phylogenies and community ecology. Annual Review of Ecology and Systematics 33:475-505.

Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D. et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26, 1463-1464.
}
\author{
Daliang Ning
}
\note{
Version 2: 2020.8.19, update help document, add example.
Version 1: 2017.3.13
}
\seealso{
\code{\link{NRI.p}}
}
\examples{
data("example.data")
comm=example.data$comm
pd=example.data$pd
mpd=mpdn(comm = comm, pd = pd, abundance.weighted = TRUE)
}
\keyword{phylogenetic}
