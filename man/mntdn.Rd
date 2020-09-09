\name{mntdn}
\alias{mntdn}
\title{
Mean nearest taxon distance (MNTD)
}
\description{
Calculate mean nearest taxon distance (MNTD) in each community in a given community matrix.
}
\usage{
mntdn(comm, pd, abundance.weighted = TRUE,
      check.name = TRUE, memory.G = 50, time.count = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are OTU ids.}
  \item{pd}{matrix, pairwise phylogenetic distance matrix.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
  \item{check.name}{logic, whether to check the OTU ids (species names) in community matrix and phylogenetic distance matrix are the same.}
  \item{memory.G}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb}
  \item{time.count}{logic, whether to count calculation time, default is FALSE.}
}
\details{
mean nearest taxon distance (MNTD) in each community, using the same algrithm as the function 'mntd' in package 'picante'.
}
\value{
result is a numeric vector with sample names
}
\references{
Webb CO, Ackerly DD, and Kembel SW. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 18:2098-2100

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
\code{\link{NTI.p}}
}
\examples{
data("example.data")
comm=example.data$comm
pd=example.data$pd
mntd=mntdn(comm=comm,pd=pd,abundance.weighted = TRUE)
}
\keyword{phylogenetic}
