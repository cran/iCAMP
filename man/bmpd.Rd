\name{bmpd}
\alias{bmpd}
\title{
Beta mean pairwise distance (betaMPD)
}
\description{
Calculates mean pairwise distance separating taxa in each pair of communities in a given community matrix.
}
\usage{
bmpd(comm, pd, abundance.weighted = TRUE, na.zero = TRUE,
     time.output = FALSE, unit.sum = NULL)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are OTU ids.}
  \item{pd}{matrix, pairwise phylogenetic distance matrix.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
  \item{na.zero}{logic. when the sum of a row (a sample) is zero in community data matrix, the relative abundance will be NAN. Sometimes, to avoid some problem in following calculation, this kind of NAN value need be set as zero. Defalt is TRUE.}
  \item{time.output}{logic, whether to count calculation time, default is FALSE.}
  \item{unit.sum}{When a beta diversity index is calculated for a bin, the taxa abundances will be divided by unit.sum to calculate the relative abundances. usually, unit.sum can be set as the sequencing depth in each sample. Default setting is NULL, means not to do this special transformation.}
}

\details{
beta mean pairwise distance.
}
\value{
Output is a distance object of pairwise betaMPD between samples.
}
\references{
Webb CO, Ackerly DD, and Kembel SW. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 18:2098-2100
}
\author{
Daliang Ning
}
\note{
Version 3: 2020.8.16, add examples.
Version 2: 2018.10.3, add unit.sum option. if unit.sum!=NULL, will calculate relative abundance according to unit.sum
Version 1: 2015.8.21.
}
\seealso{
\code{\link{bNRIn.p}}
}
\examples{
data("example.data")
comm=example.data$comm
pd=example.data$pd
bmpd.wt=bmpd(comm, pd, abundance.weighted = TRUE)
}

\keyword{phylogenetic}
