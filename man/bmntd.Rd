\name{bmntd}
\alias{bmntd}
\title{
beta mean nearest taxon distance (betaMNTD)
}
\description{
Calculates beta MNTD (beta mean nearest taxon distance, Webb et al 2008) for taxa in each pair of communities in a givern community matrix.
}
\usage{
bmntd(comm, pd, abundance.weighted = TRUE,
      exclude.conspecifics = FALSE,time.output=FALSE,
      unit.sum=NULL, spname.check = TRUE, silent = TRUE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are OTU ids.}
  \item{pd}{matrix, pairwise phylogenetic distance matrix.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
  \item{exclude.conspecifics}{logic, whether conspecific taxa in different communities be exclude from beta MNTD calculations, default is FALSE.}
  \item{time.output}{logic, whether to count calculation time, default is FALSE.}
  \item{unit.sum}{NULL or a number or a nemeric vector. When unit.sum is not NULL and a beta diversity index is calculated for a bin, the taxa abundances will be divided by unit.sum to calculate the relative abundances. usually, unit.sum can be set as the sequencing depth in each sample. Default setting is NULL, means not to do this special transformation.} 
  \item{spname.check}{logic, whether to check the species names in comm and pd.} 
  \item{silent}{logic, if FALSE, some messages will be showed if any mismatch in spcies names.}
}
\details{
beta mean nearest taxon distance for taxa in each pair of communities. Modified from 'comdistnt' in package 'picante'(Kembel et al 2010), this function includes matrix multiplication to be efficient for medium size dataset. 
}
\value{
result is a distance object of pairwise beta MNTD between samples. 
}
\references{
Webb, C.O., Ackerly, D.D. & Kembel, S.W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics, 24, 2098-2100.

Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D. et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26, 1463-1464.
}
\author{
Daliang Ning
}
\note{
Version 3: 2020.8.16, add examples.
Version 2: 2018.10.15, add unit.sum option. if unit.sum!=NULL, will calculate relative abundance according to unit.sum.
Version 1: 2015.9.23  
}
\seealso{
\code{\link{bNTIn.p}}
}
\examples{
data("example.data")
comm=example.data$comm
pd=example.data$pd
bmntd.wt=bmntd(comm, pd, abundance.weighted = TRUE,
               exclude.conspecifics = FALSE)
}

\keyword{phylogenetic}
