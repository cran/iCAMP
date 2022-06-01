\name{bmntd.big}
\alias{bmntd.big}
\title{
beta mean nearest taxon distance (betaMNTD) from big data
}
\description{
Calculates beta MNTD (beta mean nearest taxon distance, Webb et al 2008) for taxa in each pair of communities in a givern community matrix, using bigmemory (Kane et al 2013) to deal with too large dataset.
}
\usage{
bmntd.big(comm, pd.desc = "pd.desc", pd.spname, pd.wd,
          spname.check = FALSE, abundance.weighted = TRUE,
          exclude.conspecifics = FALSE, time.output = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are taxa ids.}
  \item{pd.desc}{character, the name to describe bigmemory file of phylogenetic distance matrix, default is "pd.desc".}
  \item{pd.spname}{vector, the OTU ids (species names) in exactly the same order as the phylogenetic matrix rows or columns}
  \item{pd.wd}{the path of the folder saving the phylogenetic distance matrix.}
  \item{spname.check}{logic, whether to check the OTU ids (species names) in community matrix and phylogenetic distance matrix are the same.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
  \item{exclude.conspecifics}{logic, whether conspecific taxa in different communities be exclude from beta MNTD calculations, default is FALSE.}
  \item{time.output}{logic, whether to count calculation time, default is FALSE.}
}
\details{
beta mean nearest taxon distance for taxa in each pair of communities. Improved from 'comdistnt' in package 'picante'(Kembel et al 2010). This function adds bigmemory part (Kane et al 2013) to deal with large dataset.
}
\value{
result is a distance object.
}
\references{
Webb, C.O., Ackerly, D.D. & Kembel, S.W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics, 24, 2098-2100.

Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D. et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26, 1463-1464.

Kane, M.J., Emerson, J., Weston, S. (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning (ningdaliang@ou.edu)
}
\note{
Version 4: 2020.12.5, copy from package NST to iCAMP to improve the function qpen.
Version 3: 2020.9.9, remove setwd; change dontrun to donttest and revise save.wd in help doc.
Version 2: 2020.8.22, add to NST package, update help document.
Version 1: 2017.3.13
}
\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree

# since it needs to save some file to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
wd0=getwd()
save.wd=paste0(tempdir(),"/pdbig.bmntd.big")
# you may change save.wd to the folder you want to save the pd.big output.
nworker=2 # parallel computing thread number
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
bmntd.wt=bmntd.big(comm=comm, pd.desc = pd.big$pd.file,
                   pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                   abundance.weighted = TRUE)
setwd(wd0)
}
}
\keyword{phylogenetic}
