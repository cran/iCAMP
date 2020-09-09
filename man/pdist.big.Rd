\name{pdist.big}
\alias{pdist.big}
\title{
Pairwise phylogenetic distance matrix from big tree
}
\description{
Calculates between-species phylogenetic distance matrix from a tree, using bigmemory to deal with too large dataset.
}
\usage{
pdist.big(tree, wd = getwd(), tree.asbig = FALSE,
          output = FALSE, nworker = 4, nworker.pd = nworker,
          memory.G = 50, time.count = FALSE,
          treepath.file="path.rda", pd.spname.file="pd.taxon.name.csv",
          pd.backingfile="pd.bin", pd.desc.file="pd.desc",
          tree.backingfile="treeinfo.bin", tree.desc.file="treeinfo.desc")
}
\arguments{
  \item{tree}{phylogenetic tree, an object of class "phylo".}
  \item{wd}{path of a folder to save the big phylogenetic distance matrix, default is current work directory.}
  \item{tree.asbig}{logic, whether to treat tree attributes also as big data, default is FALSE, generally no need to set as TRUE.}
  \item{output}{logic, whether to output the big phylogenetic distance matrix, default is FALSE, generally do not output it, could be too large.}
  \item{nworker}{for parallel computing the tree paths. a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{nworker.pd}{for parallel computing the phylogenetic distance matrix. default is set the same as nworker. may need to set lower than nworker if the matrix is too large.}
  \item{memory.G}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb}
  \item{time.count}{logic, whether to count calculation time, default is FALSE.}
  \item{treepath.file}{character, name of the file saving the tree.path, which is a list of all the nodes and edge lengthes from root to every tip and/or node. it should be a .rda filename.}
  \item{pd.spname.file}{character, name of the file saving the taxa IDs, which has exactly the same order as the row names (and column names) of the big phylogenetic distance matrix. it should be a .csv filename.}
  \item{pd.backingfile}{character, the root name for the file for the cache of the big phylogenetic distance matrix. it should be a .bin filename.}
  \item{pd.desc.file}{character, name of the file to hold the backingfile description for the big phylogenetic distance matrix. it should be a .desc filename.}
  \item{tree.backingfile}{character, the root name for the file for the cache of the 3-column matrix of the tree information, including edge and edge length. it should be a .bin filename.}
  \item{tree.desc.file}{character, name of the file to hold the backingfile description for the tree information matrix. it should be a .desc filename.}
}
\details{
The cophenetic distance between each pair of taxa is calculated (Sokal and Rohlf 1962). Modified from the function "cophenetic" in package "ape" (Paradis & Schliep 2018), this function can calculate pairwise distance from large phylogenetic tree quickly by parallel computing. This function uses bigmemory (Kane et al 2013) to deal with large phylogenetic distance matrix, which will not occupy memory but directly be saved at the hard disk.
}
\value{
Output is a list
  \item{tip.label}{OTU ids or species names, which is tip.label in tree file.}
  \item{pd.wd}{the folder saving the big phylogenetic distance matrix.}
  \item{pd.file}{the folder saving the big phylogenetic distance matrix.}
  \item{pd.name.file}{the file saving the tip.label information.}
}
\references{
Sokal, R. R. & Rohlf, F. J.. (1962). The comparison of dendrograms by objective methods. Taxon, 11:33-40

Paradis, E. & Schliep, K. (2018). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526-528.

Kane, M.J., Emerson, J., & Weston, S. (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning
}
\note{
Version 4: 2020.9.1, remove setwd; add options to specify the file names; change dontrun to donttest and revise save.wd in help doc.
Version 3: 2020.8.19, add example.
Version 2: 2017.3.13
Version 1: 2015.7.24
}
\examples{
data("example.data")
tree=example.data$tree
# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
  save.wd=tempdir() # please change to the folder you want to save the pd.big output.
  nworker=2 # parallel computing thread number
  pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
}
}
\keyword{phylogenetic}
