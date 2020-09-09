\name{midpoint.root.big}
\alias{midpoint.root.big}
\title{
Midpoint root a large phylogeny
}
\description{
This is modified from the function "modpoint.root" in package "phytools". To deal with a large tree, phylogenetic distance is calculated and saved by using bigmemory in advance.
}
\usage{
midpoint.root.big(tree, pd.desc, pd.spname, pd.wd, nworker = 4)
}
\arguments{
  \item{tree}{phylogenetic tree, an object of class "phylo".}
  \item{pd.desc}{the name of the file to hold the backingfile description of the phylogenetic distance matrix, it is usually "pd.desc" if using default setting in pdist.big function.}
  \item{pd.spname}{character vector, taxa id in the same rank as the big matrix of phylogenetic distances.}
  \item{pd.wd}{folder path, where the bigmemmory file of the phylogenetic distance matrix are saved.}
  \item{nworker}{integer, for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
}
\details{
iCAMP analysis need a rooted tree. If it is difficult to figure out the root, midpoint root is recommended for iCAMP analysis. Modified from the function 'midpoint.root' in package 'phytool'(Revell 2012), this function uses bigmemory (Kane et al 2013) to deal with large datasets.
}
\value{
Output is a list with two elements.
\item{tree}{The rooted tree.}
\item{max.pd}{The maximum pairwise phylogenetic distance.}
}
\references{
Farris, J. (1972) Estimating phylogenetic trees from distance matrices. American Naturalist, 106, 645-667.

Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289-290.

Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223. 

Michael J. Kane, John Emerson, Stephen Weston (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning
}
\note{
Version 3: 2020.9.1, remove setwd; change dontrun to donttest and revise save.wd in help doc.
Version 2: 2020.8.19, update help document, add example.
Version 1: 2015.12.16
}
\seealso{
\code{\link{maxbigm}}
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
  
  mroot=midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                          pd.spname = pd.big$tip.label,
                          pd.wd = pd.big$pd.wd, nworker = nworker)
}
}

\keyword{phylogenetic}% use one of  RShowDoc("KEYWORDS")
