\name{bNTI.big}
\alias{bNTI.big}
\title{
Beta nearest taxon index (betaNTI) from big data
}
\description{
To calculate pairwise beta nearest taxon index (betaNTI) by randomizing in the whole species pool or within each group. Package bigmemory (Kane et al 2013) is used to deal with large datasets.
}
\usage{
bNTI.big(comm, meta.group=NULL, pd.desc="pd.desc",
         pd.spname,pd.wd, spname.check=TRUE,
         nworker=4, memo.size.GB=50, weighted=TRUE,
         exclude.consp=FALSE,rand=1000,output.dtail=FALSE,
         RC=FALSE, trace=TRUE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a species or OTU or gene, thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. rownames are sample IDs. first column is metacommunity IDs. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples are under the same metacommunity (the same regional species pool).}
  \item{pd.desc}{the name of the file to hold the backingfile description of the phylogenetic distance matrix, it is usually "pd.desc" if using default setting in pdist.big function.}
  \item{pd.spname}{character vector, taxa id in the same rank as the big matrix of phylogenetic distances.}
  \item{pd.wd}{folder path, where the bigmemmory file of the phylogenetic distance matrix are saved.}
  \item{spname.check}{logic, whether to check the OTU ids (species names) in community matrix and phylogenetic distance matrix are the same.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memo.size.GB}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{weighted}{Logic, consider abundances or not (just presence/absence). default is TRUE.}
  \item{exclude.consp}{Logic, should conspecific taxa in different communities be exclude from MNTD calculations? default is FALSE. The same as in the function bmntd.}
  \item{rand}{integer, randomization times. default is 1000.}
  \item{output.dtail}{logic, if TRUE, the betaNTI, RC value, observed betaMNTD, all null betaMNTD values will all be output, if FALSE, only output betaNTI or RC.}
  \item{RC}{logic, whether to use modified RC merics to evaluate significance of betaMNTD insteal of betaNTI (standardized effect size).}
  \item{trace}{logic, whether to show the progress when the code is running.}
}
\details{
The beta nearest taxon index (betaNTI) is a standardized measure of the mean phylogenetic distance to the nearest taxon between samples/communities (betaMNTD) and quantifies the extent of terminal clustering, independent of deep level clustering. There are a lot of null models for randomization, but this function only use phylogeny shuffle (the same as taxa.labels in ses.mntd).

In the output of betaNTI, the diagonal are set as zero. If the randomized results are all the same, the standard deviation will be zero and betaNTI will be NAN. In this case, beta NTI will be set as zero, since the observed result is not differentiable from randomized results. If the observed betaMNTD has NA values, the corresponding betaNTI will remain NA. Modified RC (Chase 2010) is another metric to evaluate how the observed betaMNTD deviates from null expectation, which could be a better metric than standardized effect size (classic betaNTI) in some cases.
}
\value{
If output.detail=FALSE (default), a matrix of betaNTI values (if RC=FALSE) or RC values (if RC=TRUE) is returned. If output.detail=TRUE, a list is returned.

bNTI: a matrix of pairwise betaNTI values;
RC.bMNTD: a matrix of RC values based on null model test of betaMNTD;
bMNTD: observed betaMNTD values;
bMNTD.rand: a matrix of all null results;
}
\references{
Webb, C.O., Ackerly, D.D. & Kembel, S.W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics, 24, 2098-2100.

Kembel, S.W. (2009). Disentangling niche and neutral influences on community assembly: assessing the performance of community phylogenetic structure tests. Ecol Lett, 12, 949-960.

Stegen, J.C., Lin, X., Konopka, A.E. & Fredrickson, J.K. (2012). Stochastic and deterministic assembly processes in subsurface microbial communities. Isme Journal, 6, 1653-1664.

Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. & Inouye, B.D. (2011). Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere, 2, 1-11.

Kane, M.J., Emerson, J., Weston, S. (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{Daliang Ning}
\note{
Version 2: 2020.12.5, included into iCAMP package to improve the function qpen.
Version 1: 2017.7.12
}
\seealso{\code{\link{bmntd.big}},\link{qpen}}
\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree

# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer after change the path for 'save.wd'.
\donttest{
save.wd=tempdir() # please change to the folder you want to save the pd.big output.
nworker=2 # parallel computing thread number
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
rand.time=20 # usually use 1000 for real data.
bNTI=bNTI.big(comm=comm, pd.desc=pd.big$pd.file,
              pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
              spname.check=TRUE, nworker=nworker, memo.size.GB=50,
              weighted=TRUE, exclude.consp=FALSE,rand=rand.time,
              output.dtail=FALSE, RC=FALSE, trace=TRUE)
}
}
\keyword{phylogenetic}
