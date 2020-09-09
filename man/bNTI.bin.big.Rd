\name{bNTI.bin.big}
\alias{bNTI.bin.big}
\title{
Calculate beta nearest taxon index (betaNTI) for each phylogenetic bin
}
\description{
Perform null model test based on a phylogenetic beta diversity index, beta mean phylogenetic distance to the nearest taxon (betaMNTD), in each bin; calculate beta nearest taxon index (betaNTI), or modified Raup-Crick metric, or confidence level based on the comparison between observed and null betaMNTD in each bin. The package bigmemory (Kane et al 2013) is used to handle very large phylogenetic distance matrix.
}
\usage{
bNTI.bin.big(comm, pd.desc, pd.spname, pd.wd, pdid.bin,
             sp.bin, spname.check = FALSE, nworker = 4,
             memo.size.GB = 50, weighted = c(TRUE, FALSE),
             rand = 1000, output.bMNTD = c(FALSE, TRUE),
             sig.index=c("SES","Confidence","RC","bNTI"),
             unit.sum = NULL, correct.special = FALSE, detail.null = FALSE,
             special.method = c("MNTD","MPD","both"),
             ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975,
             exclude.conspecifics = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a species or OTU or gene, thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{pd.desc}{the name of the file to hold the backingfile description of the phylogenetic distance matrix, it is usually "pd.desc" if using default setting in pdist.big function.}
  \item{pd.spname}{character vector, taxa id in the same rank as the big matrix of phylogenetic distances.}
  \item{pd.wd}{folder path, where the bigmemmory file of the phylogenetic distance matrix are saved.}
  \item{pdid.bin}{list, each element is a vector of integer, indicating which rows/columns in the big phylogenetic matrix represent the taxa in a bin.}
  \item{sp.bin}{one-column matrix, rownames are taxa IDs (i.e. OTU IDs), the only column shows the bin ID of each taxon. Bin IDs are integers in the same order as the elements in the list of pdid.bin.}
  \item{spname.check}{logic, whether to check the OTU ids (species names) in community matrix and phylogenetic distance matrix are the same.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memo.size.GB}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{weighted}{Logic, consider abundances or not (just presence/absence). default is TRUE.}
  \item{rand}{integer, randomization times. default is 1000.}
  \item{output.bMNTD}{logic, if TRUE, the output will include betaMNTD.}
  \item{sig.index}{character, the index for null model significance test. SES or bNTI, standard effect size, i.e. beta nearest taxon index (betaNTI); Confidence, percentage of null values less extreme than the observed value, i.e. non-parametric one-side confidence level; RC, modified Raup-Crick index (RC) based on betaMNTD, i.e. count the number of null betaMNTD lower than observed betaMNTD plus a half of the number of null betaMNTD equal to observed betaMNTD, to get alpha, then calculate betaMNTD-based RC as (2 x alpha - 1). default is SES. If input a vector, only the first element will be used.}
   \item{unit.sum}{NULL or a number or a nemeric vector. When a beta diversity index is calculated for a bin, the taxa abundances will be divided by unit.sum to calculate the relative abundances. Usually, unit.sum can be set as the sequencing depth in each sample. Default setting is NULL, means not to do this transformation.}
  \item{correct.special}{logic, whether to correct the special cases. Default is FALSE.}
  \item{detail.null}{logic, if TRUE, the output will include all the null values. Default is FALSE.}
  \item{special.method}{When correct.special is TRUE, which method will be used to check underestimation of deterministic pattern(s) in special cases. MNTD, use null model test of mean distance to the nearest taxon; MPD, use null model test based on mean pairwise distance; both, use null model test of both MPD and MNTD. Default is MNTD.}
  \item{ses.cut}{numeric, the cutoff of significant standard effect size, default is 1.96.}
  \item{rc.cut}{numeric, the cutoff of significant modified Raup-Crick metric, default is 0.95.}
  \item{conf.cut}{numeric, the cutoff of significant one-side confidence level, default is 0.975.}
  \item{exclude.conspecifics}{Logic, should conspecific taxa in different communities be exclude from MNTD calculations? default is FALSE. The same as in the function bmntd.}
}
\details{
The beta nearest taxon index (betaNTI; Webb et al. 2008, Stegen et al. 2012) is calculated for each phylogenetic bin. betaNTI is a standardized measure of the mean phylogenetic distance to the nearest taxon between samples/communities (beta MNTD) and quantifies the extent of terminal clustering, independent of deep level clustering. Parallel computing is used to improve the speed.

The null model algorithm is "taxa shuffle" (Kembel 2009), i.e. shuffling taxa labels across the
tips of the phylogenetic tree to randomize phylogenetic relationships among species. In this function, taxa will be randomized across all bins.

In the betaNTI of each bin, the diagonal are set as zero. If the randomized results are all the same, the standard deviation will be zero and betaNTI will be NAN. In this case, betaNTI will be set as zero, since the observed result is not differentiable from randomized results.

Modified RC (Chase et al 2011) and Confidence (Ning et al 2020) are alternative significance test indexes to evaluate how the observed beta diversity index deviates from null expectation, which could be a better metric than standardized effect size (betaNTI) in some cases, e.g. null values do not follow normal distribution.
}
\value{
Output is a list with following elements:
\item{index}{list, each element is a square matrix of betaNTI values of a bin. The elements (bins) are in the same order as in the input pdid.bin.}
\item{sig.index}{character, indicates the index for null model significance test, SES (i.e. betaNTI), RC, or Confidence.}
\item{betaMNTD.obs}{Output only if output.bMNTD is TRUE. A list, each element is a square matrix of observed beta MNTD values of a bin. The elements (bins) are in the same order as in the input pdid.bin.}
\item{rand}{Output only if detail.null is TRUE. A list, each element is a matrix with null values of beta MNTD for each turnover of a bin. The elements (bins) are in the same order as in the input pdid.bin.}
\item{special.crct}{Output only if detail.null is TRUE. NULL if correct.special is FALSE. A list with three elements, corresponding to three different null model significance testing indexes, i.e. SES, RC, and Confidence. Each element is a matrix, where the value is zero if the result for a turnover of a bin does not need to correct, otherwise there will be a corrected value.}
}
\references{
Webb, C.O., Ackerly, D.D. & Kembel, S.W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics, 24, 2098-2100.

Kembel, S.W. (2009). Disentangling niche and neutral influences on community assembly: assessing the performance of community phylogenetic structure tests. Ecol Lett, 12, 949-960.

Stegen, J.C., Lin, X., Konopka, A.E. & Fredrickson, J.K. (2012). Stochastic and deterministic assembly processes in subsurface microbial communities. Isme Journal, 6, 1653-1664.

Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. & Inouye, B.D. (2011). Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere, 2, 1-11.

Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals the ecological drivers of grassland soil microbial community assembly in response to warming. bioRxiv. doi:10.1101/2020.02.22.960872 (Nature Communications, final revision).

Kane, M.J., Emerson, J., Weston, S. (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning
}
\note{
Version 7: 2020.9.1, remove setwd. change dontrun to donttest and revise save.wd in help doc.
Version 6: 2020.8.18, update help document, add example.
Version 5: 2020.8.1, change RC opiton to sig.index, add detail.null, rc.cut and conf.cut. 
Version 4: 2019.11.6, add exclude.conspecifics.
Version 3: 2018.10.15, add unit.sum, correct.special.
Version 2: 2016.3.26, add RC option.
Version 1: 2015.12.16
}

\seealso{
\code{\link{bNTIn.p}},\code{\link{bmntd}}
}
\examples{
# function 'bNTI.bin.big' is usually used in the main function, 'icamp.big',
# when setting phylo.rand.scale="across",
# means randomization across all bins in phylogenetic null model.

data("example.data")
comm=example.data$comm
tree=example.data$tree
pdid.bin=example.data$pdid.bin
sp.bin=example.data$sp.bin

# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer after change the path for 'save.wd'.
\donttest{
save.wd=tempdir() # please change to the folder you want to save the pd.big output.
nworker=2 # parallel computing thread number
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
rand.time=20 # usually use 1000 for real data.

bNTIbins=bNTI.bin.big(comm=comm, pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label,
                      pd.wd=pd.big$pd.wd, pdid.bin=pdid.bin, sp.bin=sp.bin,
                      spname.check = TRUE, nworker = nworker, memo.size.GB = 50,
                      weighted = TRUE, rand = rand.time, output.bMNTD = FALSE,
                      sig.index="SES", unit.sum = NULL, correct.special = TRUE,
                      detail.null = FALSE, special.method = "MNTD",
                      exclude.conspecifics = FALSE)
}
}

\keyword{phylogenetic}
