\name{bNTI.cm}
\alias{bNTI.cm}
\title{
Calculate beta nearest taxon index (betaNTI) with parallel computing under multiple metacommunities
}
\description{
Perform null model test based on a phylogenetic beta diversity index, beta mean phylogenetic distance to the nearest taxon (betaMNTD); calculate beta nearest taxon index (betaNTI), or modified Raup-Crick metric, or confidence level, based on the comparison between observed and null betaMNTD. Run by parallel computing. This function can deal with local communities under different metacommunities (regional pools).
}
\usage{
bNTI.cm(comm, dis, nworker = 4, memo.size.GB = 50,
        meta.group = NULL, meta.spool = NULL,
        meta.frequency = NULL, meta.ab = NULL,
        weighted = c(TRUE, FALSE), exclude.consp = FALSE,
        rand = 1000, output.bMNTD = c(FALSE, TRUE),
        sig.index = c("SES", "Confidence", "RC", "bNTI"),
        unit.sum = NULL, correct.special = FALSE,
        detail.null = FALSE, special.method = c("MNTD", "MPD", "both"),
        ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975,
        dirichlet = FALSE)
}
\arguments{
  \item{comm}{community data matrix. rownames are sample names. colnames are species names}
  \item{dis}{Phylogenetic distance matrix.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memo.size.GB}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. rownames are sample IDs. first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.spool}{a list object, each element is a character vector listing all taxa IDs in a metacommunity. The names of the elements indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default is NULL, means to use the observed taxa in comm across samples within the same metacommunity that is defined by meta.group.}
  \item{meta.frequency}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the occurrence frequency of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.frequency as occurrence frequency of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{meta.ab}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the aubndance (or relative abundance) of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.ab as average relative abundance of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{weighted}{Logic, consider abundances or not (just presence/absence). default is TRUE.}
  \item{exclude.consp}{Logic, should conspecific taxa in different communities be exclude from MNTD calculations? default is FALSE. The same as in the function bmntd.}
  \item{rand}{integer, randomization times. default is 1000.}
  \item{output.bMNTD}{logic, if TRUE, the output will include the observed betaMNTD.}
  \item{sig.index}{character, the index for null model significance test. SES or bNTI, standard effect size, i.e. beta nearest taxon index (betaNTI); Confidence, percentage of null values less extreme than the observed value, i.e. non-parametric one-side confidence level; RC, modified Raup-Crick index (RC) based on betaMNTD, i.e. count the number of null betaMNTD lower than observed betaMNTD plus a half of the number of null betaMNTD equal to observed betaMNTD, to get alpha, then calculate betaMNTD-based RC as (2 x alpha - 1). default is SES. If input a vector, only the first element will be used.}
  \item{unit.sum}{NULL or a number or a nemeric vector. When a beta diversity index is calculated for a bin, the taxa abundances will be divided by unit.sum to calculate the relative abundances. Usually, unit.sum can be set as the sequencing depth in each sample. Default setting is NULL, means not to do this transformation.}  
  \item{correct.special}{logic, whether to correct the special cases when calculating bNTI. Default is FALSE.}
  \item{detail.null}{logic, if TRUE, the output will include all the null values. Default is FALSE.}
  \item{special.method}{When correct.special is TRUE, which method will be used to check underestimation of deterministic pattern(s) in special cases. MNTD, use null model test of mean distance to the nearest taxon; MPD, use null model test based on mean pairwise distance; both, use null model test of both MPD and MNTD. Default is MNTD.}
  \item{ses.cut}{numeric, the cutoff of significant standard effect size, default is 1.96.}
  \item{rc.cut}{numeric, the cutoff of significant modified Raup-Crick metric, default is 0.95.}
  \item{conf.cut}{numeric, the cutoff of significant one-side confidence level, default is 0.975.}
  \item{dirichlet}{Logic. If TRUE, the taxonomic null model for correcting special cases will use Dirichlet distribution to generate relative abundances in randomized community matrix. default is FALSE.}
}
\details{
This function is particularly designed for samples from different metacommunities. The null model "taxa shuffle" will be done under different metacommunities, separately (and independently). All other details are the same as the function bNTIn.p. 
}
\value{
Output is a list with following elements:
\item{index}{a square matrix of betaNTI (or RC or Confidence based on betaMNTD) values.}
\item{sig.index}{character, indicates the index for null model significance test, SES (i.e. betaNTI), RC, or Confidence.}
\item{betaMNTD.obs}{Output only if output.bMNTD is TRUE. A square matrix of observed beta MNTD values.}
\item{rand}{Output only if detail.null is TRUE. A matrix with null values of beta MNTD for each turnover.}
\item{special.crct}{Output only if detail.null is TRUE. it will be NULL if correct.special is FALSE. Otherwise, it will be a list with three elements, corresponding to three different null model significance testing indexes, i.e. SES, RC, and Confidence. Each element is a square matrix, where the value is zero if the result for a turnover does not need to correct, otherwise there will be a corrected value.}
}
\references{
Webb, C.O., Ackerly, D.D. & Kembel, S.W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics, 24, 2098-2100.

Kembel, S.W. (2009). Disentangling niche and neutral influences on community assembly: assessing the performance of community phylogenetic structure tests. Ecol Lett, 12, 949-960.

Stegen, J.C., Lin, X., Konopka, A.E. & Fredrickson, J.K. (2012). Stochastic and deterministic assembly processes in subsurface microbial communities. Isme Journal, 6, 1653-1664.

Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. & Inouye, B.D. (2011). Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere, 2, 1-11.

Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.
}
\author{
Daliang Ning
}
\note{
Version 1: 2021.8.4
}
\seealso{
\code{\link{bNTIn.p}}
}
\examples{
data("example.data")
comm=example.data$comm
pd=example.data$pd

# in this example, 10 samples from one metacommunity,
# the other 10 samples from another metacommunity.
meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
rownames(meta.group)=rownames(comm)

nworker=2 # parallel computing thread number
rand.time=4 # usually use 1000 for real data.
bNTI=bNTI.cm(comm=comm, meta.group=meta.group,
             dis=pd, nworker = nworker, memo.size.GB = 50,
             weighted = TRUE, exclude.consp = FALSE, rand = rand.time,
             output.bMNTD = FALSE, sig.index = "SES", unit.sum = NULL,
             correct.special = TRUE, detail.null = FALSE,
             special.method = "MNTD")
}

\keyword{phylogenetic}
