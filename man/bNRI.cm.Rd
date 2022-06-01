\name{bNRI.cm}
\alias{bNRI.cm}
\title{
Calculate beta net relatedness index with parallel computing under multiple metacommunities
}
\description{
Perform null model test based on a phylogenetic beta diversity index, beta mean pairwise distance (betaMPD); calculate beta net relatedness index (betaNRI), or modified Raup-Crick metric, or confidence level based on the comparison between observed and null betaMPD. Run by parallel computing. This function can deal with local communities under different metacommunities (regional pools).
}
\usage{
bNRI.cm(comm, dis, nworker = 4, memo.size.GB = 50,
        meta.group = NULL, meta.spool = NULL,
        meta.frequency = NULL, meta.ab = NULL,
        weighted = c(TRUE, FALSE), rand = 1000,
        output.bMPD = c(FALSE, TRUE),
        sig.index = c("SES", "Confidence", "RC", "bNRI"),
        unit.sum = NULL, correct.special = FALSE,
        detail.null = FALSE, special.method = c("MPD", "MNTD", "both"),
        ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975,
        dirichlet = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a species or OTU or gene, thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{dis}{matrix, pairwise phylogenetic distance matrix.}
  \item{nworker}{integer, for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memo.size.GB}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. rownames are sample IDs. first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.spool}{a list object, each element is a character vector listing all taxa IDs in a metacommunity. The names of the elements indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default is NULL, means to use the observed taxa in comm across samples within the same metacommunity that is defined by meta.group.}
  \item{meta.frequency}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the occurrence frequency of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.frequency as occurrence frequency of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{meta.ab}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the aubndance (or relative abundance) of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.ab as average relative abundance of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{weighted}{logic, whether to use abundance-weighted or unweighted metrics. Default is TRUE.}
  \item{rand}{integer, randomization times. default is 1000.}
  \item{output.bMPD}{logic, if TRUE, the output will include beta mean pairwise distance (betaMPD).}
  \item{sig.index}{character, the index for null model significance test. SES or bNRI, standard effect size, i.e. beta net relatedness index (betaNRI); Confidence, percentage of null values less extreme than the observed value, i.e. non-parametric one-side confidence level; RC, modified Raup-Crick index (RC) based on betaMPD, i.e. count the number of null betaMPD lower than observed betaMPD plus a half of the number of null betaMPD equal to observed betaMPD, to get alpha, then calculate betaMPD-based RC as (2 x alpha - 1). default is SES. If input a vector, only the first element will be used.}
  \item{unit.sum}{NULL or a number or a nemeric vector. When a beta diversity index is calculated for a bin, the taxa abundances will be divided by unit.sum to calculate the relative abundances. Usually, unit.sum can be set as the sequencing depth in each sample. Default setting is NULL, means not to do this transformation.}
  \item{correct.special}{logic, whether to correct the special cases. Default is FALSE.}
  \item{detail.null}{logic, if TRUE, the output will include all the null values. Default is FALSE.}
  \item{special.method}{When correct.special is TRUE, which method will be used to check underestimation of deterministic pattern(s) in special cases. MPD, use null model test based on mean pairwise distance; MNTD, use null model test of mean nearest taxon distance; both, use null model test of both MPD and MNTD. Default is MPD.}
  \item{ses.cut}{numeric, the cutoff of significant standard effect size, default is 1.96.}
  \item{rc.cut}{numeric, the cutoff of significant modified Raup-Crick metric, default is 0.95.}
  \item{conf.cut}{numeric, the cutoff of significant one-side confidence level, default is 0.975.}
  \item{dirichlet}{Logic. If TRUE, the taxonomic null model for correcting special cases will use Dirichlet distribution to generate relative abundances in randomized community matrix. default is FALSE.}
}
\details{
This function is particularly designed for samples from different metacommunities. The null model "taxa shuffle" will be done under different metacommunities, separately (and independently). All other details are the same as the function bNRIn.p. 
}
\value{
Output is a list with following elements:
\item{index}{a square matrix of betaNRI (or RC or Confidence based on betaMPD) values.}
\item{sig.index}{character, indicates the index for null model significance test, SES (i.e. betaNRI), RC, or Confidence.}
\item{betaMPD.obs}{Output only if output.bMPD is TRUE. A square matrix of observed beta MPD values.}
\item{rand}{Output only if detail.null is TRUE. A matrix with null values of beta MPD for each turnover.}
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
\code{\link{bNRIn.p}}
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
bNRI=bNRI.cm(comm=comm, meta.group=meta.group,
             dis=pd, nworker = nworker, memo.size.GB = 50,
             weighted = TRUE, rand = rand.time, output.bMPD = FALSE, 
             sig.index = "SES", unit.sum = NULL, correct.special = TRUE,
             detail.null = FALSE, special.method = "MPD")
}

\keyword{phylogenetic}
