\name{RC.cm}
\alias{RC.cm}
\title{
Modified Raup-Crick index based on Bray-Curtis similarity under multiple metacommunities
}
\description{
The Raup-Crick based on taxonomic dissimilarity index (i.e. Bray-Curtis) is to use null models to disentangle variation in community dissimilarity from variation in alpha-diversity. This function can deal with local communities under different metacommunities (regional pools).
}
\usage{
RC.cm(comm, rand = 1000, na.zero = TRUE, nworker = 4,
      meta.group = NULL, meta.frequency = NULL, meta.ab = NULL,
      memory.G = 50, weighted = TRUE, unit.sum = NULL,
      sig.index = c("RC", "Confidence", "SES"), detail.null = FALSE,
      output.bray = FALSE, silent = FALSE, taxo.metric = "bray",
      transform.method = NULL, logbase = 2, dirichlet = FALSE)
}

\arguments{
  \item{comm}{Community data matrix. rownames are sample names. colnames are species names.}
  \item{rand}{integer, randomization times, default is 1000.}
  \item{na.zero}{logic. If community data marix has any zero-sum row (sample), Bray-Curtis index will be NA. Somtimes, this kind of NA need be set as zero to avoid some format problem in following calculation. Default is TRUE.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. rownames are sample IDs. first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.frequency}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the occurrence frequency of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.frequency as occurrence frequency of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{meta.ab}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the aubndance (or relative abundance) of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.ab as average relative abundance of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{memory.G}{numeric, to set the memory size as you need, so that calculation of big data will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{weighted}{lOgic, whether to use abundance-weighted metrics. default is TRUE}
  \item{unit.sum}{If unit.sum is set as a number or a numeric vector, the taxa abundances will be divided by unit.sum to calculate the relative abundances, and the Bray-Cuits index in each bin will become manhattan index divided by 2. usually, unit.sum can be set as the sequencing depth in each sample. Default setting is NULL, means not to do this special transformation.}
  \item{sig.index}{character, the index for null model significance test. RC, modified Raup-Crick index (RC) based on Bray-Curtis (BC), i.e. count the number of null BC lower than observed BC plus a half of the number of null BC equal to observed BC, to get alpha, then calculate RCbray as (2 x alpha - 1). SES, standard effect size; Confidence, percentage of null values less extreme than the observed value, i.e. non-parametric one-side confidence level. default is RC. If input a vector, only the first element will be used.}
  \item{detail.null}{logic, if TRUE, the output will include all the null values. Default is FALSE.}
  \item{output.bray}{logic, if TRUE, the output will include observed taxonomic dissimilarity (default is Bray-Curtis).}
  \item{silent}{logic, if FALSE, some messages will show during calculation.}
  \item{taxo.metric}{taxonomic beta diversity index, the same as 'method' in the function 'vegdist' in package 'vegan', including "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis". If taxo.metric='bray' and transform.method=NULL, RC will be calculated based on Bray-Curtis dissimilarity as recommended in original iCAMP; otherwise, unit.sum setting will be ignored.}
  \item{transform.method}{character or a defined function, to specify how to transform community matrix before calculating dissimilarity. if it is a characher, it should be a method name as in the function 'decostand' in package 'vegan', including 'total','max','freq','normalize','range','standardize','pa','chi.square','cmdscale','hellinger','log'.}
  \item{logbase}{numeric, the logarithm base used when transform.method='log'.} 
  \item{dirichlet}{Logic. If TRUE, the taxonomic null model will use Dirichlet distribution to generate relative abundances in randomized community matrix. If the input community matrix has all row sums no more than 1, the function will automatically set dirichlet=TRUE. default is FALSE.}
}
\details{
While all other details are the same as the function RC.pc, this function is particularly designed for samples from different metacommunities. The null model will randomize the commuity matrix under different metacommunities, separately (and independently).
}
\value{
Output is a list. 
\item{index}{a square matrix of RC (or SES or Confidence based on Bray-Curtis) values.}
\item{sig.index}{character, indicates the index for null model significance test, RC, Confidence, or SES.}
\item{BC.obs}{Output only if output.bray is TRUE. A square matrix of observed taxonomic dissimilarity index (default is Bray-Curtis dissimilarity) values.}
\item{rand}{Output only if detail.null is TRUE. A matrix with all null values of Bray-Curtis index for each turnover.}
}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.

Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. & Inouye, B.D. (2011). Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere, 2, 1-11.

Stegen, J.C., Lin, X., Fredrickson, J.K., Chen, X., Kennedy, D.W., Murray, C.J. et al. (2013). Quantifying community assembly processes and identifying features that impose them. ISME J, 7, 2069.

Kraft, N.J.B., Comita, L.S., Chase, J.M., Sanders, N.J., Swenson, N.G., Crist, T.O. et al. (2011). Disentangling the drivers of beta diversity along latitudinal and elevational gradients. Science, 333, 1755-1758.
}
\author{
Daliang Ning
}
\note{
Version 1: 2021.8.2
}
\seealso{
\code{\link{RC.pc}}
}
\examples{
data("example.data")
comm=example.data$comm
rand.time=20 # usually use 1000 for real data.

# in this example, 10 samples from one metacommunity,
# the other 10 samples from another metacommunity.
meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
rownames(meta.group)=rownames(comm)

nworker=2 # parallel computing thread number
RC=RC.cm(comm=comm, rand = rand.time,
         nworker = nworker, meta.group=meta.group,
         weighted = TRUE, sig.index="RC")
}

\keyword{taxonomic}
