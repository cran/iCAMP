\name{bNTI.bin.cm}
\alias{bNTI.bin.cm}
\title{
Calculate beta nearest taxon index (betaNTI) for each phylogenetic bin under multiple metacommunities
}
\description{
Perform null model test based on a phylogenetic beta diversity index, beta mean phylogenetic distance to the nearest taxon (betaMNTD), in each bin; calculate beta nearest taxon index (betaNTI), or modified Raup-Crick metric, or confidence level based on the comparison between observed and null betaMNTD in each bin. The package bigmemory (Kane et al 2013) is used to handle very large phylogenetic distance matrix. This function can deal with local communities under different metacommunities (regional pools).
}
\usage{
bNTI.bin.cm(comm, meta.group = NULL, meta.spool = NULL,
            meta.frequency = NULL, meta.ab = NULL,
            pd.desc, pd.spname, pd.wd, pdid.bin, sp.bin,
            spname.check = FALSE, nworker = 4,
            memo.size.GB = 50, weighted = c(TRUE, FALSE),
            rand = 1000, output.bMNTD = c(FALSE, TRUE),
            sig.index = c("SES", "Confidence", "RC", "bNTI"),
            unit.sum = NULL, correct.special = FALSE,
            detail.null = FALSE, special.method = c("MNTD", "MPD", "both"),
            ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975,
            exclude.conspecifics = FALSE, dirichlet = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a species or OTU or gene, thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. rownames are sample IDs. first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.spool}{a list object, each element is a character vector listing all taxa IDs in a metacommunity. The names of the elements indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default is NULL, means to use the observed taxa in comm across samples within the same metacommunity that is defined by meta.group.}
  \item{meta.frequency}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the occurrence frequency of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.frequency as occurrence frequency of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{meta.ab}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the aubndance (or relative abundance) of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.ab as average relative abundance of each taxon in comm across the samples within each metacommunity defined by meta.group.}
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
  \item{dirichlet}{Logic. If TRUE, the taxonomic null model for correcting special cases will use Dirichlet distribution to generate relative abundances in randomized community matrix. default is FALSE.}
}
\details{
This function is particularly designed for samples from different metacommunities. The null model "taxa shuffle" will be done under different metacommunities, separately (and independently). All other details are the same as the function bNTI.bin.big. 
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

Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.

Kane, M.J., Emerson, J., Weston, S. (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning
}
\note{
Version 2: 2022.4.26, fixed error when correcting special case.
Version 1: 2021.8.4
}
\seealso{
\code{\link{bNTI.bin.big}}, \code{\link{icamp.cm}}, \code{\link{bNTI.cm}}, \code{\link{bNTI.big.cm}}
}
\examples{
# function 'bNTI.bin.cm' is usually used in the main function, 'icamp.cm',
# when setting phylo.rand.scale="across",
# means randomization across all bins in phylogenetic null model.

data("example.data")
comm=example.data$comm
tree=example.data$tree
pdid.bin=example.data$pdid.bin
sp.bin=example.data$sp.bin

# in this example, 10 samples from one metacommunity,
# the other 10 samples from another metacommunity.
meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
rownames(meta.group)=rownames(comm)

# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer after change the path for 'save.wd'.
\donttest{
wd0=getwd()
save.wd=paste0(tempdir(),"/pdbig.bNTI.bin.cm")
# you may change save.wd to the folder you want to save the pd.big output.
nworker=2 # parallel computing thread number
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
rand.time=20 # usually use 1000 for real data.

bNTIbins=bNTI.bin.cm(comm=comm,meta.group=meta.group,
                     pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label,
                     pd.wd=pd.big$pd.wd, pdid.bin=pdid.bin, sp.bin=sp.bin,
                     spname.check = TRUE, nworker = nworker, memo.size.GB = 50,
                     weighted = TRUE, rand = rand.time, output.bMNTD = FALSE,
                     sig.index="SES", unit.sum = NULL, correct.special = TRUE,
                     detail.null = FALSE, special.method = "MNTD",
                     exclude.conspecifics = FALSE)
setwd(wd0)
}
}

\keyword{phylogenetic}
