\name{qpen.cm}
\alias{qpen.cm}
\title{
Quantifying assembly processes based on entire-community null model analysis under multiple metacommunities
}
\description{
The framework quantifying assembly processes based on entire-community null model analysis (Stegen et al 2013, 2015). add bigmemory function to handle large datasets. This function can deal with local communities under different metacommunities (regional pools).
}
\usage{
qpen.cm(comm, pd = NULL, pd.big.wd = NULL,
        pd.big.spname = NULL, tree = NULL,
        meta.group = NULL, meta.com = NULL,
        meta.frequency = NULL, meta.ab = NULL,
        ab.weight = TRUE, exclude.conspecifics = FALSE,
        rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
        nworker = 4, memory.G = 50, project = NA, wd = getwd(),
        output.detail = FALSE, save.bNTIRC = FALSE,
        taxo.metric = "bray", transform.method = NULL,
        logbase = 2, dirichlet = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a taxon (species or OTU or ASV), thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{pd}{a character or a matrix or NULL. If it is a character, it specifies the name of the file to hold the backingfile description of the big phylogenetic distance matrix, it is usually "pd.desc" if using default setting in pdist.big function. If it is matrix, it is the phylogenetic distance matrix with taxa IDs as colnames and rownames. The function will check the consistency of species names between phylogenetic matrix and comm. if pd is NULL, the function will calculate phylogenetic distance matrix from tree.}
  \item{pd.big.wd}{a folder path, where the bigmemmory file of the phylogenetic distance matrix are saved. If pd has given the phylogenetic matrix, pd.big.wd should be NULL. If pd is NULL and pd.big.wd is NULL either, a folder will be created to save the big phylgoenetic distance matrix.}
  \item{pd.big.spname}{character vector, taxa id in the same rank as the big matrix of phylogenetic distances. not necessary if pd has given the phylogenetic matrix.}
  \item{tree}{phylogenetic tree, an object of class "phylo".}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. Rownames are sample IDs. The first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.com}{a list object, each element is a matrix or data.frame to define abundance (or relative abundance) of taxa in a metacommunity (regional pool). The element names indicate metacommunity names, which should be consistent with the metacommunity names defined in meta.group. If there is only one metacommunity, meta.com can be a matrix or data.frame to define taxa abundance (or relative abundance) in the metacommunity. Default is NULL, means to calculate metacommunity structure from comm according to metacommunities defined in meta.group.
}
  \item{meta.frequency}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the occurrence frequency of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.frequency as occurrence frequency of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{meta.ab}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the aubndance (or relative abundance) of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.ab as average relative abundance of each taxon in comm across the samples within each metacommunity defined by meta.group.}
  \item{ab.weight}{logic, abundance-weighted or binary. default is TRUE, means abundance-weighted.}
  \item{exclude.conspecifics}{Logic, should conspecific taxa in different communities be exclude from MNTD calculations? default is FALSE. The same as in the function bmntd.}
  \item{rand.time}{integer, randomization times. default is 1000.}
  \item{sig.bNTI}{numeric, the cutoff of significant betaNTI, default is 1.96.}
  \item{sig.rc}{numeric, the cutoff of significant modified Raup-Crick metric, default is 0.95.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memory.G}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{project}{character string, the prefix of saved output files.}
  \item{wd}{a folder path, where the files will be saved.}
  \item{output.detail}{logic, if TRUE, some detailed results including all null values will be output.}
  \item{save.bNTIRC}{logic, if TRUE, a file will be saved in the folder specified by wd.}
  \item{taxo.metric}{taxonomic beta diversity index, the same as 'method' in the function 'vegdist' in package 'vegan', including "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis". If taxo.metric='bray' and transform.method=NULL, RC will be calculated based on Bray-Curtis dissimilarity as recommended in original iCAMP; otherwise, unit.sum setting will be ignored.}
  \item{transform.method}{character or a defined function, to specify how to transform community matrix before calculating taxonomic dissimilarity. Due to the definition of the phylogenetic dissimilarity (bMNTD), it is not affected by data transformation. If transform.method is a characher, it should be a method name as in the function 'decostand' in package 'vegan', including 'total','max','freq','normalize','range','standardize','pa','chi.square','cmdscale','hellinger','log'.}
  \item{logbase}{numeric, the logarithm base used when transform.method='log'.} 
  \item{dirichlet}{Logic. If TRUE, the taxonomic null model will use Dirichlet distribution to generate relative abundances in randomized community matrix. If the input community matrix has all row sums no more than 1, the function will automatically set dirichlet=TRUE. default is FALSE.}
}
\details{
This function is particularly designed for samples from different metacommunities. The null model will randomize the commuity matrix under different metacommunities, separately (and independently). All other details are the same as the function qpen. 
}
\value{
output is a list with following elements.
\item{ratio}{The overall percentage of turnovers governed by each process}
\item{result}{a matrix, showing betaMNTD, Bray-Curtis dissimilarity, betaNTI, RC, and the identified process governing each turnover. Each row represents a turnover between two samples.}
\item{pd}{phylogenetic distance matrix or the backingfile name if phylogenetic distance is saved by bigmemeory function.}
\item{bMNTD}{a square matrix of observed betaMNTD.}
\item{BC}{a square matrix of observed Bary-Curtis index.}
\item{bMNTD.rand}{a matrix showing all null values of betaMNTD.}
\item{BC.rand}{a matrix showing all null values of Bray-Curtis index.}
\item{setting}{a data.frame showing all settings for this function.}
}
\references{
Stegen, J.C., Lin, X., Fredrickson, J.K., Chen, X., Kennedy, D.W., Murray, C.J. et al. (2013). Quantifying community assembly processes and identifying features that impose them. ISME J, 7, 2069.

Stegen, J.C., Lin, X., Fredrickson, J.K. & Konopka, A.E. (2015). Estimating and mapping ecological processes influencing microbial community assembly. Front Microbiol, 6, 370.

Kane, M.J., Emerson, J., & Weston, S. (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning
}
\note{
Version 1: 2021.8.3
}
\seealso{
\code{\link{qpen}}, \code{\link{RC.cm}}, \code{\link{bNTI.cm}}, \code{\link{bNTI.big.cm}}
}
\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree

# in this example, 10 samples from one metacommunity,
# the other 10 samples from another metacommunity.
meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
rownames(meta.group)=rownames(comm)

# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
  wd0=getwd()
  nworker=2 # parallel computing thread number
  rand.time=5 # usually use 1000 for real data.
  
  # for a small dataset, phylogenetic distance matrix can be directly used.
  pd=example.data$pd
  qp=qpen.cm(comm=comm, meta.group=meta.group, pd=pd,
          rand.time=rand.time,nworker=nworker)
  
  # for a big dataset, pdist.big may be used
  save.wd=paste0(tempdir(),"/pdbig.qpen.cm")
  # please change to the folder you want to save the pd.big output.
  
  pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
  qp2=qpen.cm(comm=comm, meta.group=meta.group,
           pd=pd.big$pd.file, pd.big.wd=pd.big$pd.wd,
           pd.big.spname=pd.big$tip.label, tree=tree,
           rand.time=rand.time, nworker=nworker)
  setwd(wd0)
}
}
\keyword{Other methods}
