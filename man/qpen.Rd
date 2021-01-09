\name{qpen}
\alias{qpen}
\title{
Quantifying assembly processes based on entire-community null model analysis
}
\description{
The previous framework quantifying assembly processes based on entire-community null model analysis (Stegen et al 2013, 2015). add bigmemory function to handle large datasets.
}
\usage{
qpen(comm = NULL, pd = NULL, pd.big.wd = NULL,
     pd.big.spname = NULL, tree = NULL,
     bNTI = NULL, RC = NULL, ab.weight = TRUE,
     meta.ab = NULL, exclude.conspecifics = FALSE,
     rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
     nworker = 4, memory.G = 50, project = NA,
     wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a taxon (species or OTU or ASV), thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{pd}{a character or a matrix or NULL. If it is a character, it specifies the name of the file to hold the backingfile description of the big phylogenetic distance matrix, it is usually "pd.desc" if using default setting in pdist.big function. If it is matrix, it is the phylogenetic distance matrix with taxa IDs as colnames and rownames. The function will check the consistency of species names between phylogenetic matrix and comm. if pd is NULL, the function will calculate phylogenetic distance matrix from tree.}
  \item{pd.big.wd}{a folder path, where the bigmemmory file of the phylogenetic distance matrix are saved. If pd has given the phylogenetic matrix, pd.big.wd should be NULL. If pd is NULL and pd.big.wd is NULL either, a folder will be created to save the big phylgoenetic distance matrix.}
  \item{pd.big.spname}{character vector, taxa id in the same rank as the big matrix of phylogenetic distances. not necessary if pd has given the phylogenetic matrix.}
  \item{tree}{phylogenetic tree, an object of class "phylo".}
  \item{bNTI}{a matrix, if beta nearest taxon index (betaNTI) values are available, just directly input them as a squre matrix here.}
  \item{RC}{a matrix, if modified Raup-Crick metric based Bray-Curtis (RCbray) values are available, just directly input them as a squre matrix here.}
  \item{ab.weight}{logic, abundance-weighted or binary. default is TRUE, means abundance-weighted.}
  \item{meta.ab}{a numeric vector, to define the relative aubndance of each species in the regional pool. Default setting is NULL, means to calculate meta.ab as average relative abundance of each species across the samples.}
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
}
\details{
This method is developed by Dr. James Stegen (2013 and 2015), from which we developed iCAMP. In all pairwised comparisons between samples/communities, the non-random phylogenetic turnovers recognized by phylogeny shuffle were counted as influence of environment selection (betaNTI>1.96 or <-1.96), and the non-random taxonomic turnovers in the rest pairwised comparisons were counted as influence of dispersal limitation (RCbray>0.95) or homogenizing dispersal (RCbray<-0.95). The rest part is called undominated. Please read the references for details.

betaNTI is the standardized effect size between the observed and null values of beta mean nearest taxon distance (betaMNTD). RCbray is the modified Roup-Crick metric based on taxonomic Bray-Curtis dissimilarity index.

Bigmemory (Kane et al 2013) is used to deal with large datasets.
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
Version 6: 2020.12.5, use function bNTI.big to calculate betaNTI, which use bigmemory better.
Version 5: 2020.9.1, remove setwd; change dontrun to donttest and revise save.wd in help doc.
Version 4: 2020.8.21, update help document, add example.
Version 3: 2016.2.15, add options to use bigmemory to handle large phylogenetic distance matrixes.
Version 2: 2015.11.22
}

\seealso{
\code{\link{bNTIn.p}}, \code{\link{bNTI.big}}, \code{\link{RC.pc}}
}

\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree
# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
  nworker=2 # parallel computing thread number
  rand.time=5 # usually use 1000 for real data.
  
  # for a small dataset, phylogenetic distance matrix can be directly used.
  pd=example.data$pd
  qp=qpen(comm=comm, pd=pd,
          rand.time=rand.time,nworker=nworker)
  
  # for a big dataset, pdist.big may be used
  save.wd=tempdir() # please change to the folder you want to save the pd.big output.
  pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
  qp2=qpen(comm=comm, pd=pd.big$pd.file, pd.big.wd=pd.big$pd.wd,
           pd.big.spname=pd.big$tip.label, tree=tree,
           rand.time=rand.time, nworker=nworker)
}
}
\keyword{Other methods}% use one of  RShowDoc("KEYWORDS")

