\name{dniche}
\alias{dniche}
\title{
Calculate niche difference between species
}
\description{
Calculate niche difference between species based on each environmental variable, directly output the matrix or save the result matrix as big.matrix.
}
\usage{
dniche(env, comm,
       method = c("ab.overlap", "niche.value", "prefer.overlap"),
       nworker = 4, memory.G = 50, out.dist = FALSE,
       bigmemo = TRUE, nd.wd = getwd(), nd.spname.file="nd.names.csv",
       detail.file="ND.res.rda")
}
\arguments{
  \item{env}{matrix or data.frame, each row is a sample, each column is an environmental factor which may be important to represent the niche, thus rownames are sample IDs, and colnames are environmental factor names.}
  \item{comm}{matrix or data.frame, each row is a sample, each column is a spcies (OTU or ASV), thus rownames are sample IDs, colnames are species/OTU/ASV IDs.}
  \item{method}{methods to calculate niche difference. ab.overlap means to calculate from overlapp based on observed abundances along an environment gradient. niche.value means to calculate the difference from abundance weighted mean of each environment factor for each species. prefer.overlap is similar to ab.overlap, but the observed abundances of each species are divided by total abundance sum of the species before calculating overlapping.  If list multiple methods as a vector, only the first element will be used.}
  \item{nworker}{integer, for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memory.G}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{out.dist}{logic, if TRUE, the output niche difference matrix of each environment factor will be a distance object, otherwise will be a matrix in the output list.}
  \item{bigmemo}{logic, if TRUE, big.matrix in R package bigmemory will be used to save each niche differnece matrix as a big matrix on hard disk.}
  \item{nd.wd}{folder path, when bigmemo is TRUE, where the big matrixes are saved.}
  \item{nd.spname.file}{character, name of the file saving taxa IDs, which should be in exactly the same order as in the row names (and column names) of the big niche difference matrix, if bigmemo is TRUE. it should be a .csv file.}
  \item{detail.file}{character, name of the file saving all output information in R data format. it should be a .rda file.}
}
\details{
The method niche.value is to calculate niche difference as the absolute difference of niche values between each pair of species. The niche value of a species is calculated as abundance-weighted mean of each environmental factor as previously reported (Stegen et al 2012 ISME J). In the method ab.overlap, the abundance of each species along the gradient of an environment factor is estimated using the density function using Gaussian kernel with 512 points. Then, the niche difference between two species is calculated as the sum of absolute abundance difference at each point divided by the sum of the higher abundance at each point, like Ruzicka dissimilarity (weighted Jaccard). It is like 1 - niche overlap based on abundance profile overlap, thus called ab.overlap. The method prefer.overlap is very similar to ab.overlap, just one modification, i.e. the observed abundance of each species in each sample is divied by the total abundance of the species across all sample, to normalize the profile, before calcuating niche difference.

Bigmemory (Kane et al 2013) is used to deal with large datasets.
}
\value{
The output is a list object, with several elements.
  \item{bigmemo}{logic, to show whether big.matrix is used.}
  \item{nd}{if bigmemo is FALSE, this is a list of matrixes or distance objects showing the niche difference matrix based on each environment factor. if bigmemo is TRUE, this is a list of big matrix file names.}
  \item{nd.wd}{only appear when bigmemo is TRUE, shows the folder path where the big matrixes are saved.}
  \item{names}{only appear when bigmemo is TRUE, shows species (OTU or ASV) IDs, in the same order as rownames and colnames in the niche difference matrixes.}
  \item{method}{The method used.}
}
\references{
Stegen, J.C., Lin, X., Konopka, A.E. & Fredrickson, J.K. (2012). Stochastic and deterministic assembly processes in subsurface microbial communities. ISME J, 6, 1653-1664.

Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.

Kane, M.J., Emerson, J., Weston, S. (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning
}
\note{
Version 3: 2020.9.1, add nd.spname.file and detail.file; remove setwd; change dontrun to donttest and revise save.wd in help doc.
Version 2: 2020.8.18, add example.
Version 1: 2020.5.15
}

\seealso{
\code{\link{ps.bin}}
}
\examples{
data("example.data")
comm=example.data$comm
env=example.data$env

# if data is small, you do not need to use big.memory
niche.dif=dniche(env = env, comm = comm, method = "niche.value",
                 nworker = 1,out.dist=FALSE,bigmemo=FALSE,nd.wd = NULL)

# if data is large, you need to use big.memory
# since big.memory need to specify a certain folder,
# it is set as 'not test'.
# but you may test the code on your computer after change the path for 'save.wd'.
\donttest{
  save.wd=tempdir() # please change to the folder you want to save the big niche difference matrix.
  nworker=2 # parallel computing thread number
  niche.dif=dniche(env = env, comm = comm,
                   method = "niche.value", nworker = nworker,
                   out.dist=FALSE,bigmemo=TRUE,nd.wd = save.wd)
}
}
\keyword{phylogenetic signal}
