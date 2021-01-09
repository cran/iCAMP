\name{icamp.cate}
\alias{icamp.cate}
\title{
Summarize iCAMP result for different categories of taxa
}
\description{
This function is to calculate various statistic index to assess relative importance of each process on different categories of taxa. The categories can be defined in various ways. For example, core, consistently and occasionally rare taxa; or different phyla; or various particular functional groups.
}
\usage{
icamp.cate(icamp.bins.result, comm, cate, treat = NULL,
           silent = FALSE, between.group = FALSE)
}
\arguments{
  \item{icamp.bins.result}{list object, the output from \code{\link{icamp.bins}}}
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a taxon (a species or OTU or ASV), thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{cate}{matrix or data.frame, indicating the category of each taxon, rownames are taxa IDs. If the matrix has multiple columns, only the first column will be used.}
  \item{treat}{matrix or data.frame, indicating the group or treatment of each sample, rownames are sample IDs. Allow to input multi-column matrix, different columns represent different ways to group the samples.}
  \item{silent}{logic, if FALSE, some messages will show during calculation.}
  \item{between.group}{logic, whether to analyze between-treatment turnovers.}
}
\details{
This function simply sums up the relative abundance of taxa of a category in different bins governed by a process to summarize the relative importance of the process on the category.
}
\value{
Output is a list object.
\item{Ptuvx}{Relative importance of each process in governing each category's turnover between each pair of communities (samples).}
\item{Ptx}{Relative importance of each process in governing each category's turnovers among a group of samples.}
}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.
}
\author{
Daliang Ning
}
\note{
Version 2: 2021.1.7, add help document; fixed NAN error.
Version 1: 2020.12.9.
}
\seealso{
\code{\link{icamp.bins}}
}
\examples{
data("icamp.out")
data("example.data")
comm=example.data$comm
treatment=example.data$treat
classification=example.data$classification
rand.time=20 # usually use 1000 for real data.
# 1 # summarize each bin
icampbin=icamp.bins(icamp.detail = icamp.out, treat = treatment,
                    clas = classification, boot = TRUE,
                    rand.time = rand.time, between.group = TRUE)

# 2 # define category
cate=data.frame(type=rep("others",ncol(comm)),stringsAsFactors = FALSE)
rownames(cate)=colnames(comm)
tax.frequency=colSums(comm>0)/nrow(comm)
tax.relative.ab=colMeans(comm/rowSums(comm))
cate[which(tax.frequency>0.75 & tax.relative.ab>0.05),1]="core"
cate[which(tax.relative.ab<0.02),1]="rare"

# 3 # summarize each category
icampcate=icamp.cate(icamp.bins.result = icampbin, comm = comm, cate = cate,
                     treat = treatment, silent = FALSE, between.group = TRUE)
}

\keyword{main}
