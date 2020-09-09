\name{icamp.bins}
\alias{icamp.bins}
\title{
Summarize iCAMP result in each bin
}
\description{
This function is to calculate various statistic index to assess relative importance of each process in each bin and each turnover, and bin's contribution to each process.
}
\usage{
icamp.bins(icamp.detail, treat = NULL, clas = NULL, silent = FALSE,
          boot = FALSE, rand.time = 1000, between.group = FALSE)
}
\arguments{
  \item{icamp.detail}{list object, the output or the "detail" element of the output from \code{\link{icamp.big}}}
  \item{treat}{matrix or data.frame, indicating the group or treatment of each sample, rownames are sample IDs. Allow to input multi-column matrix, different columns represent different ways to group the samples.}
  \item{clas}{matrix or data.frame, the classification information of species (OTUs).}
  \item{silent}{Logic, whether to show messages. Default is FALSE, thus all messages will be showed.}
  \item{boot}{Logic, whether to do bootstrapping test to get significance of dominating process in each bin.}
  \item{rand.time}{integer, bootstrapping times.}
  \item{between.group}{Logic, whether to analyze between-treatment turnovers.}
}
\details{
Bin level analysis can provide insights into community assembly mechanisms. This function provides more detailed statistics with the output of the main function \code{\link{icamp.big}}.
}
\value{
Output is a list object.
\item{Wtuvk}{The dominant process in each turnover of each bin.}
\item{Ptuv}{Relative importance of each process in governing the turnovers between each pair of communities (samples).}
\item{Ptk}{Relative importance of each process in governing the turnovers of each bin among a group of samples.}
\item{Pt}{Relative importance of each process in governing the turnovers in a group of samples.}
\item{BPtk}{Bin contribution to each process, measuring the contribution of each bin to the relative importance of each process in the assembly of a group of communities.}
\item{BRPtk}{Bin relative contribution to each process, measuring the relative contribution of each bin to a certain process.}
\item{Binwt}{Output if treat is given. Bin relative abundance in each group (treatment) of samples.}
\item{Bin.TopClass}{Output if clas is given. A matrix showing the bin relative abundance; the top taxon ID, percentage in bin, and classification; the most abundant name at each phylogeny level in the bin.}
\item{Class.Bin}{Output if clas is given. A matrix showing the bin ID and classification information for each taxon.}
}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals the ecological drivers of grassland soil microbial community assembly in response to warming. bioRxiv. doi:10.1101/2020.02.22.960872 (Nature Communications, final revision).
}
\author{
Daliang Ning
}
\note{
Version 2: 2020.8.19, update help document, add example.
Version 1: 2019.12.11
}
\seealso{
\code{\link{icamp.big}}
}
\examples{
data("icamp.out")
data("example.data")
treatment=example.data$treat
classification=example.data$classification
rand.time=20 # usually use 1000 for real data.
icampbin=icamp.bins(icamp.detail = icamp.out, treat = treatment,
                    clas = classification, boot = TRUE,
                    rand.time = rand.time, between.group = TRUE)
}

\keyword{main}
