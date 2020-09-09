\name{maxbigm}
\alias{maxbigm}
\title{
Find maximum value in a big matrix
}
\description{
Return the maxium value and its (their) location(s) in a big matrix.
}
\usage{
maxbigm(m.desc, m.wd, nworker = 1, rm.na = TRUE, size.limit = 10000 * 10000)
}
\arguments{
  \item{m.desc}{the name of the file to hold the backingfile description of the big matrix.}
  \item{m.wd}{the path of the folder holding the big matrix file.}
  \item{nworker}{integer, for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{rm.na}{logic, whether to remove NA. Default is TRUE.}
  \item{size.limit}{the matrix size which your current computer memory can easily handle at each time.}
}
\details{
A tool to figure out the maximum value in the big phylogenetic distance matrix.
}
\value{
Output is a list of two elements.
\item{max.value}{Numeric, the maximum value.}
\item{row.col}{Matrix, the row(s) and column(s), i.e. the location(s), of the maximum value in the big matrix.}
}
\references{
Michael J. Kane, John Emerson, Stephen Weston (2013). Scalable Strategies for Computing with Massive Data. Journal of Statistical Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.
}
\author{
Daliang Ning
}
\note{
Version 3: 2020.9.1, remove setwd; change dontrun to donttest and revise save.wd in help doc.
Version 2: 2020.8.19, add example.
Version 1: 2015.12.16
}
\seealso{
\code{\link{midpoint.root.big}}
}
\examples{
# this example shows how to find maximum value
# in a big phylogenetic distance matrix.
data("example.data")
tree=example.data$tree
# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
  save.wd=tempdir() # please change to the folder you want to save the pd.big output.
  nworker=2 # parallel computing thread number
  pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
  
  maxb=maxbigm(m.desc = pd.big$pd.file, m.wd = pd.big$pd.wd,
               nworker = nworker, rm.na = TRUE)
}
}
\keyword{tool}
