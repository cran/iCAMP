\name{match.2col}
\alias{match.2col}
\title{
Check the consistency of the first two columns of different matrixes 
}
\description{
This function is usually used to check the consistency of samples names in different pairwise comparison matrixes.
}
\usage{
match.2col(check.list, name.check = NULL, rerank = TRUE, silent = FALSE)
}
\arguments{
  \item{check.list}{List, each element is a matrix. It must be set in a format like "check.list=list(A=A,B=B)".
  The first two columns of the matrixes will be compared and matched with each other. }
  \item{name.check}{matrix, the first two columns will be used as a standard. The pairs not appear in this matrix will be removed from all matrixes.}
  \item{rerank}{Logic, make the first two columns in all matrixes in the same rank or not. Default is TRUE.}
  \item{silent}{Logic, whether to show messages. Default is FALSE, thus all messages will be showed.}
}
\details{
A tool to match IDs.
}
\value{
Return a list object, new matrixes with the same first two columns. Some messages will return if some names are removed or all names matches very well.
}
\author{
Daliang Ning
}
\note{
Version 2: 2020.8.19, add example.
Version 1: 2018.10.20
}
\examples{
# here two simple matrixes are generated and the pairwise comparison IDs not matched are removed.
A=1:5
names(A)=paste0("S",1:5)
B=1:6
names(B)=paste0("S",1:6)
DA3c=dist.3col(dist(A))
DB3c=dist.3col(dist(B))

checkid=match.2col(check.list = list(DA3c=DA3c,DB3c=DB3c))
DA3cnew=checkid$DA3c
DB3cnew=checkid$DB3c
}

\keyword{tool}
