\name{qpen.test}
\alias{qpen.test}
\title{
Summary and comparison of QPEN results based on bootstrapping
}
\description{
Bootstrapping analysis of the results from QPEN (quantifying assembly processes based on entire-community null model analysis, Stegen et al 2013, 2015), to estimate the mean and variation of each index and each process influence in each group, and calculate the significance of the difference between groups.
}
\usage{
qpen.test(qpen.result, treat, rand.time = 1000,
          between.group = FALSE, out.detail = TRUE,
          silent = FALSE)
}
\arguments{
  \item{qpen.result}{the output of the function 'qpen', or the element 'result' (data.frame) in qpen output.}
  \item{treat}{matrix or data.frame, each column indicates the group or treatment of each sample, rownames are sample IDs.}
  \item{rand.time}{integer, bootstrapping times. default is 1000.}
  \item{between.group}{logic. if True, the turnovers between each pair of treatments will also be calculated as a group.}
  \item{out.detail}{logic. if True, the 'qpen' results and the bootstrapping results in each group will also be output.}
  \item{silent}{logic. if FALSE, some messages will show during calculation.}
}
\details{
Basically use bootstrapping of samples to estimate the variation of each index's mean and relative importance of each process in each group, as well as the effect size and signficance of the difference between different groups.
}
\value{
Output is a list.
  \item{obs.summary}{the mean, standarded deviation, quartile, and boxplot elements for each observed index (e.g. bMNTD, bNTI, etc.) in each group.}
  \item{boot.summary}{the mean, standarded deviation, quartile, and boxplot elements of the average level of each index and the estimated relative importance of each process in each group.}
  \item{compare}{The relative difference, Cohen's d, and P value of the difference of each index or each process importance between different groups.}
  \item{group.results.detail}{the qpen results in each group.}
  \item{boot.detail}{the average value of each index or estimated relative importanc of each process in each group in each time of bootstrapping.}
}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.

Stegen, J.C., Lin, X., Fredrickson, J.K., Chen, X., Kennedy, D.W., Murray, C.J. et al. (2013). Quantifying community assembly processes and identifying features that impose them. ISME J, 7, 2069.

Stegen, J.C., Lin, X., Fredrickson, J.K. & Konopka, A.E. (2015). Estimating and mapping ecological processes influencing microbial community assembly. Front Microbiol, 6, 370.
}
\author{
Daliang Ning
}
\note{
Version 2: 2021.6.9 speed up transformation for a huge number of comparisons by using package data.table.
Version 1: 2021.4.15 include the function into package iCAMP.
}
\seealso{
\code{\link{qpen}}
}
\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree
treat=example.data$treat
# since pdist.big need to save output to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
  wd0=getwd()
  nworker=2 # parallel computing thread number
  rand.time=5 # usually use 1000 for real data.
  
  # for a big dataset, pdist.big may be used
  save.wd=paste0(tempdir(),"/pdbig.qpen.test")
  # please change to the folder you want to save the pd.big output.
  
  pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
  qp2=qpen(comm=comm, pd=pd.big$pd.file, pd.big.wd=pd.big$pd.wd,
           pd.big.spname=pd.big$tip.label, tree=tree,
           rand.time=rand.time, nworker=nworker)

  qptest=qpen.test(qpen.result=qp2, treat=treat)
  setwd(wd0)
}
}
\keyword{Other methods}
