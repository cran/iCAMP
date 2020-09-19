\name{iCAMP-package}
\alias{iCAMP-package}
\alias{iCAMP}
\docType{package}
\title{
Infer Community Assembly Mechanisms by Phylogenetic-bin-based null model analysis
}
\description{
This package is to implement a general framework to quantitatively infer Community Assembly Mechanisms by Phylogenetic-bin-based null model analysis, abbreviated as iCAMP (Ning et al 2020). It can quantitatively assess the relative importance of different community assembly processes, such as selection, dispersal, and drift, for both communities and each phylogenetic group ('bin'). Each bin usually consists of different taxa from a family or an order. The package also provides functions to implement some other published methods, including neutral taxa percentage (Burns et al 2016) based on neutral theory model (Sloan et al 2006) and quantifying assembly processes based on entire-community null models (Stegen et al 2013). It also includes quite a few handy functions, particularly for big datasets, such as phylogenetic and taxonomic null model analysis at both community and bin levels, between-taxa niche difference and phylogenetic distance calculation, phylogenetic signal test within phylogenetic groups, midpoint root of big trees, etc. URL: https://github.com/DaliangNing/iCAMP1

Version 1.2.4: the first formal version of iCAMP for CRAN.
Version 1.2.5: correct typo in description and fix the error of memory.limit issue.
Version 1.2.6: revise the help document of qpen to include an example for big datasets.
Version 1.2.7: remove setwd in functions; add options to specify file names; change dontrun to donttest and revise save.wd in some help documents.
Version 1.2.8: revise dniche to avoid unnecessary file.
Version 1.2.9: update iCAMP paper newly published on Nature Communications and the GitHub link.
}
\details{
\tabular{ll}{
Package: \tab iCAMP\cr
Type: \tab Package\cr
Version: \tab 1.2.9\cr
Date: \tab 2020-09-19\cr
License: \tab GPL-2\cr
}
}
\author{
Daliang Ning <ningdaliang@ou.edu>
}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.

Burns, A.R., Stephens, W.Z., Stagaman, K., Wong, S., Rawls, J.F., Guillemin, K. et al. (2016). Contribution of neutral processes to the assembly of gut microbial communities in the zebrafish over host development. Isme Journal, 10, 655-664.

Sloan, W.T., Lunn, M., Woodcock, S., Head, I.M., Nee, S. & Curtis, T.P. (2006). Quantifying the roles of immigration and chance in shaping prokaryote community structure. Environmental Microbiology, 8, 732-740.

Stegen, J.C., Lin, X., Fredrickson, J.K., Chen, X., Kennedy, D.W., Murray, C.J. et al. (2013). Quantifying community assembly processes and identifying features that impose them. Isme Journal, 7, 2069-2079.
}
\keyword{ package }
\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree
# since need to save some outputs to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the path for 'save.wd'.
\donttest{
  save.wd=tempdir() # please change to the folder you want to save the pd.big output.
  pd.wd=paste0(save.wd,"/pdbig")
  nworker=4 # parallel computing thread number
  rand.time=20 # usually use 1000 for real data.
  
  bin.size.limit=5 # for real data, usually use a proper number
  # according to phylogenetic signal test or try some settings
  # then choose the reasonable stochasticity level.
  # our experience is 12, or 24, or 48.
  # but for this example dataset which is too small, have to use 5.
  
  setwd(save.wd)
  icamp.out=icamp.big(comm=comm,tree=tree,pd.wd=pd.wd,
                      rand=rand.time, nworker=nworker,
                      bin.size.limit=bin.size.limit)
}
}
