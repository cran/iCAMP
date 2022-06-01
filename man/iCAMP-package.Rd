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
Version 1.2.10: fix minor bug when output.wd is NULL in icamp.big.
Version 1.2.11: fix minor bug in icamp.big when comm is data.frame.
Version 1.3.1: add bNTI.big and bMNTD.big, and revise qpen to handle big datasets better.
Version 1.3.2: revise icamp.bins to fix error when an input taxonomy name has unrecognizable character; revise icamp.boot to fix error when there is no outlier.
Version 1.3.3: add icamp.cate to summary for each category of taxa, e.g. core versus rare taxa.
Version 1.3.4: typo and format.
Version 1.3.5: revise icamp.big to correct error when using strict bin IDs when omit small bins.
Version 1.4.1: add function 'qpen.test' for bootstrapping test on 'qpen' results.
Version 1.4.2: add options in icamp.big, RC.pc, and RC.bin.bigc to allow relative abundances (value < 1) in community matrix, community data transformation, and use of other taxonomic dissimilarity indexes.
Version 1.4.3: debug to allow input community matrix only has two samples. Also provide a temporary solution for the failure of makeCluster in some OS.
Version 1.4.4: debug ps.bin and icamp.cate to avoid error in special cases.
Version 1.4.5: add option to taxa.binphy.big and icamp.big to handle trees with single edge from root.
Version 1.4.6: debug icamp.big, fix 'differing number of rows' issue in version 1.4.2 to 1.4.5.
Version 1.4.7: speed up qpen.test when there are numerous between-group comparisons.
Version 1.4.8: fix a bug in function taxa.binphy.big.
Version 1.4.9: internal version.
Version 1.4.10: fix a potential bug in function maxbigm.
Version 1.4.11: debug for function icamp.boot.
Version 1.5.1: add functions qpen.cm, RC.cm, bNTI.cm, and bNTI.big.cm, to deal with samples from multiple metacommunities.
Version 1.5.2: add functions icamp.cm, NTI.cm, NRI.cm, bNRI.cm, bNTI.bin.cm, bNRI.bin.cm, and RC.bin.cm, to deal with samples from multiple metacommunities.
Version 1.5.3(20210924): add function pdist.p to calculate phylogenetic distance for relatively small datasets.
Version 1.5.4(20211209): correct 'paste' error in functions bNRI.bin.big and bNRI.bin.cm.
Version 1.5.5(20220210): add icamp.cm2 function to allow different metacommunity settings for taxonomic and phylogenetic null models.
Version 1.5.6(20220410): fix error and warnings from package check.
Version 1.5.7(20220410): fix notes from package check.
Version 1.5.8(20220421): fix error when nworker=1 in several functions.
Version 1.5.9(20220421): fix error in function bNTI.bin.cm.
Version 1.5.10(20220425): correct parallel thread number in examples.
Version 1.5.11(20220426): fix error in special cases in fucntion bNTI.bin.cm.
Version 1.5.12(20220529): fix warnings due to working directory issues of some examples in help documents.
}
\details{
\tabular{ll}{
Package: \tab iCAMP\cr
Type: \tab Package\cr
Version: \tab 1.5.12\cr
Date: \tab 2022-5-29\cr
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
  wd0=getwd() # please change to the folder you want to save the pd.big output.
  save.wd=paste0(tempdir(),"/pdbig")
  nworker=2 # parallel computing thread number
  rand.time=20 # usually use 1000 for real data.
  
  bin.size.limit=5 # for real data, usually use a proper number
  # according to phylogenetic signal test or try some settings
  # then choose the reasonable stochasticity level.
  # our experience is 12, or 24, or 48.
  # but for this example dataset which is too small, have to use 5.
  
  icamp.out=icamp.big(comm=comm,tree=tree,pd.wd=save.wd,
                      rand=rand.time, nworker=nworker,
                      bin.size.limit=bin.size.limit)
  setwd(wd0)
}
}
