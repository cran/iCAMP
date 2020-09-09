\name{change.sigindex}
\alias{change.sigindex}
\title{
Change significance index option in iCAMP analysis
}
\description{
This function is to change the method to calculate significance between null and observed dissimilarity and/or change the significance threshold values.
}
\usage{
change.sigindex(icamp.output, sig.index = c("Confidence", "SES.RC", "SES", "RC"),
                detail.save = TRUE, detail.null = FALSE,
                ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975)
}
\arguments{
  \item{icamp.output}{list, the exact output of the function icamp.big in which detail.null must be TRUE, to save all null values.}
  \item{sig.index}{character, Confidence means to directly count the percentage of null values higher/lower than observed value; SES.RC means to use Standard Effect Size (e.g. betaNRI, betaNTI) for phylogenetic beta diversity and use modified Raup-Crick for taxonomic beta diversity, which is typical practice in the previous method; SES means to use Standard Effect Size for both phylogenetic and taxonomic beta diversity; RC means to use modified Raup-Crick for both phylogenetic and taxonomic beta diversity.}
  \item{detail.save}{logic, whether to output the details, including binning information, significance indexes, bin abundances, and some key parameter settings for iCAMP analysis. Default is TRUE}
  \item{detail.null}{logic, whether to output all observed and null values of beta diversity indexes. Default is FALSE.}
  \item{ses.cut}{numeric, the cutoff of significant standard effect size, default is 1.96.}
  \item{rc.cut}{numeric, the cutoff of significant modified Raup-Crick index value, default is 0.95.}
  \item{conf.cut}{numeric, the cutoff of significant confidence level (one-tail), default is 0.975.}
}
\details{
This function is to re-calculate significance using another index or a different threshold value using previously saved null model values. Since the null values are directly extracted from previous icamp.big results, it can skip the most time-consuming step (randomization) and quickly complete calculation.

The default threshold values of Confidence (0.975), SES (1.96), and RC (0.95) mean to capture the 0.95 two-tail confidence level (P=0.05). But, SES need to assume the null values follow normal distribution. RC counts in a half of the special cases that null values are equal to observed values, which is good for obtaining a symmetric metric but theoretically has risk to misestimate significance level (but very slight). Thus, Confidence is preferred as long as the 1000-time randomization is representative.
}
\value{
The output will be the same as icamp.big.
}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals the ecological drivers of grassland soil microbial community assembly in response to warming. bioRxiv. doi:10.1101/2020.02.22.960872 (Nature Communications, final revision).
}
\author{
Daliang Ning
}
\note{
Version 2: 2020.8.18, update help document, add example.
Version 1: 2020.8.1
}

\seealso{
\code{\link{icamp.big}}, \code{\link{null.norm}}
}
\examples{
data("icamp.out")
icamp.out.new=change.sigindex(icamp.output=icamp.out, sig.index = "Confidence")
}
\keyword{test}
