\name{icamp.boot}
\alias{icamp.boot}
\title{
Bootstrapping analysis of icamp results
}
\description{
Use bootstrapping to estimate the variation of relative importance of each process in each group, and compare the difference between groups.
}
\usage{
icamp.boot(icamp.result, treat, rand.time = 1000, compare = TRUE,
           silent = FALSE, between.group = FALSE, ST.estimation = FALSE)
}
\arguments{
  \item{icamp.result}{data.frame object, from the output of \code{\link{icamp.big}}. the first two columns are sample IDs, the third to seventh columns are the relative importance of the five ecological processes.}
  \item{treat}{matrix or data.frame, a one-column (n x 1) matrix indicating the group or treatment of each sample, rownames are sample IDs. if input a n x m matrix, only the first column is used.}
  \item{rand.time}{integer, bootstrapping times. default is 1000.}
  \item{compare}{logic, whether to compare icamp reults between different groups.}
  \item{silent}{logic, if FALSE, some messages will show during calculation.}
  \item{between.group}{logic, whether to analyze between-treatment turnovers.}
  \item{ST.estimation}{logic, whether to estimate stochasticity as the total relative importance of dispersal and drift.}
}
\details{
Bootstrapping is implemented by random draw samples with replacement, to estimate the variation of relative importance of each process in each group, and calculate the relative difference, effect size, and significance of the difference between each two groups.}
\value{
Output is a list with three elements.
\item{summary}{data.frame, summary of each group.

Group: group name from the input "treat".

Process: process name from the icamp.result.

Observed: the mean relative importance of each process in each group.

Mean, Stdev, Min, Quartile25, Median, Quartile75, and Max: mean, standard deviation, minimum, 25 percent-quantile, median, 75 percent-quantile, and maximum of bootstrapping results, respectively.

Lower.whisker, Lower.hinge, Mediean.1, Higher.hinge, Higher.whisker, Outerlier1...: boxplot elements.}

\item{compare}{data.frame, summary of comaprison between each two groups. First two columns are group names. From the third column, different indexes for comparison are showed, including Cohen's d (Cohen.d), effect size magnitude according to Cohen's d (Effect.Size), and P value from bootstrapping test (P.value).}

\item{boot.detail}{a list of matrixes, each matrix corresponds to a group, showing detailed bootstrapping results in each random draw.}

}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.
}
\author{
Daliang Ning
}
\note{
Version 4: 2021.7.1, fix a bug leading to zero cohen's d.
Version 3: 2021.1.5, fix error when there is no outlier.
Version 2: 2020.8.19, update help document, add example.
Version 1: 2019.11.14
}
\seealso{
\code{\link{icamp.big}}
}
\examples{
data("icamp.out")
data("example.data")
treatment=example.data$treat
rand.time=20 # usually use 1000 for real data.
icampbt=icamp.boot(icamp.result = icamp.out$bNRIiRCa, treat = treatment,
                   rand.time = rand.time, compare = TRUE,
                   between.group = TRUE, ST.estimation = TRUE)
}
\keyword{main}
