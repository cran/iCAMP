\name{cohend}
\alias{cohend}
\title{
Cohen's d effect size
}
\description{
This function is to calculate the popular effect size index Cohen's d.
}
\usage{
cohend(treat, control, paired = FALSE)
}
\arguments{
  \item{treat}{a numeric vector. treatment group.}
  \item{control}{a numberic vector. control group.}
  \item{paired}{logic. Whether the samples in treatment and control groups are paired. default is FALSE.}
}
\details{
This function computes the value of Cohen's d statistics (Cohen 1988).
The effect size magnitude is performed using the thresholds proposed by Cohen (1992), i.e. |d|<0.2 "negligible", 0.2<=|d|<0.5 "small", 0.5<=|d|<0.8 "medium", |d|>=0.8 "large".
The variance of the d is calculate using the conversion formula reportead at page 238 of Cooper et al. (2009):
((n1+n2)/(n1*n2) + .5*d^2/df) * ((n1+n2)/df)
Its square root is output as standard deviation of d.
}
\value{
A list of values will be returned
  \item{d}{Cohen's d value, (mean(treat)-mean(control))/sd}
  \item{sd}{standard deviation of d}
  \item{magnitude}{a qualitative assessment of the magnitude of effect size}
  \item{paired}{whether the samples are paired}
}
\references{
Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). New York:Academic Press

Cohen, J. (1992). A power primer. Psychological Bulletin, 112, 155-159.

The Handbook of Research Synthesis and Meta-Analysis (Cooper, Hedges, & Valentine, 2009)

}
\author{
Daliang Ning
}
\note{
version 1: 2016.2.12
}
\examples{
x=c(1,5,8)
y=c(2,6,10)
cohend(x,y)
cohend(x,y,paired=TRUE)
}
\keyword{ effect size }
