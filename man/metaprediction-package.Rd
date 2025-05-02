\name{metaprediction-package}
\alias{metaprediction-package}
\alias{metaprediction}
\docType{package}
\title{
  Frequentist Predictive Distributions and Intervals for Meta-Analysis
}
\description{
  Computation of Predictive Distributions and Prediction Intervals utilizing
  the Edgington combined p-value function inducing a confidence distribution of 
  the average effect and accounting for uncertainty about the between-study heterogeneity
  estimate by use of the generalized heterogeneity statistic.
}
\details{
  PredDist can be used to compute both predictive distributions as well as prediction intervals
  in a fully Frequentist framework.
}
\author{
David Kronthaler

Maintainer: David Kronthaler <david.kronthaler@uzh.ch>
}
\references{
  This optional section can contain literature or other references for
  background information.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
es <- c(0.17,  1.20,  1.10, -0.0019, -2.33)
se <- c(0.52, 0.93, 0.63, 0.3, 0.28)
PredDist(es = es, se = se, method = "FullCD")
}
