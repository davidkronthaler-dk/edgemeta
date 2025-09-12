\name{edgemeta-package}
\alias{edgemeta-package}
\alias{edgemeta}
\docType{package}
\title{
  Edgington's Method for Random-Effects Meta-Analysis: Estimation and Prediction
}
\description{
  Provides methods for estimation and prediction in random-effects meta-analysis using Edgington's method
  for combination of p-value functions, or equivalently, confidence distributions.
}
\details{
  \code{remaeffect} can be used to estimate the average effect.

  \code{PredDist} can be used to compute predictive distributions as well as prediction intervals.

  \code{crps} provides a Monte Carlo approximation to compute continuous ranked probability scores for computational efficiency.
  
  \code{cdtau2} provides a confidence density of the heterogeneity paramter \eqn{\tau^2} derived from the generalized heterogeneity statistic.
}
\author{
David Kronthaler

Maintainer: David Kronthaler <david.kronthaler@uzh.ch>
}
\references{
Viechtbauer, W. (2007). \emph{Confidence intervals for the amount of heterogeneity in meta‐analysis}. Statistics in medicine, 26(1), 37-52. https://doi.org/10.1002/sim.2514

Held, L., Hofmann, F., & Pawel, S. (2025). \emph{A comparison of combined p-value functions for meta-analysis}. doi:10.1017/rsm.2025.26
}
\keyword{ package }
\seealso{
 \code{\link{PredDist}}
  \code{\link{remaeffect}}
}
