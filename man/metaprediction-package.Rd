\name{metaprediction-package}
\alias{metaprediction-package}
\alias{metaprediction}
\docType{package}
\title{
  Frequentist Predictive Distributions and Intervals for Random-Effects Meta-Analysis
}
\description{
  Provides methods to compute predictive distributions and prediction intervals in random-effects meta-analysis using confidence distributions. The confidence distribution for the average effect is derived from the Edgington combined p-value function (Held et al., 2025), while the confidence distribution for the heterogeneity parameter is based on the generalized heterogeneity statistic (Viechtbauer, 2007). The methods return the full predictive distribution as Monte Carlo samples, enabling the computation of any summary statistic.

  The package also includes an estimator for the average effect based on the Edgington combined p-value function, which accounts for the uncertainty in the estimation of heterogeneity. This estimator is obtained by repeatedly constructing the confidence distribution of the average effect across values of the heterogeneity parameter, sampled according to its own confidence distribution.
  
}
\details{
  \code{PredDist} can be used to compute predictive distributions as well as prediction intervals in a fully frequentist framework.
  
  \code{metaeffect} can be used to estimate the average effect.
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
