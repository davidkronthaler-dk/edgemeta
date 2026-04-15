
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- README.Rmd needs to be rendered to MD with devtools::build_readme() -->

# edgemeta

<!-- badges: start -->

<!-- badges: end -->

The `edgemeta` package provides methods for frequentist random-effects
meta-analysis using Edgington’s method. It contains methods for
estimation of the average effect and prediction of future study effects,
constructed using confidence distributions.

The methods are implemented in the C++ programming language to ensure
computational efficiency and are made available in the R environment
using the Rcpp interface.

## Estimation

A combined $p$-value function, or equivalently, confidence distribution,
is constructed by combining one-sided $p$-value functions from
individual studies using Edgington’s $p$-value combination method. The
combined $p$-value function of the average effect is a confidence
distribution of the average effect conditional on the heterogeneity
parameter. To account for heterogeneity estimation uncertainty, it is
then marginalized over the confidence distribution of the heterogeneity
parameter. Computation can be performed using a Monte Carlo algorithm
generating independent samples, or using deterministic global adaptive
quadrature integration. The latter is not recommended for few studies
(3-5), since under these marginal distributions are slightly too wide.
For more than five studies, both methods typically produce very similar
results.

## Prediction

The methods assumes that future effects and effects from observed
studies are exchangeable and distributed according to a normal
distribution. This joint effect distribution is maginalized over the
confidence distribution from Edgington’s method and over the confidence
distribution of the heterogeneity parameter, thereby incorporating
uncertainty about parameter estimation in both parameters. The
predictive distributions are computed using a Monte Carlo sampling
algorithm. Implemented are three forms of the predictive distribution: -
PCD-full: Uses the full marginalization scheme and is therein the most
accurate (recommended). - PCD-simplified: Uses a simplified approach,
but accounts for uncertainty about both parameters. - PCD-fixed: Does
not account for heterogeneity estimation uncertainty and is typically
too narrow except for a large number of studies.

## Installation

You can install the `edgemeta` package from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("davidkronthaler-dk/edgemeta")
# or
# install.packages("devtools")
devtools::install_github("davidkronthaler-dk/edgemeta")
```

## Example

Assume observed data (estimates and standard errors) from *k = 5*
independent studies. The average effect can be estimated using the
`remaeffect` function:

``` r
library(edgemeta)

# Estimates and standard errors from 5 studies
es <- c(1.17,  2.20,  1.10, -0.0019, -1.33) 
se <- c(0.52, 0.93, 0.63, 0.3, 0.28)

# Average effect and confidence interval
remaeffect(es = es, se = se, method = "MC") 
#> 
#> CD-Edgington Random-Effects Meta-Analysis
#> 
#> Details:
#> Monte Carlo Algorithm (stochastic, independent samples)
#> Number of Monte Carlo samples: 100,000 
#> 
#> Number of studies: 5 
#> Average effect: 0.630 
#> 95% Confidence interval from -0.901 to 2.127 
#> Two-sided p-value against H0: mu = 0 is 0.3579 
#> 
#> Summary of confidence distribution of the average effect:
#>        2.5%     25.0%    Median      Mean    75.0%    97.5%
#>  -0.9012554 0.1774378 0.6353253 0.6302936 1.091239 2.127029

# or: remaeffect(es = es, se = se, method = "GAQ")
# or: remaeffect(es = es, se = se, method = "MHEU")

# Further information: 'help(remaeffect)'
```

A frequentist predictive distribution can be computed using:

``` r
# Computation of predictive distribution
pd <- PredDist(es = es, se = se, method = "PCD-full")
#> 
#> Predictive Distribution for Random-Effects Meta-Analysis
#> 
#> Number of studies: 5 
#> Method: PCD-full 
#> Number of Monte Carlo samples: 100,000 
#> 
#> 95% prediction interval from -3.238 to 4.523 
#> 
#> Summary of predictive distribution:
#>    2.5%  25.0% Median  Mean 75.0% 97.5%
#>  -3.238 -0.386  0.637 0.629 1.645 4.523

# or: PredDist(es = es, se = se, method = "PCD-simplified")
# or: PredDist(es = es, se = se, method = "PCD-fixed")

# 95% prediction interval
pd$PI
#>      2.5%     97.5% 
#> -3.238030  4.522863

# Plot the predictive distribution
# plot(pd, param = "theta_new", breaks = 200, xlim = c(-7, 7))

# Futher information: 'help(PredDist)'
```
