
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- README.Rmd needs to be rendered to MD with devtools::build_readme() -->

# edgemeta

<!-- badges: start -->

<!-- badges: end -->

The `edgemeta` package provides methods for frequentist random-effects
meta-analysis using Edgington’s method. It contains methods for
estimation of the average effect and prediction of future effects,
constructed using confidence distributions.

The methods are implemented in the C++ programming language to ensure
computational efficiency and are made available in the R environment
using the Rcpp interface.

## Estimation

A combined $p$-value function, or equivalently, confidence distribution,
is constructed by combining one-sided $p$-value functions from
individual studies. The confidence distribution of the average effect is
then marginalized over the confidence distribution of the heterogeneity
parameter to account for heterogeneity estimation uncertainty.
Computation can be performed using a Monte Carlo algorithm generating
independent samples, or using deterministic global adaptive quadrature
integration. The latter is not recommended for few studies (3-5), since
under these marginal distributions are slightly too wide. For mote than
five studies, both methods typically produce very similar results.

## Prediction

The methods assumes that future effects and effects from observed
studies are exchangeable and distributed according to a normal
distribution. This joint effect distribution is maginalized over the
confidence distribution from Edgington’s method and over the confidence
distribution of the heterogeneity parameter, thereby incorporating
uncertainty about parameter estimation in both parameters. The
predictive distributions are computed using a Monte Carlo sampling
algorithm. Implemented are three forms of the predictive distribution: -
Full CD: Uses the full marginalization scheme and is therein the most
accurate (recommended). - Simplified CD: Uses a simplified approach, but
accounts for uncertainty about both parameters. - Fixed $\hat{\tau^2}$:
Does not account for heterogeneity estimation uncertainty and is
typically too narrow except for a large number of studies.

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

The `PredDist` function constructs an object of class `metaprediction`.
Generic S3 methods (such as print, summary, or plot) are applicable. The
function `PredDist` requires minimal input: a numeric vector of effect
estimates, a corresponding vector of standard errors, and a string
specifying the desired method. In case exact reproducibility of the
result is required, a seed can be additionally set to ensure the latter.

Assume observed data (estimates and standard errors) from *k = 5*
independent studies. The average effect can be estimated using the
`remaeffect` function:

``` r
library(edgemeta)

# Estimates and standard errors from 5 studies
es <- c(1.17,  2.20,  1.10, -0.0019, -1.33) 
se <- c(0.52, 0.93, 0.63, 0.3, 0.28)

# Average effect and confidence interval
me <- remaeffect(es = es, se = se, method = "MC") 
#> 
#> CD-Edgington Random-Effects Meta-Analysis
#> 
#> Details:
#> Monte Carlo Algorithm (stochastic, independent samples)
#> Number of Monte Carlo samples: 100,000 
#> 
#> Number of studies: 5 
#> Average effect: 0.637 
#> 95% Confidence interval from -0.902 to 2.165 
#> Two-sided p-value against H0: mu = 0 is 0.35686 
#> 
#> Summary of confidence distribution of the average effect:
#>        2.5%     25.0%    Median      Mean   75.0%    97.5%
#>  -0.9020023 0.1783121 0.6423458 0.6370906 1.09678 2.164598

# or me <- remaeffect(es = es, se = se, method = "GAQ")
# or me <- remaeffect(es = es, se = se, method = "MHEU")

# Further information: 'help(remaeffect)'
```

A frequentist predictive distribution can be computed using:

``` r
# Computation of predictive distribution
pd <- PredDist(es = es, se = se, method = "FullCD")
#> 
#> Predictive and Confidence Distributions for Random-Effects Meta-Analysis
#> 
#> Number of studies: 5 
#> Method: FullCD 
#> Number of Monte Carlo samples: 100,000 
#> 
#> 95% prediction interval from -3.239 to 4.501 
#> 
#> Summary of predictive distribution (PD):
#>                2.5%  25.0% Median  Mean 75.0% 97.5%
#> PD theta_new -3.239 -0.374   0.64 0.633 1.641 4.501
#> 
#> Confidence calculations:
#> Confidence of `theta_new` lying between 0 and Inf:
#> 0.66758
#> Confidence of `theta_new` lying between -Inf and 0:
#> 0.33242
#> 
#> Summary of confidence distributions (CD):
#>           2.5% 25.0% Median  Mean 75.0%  97.5%
#> CD mu   -0.901 0.179  0.643 0.636 1.094  2.142
#> CD tau2  0.365 1.027  1.815 3.226 3.409 14.593

# 95% prediction interval
pd$PI
#>      2.5%     97.5% 
#> -3.239477  4.500845

# Plot the predictive distribution
plot(pd, param = "theta_new", breaks = 200, xlim = c(-7, 7))
```

<img src="man/figures/README-example predictive distribution-1.png" width="100%" />

``` r

# Futher information: 'help(PredDist)'
```
