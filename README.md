
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- README.Rmd needs to be rendered to MD with devtools::build_readme() -->

# metaprediction

<!-- badges: start -->
<!-- badges: end -->

The `metaprediction` package provides methods to construct frequentist
predictive distributions in random-effects meta-analysis using
confidence distributions.

The `metaprediction` package introduces the S3 class `metaprediction`.
The methods are implemented in the C++ programming language to ensure
computational efficiency and are made available in the R environment
using the Rcpp interface.

The methodology assume that future effects and effects from observed
studies are exchangeable and distributed according to a normal
distribution. The Edgington combined $p$-value function is used to
construct a confidence distribution for the average effect, and the
generalized heterogeneity statistic is used to construct a confidence
distribution for the heterogeneity parameter.

The considered predictive distributions are then constructed by
marginalizing the normal effect distribution over the latter confidence
distributions to account for uncertainty in the estimation. The
predictive distributions are computed using a Monte Carlo sampling
algorithm, also yielding an estimator for the average effect, which
accounts for uncertainty in the heterogeneity estimate.

## Installation

You can install the development version of `metaprediction` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
remotes::install_github("davidkronthaler-dk/metaprediction")
# or
# install.packages("devtools")
devtools::install_github("davidkronthaler-dk/metaprediction")
```

## Example

The `PredDist` function constructs an object of class `metaprediction`.
Generic S3 methods (such as print, summary, or plot) are applicable,
while all internal helper functions are hidden from the user. The
function `PredDist` requires minimal input: a numeric vector of effect
estimates, a corresponding vector of standard errors, and a string
specifying the desired method. In case exact reproducibility of the
result is required, a seed can be additionally set to ensure the latter.

Assume observed data (estimates and standard errors) from *k = 3*
independent studies. A frequentist predictive distribution can be
computed as:

``` r
library(metaprediction)

# Estimates and standard errors from 5 studies
es <- c(1.17,  2.20,  1.10, -0.0019, -1.33) 
se <- c(0.52, 0.93, 0.63, 0.3, 0.28)

# Computation of predictive distribution
pd <- PredDist(es = es, se = se, method = "FullCD")

# 95% prediction interval
pd$PI
#>      2.5%     97.5% 
#> -3.309570  4.474812

# Plot the predictive distribution
plot(pd, param = "theta_new", breaks = 200, xlim = c(-7, 7))
```

<img src="man/figures/README-example predictive distribution-1.png" width="100%" />

``` r

# Futher information: 'help(PredDist)'
```

The `metaeffect` function can be used to estimate the average effect
$\mu$ and its confidence interval, which incorporates the uncertainty in
the estimation of the between-study heterogeneity:

``` r
# Average effect and confidence interval
me <- metaeffect(es = es, se = se)
#> 
#> ===== Random-Effects Meta-Analysis (Confidence Distributions) =====
#> 
#> Number of studies             : 5 
#> Number of Monte Carlo samples : 100,000 
#> Confidence level              : 95% 
#> 
#> Average effect                : 0.634 
#> Confidence interval           : [ -0.901, 2.121 ]
#> 
#> Summary of confidence distribution of average effect
#>        2.5%     25.0%    Median      Mean    75.0%    97.5%
#>  -0.9010965 0.1800893 0.6418394 0.6335652 1.095079 2.120933
#> 
#> ===================================================================
me$estimate
#> [1] 0.6335652
me$CI
#>       2.5%      97.5% 
#> -0.9010965  2.1209334

# Further information: 'help(metaeffect)'
```
