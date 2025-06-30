
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- README.Rmd needs to be rendered to MD with devtools::build_readme() -->

# metaprediction

<!-- badges: start -->
<!-- badges: end -->

The `metaprediction` package provides methods to construct frequentist
predictive distributions in random-effects meta-analysis using
confidence distributions.

The package introduces the S3 class `metaprediction`. The methods are
implemented in the C++ programming language to ensure computational
efficiency and are made available in the R environment using the Rcpp
interface.

The methods assumes that future effects and effects from observed
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
Generic S3 methods (such as print, summary, or plot) are applicable. The
function `PredDist` requires minimal input: a numeric vector of effect
estimates, a corresponding vector of standard errors, and a string
specifying the desired method. In case exact reproducibility of the
result is required, a seed can be additionally set to ensure the latter.

Assume observed data (estimates and standard errors) from *k = 5*
independent studies. A frequentist predictive distribution can be
computed using:

``` r
library(metaprediction)

# Estimates and standard errors from 5 studies
es <- c(1.17,  2.20,  1.10, -0.0019, -1.33) 
se <- c(0.52, 0.93, 0.63, 0.3, 0.28)

# Computation of predictive distribution
pd <- PredDist(es = es, se = se, method = "FullCD")
#> 
#> =================== MetaPrediction Summary ==================
#> 
#> Number of studies: 5 
#> Method:  FullCD 
#> Number of Monte Carlo samples: 100000 
#> 
#> 95 % prediction interval from -3.285 to 4.502 
#> 
#> Summary of predictive distribution (PD) and confidence distributions (CD)
#>                    2.5%      25.0%    Median      Mean    75.0%    97.5%
#> PD theta_new -3.2845101 -0.3922720 0.6333064 0.6252145 1.639966  4.50216
#> CD mu        -0.8879358  0.1781743 0.6365155 0.6321113 1.089866  2.12975
#> CD tau2       0.3672642  1.0246826 1.8132607 3.2422750 3.413068 14.74459
#> 
#> Confidence calculations:
#> Confidence of `theta_new` lying between 0 and Inf :
#> 0.66492
#> Confidence of `theta_new` lying between -Inf and 0 :
#> 0.33508
#> 
#> =============================================================

# 95% prediction interval
pd$PI
#>     2.5%    97.5% 
#> -3.28451  4.50216

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
#> Average effect                : 0.632 
#> Confidence interval           : [ -0.903, 2.133 ]
#> 
#> Summary of confidence distribution of average effect
#>        2.5%     25.0%    Median      Mean   75.0% 97.5%
#>  -0.9033486 0.1768678 0.6374315 0.6315653 1.09519 2.133
#> 
#> ===================================================================
me$estimate
#> [1] 0.6315653
me$CI
#>       2.5%      97.5% 
#> -0.9033486  2.1330003

# Further information: 'help(metaeffect)'
```
