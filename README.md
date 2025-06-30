
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
# install.packages("remotes")
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
#> Predictive & Confidence Distributions for Random-Effects Meta-Analysis
#> 
#> Number of studies: 5 
#> Method: FullCD 
#> Number of Monte Carlo samples: 100,000 
#> 
#> 95% prediction interval from -3.24 to 4.528 
#> 
#> Summary of predictive distribution (PD):
#>               2.5%  25.0% Median  Mean 75.0% 97.5%
#> PD theta_new -3.24 -0.373   0.63 0.633 1.644 4.528
#> 
#> Confidence calculations:
#> Confidence of `theta_new` lying between 0 and Inf:
#> 0.66746
#> Confidence of `theta_new` lying between -Inf and 0:
#> 0.33254
#> 
#> Summary of confidence distributions (CD):
#>           2.5% 25.0% Median  Mean 75.0%  97.5%
#> CD mu   -0.899 0.178  0.640 0.635 1.095  2.166
#> CD tau2  0.365 1.021  1.812 3.246 3.416 14.554

# 95% prediction interval
pd$PI
#>      2.5%     97.5% 
#> -3.239617  4.528168

# Plot the predictive distribution
plot(pd, param = "theta_new", breaks = 200, xlim = c(-7, 7))
```

<img src="man/figures/README-example predictive distribution-1.png" width="100%" />

``` r

# Futher information: 'help(PredDist)'
```

The `remaeffect` function can be used to estimate the average effect
$\mu$ and its confidence interval, which incorporates the uncertainty in
the estimation of the between-study heterogeneity:

``` r
# Average effect and confidence interval
me <- remaeffect(es = es, se = se)
#> 
#> Random-Effects Meta-Analysis using Confidence Distributions
#> 
#> Number of studies: 5 
#> Number of Monte Carlo samples: 100,000 
#> 
#> Average effect: 0.633 
#> 95% Confidence interval from -0.902 to 2.143 
#> 
#> Summary of confidence distribution of the average effect:
#>        2.5%     25.0%    Median      Mean    75.0%    97.5%
#>  -0.9024812 0.1778291 0.6416051 0.6333235 1.091761 2.143147
me$estimate
#> [1] 0.6333235
me$CI
#>       2.5%      97.5% 
#> -0.9024812  2.1431473

# Further information: 'help(remaeffect)'
```
