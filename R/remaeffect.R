#' Estimate the Average Effect in Random-Effects Meta-Analysis via Edgington's Confidence Distribution
#'
#' Confidence distributions for the average effect \eqn{\mu}, or equivalently, \eqn{p}-value functions, from individual studies
#' are combined using Edgington's method. The resulting combined confidence distribution of \eqn{\mu}, 
#' which is conditional on the heterogeneity parameter \eqn{\tau^2}, is integrated over an approximate confidence distribution of the latter
#' to account for heterogeneity estimation uncertainty.
#'
#' @param es Numeric vector of effect estimates from individual studies (length >= 2).
#' @param se Numeric vector of standard errors corresponding to each effect estimate (length >= 2).
#' @param method Either "MC" for Monte Carlo sampling algorithm or "GAQ" for global adaptive quadrature integration (compare details).
#' @param w Study-specific weights.
#' @param level.ci Confidence level for the interval estimate (default is 0.95).
#' @param B Number of Monte Carlo samples used to estimate the confidence distributions (default is 100,000). Relevant when selecting \code{method = "MC"}.
#' @param seed Optional integer to ensure reproducibility of the random sampling. Relevant when selecting \code{method = "MC"}.
#' @param mu0 Parameter value under the null-hypothesis against which two-sided p-value is computed.
#'
#' @return A list, including not all, but some of these, depending on the method:
#' \describe{
#'   \item{estimate}{Point estimate of the average effect \eqn{\mu}.}
#'   \item{CI}{Confidence interval for the average effect \eqn{\mu}.}
#'   \item{pval}{Two-sided p-value against H0: \eqn{\mu} = \eqn{\mu_0}.}
#'   \item{cd_mu}{Vector of samples from the confidence distribution of the average effect \eqn{\mu}. Obtained for method \code{MC}}.
#'   \item{cd_tau2}{Vector of samples from the confidence distribution of the between-study heterogeneity \eqn{\tau^2}. Obtained for method \code{MC}}.
#'   \item{fcd}{Marginalized confidence density function of \eqn{\mu}. Obtained for method \code{GAQ}.}
#' }
#'
#' @author David Kronthaler
#'
#' @details
#' This function performs a random-effects meta-analysis using Edgington's confidence distribution: study-specific confidence distributions, respectively, 
#' one-sided p-value functions for the alternative "greater" are constructed from normal pivots unter the normal random-effects model. These 
#' are then combined using Edgington's method, yielding a combined confidence distribution, respectively  
#' combined p-value function, of the average effect \eqn{\mu}. Then, the two approaches proceed as:
#' \describe{
#'      \item{\code{MC}:}{Generate samples from the approximate confidence distribution of the heterogeneity parameter \eqn{\tau^2},
#'      implied by the generalized heterogeneity statistic (Viechtbauer, 2006). For each sampled \eqn{\tau^{2*}},
#'      generate one \eqn{\mu^*} from Edgington's confidence distribution. Repeatedly sampling heterogeneity  \eqn{\tau^{2*}} and
#'      conditional sampling \eqn{\mu^*}, a confidence distribution of \eqn{\mu}, incorporating uncertainty about the
#'      heterogeneity parameter, is constructed, from which point estimates and confidence intervals are computed.}
#'      \item{\code{GAQ}:}{The confidence density of \eqn{\tau^2}, implied by the generalized heterogeneity statistic,
#'      is constructed by change of variables. Global adaptive quadrature integration is used to integrate Edgington's confidence density
#'      over the confidence density of \eqn{\tau^2}. This typically corresponds to the MC confidence distribution, except if very few
#'      studies are available. In this case, GAQ integration may produce slightly too wide confidence distributions and intervals.}
#' }
#' 
#' Additionally, study-specific weights can be incorporated, for example to weight studies by
#' the inverse of their standard errors or variances, or to downweight studies considered to 
#' be at higher risk of bias. When weights other than one are specified, Edgington's weighted
#' \eqn{p}-value function is used to construct the confidence distribution for the average effect.
#'
#' @references
#' Viechtbauer, W. (2007). *Confidence intervals for the amount of heterogeneity in meta‐analysi*s. Statistics in medicine, 26(1), 37-52. https://doi.org/10.1002/sim.2514
#'
#' Held, L., Hofmann, F., & Pawel, S. (2025). *A comparison of combined p-value functions for meta-analysis*. doi:10.1017/rsm.2025.26

#'
#' @examples
#' es <- c(0.17,  1.20,  1.10, -0.0019, -2.33)
#' se <- c(0.52, 0.93, 0.63, 0.3, 0.28)
#' remaeffect(es = es, se = se, method = "MC")
#'
#'
#'
remaeffect <- function(es,
                       se,
                       method = "MC",
                       w = rep(1, length(es)),
                       level.ci = 0.95,
                       B = 100000L,
                       seed = NULL,
                       mu0 = 0) {
  
  # validate input
  vd_remaeffect(es = es, se = se, w = w, method = method, level.ci = level.ci,
                B = B, mu0 = mu0)
  
  out <- switch(
    method,
    "MC" = reffMC(es = es, se = se, w = w, level.ci = level.ci, B = B,
                  seed = seed, mu0 = mu0),
    "GAQ" = reffAQ(es = es, se = se, w = w, level.ci = level.ci, mu0 = mu0)
  )
  
  out$method <- method
  out$level.ci <- level.ci
  out$k <- length(es)
  out$weights <- w
  out$mu0 <- mu0
  out$es  <- es
  out$se  <- se
  if (method == "MC") {
    out$B <- B
  }
  class(out) <- "remaeffect"
  
  out
}


# Estimate average effect (mu) using Monte Carlo algorithm
reffMC <- function(es, se, w, level.ci, B, seed, mu0) {
    
    # Reproducibility
    if (!is.null(seed)) {
      if (!is.wholepositivenumber(seed)) {
        warning("Seed must be a valid scalar integer.")
      }
      set.seed(seed)
    }
    
    # Sample heterogeneity (tau2*) from its confidence distribution
    ma <- run_metagen(es = es, se = se, mtau2 = "REML")
    s_tau2 <- samptau2(B = B, es = es, se = se, upper = ma$tau2+100*ma$se.tau2)
    
    # Sample average effect (mu*) from Edgingtons confidence distribution
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se, w = w)
    
    # Point estimate
    hmu <- base::mean(s_mu, na.rm = T)
    ci <- stats::quantile(x = s_mu, p = (1 + level.ci * c(-1,1)) / 2, na.rm = T)
    
    # ---------------------------
    # could include HCD intervals
    # ---------------------------
    
    # Two-sided p-value at mu = mu0
    p2s <- 2 * min(mean(s_mu <= mu0), mean(s_mu >= mu0))
    
    # Return
    invisible(list(
      estimate = hmu,
      CI = ci,
      pval = p2s,
      cd_mu = s_mu,
      cd_tau2 = s_tau2
    ))
  }

# Estimate average effect (mu) with global adaptive quadrature integration
reffAQ <- function(es, se, w, level.ci, mu0) {
  
  # Marginal confidence distribution function of mu 
  ma <- run_metagen(es = es, se = se, mtau2 = "REML")
  utau2 <- ma$tau2 + 100 * ma$se.tau2
  cdfmu <- reff(es = es, se = se, w = w, utau2 = utau2, grid_step = 0.01)
  
  # Confidence interval (by inversion of CDF)
  quant <- function(p) {
    stats::approx(cdfmu[, 2], cdfmu[, 1], xout = p, rule = 2)$y
  }
  ci <- quant(c((1 - level.ci) / 2, (1 + level.ci) / 2))
  
  # Confidence density function
  fcd <- function(mu) marCD(mu = mu, es = es, se = se, w = w,  utau2 = utau2)
  
  # Point estimate
  dx <- diff(cdfmu[, 1]) # equi-spaced distance
  df <- diff(cdfmu[, 2]) # diff in F() from i to i+1
  x_mid <- cdfmu[, 1][-length(cdfmu[, 1])] + dx / 2 # midpoints
  hmu <-  sum(x_mid * df) # Numerical approximation
  
  # P-value against H0: mu = 0
  p1sf <- stats::approxfun(cdfmu[, 1], cdfmu[, 2], yleft = 0, yright = 1)
  p1s <- p1sf(mu0)                           # one-sided p-value
  p2s <- ifelse(p1s <= 0.5, p1s*2, 2*(1-p1s))   # two-sided p-value
  
  # Return
  invisible(list(
    estimate = hmu,
    CI = ci,
    pval = p2s,
    fcd = fcd
  ))
}

#' @export
print.remaeffect <- function(x, ...) {
  cat("\nCD-Edgington Random-Effects Meta-Analysis\n\n")
  cat("Method:", x$method, "\n")
  weight_type <- "custom"
  if (isTRUE(all.equal(x$w, 1/x$se^2, tolerance = 1e-4, check.attributes = FALSE))) {
    weight_type <- "inverse squared standard errors (1/se^2)"
  } else if (isTRUE(all.equal(x$w, 1/x$se, tolerance = 1e-4, check.attributes = FALSE))) {
    weight_type <- "inverse standard errors (1/se)"
  } else if (all(abs(x$w - 1) < 1e-4)) {
    weight_type <- "unweighted"
  }
  base::cat("Weights:", weight_type, "\n")
  if (!is.null(x$B)) {
    cat("Number of Monte Carlo samples:", format(x$B, big.mark = ","), "\n")
  }
  cat("Number of studies:", x$k, "\n\n")
  cat("Average effect:", sprintf("%.3f", x$estimate), "\n")
  cat(paste0(x$level.ci * 100, "%"), "confidence interval from",
      paste(sprintf("%.3f", x$CI), collapse = " to "), "\n"
  )
  cat("Two-sided p-value against H0: mu =", x$mu0, "is", round(x$pval, 5), "\n")
  invisible(x)
}








