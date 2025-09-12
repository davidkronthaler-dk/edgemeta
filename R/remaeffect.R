#' Estimate the Average Effect in Random-Effects Meta-Analysis via Edgington's Confidence Distribution
#'
#' Estimates the average effect in random-effects meta-analysis using Edgingtons Confidence Distribution.
#' Confidence distributions for the average effect \eqn{\mu}, or equivalently, \eqn{p}-value functions, from individual studies
#' are combined using Edgington's approach. The resulting combined confidence distribution of \eqn{\mu}, 
#' which is conditional on the heterogeneity parameter \eqn{\tau^2}, is marginalized over a confidence distribution of the latter
#' to account for estimation uncertainty. Alternatively, the approach can be performed without such marginalization using a fixed
#' heterogeneity estimate.
#'
#' @param es Numeric vector of effect estimates from individual studies (length >= 2).
#' @param se Numeric vector of standard errors corresponding to each effect estimate (length >= 2).
#' @param method Either "MC" for Monte Carlo sampling algorithm, "GAQ" for global adaptive quadrature integration, or "NHEU" for the unadjusted (unmarginalized) approach.
#' @param level.ci Confidence level for the interval estimate (default is 0.95).
#' @param n_samples Number of Monte Carlo samples used to estimate the confidence distributions (default is 100,000).
#' @param seed Optional integer to ensure reproducibility of the random sampling, relevant in case method 'MC' is chosen.
#' @param mu0 Parameter value under the null-hypothesis against which two-sided p-value is computed.
#'
#' @return A list, including not all, but some of these, depending on the method:
#' \describe{
#'   \item{estimate}{Point estimate of the average effect \eqn{\mu}.}
#'   \item{CI}{Confidence interval for the average effect \eqn{\mu}.}
#'   \item{pval}{Two-sided p-value against H0: \eqn{\mu} = \eqn{\mu_0}.}
#'   \item{cd_mu}{Vector of samples from the confidence distribution of the average effect \eqn{\mu}.}
#'   \item{cd_tau2}{Vector of samples from the confidence distribution of the between-study heterogeneity \eqn{\tau^2}.}
#'   \item{fcd}{Obtained for method 'GAQ': Marginalized confidence density function of \eqn{\mu}.}
#' }
#'
#' @author David Kronthaler
#'
#' @details
#' This function performs a random-effects meta-analysis using confidence distributions.
#' It generates samples from the confidence distribution of the between-study heterogeneity \eqn{\tau^2},
#' derived from the generalized heterogeneity statistic (Viechtbauer, 2006). For each sampled value of \eqn{\tau^2},
#' the function then computes a corresponding confidence distribution for the average
#' effect \eqn{\mu}, induced by the Edgington combined p-value function (Held et al., 2025), and generates a sample from
#' this distribution.
#'
#' By repeatedly sampling \eqn{\tau^2} and recalculating the confidence distribution of \eqn{\mu},
#' the method incorporates the uncertainty in heterogeneity estimation into the final
#' average effect estimate. The function returns the point estimate, its confidence
#' interval, and the full sample of values drawn from the confidence distributions.
#'
#' Since the method accounts for the uncertainty in \eqn{\tau^2}, it always returns a
#' confidence interval that adjusts for potential between-study heterogeneity,
#' that is, the framework of a fixed-effects meta-analysis is not applicable.
#'
#' @references
#' Viechtbauer, W. (2007). *Confidence intervals for the amount of heterogeneity in meta‐analysi*s. Statistics in medicine, 26(1), 37-52. https://doi.org/10.1002/sim.2514
#'
#' Held, L., Hofmann, F., & Pawel, S. (2025). *A comparison of combined p-value functions for meta-analysis*. doi:10.1017/rsm.2025.26

#'
#' @examples
#' es <- c(0.17,  1.20,  1.10, -0.0019, -2.33)
#' se <- c(0.52, 0.93, 0.63, 0.3, 0.28)
#' remaeffect(es = es, se = se)
#'
#'
#'
remaeffect <- function(es,
                       se,
                       method = c("MC", "GAQ", "NHEU"),
                       level.ci = 0.95,
                       n_samples = 100000L,
                       seed = NULL,
                       mu0 = 0) {
  method <- match.arg(method)
  
  validate_input(
    es = es,
    se = se,
    method = "estimate",
    lpi = level.ci,
    ns = n_samples,
    mtau2 = "estimate"
  )
  
  if (method == "MC") {
    r <- reffMC(
      es = es,
      se = se,
      level.ci = level.ci,
      n_samples = n_samples,
      seed = seed,
      mu0 = mu0
    )
  } else if (method == "GAQ") {
    r <- reffAQ(
      es = es,
      se = se,
      level.ci = level.ci,
      mu0 = mu0
    )
  } else if (method == "NHEU") {
    r <- reffNHEU(
      es = es,
      se = se,
      level.ci = level.ci,
      mu0 = mu0
    )
  }
  
  invisible(r)
}


# Estimate 'mu' using Monte Carlo algorithm
reffMC <-
  function(es,
           se,
           level.ci = 0.95,
           n_samples = 100000L,
           seed = NULL,
           mu0) {
    
    # Reproducibility under MC
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    # Initial tau2 for definition of grid
    ma <- run_metagen(es = es, se = se, mtau2 = "REML")
    
    # Sampling tau2 from its confidence distribution
    s_tau2 <- samptau2(
      ns = n_samples,
      es = es,
      se = se,
      upper = ma$tau2 + 100 * ma$se.tau2
    )
    
    # Sampling mu from Edgingtons confidence distribution
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se)
    
    # Point estimate
    hmu <- base::mean(s_mu, na.rm = T)
    ci <- stats::quantile(x = s_mu,
                          p = (1 + level.ci * c(-1, 1)) / 2,
                          na.rm = T)
    
    # P-value at mu = mu0
    p1s <- mean(s_mu < mu0) # one-sided
    p2s <- p1tp2(p1s)       # two-sided
    
    # Visual output
    cat("\nRandom-Effects Meta-Analysis using Confidence Distributions\n\n")
    cat("Number of studies:", length(es), "\n")
    cat("Number of Monte Carlo samples:",
        format(n_samples, big.mark = ","),
        "\n\n")
    cat("Average effect:", sprintf("%.3f", hmu), "\n")
    cat(
      paste0(level.ci * 100, "%"),
      "Confidence interval from",
      paste(sprintf("%.3f", ci), collapse = " to "),
      "\n"
    )
    cat("Two-sided p-value against H0: mu =",
        mu0,
        "is",
        round(p2s, 5),
        "\n\n")
    cat("Summary of confidence distribution of the average effect:\n")
    probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    qs <- stats::quantile(s_mu, probs, na.rm = TRUE)
    m <- mean(s_mu, na.rm = TRUE)
    stats <- c(qs[1:2], Median = qs[3], Mean = m, qs[4:5])
    names(stats) <- c(
      sprintf("%.1f%%", probs[1] * 100),
      sprintf("%.1f%%", probs[2] * 100),
      "Median",
      "Mean",
      sprintf("%.1f%%", probs[4] * 100),
      sprintf("%.1f%%", probs[5] * 100)
    )
    print(as.data.frame(t(stats)), row.names = FALSE)
    
    # Return
    invisible(list(
      estimate = hmu,
      CI = ci,
      pval = p2s,
      cd_mu = s_mu,
      cd_tau2 = s_tau2
    ))
  }

# Estimate 'mu' using deterministic global adaptive quadrature integration
reffAQ <- function(es, se, level.ci = 0.95, mu0) {
  
  # Upper integration bound of tau2
  ma <- meta::metagen(es, se, method.tau = "REML")
  utau2 <- ma$tau2 + 100 * ma$se.tau2
  
  # Marginal confidence distribution function of mu
  cdfmu <- reff(es, se, utau2) 
  
  # Confidence interval (by inversion of CDF)
  quant <- function(p)
    stats::approx(cdfmu[, 2], cdfmu[, 1], xout = p, rule = 2)$y
  ci <- quant(c((1 - level.ci) / 2, (1 + level.ci) / 2))
  
  # Confidence density function
  fcd <- function(mu)
    marCD(
      mu = mu,
      es = es,
      se = se,
      utau2 = utau2
    )
  
  # Point estimate
  dx <- diff(cdfmu[, 1]) # equi-spaced distance
  df <- diff(cdfmu[, 2]) # diff in F() from i to i+1
  x_mid <- cdfmu[, 1][-length(cdfmu[, 1])] + dx / 2 # midpoints
  hmu <-  sum(x_mid * df) # Numerical approximation
  
  # P-value against H0: mu = 0
  p1sf <- stats::approxfun(cdfmu[, 1], cdfmu[, 2])
  p1s <- p1sf(mu0)    # one-sided p-value
  p2s <- p1tp2(p1s)   # two-sided p-value
  
  # Visual output
  cat("\nRandom-Effects Meta-Analysis using Confidence Distributions\n\n")
  cat("Number of studies:", length(es), "\n")
  cat("Average effect:", sprintf("%.3f", hmu), "\n")
  cat(
    paste0(level.ci * 100, "%"),
    "Confidence interval from",
    paste(sprintf("%.3f", ci), collapse = " to "),
    "\n"
  )
  cat("Two-sided p-value against H0: mu =",
      mu0,
      "is",
      round(p2s, 5),
      "\n\n")
  
  # Return
  invisible(list(
    estimate = hmu,
    CI = ci,
    pval = p2s,
    fcd = fcd
  ))
}

# Estimate 'mu' without adjustment for heterogenetiy estimation uncertainty
reffNHEU <- function(es, se, level.ci, mu0) {
  hmu <- opti_num(es = es, se = se, point = TRUE, ci = TRUE, levelci = level.ci)
  ci <- c(hmu$cilower, hmu$ciupper)
  p1s <- pfctedge(h0 = mu0, es = es, se = se)
  p2s <- p1tp2(p1s)
  
  # Visual output
  cat("\nRandom-Effects Meta-Analysis using Edgington's Confidence Distribution\n\n")
  cat("Details: Heterogenetiy estimation uncertainty is NOT accounted for.\n")
  cat("Consider method = 'MC' or 'GAQ' to incorporate this.\n\n")
  cat("Number of studies:", length(es), "\n")
  cat("Average effect:", sprintf("%.3f", hmu$point), "\n")
  cat(
    paste0(level.ci * 100, "%"),
    "Confidence interval from",
    paste(sprintf("%.3f", ci), collapse = " to "),
    "\n"
  )
  cat("Two-sided p-value against H0: mu =",
      mu0,
      "is",
      round(p2s, 5),
      "\n\n")
  
  # Return
  invisible(list(
    estimate = hmu$point,
    CI = ci,
    pval = p2s
  ))
}


# Transform one-sided to two-sided p-value
p1tp2 <- function(p1) {
  if (p1 <= 0.5) {
    return(p1 * 2)
  } else if (p1 > 0.5) {
    return(2 * (1 - p1))
  }
}












