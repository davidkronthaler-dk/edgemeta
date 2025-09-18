#' Estimate the Average Effect in Random-Effects Meta-Analysis via Edgington's Confidence Distribution
#'
#' Confidence distributions for the average effect \eqn{\mu}, or equivalently, \eqn{p}-value functions, from individual studies
#' are combined using Edgington's approach. The resulting combined confidence distribution of \eqn{\mu}, 
#' which is conditional on the heterogeneity parameter \eqn{\tau^2}, is integrated over an approximate confidence distribution of the latter
#' to account for heterogeneity estimation uncertainty. Alternatively, the approach can be performed without such marginalization using a fixed
#' heterogeneity estimate.
#'
#' @param es Numeric vector of effect estimates from individual studies (length >= 2).
#' @param se Numeric vector of standard errors corresponding to each effect estimate (length >= 2).
#' @param method Either "MC" for Monte Carlo sampling algorithm, "GAQ" for global adaptive quadrature integration, or "NHEU" for the unadjusted (unmarginalized) approach (compare details).
#' @param level.ci Confidence level for the interval estimate (default is 0.95).
#' @param n_samples Number of Monte Carlo samples used to estimate the confidence distributions (default is 100,000). Relevant when selecting \code{method = "MC"}.
#' @param seed Optional integer to ensure reproducibility of the random sampling. Relevant when selecting \code{method = "MC"}.
#' @param mu0 Parameter value under the null-hypothesis against which two-sided p-value is computed.
#' @param method.tau2 Method to compute the heterogeneity estimate. Relevant when selecting \code{method = "NHEU"}.
#'
#' @return A list, including not all, but some of these, depending on the method:
#' \describe{
#'   \item{estimate}{Point estimate of the average effect \eqn{\mu}.}
#'   \item{CI}{Confidence interval for the average effect \eqn{\mu}.}
#'   \item{pval}{Two-sided p-value against H0: \eqn{\mu} = \eqn{\mu_0}.}
#'   \item{cd_mu}{Vector of samples from the confidence distribution of the average effect \eqn{\mu}. Obtained for method \code{MC}}
#'   \item{cd_tau2}{Vector of samples from the confidence distribution of the between-study heterogeneity \eqn{\tau^2}. Obtained for method \code{MC}}
#'   \item{fcd}{Marginalized confidence density function of \eqn{\mu}. Obtained for method \code{GAQ}.}
#' }
#'
#' @author David Kronthaler
#'
#' @details
#' This function performs a random-effects meta-analysis using Edgington's confidence distribution: study-specific confidence distributions, respectively, 
#' one-sided p-value functions for the alternative "greater" are constructed from normal pivots unter the normal random-effects model. These 
#' are then combined using Edgington's method, yielding a combined confidence distribution, respectively  
#' combined p-value function, of the average effect \eqn{\mu}. Then, the three approaches proceed as:
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
#'      \item{\code{NHEU}:}{Derives point estimates and intervals directly from Edgington's approach wihtout adjusting for uncertainty
#'      in the heterogeneity estimate (Held et al, 2025).}
#' }
#' 
#' The methods \code{MC} and \code{GAQ} always produces confidence intervals reflecting potential heterogeneity, and are applicable only 
#' under a random-effects framework.
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
                       method = c("MC", "GAQ", "NHEU"),
                       level.ci = 0.95,
                       n_samples = 100000L,
                       seed = NULL,
                       mu0 = 0,
                       method.tau2 = "REML") {
  
  # validate input
  vd_remaeffect(
    es = es,
    se = se,
    method = method,
    level.ci = level.ci,
    n_samples = n_samples,
    mu0 = mu0,
    method.tau2 = method.tau2
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
      mu0 = mu0,
      mtau2 = method.tau2
    )
  }
  
  invisible(r)
}


# Estimate 'mu' using Monte Carlo algorithm
reffMC <-
  function(es,
           se,
           level.ci,
           n_samples,
           seed,
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
    
    # -------------------------------------------------------------------
    # could include alternatively
    # pimeta::pima(y = es, se = se, method = "boot", B = n_samples)$rnd
    # to sample from the exact distribution of tau2 (also for 'PredDist')
    # -------------------------------------------------------------------
    
    # Sampling from Edgingtons confidence distribution
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se)
    
    # Point estimate
    hmu <- base::mean(s_mu, na.rm = T)
    ci <- stats::quantile(x = s_mu,
                          p = (1 + level.ci * c(-1, 1)) / 2,
                          na.rm = T)
    
    # ---------------------------
    # could include HCD intervals
    # ---------------------------
    
    # Two-sided p-value at mu = mu0
    p2s <- 2 * min(mean(s_mu <= mu0), mean(s_mu >= mu0))
    
    # Visual output
    cat("\nCD-Edgington Random-Effects Meta-Analysis\n\n")
    cat("Details:\nMonte Carlo Algorithm (stochastic, independent samples)\n")
    cat("Number of Monte Carlo samples:",
        format(n_samples, big.mark = ","),
        "\n\n")
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
reffAQ <- function(es,
                   se, 
                   level.ci, 
                   mu0) {
  
  # Upper integration bound of tau2
  ma <- run_metagen(es = es, se = se, mtau2 = "REML")
  utau2 <- ma$tau2 + 100 * ma$se.tau2
  
  # Marginal confidence distribution function of mu 
  cdfmu <- reff(es, se, utau2) 
  
  # Confidence interval (by inversion of CDF)
  quant <- function(p) {
    stats::approx(cdfmu[, 2], cdfmu[, 1], xout = p, rule = 2)$y
  }
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
  p1sf <- stats::approxfun(cdfmu[, 1], cdfmu[, 2], yleft = 0, yright = 1)
  p1s <- p1sf(mu0)    # one-sided p-value
  p2s <- p1tp2(p1s)   # two-sided p-value
  
  # Visual output
  cat("\nCD-Edgington Random-Effects Meta-Analysis\n\n")
  cat("Details: Global Adaptive Quadrature Integration\n\n")
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

# Estimate 'mu' without adjustment for heterogeneity estimation uncertainty
reffNHEU <- function(es,
                     se, 
                     level.ci,
                     mu0, 
                     mtau2) {
  
  tau2 <- run_metagen(es = es, se = se, mtau2 = mtau2)$tau2
  hmu <- opti_num(es = es, se = sqrt( se ^ 2 + tau2), 
                  point = TRUE, ci = TRUE, levelci = level.ci)
  ci <- c(hmu$cilower, hmu$ciupper)
  p1s <- pfctedge(h0 = mu0, es = es, se = se)
  p2s <- p1tp2(p1s)
  
  cat("\nRandom-Effects Meta-Analysis using Edgington's Method\n\n")
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












