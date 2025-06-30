#' Estimate the Average Effect in Random-Effects Meta-Analysis via Confidence Distributions
#'
#' Estimates the average effect in random-effects meta-analysis by generating samples 
#' from the confidence distribution of the average effect \eqn{\mu}. The confidence
#' distribution of \eqn{\mu} is adjusted for the uncertainty in the estimation of the
#' between-study heterogeneity.
#'
#' @param es Numeric vector of effect estimates from individual studies (length >= 2).
#' @param se Numeric vector of standard errors corresponding to each effect estimate (length >= 2).
#' @param level.ci Confidence level for the interval estimate (default is 0.95).
#' @param n_samples Number of Monte Carlo samples used to estimate the confidence distributions (default is 100,000).
#' @param seed Optional integer to ensure reproducibility of the random sampling.
#'
#' @return A list containing:
#' \describe{
#'   \item{estimate}{Point estimate of the average effect \eqn{\mu}.}
#'   \item{CI}{Confidence interval for the average effect \eqn{\mu}.}
#'   \item{cd_mu}{Vector of samples from the confidence distribution of the average effect \eqn{\mu}.}
#'   \item{cd_tau2}{Vector of samples from the confidence distribution of the between-study heterogeneity \eqn{\tau^2}.}
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
remaeffect <-
  function(es,
           se,
           level.ci = 0.95,
           n_samples = 100000L,
           seed = NULL) {
    # reproducibility
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    # validate
    validate_input(
      es = es,
      se = se,
      method = "estimate",
      lpi = level.ci,
      ns = n_samples,
      mtau2 = "estimate"
    )
    
    # Initial guess for tau2
    ma <- run_metagen(es = es, se = se, mtau2 = "REML")
    
    # Samples of tau2
    s_tau2 <- samptau2(
      ns = n_samples,
      es = es,
      se = se,
      upper = ma$tau2 + 100 * ma$se.tau2
    )
    # Samples of mu
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se)
    
    # Point estimate and confidence interval
    hmu <- base::mean(s_mu, na.rm = T)
    ci <-
      stats::quantile(x = s_mu,
                      p = (1 + level.ci * c(-1, 1)) / 2,
                      na.rm = T)
    
    # Visual output
    cat("\nRandom-Effects Meta-Analysis using Confidence Distributions\n\n")
    cat("Number of studies:", length(es), "\n")
    cat("Number of Monte Carlo samples:", format(n_samples, big.mark = ","),"\n\n")
    cat("Average effect:", sprintf("%.3f", hmu), "\n")
    cat(paste0(level.ci * 100, "%"), "Confidence interval from",
        paste(sprintf("%.3f", ci), collapse = " to "), "\n\n")
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
      cd_mu = s_mu,
      cd_tau2 = s_tau2
    ))
  }