#' Frequentist Predictive Distributions and Prediction Intervals for Meta-Analysis
#'
#' @param es Estimates from individual studies (numeric vector, length >= 2).
#' @param se Standard errors from individual studies (numeric vector, length >= 2).
#' @param method One of "FullCD", "SimplifiedCD", "FixedTau2". Check details for information.
#' @param level.pi Level of the prediction interval computed (numeric, between 0 and 1).
#' @param n_samples In case method is "FullCD" or "SimplifiedCD", this determines the number of samples generated for the computation of the predictive distribution. For method "FixedTau", this determines the number of samples generated from the deterministic predictive distribution function.
#' @param theta_new In case method is "FixedTau2", these are the points for which the predictive density is evaluated and returned.
#' @param method.tau2 In case method is "FixedTau2" or "SimplifiedCD", this determines the method of estimating tau2. Check 'help(meta)' for information.
#' @param subdivisions Number of subdivisions used for integration (relevant for method "FixedTau2)"
#' @param ... Additional arguments handed to 'stats::integrate' (relevant for method "FixedTau2).
#'
#' @return Return a prediction interval.
#' @export
#'
#' @examples
#' es <- c(0.17,  1.20,  1.10, -0.0019, -2.33)
#' se <- c(0.52, 0.93, 0.63, 0.3, 0.28)
#' PredDist(es = es, se = se, method = "FullCD")
PredDist <-
  function(es,
           se,
           method = c("FullCD", "SimplifiedCD", "FixedTau2"),
           level.pi = 0.95,
           n_samples = 100000L,
           theta_new = NULL,
           method.tau2 = "REML",
           subdivisions = 300L,
           ...) {
    # validate input
    validate_input(
      es = es,
      se = se,
      method = method,
      lpi = level.pi,
      ns = n_samples,
      theta_new = theta_new,
      mtau2 = method.tau2,
      subdivisions = subdivisions
    )
    
    # computation
    if (method == "FixedTau2") {
      rt = pd_cd(
        es = es,
        se = se,
        mtau2 = method.tau2,
        lpi = level.pi,
        theta_new = theta_new,
        subdivisions = subdivisions,
        ns = n_samples,
        ...
      )
    } else if (method == "SimplifiedCD" | method == "FullCD") {
      rt = pd_cd_tau2(
        es = es,
        se = se,
        lpi = level.pi,
        method = method,
        ns = n_samples,
        mtau2 = method.tau2
      )
    }
    
    # Return
    return(rt)
    
  }

## PD and PI based on confidence density, fixed tau2
##------------------------------------------------------------------------------
pd_cd <- function(es, se, mtau2 = "REML", lpi = 0.95, theta_new = NULL,
                  ns = 100000L, subdivisions = 100L, ...){

  # Estimate tau2
  ma <- meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2)
  tau2 = ma$tau2

  # Adjust standard errors with tau2
  sea <- base::sqrt(se ^ 2 + tau2)

  # If estimated tau2 is zero, return the confidence interval for mu
  if (tau2 == 0L) {
    opt <- opti_num(es, sea)
    warning("Tau2 is estimated to be zero (no between-study heterogeneity).
The confidence interval and confidence density for the pooled effect are returned.")
    return(list(CI = c(opt$cilower, opt$ciupper),
                fCD = function(mu) CD_cpp(mu, es = es, se = se)))
  }

  # Approximate the confidence density of mu
  lim <- ma$TE.random + c(-10, 10) * ma$seTE.random
  theta_grid <- base::seq(from = lim[1], to = lim[2], l = 500)
  cd_apprx <- stats::approxfun(x = theta_grid,
                               y = CD_cpp(h0 = theta_grid, es = es, se = sea),
                               rule = 2)

  # Joint density of theta_new and mu
  pd_intern <- base::Vectorize(function(tn_int) {
    result <- base::tryCatch({
      stats::integrate(f = function(mu) {
        return(stats::dnorm(tn_int, mu, sqrt(ma$tau2)) * cd_apprx(mu))
      }, subdivisions = subdivisions, lower = lim[1], upper = lim[2], ...)$value
    }, error = function(e) {
      return(NA)
    })
  }, vectorize.args = "tn_int")

  # Approximate the predictive density
  pd_improper <- stats::approxfun(x = theta_grid, y = pd_intern(theta_grid),
                                  yleft = 0, yright = 0)

  # Normalize the predictive density after approximation
  c <- stats::integrate(f = pd_improper, lower = lim[1], upper = lim[2],
                        subdivisions = subdivisions, ...)$value
  pd <- function(x) pd_improper(x)/c

  # Prediction interval
  pi <- pinum(density = pd, lpi = lpi, lower = lim[1], upper = lim[2])
  
  # Generate samples from the predictive distribution
  if (ns > 0) {
    s_mu <- samplemusimple(n_samples = ns, tau2 = tau2, es = es, se = se)
    s_tn <- base::suppressWarnings(stats::rnorm(n = ns, mean = s_mu, sd = sqrt(tau2)))
    s <- base::cbind(s_mu, s_tn)
    base::colnames(s) <- c("mu", "theta_new")
  } else {
    s <- NULL
  }

  # Return
  return(list(PI = pi, PD.theta_new = pd(theta_new), fPD = pd, samples = s))
}


## Predictive distribution generated using sampling algorithm
##------------------------------------------------------------------------------
pd_cd_tau2 <- function(es, se, lpi = 0.95, method = c("SimplifiedCD", "FullCD"),
                       ns = 100000L, mtau2 = "REML") {

  # Initial guess for tau2
  ma <- meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2)

  # Samples of tau2
  s_tau2 <- samptau2(ns = ns, es = es, se = se, 
                     upper = ma$tau2 + 100 * ma$se.tau2)
  # Samples of mu 
  if (method == "SimplifiedCD") {
    s_mu <- samplemusimple(n_samples = ns, tau2 = ma$tau2, es = es, se = se)
  } else if (method == "FullCD") {
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se)
  }
  
  # Samples of tn
  s_tn <- base::suppressWarnings(stats::rnorm(n = ns, mean = s_mu, sd = sqrt(s_tau2)))

  # Combine 
  s <- base::cbind(s_tau2, s_mu, s_tn)
  base::colnames(s) <- c("tau2", "mu", "theta_new")

  # Prediction interval
  PI <- stats::quantile(x = s_tn, p = (1 + lpi * c(-1, 1)) / 2, na.rm = T)
  
  # Return
  return(list(PI = PI, samples = s))

}

## Numerically evaluate density to obtain prediction intervals
##------------------------------------------------------------------------------

pinum <- function(density, lpi = 0.95, lower, upper) {

  # quantiles
  probs <- c((1 - lpi) / 2, (1 + lpi) / 2)

  # Function to find bounds
  fb <- function(p) {
    integrand <- function(u) {
      base::abs(stats::integrate(f = density, lower = lower, upper = u)$value - p)
    }

    # Integration with fallback option
    b <- tryCatch(
      stats::optimize(f = integrand, lower = lower, upper = upper)$minimum,
      error = function(e) {
        sfb <- function(u) {
          base::abs(stats::integrate(density, lower = lower, upper = u,
                                     subdivisions = 1000)$value - p)
        }
        stats::optimize(f = sfb, lower = lower * 0.9, upper = upper * 0.9)$minimum
      }
    )
    return(b)
  }

  # Compute lower and upper bounds
  lwr <- fb(probs[1])
  upr <- fb(probs[2])

  # Return
  invisible(c(lwr, upr))
}












