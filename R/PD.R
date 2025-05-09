#' Frequentist Predictive Distributions and Prediction Intervals for Meta-Analysis
#'
#' @param es Estimates from individual studies (numeric vector, length >= 2).
#' @param se Standard errors from individual studies (numeric vector, length >= 2).
#' @param method One of "FullCD", "SimplifiedCD", "FixedTau2". Check details for information.
#' @param level.pi Level of the prediction interval computed (numeric, between 0 and 1).
#' @param n_samples In case method is "FullCD", this determines the number of samples generated for the computation of the predictive distribution.
#' @param theta_new In case method is "FixedTau2" or "SimplifiedCD", these are the points for which the predictive density is evaluated and returned.
#' @param method.tau2 In case method is "FixedTau2", this determines the method of estimating tau2. Check 'help(meta)' for information.
#' @param subdivisions Number of subdivisions used for integration.
#' @param ... Additional arguments handed to 'stats::integrate'.
#'
#' @return Return a prediction interval. In case method is "FullCD" returns a matrix containing samples of theta.new, mu and tau2. In case method is "SimplifiedCD" or "FixedTau2" returns additionally the predictive density in form of a function.
#' @export
#'
#' @examples
#' es <- c(0.17,  1.20,  1.10, -0.0019, -2.33)
#' se <- c(0.52, 0.93, 0.63, 0.3, 0.28)
#' PredDist(es = es, se = se, method = "FullCD")
PredDist <- function(es, se, method = c("FullCD", "SimplifiedCD", "FixedTau2"),
                     level.pi = 0.95, n_samples = 100000L, theta_new = NULL,
                     method.tau2 = "REML", subdivisions = 300L, ...) {

  # validate input
  validate_input(es = es, se = se, method = method, lpi = level.pi,
                 ns = n_samples, theta_new = theta_new,
                 mtau2 = method.tau2, subdivisions = subdivisions)

  # computation
  if (method == "FixedTau2") {
    rt = pd_cd(es = es, se = se, mtau2 = method.tau2, lpi = level.pi,
                theta_new = theta_new, subdivisions = subdivisions, ...)
  } else if (method == "SimplifiedCD") {
    rt = pd_cd_tau2(es = es, se = se, lpi = level.pi, theta_new = theta_new,
                    method = method, subdivisions = subdivisions, ...)
  } else if (method == "FullCD") {
    rt = pd_cd_tau2(es = es, se = se, lpi = level.pi, method = method,
                    ns = n_samples)
  }

  # Return
  return(rt)

}


## PD and PI based on confidence density, fixed tau2
##------------------------------------------------------------------------------
pd_cd <- function(es, se, mtau2 = "REML", lpi = 0.95,  theta_new = NULL,
                  ns = 100000L, subdivisions = 100L, ...){

  # Compute tau2 with the 'meta' package
  ma <- meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2)
  tau2 = ma$tau2

  # Adjust standard errors with tau2
  sea <- base::sqrt(se^2 + tau2)

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
  
  # Generate samples from the Predictive distribution
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


## Function to compute predictive distribution and PI for CD and tau2 density
##------------------------------------------------------------------------------
pd_cd_tau2 <- function(es, se, lpi = 0.95, theta_new = NULL,
                       method = c("SimplifiedCD", "FullCD"),
                       ns = 100000L, subdivisions = 300L) {

  # Compute a grid for mu and tau2
  ma <- meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = "REML")
  lim <- ma$TE.random + c(-12, 12) * ma$seTE.random
  gmu <- base::seq(from = lim[1], to = lim[2], length.out = 100)
  gtau2 <- base::seq(from = 0, to = ma$tau2 + 10 * ma$se.tau2, length.out = 100)

  # Prediction interval using Adaptive Quadrature (AQ) integration
  if (method == "SimplifiedCD") {

    # Compute the predictive distribution
    pd <- aq_pd(max.gtau2 = gtau2[length(gtau2)], gmu = gmu, lim = lim,
                es = es, se = se, es.tau2 = ma$tau2,
                subdivisions = subdivisions)

    # Prediction interval
    pi <- pinum(density = pd, lpi = lpi, lower = lim[1], upper = lim[2])
    
    # Generate samples from the Predictive distribution
    if (ns > 0) {
      s_tau2 <- samptau2(ns = ns, es = es, se = se,
                         upper = ma$tau2 + 100 * ma$se.tau2)
      s_mu <- samplemusimple(n_samples = ns, tau2 = ma$tau2, es = es, se = se)
      s_tn <- base::suppressWarnings(stats::rnorm(n = ns, mean = s_mu, sd = sqrt(s_tau2)))
      s <- base::cbind(s_tau2, s_mu, s_tn)
      base::colnames(s) <- c("tau2", "mu", "theta_new")
    } else {
      s <- NULL
    }

    # Return
    return(list(PI = pi, PD.theta_new = pd(theta_new), fPD = pd, samples = s))

  # Prediction interval using sampling (MC and MCMC)
  } else if (method == "FullCD") {

    # Generate samples of tau2, mu, theta_new
    s_tau2 <- samptau2(ns = ns, es = es, se = se,
                              upper = ma$tau2 + 100 * ma$se.tau2)
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se)
    s_tn <- base::suppressWarnings(stats::rnorm(n = ns, mean = s_mu, sd = sqrt(s_tau2)))
    s <- base::cbind(s_tau2, s_mu, s_tn)
    base::colnames(s) <- c("tau2", "mu", "theta_new")
    
    # Prediction interval
    PI <- stats::quantile(x = s_tn, p = (1 + lpi * c(-1, 1))/2, na.rm = T)

    # Return
    return(list(PI = PI, samples = s))

  }
}

## Adaptive Quadrature: Predictive distribution
##------------------------------------------------------------------------------
aq_pd <- function(max.gtau2, gmu, es, se, es.tau2, lim,
                  subdivisions) {

  # Joint density of theta_new, mu, tau2
  joint3 <- function(tn, mu, tau2) {
    return(stats::dnorm(tn, mu, base::sqrt(tau2)) *  # P(theta_new | mu, tau2)
             CD_cpp(mu, es, base::sqrt(se^2 + es.tau2)) * # Marginal P(mu)
        ftau2(es, se, tau2)) # Marginal P(tau2)
  }

  # Marginalize over f(tau2)
  joint2 <- base::Vectorize(function(tn, mu) {
    tryCatch({
      stats::integrate(function(tau2_intern) joint3(tn, mu, tau2_intern),
                subdivisions = subdivisions,
                lower = .Machine$double.eps,
                upper = max.gtau2)$value
    }, error = function(e) { return(NA) })
  }, vectorize.args = c("tn", "mu"))

  # Marginalize over f(mu)
  marginal_tn_improper <- base::Vectorize(function(tn) {
    base::tryCatch({
      stats::integrate(function(mu) joint2(tn, mu),
                subdivisions = subdivisions,
                lower = lim[1],
                upper = lim[2])$value
    }, error = function(e) { return(NA) })
  }, vectorize.args = "tn")

  # Approximate the predictive density
  marginal_tn_approx <- stats::approxfun(gmu, marginal_tn_improper(gmu),
                                  yleft = 0, yright = 0)

  # Ensure the pred. density integrates to 1
  c <- base::tryCatch({
    stats::integrate(f = marginal_tn_approx, lower = lim[1], upper = lim[2],
                     subdivisions = subdivisions)$value},
                    error = function(e) { return(1) })

  # Normalized predictive density
  marginal_tn <- function(theta_new) { marginal_tn_approx(theta_new)/c }

  # Return
  return(marginal_tn)
}


## Numerically evaluate density to obtain prediction (or confidence) intervals
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












