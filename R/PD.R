#' Frequentist Predictive Distributions and Prediction Intervals for Meta-Analysis
#'
#' @param es Estimates from individual studies (numeric vector, length >= 2).
#' @param se Standard errors from individual studies (numeric vector, length >= 2).
#' @param method One of "FullCD", "SimplifiedCD", "FixedTau2". Check details for information.
#' @param level.pi Level of the prediction interval computed (numeric, between 0 and 1).
#' @param n_samples In case method is "FullCD" or "SimplifiedCD", this determines the number of samples generated for the computation of the predictive distribution. For method "FixedTau", this determines the number of samples generated from the deterministic predictive distribution function.
#' @param theta_new In case method is "FixedTau2", these are the points for which the predictive density is evaluated and returned.
#' @param method.tau2 In case method is "FixedTau2" or "SimplifiedCD", this determines the method of estimating tau2. Check 'help(meta)' for information.
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
           method.tau2 = "REML") {
    # validate input
    validate_input(
      es = es,
      se = se,
      method = method,
      lpi = level.pi,
      ns = n_samples,
      mtau2 = method.tau2
    )
    
    # computation
    if (method == "FixedTau2") {
      rt = pd_cd(
        es = es,
        se = se,
        mtau2 = method.tau2,
        lpi = level.pi,
        ns = n_samples
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
pd_cd <- function(es, se, mtau2 = "REML", lpi = 0.95, ns = 100000L){

  # Estimate tau2
  ma <- meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2)
  tau2 = ma$tau2

  # If estimated tau2 is zero, return the confidence interval for mu
  if (tau2 == 0L) {
    opt <- opti_num(es, sqrt(se ^ 2 + tau2))
    warning("Tau2 is estimated to be zero (no between-study heterogeneity).
The confidence interval and confidence density for the pooled effect are returned.")
    return(list(CI = c(opt$cilower, opt$ciupper),
                fCD = function(mu) CD_cpp(mu, es = es, se = se)))
  }
  
  # Generate samples of mu
  s_mu <- samplemusimple(n_samples = ns, tau2 = tau2, es = es, se = se)
  
  # Generate samples of theta_new
  s_tn <- base::suppressWarnings(stats::rnorm(n = ns, mean = s_mu, sd = sqrt(tau2)))
  
  # Combine
  s <- base::cbind(s_mu, s_tn)
  base::colnames(s) <- c("mu", "theta_new")
  
  # Prediction interval
  pi <- stats::quantile(x = s_tn, p = (1 + lpi * c(-1, 1)) / 2, na.rm = T)

  # Return
  return(list(PI = pi, samples = s))
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


