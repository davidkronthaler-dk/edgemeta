#' Frequentist Predictive Distributions and Intervals for Random-Effects Meta-Analysis using Edgington's Method
#' 
#' Computes frequentist predictive distributions and prediction intervals for future effects
#' \eqn{\theta_{new}} in random-effects meta-analysis. Predictive distributions are adjusted for estimation uncertainty
#' through Edgington's confidence distributions of the average effect (\eqn{\mu}) and through an approximate confidence distribution of 
#' the between-study heterogeneity (\eqn{\tau^2}).
#' Computation is performed using a Monte Carlo sampling algorithm generating independent samples.
#'
#' @param es Numeric vector of effect estimates from individual studies (length >= 2).
#' @param se Numeric vector of standard errors corresponding to each effect estimate (length >= 2).
#' @param method Either "FullCD" (recommended for practical application), "SimplifiedCD" or "FixedTau2". Check details for information. 
#' @param level.pi Coverage level of the equi-tailed prediction interval (numeric, between 0 and 1).
#' @param n_samples Number of Monte Carlo samples used to construct the predictive distribution (default is 100,000).
#' @param method.tau2 In case method is "FixedTau2" or "SimplifiedCD", this determines the method of estimating the between-study heterogeneity. 
#' Check 'help(meta)' for information on estimation methods (default is "REML").
#' @param seed Optional integer to ensure reproducibility of the random sampling.
#' @return A list containing:
#' \describe{
#'   \item{PI}{Prediction interval for a future effect \eqn{\theta_{new}}.}
#'   \item{samples}{A matrix containing samples \eqn{\theta_{new}^*} from the predictive distribution ('theta_new'), samples \eqn{\mu^*} from the confidence distribution of the average effect ('mu') and samples \eqn{\tau^{2*}} from the confidence distribution of the between-study heterogeneity ('tau2').
#'   The latter are not included in case \code{method = "FixedTau2"} is selected.}
#' }

#' @author David Kronthaler
#' 
#' @details
#' Predictive distributions are constructed by integrating the assumed normal effect distribution \eqn{\theta_{new}} ~ N(\eqn{\mu}, \eqn{\tau^2}) over Edgington's confidence distribution
#' to account for uncertainty in \eqn{\mu} and over an approximate confidence distribution of the heterogeneity parameter to account for uncertainty in \eqn{\tau^2}. 
#' The approximate confidence distribution of \eqn{\tau^2} is implied by the generalized heterogeneity statistic (Viechtbauer, 2006).
#' 
#'  The function supports three methods:
#' 
#' - **"FullCD"**: Samples \eqn{\tau^{2*}} are drawn from the approximate confidence distribution of \eqn{\tau^2}. For each drawn
#'  \eqn{\tau^{2*} }, a corresponding \eqn{\mu^*} is conditionally drawn from Edgington's confidence distribution (, which is conditional
#'  on the heterogeneity parameter). For each pair (\eqn{\mu^*}, \eqn{\tau^{2*}}), a future effect \eqn{\theta_{new}^*} is drawn from
#'  a N(\eqn{\mu^*}, \eqn{\tau^{2*}}). This method fully accounts for uncertainty in both parameters and the dependence of Edgington's
#'  confidence distribution on the heterogeneity parameter.
#' 
#' - **"SimplifiedCD"**: This method generates samples \eqn{\tau^{2*}} from the confidence distribution of 
#'   \eqn{\tau^2}, but generates samples \eqn{\mu^*} using a simplified
#'   approach. Specifically, samples \eqn{\mu^*} are drawn conditional on a fixed \eqn{\hat{\tau}^2}.
#'   For each pair (\eqn{\mu^*}, \eqn{\tau^{2*}}), a future effect \eqn{\theta_{new}^*} is drawn from
#'   a N(\eqn{\mu^*}, \eqn{\tau^{2*}}).
#'   
#' 
#' - **"FixedTau2"**: This method used a fixed estimate \eqn{\hat{\tau}^2}, and generates samples \eqn{\mu^*}
#'   conditional on this fixed estimate. Samples \eqn{\theta_{new}^*} are generated from a N(\eqn{\mu^*}, \eqn{\hat{\tau}^2}).
#' 
#' 
#' The empirical distribution of the sampled \eqn{\theta_{new}^*} values serves as the estimated predictive 
#' distribution of future effects. This distribution can be used to compute prediction intervals, 
#' summarize predictive uncertainty, or generate visualizations.
#' 
#' @references
#' Viechtbauer, W. (2007). *Confidence intervals for the amount of heterogeneity in meta‐analysis*. Statistics in medicine, 26(1), 37-52. https://doi.org/10.1002/sim.2514
#' 
#' Held, L., Hofmann, F., & Pawel, S. (2025). *A comparison of combined p-value functions for meta-analysis*. doi:10.1017/rsm.2025.26

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
           method.tau2 = "REML",
           seed = NULL) {
    
    # validate input
    vd_PredDist(es = es,
                se = se,
                method = method,
                lpi = level.pi,
                ns = n_samples,
                mtau2 = method.tau2)
    
    # reproducibility, if required
    if (!is.null(seed)) set.seed(seed)
    
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
    
    class(rt) <- "metaprediction"
    attr(rt, "method") <- method
    attr(rt, "n_samples") <- n_samples
    attr(rt, "level_pi") <- level.pi
    attr(rt, "k") <- length(es)
    print(rt)
    
    invisible(rt)
  }

# Predictive distribution (fixed tau2)
pd_cd <- function(es,
                  se, 
                  mtau2, 
                  lpi, 
                  ns){

  # Estimate tau2
  ma <- run_metagen(es = es, se = se, mtau2 = mtau2)
  tau2 <- ma$tau2

  # If estimated tau2 is zero, return the confidence interval for mu
  if (tau2 == 0L) {
    opt <- opti_num(es, se)
    warning("Tau2 is estimated to be zero (no between-study heterogeneity). 
             The pooled effect estimate and confidence interval are returned.")
    return(list(estimate = opt$point, CI = c(opt$cilower, opt$ciupper)))
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

# Predictive distribution, adjusted for uncertainty in tau2
pd_cd_tau2 <- function(es, 
                       se, 
                       lpi, 
                       method,
                       ns,
                       mtau2) {

  # Initial tau2
  ma <- run_metagen(es = es, se = se, mtau2 = mtau2)

  # Samples of tau2
  s_tau2 <- samptau2(ns = ns, es = es, se = se, 
                     upper = ma$tau2 + 100 * ma$se.tau2)
  # Samples of mu 
  if (method == "SimplifiedCD") {
    s_mu <- samplemusimple(n_samples = ns, tau2 = ma$tau2, es = es, se = se)
  } else if (method == "FullCD") {
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se)
  }
  
  # Samples of theta_new
  s_tn <- base::suppressWarnings(stats::rnorm(n = ns, mean = s_mu, sd = sqrt(s_tau2)))

  # Combine 
  s <- base::cbind(s_tau2, s_mu, s_tn)
  base::colnames(s) <- c("tau2", "mu", "theta_new")

  # Prediction interval
  PI <- stats::quantile(x = s_tn, p = (1 + lpi * c(-1, 1)) / 2, na.rm = T)
  
  # Return
  return(list(PI = PI, samples = s))
}


# Estimation of between-study heterogeneity
run_metagen <- function(es, 
                        se,
                        mtau2) {
  
  r <- tryCatch( # default settings
    meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2),
    error = function(e) NULL
  )
  
  if (is.null(r)) { # increase maxiter
    r <- tryCatch( 
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(maxiter = 10000)),
      error = function(e) NULL
    )
  }
  
  if (is.null(r)) { # decrease step size
    r <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(stepadj = 0.5)),
      error = function(e) NULL
    )
  }
  
  if (is.null(r)) { # increase maxiter and decrease step size
    r <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(maxiter = 10000, stepadj = 0.5)),
      error = function(e) NULL
    )
  }
  
  if (is.null(r)) { # increase maxiter and decrease step size
    r <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(maxiter = 100000, stepadj = 0.25)),
      error = function(e) NULL
    )
  }
  
  return(r)
}