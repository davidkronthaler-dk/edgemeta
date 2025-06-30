#' Frequentist Predictive Distributions and Intervals for Random-Effects Meta-Analysis
#' 
#' This function computes frequentist predictive distributions and prediction intervals for a future effect \eqn{\theta_{new}} in random-effects meta-analysis. Predictive distributions are constructed from confidence distributions of the average effect (\eqn{\mu}) and the between-study heterogeneity (\eqn{\tau^2}), using a Monte Carlo sampling algorithm. It supports three methods: "FullCD" (recommended), "SimplifiedCD", and "FixedTau2".
#'
#' @param es Numeric vector of effect estimates from individual studies (length >= 2).
#' @param se Numeric vector of standard errors corresponding to each effect estimate (length >= 2).
#' @param method Either "FullCD" (recommended for practical application), "SimplifiedCD" or "FixedTau2". Check details for information. 
#' @param level.pi Coverage level of the prediction interval computed (numeric, between 0 and 1).
#' @param n_samples Number of Monte Carlo samples used to construct the predictive distribution (default is 100,000).
#' @param method.tau2 In case method is "FixedTau2" or "SimplifiedCD", this determines the method of estimating the between-study heterogeneity. Check 'help(meta)' for information on estimation methods (default is "REML").
#' @param seed Optional integer to ensure reproducibility of the random sampling.
#' @return A list containing:
#' \describe{
#'   \item{PI}{Prediction interval for a future effect \eqn{\theta_{new}}.}
#'   \item{samples}{A matrix containing samples \eqn{\theta_{new}} from the predictive distribution ('theta_new'), samples \eqn{\mu} from the confidence distribution of the average effect ('mu') and samples \eqn{\tau^2} from the confidence distribution of the between-study heterogeneity ('tau2').}
#' }

#' @author David Kronthaler
#' 
#' @details
#' Predictive distributions are constructed from the confidence distribution of the average effect 
#' \eqn{\mu} and from the confidence distribution of heterogeneity parameter \eqn{\tau^2}. The confidence distribution of the 
#' average effect \eqn{\mu} is obtained from the Edgington combined \eqn{p}-value function (Held et al., 2025).
#' The confidence distribution of the between-study heterogeneity parameter \eqn{\tau^2} is 
#' derived from the generalized heterogeneity statistic (Viechtbauer, 2006).
#' 
#'  The function supports three methods:
#' 
#' - **"FullCD"**: This method generates samples from the confidence distribution of 
#'   \eqn{\tau^2}, and for each sampled \eqn{\tau^2}, it computes the corresponding confidence 
#'   distribution of \eqn{\mu}. It then generates a sample of the future effect \eqn{\theta_{new}} 
#'   for each (\eqn{\tau^2}, \eqn{\mu}) pair. This is the most comprehensive method, as it fully 
#'   accounts for uncertainty in both parameters and the dependence of the Edgington combined \eqn{p}-value function on the heterogeneity parameter.
#' 
#' - **"SimplifiedCD"**: This method generates samples from the confidence distribution of 
#'   \eqn{\tau^2}, but computes the confidence distribution of \eqn{\mu} using a simplified 
#'   approach. Specifically, it uses a fixed \eqn{\tau^2} to compute the Edgington combined 
#'   \eqn{p}-value function, from which samples of \eqn{\mu} are drawn. It then generates a sample of
#'    the future effect \eqn{\theta_{new}} for each (\eqn{\tau^2}, \eqn{\mu}) pair
#' 
#' - **"FixedTau2"**: This method assumes a fixed value for \eqn{\tau^2}, and uses it to compute 
#'   the confidence distribution of \eqn{\mu}. It does 
#'   not account for uncertainty in the estimation of \eqn{\tau^2}.
#' 
#' 
#' The empirical distribution of the sampled \eqn{\theta_{new}} values serves as the estimated predictive 
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
    validate_input(
      es = es,
      se = se,
      method = method,
      lpi = level.pi,
      ns = n_samples,
      mtau2 = method.tau2
    )
    
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

## Predictive distribution, fixed tau2
##------------------------------------------------------------------------------
pd_cd <- function(es, se, mtau2 = "REML", lpi = 0.95, ns = 100000L){

  # Estimate tau2
  ma <- run_metagen(es = es, se = se, mtau2 = mtau2)
  tau2 <- ma$tau2

  # If estimated tau2 is zero, return the confidence interval for mu
  if (tau2 == 0L) {
    opt <- opti_num(es, sqrt(se ^ 2 + tau2))
    warning("Tau2 is estimated to be zero (no between-study heterogeneity).
The confidence interval and confidence density for the pooled effect are returned.")
    return(list(estimate = opt$point,
                CI = c(opt$cilower, opt$ciupper),
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


## Predictive distribution, adjusted for uncertainty in tau2
##------------------------------------------------------------------------------
pd_cd_tau2 <- function(es, se, lpi = 0.95, method = c("SimplifiedCD", "FullCD"),
                       ns = 100000L, mtau2 = "REML") {

  # Initial guess for tau2
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

## Estimation of between-study heterogeneity
## -----------------------------------------------------------------------------
run_metagen <- function(es, se, mtau2 = "REML") {
  # Attempt 1: default settings
  result <- tryCatch(
    meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2),
    error = function(e) NULL
  )
  
  # Attempt 2: increase max iterations
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(maxiter = 10000)),
      error = function(e) NULL
    )
  }
  
  # Attempt 3: decrease step size
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(stepadj = 0.5)),
      error = function(e) NULL
    )
  }
  
  # Attempt 4: increase both maxiter and decrease step size
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(maxiter = 10000, stepadj = 0.5)),
      error = function(e) NULL
    )
  }
  
  # Attempt 5: increase both maxiter and decrease step size
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = mtau2,
                    control = list(maxiter = 100000, stepadj = 0.25)),
      error = function(e) NULL
    )
  }
  
  return(result)
}