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
#' @param method Either "PCD-full" (recommended for practical application), "PCD-simplified" or "PCD-fixed". Check details for information. 
#' @param w Study-specific weights.
#' @param level.pi Coverage level of the equi-tailed prediction interval (numeric, between 0 and 1).
#' @param B Number of Monte Carlo samples used to construct the predictive distribution (default is 100,000).
#' @param method.tau2 In case method is "PCD-fixed" or "PCD-simplified", this determines the method of estimating between-study heterogeneity. 
#' Check 'help(meta)' for information on estimation methods (default is "REML").
#' @param seed Optional integer to ensure reproducibility of the random sampling.
#' @return A list containing:
#' \describe{
#'   \item{PI}{Prediction interval for a future effect \eqn{\theta_{new}}.}
#'   \item{samples}{A matrix containing samples \eqn{\theta_{new}^*} from the predictive distribution ('theta_new'), samples \eqn{\mu^*} from the confidence distribution of the average effect ('mu') and samples \eqn{\tau^{2*}} from the confidence distribution of the between-study heterogeneity ('tau2').
#'   The latter are not included in case \code{method = "PCD-fixed"} is selected.}
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
#' - **"PCD-full"**: Samples \eqn{\tau^{2*}} are drawn from the approximate confidence distribution of \eqn{\tau^2}. For each drawn
#'  \eqn{\tau^{2*} }, a corresponding \eqn{\mu^*} is conditionally drawn from Edgington's confidence distribution (, which is conditional
#'  on the heterogeneity parameter). For each pair (\eqn{\mu^*}, \eqn{\tau^{2*}}), a future effect \eqn{\theta_{new}^*} is drawn from
#'  a N(\eqn{\mu^*}, \eqn{\tau^{2*}}). This method fully accounts for uncertainty in both parameters and the dependence of Edgington's
#'  confidence distribution on the heterogeneity parameter.
#' 
#' - **"PCD-simplified"**: This method generates samples \eqn{\tau^{2*}} from the confidence distribution of 
#'   \eqn{\tau^2}, but generates samples \eqn{\mu^*} using a simplified
#'   approach. Specifically, samples \eqn{\mu^*} are drawn conditional on a fixed \eqn{\hat{\tau}^2}.
#'   For each pair (\eqn{\mu^*}, \eqn{\tau^{2*}}), a future effect \eqn{\theta_{new}^*} is drawn from
#'   a N(\eqn{\mu^*}, \eqn{\tau^{2*}}).
#'   
#' 
#' - **"PCD-fixed"**: This method used a fixed estimate \eqn{\hat{\tau}^2}, and generates samples \eqn{\mu^*}
#'   conditional on this fixed estimate. Samples \eqn{\theta_{new}^*} are generated from a N(\eqn{\mu^*}, \eqn{\hat{\tau}^2}).
#' 
#' 
#' The empirical distribution of the sampled \eqn{\theta_{new}^*} values serves as the estimated predictive 
#' distribution of future effects. This distribution can be used to compute prediction intervals, 
#' summarize predictive uncertainty, or generate visualizations.
#' 
#' Additionally, study-specific weights can be incorporated, for example to weight studies by
#' the inverse of their standard errors or variances, or to downweight studies considered to 
#' be at higher risk of bias. When weights other than one are specified, Edgington's weighted
#' \eqn{p}-value function is used to construct the confidence distribution for the average effect.
#' 
#' @references
#' Viechtbauer, W. (2007). *Confidence intervals for the amount of heterogeneity in meta‐analysis*. Statistics in medicine, 26(1), 37-52. https://doi.org/10.1002/sim.2514
#' 
#' Held, L., Hofmann, F., & Pawel, S. (2025). *A comparison of combined p-value functions for meta-analysis*. doi:10.1017/rsm.2025.26

#' @export
#'
#' @examples
#' # Effect estimates and standard errors
#' es <- c(0.69, -0.22, -0.52, -0.77,  1.38, -0.34, -0.09)
#' se <- c(1.15, 0.25, 0.15, 0.42, 0.93, 0.32, 0.58)
#' 
#' # Predictive distribution PCD-full
#' PredDist(es = es, se = se, method = "PCD-full")
#' 
#' # Predictive distribution PCD-full weighting by inverse squared standard errors
#' PredDist(es = es, se = se, method = "PCD-full", w = 1/se^2)
PredDist <-
  function(es,
           se,
           method = "PCD-full",
           w = rep(1, length(es)),
           level.pi = 0.95,
           B = 100000L,
           method.tau2 = "REML",
           seed = NULL) {
    
    # validate input
    vd_PredDist(es = es, se = se, w  = w, method = method, lpi = level.pi,
                B = B, mtau2 = method.tau2)
    
    # Reproducibility under MC
    if (!is.null(seed)) {
      if (!is.wholepositivenumber(seed)) {
        warning("Seed must be a valid scalar integer.")
      }
      set.seed(seed)
    }
    
    # computation
    out <- switch(
      method,
      "PCD-fixed" = pd_cd(es = es, se = se, w = w, mtau2 = method.tau2,
                          lpi = level.pi, B = B),
      "PCD-simplified" =,
      "PCD-full" = pd_cd_tau2(es = es, se = se, w = w, lpi = level.pi,
                              method = method, B = B, mtau2 = method.tau2)
    )
    
    out$method   <- method
    out$B        <- B
    out$level_pi <- level.pi
    out$k        <- length(es)
    out$es       <- out$es
    out$se       <- se
    out$w        <- w
    class(out) <- "edgemeta"
    
    return(out)
  }

# Predictive distribution (fixed tau2)
pd_cd <- function(es,
                  se, 
                  w,
                  mtau2, 
                  lpi, 
                  B){

  # Estimate tau2
  ma <- run_metagen(es = es, se = se, mtau2 = mtau2)
  tau2 <- ma$tau2

  if (tau2 == 0L) {
    warning(
      paste0(
        "Estimated heterogeneity (tau2) is zero. ",
        "Consider using `PCD-simplified` or `PCD-full` ",
        "to account for uncertainty in the heterogeneity estimate. ",
        "Returning the ", 100 * lpi, 
        "% confidence interval for the overall effect (mu) instead."
      ),
      call. = FALSE
    )
    return(remaeffect(es = es, se = se, method = "MC", w = w, level.ci = lpi))
  }
  
  # Generate samples of mu
  s_mu <- samplemusimple(B = B, tau2 = tau2, es = es, se = se, w = w)

  # Generate samples of theta_new
  s_tn <- base::suppressWarnings(stats::rnorm(n = B, mean = s_mu, sd = sqrt(tau2)))
  
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
                       w,
                       lpi, 
                       method,
                       B,
                       mtau2) {

  # Initial tau2
  ma <- run_metagen(es = es, se = se, mtau2 = mtau2)

  # Samples of tau2
  s_tau2 <- samptau2(B = B, es = es, se = se, 
                     upper = ma$tau2 + 100 * ma$se.tau2)
  # Samples of mu 
  if (method == "PCD-simplified") {
    s_mu <- samplemusimple(B = B, tau2 = ma$tau2, es = es, se = se, w = w)
  } else if (method == "PCD-full") {
    s_mu <- samplemu(s_tau2 = s_tau2, es = es, se = se, w = w)
  }
  
  # Samples of theta_new
  s_tn <- base::suppressWarnings(stats::rnorm(n = B, mean = s_mu, sd = sqrt(s_tau2)))

  # Combine 
  s <- base::cbind(s_tau2, s_mu, s_tn)
  base::colnames(s) <- c("tau2", "mu", "theta_new")

  # Prediction interval
  PI <- stats::quantile(x = s_tn, p = (1 + lpi * c(-1, 1)) / 2, na.rm = T)
  
  # Return
  return(list(PI = PI, samples = s))
}


# Estimation of between-study heterogeneity
run_metagen <- function(es, se, mtau2) {
  
  controls <- list(
    NULL,
    list(maxiter = 10000),
    list(stepadj = 0.5),
    list(maxiter = 10000, stepadj = 0.5),
    list(maxiter = 100000, stepadj = 0.25)
  )
  
  for (ctrl in controls) {
    
    r <- tryCatch(
      meta::metagen(
        TE = es,
        seTE = se,
        random = TRUE,
        method.tau = mtau2,
        control = ctrl
      ),
      error = function(e) NULL
    )
    
    if (!is.null(r)) {
      return(r)
    }
  }
  NULL
}


# For RNG seeds
is.wholepositivenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (!is.numeric(x)) return(FALSE)
  return(is.finite(x) & abs(x - round(x)) < tol & x > 0)
}


#' Plot Predictive and Confidence Distributions from edgemeta Class Objects
#' 
#' Displays the predictive distribution for \code{param = "theta_new"} (future effect from random effects).
#' For \code{param = "mu"} and \code{param = "tau2"}, the corresponding confidence distributions are shown.
#' The latter is only valid if the estimation method in \code{PredDist} was not set to \code{"PCD-fixed"}.
#'
#' @param x An object of class \code{edgemeta}.
#' @param param The parameter to display. One of \code{theta_new} (predictive distribution), \code{mu} (confidence distribution), or \code{tau2} (confidence distribution).
#' @param ... Additional arguments passed to \code{graphics::hist()}, including optional \code{main}, \code{xlab}, and \code{ylab} to override default plot title and axis labels.
#'
#' @author David Kronthaler
#'
#' @export
#'
#' @examples
#' es <- c(0.17, 1.20, 1.10, -0.0019, -2.33)
#' se <- c(0.52, 0.93, 0.63, 0.3, 0.28)
#' 
#' # Compute predictive distribution
#' pd <- PredDist(es = es, se = se, method = "PCD-full")
#' 
#' # Plot predictive distribution
#' plot(pd, param = "theta_new", breaks = 100)
#' 
#' # Plot confidence distribution of mu
#' plot(pd, param = "mu", xlab = "mu", ylab = "c(mu)")
#' @export
plot.edgemeta <- function(x, param = c("theta_new", "mu", "tau2"), ...) {
  param <- base::match.arg(param)
  graphics::par(mar = c(5, 4.5, 4, 2))
  dots <- list(...)
  default_xlab <- switch(param,
                         "theta_new" = expression(theta[new]),
                         "mu" = expression(mu),
                         "tau2" = expression(tau^2))
  default_ylab <- switch(param,
                         "theta_new" = expression(f(theta[new])),
                         "mu" = expression(c(mu)),
                         "tau2" = expression(c(tau^2)))
  default_main <- switch(param,
                         "theta_new" = "Predictive distribution",
                         "mu" = "Confidence distribution",
                         "tau2" = "Confidence distribution")
  xlab <- if ("xlab" %in% names(dots)) dots$xlab else default_xlab
  ylab <- if ("ylab" %in% names(dots)) dots$ylab else default_ylab
  main <- if ("main" %in% names(dots)) dots$main else default_main
  col <- switch(param,
                "theta_new" = "theta_new",
                "mu" = "mu",
                "tau2" = "tau2")
  if (param == "tau2" && identical(x$method, "PCD-fixed")) {
    warning("Distribution of tau2 not available for method 'PCD-fixed'")
    return(invisible(NULL))
  }
  if (!col %in% colnames(x$samples)) {
    stop(sprintf("Column '%s' not found in x$samples", col))
  }
  hist_args <- c(
    list(x = x$samples[, col],
         xlab = xlab,
         ylab = ylab,
         main = main,
         probability = TRUE),
    dots[!names(dots) %in% c("xlab", "ylab", "main")]
  )
  do.call(graphics::hist, hist_args)
}

#' Confidence Probabilities for Future Effects in Random-Effects Meta-Analysis
#'
#' Computes the confidence probability that a future effect \code{theta_new}
#' lies within a specified interval, based on the predictive distribution from a 
#' random-effects meta-analysis.
#'
#' @param obj An object of class \code{edgemeta}.
#' @param lower Lower bound of the interval (default is \code{-Inf}).
#' @param upper Upper bound of the interval (default is \code{Inf}).
#' 
#' @author David Kronthaler
#'
#' @return Returns the confidence probability (a numeric value between 0 and 1).
#'
#' @seealso \code{\link{PredDist}}
#' @export
conf <- function(obj, lower, upper) {
  base::UseMethod("conf")
}

#' @export
conf.edgemeta <- function(obj, lower = -Inf, upper = Inf) {
  p <- base::mean(obj$samples[,"theta_new"] <= upper &
                    obj$samples[,"theta_new"] >= lower,
                  na.rm = TRUE)
  base::cat("Confidence of `theta_new` lying between ", lower, " and ", upper, ":\n", sep = "")
  base::cat(sprintf("%s\n", p))
  invisible(p)
}

#' @export
print.edgemeta <- function(x, lower = NULL, upper = NULL, ...) {
  if (!is.null(x$samples) && nrow(x$samples) > 0) {
    s <- summary_edgemeta(x)
    base::class(s) <- "data.frame"
    param_labels <- base::colnames(x$samples)
    base::rownames(s) <- rev(base::paste(c("CD", "CD", "PD")[base::seq_along(param_labels)], param_labels))
    
    base::cat("\nPredictive Distribution for Random-Effects Meta-Analysis\n\n")
    
    base::cat("Number of studies:", x$k, "\n")
    base::cat("Method:", x$method, "\n")
    
    weight_type <- "custom"
    if (isTRUE(all.equal(x$w, 1/x$se^2, tolerance = 1e-4, check.attributes = FALSE))) {
      weight_type <- "inverse squared standard errors (1/se^2)"
    } else if (isTRUE(all.equal(x$w, 1/se, tolerance = 1e-4, check.attributes = FALSE))) {
      weight_type <- "inverse standard errors (1/se)"
    } else if (all(abs(x$w - 1) < 1e-4)) {
      weight_type <- "unweighted"
    }
    base::cat("Weights:", weight_type, "\n")
    base::cat("Number of Monte Carlo samples:", 
              format(x$B, big.mark = ",", scientific = FALSE), "\n\n")
    base::cat(paste0(x$level_pi * 100, "%"), "prediction interval from", round(x$PI[1], 3),
              "to", round(x$PI[2], 3), "\n")
    
    base::cat("\nSummary of predictive distribution:\n")
    pd <- s[1, , drop = FALSE]
    rownames(pd) <- NULL
    base::print(round(pd, 3), row.names = FALSE)
    
    if (!is.null(lower) || !is.null(upper)) {
      conf.edgemeta(x, lower, upper)
    }
    
  }
}

summary_edgemeta <- function(object, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
  method <- object$method
  stats <- base::t(base::apply(object$samples, 2, function(x) {
    m <- base::mean(x, na.rm = TRUE)
    qs <- stats::quantile(x, probs = probs, na.rm = TRUE)
    c(qs[1:2], Median = qs[3], Mean = m, qs[4:5])
  }))
  base::colnames(stats) <- c(
    base::sprintf("%.1f%%", probs[1] * 100),
    base::sprintf("%.1f%%", probs[2] * 100),
    "Median", "Mean",
    base::sprintf("%.1f%%", probs[4] * 100),
    base::sprintf("%.1f%%", probs[5] * 100)
  )
  base::rownames(stats) <- base::colnames(object$samples)
  stats <- stats[base::nrow(stats):1, ]
  res <- base::as.data.frame(stats)
  base::attr(res, "method") <- method
  base::class(res) <- c("summary.edgemeta", "data.frame")
  return(res)
}