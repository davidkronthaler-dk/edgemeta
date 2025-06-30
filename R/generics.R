#' Plotting Predictive and Confidence Distributions from 'metaprediction' Class Objects
#'
#' @param obj Object of class 'metaprediction', obtained from the 'PredDist' function.
#' @param param The parameter to display. One of \code{"theta_new"}, \code{"mu"}, or \code{"tau2"}.
#' @param ... Additional arguments passed to \code{graphics::hist()}, including optional \code{main}, \code{xlab}, and \code{ylab} to override default plot title and axis labels.
#'
#' @details
#' Displays the predictive distribution for \code{param = "theta_new"} (future effect from random effects).
#' For \code{"mu"} and \code{"tau2"}, the corresponding confidence distributions are shown.
#' \code{"tau2"} is only valid if the estimation method was not set to \code{"FixedTau2"}.
#'
#' @author David Kronthaler
#'
#' @export
#'
#' @examples
#' es <- c(0.17, 1.20, 1.10, -0.0019, -2.33)
#' se <- c(0.52, 0.93, 0.63, 0.3, 0.28)
#' pd <- PredDist(es = es, se = se, method = "FullCD")
#' plot(pd, param = "theta_new", breaks = 100)
#' plot(pd, param = "mu", main = "CD for Mu", xlab = "Mu values", ylab = "Density")
#' @export
plot.metaprediction <- function(obj, param = c("theta_new", "mu", "tau2"), ...) {
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
  if (!col %in% colnames(obj$samples)) {
    stop(sprintf("Column '%s' not found in obj$samples", col))
  }
  if (param == "tau2" && identical(attr(obj, "method"), "FixedTau2")) {
    warning("Distribution of Tau2 not available for method 'FixedTau2'")
    return(invisible(NULL))
  }
  hist_args <- c(
    list(x = obj$samples[, col],
         xlab = xlab,
         ylab = ylab,
         main = main,
         probability = TRUE),
    dots[!names(dots) %in% c("xlab", "ylab", "main")]
  )
  do.call(graphics::hist, hist_args)
}

#' @export
summary.metaprediction <- function(obj, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
  method <- base::attr(obj, "method")
  stats <- base::t(base::apply(obj$samples, 2, function(x) {
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
  base::rownames(stats) <- base::colnames(obj$samples)
  stats <- stats[base::nrow(stats):1, ]
  res <- base::as.data.frame(stats)
  base::attr(res, "method") <- method
  base::class(res) <- c("summary.metaprediction", "data.frame")
  return(res)
}

#' @export
print.summary.metaprediction <- function(x, ...) {
  base::cat("Summary of metaprediction object\n")
  base::cat("Method:", attr(x, "method"), "\n\n")
  base::NextMethod("print")
}

#' Confidence Probability for Future Effect in Meta-Analysis
#'
#' Computes the confidence probability that a future outcome \code{theta_new}
#' lies within a specified interval, based on the predictive distribution from a 
#' random-effects meta-analysis.
#'
#' @param obj An object of class \code{metaprediction}.
#' @param lower Lower bound of the interval (default is \code{-Inf}).
#' @param upper Upper bound of the interval (default is \code{Inf}).
#'
#' @return Invisibly returns the confidence probability (a numeric value between 0 and 1).
#' Prints the result to the console.
#'
#' @seealso \code{\link{PredDist}}
#' @export
conf <- function(obj, ...) {
  base::UseMethod("conf")
}

#' @export
conf.metaprediction <- function(obj, lower = -Inf, upper = Inf) {
  p <- base::mean(obj$samples[,"theta_new"] <= upper &
              obj$samples[,"theta_new"] >= lower,
            na.rm = TRUE)
  base::cat("Confidence of `theta_new` lying between", lower, "and", upper, ":")
  base::cat(sprintf("%s\n", p))
  invisible(p)
}

#' @export
print.metaprediction <- function(obj, lower = NULL, upper = NULL, ...) {
  if (!is.null(obj$samples) && nrow(obj$samples) > 0) {
    s <- summary.metaprediction(obj)
    base::class(s) <- "data.frame"
    param_labels <- base::colnames(obj$samples)
    base::rownames(s) <- rev(base::paste(c("CD", "CD", "PD")[base::seq_along(param_labels)], param_labels))
    
    base::cat("\n=================== MetaPrediction Summary ==================\n\n")
    
    base::cat("Number of studies:", attr(obj, "k"), "\n")
    base::cat("Method: ", attr(obj, "method"), "\n")
    base::cat("Number of Monte Carlo samples:", attr(obj, "n_samples"), "\n\n")
    base::cat(attr(obj, "level_pi") * 100, "% prediction interval from", round(obj$PI[1], 3),
              "to", round(obj$PI[2], 3), "\n")
    
    base::cat("\nSummary of predictive distribution (PD) and confidence distributions (CD)\n")
    base::print(s)
    
    base::cat("\nConfidence calculations:\n")
    conf.metaprediction(obj, 0, Inf)
    conf.metaprediction(obj, -Inf, 0)
    
    if (!is.null(lower) || !is.null(upper)) {
      conf.metaprediction(obj, lower, upper)
    }
    
    base::cat("\n=============================================================\n")
  }
}

