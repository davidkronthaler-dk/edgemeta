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
  if (param == "tau2" && identical(attr(x, "method"), "PCD-fixed")) {
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
    
    base::cat("Number of studies:", attr(x, "k"), "\n")
    base::cat("Method:", attr(x, "method"), "\n")
    base::cat("Number of Monte Carlo samples:", 
              format(attr(x, "B"), big.mark = ",", scientific = FALSE), "\n\n")
    base::cat(paste0(attr(x, "level_pi") * 100, "%"), "prediction interval from", round(x$PI[1], 3),
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
  method <- base::attr(object, "method")
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
