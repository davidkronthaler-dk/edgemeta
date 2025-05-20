#' Plotting Predictive and Confidence Distributions of 'metaprediction' Class Objects
#'
#' @param obj Object of class 'metaprediction', obtained from the 'PredDist' function.
#' @param param The parameter to display. One of \code{"theta_new"}, \code{"mu"}, or \code{"tau2"}.
#' @param ... Additional arguments passed to \code{graphics::hist()}, including optional \code{main}, \code{xlab}, and \code{ylab} to override default plot title and axis labels.
#'
#' @details
#' Displays the predictive distribution for \code{param = "theta_new"} (future effect from random effects).
#' For \code{"mu"} and \code{"tau2"}, the corresponding confidence distributions (CDs) are shown.
#' \code{"tau2"} is only valid if the estimation method was not set to \code{"FixedTau2"}.
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
                         "mu" = expression(CD(mu)),
                         "tau2" = expression(CD(tau^2)))
  default_main <- switch(param,
                         "theta_new" = "Predictive distribution",
                         "mu" = "Confidence density",
                         "tau2" = "Confidence density")
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

#' @export
prob <- function(obj, ...) {
  base::UseMethod("prob")
}

#' @export
prob.metaprediction <- function(obj, lower = -Inf, upper = Inf) {
  p <- base::mean(obj$samples[,"theta_new"] <= upper &
              obj$samples[,"theta_new"] >= lower,
            na.rm = TRUE)
  base::cat("Probability of `theta_new` lying between", lower, "and", upper, ":\n")
  base::cat(sprintf("%s\n", p))
  invisible(p)
}

#' @export
print.metaprediction <- function(obj, lower = NULL, upper = NULL, ...) {
  s <- summary.metaprediction(obj)
  base::class(s) <- "data.frame"
  param_labels <- base::colnames(obj$samples)
  base::rownames(s) <- base::paste(c("PD", "CD", "CD")[base::seq_along(param_labels)], param_labels)
  base::cat("\n============================== MetaPrediction Summary ==============================\n")
  base::cat("\nSummary of predictive distribution (PD) and confidence distributions (CD)\n\n")
  base::print(s)
  base::cat("\nProbability calculations:\n")
  prob.metaprediction(obj, 0, )
  base::cat("\n")
  prob.metaprediction(obj, Inf, 0)
  if (!base::is.null(lower) || !base::is.null(upper)) {
    prob.metaprediction(obj, lower, upper)
  }
  base::cat("\n=======================================================================================\n")
}

