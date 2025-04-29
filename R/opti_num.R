opti_num <- function(es, se, point = TRUE, ci = TRUE, level.ci = 0.95) {

  # Compute the point estimate
  if (point) {

    point_est <- stats::uniroot(f = function(x) pfct_edge_cpp(x, es, se) - 0.5,
                                interval = c(base::min(es), base::max(es)))$root
  } else {
    point_est <- NA
  }

  # Function to compute the CI bounds
  if (ci) {
    interval_ci <- c(base::min(es) - 2 * base::max(se), base::max(es) + 2 * base::max(se))
    alpha <- 1 - level.ci
    ci_lower <- stats::uniroot(f = function(x) pfct_edge_cpp(x, es, se) - alpha / 2,
                               interval = interval_ci, extendInt = "yes")$root
    ci_upper <- stats::uniroot(f = function(x) pfct_edge_cpp(x, es, se) - (1 - alpha / 2),
                               interval = interval_ci, extendInt  = "yes")$root
  } else {
    ci_lower <- NA
    ci_upper <- NA
  }

  return(list(point_estimate = point_est, ci_lower = ci_lower, ci_upper = ci_upper))
}
