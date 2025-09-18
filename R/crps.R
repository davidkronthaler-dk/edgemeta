#' Continuous Ranked Probability Score (CRPS)
#'
#' Computes the Continuous Ranked Probability Score (CRPS) for a predictive distribution using a Monte Carlo (MC) approximation.
#' The MC approximation follows the approach for the computation of the *energy score* proposed by Gneiting et al. (2008).
#' To ensure computational efficiency, especially when handling large numbers of predictive samples or multiple future truths,
#' the Monte Carlo approximation is implemented with a C++ backend.
#' If multiple future truths are provided, the function returns the mean CRPS across all provided future truths.
#'
#' @param s A numeric vector of samples from the predictive distribution.
#' @param tn A numeric vector of future truths.
#'
#' @return A single numeric value representing the mean CRPS across all future truths.
#' 
#' 
#' @details
#' A likely misconception is that supplying the same vector for both the predictive samples
#' \code{s} and the truths \code{tn} should result in a CRPS of zero. This is not the case,
#' unless the vector contains only a single unique element. The reason is that each observed
#' value in \code{tn} is evaluated against the *entire* predictive distribution in \code{s}.
#' By contrast, if a truth value is evaluated against a predictive distribution that
#' places all of its mass at that same value (i.e., a point mass), the CRPS is indeed zero.
#'
#' @author David Kronthaler
#'
#' @references
#' Gneiting, T., Stanberry, L. I., Grimit, E. P., Held, L., & Johnson, N. A. (2008).
#' Assessing probabilistic forecasts of multivariate quantities, with an application
#' to ensemble predictions of surface winds. Test, 17, 211–235.
#' 10.1007/s11749-008-0114-x
#' 
#' @examples
#' # True data-generating distribution
#' truth <- rnorm(100, 0, 3)
#' pd <- rnorm(1000, 0.2, 2.8)
#'
#' # Compute CRPS
#' crps(pd, truth)
#'
#' @export
crps <- function(s, tn) {
  vd_crps(s = s, tn = tn)
  r <- crpsCPP(s = s, t = tn)
  return(r)
}

