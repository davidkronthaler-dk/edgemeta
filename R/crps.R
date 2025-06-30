#' Continuous Ranked Probability Score (CRPS) for Sampled Predictive Distributions
#'
#' This function computes the Continuous Ranked Probability Score (CRPS) for a predictive distribution,
#' following the definition by Gneiting et al. (2008). It uses a C++ backend for computational efficiency.
#'
#' @param s A numeric vector of samples from the predictive distribution.
#' @param tn A numeric vector of values from the true distribution of future effects.
#'
#' @return A single numeric value representing the mean CRPS across all observations.
#'
#' @author David Kronthaler
#' 
#' @name crps
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
crps
