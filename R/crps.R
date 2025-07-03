#' Continuous Ranked Probability Score (CRPS)
#'
#' Computes the Continuous Ranked Probability Score (CRPS) for a predictive distribution using a Monte Carlo (MC) approximation.
#' The MC approximation follows the approach for the computation of the *energy score* proposed by Gneiting et al. (2008).
#' To ensure computational efficiency, especially when handling large numbers of predictive samples or multiple future effect samples,
#' the Monte Carlo approximation is implemented with a C++ backend.
#' If multiple future effects are provided, the function returns the mean CRPS across all provided samples.
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
