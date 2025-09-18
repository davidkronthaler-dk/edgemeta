#' Confidence density (approximate) of heterogeneity parameter in random-effects meta-analysis
#'
#' The generalized heterogeneity statistic \eqn{Q(\tau^2)} (Viechtbauer, 2006) implies an approximate confidence
#' distribution of the heterogeneity parameter \eqn{\tau^2} under the normal random-effects model.
#' Applying change of variables to \eqn{Q(\tau^2)} the approximate confidence density is obtained.
#'
#' @param tau2 Value of tau2 at which to evaluate the confidence density (numeric, scalar or vector).
#' @param es Numeric vector of effect estimates from individual studies (length >= 2).
#' @param se Numeric vector of standard errors corresponding to each effect estimate (length >= 2).
#'
#' @return The value of the confidence density.
#'
#' @author David Kronthaler
#'
#' @references
#' Viechtbauer, W. (2006). *Confidence intervals for the amount of heterogeneity in meta‐analysis*. Statistics in medicine, 26(1), 37-52. https://doi.org/10.1002/sim.2514
#'
#' @export
#'
#' @examples
#' es <- c(0.17,  1.20,  1.10, -0.0019, -2.33, -1.1, -0.98)
#' se <- c(0.52, 0.93, 0.63, 0.3, 0.28, 0.45, 0.6)
#' curve(cdtau2(x, es, se), 0, 7.5)
cdtau2 <- function(tau2, es, se){
  vd_es_se(es = es, se = se)
  vd_tau2(tau2)
  cd <- ftau2(es = es, se = se, tau2 = tau2)
  return(cd)
}
