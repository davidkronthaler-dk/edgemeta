vd_PredDist <- function(es, 
                        se, 
                        method = c("FullCD", "SimplifiedCD", "FixedTau2"),
                        lpi, 
                        ns, 
                        mtau2 = c("REML", "PM", "DL", "ML", "HS", "SJ", "HE", "EB")
                        ) {
  
  if (!requireNamespace("meta", quietly = TRUE)) {
    stop("Package 'meta' is required. Please install it using 'install.packages(\"meta\")'.")
  } 
  
  method <- match.arg(method)
  mtau2 <- match.arg(mtau2)
  vd_es_se(es = es, se = se)
  vd_l(lpi)
  vd_ns(ns)
}

vd_remaeffect <- function(
  es,
  se,
  method = c("MC", "GAQ", "NHEU"),
  level.ci,
  n_samples,
  mu0,
  method.tau2 = c("REML", "PM", "DL", "ML", "HS", "SJ", "HE", "EB", "estimate")) {
  
  if (!requireNamespace("meta", quietly = TRUE)) {
    stop("Package 'meta' is required. Please install it using 'install.packages(\"meta\")'.")
  } 
  
  method <- match.arg(method)
  method.tau2 <- match.arg(method.tau2)
  vd_es_se(es = es, se = se)
  vd_l(level.ci)
  vd_ns(n_samples)
  vd_mu0(mu0)
}

vd_crps <- function(s, tn) {
  if (!is.numeric(s)) {
    stop("Argument 's' must be a numeric vector of predictive samples.")
  }
  if (!is.numeric(tn)) {
    stop("Argument 'tn' must be a numeric vector of truth values.")
  }
  if (length(s) == 0) {
    stop("Argument 's' must not be empty.")
  }
  if (length(tn) == 0) {
    stop("Argument 'tn' must not be empty.")
  }
  if (anyNA(s) || any(!is.finite(s))) {
    stop("Argument 's' must not contain NA, NaN, or Inf values.")
  }
  if (anyNA(tn) || any(!is.finite(tn))) {
    stop("Argument 'tn' must not contain NA, NaN, or Inf values.")
  }
}

vd_es_se <- function(es, se) {
  if (missing(es) || missing(se)) {
    stop("Estimates ('es') and standard errors ('se') must be provided.", call. = FALSE)
  }
  
  if (length(es) != length(se) || length(es) <= 1) {
    stop("Estimates ('es') and standard errors ('se') must be numeric vectors of the same length (at least 2).", call. = FALSE)
  }
  
  if (!is.numeric(es) || !is.vector(es) || anyNA(es)) {
    stop("Estimates ('es') must be a numeric vector without missing values.", call. = FALSE)
  }
  
  if (!is.numeric(se) || !is.vector(se) || anyNA(se)) {
    stop("Standard errors ('se') must be a numeric vector without missing values.", call. = FALSE)
  } else if (any(se <= 0)) {
    stop("All standard errors ('se') must be greater than 0.", call. = FALSE)
  }
}

vd_l <- function(l){
  if (!is.numeric(l) || length(l) != 1 || l <= 0 || l >= 1) {
    stop("Level of the confidence / prediction interval must be a scalar in (0, 1).", call. = FALSE)
  }
}

vd_ns <- function(ns) {
  if (!is.numeric(ns) || length(ns) != 1 || ns %% 1 != 0 || ns <= 0) {
    stop("Number of samples ('n_samples') must be a positive integer.", call. = FALSE)
  }
}

vd_tau2 <- function(tau2) {
  if (missing(tau2)) stop("tau2 must be provided.")
  if (!is.numeric(tau2)) stop("tau2 must be numeric.")
  if (any(is.na(tau2))) stop("tau2 cannot contain NA values.")
  if (any(is.nan(tau2))) stop("tau2 cannot contain NaN values.")
  if (any(is.infinite(tau2))) stop("tau2 cannot contain Inf values.")
}

vd_mu0 <- function(mu0) {
  if (!is.numeric(mu0)) stop("mu0 must be numeric.")
  if (is.na(mu0)) stop("mu0 cannot contain NA values.")
  if (is.nan(mu0)) stop("mu0 cannot contain NaN values.")
  if (length(mu0) != 1) stop("mu0 must be a scalar.")
}





