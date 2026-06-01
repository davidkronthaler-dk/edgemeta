vd_PredDist <- function(es, 
                        se, 
                        w,
                        method = c("PCD-full", "PCD-simplified", "PCD-fixed"),
                        lpi, 
                        B, 
                        mtau2 = c("REML", "PM", "DL", "ML", "HS", "SJ", "HE", "EB")
                        ) {
  
  if (!requireNamespace("meta", quietly = TRUE)) {
    stop("Package 'meta' is required. Please install it using 'install.packages(\"meta\")'.")
  } 
  
  if (is.null(method) || is.na(method) || is.nan(method)) stop("'method' must not be NULL or NA.")
  if (!is.character(method)) stop("'method' must be a character")
  method <- match.arg(method)
  if (is.null(mtau2) || is.na(mtau2) || is.nan(mtau2)) stop("'method.tau2' must not be NULL or NA.")
  if (!is.character(mtau2)) stop("'method.tau2' must be a character")
  mtau2 <- match.arg(mtau2)
  vd_es_se_w(es = es, se = se, w = w)
  vd_l(lpi)
  vd_B(B)
}

vd_remaeffect <- function(
  es,
  se,
  method = c("MC", "GAQ"),
  w,
  level.ci,
  B,
  mu0) {
  
  if (!requireNamespace("meta", quietly = TRUE)) {
    stop("Package 'meta' is required. Please install it using 'install.packages(\"meta\")'.")
  } 
  
  if (is.null(method) || is.na(method) || is.nan(method)) stop("'method' must not be NULL or NA.")
  if (!is.character(method)) stop("'method' must be a character")
  method <- match.arg(method)
  vd_es_se_w(es = es, se = se, w = w)
  vd_l(level.ci)
  vd_B(B)
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

vd_es_se_w <- function(es, se, w) {
  if (missing(es) || missing(se) || missing(w)) {
    stop("Estimates ('es'), standard errors ('se'), and weights ('w') must be provided.",
         call. = FALSE)
  }
  
  if (length(es) != length(se) ||
      length(es) != length(w) ||
      length(es) <= 1) {
    stop(paste0(
      "Estimates ('es'), standard errors ('se'), ",
      "and weights ('w') must be numeric vectors ",
      "of the same length (at least 2)."
    ), call. = FALSE)
  }
  
  if (!is.numeric(es) || !is.vector(es) ||
      anyNA(es) || any(!is.finite(es))) {
    stop("Estimates ('es') must be a numeric, finite vector without missing values.",
         call. = FALSE)
  }
  
  if (!is.numeric(se) || !is.vector(se) ||
      anyNA(se) || any(!is.finite(se))) {
    stop("Standard errors ('se') must be a numeric, finite vector without missing values.",
         call. = FALSE)
  } else if (any(se <= 0)) {
    stop("All standard errors ('se') must be greater than 0.",
         call. = FALSE)
  }
  
  if (!is.numeric(w) || !is.vector(w) ||
      anyNA(w) || any(!is.finite(w))) {
    stop("Weights ('w') must be a numeric, finite vector without missing values.",
         call. = FALSE)
  } else if (any(w <= 0)) {
    stop("All weights ('w') must be greater than 0.",
         call. = FALSE)
  }
}

vd_l <- function(l){
  if (!is.numeric(l) || length(l) != 1 || l <= 0 || l >= 1) {
    stop("Level of the confidence / prediction interval must be a scalar in (0, 1).", call. = FALSE)
  }
}

vd_B <- function(B) {
  if (!is.numeric(B) || length(B) != 1 || B %% 1 != 0 || B <= 0 || !is.finite(B)) {
    stop("Number of samples ('B') must be a finite, positive integer.", call. = FALSE)
  }
  if (B > 1e7 & B < 1e8) cat("\n'B' is set very high. Are you sure you require that many samples? Proceeding anyway:\n\n")
  if (B > 1e8) cat("\nThat many samples may be computationally stressful. Proceeding anyway:\n\n")
}

vd_tau2 <- function(tau2) {
  if (missing(tau2)) stop("tau2 must be provided.")
  if (!is.numeric(tau2) || !is.vector(tau2)) stop("tau2 must be numeric vector or scalar.")
  if (any(is.na(tau2)) || any(is.nan(tau2)) || any(is.infinite(tau2))) stop("tau2 cannot contain NA, NaN or Inf values.")
}

vd_mu0 <- function(mu0) {
  if (length(mu0) != 1) stop("mu0 must be a scalar.")
  if (!is.numeric(mu0)) stop("mu0 must be numeric.")
  if (is.na(mu0)) stop("mu0 cannot contain NA values.")
  if (is.nan(mu0)) stop("mu0 cannot contain NaN values.")
}





