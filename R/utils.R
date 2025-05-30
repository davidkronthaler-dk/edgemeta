# Validate input
validate_input <- function(es, se, method = c("FixedTau2", "SimplifiedCD", "FullCD", "estimate"),
                           lpi, ns, 
                           mtau2 = c("REML", "PM", "DL", "ML", "HS", "SJ", "HE", "EB", "estimate")) {
  
  # Check if meta-package is available
  if (!requireNamespace("meta", quietly = TRUE)) {
    stop("Package 'meta' is required. Please install it using 'install.packages(\"meta\")'.")
  } 
  
  # Validate character arguments
  method <- match.arg(method)
  mtau2 <- match.arg(mtau2)
  
  # Check if required inputs are present
  if (missing(es) || missing(se)) {
    stop("Estimates ('es') and standard errors ('se') must be provided.", call. = FALSE)
  }
  
  # Check if estimates and standard errors match
  if (length(es) != length(se) || length(es) <= 1) {
    stop("Estimates ('es') and standard errors ('se') must be numeric vectors of the same length (at least 2).", call. = FALSE)
  }
  
  # Validate estimates
  if (!is.numeric(es) || !is.vector(es) || anyNA(es)) {
    stop("Estimates ('es') must be a numeric vector without missing values.", call. = FALSE)
  }
  
  # Validate standard errors
  if (!is.numeric(se) || !is.vector(se) || anyNA(se)) {
    stop("Standard errors ('se') must be a numeric vector without missing values.", call. = FALSE)
  } else if (any(se <= 0)) {
    stop("All standard errors ('se') must be greater than 0.", call. = FALSE)
  }
  
  # Validate level.pi
  if (!is.numeric(lpi) || length(lpi) != 1 || lpi <= 0 || lpi >= 1) {
    stop("Level of the prediction interval ('level.pi') must be a scalar in (0, 1).", call. = FALSE)
  }
  
  # Validate Number of samples
  if (!is.numeric(ns) || length(ns) != 1 || ns %% 1 != 0 || ns <= 0) {
    stop("Number of samples ('n_samples') must be a positive integer.", call. = FALSE)
  }

}










