# Internal helper: .assert_numeric_vector
# purpose: validate finite numeric vector inputs used throughout package internals
# inputs:
#   x    = object to validate
#   name = character, input name for error messages
# outputs:
#   NULL (invisibly) on success; throws error on invalid input
.assert_numeric_vector <- function(x, name) {
  if (!is.numeric(x) || length(x) < 1L || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a non-empty finite numeric vector.", name), call. = FALSE)
  }
}

# Internal helper: .logsumexp
# purpose: stable computation of log(sum(exp(x)))
# inputs:
#   x = numeric vector
# outputs:
#   numeric scalar log-sum-exp
.logsumexp <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Helper: increments_to_tsm
# purpose: construct a TSM-style path from one-step increments
# inputs:
#   increments  = numeric length-N vector of increments (or log-increments if log=TRUE)
#   initial     = numeric scalar initial wealth (default 1)
#   running_max = logical; if TRUE, return running max of cumulative product
#   log         = logical; if TRUE, interpret `increments` as log-increments and accumulate by summation
# outputs:
#   numeric length-N vector with cumulative path
increments_to_tsm <- function(increments, initial = 1, running_max = FALSE, log = FALSE) {
  .assert_numeric_vector(increments, "increments")
  if(!log && any(increments < 0)){stop("increments are negative but log is FALSE!")}
  if (length(initial) != 1L || !is.finite(initial) || initial <= 0) {
    stop("`initial` must be a positive finite scalar.", call. = FALSE)
  }
  if (log) {
    path <- log(initial) + cumsum(increments)
  } else {
    path <- initial * cumprod(pmax(increments, .Machine$double.eps))
  }
  if (isTRUE(running_max)) {
    return(cummax(path))
  }
  path
}

# Internal helper: .geometric_spending
# purpose: create the first n terms of a geometric spending schedule
# inputs:
#   n = positive integer horizon
#   p = geometric parameter in (0,1), default 0.01 (mean 100)
# outputs:
#   numeric length-n vector with entries p*(1-p)^(t-1) (no renormalization)
.geometric_spending <- function(n, p = 0.01) {
  idx <- seq_len(n)
  p * (1 - p)^(idx - 1L)
}

# Internal helper: .composition_grid
# purpose: enumerate integer compositions of m across K dimensions
# inputs:
#   K = positive integer, dimension
#   m = nonnegative integer, total mass to allocate
# outputs:
#   numeric matrix where each row sums to m (used for simplex grid construction)
.composition_grid <- function(K, m) {
  if (K == 1L) {
    return(matrix(m, nrow = 1L))
  }
  out <- list()
  recurse <- function(prefix, left, dim_left) {
    if (dim_left == 1L) {
      out[[length(out) + 1L]] <<- c(prefix, left)
      return(invisible(NULL))
    }
    for (v in 0:left) {
      recurse(c(prefix, v), left - v, dim_left - 1L)
    }
    invisible(NULL)
  }
  recurse(numeric(0), m, K)
  do.call(rbind, out)
}
