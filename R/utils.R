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

# Internal helper: .lag_mean_sd
# purpose: compute lagged running mean and biased SD of a sequence
# inputs:
#   x    = numeric vector of observations
#   mu_0 = numeric scalar, default mean inserted at position 1 (before any data)
#   sd_0 = numeric scalar, default SD inserted at position 1 (before any data)
# outputs:
#   list with elements:
#     mean = numeric length-n vector of lagged running means
#     sd   = numeric length-n vector of lagged running biased SDs
.lag_mean_sd <- function(x, mu_0 = 0.5, sd_0 = 0.25) {
  n <- length(x)
  if (n == 0L) return(list(mean = numeric(0), sd = numeric(0)))
  t       <- seq_len(n)
  mu_hat  <- cumsum(x) / t
  sum_sq  <- cumsum(x^2)
  var_hat <- pmax(sum_sq / t - mu_hat^2, 0)
  lag_mean <- c(mu_0, mu_hat[-n])
  lag_sd   <- c(sd_0, sqrt(var_hat[-n]))
  list(mean = lag_mean, sd = lag_sd)
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
