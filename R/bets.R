# Class: BettingStrategy
# purpose: virtual base class for all AGRAPA-style one-step betting strategies
#          used by BoundedModel.  Subclasses implement compute_bet().
setClass("BettingStrategy", contains = "VIRTUAL")

# Class: AGRAPABet
# purpose: adaptive betting strategy from Waudby-Smith & Ramdas (2024).
#          At time t the bet uses the mean and biased SD of the most recent
#          `window` non-NA observations X_{t-1},...,X_{max(1,t-window)}
#          via the Kelly-inspired AGRAPA formula.
# slots:
#   c      = truncation factor; bet is capped at c / (eta + eps)
#   sd_min = floor on the lagged SD to prevent division by zero
#   eps    = small constant added to eta in the cap denominator
#   window = positive integer (or Inf) number of most-recent time steps used;
#            Inf means use the full history (no windowing)
setClass(
  "AGRAPABet",
  contains = "BettingStrategy",
  slots    = c(c = "numeric", sd_min = "numeric", eps = "numeric",
               window = "numeric")
)

# Class: FixedBet
# purpose: constant betting strategy — always places the same bet lambda,
#          regardless of history.  Useful as a baseline.
# slots:
#   lambda = non-negative scalar bet
setClass(
  "FixedBet",
  contains = "BettingStrategy",
  slots    = c(lambda = "numeric")
)

# Constructor: AGRAPABet
# inputs:
#   c      = truncation factor in (0, 1], default 0.75
#   sd_min = SD floor, default 0.01
#   eps    = cap denominator nudge, default 1e-5
#   window = positive integer or Inf; number of most-recent time steps to use
#            when estimating the lagged mean and SD.  Default 30.
#            Inf disables windowing (full history, original AGRAPA behaviour).
# outputs:
#   AGRAPABet object
AGRAPABet <- function(c = 0.75, sd_min = 0.01, eps = 1e-5, window = 30) {
  if (length(c) != 1L || !is.finite(c) || c <= 0 || c > 1)
    stop("`c` must be a scalar in (0, 1].", call. = FALSE)
  if (length(sd_min) != 1L || sd_min <= 0)
    stop("`sd_min` must be a positive scalar.", call. = FALSE)
  if (length(eps) != 1L || eps < 0)
    stop("`eps` must be a non-negative scalar.", call. = FALSE)
  if (length(window) != 1L || (!is.infinite(window) && window < 1))
    stop("`window` must be a positive integer or Inf.", call. = FALSE)
  new("AGRAPABet", c = as.numeric(c), sd_min = as.numeric(sd_min),
      eps = as.numeric(eps), window = as.numeric(window))
}

# Constructor: FixedBet
# inputs:
#   lambda = non-negative scalar bet
# outputs:
#   FixedBet object
FixedBet <- function(lambda) {
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
    stop("`lambda` must be a non-negative finite scalar.", call. = FALSE)
  new("FixedBet", lambda = as.numeric(lambda))
}

# Method: compute_bet for AGRAPABet
# Uses the biased mean and SD of the most recent `window` non-NA observations.
# Returns 0 when: history is empty, or windowed mean <= eta (one-sided test).
setMethod("compute_bet", "AGRAPABet", function(strategy, history, eta) {
  h <- if (is.null(history)) numeric(0) else as.numeric(history[!is.na(history)])
  n_h <- length(h)
  if (n_h == 0L) return(0)
  if (is.finite(strategy@window) && n_h > strategy@window)
    h <- tail(h, as.integer(strategy@window))
  n_h     <- length(h)
  mu      <- mean(h)
  sd_hat  <- if (n_h > 1L) sqrt(mean((h - mu)^2)) else 0
  sd_hat  <- max(sd_hat, strategy@sd_min)
  lam_raw <- (mu - eta) / (sd_hat^2 + (mu - eta)^2)
  max(0, min(lam_raw, strategy@c / (eta + strategy@eps)))
})

# Method: compute_bet for FixedBet
# Always returns the stored lambda, ignoring history.
setMethod("compute_bet", "FixedBet", function(strategy, history, eta) {
  strategy@lambda
})

# Function: compute_kelly_optimal_bet
# purpose: find the Kelly-optimal bet lambda* that maximises expected log-growth
#          E_p[log(1 + lambda*(X - eta))] over the feasible interval (0, 1/eta).
#          Equivalently, finds the root of the derivative
#            f(lambda) = sum_i p_i * (l_i - eta) / (1 + lambda*(l_i - eta)) = 0.
#          If the weighted mean is <= eta (no evidence against the null), returns 0.
#          If all locations l_i >= eta (no downside), returns the cap 1/eta.
# inputs:
#   l_i = numeric vector of observation locations in [0, 1]
#   eta = numeric scalar null mean upper bound in (0, 1)
#   p_i = numeric non-negative weight vector summing to 1 (default: uniform 1/n)
# outputs:
#   numeric scalar lambda* in [0, 1/eta)
compute_kelly_optimal_bet <- function(l_i, eta, p_i = NULL) {
  n <- length(l_i)
  if (!is.numeric(l_i) || n < 1L)
    stop("`l_i` must be a non-empty numeric vector.", call. = FALSE)
  if (!is.numeric(eta) || length(eta) != 1L || eta <= 0 || eta >= 1)
    stop("`eta` must be a scalar in (0, 1).", call. = FALSE)
  if (is.null(p_i)) {
    p_i <- rep(1 / n, n)
  } else {
    if (!is.numeric(p_i) || length(p_i) != n || any(p_i < 0))
      stop("`p_i` must be a non-negative numeric vector of the same length as `l_i`.",
           call. = FALSE)
    p_i <- p_i / sum(p_i)
  }

  # If weighted mean <= eta, optimal bet is 0 (no evidence against null)
  if (sum(p_i * l_i) <= eta) return(0)

  # f(lambda) = E_p[(X - eta) / (1 + lambda*(X - eta))]
  # Strictly decreasing in lambda when E[X] > eta.
  f <- function(lam) sum(p_i * (l_i - eta) / (1 + lam * (l_i - eta)))

  # Feasibility constraint from x = 0 on [0,1]: 1 - lambda*eta >= 0 => lambda <= 1/eta.
  # Stay strictly below to avoid singularity when any l_i = 0.
  lam_cap <- (1 - sqrt(.Machine$double.eps)) / eta

  f_cap <- tryCatch(f(lam_cap), error = function(e) NA_real_)

  # If f is still non-negative at the cap (all l_i >= eta), return the cap.
  if (is.na(f_cap) || f_cap >= 0) return(lam_cap)

  # Otherwise find the unique root in (0, lam_cap).
  uniroot(f, lower = 0, upper = lam_cap,
          tol = .Machine$double.eps^0.5)$root
}
