# Class: TestSupermartingale
# purpose: base class for objects that convert data streams into evidence processes
# slots:
#   model   = Model subclass object used for likelihood increments
#   initial = numeric scalar, initial wealth/evidence
setClass("TestSupermartingale", slots = c(model = "Model", initial = "numeric"))

# Class: TSM
# purpose: generic TSM class built from a GaussianVARModel, BernoulliModel, or any Model subclass
setClass("TSM", contains = "TestSupermartingale")

# Backward-compatible alias class
setClass("SimpleVsSimpleTSM", contains = "TSM")

# Constructor: TSM
# inputs:
#   model   = Model subclass instance
#   initial = numeric scalar > 0
# outputs:
#   TSM object
TSM <- function(model, initial = 1) {
  stopifnot(is(model, "Model"), length(initial) == 1L, initial > 0)
  new("TSM", model = model, initial = initial)
}

# Backward-compatible constructor alias
SimpleVsSimpleTSM <- function(model, initial = 1) TSM(model = model, initial = initial)

# Class: ConformalTSM (placeholder — not yet implemented)
# purpose: placeholder for conformal martingale evidence processes
setClass("ConformalTSM",
         slots = c(calibration = "numeric", score_fn = "function", initial = "numeric"))

ConformalTSM <- function(...) stop("ConformalTSM is not yet implemented.", call. = FALSE)

# ---- vectorised fast-path helpers ------------------------------------------

# Fast path: BoundedModel + AGRAPABet — O(N) via sliding prefix sums.
# Supports hard windowing: at time t the lagged statistics use only the
# most recent `window` time steps (not the full history).  Inf = no windowing.
# NA handling: missing observations are excluded from the lagged statistics
# and contribute an increment of 1 (no evidence update at that step).
.bounded_increments_fast <- function(model, x, log) {
  bets <- model@bets   # AGRAPABet
  x    <- as.numeric(x)
  .assert_numeric_vector(x, "x")
  n   <- length(x)
  eta <- model@eta
  win <- if (is.finite(bets@window)) as.integer(bets@window) else n

  not_na  <- !is.na(x)
  x_obs   <- ifelse(not_na, x,    0)   # NA → 0 so cumsum stays finite
  x2_obs  <- ifelse(not_na, x^2,  0)

  # Build length-(n+1) prefix arrays: prefix[k] = aggregate over x[1..k-1].
  # prefix[1] = 0 (empty history), prefix[t] = aggregate over x[1..t-1].
  full_n  <- c(0L, cumsum(as.integer(not_na)))
  full_s  <- c(0,  cumsum(x_obs))
  full_s2 <- c(0,  cumsum(x2_obs))

  # At time t the window covers x[max(1, t-W)..t-1] in the time-step sense.
  # Using prefix arrays: lag quantity = prefix[t] - prefix[max(1, t-win)].
  t_idx  <- seq_len(n)
  start  <- pmax(1L, t_idx - win)   # first prefix index to subtract

  lag_n  <- full_n[t_idx]  - full_n[start]
  lag_s  <- full_s[t_idx]  - full_s[start]
  lag_s2 <- full_s2[t_idx] - full_s2[start]

  # Lagged mean and biased SD (safe at lag_n=0: lam forced to 0 below)
  lag_mu <- lag_s / pmax(lag_n, 1L)
  lag_sd <- pmax(sqrt(pmax(lag_s2 / pmax(lag_n, 1L) - lag_mu^2, 0)),
                 bets@sd_min)

  lam_raw          <- (lag_mu - eta) / (lag_sd^2 + (lag_mu - eta)^2)
  lam              <- pmax(0, pmin(lam_raw, bets@c / (eta + bets@eps)))
  lam[lag_n == 0L] <- 0   # empty window → no bet

  # Missing observation at t → increment 1 (no evidence update at that step)
  inc <- ifelse(not_na,
                pmax(1 + lam * (x - eta), .Machine$double.eps),
                1)
  if (log) base::log(inc) else inc
}

# Fast path: BoundedModel + FixedBet — fully vectorised, O(N).
# Lambda is constant so each increment is computed in one expression, no loop.
.fixedbet_increments_fast <- function(model, x, log) {
  lam <- model@bets@lambda
  eta <- model@eta
  x   <- as.numeric(x)
  .assert_numeric_vector(x, "x")
  inc <- ifelse(!is.na(x),
                pmax(1 + lam * (x - eta), .Machine$double.eps),
                1)
  if (log) base::log(inc) else inc
}

# Fast path: BoundedModel + EWMABet — O(N) sequential EWMA update.
# EWMA recursion (predictable: uses X_{t-1} to update before observing X_t):
#   mu_t = (1-rho)*mu_{t-1} + rho*X_{t-1}
#   v_t  = (1-rho)*v_{t-1}  + rho*(X_{t-1} - mu_{t-1})^2
# NA observations are skipped: EWMA state is not updated on NA.
.ewma_increments_fast <- function(model, x, log) {
  bets   <- model@bets   # EWMABet
  x      <- as.numeric(x)
  .assert_numeric_vector(x, "x")
  n      <- length(x)
  eta    <- model@eta          # scalar or length-n vector
  is_vec <- length(eta) > 1L
  rho    <- bets@rho
  mu     <- bets@mu_init
  v      <- bets@v_init

  inc <- numeric(n)
  for (t in seq_len(n)) {
    xt    <- x[t]
    eta_t <- if (is_vec) eta[t] else eta
    # Bet using current predictable (lagged) EWMA estimates
    d      <- mu - eta_t
    sd_hat <- max(sqrt(v), bets@sd_min)
    lam    <- max(0, min(d / (sd_hat^2 + d^2), bets@c / (eta_t + bets@eps)))
    if (!is.na(xt)) {
      inc[t] <- max(1 + lam * (xt - eta_t), .Machine$double.eps)
      # Update EWMA state: new estimate incorporates X_t for the NEXT step
      v  <- (1 - rho) * v  + rho * (xt - mu)^2
      mu <- (1 - rho) * mu + rho * xt
    } else {
      inc[t] <- 1   # missing → no evidence, EWMA state unchanged
    }
  }
  if (log) base::log(inc) else inc
}

# Fast path: GaussianVARModel — O(N * K^2) vectorised log-LR computation.
# Builds the full N-by-K conditional mean matrix for pre and post in one pass,
# then computes all N Mahalanobis distances at once.
.gvar_increments_fast <- function(model, x, log) {
  if (is.null(dim(x))) x <- matrix(as.numeric(x), ncol = 1L)
  .assert_numeric_vector(as.numeric(x), "x")
  N  <- nrow(x)
  K  <- ncol(x)
  p  <- length(model@Phi_pre)
  I_K <- if (K == 1L) matrix(1) else diag(K)

  # Intercept vectors:  nu = (I - sum Phi_j) * mean  (for IID: nu = mean)
  nu_pre  <- if (p > 0L) as.numeric((I_K - Reduce("+", model@Phi_pre))  %*% model@mean_pre)
             else model@mean_pre
  nu_post <- if (p > 0L) as.numeric((I_K - Reduce("+", model@Phi_post)) %*% model@mean_post)
             else model@mean_post

  # Initialise conditional mean matrices as N x K broadcasts of the intercept
  cond_pre  <- matrix(nu_pre,  nrow = N, ncol = K, byrow = TRUE)
  cond_post <- matrix(nu_post, nrow = N, ncol = K, byrow = TRUE)

  if (p > 0L) {
    # Prepend p rows of x0 for lag initialisation
    x_pad <- rbind(matrix(model@x0, nrow = p, ncol = K, byrow = TRUE), x)
    for (j in seq_len(p)) {
      # Rows of x_pad corresponding to lag j at times t = 1..N
      lag_mat    <- x_pad[(p - j + 1L):(p + N - j), , drop = FALSE]  # N x K
      cond_pre   <- cond_pre  + lag_mat %*% t(model@Phi_pre[[j]])
      cond_post  <- cond_post + lag_mat %*% t(model@Phi_post[[j]])
    }
  }

  # Vectorised log MVN density: log_const - 0.5 * rowSums((resid %*% Sigma_inv) * resid)
  .log_mvn_rows <- function(x_mat, mean_mat, Sigma) {
    log_const <- -0.5 * (K * log(2 * pi) +
                   as.numeric(determinant(Sigma, logarithm = TRUE)$modulus))
    S_inv <- solve(Sigma)
    resid <- x_mat - mean_mat
    log_const - 0.5 * rowSums((resid %*% S_inv) * resid)
  }

  log_inc <- .log_mvn_rows(x, cond_post, model@Sigma_post) -
             .log_mvn_rows(x, cond_pre,  model@Sigma_pre)
  if (log) log_inc else pmax(exp(log_inc), .Machine$double.eps)
}

# ---- compute_increments for TSM --------------------------------------------

# Method: compute_increments for TSM
# inputs:
#   object = TSM object
#   x      = numeric vector (univariate stream) or numeric N-by-K matrix (multivariate stream)
#   log    = logical; if TRUE return log-increments
# outputs:
#   numeric length-N vector with one-step increments (or log-increments)
setMethod("compute_increments", "TSM", function(object, x, log = FALSE) {
  model <- object@model

  # Vectorised fast paths (bypass the generic observation-by-observation loop)
  if (is(model, "GaussianVARModel"))
    return(.gvar_increments_fast(model, x, log))
  if (is(model, "BoundedModel") && is.null(dim(x))) {
    # Vector eta: validate length matches x before dispatching.
    if (length(model@eta) > 1L && length(model@eta) != length(x))
      stop("length(eta) must equal length(x) when eta is a vector.",
           call. = FALSE)
    if (is(model@bets, "AGRAPABet"))
      return(.bounded_increments_fast(model, x, log))
    if (is(model@bets, "EWMABet"))
      return(.ewma_increments_fast(model, x, log))
    if (is(model@bets, "FixedBet"))
      return(.fixedbet_increments_fast(model, x, log))
  }

  # Generic loop — used for BernoulliModel, BoundedModel with non-AGRAPA bets,
  # and any future Model subclasses.
  # History is accumulated one observation at a time so that likelihood_increment
  # has access to all past data (supports arbitrary dependence on lagged values).
  # Missing observations (NA) yield increment 1 and do not update history.
  if (is.null(dim(x))) {
    .assert_numeric_vector(as.numeric(x), "x")
    x       <- as.numeric(x)
    n       <- length(x)
    out     <- rep(NA_real_, n)
    history <- numeric(0)

    for (t in seq_len(n)) {
      if (is.na(x[t])) {
        out[t] <- if (log) 0 else 1   # missing → no evidence
        next
      }
      out[t]  <- likelihood_increment(model, x = x[t], history = history, log = log)
      history <- c(history, x[t])
    }
    return(out)
  }

  if (!is.matrix(x) || !is.numeric(x))
    stop("`x` must be a numeric vector or matrix.", call. = FALSE)
  if (!is(model, "MultivariateModel"))
    stop("Matrix input requires a model inheriting from `MultivariateModel`.", call. = FALSE)

  n       <- nrow(x)
  out     <- numeric(n)
  history <- matrix(numeric(0), nrow = 0, ncol = ncol(x))

  for (t in seq_len(n)) {
    out[t]  <- likelihood_increment(model, x = x[t, ], history = history, log = log)
    history <- rbind(history, x[t, , drop = FALSE])
  }
  out
})

# Method: compute_increments for ConformalTSM (not yet implemented)
setMethod("compute_increments", "ConformalTSM", function(object, x, log = FALSE) {
  stop("ConformalTSM is not yet implemented.", call. = FALSE)
})

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
  if (!log && any(increments < 0)) stop("increments are negative but log is FALSE!")
  if (log && !any(increments < 0)) warning("log is TRUE but no increments are negative; did you pass log-increments?")
  if (length(initial) != 1L || !is.finite(initial) || initial <= 0)
    stop("`initial` must be a positive finite scalar.", call. = FALSE)
  path <- if (log) log(initial) + cumsum(increments) else initial * cumprod(pmax(increments, .Machine$double.eps))
  if (isTRUE(running_max)) return(cummax(path))
  path
}

# Method: compute_tsm for TSM
# inputs:
#   object = TSM object
#   x      = numeric vector or matrix sequence
#   log    = logical; if TRUE return log-TSM path
# outputs:
#   numeric length-N vector with cumulative TSM path (or log-TSM path)
setMethod("compute_tsm", "TSM", function(object, x, log = FALSE) {
  inc <- compute_increments(object, x, log = log)
  increments_to_tsm(inc, initial = object@initial, running_max = FALSE, log = log)
})
