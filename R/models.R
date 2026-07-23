# Class: Model
# purpose: base class for pre/post change models used in likelihood increments
setClass("Model", slots = c(name = "character"))

# Class: UnivariateModel
# purpose: base class for models with scalar observations at each time t
setClass("UnivariateModel", contains = "Model")

# Class: MultivariateModel
# purpose: base class for models with vector observations at each time t
setClass("MultivariateModel", contains = "Model")

# Class: GaussianVARModel
# purpose: unified Gaussian VAR(p) model covering univariate IID (K=1, Phi=list()),
#          univariate AR(p) (K=1, Phi=list of 1x1 matrices), and vector VAR(p) (K>1).
#          Pre and post regimes share Phi and Sigma by default; the primary use-case
#          is detecting a shift in the long-run mean vector.
# slots:
#   Phi_pre, Phi_post   = list of p K-by-K AR coefficient matrices (empty list = IID)
#   Sigma_pre, Sigma_post = K-by-K innovation covariance matrices
#   mean_pre, mean_post   = length-K long-run mean vectors
#   x0                    = length-K initialisation vector used when history is too short
setClass(
  "GaussianVARModel",
  contains = "MultivariateModel",
  slots = c(
    Phi_pre   = "list",   Sigma_pre  = "matrix", mean_pre  = "numeric",
    Phi_post  = "list",   Sigma_post = "matrix", mean_post = "numeric",
    x0        = "numeric"
  )
)

# Class: BernoulliModel
# purpose: Bernoulli pre/post model
setClass("BernoulliModel", contains = "UnivariateModel",
         slots = c(p_pre = "numeric", p_post = "numeric"))

# Class: BoundedModel
# purpose: test martingale model for [0,1]-bounded univariate observations.
#          The one-sided null is E[X_t] <= eta_t; alternative is E[X_t] > eta_t.
#          Increment at time t: 1 + lambda_t * (X_t - eta_t), where lambda_t is
#          determined by the supplied BettingStrategy (default: AGRAPABet).
# slots:
#   eta  = null mean upper bound: scalar in (0,1) (IID null, length-1) or
#          a pre-computed length-N vector (time-varying null, e.g. from ARBoundedModel)
#   bets = BettingStrategy object that maps history -> lambda_t
setClass(
  "BoundedModel",
  contains = "UnivariateModel",
  slots = c(
    eta  = "numeric",
    bets = "BettingStrategy"
  )
)

# ---- internal helpers ----

# helper: VAR(p) stationarity check via companion matrix eigenvalues
# A VAR(p) is stationary iff all eigenvalues of the Kp-by-Kp companion matrix
# lie strictly inside the complex unit circle.
.var_is_stationary <- function(Phi_list) {
  p <- length(Phi_list)
  if (p == 0L) return(TRUE)
  K    <- nrow(Phi_list[[1]])
  top  <- do.call(cbind, Phi_list)
  comp <- if (p > 1L) {
    eye  <- diag(K * (p - 1L))
    zero <- matrix(0, nrow = K * (p - 1L), ncol = K)
    rbind(top, cbind(eye, zero))
  } else {
    top
  }
  all(Mod(eigen(comp, only.values = TRUE)$values) < 1 - 1e-8)
}

# helper: convert a vector or row to a 1-by-k matrix; verify dimension if k supplied
.mv_as_matrix <- function(x, k = NULL) {
  if (is.null(dim(x))) {
    x_vec <- as.numeric(x)
    if (!is.null(k) && length(x_vec) != k)
      stop("Observation dimension does not match model dimension.", call. = FALSE)
    return(matrix(x_vec, nrow = 1L))
  }
  x_mat <- as.matrix(x)
  if (!is.null(k) && ncol(x_mat) != k)
    stop("Observation dimension does not match model dimension.", call. = FALSE)
  x_mat
}

# helper: multivariate Gaussian log-density; returns a scalar for vector x
.mv_log_density <- function(x, mean, Sigma) {
  x_mat <- .mv_as_matrix(x, k = length(mean))
  out   <- mvtnorm::dmvnorm(x_mat, mean = mean, sigma = Sigma, log = TRUE)
  if (is.null(dim(x))) out[1L] else out
}

# helper: normalise the history argument to a 0-row or n-row K-column matrix
.gvar_normalize_history <- function(history, K) {
  if (is.null(history) || length(history) == 0L)
    return(matrix(numeric(0), nrow = 0L, ncol = K))
  if (is.null(dim(history)))
    return(matrix(as.numeric(history), ncol = K))
  as.matrix(history)
}

# helper: VAR(p) conditional mean given Phi list, long-run mean, x0, and history matrix.
# Pads with x0 rows when history is shorter than the lag order.
.gvar_cond_mean <- function(Phi, mean_vec, x0, hist_mat) {
  p <- length(Phi)
  K <- length(mean_vec)
  if (p == 0L) return(mean_vec)
  I_K <- if (K == 1L) matrix(1) else diag(K)
  nu  <- as.numeric((I_K - Reduce("+", Phi)) %*% mean_vec)
  n_h <- nrow(hist_mat)
  if (n_h >= p) {
    lag_mat <- hist_mat[(n_h - p + 1L):n_h, , drop = FALSE]
    lag_mat <- lag_mat[seq(p, 1L), , drop = FALSE]   # row 1 = lag-1, row p = lag-p
  } else {
    pad <- matrix(rep(x0, p - n_h), nrow = p - n_h, ncol = K, byrow = TRUE)
    lag_mat <- if (n_h > 0L) rbind(hist_mat[seq(n_h, 1L), , drop = FALSE], pad) else pad
  }
  nu + Reduce("+", lapply(seq_len(p), function(l) as.numeric(Phi[[l]] %*% lag_mat[l, ])))
}

# ---- constructors ----

# Constructor: GaussianVARModel
# inputs:
#   Phi_pre, Phi_post     = list of p K-by-K AR coefficient matrices (lag 1 ... lag p);
#                           a bare K-by-K matrix is coerced to a length-1 list;
#                           a plain numeric vector is treated as K=1 AR coefficients and
#                           each element is wrapped in a 1-by-1 matrix;
#                           use list() for IID (VAR(0));
#                           Phi_post defaults to Phi_pre
#   Sigma_pre, Sigma_post = K-by-K innovation covariance matrices;
#                           Sigma_post defaults to Sigma_pre
#   mean_pre, mean_post   = length-K long-run mean vectors
#   x0                    = length-K initialisation vector (default zero vector)
#   name                  = model label
# outputs:
#   GaussianVARModel object
GaussianVARModel <- function(Phi_pre, Sigma_pre, mean_pre,
                        mean_post,
                        Phi_post   = Phi_pre,
                        Sigma_post = Sigma_pre,
                        x0         = NULL,
                        name       = "gaussian-var") {
  # Coerce univariate AR coefficient vector to list of 1-by-1 matrices
  if (is.numeric(Phi_pre)  && !is.matrix(Phi_pre))
    Phi_pre  <- lapply(as.numeric(Phi_pre),  matrix)
  if (is.numeric(Phi_post) && !is.matrix(Phi_post))
    Phi_post <- lapply(as.numeric(Phi_post), matrix)
  # Coerce bare matrix to single-element list
  if (is.matrix(Phi_pre))  Phi_pre  <- list(Phi_pre)
  if (is.matrix(Phi_post)) Phi_post <- list(Phi_post)

  mean_pre  <- as.numeric(mean_pre)
  mean_post <- as.numeric(mean_post)
  K         <- length(mean_pre)

  if (length(mean_post) != K)
    stop("`mean_pre` and `mean_post` must have the same length.", call. = FALSE)
  if (!is.matrix(Sigma_pre)  || nrow(Sigma_pre)  != K || ncol(Sigma_pre)  != K)
    stop("`Sigma_pre` must be a K-by-K matrix.", call. = FALSE)
  if (!is.matrix(Sigma_post) || nrow(Sigma_post) != K || ncol(Sigma_post) != K)
    stop("`Sigma_post` must be a K-by-K matrix.", call. = FALSE)

  for (i in seq_along(Phi_pre))
    if (!is.matrix(Phi_pre[[i]]) || nrow(Phi_pre[[i]]) != K || ncol(Phi_pre[[i]]) != K)
      stop("Each element of `Phi_pre` must be a K-by-K matrix.", call. = FALSE)
  for (i in seq_along(Phi_post))
    if (!is.matrix(Phi_post[[i]]) || nrow(Phi_post[[i]]) != K || ncol(Phi_post[[i]]) != K)
      stop("Each element of `Phi_post` must be a K-by-K matrix.", call. = FALSE)

  if (length(Phi_pre)  > 0L && !.var_is_stationary(Phi_pre))
    stop("Pre-change VAR model is not stationary: all companion-matrix eigenvalues must lie strictly inside the unit circle.", call. = FALSE)
  if (length(Phi_post) > 0L && !.var_is_stationary(Phi_post))
    stop("Post-change VAR model is not stationary: all companion-matrix eigenvalues must lie strictly inside the unit circle.", call. = FALSE)

  if (is.null(x0)) x0 <- rep(0, K)
  if (length(x0) != K)
    stop("`x0` must have length K.", call. = FALSE)

  new("GaussianVARModel",
      name       = name,
      Phi_pre    = Phi_pre,
      Sigma_pre  = Sigma_pre,
      mean_pre   = mean_pre,
      Phi_post   = Phi_post,
      Sigma_post = Sigma_post,
      mean_post  = mean_post,
      x0         = as.numeric(x0))
}

# Constructor: GaussianModel
# purpose: convenience wrapper — K=1 IID Gaussian; returns a GaussianVARModel with Phi = list()
# inputs:
#   mean_pre, mean_post = scalar long-run means
#   sd_pre, sd_post     = positive innovation standard deviations (default 1)
#   name                = model label
# outputs:
#   GaussianVARModel object (K = 1, VAR(0))
GaussianModel <- function(mean_pre, mean_post, sd_pre = 1, sd_post = 1, name = "gaussian") {
  if (is.null(sd_pre)  || length(sd_pre)  == 0L) sd_pre  <- 1
  if (is.null(sd_post) || length(sd_post) == 0L) sd_post <- 1
  GaussianVARModel(
    Phi_pre    = list(),
    Sigma_pre  = matrix(sd_pre^2),
    mean_pre   = as.numeric(mean_pre),
    Phi_post   = list(),
    Sigma_post = matrix(sd_post^2),
    mean_post  = as.numeric(mean_post),
    name       = name
  )
}

# Constructor: MultivariateGaussianModel
# purpose: convenience wrapper — IID multivariate Gaussian; returns GaussianVARModel with Phi = list()
# inputs:
#   mu_pre, mu_post       = length-K long-run mean vectors
#   Sigma_pre, Sigma_post = K-by-K covariance matrices (default identity)
#   name                  = model label
# outputs:
#   GaussianVARModel object (K >= 1, VAR(0))
MultivariateGaussianModel <- function(mu_pre, Sigma_pre = NULL, mu_post, Sigma_post = NULL,
                                      name = "mv-gaussian") {
  K <- length(mu_pre)
  if (is.null(Sigma_pre))  Sigma_pre  <- diag(K)
  if (is.null(Sigma_post)) Sigma_post <- diag(K)
  GaussianVARModel(
    Phi_pre    = list(),
    Sigma_pre  = Sigma_pre,
    mean_pre   = as.numeric(mu_pre),
    Phi_post   = list(),
    Sigma_post = Sigma_post,
    mean_post  = as.numeric(mu_post),
    name       = name
  )
}

# Constructor: BernoulliModel
# inputs:
#   p_pre, p_post = scalar probabilities in (0, 1)
#   name          = model label
# outputs:
#   BernoulliModel object
BernoulliModel <- function(p_pre, p_post, name = "bernoulli") {
  stopifnot(length(p_pre) == 1L, length(p_post) == 1L,
            p_pre > 0, p_pre < 1, p_post > 0, p_post < 1)
  new("BernoulliModel", name = name, p_pre = p_pre, p_post = p_post)
}

# Constructor: BoundedModel
# purpose: builds a non-negative-data test martingale model with a pluggable betting strategy.
#          The increment 1 + lambda*(X - eta) is non-negative whenever X >= 0 and
#          lambda <= 1/eta, which the AGRAPA cap (c/eta, c < 1) enforces automatically.
#          eta need only be a positive upper bound on E_null[X]; it does not need to be < 1.
# inputs:
#   eta  = null mean upper bound: a positive scalar or a length-N positive numeric vector
#          of time-varying conditional null means.
#          Classic use case: X in [0,1] data with eta in (0,1).
#          Extended use case: X >= 0 data (e.g. exp-residuals) with eta in (0, Inf).
#          Use ARBoundedModel() to construct the eta vector for AR(p) pre-change data.
#   bets = BettingStrategy object (default: AGRAPABet())
#   name = character model label
# outputs:
#   BoundedModel object
BoundedModel <- function(eta, bets = AGRAPABet(), name = "bounded") {
  eta <- as.numeric(eta)
  if (length(eta) == 0L || !all(is.finite(eta)) || any(eta <= 0))
    stop("`eta` must be a non-empty positive finite numeric scalar or vector.", call. = FALSE)
  if (!is(bets, "BettingStrategy"))
    stop("`bets` must be a BettingStrategy object.", call. = FALSE)
  new("BoundedModel", name = name, eta = eta, bets = bets)
}

# ---- model_density methods ----

# Method: model_density for GaussianVARModel
# inputs:
#   object  = GaussianVARModel object
#   x       = length-K numeric vector (single observation; scalar acceptable when K=1)
#   regime  = "pre" or "post"
#   history = optional numeric N_hist-by-K matrix (or length-N vector for K=1) of past
#             observations; rows shorter than the lag order are padded with x0
# outputs:
#   numeric scalar Gaussian conditional density
setMethod("model_density", "GaussianVARModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  Phi    <- if (regime == "pre") object@Phi_pre  else object@Phi_post
  Sigma  <- if (regime == "pre") object@Sigma_pre else object@Sigma_post
  m      <- if (regime == "pre") object@mean_pre  else object@mean_post
  K      <- length(m)

  hist_mat  <- .gvar_normalize_history(history, K)
  cond_mean <- .gvar_cond_mean(Phi, m, object@x0, hist_mat)
  exp(.mv_log_density(as.numeric(x), mean = cond_mean, Sigma = Sigma))
})

# Method: model_density for BernoulliModel
# inputs:
#   object  = BernoulliModel
#   x       = 0/1 observation(s)
#   regime  = "pre" or "post"
#   history = unused (kept for generic compatibility)
# outputs:
#   numeric Bernoulli probability mass value(s)
setMethod("model_density", "BernoulliModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  prob   <- if (regime == "pre") object@p_pre else object@p_post
  stats::dbinom(x, size = 1L, prob = prob)
})

# ---- likelihood_increment methods ----

# Method: likelihood_increment for GaussianVARModel
# inputs:
#   object  = GaussianVARModel object
#   x       = length-K numeric vector (single observation; scalar acceptable when K=1)
#   history = optional history: numeric vector (K=1) or matrix (K>1) of past observations
#   log     = logical; return log-increment when TRUE
# outputs:
#   numeric scalar increment (or log-increment)
setMethod("likelihood_increment", "GaussianVARModel", function(object, x, history = NULL, log = FALSE) {
  K        <- length(object@mean_pre)
  hist_mat <- .gvar_normalize_history(history, K)
  x_vec    <- as.numeric(x)

  log_pre  <- .mv_log_density(x_vec,
    mean  = .gvar_cond_mean(object@Phi_pre,  object@mean_pre,  object@x0, hist_mat),
    Sigma = object@Sigma_pre)
  log_post <- .mv_log_density(x_vec,
    mean  = .gvar_cond_mean(object@Phi_post, object@mean_post, object@x0, hist_mat),
    Sigma = object@Sigma_post)

  log_inc <- log_post - log_pre
  if (log) log_inc else pmax(exp(log_inc), .Machine$double.eps)
})

# Method: likelihood_increment for Model (generic fallback via density ratio)
# inputs:
#   object  = Model subclass with model_density implemented
#   x       = scalar/vector observation(s)
#   history = optional history passed to model_density
#   log     = logical; return log-increment when TRUE
# outputs:
#   numeric increment(s)
setMethod("likelihood_increment", "Model", function(object, x, history = NULL, log = FALSE) {
  pre  <- model_density(object, x = x, regime = "pre",  history = history)
  post <- model_density(object, x = x, regime = "post", history = history)
  if (log) {
    return(log(pmax(post, .Machine$double.eps)) - log(pmax(pre, .Machine$double.eps)))
  }
  pmax(post / pmax(pre, .Machine$double.eps), .Machine$double.eps)
})

# Method: likelihood_increment for BoundedModel
# purpose: betting increment 1 + lambda_t * (x_t - eta), where lambda_t is
#          determined by the model's BettingStrategy.  Missing observations
#          (NA) yield an increment of 1 (no evidence update).
# inputs:
#   object  = BoundedModel object
#   x       = scalar observation in [0, 1] (or NA for missing)
#   history = numeric vector of past observations X_1,...,X_{t-1} (or NULL)
#   log     = logical; return log-increment when TRUE
# outputs:
#   numeric scalar increment (or log-increment)
setMethod("likelihood_increment", "BoundedModel", function(object, x, history = NULL, log = FALSE) {
  if (is.na(x)) return(if (log) 0 else 1)
  # For a time-varying eta vector, infer the current time from history length.
  eta <- if (length(object@eta) == 1L) object@eta
         else object@eta[length(history) + 1L]
  lam <- compute_bet(object@bets, history, eta)
  inc <- max(1 + lam * (as.numeric(x) - eta), .Machine$double.eps)
  if (log) base::log(inc) else inc
})

# Constructor: ARBoundedModel
# purpose: builds a BoundedModel with a time-varying null mean eta_t equal to
#          the AR(p) conditional pre-change mean E[X_t | X_{t-1},...,X_{t-p}].
#          This makes the resulting BoundedModel TSM valid under AR(p) pre-change
#          data: under H0, E[m_t | F_{t-1}] = 1 exactly (true martingale), while
#          under H1 the increment has positive expected value, enabling detection.
#
#          The function computes eta_t from the observed stream x and the AR
#          parameters (phi, mu), then calls BoundedModel(eta = eta_t, bets).
#          Pass the same x to both ARBoundedModel() and compute_increments().
#
# inputs:
#   phi  = numeric vector of AR coefficients (phi_1, ..., phi_p); p = length(phi)
#   mu   = numeric scalar unconditional pre-change mean in (0, 1)
#   x    = numeric vector of observed data (pre- and post-change); length N
#   bets = BettingStrategy (default: EWMABet with mu_init = mu)
#   name = character model label
# outputs:
#   BoundedModel with eta slot set to the length-N vector of conditional null means
ARBoundedModel <- function(phi, mu, x,
                           bets = EWMABet(rho = 0.1, mu_init = mu),
                           name = "AR-bounded") {
  phi <- as.numeric(phi)
  mu  <- as.numeric(mu)
  x   <- as.numeric(x)
  if (!all(is.finite(phi)))
    stop("`phi` must be a finite numeric vector.", call. = FALSE)
  if (length(mu) != 1L || !is.finite(mu) || mu <= 0 || mu >= 1)
    stop("`mu` must be a finite scalar in (0, 1).", call. = FALSE)
  if (length(x) == 0L)
    stop("`x` must be a non-empty numeric vector.", call. = FALSE)

  n <- length(x)
  p <- length(phi)

  # eta_t = mu + sum_{j=1}^{min(p, t-1)} phi_j * (x[t-j] - mu)
  # At t <= p, only the available lags are used; at t = 1, eta_1 = mu.
  eta_t <- numeric(n)
  for (t in seq_len(n)) {
    avail    <- min(p, t - 1L)
    eta_t[t] <- mu + if (avail == 0L) 0
                     else sum(phi[seq_len(avail)] *
                              (x[seq.int(t - 1L, t - avail)] - mu))
  }

  # Check that all conditional null means are in (0, 1), as required by
  # BoundedModel.  Values outside this range indicate a mis-specified AR model
  # (e.g., non-stationary phi, or extreme observations outside [0, 1]).
  if (any(!is.finite(eta_t)))
    stop("computed eta_t contains non-finite values; check phi and x.",
         call. = FALSE)
  if (any(eta_t <= 0) || any(eta_t >= 1))
    stop(
      "computed eta_t values must lie in (0, 1).  ",
      sprintf("%d value(s) outside range (min=%.4g, max=%.4g).  ",
              sum(eta_t <= 0 | eta_t >= 1), min(eta_t), max(eta_t)),
      "Ensure the AR model is stationary and all observations are in (0, 1).",
      call. = FALSE)

  BoundedModel(eta = eta_t, bets = bets, name = name)
}
