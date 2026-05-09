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
# purpose: test martingale model for [0,1]-bounded univariate observations using
#          AGRAPA bets (Waudby-Smith and Ramdas, 2024).  The one-sided null is
#          E[X_t] <= eta; post-change alternative is E[X_t] > eta.
#          Increment at time t: 1 + lambda_t * (X_t - eta) where lambda_t is
#          the AGRAPA bet computed from the lagged history X_1,...,X_{t-1}.
# slots:
#   eta    = null mean upper bound in (0, 1)
#   c      = AGRAPA truncation factor; bet is capped at c / eta
#   sd_min = floor on the lagged SD estimate (prevents division by zero)
#   eps    = small constant added to eta in the cap denominator
setClass(
  "BoundedModel",
  contains = "UnivariateModel",
  slots = c(
    eta    = "numeric",
    c      = "numeric",
    sd_min = "numeric",
    eps    = "numeric"
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
# purpose: builds a [0,1]-bounded test martingale model with AGRAPA bets
# inputs:
#   eta    = numeric scalar null mean upper bound in (0, 1)
#   c      = numeric scalar AGRAPA truncation factor in (0, 1], default 0.75
#   sd_min = numeric scalar floor for the lagged SD estimate, default 0.01
#   eps    = numeric scalar added to eta in the cap denominator, default 1e-5
#   name   = character model label
# outputs:
#   BoundedModel object
BoundedModel <- function(eta, c = 0.75, sd_min = 0.01, eps = 1e-5, name = "bounded") {
  if (length(eta) != 1L || !is.finite(eta) || eta <= 0 || eta >= 1)
    stop("`eta` must be a scalar in (0, 1).", call. = FALSE)
  if (length(c) != 1L || !is.finite(c) || c <= 0 || c > 1)
    stop("`c` must be a scalar in (0, 1].", call. = FALSE)
  if (length(sd_min) != 1L || sd_min <= 0)
    stop("`sd_min` must be a positive scalar.", call. = FALSE)
  if (length(eps) != 1L || eps < 0)
    stop("`eps` must be a non-negative scalar.", call. = FALSE)
  new("BoundedModel", name = name,
      eta = as.numeric(eta), c = as.numeric(c),
      sd_min = as.numeric(sd_min), eps = as.numeric(eps))
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
# purpose: AGRAPA betting increment — 1 + lambda_t * (x_t - eta) where lambda_t
#          is computed from the lagged history using the AGRAPA formula.
#          When history is empty, lambda_t = 0 and the increment is 1.
# inputs:
#   object  = BoundedModel object
#   x       = scalar observation in [0, 1]
#   history = numeric vector of past observations X_1, ..., X_{t-1} (or NULL)
#   log     = logical; return log-increment when TRUE
# outputs:
#   numeric scalar increment (or log-increment)
setMethod("likelihood_increment", "BoundedModel", function(object, x, history = NULL, log = FALSE) {
  eta <- object@eta
  n_h <- if (is.null(history)) 0L else length(as.numeric(history))
  if (n_h == 0L) {
    lam <- 0
  } else {
    hist_vec <- as.numeric(history)
    lag_mu   <- mean(hist_vec)
    lag_sd   <- if (n_h > 1L) sqrt(mean((hist_vec - lag_mu)^2)) else 0
    lag_sd   <- max(lag_sd, object@sd_min)
    lam_raw  <- (lag_mu - eta) / (lag_sd^2 + (lag_mu - eta)^2)
    cap      <- object@c / (eta + object@eps)
    lam      <- max(0, min(lam_raw, cap))
  }
  inc <- max(1 + lam * (as.numeric(x) - eta), .Machine$double.eps)
  if (log) base::log(inc) else inc
})
