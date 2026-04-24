# Class: Model
# purpose: base class for pre/post change models used in likelihood increments
setClass("Model", slots = c(name = "character"))

# Class: UnivariateModel
# purpose: base class for models with scalar observations at each time t
setClass("UnivariateModel", contains = "Model")

# Class: MultivariateModel
# purpose: base class for models with vector observations at each time t
setClass("MultivariateModel", contains = "Model")

# Class: GaussianModel
# purpose: unified univariate Gaussian model supporting simple/composite pre/post means
# slots:
#   mean_pre, mean_post = scalar (simple) or length-2 interval c(a,b) (composite)
#   sd_pre, sd_post     = scalar (>0) or empty (estimate via MLE)
#   method              = "predictable" or "mixture" for composite post-change mean
#   grid_size           = grid points for mixture approximation
#   prior_weights       = optional fixed mixture weights over post-change grid
#   min_sd              = floor for sd estimates
setClass(
  "GaussianModel",
  contains = "UnivariateModel",
  slots = c(
    mean_pre = "ANY",
    mean_post = "ANY",
    sd_pre = "numeric",
    sd_post = "numeric",
    method = "character",
    grid_size = "numeric",
    prior_weights = "numeric",
    min_sd = "numeric"
  )
)

# Class: MultivariateGaussianModel
# purpose: unified multivariate Gaussian model supporting simple/composite pre/post means
# slots:
#   mu_pre, mu_post       = vector length K (simple) or K-by-2 box matrix (composite)
#   Sigma_pre, Sigma_post = K-by-K covariance matrix or NULL (estimate via MLE)
#   method                = "predictable" or "mixture" for composite post-change mean box
#   grid_size             = points per axis for mixture grid approximation
#   prior_weights         = optional fixed mixture weights over mean grid
#   ridge                 = small ridge for numerical PD covariance handling
#   max_grid_points       = cap on mixture grid size
setClass(
  "MultivariateGaussianModel",
  contains = "MultivariateModel",
  slots = c(
    mu_pre = "ANY",
    mu_post = "ANY",
    Sigma_pre = "ANY",
    Sigma_post = "ANY",
    method = "character",
    grid_size = "numeric",
    prior_weights = "numeric",
    ridge = "numeric",
    max_grid_points = "numeric"
  )
)

# Class: BernoulliModel
# purpose: Bernoulli pre/post model
setClass("BernoulliModel", contains = "UnivariateModel", slots = c(p_pre = "numeric", p_post = "numeric"))

# Class: AR1Model
# purpose: Gaussian AR(1) model, parameterized by long-run mean and AR coefficient
setClass(
  "AR1Model",
  contains = "UnivariateModel",
  slots = c(
    phi_pre = "numeric", sigma_pre = "numeric", mu_pre = "numeric",
    phi_post = "numeric", sigma_post = "numeric", mu_post = "numeric", x0 = "numeric"
  )
)

# Class: ARpModel
# purpose: Gaussian AR(p) model with potentially different pre/post orders;
#          parameterized by long-run mean, AR coefficients, and innovation sd
# slots:
#   phi_pre, phi_post     = numeric vectors of AR coefficients (length = order p)
#   sigma_pre, sigma_post = positive innovation sd scalars
#   mean_pre, mean_post   = long-run (unconditional) means
#   x0                    = numeric scalar used as the initial lag when history is too short
setClass(
  "ARpModel",
  contains = "UnivariateModel",
  slots = c(
    phi_pre   = "numeric", sigma_pre  = "numeric", mean_pre  = "numeric",
    phi_post  = "numeric", sigma_post = "numeric", mean_post = "numeric",
    x0        = "numeric"
  )
)

# ---- helper functions ----

# helper: identify scalar numeric model specifications
.is_scalar_spec <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

# helper: identify univariate interval specifications c(a, b), a < b
.is_interval_spec <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x)) && x[1] < x[2]

# helper: identify K-by-2 box constraints for multivariate means
.is_matrix_box_spec <- function(x) {
  is.matrix(x) && is.numeric(x) && ncol(x) == 2L && all(is.finite(x)) && all(x[, 1] < x[, 2])
}

# helper: clamp scalar x into [lo, hi]
.clamp <- function(x, lo, hi) min(max(x, lo), hi)

# helper: check whether an AR coefficient vector defines a stationary process
# Stationarity iff all roots of the characteristic polynomial 1 - phi_1*z - ... - phi_p*z^p
# lie strictly outside the complex unit circle, equivalently all roots of
# z^p - phi_1*z^(p-1) - ... - phi_p have modulus > 1.
.ar_is_stationary <- function(phi) {
  p <- length(phi)
  if (p == 0L) return(TRUE)
  # Companion form: characteristic polynomial coefficients (highest degree first)
  char_poly <- c(1, -phi)
  roots <- polyroot(char_poly)
  all(Mod(roots) > 1 + 1e-8)
}

# helper: convert long-run mean m to AR intercept
# For AR(p): E[X] = mu / (1 - sum(phi))  =>  mu = m * (1 - sum(phi))
.ar_intercept_from_mean <- function(m, phi) {
  m * (1 - sum(phi))
}

# helper: convert AR intercept mu to long-run mean
.ar_mean_from_intercept <- function(mu, phi) {
  mu / (1 - sum(phi))
}

# helper: project vector v into a K-by-2 box
.project_box <- function(v, box) pmin(pmax(v, box[, 1]), box[, 2])

# helper: geometric grid on an interval for mixture Riemann sums
.gaussian_grid <- function(interval, grid_size) {
  seq(interval[1], interval[2], length.out = as.integer(grid_size))
}

# helper: Cartesian product grid on a K-dimensional box for mixture Riemann sums
.mv_box_grid <- function(box, grid_size, max_grid_points) {
  axes <- lapply(seq_len(nrow(box)), function(j) {
    seq(box[j, 1], box[j, 2], length.out = as.integer(grid_size))
  })

  grid <- as.matrix(expand.grid(axes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  if (nrow(grid) > max_grid_points) {
    idx <- unique(round(seq(1, nrow(grid), length.out = max_grid_points)))
    grid <- grid[idx, , drop = FALSE]
  }
  grid
}

# helper: normalize fixed mixture weights; defaults to flat weights
.normalize_weights <- function(weights, n_points) {
  if (length(weights) == 0L) {
    return(rep(1 / n_points, n_points))
  }
  if (length(weights) != n_points) {
    stop("`prior_weights` must be empty or have one entry per grid point.", call. = FALSE)
  }
  if (any(weights < 0) || sum(weights) <= 0) {
    stop("`prior_weights` must be nonnegative and sum to a positive value.", call. = FALSE)
  }
  weights / sum(weights)
}

# helper: convert vector/row inputs to an n-by-k matrix and preserve scalar-row behavior
.mv_as_matrix <- function(x, k = NULL) {
  if (is.null(dim(x))) {
    x_vec <- as.numeric(x)
    if (!is.null(k) && length(x_vec) != k) {
      stop("Observation dimension does not match model dimension.", call. = FALSE)
    }
    return(matrix(x_vec, nrow = 1L))
  }

  x_mat <- as.matrix(x)
  if (!is.null(k) && ncol(x_mat) != k) {
    stop("Observation dimension does not match model dimension.", call. = FALSE)
  }
  x_mat
}

# helper: thin wrapper around mvtnorm::dmvnorm with scalar return for vector input
.mv_log_density <- function(x, mean, Sigma) {
  x_mat <- .mv_as_matrix(x, k = length(mean))
  out <- mvtnorm::dmvnorm(x_mat, mean = mean, sigma = Sigma, log = TRUE)
  if (is.null(dim(x))) out[1L] else out
}

# helper: ensure covariance matrix is symmetric PD by adding ridge if needed
.ensure_pd <- function(S, ridge, fallback) {
  if (!is.matrix(S) || any(!is.finite(S))) {
    return(fallback)
  }

  S <- (S + t(S)) / 2
  p <- nrow(S)
  for (k in 0:8) {
    S_try <- S + diag(ridge * 10^k, p)
    if (!inherits(try(chol(S_try), silent = TRUE), "try-error")) {
      return(S_try)
    }
  }
  fallback
}

# helper: Gaussian mean MLE with optional interval truncation
.gaussian_mean_mle <- function(x, mean_spec) {
  if (.is_scalar_spec(mean_spec)) {
    return(as.numeric(mean_spec))
  }

  interval <- as.numeric(mean_spec)
  if (length(x) == 0L) {
    return(mean(interval))
  }
  .clamp(mean(x), interval[1], interval[2])
}

# helper: Gaussian sd MLE (dividing by n) around provided mean
.gaussian_sd_mle <- function(x, mu, min_sd, fallback_sd = 1) {
  fallback_sd <- max(as.numeric(fallback_sd), min_sd)
  if (length(x) == 0L) {
    return(fallback_sd)
  }

  s <- sqrt(mean((x - mu)^2))
  if (!is.finite(s) || s <= 0) {
    return(fallback_sd)
  }
  max(s, min_sd)
}

# helper: lagged Gaussian sd MLE for post-change unknown sd
.gaussian_post_sd_lagged <- function(history, model) {
  if (length(model@sd_post) == 1L) {
    return(model@sd_post)
  }

  fallback <- if (length(model@sd_pre) == 1L) model@sd_pre else 1
  if (length(history) == 0L) {
    return(max(fallback, model@min_sd))
  }

  mu_hist <- mean(history)
  .gaussian_sd_mle(history, mu = mu_hist, min_sd = model@min_sd, fallback_sd = fallback)
}

# helper: multivariate mean MLE with optional box truncation
.mv_mean_mle <- function(x, mu_spec) {
  if (is.numeric(mu_spec) && is.null(dim(mu_spec))) {
    return(as.numeric(mu_spec))
  }

  box <- mu_spec
  if (nrow(x) == 0L) {
    return(rowMeans(box))
  }
  .project_box(colMeans(x), box)
}

# helper: multivariate covariance MLE (dividing by n) around provided mean
.mv_cov_mle <- function(x, mu, ridge, fallback) {
  if (nrow(x) == 0L) {
    return(fallback)
  }

  xc <- sweep(x, 2, as.numeric(mu), "-")
  S <- crossprod(xc) / nrow(x)
  .ensure_pd(S, ridge = ridge, fallback = fallback)
}

# helper: lagged multivariate covariance MLE for post-change unknown Sigma
.mv_post_cov_lagged <- function(history, model, k) {
  if (!is.null(model@Sigma_post)) {
    return(model@Sigma_post)
  }

  fallback <- if (!is.null(model@Sigma_pre)) model@Sigma_pre else diag(1, k)
  if (nrow(history) == 0L) {
    return(fallback)
  }

  mu_hist <- colMeans(history)
  .mv_cov_mle(history, mu = mu_hist, ridge = model@ridge, fallback = fallback)
}

# helper: check whether two K-by-2 boxes overlap
.boxes_overlap <- function(b1, b2) {
  all(pmax(b1[, 1], b2[, 1]) <= pmin(b1[, 2], b2[, 2]))
}

# helper: cumulative Gaussian log-likelihood path under pre-change model
.gaussian_pre_loglik_path <- function(model, x) {
  n <- length(x)
  if (n == 0L) return(numeric(0))

  simple_pre <- .is_scalar_spec(model@mean_pre) && length(model@sd_pre) == 1L
  if (simple_pre) {
    return(cumsum(stats::dnorm(x, mean = as.numeric(model@mean_pre), sd = model@sd_pre, log = TRUE)))
  }

  out <- numeric(n)
  fallback_sd <- if (length(model@sd_post) == 1L) model@sd_post else 1

  for (t in seq_len(n)) {
    xt <- x[seq_len(t)]
    mu0 <- .gaussian_mean_mle(xt, model@mean_pre)
    sd0 <- if (length(model@sd_pre) == 1L) {
      model@sd_pre
    } else {
      .gaussian_sd_mle(xt, mu = mu0, min_sd = model@min_sd, fallback_sd = fallback_sd)
    }
    out[t] <- sum(stats::dnorm(xt, mean = mu0, sd = sd0, log = TRUE))
  }

  out
}

# helper: cumulative Gaussian log-likelihood path under post-change model
.gaussian_post_loglik_path <- function(model, x) {
  n <- length(x)
  if (n == 0L) return(numeric(0))

  post_is_interval <- .is_interval_spec(model@mean_post)
  is_mixture_post <- post_is_interval && identical(model@method, "mixture")

  if (!is_mixture_post) {
    out <- numeric(n)
    running <- 0

    for (t in seq_len(n)) {
      history <- if (t == 1L) numeric(0) else x[seq_len(t - 1L)]
      mu1 <- .gaussian_mean_mle(history, model@mean_post)
      sd1 <- .gaussian_post_sd_lagged(history, model)
      running <- running + stats::dnorm(x[t], mean = mu1, sd = sd1, log = TRUE)
      out[t] <- running
    }

    return(out)
  }

  grid <- .gaussian_grid(as.numeric(model@mean_post), model@grid_size)
  w <- .normalize_weights(model@prior_weights, length(grid))
  logw <- ifelse(w > 0, log(w), -Inf)

  out <- numeric(n)
  log_acc <- rep(0, length(grid))

  for (t in seq_len(n)) {
    history <- if (t == 1L) numeric(0) else x[seq_len(t - 1L)]
    sd1 <- .gaussian_post_sd_lagged(history, model)
    log_acc <- log_acc + stats::dnorm(x[t], mean = grid, sd = sd1, log = TRUE)
    out[t] <- .logsumexp(logw + log_acc)
  }

  out
}

# helper: cumulative MV Gaussian log-likelihood path under pre-change model
.mv_pre_loglik_path <- function(model, x) {
  n <- nrow(x)
  if (n == 0L) return(numeric(0))

  simple_pre <- (is.numeric(model@mu_pre) && is.null(dim(model@mu_pre))) && !is.null(model@Sigma_pre)
  if (simple_pre) {
    return(cumsum(mvtnorm::dmvnorm(x, mean = as.numeric(model@mu_pre), sigma = model@Sigma_pre, log = TRUE)))
  }

  out <- numeric(n)
  k <- ncol(x)
  fallback <- if (!is.null(model@Sigma_post)) model@Sigma_post else diag(1, k)

  for (t in seq_len(n)) {
    xt <- x[seq_len(t), , drop = FALSE]
    mu0 <- .mv_mean_mle(xt, model@mu_pre)
    Sigma0 <- if (!is.null(model@Sigma_pre)) {
      model@Sigma_pre
    } else {
      .mv_cov_mle(xt, mu = mu0, ridge = model@ridge, fallback = fallback)
    }
    out[t] <- sum(mvtnorm::dmvnorm(xt, mean = mu0, sigma = Sigma0, log = TRUE))
  }

  out
}

# helper: cumulative MV Gaussian log-likelihood path under post-change model
.mv_post_loglik_path <- function(model, x) {
  n <- nrow(x)
  if (n == 0L) return(numeric(0))

  k <- ncol(x)
  post_is_box <- .is_matrix_box_spec(model@mu_post)
  is_mixture_post <- post_is_box && identical(model@method, "mixture")

  if (!is_mixture_post) {
    out <- numeric(n)
    running <- 0

    for (t in seq_len(n)) {
      history <- if (t == 1L) {
        matrix(numeric(0), nrow = 0, ncol = k)
      } else {
        x[seq_len(t - 1L), , drop = FALSE]
      }

      mu1 <- .mv_mean_mle(history, model@mu_post)
      Sigma1 <- .mv_post_cov_lagged(history, model, k)
      running <- running + mvtnorm::dmvnorm(x[t, , drop = FALSE], mean = mu1, sigma = Sigma1, log = TRUE)
      out[t] <- running
    }

    return(out)
  }

  grid <- .mv_box_grid(model@mu_post, model@grid_size, model@max_grid_points)
  g <- nrow(grid)
  w <- .normalize_weights(model@prior_weights, g)
  logw <- ifelse(w > 0, log(w), -Inf)

  out <- numeric(n)
  log_acc <- rep(0, g)

  for (t in seq_len(n)) {
    history <- if (t == 1L) {
      matrix(numeric(0), nrow = 0, ncol = k)
    } else {
      x[seq_len(t - 1L), , drop = FALSE]
    }

    Sigma1 <- .mv_post_cov_lagged(history, model, k)
    log_t <- vapply(seq_len(g), function(i) {
      mvtnorm::dmvnorm(x[t, , drop = FALSE], mean = grid[i, ], sigma = Sigma1, log = TRUE)
    }, numeric(1L))

    log_acc <- log_acc + log_t
    out[t] <- .logsumexp(logw + log_acc)
  }

  out
}

# helper: cumulative Gaussian log-likelihood ratio path
.gaussian_log_lr_cum_path <- function(model, x) {
  .gaussian_post_loglik_path(model, x) - .gaussian_pre_loglik_path(model, x)
}

# helper: cumulative multivariate Gaussian log-likelihood ratio path
.mv_log_lr_cum_path <- function(model, x) {
  .mv_post_loglik_path(model, x) - .mv_pre_loglik_path(model, x)
}

# helper: fast-path wrapper used by TSM for Gaussian mixture-post models
.gaussian_mix_increments_fast <- function(model, x, log = FALSE) {
  if (length(x) == 0L) return(numeric(0))
  log_z <- .gaussian_log_lr_cum_path(model, x)
  log_inc <- c(log_z[1L], diff(log_z))
  if (log) return(log_inc)
  pmax(exp(log_inc), .Machine$double.eps)
}

# helper: fast-path wrapper used by TSM for MV Gaussian mixture-post models
.mv_gaussian_mix_increments_fast <- function(model, x, log = FALSE) {
  if (nrow(x) == 0L) return(numeric(0))
  log_z <- .mv_log_lr_cum_path(model, x)
  log_inc <- c(log_z[1L], diff(log_z))
  if (log) return(log_inc)
  pmax(exp(log_inc), .Machine$double.eps)
}

# ---- constructors ----

# Constructor: GaussianModel
# inputs:
#   mean_pre, mean_post = scalar or interval c(a,b)
#   sd_pre, sd_post     = positive scalar or empty (estimate by MLE)
#   method              = "predictable" or "mixture"
#   grid_size           = integer >= 2, used for mixture over interval means
#   prior_weights       = optional fixed mixture weights on post-change grid
#   min_sd              = positive sd floor for numerical stability
#   name                = model label
# outputs:
#   GaussianModel object
GaussianModel <- function(mean_pre,
                          sd_pre = NULL,
                          mean_post,
                          sd_post = NULL,
                          method = c("predictable", "mixture"),
                          grid_size = 101,
                          prior_weights = numeric(0),
                          min_sd = 1e-8,
                          name = "gaussian") {
  method <- match.arg(method)

  if (!(.is_scalar_spec(mean_pre) || .is_interval_spec(mean_pre))) {
    stop("`mean_pre` must be a scalar or interval c(a,b).", call. = FALSE)
  }
  if (!(.is_scalar_spec(mean_post) || .is_interval_spec(mean_post))) {
    stop("`mean_post` must be a scalar or interval c(a,b).", call. = FALSE)
  }

  if (length(sd_pre) > 1L || (length(sd_pre) == 1L && (!is.finite(sd_pre) || sd_pre <= 0))) {
    stop("`sd_pre` must be empty or a positive scalar.", call. = FALSE)
  }
  if (length(sd_post) > 1L || (length(sd_post) == 1L && (!is.finite(sd_post) || sd_post <= 0))) {
    stop("`sd_post` must be empty or a positive scalar.", call. = FALSE)
  }

  if (length(grid_size) != 1L || grid_size < 2 || grid_size != as.integer(grid_size)) {
    stop("`grid_size` must be integer >= 2.", call. = FALSE)
  }
  if (length(min_sd) != 1L || min_sd <= 0) {
    stop("`min_sd` must be a positive scalar.", call. = FALSE)
  }

  new(
    "GaussianModel",
    name = name,
    mean_pre = mean_pre,
    mean_post = mean_post,
    sd_pre = sd_pre,
    sd_post = sd_post,
    method = method,
    grid_size = as.integer(grid_size),
    prior_weights = as.numeric(prior_weights),
    min_sd = min_sd
  )
}

# Constructor: MultivariateGaussianModel
# inputs:
#   mu_pre, mu_post       = vector (simple) or K-by-2 box matrix (composite)
#   Sigma_pre, Sigma_post = K-by-K covariance matrix or NULL (estimate by MLE)
#   method                = "predictable" or "mixture"
#   grid_size             = integer >= 2 points per axis for mixture grid
#   prior_weights         = optional fixed mixture weights on mean grid
#   ridge                 = positive scalar ridge for covariance regularization
#   max_grid_points       = positive integer cap on mixture grid size
#   name                  = model label
# outputs:
#   MultivariateGaussianModel object
MultivariateGaussianModel <- function(mu_pre,
                                      Sigma_pre = NULL,
                                      mu_post,
                                      Sigma_post = NULL,
                                      method = c("predictable", "mixture"),
                                      grid_size = 5,
                                      prior_weights = numeric(0),
                                      ridge = 1e-6,
                                      max_grid_points = 5000,
                                      name = "mv-gaussian") {
  method <- match.arg(method)

  pre_is_vec <- is.numeric(mu_pre) && is.null(dim(mu_pre))
  pre_is_box <- .is_matrix_box_spec(mu_pre)
  post_is_vec <- is.numeric(mu_post) && is.null(dim(mu_post))
  post_is_box <- .is_matrix_box_spec(mu_post)

  if (!(pre_is_vec || pre_is_box) || !(post_is_vec || post_is_box)) {
    stop("`mu_pre`/`mu_post` must each be a vector (simple) or K-by-2 box matrix (composite).", call. = FALSE)
  }

  K_pre <- if (pre_is_vec) length(mu_pre) else nrow(mu_pre)
  K_post <- if (post_is_vec) length(mu_post) else nrow(mu_post)
  if (K_pre != K_post) {
    stop("Pre/post mean dimensions must match.", call. = FALSE)
  }
  K <- K_pre

  if (pre_is_box && post_is_box && .boxes_overlap(mu_pre, mu_post)) {
    stop("Pre-change and post-change mean boxes must not overlap.", call. = FALSE)
  }

  if (!is.null(Sigma_pre)) {
    if (!is.matrix(Sigma_pre) || nrow(Sigma_pre) != K || ncol(Sigma_pre) != K) {
      stop("`Sigma_pre` must be K-by-K.", call. = FALSE)
    }
    if (inherits(try(chol((Sigma_pre + t(Sigma_pre)) / 2), silent = TRUE), "try-error")) {
      stop("`Sigma_pre` must be positive definite.", call. = FALSE)
    }
  }

  if (!is.null(Sigma_post)) {
    if (!is.matrix(Sigma_post) || nrow(Sigma_post) != K || ncol(Sigma_post) != K) {
      stop("`Sigma_post` must be K-by-K.", call. = FALSE)
    }
    if (inherits(try(chol((Sigma_post + t(Sigma_post)) / 2), silent = TRUE), "try-error")) {
      stop("`Sigma_post` must be positive definite.", call. = FALSE)
    }
  }

  if (length(grid_size) != 1L || grid_size < 2 || grid_size != as.integer(grid_size)) {
    stop("`grid_size` must be integer >= 2.", call. = FALSE)
  }
  if (length(ridge) != 1L || ridge <= 0) {
    stop("`ridge` must be a positive scalar.", call. = FALSE)
  }
  if (length(max_grid_points) != 1L || max_grid_points < 1 || max_grid_points != as.integer(max_grid_points)) {
    stop("`max_grid_points` must be a positive integer.", call. = FALSE)
  }

  new(
    "MultivariateGaussianModel",
    name = name,
    mu_pre = mu_pre,
    mu_post = mu_post,
    Sigma_pre = Sigma_pre,
    Sigma_post = Sigma_post,
    method = method,
    grid_size = as.integer(grid_size),
    prior_weights = as.numeric(prior_weights),
    ridge = ridge,
    max_grid_points = as.integer(max_grid_points)
  )
}

# Constructor: BernoulliModel
# inputs:
#   p_pre, p_post = scalar probabilities in (0,1)
#   name          = model label
# outputs:
#   BernoulliModel object
BernoulliModel <- function(p_pre, p_post, name = "bernoulli") {
  stopifnot(length(p_pre) == 1L, length(p_post) == 1L, p_pre > 0, p_pre < 1, p_post > 0, p_post < 1)
  new("BernoulliModel", name = name, p_pre = p_pre, p_post = p_post)
}

# Constructor: AR1Model
# inputs:
#   phi_pre, phi_post     = AR(1) coefficients in (-1, 1)
#   sigma_pre, sigma_post = positive innovation sd
#   mean_pre, mean_post   = long-run (unconditional) means of the pre/post processes
#   x0                    = initial value when no history is available
#   name                  = model label
# outputs:
#   AR1Model object
AR1Model <- function(phi_pre, sigma_pre, mean_pre = 0, phi_post, sigma_post, mean_post = 0, x0 = 0, name = "ar1") {
  stopifnot(
    length(phi_pre)   == 1L, length(phi_post)  == 1L,
    abs(phi_pre) < 1,        abs(phi_post) < 1,
    sigma_pre > 0,           sigma_post > 0,
    length(mean_pre)  == 1L, length(mean_post) == 1L,
    length(x0) == 1L
  )

  new(
    "AR1Model",
    name       = name,
    phi_pre    = phi_pre,
    sigma_pre  = sigma_pre,
    mu_pre     = .ar_intercept_from_mean(mean_pre, phi_pre),
    phi_post   = phi_post,
    sigma_post = sigma_post,
    mu_post    = .ar_intercept_from_mean(mean_post, phi_post),
    x0         = x0
  )
}

# Constructor: ARpModel
# inputs:
#   phi_pre, phi_post     = numeric vectors of AR coefficients (order p = length of vector)
#   sigma_pre, sigma_post = positive innovation sd scalars
#   mean_pre, mean_post   = long-run (unconditional) means
#   x0                    = scalar used as the initial lag when history is too short
#   name                  = model label
# outputs:
#   ARpModel object
ARpModel <- function(phi_pre, sigma_pre, mean_pre = 0,
                     phi_post, sigma_post, mean_post = 0,
                     x0 = 0, name = "arp") {
  phi_pre  <- as.numeric(phi_pre)
  phi_post <- as.numeric(phi_post)
  stopifnot(
    length(phi_pre) >= 1L, length(phi_post) >= 1L,
    sigma_pre > 0, sigma_post > 0,
    length(mean_pre) == 1L, length(mean_post) == 1L,
    length(x0) == 1L
  )
  if (!.ar_is_stationary(phi_pre)) {
    stop("Pre-change AR model is not stationary: all roots of the characteristic polynomial must lie outside the unit circle.", call. = FALSE)
  }
  if (!.ar_is_stationary(phi_post)) {
    stop("Post-change AR model is not stationary: all roots of the characteristic polynomial must lie outside the unit circle.", call. = FALSE)
  }

  new(
    "ARpModel",
    name       = name,
    phi_pre    = phi_pre,
    sigma_pre  = sigma_pre,
    mean_pre   = as.numeric(mean_pre),
    phi_post   = phi_post,
    sigma_post = sigma_post,
    mean_post  = as.numeric(mean_post),
    x0         = as.numeric(x0)
  )
}

# ---- density methods ----

# Method: model_density for GaussianModel
# inputs:
#   object  = GaussianModel
#   x       = scalar/vector observation(s)
#   regime  = "pre" or "post"
#   history = unused (kept for generic compatibility)
# outputs:
#   numeric density value(s)
setMethod("model_density", "GaussianModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)

  if (regime == "pre") {
    mu <- if (.is_scalar_spec(object@mean_pre)) as.numeric(object@mean_pre) else mean(as.numeric(object@mean_pre))
    sd <- if (length(object@sd_pre) == 1L) object@sd_pre else 1
    return(stats::dnorm(x, mean = mu, sd = sd))
  }

  mu <- if (.is_scalar_spec(object@mean_post)) as.numeric(object@mean_post) else mean(as.numeric(object@mean_post))
  sd <- if (length(object@sd_post) == 1L) object@sd_post else 1
  stats::dnorm(x, mean = mu, sd = sd)
})

# Method: model_density for MultivariateGaussianModel
# inputs:
#   object  = MultivariateGaussianModel
#   x       = numeric vector/matrix observation(s)
#   regime  = "pre" or "post"
#   history = unused (kept for generic compatibility)
# outputs:
#   numeric density value(s)
setMethod("model_density", "MultivariateGaussianModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)

  if (regime == "pre") {
    mu <- if (is.numeric(object@mu_pre) && is.null(dim(object@mu_pre))) as.numeric(object@mu_pre) else rowMeans(object@mu_pre)
    Sigma <- if (!is.null(object@Sigma_pre)) object@Sigma_pre else diag(1, length(mu))
    return(exp(.mv_log_density(x, mean = mu, Sigma = Sigma)))
  }

  mu <- if (is.numeric(object@mu_post) && is.null(dim(object@mu_post))) as.numeric(object@mu_post) else rowMeans(object@mu_post)
  Sigma <- if (!is.null(object@Sigma_post)) object@Sigma_post else diag(1, length(mu))
  exp(.mv_log_density(x, mean = mu, Sigma = Sigma))
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
  prob <- if (regime == "pre") object@p_pre else object@p_post
  stats::dbinom(x, size = 1L, prob = prob)
})

# Method: model_density for AR1Model
# inputs:
#   object  = AR1Model
#   x       = scalar/vector observation(s)
#   regime  = "pre" or "post"
#   history = optional numeric history; last element is used as X_{t-1}
# outputs:
#   numeric Gaussian conditional density value(s)
setMethod("model_density", "AR1Model", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  prev <- object@x0
  if (!is.null(history) && length(history) >= 1L) {
    prev <- history[length(history)]
  }

  if (regime == "pre") {
    return(stats::dnorm(x, mean = object@mu_pre + object@phi_pre * prev, sd = object@sigma_pre))
  }
  stats::dnorm(x, mean = object@mu_post + object@phi_post * prev, sd = object@sigma_post)
})

# Method: model_density for ARpModel
# inputs:
#   object  = ARpModel
#   x       = scalar observation
#   regime  = "pre" or "post"
#   history = numeric vector of past observations (most recent last); if shorter
#             than order p, earlier lags are filled with x0
# outputs:
#   numeric Gaussian conditional density value
setMethod("model_density", "ARpModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)

  phi    <- if (regime == "pre") object@phi_pre  else object@phi_post
  sigma  <- if (regime == "pre") object@sigma_pre else object@sigma_post
  m      <- if (regime == "pre") object@mean_pre  else object@mean_post
  mu     <- .ar_intercept_from_mean(m, phi)   # conditional intercept
  p      <- length(phi)

  # Build lag vector of length p, padding with x0 when history is too short
  hist <- if (is.null(history)) numeric(0) else as.numeric(history)
  n_hist <- length(hist)
  lags <- if (n_hist >= p) {
    rev(hist[(n_hist - p + 1L):n_hist])   # [X_{t-1}, X_{t-2}, ..., X_{t-p}]
  } else {
    c(rev(hist), rep(object@x0, p - n_hist))
  }

  cond_mean <- mu + sum(phi * lags)
  stats::dnorm(x, mean = cond_mean, sd = sigma)
})

# ---- likelihood increment methods ----

# Method: likelihood_increment for GaussianModel
# inputs:
#   object  = GaussianModel
#   x       = scalar/vector observation(s)
#   history = optional numeric history observed before x
#   log     = logical; return log-increment when TRUE
# outputs:
#   numeric increment(s)
setMethod("likelihood_increment", "GaussianModel", function(object, x, history = NULL, log = FALSE) {
  simple_fast <- .is_scalar_spec(object@mean_pre) && .is_scalar_spec(object@mean_post) &&
    length(object@sd_pre) == 1L && length(object@sd_post) == 1L

  if (simple_fast) {
    out_log <- stats::dnorm(x, mean = as.numeric(object@mean_post), sd = object@sd_post, log = TRUE) -
      stats::dnorm(x, mean = as.numeric(object@mean_pre), sd = object@sd_pre, log = TRUE)
    return(if (log) out_log else pmax(exp(out_log), .Machine$double.eps))
  }

  history <- if (is.null(history)) numeric(0) else as.numeric(history)
  x_new <- as.numeric(x)
  xx <- c(history, x_new)

  if (length(xx) == 0L) {
    return(numeric(0))
  }

  log_z <- .gaussian_log_lr_cum_path(object, xx)

  if (length(history) == 0L) {
    log_inc <- c(log_z[1L], diff(log_z))
  } else {
    seg <- log_z[(length(history) + 1L):length(log_z)]
    base <- log_z[length(history)]
    log_inc <- c(seg[1L] - base, diff(seg))
  }

  if (length(x_new) == 1L) {
    log_inc <- log_inc[1L]
  }

  if (log) log_inc else pmax(exp(log_inc), .Machine$double.eps)
})

# Method: likelihood_increment for MultivariateGaussianModel
# inputs:
#   object  = MultivariateGaussianModel
#   x       = numeric vector (single observation) or matrix (multiple observations)
#   history = optional matrix history observed before x
#   log     = logical; return log-increment when TRUE
# outputs:
#   numeric increment(s)
setMethod("likelihood_increment", "MultivariateGaussianModel", function(object, x, history = NULL, log = FALSE) {
  simple_fast <- (is.numeric(object@mu_pre) && is.null(dim(object@mu_pre))) &&
    (is.numeric(object@mu_post) && is.null(dim(object@mu_post))) &&
    !is.null(object@Sigma_pre) && !is.null(object@Sigma_post)

  if (simple_fast) {
    out_log <- .mv_log_density(x, mean = as.numeric(object@mu_post), Sigma = object@Sigma_post) -
      .mv_log_density(x, mean = as.numeric(object@mu_pre), Sigma = object@Sigma_pre)
    return(if (log) out_log else pmax(exp(out_log), .Machine$double.eps))
  }

  k <- if (is.numeric(object@mu_pre) && is.null(dim(object@mu_pre))) {
    length(object@mu_pre)
  } else {
    nrow(object@mu_pre)
  }

  x_new <- .mv_as_matrix(x, k = k)
  hist_mat <- if (is.null(history) || length(history) == 0L) {
    matrix(numeric(0), nrow = 0, ncol = k)
  } else {
    .mv_as_matrix(history, k = k)
  }

  xx <- rbind(hist_mat, x_new)
  if (nrow(xx) == 0L) {
    return(numeric(0))
  }

  log_z <- .mv_log_lr_cum_path(object, xx)

  if (nrow(hist_mat) == 0L) {
    log_inc <- c(log_z[1L], diff(log_z))
  } else {
    seg <- log_z[(nrow(hist_mat) + 1L):length(log_z)]
    base <- log_z[nrow(hist_mat)]
    log_inc <- c(seg[1L] - base, diff(seg))
  }

  if (nrow(x_new) == 1L) {
    log_inc <- log_inc[1L]
  }

  if (log) log_inc else pmax(exp(log_inc), .Machine$double.eps)
})

# Method: likelihood_increment for generic Model fallback
# inputs:
#   object  = Model subclass with model_density implemented
#   x       = scalar/vector observation(s)
#   history = optional history passed to model_density
#   log     = logical; return log-increment when TRUE
# outputs:
#   numeric increment(s)
setMethod("likelihood_increment", "Model", function(object, x, history = NULL, log = FALSE) {
  pre <- model_density(object, x = x, regime = "pre", history = history)
  post <- model_density(object, x = x, regime = "post", history = history)

  if (log) {
    return(log(pmax(post, .Machine$double.eps)) - log(pmax(pre, .Machine$double.eps)))
  }

  pmax(post / pmax(pre, .Machine$double.eps), .Machine$double.eps)
})
