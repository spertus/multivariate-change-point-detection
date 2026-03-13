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
# purpose: unified univariate Gaussian model supporting simple/composite pre/post means and optional unknown sd
# slots:
#   mean_pre, mean_post = scalar (simple) or length-2 interval c(a,b) (composite)
#   sd_pre, sd_post     = scalar (>0) or empty (estimate from data)
#   method              = for composite post mean: "predictable" or "mixture"
#   update_window       = update frequency for predictable estimators
#   grid_size           = grid points for mixture approximation
#   prior_weights       = optional weights for post-mean grid (empty => uniform)
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
    update_window = "numeric",
    grid_size = "numeric",
    prior_weights = "numeric",
    min_sd = "numeric"
  )
)

# Class: MultivariateGaussianModel
# purpose: unified multivariate Gaussian model supporting simple/composite pre/post means and optional unknown covariance
# slots:
#   mu_pre, mu_post     = vector length K (simple) or K-by-2 box matrix (composite)
#   Sigma_pre, Sigma_post = K-by-K covariance matrix or NULL (estimate from data)
#   method              = for composite post mean box: "predictable" or "mixture"
#   update_window       = update frequency for predictable estimators
#   grid_size           = points per axis for mixture grid approximation
#   prior_weights       = optional weights for mean grid (empty => uniform)
#   ridge               = ridge added to covariance estimates for numerical stability
#   max_grid_points     = cap on mixture grid size
setClass(
  "MultivariateGaussianModel",
  contains = "MultivariateModel",
  slots = c(
    mu_pre = "ANY",
    mu_post = "ANY",
    Sigma_pre = "ANY",
    Sigma_post = "ANY",
    method = "character",
    update_window = "numeric",
    grid_size = "numeric",
    prior_weights = "numeric",
    ridge = "numeric",
    max_grid_points = "numeric"
  )
)

# Class: BernoulliModel
setClass("BernoulliModel", contains = "UnivariateModel", slots = c(p_pre = "numeric", p_post = "numeric"))

# Class: AR1Model
setClass(
  "AR1Model",
  contains = "UnivariateModel",
  slots = c(phi_pre = "numeric", sigma_pre = "numeric", mu_pre = "numeric",
            phi_post = "numeric", sigma_post = "numeric", mu_post = "numeric", x0 = "numeric")
)

# ---- helpers ----

.is_scalar_spec <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
.is_interval_spec <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x)) && x[1] < x[2]
.is_matrix_box_spec <- function(x) is.matrix(x) && is.numeric(x) && ncol(x) == 2L && all(is.finite(x)) && all(x[, 1] < x[, 2])

.clamp <- function(x, lo, hi) min(max(x, lo), hi) # clamp a scalar in an interval
.project_box <- function(v, box) pmin(pmax(v, box[, 1]), box[, 2]) # clamp a vector in a box

.logsumexp <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.gaussian_grid <- function(interval, grid_size) seq(interval[1], interval[2], length.out = as.integer(grid_size))

.mv_box_grid <- function(box, grid_size, max_grid_points) {
  axes <- lapply(seq_len(nrow(box)), function(j) seq(box[j, 1], box[j, 2], length.out = as.integer(grid_size)))
  g <- as.matrix(expand.grid(axes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  if (nrow(g) > max_grid_points) {
    idx <- unique(round(seq(1, nrow(g), length.out = max_grid_points)))
    g <- g[idx, , drop = FALSE]
  }
  g
}

.gaussian_sd_mle <- function(x, mu, min_sd = 1e-8) {
  if (length(x) == 0L) return(min_sd)
  s <- sqrt(mean((x - mu)^2))
  if (!is.finite(s) || s <= 0) min_sd else max(s, min_sd)
}

.predictable_mean_estimate <- function(history, update_window, lower, upper, init_mean) {
  t <- length(history) + 1L
  n_used <- ((t - 1L) %/% update_window) * update_window
  if (n_used <= 0L) return(.clamp(init_mean, lower, upper))
  .clamp(mean(history[seq_len(n_used)]), lower, upper)
}

.lagged_sd_seq <- function(x, update_window, fallback_sd, min_sd) {
  n <- length(x)
  out <- rep(fallback_sd, n)
  for (i in seq_len(n)) {
    n_used <- ((i - 1L) %/% update_window) * update_window
    if (n_used <= 1L) {
      out[i] <- max(fallback_sd, min_sd)
    } else {
      s <- stats::sd(x[seq_len(n_used)])
      out[i] <- if (is.finite(s) && s > 0) max(s, min_sd) else max(fallback_sd, min_sd)
    }
  }
  out
}

.mvnorm_log_density <- function(x, mean, Sigma) {
  x_mat <- if (is.null(dim(x))) matrix(as.numeric(x), nrow = 1L) else as.matrix(x)
  p <- length(mean)
  if (ncol(x_mat) != p) stop("Observation dimension does not match mean/covariance dimension.", call. = FALSE)

  R <- chol(Sigma)
  log_det <- 2 * sum(log(diag(R)))
  c0 <- -0.5 * (p * log(2 * pi) + log_det)

  out <- numeric(nrow(x_mat))
  for (i in seq_len(nrow(x_mat))) {
    d <- x_mat[i, ] - mean
    z <- forwardsolve(t(R), d)
    out[i] <- c0 - 0.5 * sum(z^2)
  }
  if (is.null(dim(x))) out[1L] else out
}

.mvnorm_density <- function(x, mean, Sigma) exp(.mvnorm_log_density(x, mean, Sigma))

.mv_cov_mle <- function(x, mu, ridge = 1e-6) {
  p <- length(mu)
  if (is.null(dim(x))) x <- matrix(x, ncol = p)
  if (nrow(x) <= 1L) return(diag(ridge, p))

  xc <- x - matrix(rep(mu, each = nrow(x)), nrow = nrow(x))
  S <- crossprod(xc) / nrow(xc)
  S <- (S + t(S)) / 2
  for (k in 0:6) {
    S_try <- S + diag(ridge * 10^k, p)
    if (!inherits(try(chol(S_try), silent = TRUE), "try-error")) return(S_try)
  }
  diag(ridge, p)
}

.lagged_cov_seq <- function(x, update_window, fallback, ridge) {
  n <- nrow(x)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    n_used <- ((i - 1L) %/% update_window) * update_window
    if (n_used <= 1L) out[[i]] <- fallback else out[[i]] <- .mv_cov_mle(x[seq_len(n_used), , drop = FALSE], colMeans(x[seq_len(n_used), , drop = FALSE]), ridge)
  }
  out
}

.boxes_overlap <- function(b1, b2) {
  # overlap is non-empty if all coordinates overlap
  all(pmax(b1[, 1], b2[, 1]) <= pmin(b1[, 2], b2[, 2]))
}

# ---- constructors ----

GaussianModel <- function(mean_pre,
                          sd_pre = NULL,
                          mean_post,
                          sd_post = NULL,
                          method = c("predictable", "mixture"),
                          update_window = 20,
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

  if (length(update_window) != 1L || update_window < 1 || update_window != as.integer(update_window)) {
    stop("`update_window` must be a positive integer.", call. = FALSE)
  }
  if (length(grid_size) != 1L || grid_size < 2 || grid_size != as.integer(grid_size)) {
    stop("`grid_size` must be integer >= 2.", call. = FALSE)
  }
  if (length(min_sd) != 1L || min_sd <= 0) {
    stop("`min_sd` must be a positive scalar.", call. = FALSE)
  }

  new("GaussianModel",
      name = name,
      mean_pre = mean_pre,
      mean_post = mean_post,
      sd_pre = sd_pre,
      sd_post = sd_post,
      method = method,
      update_window = as.integer(update_window),
      grid_size = as.integer(grid_size),
      prior_weights = as.numeric(prior_weights),
      min_sd = min_sd)
}

MultivariateGaussianModel <- function(mu_pre,
                                      Sigma_pre = NULL,
                                      mu_post,
                                      Sigma_post = NULL,
                                      method = c("predictable", "mixture"),
                                      update_window = 20,
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
  if (K_pre != K_post) stop("Pre/post mean dimensions must match.", call. = FALSE)
  K <- K_pre

  if (pre_is_box && post_is_box && .boxes_overlap(mu_pre, mu_post)) {
    stop("Pre-change and post-change mean boxes must not overlap.", call. = FALSE)
  }

  if (!is.null(Sigma_pre)) {
    if (!is.matrix(Sigma_pre) || nrow(Sigma_pre) != K || ncol(Sigma_pre) != K) stop("`Sigma_pre` must be K-by-K.", call. = FALSE)
    if (inherits(try(chol((Sigma_pre + t(Sigma_pre)) / 2), silent = TRUE), "try-error")) stop("`Sigma_pre` must be positive definite.", call. = FALSE)
  }
  if (!is.null(Sigma_post)) {
    if (!is.matrix(Sigma_post) || nrow(Sigma_post) != K || ncol(Sigma_post) != K) stop("`Sigma_post` must be K-by-K.", call. = FALSE)
    if (inherits(try(chol((Sigma_post + t(Sigma_post)) / 2), silent = TRUE), "try-error")) stop("`Sigma_post` must be positive definite.", call. = FALSE)
  }

  if (length(update_window) != 1L || update_window < 1 || update_window != as.integer(update_window)) stop("`update_window` must be positive integer.", call. = FALSE)
  if (length(grid_size) != 1L || grid_size < 2 || grid_size != as.integer(grid_size)) stop("`grid_size` must be integer >= 2.", call. = FALSE)
  if (length(ridge) != 1L || ridge <= 0) stop("`ridge` must be positive scalar.", call. = FALSE)

  new("MultivariateGaussianModel",
      name = name,
      mu_pre = mu_pre,
      mu_post = mu_post,
      Sigma_pre = Sigma_pre,
      Sigma_post = Sigma_post,
      method = method,
      update_window = as.integer(update_window),
      grid_size = as.integer(grid_size),
      prior_weights = as.numeric(prior_weights),
      ridge = ridge,
      max_grid_points = as.integer(max_grid_points))
}

BernoulliModel <- function(p_pre, p_post, name = "bernoulli") {
  stopifnot(length(p_pre) == 1L, length(p_post) == 1L, p_pre > 0, p_pre < 1, p_post > 0, p_post < 1)
  new("BernoulliModel", name = name, p_pre = p_pre, p_post = p_post)
}

AR1Model <- function(phi_pre, sigma_pre, mu_0 = 0, phi_post, sigma_post, mu_1 = 0, x0 = 0, name = "ar1") {
  stopifnot(length(phi_pre) == 1L, length(phi_post) == 1L, abs(phi_pre) < 1, abs(phi_post) < 1,
            sigma_pre > 0, sigma_post > 0, length(mu_0) == 1L, length(mu_1) == 1L, length(x0) == 1L)
  new("AR1Model", name = name, phi_pre = phi_pre, sigma_pre = sigma_pre, mu_pre = mu_0,
      phi_post = phi_post, sigma_post = sigma_post, mu_post = mu_1, x0 = x0)
}

# ---- path builders for composite models ----
.gaussian_mean_seq_post <- function(model, x) {
  n <- length(x)
  if (.is_scalar_spec(model@mean_post)) return(rep(as.numeric(model@mean_post), n))

  # composite post interval
  interval <- as.numeric(model@mean_post)
  if (model@method == "predictable") {
    out <- numeric(n)
    init <- mean(interval)
    u <- as.integer(model@update_window)
    for (i in seq_len(n)) {
      n_used <- ((i - 1L) %/% u) * u
      if (n_used <= 0L) out[i] <- init else out[i] <- .clamp(mean(x[seq_len(n_used)]), interval[1], interval[2])
    }
    return(out)
  }

  # mixture handled separately
  rep(NA_real_, n)
}

.gaussian_pre_params_t <- function(model, x, t) {
  if (.is_scalar_spec(model@mean_pre)) {
    mu <- as.numeric(model@mean_pre)
  } else {
    iv <- as.numeric(model@mean_pre)
    mu <- .clamp(mean(x[seq_len(t)]), iv[1], iv[2])
  }

  if (length(model@sd_pre) == 1L) {
    sd <- model@sd_pre
  } else {
    sd <- .gaussian_sd_mle(x[seq_len(t)], mu = mu, min_sd = model@min_sd)
  }
  c(mu = mu, sd = sd)
}

.gaussian_log_lr_cum_path <- function(model, x) {
  n <- length(x)
  logZ <- numeric(n)

  # pre denominator reset needed when pre mean composite OR pre sd unknown
  pre_reset <- .is_interval_spec(model@mean_pre) || length(model@sd_pre) == 0L

  sd1_fallback <- if (length(model@sd_post) == 1L) model@sd_post else if (length(model@sd_pre) == 1L) model@sd_pre else 1
  sd1_seq <- if (length(model@sd_post) == 1L) rep(model@sd_post, n) else .lagged_sd_seq(x, as.integer(model@update_window), sd1_fallback, model@min_sd)

  if (!.is_interval_spec(model@mean_post) || model@method == "predictable") {
    mu1_seq <- .gaussian_mean_seq_post(model, x)
    for (t in seq_len(n)) {
      pre_par <- .gaussian_pre_params_t(model, x, t)
      log_den <- if (pre_reset) {
        sum(stats::dnorm(x[seq_len(t)], mean = pre_par["mu"], sd = pre_par["sd"], log = TRUE))
      } else {
        sum(stats::dnorm(x[seq_len(t)], mean = as.numeric(model@mean_pre), sd = model@sd_pre, log = TRUE))
      }
      log_num <- sum(stats::dnorm(x[seq_len(t)], mean = mu1_seq[seq_len(t)], sd = sd1_seq[seq_len(t)], log = TRUE))
      logZ[t] <- log_num - log_den
    }
    return(logZ)
  }

  # mixture post over interval
  grid <- .gaussian_grid(as.numeric(model@mean_post), model@grid_size)
  w <- if (length(model@prior_weights) == length(grid)) model@prior_weights else rep(1 / length(grid), length(grid))
  w <- w / sum(w)
  logw <- log(pmax(w, .Machine$double.eps))

  for (t in seq_len(n)) {
    pre_par <- .gaussian_pre_params_t(model, x, t)
    log_den <- if (pre_reset) {
      sum(stats::dnorm(x[seq_len(t)], mean = pre_par["mu"], sd = pre_par["sd"], log = TRUE))
    } else {
      sum(stats::dnorm(x[seq_len(t)], mean = as.numeric(model@mean_pre), sd = model@sd_pre, log = TRUE))
    }

    loglikes <- vapply(grid, function(mu) {
      sum(stats::dnorm(x[seq_len(t)], mean = mu, sd = sd1_seq[seq_len(t)], log = TRUE))
    }, numeric(1L))
    log_num <- .logsumexp(logw + loglikes)
    logZ[t] <- log_num - log_den
  }

  logZ
}

.mv_mean_seq_post <- function(model, x) {
  n <- nrow(x)
  if (is.numeric(model@mu_post) && is.null(dim(model@mu_post))) {
    return(matrix(rep(as.numeric(model@mu_post), each = n), nrow = n))
  }

  # composite post box + predictable
  box <- model@mu_post
  out <- matrix(0, nrow = n, ncol = nrow(box))
  init <- rowMeans(box)
  u <- as.integer(model@update_window)
  for (i in seq_len(n)) {
    n_used <- ((i - 1L) %/% u) * u
    if (n_used <= 0L) out[i, ] <- init else out[i, ] <- .project_box(colMeans(x[seq_len(n_used), , drop = FALSE]), box)
  }
  out
}

.mv_pre_params_t <- function(model, x, t) {
  if (is.numeric(model@mu_pre) && is.null(dim(model@mu_pre))) {
    mu <- as.numeric(model@mu_pre)
  } else {
    mu <- .project_box(colMeans(x[seq_len(t), , drop = FALSE]), model@mu_pre)
  }

  if (!is.null(model@Sigma_pre)) {
    Sigma <- model@Sigma_pre
  } else {
    Sigma <- .mv_cov_mle(x[seq_len(t), , drop = FALSE], mu = mu, ridge = model@ridge)
  }
  list(mu = mu, Sigma = Sigma)
}

.mv_log_lr_cum_path <- function(model, x) {
  n <- nrow(x)
  logZ <- numeric(n)
  K <- ncol(x)

  pre_reset <- .is_matrix_box_spec(model@mu_pre) || is.null(model@Sigma_pre)

  Sigma1_fallback <- if (!is.null(model@Sigma_post)) model@Sigma_post else if (!is.null(model@Sigma_pre)) model@Sigma_pre else diag(1, K)
  Sigma1_seq <- if (!is.null(model@Sigma_post)) rep(list(model@Sigma_post), n) else .lagged_cov_seq(x, as.integer(model@update_window), Sigma1_fallback, model@ridge)

  post_box <- .is_matrix_box_spec(model@mu_post)
  if (!post_box || model@method == "predictable") {
    mu1_seq <- .mv_mean_seq_post(model, x)
    for (t in seq_len(n)) {
      pre_par <- .mv_pre_params_t(model, x, t)
      log_den <- if (pre_reset) {
        sum(.mvnorm_log_density(x[seq_len(t), , drop = FALSE], mean = pre_par$mu, Sigma = pre_par$Sigma))
      } else {
        sum(.mvnorm_log_density(x[seq_len(t), , drop = FALSE], mean = as.numeric(model@mu_pre), Sigma = model@Sigma_pre))
      }

      log_num <- 0
      for (i in seq_len(t)) {
        log_num <- log_num + .mvnorm_log_density(x[i, ], mean = mu1_seq[i, ], Sigma = Sigma1_seq[[i]])
      }
      logZ[t] <- log_num - log_den
    }
    return(logZ)
  }

  # mixture post box
  grid <- .mv_box_grid(model@mu_post, model@grid_size, model@max_grid_points)
  G <- nrow(grid)
  w <- if (length(model@prior_weights) == G) model@prior_weights else rep(1 / G, G)
  w <- w / sum(w)
  logw <- log(pmax(w, .Machine$double.eps))

  for (t in seq_len(n)) {
    pre_par <- .mv_pre_params_t(model, x, t)
    log_den <- if (pre_reset) {
      sum(.mvnorm_log_density(x[seq_len(t), , drop = FALSE], mean = pre_par$mu, Sigma = pre_par$Sigma))
    } else {
      sum(.mvnorm_log_density(x[seq_len(t), , drop = FALSE], mean = as.numeric(model@mu_pre), Sigma = model@Sigma_pre))
    }

    loglikes <- numeric(G)
    for (g in seq_len(G)) {
      ll <- 0
      for (i in seq_len(t)) {
        ll <- ll + .mvnorm_log_density(x[i, ], mean = grid[g, ], Sigma = Sigma1_seq[[i]])
      }
      loglikes[g] <- ll
    }
    log_num <- .logsumexp(logw + loglikes)
    logZ[t] <- log_num - log_den
  }

  logZ
}

# ---- density methods ----

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

setMethod("model_density", "MultivariateGaussianModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  if (regime == "pre") {
    mu <- if (is.numeric(object@mu_pre) && is.null(dim(object@mu_pre))) as.numeric(object@mu_pre) else rowMeans(object@mu_pre)
    Sigma <- if (!is.null(object@Sigma_pre)) object@Sigma_pre else diag(1, length(mu))
    return(.mvnorm_density(x, mean = mu, Sigma = Sigma))
  }

  mu <- if (is.numeric(object@mu_post) && is.null(dim(object@mu_post))) as.numeric(object@mu_post) else rowMeans(object@mu_post)
  Sigma <- if (!is.null(object@Sigma_post)) object@Sigma_post else diag(1, length(mu))
  .mvnorm_density(x, mean = mu, Sigma = Sigma)
})

setMethod("model_density", "BernoulliModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  prob <- if (regime == "pre") object@p_pre else object@p_post
  stats::dbinom(x, size = 1L, prob = prob)
})

setMethod("model_density", "AR1Model", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  prev <- object@x0
  if (!is.null(history) && length(history) >= 1L) prev <- history[length(history)]
  if (regime == "pre") {
    stats::dnorm(x, mean = object@mu_pre + object@phi_pre * prev, sd = object@sigma_pre)
  } else {
    stats::dnorm(x, mean = object@mu_post + object@phi_post * prev, sd = object@sigma_post)
  }
})

# ---- likelihood increment methods ----

setMethod("likelihood_increment", "GaussianModel", function(object, x, history = NULL, log = FALSE) {
  # check if the model is simple and quickly compute if so
  simple_fast <- .is_scalar_spec(object@mean_pre) && .is_scalar_spec(object@mean_post) &&
    length(object@sd_pre) == 1L && length(object@sd_post) == 1L
  if (simple_fast) {
    out_log <- stats::dnorm(x, mean = as.numeric(object@mean_post), sd = object@sd_post, log = TRUE) -
      stats::dnorm(x, mean = as.numeric(object@mean_pre), sd = object@sd_pre, log = TRUE)
    return(if (log) out_log else pmax(exp(out_log), .Machine$double.eps))
  }
  
  # update the likelihood increments based on the history
  if ((is.null(history) || length(history) == 0L) && is.numeric(x) && length(x) > 1L && is.null(dim(x))) {
    out <- numeric(length(x))
    h <- numeric(0)
    for (i in seq_along(x)) {
      out[i] <- likelihood_increment(object, x = x[i], history = h, log = log)
      h <- c(h, x[i])
    }
    return(out)
  }

  history <- if (is.null(history)) numeric(0) else as.numeric(history)
  xx <- c(history, x)
  logZ <- .gaussian_log_lr_cum_path(object, xx)
  out_log <- if (length(logZ) == 1L) logZ[1L] else logZ[length(logZ)] - logZ[length(logZ) - 1L]
  if (log) out_log else pmax(exp(out_log), .Machine$double.eps)
})

setMethod("likelihood_increment", "MultivariateGaussianModel", function(object, x, history = NULL, log = FALSE) {
  # check if the model is simple and quickly compute increments if so 
  simple_fast <- (is.numeric(object@mu_pre) && is.null(dim(object@mu_pre))) &&
    (is.numeric(object@mu_post) && is.null(dim(object@mu_post))) &&
    !is.null(object@Sigma_pre) && !is.null(object@Sigma_post) 
  if (simple_fast) {
    out_log <- .mvnorm_log_density(x, mean = as.numeric(object@mu_post), Sigma = object@Sigma_post) -
      .mvnorm_log_density(x, mean = as.numeric(object@mu_pre), Sigma = object@Sigma_pre)
    return(if (log) out_log else pmax(exp(out_log), .Machine$double.eps))
  }

  # compute increments under composite model
  x_row <- matrix(as.numeric(x), nrow = 1L)
  hh <- if (is.null(history) || length(history) == 0L) matrix(numeric(0), nrow = 0, ncol = ncol(x_row)) else as.matrix(history)
  xx <- rbind(hh, x_row)
  logZ <- .mv_log_lr_cum_path(object, xx)
  out_log <- if (length(logZ) == 1L) logZ[1L] else logZ[length(logZ)] - logZ[length(logZ) - 1L]
  if (log) out_log else pmax(exp(out_log), .Machine$double.eps)
})

setMethod("likelihood_increment", "Model", function(object, x, history = NULL, log = FALSE) {
  pre <- model_density(object, x = x, regime = "pre", history = history)
  post <- model_density(object, x = x, regime = "post", history = history)
  if (log) return(log(pmax(post, .Machine$double.eps)) - log(pmax(pre, .Machine$double.eps)))
  pmax(post / pmax(pre, .Machine$double.eps), .Machine$double.eps)
})
