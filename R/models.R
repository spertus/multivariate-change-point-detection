# Class: Model
# purpose: base class for pre/post change models used in likelihood increments
# slots:
#   name = character label
setClass("Model", slots = c(name = "character"))

# Class: UnivariateModel
# purpose: base class for models with scalar observations at each time t
setClass("UnivariateModel", contains = "Model")

# Class: MultivariateModel
# purpose: base class for models with vector observations at each time t
setClass("MultivariateModel", contains = "Model")

# Class: GaussianModel
# purpose: simple pre/post univariate Gaussian model
# slots:
#   mean_pre, sd_pre   = pre-change Gaussian parameters
#   mean_post, sd_post = post-change Gaussian parameters
setClass(
  "GaussianModel",
  contains = "UnivariateModel",
  slots = c(mean_pre = "numeric", sd_pre = "numeric", mean_post = "numeric", sd_post = "numeric")
)

# Class: BernoulliModel
# purpose: simple pre/post univariate Bernoulli model
# slots:
#   p_pre  = pre-change Bernoulli success probability
#   p_post = post-change Bernoulli success probability
setClass(
  "BernoulliModel",
  contains = "UnivariateModel",
  slots = c(p_pre = "numeric", p_post = "numeric")
)

# Class: AR1Model
# purpose: simple pre/post univariate Gaussian AR(1) model with intercept parameterization
# model form:
#   X_t | X_{t-1} ~ N(mu + phi * X_{t-1}, sigma)
# slots:
#   phi_pre, sigma_pre, mu_pre   = pre-change AR(1) parameters
#   phi_post, sigma_post, mu_post = post-change AR(1) parameters
#   x0 = numeric scalar initial lag value used when history is empty
setClass(
  "AR1Model",
  contains = "UnivariateModel",
  slots = c(
    phi_pre = "numeric",
    sigma_pre = "numeric",
    mu_pre = "numeric",
    phi_post = "numeric",
    sigma_post = "numeric",
    mu_post = "numeric",
    x0 = "numeric"
  )
)

# Class: MultivariateGaussianModel
# purpose: simple pre/post multivariate Gaussian model for K-vectors
# slots:
#   mu_pre, Sigma_pre   = pre-change mean vector and covariance matrix
#   mu_post, Sigma_post = post-change mean vector and covariance matrix
setClass(
  "MultivariateGaussianModel",
  contains = "MultivariateModel",
  slots = c(mu_pre = "numeric", Sigma_pre = "matrix", mu_post = "numeric", Sigma_post = "matrix")
)

# Constructor: GaussianModel
# inputs:
#   mean_pre, sd_pre   = pre-change mean and sd (sd_pre > 0)
#   mean_post, sd_post = post-change mean and sd (sd_post > 0)
#   name               = character label
# outputs:
#   GaussianModel object
GaussianModel <- function(mean_pre, sd_pre, mean_post, sd_post, name = "gaussian") {
  stopifnot(length(sd_pre) == 1L, length(sd_post) == 1L, sd_pre > 0, sd_post > 0)
  new(
    "GaussianModel",
    name = name,
    mean_pre = mean_pre,
    sd_pre = sd_pre,
    mean_post = mean_post,
    sd_post = sd_post
  )
}

# Constructor: BernoulliModel
# inputs:
#   p_pre  = pre-change success probability in (0,1)
#   p_post = post-change success probability in (0,1)
#   name   = character label
# outputs:
#   BernoulliModel object
BernoulliModel <- function(p_pre, p_post, name = "bernoulli") {
  stopifnot(length(p_pre) == 1L, length(p_post) == 1L, p_pre > 0, p_pre < 1, p_post > 0, p_post < 1)
  new("BernoulliModel",
    name = name,
    p_pre = p_pre,
    p_post = p_post
  )
}

# Constructor: AR1Model
# inputs:
#   phi_pre, sigma_pre, mu_0   = pre-change AR(1) parameters
#   phi_post, sigma_post, mu_1 = post-change AR(1) parameters
#   x0                         = initial lag value when no history is provided
#   name                       = character label
# outputs:
#   AR1Model object
AR1Model <- function(phi_pre, sigma_pre, mu_0 = 0, phi_post, sigma_post, mu_1 = 0, x0 = 0, name = "ar1") {
  stopifnot(length(phi_pre) == 1L, length(phi_post) == 1L, abs(phi_pre) < 1, abs(phi_post) < 1,
            sigma_pre > 0, sigma_post > 0, length(mu_0) == 1L, length(mu_1) == 1L, length(x0) == 1L)
  new(
    "AR1Model",
    name = name,
    phi_pre = phi_pre,
    sigma_pre = sigma_pre,
    mu_pre = mu_0,
    phi_post = phi_post,
    sigma_post = sigma_post,
    mu_post = mu_1,
    x0 = x0
  )
}

# Constructor: MultivariateGaussianModel
# inputs:
#   mu_pre, Sigma_pre   = pre-change mean vector and positive-definite covariance matrix
#   mu_post, Sigma_post = post-change mean vector and positive-definite covariance matrix
#   name                = character label
# outputs:
#   MultivariateGaussianModel object
MultivariateGaussianModel <- function(mu_pre, Sigma_pre, mu_post, Sigma_post, name = "mv-gaussian") {
  mu_pre <- as.numeric(mu_pre)
  mu_post <- as.numeric(mu_post)

  if (!is.matrix(Sigma_pre) || !is.matrix(Sigma_post)) {
    stop("`Sigma_pre` and `Sigma_post` must be matrices.", call. = FALSE)
  }
  if (nrow(Sigma_pre) != ncol(Sigma_pre) || nrow(Sigma_post) != ncol(Sigma_post)) {
    stop("Covariance matrices must be square.", call. = FALSE)
  }
  if (length(mu_pre) != nrow(Sigma_pre) || length(mu_post) != nrow(Sigma_post)) {
    stop("Mean vectors must match covariance dimensions.", call. = FALSE)
  }
  if (length(mu_pre) != length(mu_post)) {
    stop("Pre/post mean vectors must have same length.", call. = FALSE)
  }
  if (nrow(Sigma_pre) != nrow(Sigma_post)) {
    stop("Pre/post covariance matrices must have same dimensions.", call. = FALSE)
  }

  if (inherits(try(chol(Sigma_pre), silent = TRUE), "try-error") ||
      inherits(try(chol(Sigma_post), silent = TRUE), "try-error")) {
    stop("Covariance matrices must be symmetric positive definite.", call. = FALSE)
  }

  new(
    "MultivariateGaussianModel",
    name = name,
    mu_pre = mu_pre,
    Sigma_pre = Sigma_pre,
    mu_post = mu_post,
    Sigma_post = Sigma_post
  )
}

# Internal helper: .ar1_cond_mean
# purpose: evaluate conditional AR(1) mean under intercept parameterization
# inputs:
#   prev = numeric lagged observation X_{t-1}
#   phi  = numeric AR coefficient
#   mu   = numeric intercept
# outputs:
#   numeric conditional mean mu + phi * prev
.ar1_cond_mean <- function(prev, phi, mu) {
  mu + phi * prev
}

# Internal helper: .mvnorm_density
# purpose: evaluate multivariate Gaussian density for one vector or matrix of row-vectors
# inputs:
#   x     = numeric vector length p, or numeric matrix n-by-p
#   mean  = numeric vector length p
#   Sigma = p-by-p positive-definite covariance matrix
# outputs:
#   numeric scalar (vector input) or numeric length-n vector (matrix input)
.mvnorm_density <- function(x, mean, Sigma) {
  x_mat <- if (is.null(dim(x))) matrix(as.numeric(x), nrow = 1L) else as.matrix(x)
  p <- length(mean)

  if (ncol(x_mat) != p) {
    stop("Observation dimension does not match mean/covariance dimension.", call. = FALSE)
  }

  R <- chol(Sigma)
  log_det <- 2 * sum(log(diag(R)))
  log_const <- -0.5 * (p * log(2 * pi) + log_det)

  out <- numeric(nrow(x_mat))
  for (i in seq_len(nrow(x_mat))) {
    d <- x_mat[i, ] - mean
    z <- forwardsolve(t(R), d)
    quad <- sum(z^2)
    out[i] <- exp(log_const - 0.5 * quad)
  }

  if (is.null(dim(x))) out[1L] else out
}

# Method: model_density for GaussianModel
# inputs:
#   object  = GaussianModel object
#   x       = numeric scalar/vector
#   regime  = "pre" or "post"
#   history = unused
# outputs:
#   numeric density value(s) under requested regime
setMethod("model_density", "GaussianModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  if (regime == "pre") {
    stats::dnorm(x, mean = object@mean_pre, sd = object@sd_pre)
  } else {
    stats::dnorm(x, mean = object@mean_post, sd = object@sd_post)
  }
})

# Method: model_density for BernoulliModel
# inputs:
#   object  = BernoulliModel object
#   x       = numeric scalar/vector taking values in {0,1}
#   regime  = "pre" or "post"
#   history = unused
# outputs:
#   numeric mass value(s) under requested regime
setMethod("model_density", "BernoulliModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  prob <- if (regime == "pre") object@p_pre else object@p_post
  stats::dbinom(x, size = 1L, prob = prob)
})

# Method: model_density for AR1Model
# inputs:
#   object  = AR1Model object
#   x       = numeric scalar/vector at current time
#   regime  = "pre" or "post"
#   history = optional numeric vector; last value used as X_{t-1}
# outputs:
#   numeric conditional density value(s) under requested regime
setMethod("model_density", "AR1Model", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  prev <- object@x0
  if (!is.null(history) && length(history) >= 1L) {
    prev <- history[length(history)]
  }

  if (regime == "pre") {
    mean_t <- .ar1_cond_mean(prev = prev, phi = object@phi_pre, mu = object@mu_pre)
    stats::dnorm(x, mean = mean_t, sd = object@sigma_pre)
  } else {
    mean_t <- .ar1_cond_mean(prev = prev, phi = object@phi_post, mu = object@mu_post)
    stats::dnorm(x, mean = mean_t, sd = object@sigma_post)
  }
})

# Method: model_density for MultivariateGaussianModel
# inputs:
#   object  = MultivariateGaussianModel object
#   x       = numeric vector length K, or matrix n-by-K
#   regime  = "pre" or "post"
#   history = unused
# outputs:
#   numeric joint density value(s) under requested regime
setMethod("model_density", "MultivariateGaussianModel", function(object, x, regime = c("pre", "post"), history = NULL) {
  regime <- match.arg(regime)
  if (regime == "pre") {
    .mvnorm_density(x, mean = object@mu_pre, Sigma = object@Sigma_pre)
  } else {
    .mvnorm_density(x, mean = object@mu_post, Sigma = object@Sigma_post)
  }
})

# Method: likelihood_increment for Model
# inputs:
#   object  = Model subclass object
#   x       = current observation(s): scalar/vector for univariate; vector/row-matrix for multivariate
#   history = optional history object used by conditional models
# outputs:
#   numeric positive likelihood-ratio increments f_post(x)/f_pre(x)
setMethod("likelihood_increment", "Model", function(object, x, history = NULL) {
  pre <- model_density(object, x = x, regime = "pre", history = history)
  post <- model_density(object, x = x, regime = "post", history = history)
  pmax(post / pmax(pre, .Machine$double.eps), .Machine$double.eps)
})
