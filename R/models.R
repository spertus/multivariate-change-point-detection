# Class: Model
# purpose: base class for pre/post change models used in likelihood increments
# slots:
#   name = character label
setClass("Model", slots = c(name = "character"))

# Class: GaussianModel
# purpose: simple pre/post Gaussian model
# slots:
#   mean_pre, sd_pre   = pre-change Gaussian parameters
#   mean_post, sd_post = post-change Gaussian parameters
setClass(
  "GaussianModel",
  contains = "Model",
  slots = c(mean_pre = "numeric", sd_pre = "numeric", mean_post = "numeric", sd_post = "numeric")
)

# Class: BernoulliModel
# purpose: simple pre/post Bernoulli model
# slots:
#   p_pre  = pre-change Bernoulli success probability
#   p_post = post-change Bernoulli success probability
setClass(
  "BernoulliModel",
  contains = "Model",
  slots = c(p_pre = "numeric", p_post = "numeric")
)

# Class: AR1Model
# purpose: simple pre/post Gaussian AR(1) model with intercept parameterization
# model form:
#   X_t | X_{t-1} ~ N(mu + phi * X_{t-1}, sigma)
# slots:
#   phi_pre, sigma_pre, mu_pre   = pre-change AR(1) parameters
#   phi_post, sigma_post, mu_post = post-change AR(1) parameters
#   x0 = numeric scalar initial lag value used when history is empty
setClass(
  "AR1Model",
  contains = "Model",
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

# Method: likelihood_increment for Model
# inputs:
#   object  = Model subclass object
#   x       = numeric scalar/vector, current observation(s)
#   history = optional numeric vector of prior observations
# outputs:
#   numeric positive likelihood-ratio increments f_post(x)/f_pre(x)
setMethod("likelihood_increment", "Model", function(object, x, history = NULL) {
  pre <- model_density(object, x = x, regime = "pre", history = history)
  post <- model_density(object, x = x, regime = "post", history = history)
  pmax(post / pmax(pre, .Machine$double.eps), .Machine$double.eps)
})
