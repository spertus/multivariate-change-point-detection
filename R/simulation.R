# Class: DGP
# purpose: wrapper for synthetic data generation under a change-point model
# slots:
#   generator   = function implementing data generation
#   pre_params  = list of pre-change generator parameters
#   post_params = list of post-change generator parameters
#   nu          = numeric scalar change-point (Inf means no change)
#   name        = character label
setClass(
  "DGP",
  slots = c(generator = "function", pre_params = "list", post_params = "list", nu = "numeric", name = "character")
)

# Constructor: DGP
# inputs:
#   generator   = function(N, K, nu, pre_params, post_params)
#   pre_params  = named list for pre-change generator settings
#   post_params = named list for post-change generator settings
#   nu          = numeric scalar change-point (default Inf)
#   name        = character label
# outputs:
#   DGP object
DGP <- function(generator, pre_params = list(), post_params = list(), nu = Inf, name = "dgp") {
  stopifnot(is.function(generator), length(nu) == 1L)
  new("DGP", generator = generator, pre_params = pre_params, post_params = post_params, nu = nu, name = name)
}

# Method: generate_stream for DGP
# inputs:
#   object = DGP object
#   N      = integer horizon
#   K      = integer number of streams
# outputs:
#   numeric vector (K=1) or numeric N-by-K matrix (K>1)
setMethod("generate_stream", "DGP", function(object, N, K = 1L) {
  object@generator(
    N = as.integer(N),
    K = as.integer(K),
    nu = object@nu,
    pre_params = object@pre_params,
    post_params = object@post_params
  )
})

# Function: expand_dgp_grid
# purpose: create a list of DGP objects by varying selected pre-change parameters over a grid
# inputs:
#   template_dgp = DGP object used as template
#   param_grid   = named list of vectors to expand with expand.grid
# outputs:
#   list of DGP objects, one per row of expanded parameter grid
expand_dgp_grid <- function(template_dgp, param_grid) {
  stopifnot(is(template_dgp, "DGP"), is.list(param_grid), length(param_grid) > 0)
  keys <- names(param_grid)
  if (is.null(keys) || any(keys == "")) {
    stop("`param_grid` entries must be named.", call. = FALSE)
  }

  grid <- expand.grid(param_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    pre <- template_dgp@pre_params
    for (nm in keys) {
      pre[[nm]] <- grid[i, nm][[1L]]
    }
    out[[i]] <- DGP(
      generator = template_dgp@generator,
      pre_params = pre,
      post_params = template_dgp@post_params,
      nu = template_dgp@nu,
      name = sprintf("%s_%d", template_dgp@name, i)
    )
  }

  out
}

# Function: default_gaussian_dgp
# purpose: generate IID Gaussian streams with a single shared change-point
# inputs:
#   N           = integer horizon
#   K           = integer number of streams
#   nu          = numeric change-point; if finite and < N, post-change starts at nu+1
#   pre_params  = list(mean, sd) for pre-change Gaussian draws
#   post_params = list(mean, sd) for post-change Gaussian draws
# outputs:
#   numeric vector (K=1) or numeric N-by-K matrix (K>1)
default_gaussian_dgp <- function(N, K = 1L, nu = Inf, pre_params = list(mean = 0, sd = 1), post_params = list(mean = 1, sd = 1)) {
  out <- matrix(0, nrow = N, ncol = K)
  for (j in seq_len(K)) {
    pre <- stats::rnorm(N, mean = pre_params$mean, sd = pre_params$sd)
    post <- stats::rnorm(N, mean = post_params$mean, sd = post_params$sd)
    if (is.finite(nu) && nu < N) {
      out[, j] <- c(pre[1:nu], post[(nu + 1):N])
    } else {
      out[, j] <- pre
    }
  }
  if (K == 1L) {
    return(as.vector(out[, 1L]))
  }
  out
}

# Function: default_multivariate_gaussian_dgp
# purpose: generate IID multivariate Gaussian vectors with a shared change-point
# inputs:
#   N           = integer horizon
#   nu          = numeric change-point; if finite and < N, post-change starts at nu+1
#   pre_params  = list(mu, Sigma) for pre-change MVN draws
#   post_params = list(mu, Sigma) for post-change MVN draws
# outputs:
#   numeric N-by-K matrix
default_multivariate_gaussian_dgp <- function(N,
                                              K = 2L,
                                              nu = Inf,
                                              pre_params = list(mu = rep(0, K), Sigma = diag(K)),
                                              post_params = list(mu = rep(1, K), Sigma = diag(K))) {
  mu_pre_raw <- if (!is.null(pre_params$mu)) pre_params$mu else pre_params$mean
  mu_post_raw <- if (!is.null(post_params$mu)) post_params$mu else post_params$mean
  if (is.null(mu_pre_raw) || is.null(mu_post_raw)) {
    stop("`pre_params` and `post_params` must include `mu` (or alias `mean`).", call. = FALSE)
  }
  mu_pre <- as.numeric(mu_pre_raw)
  mu_post <- as.numeric(mu_post_raw)
  Sigma_pre <- pre_params$Sigma
  Sigma_post <- post_params$Sigma

  if (!is.matrix(Sigma_pre) || !is.matrix(Sigma_post)) {
    stop("`pre_params$Sigma` and `post_params$Sigma` must be matrices.", call. = FALSE)
  }
  if (nrow(Sigma_pre) != ncol(Sigma_pre) || nrow(Sigma_post) != ncol(Sigma_post)) {
    stop("Covariance matrices must be square.", call. = FALSE)
  }
  if (length(mu_pre) != nrow(Sigma_pre) || length(mu_post) != nrow(Sigma_post)) {
    stop("Mean vectors must match covariance dimensions.", call. = FALSE)
  }
  if (length(mu_pre) != K || length(mu_post) != K) {
    stop("`K` must match the dimension of provided means/covariances.", call. = FALSE)
  }

  if (inherits(try(chol(Sigma_pre), silent = TRUE), "try-error") ||
      inherits(try(chol(Sigma_post), silent = TRUE), "try-error")) {
    stop("Covariance matrices must be symmetric positive definite.", call. = FALSE)
  }

  pre <- mvtnorm::rmvnorm(N, mean = mu_pre, sigma = Sigma_pre)
  post <- mvtnorm::rmvnorm(N, mean = mu_post, sigma = Sigma_post)

  if (is.finite(nu) && nu < N) {
    rbind(pre[1:nu, , drop = FALSE], post[(nu + 1):N, , drop = FALSE])
  } else {
    pre
  }
}

# Internal helper: .run_single
# purpose: run a detector/TSM pipeline on one generated dataset (univariate or multivariate)
# inputs:
#   detector = Detector subclass object
#   tsm      = TSM object for univariate data, or list of length K for multivariate data
#   x        = numeric vector (N) or matrix (N-by-K)
#   combiner = optional Combiner object; defaults to flat average in multistream case
#   weights  = optional numeric weights for weighted combiners
#   log      = logical; if TRUE run increment generation/combining/detection on log scale
# outputs:
#   list returned by run_detector()
.run_single <- function(detector, tsm, x, combiner = NULL, weights = NULL, log = FALSE) {
  if (is.null(dim(x))) {
    increments <- compute_increments(tsm, x, log = log)
    return(run_detector(detector, increments, log = log))
  }

  if (!is.list(tsm) && is(tsm@model, "MultivariateModel")) {
    increments <- compute_increments(tsm, x, log = log)
    return(run_detector(detector, increments, log = log))
  }

  k <- ncol(x)
  if (!is.list(tsm)) {
    tsm <- replicate(k, tsm, simplify = FALSE)
  }
  if (length(tsm) != k) {
    stop("For multivariate data, `tsm` must be a list with one element per stream.", call. = FALSE)
  }

  marg <- sapply(seq_len(k), function(j) compute_increments(tsm[[j]], x[, j], log = log))
  if (is.null(combiner)) {
    combiner <- AverageCombiner()
  }
  combined_increments <- combine_streams(combiner, marg, weights = weights, log = log)
  run_detector(detector, combined_increments, log = log)
}

# Function: run_simulation
# purpose: estimate detector operating characteristics by Monte Carlo simulation
# inputs:
#   detector = Detector object or list of Detector objects
#   tsm      = TSM object (or list used by .run_single for multistream)
#   dgp      = DGP object or list of DGP objects
#   n_rep    = integer number of Monte Carlo replications per detector-DGP pair
#   N        = integer horizon (length of streams) per rep
#   K        = integer number of streams generated by DGP
#   combiner = optional Combiner object for multistream combination
#   weights  = optional numeric weights passed to combiner
#   seed     = integer random seed
#   log      = logical; if TRUE run detector pipeline using log-increments
# outputs:
#   data.frame with one row per detector-DGP pair and columns:
#     detector, dgp, n_rep, horizon, nu, false_alarm_prob, ARL, ADD
run_simulation <- function(detector, tsm, dgp, n_rep = 200, N = 500, K = 1L, combiner = NULL, weights = NULL, seed = 1L, log = FALSE) {
  if (!is.list(detector)) detector <- list(detector)
  if (!is.list(dgp)) dgp <- list(dgp)
  set.seed(seed)

  rows <- list()
  idx <- 1L

  for (di in seq_along(detector)) {
    for (gi in seq_along(dgp)) {
      stop_times <- rep(Inf, n_rep)
      nu <- dgp[[gi]]@nu

      for (r in seq_len(n_rep)) {
        x <- generate_stream(dgp[[gi]], N = N, K = K)
        out <- .run_single(detector[[di]], tsm, x, combiner = combiner, weights = weights, log = log)
        stop_times[r] <- out$stopping_time
      }

      fa <- mean(is.finite(stop_times) & (!is.finite(nu) | stop_times <= nu)) # false alarm
      arl <- mean(pmin(stop_times, N+1)) # average run length
      cad <- if (is.finite(nu)) mean(pmax(stop_times - nu, 0), na.rm = TRUE) else NA_real_ # conditional average delay

      rows[[idx]] <- data.frame(
        detector = detector[[di]]@name,
        dgp = dgp[[gi]]@name,
        n_rep = n_rep,
        horizon = N,
        nu = nu,
        false_alarm_prob = fa,
        ARL = arl,
        ADD = cad,
        CAD = cad,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}
