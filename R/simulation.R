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

# Internal helper: .gvar_simulate
# purpose: simulate N observations from a GaussianVARModel with a change at time nu.
#          Handles VAR(0) (IID) and VAR(p) cases uniformly.
# inputs:
#   model = GaussianVARModel
#   N     = integer number of observations to return
#   nu    = numeric change-point; Inf means no change (all pre-change)
# outputs:
#   numeric vector (K=1) or numeric N-by-K matrix (K>1)
.gvar_simulate <- function(model, N, nu) {
  K  <- length(model@mean_pre)
  p  <- length(model@Phi_pre)

  if (p == 0L) {
    pre  <- mvtnorm::rmvnorm(N, mean = model@mean_pre,  sigma = model@Sigma_pre)
    post <- mvtnorm::rmvnorm(N, mean = model@mean_post, sigma = model@Sigma_post)
    x <- if (is.finite(nu) && nu < N) {
      nu_i <- as.integer(nu)
      rbind(pre[seq_len(nu_i),     , drop = FALSE],
            post[seq(nu_i + 1L, N), , drop = FALSE])
    } else {
      pre
    }
  } else {
    I_K <- if (K == 1L) matrix(1) else diag(K)
    intercept_pre  <- as.numeric((I_K - Reduce("+", model@Phi_pre))  %*% model@mean_pre)
    intercept_post <- as.numeric((I_K - Reduce("+", model@Phi_post)) %*% model@mean_post)

    N_pad  <- N + p
    x      <- matrix(0, nrow = N_pad, ncol = K)
    for (i in seq_len(p)) x[i, ] <- model@mean_pre

    # Pre-draw all innovations: one Cholesky per regime, not per time step.
    nu_i   <- if (is.finite(nu)) min(as.integer(nu), N) else N
    n_pre  <- nu_i
    n_post <- N - nu_i
    innov  <- matrix(0, nrow = N_pad, ncol = K)
    innov[seq(p + 1L, p + n_pre), ] <-
      matrix(rnorm(n_pre * K), n_pre, K) %*% chol(model@Sigma_pre)
    if (n_post > 0L)
      innov[seq(p + n_pre + 1L, N_pad), ] <-
        matrix(rnorm(n_post * K), n_post, K) %*% chol(model@Sigma_post)

    for (t in seq(p + 1L, N_pad)) {
      obs_idx   <- t - p
      is_post   <- obs_idx > nu_i
      intercept <- if (is_post) intercept_post else intercept_pre
      Phi_use   <- if (is_post) model@Phi_post  else model@Phi_pre
      cond_mu   <- intercept
      for (j in seq_len(p))
        cond_mu <- cond_mu + as.numeric(Phi_use[[j]] %*% x[t - j, ])
      x[t, ] <- cond_mu + innov[t, ]
    }
    x <- x[seq(p + 1L, N_pad), , drop = FALSE]
  }

  if (K == 1L) as.vector(x[, 1L]) else x
}

# Constructor: DGP
# inputs:
#   generator   = GaussianVARModel (preferred) or function(N, K, nu, pre_params, post_params)
#   pre_params  = named list for generator settings (only used with function-based generators)
#   post_params = named list for generator settings (only used with function-based generators)
#   nu          = numeric scalar change-point (default Inf)
#   name        = character label
# outputs:
#   DGP object
DGP <- function(generator, pre_params = list(), post_params = list(), nu = Inf, name = "dgp") {
  if (is(generator, "GaussianVARModel")) {
    model  <- generator
    gen_fn <- function(N, K, nu, pre_params, post_params) .gvar_simulate(model, N, nu)
    return(new("DGP", generator = gen_fn, pre_params = list(), post_params = list(),
               nu = as.numeric(nu), name = as.character(name)))
  }
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
