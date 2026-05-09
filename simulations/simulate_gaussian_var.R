## simulate_gaussian_var.R
## Simulation study: change-point detection in Gaussian VAR(p) processes.
##
## Design follows simulation_design.tex:
##   - K streams, VAR(p), diagonal coefficient matrices (no cross-stream causality)
##   - Phi identical across streams: phi_{kj} = phi_j for all k
##   - Sigma unchanged at change-point; long-run mean shifts mu_0 -> mu_1 at nu
##   - Pre-change model P is simple (known mu_0, {phi_j})
##   - Post-change model Q is composite; three increment strategies:
##       (i)  oracle      : plug in true mu_1
##       (ii) misspecified: plug in mu_1 / 2
##       (iii) adaptive   : lagged running sample mean (predictable plug-in)
##   - Combiners: average, product (independent only), universal portfolio
##   - Performance metric: CAD (conditional average delay given alarm before horizon)
##
## Usage:
##   source("simulations/simulate_gaussian_var.R")
##   results <- run_var_simulations(n_rep = 500, seed = 1)

devtools::load_all()
library(dplyr)
library(tidyr)

# ---- simulation grid -------------------------------------------------------

sim_grid <- expand.grid(
  p           = c(0, 1, 2, 3),          # AR order (0 = iid)
  K           = c(1, 2, 10, 100),       # number of streams
  nu          = c(10, 50, 100, 1000),   # change-point location
  delta_norm  = seq(0.1, 2, length.out = 20),  # ||mu_1||_2
  sparse      = c(TRUE, FALSE),         # sparse (1 stream) vs dense (all K streams) change
  alpha       = c(0.001, 0.01, 0.1),   # control level
  criterion   = c("ARL", "PFA"),        # error criterion
  independent = c(TRUE, FALSE),         # diagonal vs. exchangeable Sigma
  combiner    = c("average", "product", "universal_portfolio"),
  stringsAsFactors = FALSE
) %>%
  # product combining only valid when streams are independent
  dplyr::filter(!(combiner == "product" & !independent)) |>
  # univariate: combining is trivial; keep only average for K=1
  dplyr::filter(
    !(K == 1 & combiner %in% c("product", "universal_portfolio"))
  ) |>
  # sample size: max(500, 5 * nu)
  mutate(N = pmax(500L, 5L * nu))

cat("Simulation grid size:", nrow(sim_grid), "scenarios\n")

# ---- fixed AR coefficients (identical across streams) ----------------------
# phi_j for orders 0..3; stationary by design
AR_COEFS <- list(
  `0` = numeric(0),
  `1` = 0.5,
  `2` = c(0.4, 0.2),
  `3` = c(0.5, 0.2, 0.1)
)

# ---- helpers ---------------------------------------------------------------

# Build a diagonal K-by-K coefficient matrix from a scalar phi_j
.phi_mat <- function(phi_j, K) diag(phi_j, K)

# Build the list of p K-by-K coefficient matrices for a VAR(p)
.make_Phi <- function(phi_vec, K) {
  if (length(phi_vec) == 0L) return(list())
  lapply(phi_vec, .phi_mat, K = K)
}

# Build Sigma: diagonal (independent) or exchangeable with rho = 0.5
.make_Sigma <- function(K, independent) {
  if (independent || K == 1L) return(diag(K))
  rho <- 0.5
  (1 - rho) * diag(K) + rho * matrix(1, K, K)
}

# Build the true post-change mean vector given norm delta and sparsity
.make_mu1 <- function(K, delta_norm, sparse) {
  if (K == 1L) return(delta_norm)
  if (sparse) {
    v <- rep(0, K); v[1] <- delta_norm; v
  } else {
    rep(delta_norm / sqrt(K), K)
  }
}

# Simulate one VAR(p) path of length N with change at nu
.simulate_var <- function(Phi_list, Sigma, c_pre, c_post, nu, N, K, x0) {
  p   <- length(Phi_list)
  chol_Sigma <- chol(Sigma)
  x   <- matrix(0, nrow = N + p, ncol = K)
  x[seq_len(p), ] <- matrix(rep(x0, p), nrow = p, byrow = TRUE)

  for (t in seq_len(N)) {
    idx   <- t + p
    c_t   <- if (t <= nu) c_pre else c_post
    ar    <- if (p > 0L) {
      Reduce("+", lapply(seq_len(p), function(l) as.numeric(Phi_list[[l]] %*% x[idx - l, ])))
    } else {
      rep(0, K)
    }
    eps   <- as.numeric(matrix(rnorm(K), nrow = 1) %*% chol_Sigma)
    x[idx, ] <- c_t + ar + eps
  }
  x[(p + 1L):(N + p), , drop = FALSE]
}

# Vectorised log-increments for one univariate AR(p) stream.
#
# For the Gaussian AR(p) model with unit innovation variance and diagonal Phi,
# the log-likelihood ratio at time t reduces to a closed form:
#
#   log LR_t = (mu1_t - mu0) * (x_t - ar_t) - 0.5*(mu1_t^2 - mu0^2)
#
# where ar_t = sum_j phi_j * x_{t-j}  (the AR contribution, same under both
# regimes since Phi is shared), and mu{0,1}_t are the pre/post long-run means
# (scalars), converted to intercepts via c = m*(1 - sum(phi)).
#
# This is O(N) per stream, replacing the previous O(N^2) growing-history loop.
.ar_log_lr_vectorised <- function(xk, phi_vec, mu0, mu1_vec) {
  N <- length(xk)
  p <- length(phi_vec)

  # AR residual: x_t - sum_j phi_j * x_{t-j}
  # pad the front with 0 (pre-sample values = 0 = mu0)
  if (p == 0L) {
    ar_contrib <- rep(0, N)
  } else {
    xk_padded <- c(rep(0, p), xk)
    # embed() gives an N x (p+1) matrix; columns 2..p+1 are lags 1..p
    lag_mat    <- embed(xk_padded, p + 1L)[, -1L, drop = FALSE]
    ar_contrib <- as.numeric(lag_mat %*% phi_vec)
  }

  # Innovation under each regime: e_t = x_t - c_k - ar_t
  # where c_k = m_k * (1 - sum(phi)) = intercept for long-run mean m_k
  sum_phi <- sum(phi_vec)
  c0      <- mu0     * (1 - sum_phi)
  c1_vec  <- mu1_vec * (1 - sum_phi)   # length-N vector (adaptive) or scalar

  innov0 <- xk - c0      - ar_contrib
  innov1 <- xk - c1_vec  - ar_contrib

  # log LR = log N(innov1; 0,1) - log N(innov0; 0,1)
  #        = -0.5*(innov1^2 - innov0^2)
  -0.5 * (innov1^2 - innov0^2)
}

# Compute combined log-increments for one strategy (vectorised, O(N) per stream)
# strategy: "oracle" | "misspecified" | "adaptive"
.compute_log_increments <- function(x, model_pre, mu1, strategy, K, Phi_list, Sigma) {
  N       <- nrow(x)
  phi_vec <- if (length(Phi_list) == 0L) numeric(0L) else
               sapply(Phi_list, function(A) A[1L, 1L])
  mu0     <- 0   # pre-change long-run mean (scalar, same for all streams)

  # per-stream post-change long-run mean (scalar for oracle/misspec; varies for adaptive)
  marg_inc <- matrix(0, nrow = N, ncol = K)

  for (k in seq_len(K)) {
    xk <- x[, k]

    mu1_vec <- switch(strategy,
      oracle       = mu1[k],
      misspecified = mu1[k] / 2,
      adaptive     = c(0, cumsum(xk[-N]) / seq_len(N - 1L))  # lagged running mean, O(N)
    )

    marg_inc[, k] <- .ar_log_lr_vectorised(xk, phi_vec, mu0, mu1_vec)
  }

  if (K == 1L) as.numeric(marg_inc) else marg_inc
}

# Run the detector on (possibly combined) log-increments
.run_sr <- function(log_inc_or_matrix, alpha, criterion, combiner_name, K) {
  det <- ShiryaevRobertsDetector(alpha = alpha, criterion = criterion)

  if (K == 1L || is.numeric(log_inc_or_matrix) && is.null(dim(log_inc_or_matrix))) {
    inc_vec <- as.numeric(log_inc_or_matrix)
  } else {
    comb <- switch(combiner_name,
      average             = AverageCombiner(),
      product             = ProductCombiner(),
      universal_portfolio = UniversalPortfolioCombiner(mode = "sparse")
    )
    inc_vec <- combine_streams(comb, log_inc_or_matrix, log = TRUE)
  }

  run_detector(det, inc_vec, log = TRUE)
}

# ---- single-scenario Monte Carlo -------------------------------------------

#' Run Monte Carlo for one row of sim_grid
#' @param sc  one-row data.frame from sim_grid
#' @param n_rep number of Monte Carlo replications
#' @return data.frame with summary statistics (CAD, alarm_prob)
run_scenario <- function(sc, n_rep = 500) {
  p     <- sc$p
  K     <- sc$K
  nu    <- sc$nu
  N     <- sc$N
  delta <- sc$delta_norm
  sparse     <- sc$sparse
  alpha      <- sc$alpha
  criterion  <- sc$criterion
  independent <- sc$independent
  combiner   <- sc$combiner

  phi_vec  <- AR_COEFS[[as.character(p)]]
  Phi_list <- .make_Phi(phi_vec, K)
  Sigma    <- .make_Sigma(K, independent)
  mu0      <- rep(0, K)
  mu1      <- .make_mu1(K, delta, sparse)
  c_pre    <- if (length(Phi_list) == 0L) mu0 else .var_intercept_from_mean(mu0, Phi_list)
  c_post   <- if (length(Phi_list) == 0L) mu1 else .var_intercept_from_mean(mu1, Phi_list)

  # Pre-change model object (for marginal density in adaptive strategy)
  if (K == 1L) {
    model_pre <- if (p == 0L) {
      GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 0, sd_post = 1)
    } else {
      ARpModel(phi_pre = phi_vec, sigma_pre = 1, mean_pre = 0,
               phi_post = phi_vec, sigma_post = 1, mean_post = 0)
    }
  } else {
    model_pre <- NULL  # not used directly for K > 1
  }

  results <- replicate(n_rep, {
    x <- .simulate_var(Phi_list, Sigma, c_pre, c_post, nu, N, K, x0 = mu0)

    run_one <- function(strategy) {
      raw <- .compute_log_increments(x, model_pre, mu1, strategy, K, Phi_list, Sigma)
      out <- .run_sr(raw, alpha, criterion, combiner, K)
      st  <- out$stopping_time
      cad   <- if (is.finite(st) && st > nu) st - nu else NA_real_
      alarm <- as.numeric(is.finite(st))
      c(cad, alarm)
    }
    r_o <- run_one("oracle")
    r_m <- run_one("misspecified")
    r_a <- run_one("adaptive")

    c(cad_oracle = r_o[1], cad_misspec = r_m[1], cad_adaptive = r_a[1],
      alarm_oracle = r_o[2], alarm_misspec = r_m[2], alarm_adaptive = r_a[2])
  }, simplify = "matrix")

  data.frame(
    p = p, K = K, nu = nu, N = N, delta_norm = delta,
    sparse = sparse, alpha = alpha, criterion = criterion,
    independent = independent, combiner = combiner,
    CAD_oracle      = mean(results["cad_oracle",    ], na.rm = TRUE),
    CAD_misspec     = mean(results["cad_misspec",   ], na.rm = TRUE),
    CAD_adaptive    = mean(results["cad_adaptive",  ], na.rm = TRUE),
    alarm_prob_oracle   = mean(results["alarm_oracle",  ]),
    alarm_prob_misspec  = mean(results["alarm_misspec", ]),
    alarm_prob_adaptive = mean(results["alarm_adaptive",])
  )
}

# ---- main entry point ------------------------------------------------------

#' Run all simulation scenarios
#' @param n_rep   Monte Carlo replications per scenario
#' @param seed    RNG seed
#' @param grid    data.frame of scenarios (defaults to full sim_grid)
#' @param outfile path to save results as .rds (NULL to skip)
#' @return data.frame of results (one row per scenario)
run_var_simulations <- function(n_rep   = 500,
                                seed    = 1,
                                grid    = sim_grid,
                                outfile = "simulations/var_simulation_results.rds") {
  set.seed(seed)
  n_scenarios <- nrow(grid)
  cat("Running", n_scenarios, "scenarios x", n_rep, "replications each.\n")

  results <- vector("list", n_scenarios)
  for (i in seq_len(n_scenarios)) {
    if (i %% 50 == 0L || i == 1L)
      cat(sprintf("  Scenario %d / %d\n", i, n_scenarios))
    results[[i]] <- tryCatch(
      run_scenario(grid[i, ], n_rep = n_rep),
      error = function(e) {
        warning(sprintf("Scenario %d failed: %s", i, conditionMessage(e)))
        NULL
      }
    )
  }

  out <- bind_rows(results)
  if (!is.null(outfile)) {
    saveRDS(out, outfile)
    cat("Results saved to", outfile, "\n")
  }
  out
}

# ---- runtime benchmark -----------------------------------------------------
#
# The most expensive scenario class: K=100, adaptive strategy, UP combiner,
# large nu => large N = max(500, 5*nu).
#
# The adaptive increment computation is O(N) per stream via .ar_log_lr_vectorised(),
# which uses embed() for the lag matrix and cumsum() for the running mean.
# This replaces the former O(N^2) growing-history loop.
#
# Usage:
#   bm <- benchmark_runtime()

benchmark_runtime <- function(n_rep_full = 500, seed = 42) {
  set.seed(seed)

  # --- benchmark adaptive increments (now O(N) via .ar_log_lr_vectorised)
  N_small  <- 200L
  phi_vec  <- AR_COEFS[["3"]]
  x_small  <- matrix(rnorm(N_small * 10L), N_small, 10L)

  t_adaptive <- system.time({
    for (k in seq_len(10L)) {
      xk <- x_small[, k]
      mu1_adap <- c(0, cumsum(xk[-N_small]) / seq_len(N_small - 1L))
      .ar_log_lr_vectorised(xk, phi_vec, mu0 = 0, mu1_vec = mu1_adap)
    }
  })["elapsed"] / 10L

  # O(N) per stream: extrapolate linearly to N=5000, then x K=100
  N_hard              <- 5000L
  K_hard              <- 100L
  sec_adaptive_hard   <- t_adaptive * (N_hard / N_small) * K_hard

  # --- benchmark UP combiner (fast: capped at max_points=300 portfolios)
  K_small  <- 5L;  N_up <- 200L
  marg_small <- matrix(rnorm(N_up * K_small), N_up, K_small)
  up <- UniversalPortfolioCombiner(resolution = 6)
  t_up <- system.time(for (i in 1:5) combine_streams(up, marg_small, log = TRUE))["elapsed"] / 5
  # Scales linearly in N; portfolio count is capped so K-scaling is negligible
  sec_up_hard <- t_up * (N_hard / N_up)

  sec_per_rep_hard <- sec_adaptive_hard + sec_up_hard

  n_hard <- nrow(sim_grid[
    sim_grid$K == K_hard & sim_grid$combiner == "universal_portfolio", ])

  # --- benchmark easy scenarios (oracle/misspec, small K)
  x_easy  <- rnorm(500)
  m_easy  <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  t_easy  <- system.time(for (i in 1:20) {
    inc <- compute_increments(TSM(m_easy), x_easy, log = TRUE)
    run_detector(ShiryaevRobertsDetector(alpha = 0.001), inc, log = TRUE)
  })["elapsed"] / 20

  n_easy <- nrow(sim_grid) - n_hard

  # --- extrapolate ---
  est_hard_hrs  <- (sec_per_rep_hard * n_rep_full * n_hard)  / 3600
  est_easy_hrs  <- (t_easy           * n_rep_full * n_easy)  / 3600
  est_total_hrs <- est_hard_hrs + est_easy_hrs

  cat(sprintf("\n=== Runtime estimate (single core, %d reps/scenario) ===\n",
              n_rep_full))
  cat(sprintf(
    "  Adaptive benchmark: N=%-4d, K=10  => %.4f sec/stream  (O(N) vectorised)\n",
    N_small, t_adaptive))
  cat(sprintf(
    "  Extrapolated to K=100, N=5000:   ~%.3f sec/rep (~%.4f hrs/rep)\n",
    sec_per_rep_hard, sec_per_rep_hard / 3600))
  cat(sprintf(
    "  Hard scenarios (K=100, UP):     %5d scenarios  => ~%.0f hrs total\n",
    n_hard, est_hard_hrs))
  cat(sprintf(
    "  Other scenarios:                %5d scenarios  => ~%.1f hrs total\n",
    n_easy, est_easy_hrs))
  cat(sprintf(
    "  TOTAL (single core):  ~%.0f hrs  (~%.1f days)\n\n",
    est_total_hrs, est_total_hrs / 24))
  cat(paste(
    "  NOTE: Adaptive increments are O(N) per stream (vectorised via embed()+cumsum()).",
    "  Runtime is dominated by the UP combiner (discrete portfolio grid).",
    sep = "\n  "
  ), "\n")

  invisible(list(
    sec_per_rep_hard  = sec_per_rep_hard,
    sec_per_rep_easy  = t_easy,
    n_hard            = n_hard,
    n_easy            = n_easy,
    n_rep_full        = n_rep_full,
    est_hard_hrs      = est_hard_hrs,
    est_easy_hrs      = est_easy_hrs,
    est_total_hrs     = est_total_hrs,
    est_total_days    = est_total_hrs / 24
  ))
}
