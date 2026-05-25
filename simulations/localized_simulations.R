# localized_simulations.R
# Simulation study for LocalizedDetector (Theorem 3, Spertus et al.)
#
# A single LocalizedDetector with uniform Bonferroni allowance (gamma_k = 1/K)
# is run on K-stream VAR(p) data drawn from the same DGP as joint_simulations.R.
# Only ARL control is evaluated (PFA spending schedules will be added later).
#
# Per-stream BoundedModel TSMs provide the increments; the LocalizedDetector
# then monitors all K streams simultaneously and declares an alarm when any
# stream's S-R statistic crosses the shared threshold 1/alpha.
#
# Output columns:
#   Condition: p, K, nu, delta, sparse, alpha, independent
#   Global:    ARL, CAD, PFA, power
#   Per-stream (k = 1..K_MAX): stream_k_ARL, stream_k_CAD
#              (NA for k > K in any given condition)
#
# ── Timing estimate (single core) ───────────────────────────────────────────
#   Quick mode (RUN_FULL = FALSE): ~6 sec  [N=1000, 100 reps, 10 conditions]
#   Full  mode (RUN_FULL = TRUE) : ~3 hrs  [N=3000, 500 reps, 960 conditions]
#
# Full grid: p(2) x K(3) x nu(2) x delta(20) x sparse(2) x alpha(1)
#            x independent(2) = 960 conditions, each with 1 detector x 500 reps.
# Use parallel::mclapply or future.apply on the outer (condition) loop to
# reduce wall time.
# ────────────────────────────────────────────────────────────────────────────

library(multichangepoints)

# ── 0. Mode switch ──────────────────────────────────────────────────────────

RUN_FULL <- FALSE          # set TRUE for the full grid

# ── 1. Fixed simulation constants ───────────────────────────────────────────

MU_PRE    <- 0.3
SIGMA_INN <- 0.06
RHO_CORR  <- 0.30

AR_COEFS <- list(
  "0" = list(),
  "1" = list(0.4),
  "2" = list(0.3, 0.1)
)

K_MAX <- 10L   # maximum K across all conditions; sets the width of stream columns

# ── 2. Parameter grids ──────────────────────────────────────────────────────
# Criterion is always ARL; PFA spending schedules will be added in a future script.
# delta_norm is log-spaced for better resolution at small signals.

quick_params <- list(
  p           = 0L,
  K           = c(1L, 2L),
  nu          = 100L,
  delta_norm  = c(0, exp(seq(log(0.01), log(0.2), length.out = 4L))),
  sparse      = FALSE,
  alpha       = 0.0001,
  independent = TRUE
)

full_params <- list(
  p           = c(0L, 1L),
  K           = c(1L, 2L, 10L),
  nu          = c(10L, 1000L),
  delta_norm  = c(0, exp(seq(log(0.01), log(0.2), length.out = 19L))),
  sparse      = c(FALSE, TRUE),
  alpha       = 0.0001,
  independent = c(TRUE, FALSE)
)

param_grid <- do.call(expand.grid,
  c(if (RUN_FULL) full_params else quick_params,
    list(stringsAsFactors = FALSE)))

N_OBS <- if (RUN_FULL) 3000L else 1000L
N_REP <- if (RUN_FULL) 500L  else 100L

# ── 3. Single-replicate runner ───────────────────────────────────────────────

# Returns global and per-stream first alarm times, false-alarm flags, and
# post-change delays.
#
# Implementation: K independent ShiryaevRobertsDetectors at alpha/K each
# (threshold = K/alpha per stream) with multiple_alarms = TRUE.  This is
# mathematically equivalent to LocalizedDetector with uniform 1/K Bonferroni
# allowance and threshold 1/alpha (the recursion scales by K).  Using
# multiple_alarms = TRUE lets each stream reset after a false alarm and keep
# accumulating evidence, so power reflects whether any stream detects the
# change within [nu+1, N] rather than being killed by pre-change alarms.
.run_rep_localized <- function(N, K, nu, mu0, mu1, phi_scalars, sigma_mat,
                               eta_vec, alpha) {
  Phi_lst   <- lapply(phi_scalars, function(ph) ph * diag(K))
  dgp_model <- GaussianVARModel(
    Phi_pre   = Phi_lst, Sigma_pre  = sigma_mat, mean_pre  = mu0,
    Phi_post  = Phi_lst, Sigma_post = sigma_mat, mean_post = mu1
  )
  x <- generate_stream(DGP(dgp_model, nu = nu), N = N)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)

  # BoundedModel TSM increments: one column per stream
  inc_mat <- matrix(NA_real_, nrow = N, ncol = K)
  for (k in seq_len(K))
    inc_mat[, k] <- compute_increments(TSM(BoundedModel(eta = eta_vec[k])),
                                       x[, k], log = TRUE)

  # K independent S-R detectors at Bonferroni-corrected level alpha/K,
  # each resetting after alarms so post-change detection is not blocked.
  alpha_k <- alpha / K
  stream_alarm_times <- vector("list", K)
  for (k in seq_len(K)) {
    det_k <- ShiryaevRobertsDetector(alpha = alpha_k, criterion = "ARL",
                                     multiple_alarms = TRUE)
    out_k <- run_detector(det_k, inc_mat[, k], log = TRUE)
    stream_alarm_times[[k]] <- out_k$alarm_times
  }

  # Global alarm times = union across streams; first alarm for ARL
  all_global <- sort(unique(unlist(stream_alarm_times)))
  first_global <- if (length(all_global) > 0L) all_global[1L] else Inf

  # First post-change alarm for power and CAD
  post_global  <- all_global[all_global > nu]
  tau_detect   <- if (length(post_global) > 0L) post_global[1L] else Inf

  # Per-stream first alarm and first post-change alarm
  stream_first <- vapply(stream_alarm_times, function(al)
    if (length(al) > 0L) al[1L] else Inf, numeric(1))
  stream_post  <- vapply(stream_alarm_times, function(al) {
    post <- al[al > nu]; if (length(post) > 0L) post[1L] else Inf
  }, numeric(1))
  stream_delays <- ifelse(is.finite(stream_post), stream_post - nu, NA_real_)
  stream_fa     <- vapply(stream_alarm_times,
                          function(al) any(al <= nu), logical(1))

  list(
    global_tau    = first_global,
    global_fa     = any(all_global <= nu),
    global_delay  = if (is.finite(tau_detect)) tau_detect - nu else NA_real_,
    stream_taus   = stream_first,
    stream_fa     = stream_fa,
    stream_delays = stream_delays
  )
}

# ── 4. Main simulation loop ──────────────────────────────────────────────────

# Pre-build the per-stream column names (K_MAX streams, ARL + CAD each)
stream_arl_names <- paste0("stream_", seq_len(K_MAX), "_ARL")
stream_cad_names <- paste0("stream_", seq_len(K_MAX), "_CAD")

results <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  row <- param_grid[i, ]
  K   <- as.integer(row$K)
  p   <- as.integer(row$p)
  nu  <- as.integer(row$nu)

  phi_sc <- AR_COEFS[[as.character(p)]]
  mu0    <- rep(MU_PRE, K)
  eta    <- rep(MU_PRE, K)

  mu1 <- if (isTRUE(row$sparse)) {
    m <- mu0; m[1L] <- mu0[1L] + row$delta_norm; m
  } else {
    mu0 + row$delta_norm / sqrt(K)
  }

  sig2 <- SIGMA_INN^2
  Sig  <- if (isTRUE(row$independent) || K == 1L) {
    diag(sig2, K)
  } else {
    S <- diag(sig2, K)
    S[upper.tri(S)] <- RHO_CORR * sig2
    S[lower.tri(S)] <- RHO_CORR * sig2
    S
  }

  reps <- vector("list", N_REP)
  for (r in seq_len(N_REP)) {
    set.seed(20000L * i + r)
    reps[[r]] <- tryCatch(
      .run_rep_localized(
        N           = N_OBS,
        K           = K,
        nu          = nu,
        mu0         = mu0,
        mu1         = mu1,
        phi_scalars = phi_sc,
        sigma_mat   = Sig,
        eta_vec     = eta,
        alpha       = row$alpha
      ),
      error = function(e) {
        list(global_tau    = NA_real_,
             global_fa     = NA,
             global_delay  = NA_real_,
             stream_taus   = rep(NA_real_, K),
             stream_fa     = rep(NA,       K),
             stream_delays = rep(NA_real_, K))
      }
    )
  }

  # ── Aggregate global results ──
  g_stops  <- vapply(reps, `[[`, numeric(1), "global_tau")
  g_delays <- vapply(reps, `[[`, numeric(1), "global_delay")
  g_fa     <- vapply(reps, `[[`, logical(1), "global_fa")

  # ── Aggregate per-stream results (wide format, NA for k > K) ──
  s_arl <- setNames(rep(NA_real_, K_MAX), stream_arl_names)
  s_cad <- setNames(rep(NA_real_, K_MAX), stream_cad_names)
  for (k in seq_len(K)) {
    k_stops  <- vapply(reps, function(r) r$stream_taus[k],   numeric(1))
    k_delays <- vapply(reps, function(r) r$stream_delays[k], numeric(1))
    s_arl[k] <- mean(k_stops)
    s_cad[k] <- mean(k_delays, na.rm = TRUE)
  }

  cond_row <- data.frame(
    p           = p,
    K           = K,
    nu          = nu,
    delta       = row$delta_norm,
    sparse      = row$sparse,
    alpha       = row$alpha,
    independent = row$independent,
    ARL         = mean(g_stops),
    CAD         = mean(g_delays, na.rm = TRUE),
    PFA         = mean(g_fa,     na.rm = TRUE),
    power       = mean(!is.na(g_delays)),
    stringsAsFactors = FALSE
  )
  cond_row <- cbind(cond_row, as.data.frame(t(s_arl)), as.data.frame(t(s_cad)))

  results[[i]] <- cond_row

  if (i %% 5L == 0L)
    message(sprintf("[%d / %d] p=%d  K=%d  nu=%d  delta=%.3f",
                    i, nrow(param_grid), p, K, nu, row$delta_norm))
}

results <- do.call(rbind, results)
rownames(results) <- NULL

# ── 5. Save ──────────────────────────────────────────────────────────────────

.pkg_root <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "DESCRIPTION")) && d != dirname(d)) d <- dirname(d)
  d
})
out_dir <- file.path(.pkg_root, "simulations", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(results, file.path(out_dir, "localized_sim_results.csv"), row.names = FALSE)
message("Saved to ", file.path(out_dir, "localized_sim_results.csv"))
