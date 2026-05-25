# joint_simulations.R
# Simulation study following Section 6 of Spertus et al.
#
# Three detector types are compared on K-stream VAR(p) data:
#   oracle    — GaussianVARModel with the true pre/post parameters
#   misspec   — GaussianVARModel with post-mean shifted up by 0.1 (mu1 replaced by mu1 + 0.1)
#   bounded   — BoundedModel per stream, combined via average or product
#
# Innovation SD is set to sigma_inn = 0.06 so that the stationary distribution
# keeps data approximately in [0.06, 0.53] ⊂ [0,1], satisfying the bounded-TSM
# assumption.  (Worst-case AR(2) stationary SD ≈ 0.08; 3-sigma range from
# mu_pre = 0.3 is [0.06, 0.54].)
#
# ── Timing estimate (single core) ───────────────────────────────────────────
#   Quick mode (RUN_FULL = FALSE): ~16 sec  [N=1000, 100 reps, 10 conditions]
#   Full  mode (RUN_FULL = TRUE) : ~1.5 hrs [N=3000, 300 reps, 640 conditions]
#
# Full grid: p(2) x K(2) x nu(2) x delta(20) x sparse(2) x alpha(1)
#            x criterion(1) x independent(2) = 640 conditions,
#            each with 3–4 detector variants x 300 reps.
# Use parallel::mclapply or future.apply on the outer (condition) loop to
# reduce wall time.
# ────────────────────────────────────────────────────────────────────────────

library(multichangepoints)

# ── 0. Mode switch ──────────────────────────────────────────────────────────

RUN_FULL <- TRUE           # set TRUE for the full Section 6 grid

# ── 1. Fixed simulation constants ───────────────────────────────────────────

MU_PRE    <- 0.3           # pre-change long-run mean (same for all streams)
SIGMA_INN <- 0.06          # innovation SD; keeps stationary data in ~[0.06, 0.54]
RHO_CORR  <- 0.30          # cross-stream correlation coefficient for the correlated case
SPEND_P   <- 0.01          # geometric spending parameter for PFA criterion

# AR coefficients (no cross-stream causality; diagonal Phi matrices)
# phi vectors chosen so sum(phi) < 1 and stationary SD * 3 < MU_PRE
AR_COEFS <- list(
  "0" = list(),                      # IID
  "1" = list(0.4),                   # AR(1)
  "2" = list(0.3, 0.1)               # AR(2)
)

# ── 2. Parameter grids ──────────────────────────────────────────────────────
# delta_norm is on the log scale so that more resolution is allocated near 0.01
# (small signals) where power curves are steepest.

quick_params <- list(
  p           = 0L,
  K           = c(1L, 2L),
  nu          = 100L,
  delta_norm  = c(0, exp(seq(log(0.01), log(0.2), length.out = 4L))),
  sparse      = FALSE,
  alpha       = 0.0001,
  criterion   = "ARL",
  independent = TRUE
)

full_params <- list(
  p           = c(0L, 1L),
  K           = c(2L, 10L),
  nu          = c(10L, 1000L),
  delta_norm  = c(0, exp(seq(log(0.01), log(0.2), length.out = 19L))),
  sparse      = c(FALSE, TRUE),
  alpha       = 0.0001,
  criterion   = "ARL",
  independent = c(TRUE, FALSE)
)

param_grid <- do.call(expand.grid,
  c(if (RUN_FULL) full_params else quick_params,
    list(stringsAsFactors = FALSE)))

N_OBS <- if (RUN_FULL) 3000L else 1000L
N_REP <- if (RUN_FULL) 300L  else 100L

# ── 3. Single-replicate runner ───────────────────────────────────────────────

# Returns a named list: stopping_time (first alarm), false_alarm, delay
# The detector uses multiple_alarms=TRUE so it resets after false alarms and
# keeps accumulating evidence. stopping_time is the first alarm (for ARL),
# false_alarm flags any alarm before nu, and delay is the gap from nu to the
# first alarm after nu (NA if no post-change alarm occurs within N steps).
.run_rep <- function(N, K, nu, mu0, mu1, phi_scalars, sigma_mat, eta_vec,
                     alpha, criterion, detector_type, combine_method) {
  Phi_lst   <- lapply(phi_scalars, function(ph) ph * diag(K))
  dgp_model <- GaussianVARModel(
    Phi_pre   = Phi_lst, Sigma_pre  = sigma_mat, mean_pre  = mu0,
    Phi_post  = Phi_lst, Sigma_post = sigma_mat, mean_post = mu1
  )
  x <- generate_stream(DGP(dgp_model, nu = nu), N = N)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)

  # ── compute per-stream or joint increments ──
  if (detector_type %in% c("oracle", "misspec")) {
    mu_use <- if (detector_type == "oracle") mu1 else mu1 + 0.1
    m <- GaussianVARModel(
      Phi_pre   = Phi_lst, Sigma_pre  = sigma_mat, mean_pre  = mu0,
      Phi_post  = Phi_lst, Sigma_post = sigma_mat, mean_post = mu_use
    )
    inc <- compute_increments(TSM(m), x, log = TRUE)

  } else {
    # BoundedModel: K separate streams, then combine (average or product only)
    inc_mat <- matrix(NA_real_, nrow = N, ncol = K)
    for (k in seq_len(K))
      inc_mat[, k] <- compute_increments(TSM(BoundedModel(eta = eta_vec[k])),
                                         x[, k], log = TRUE)

    inc <- if (K == 1L) {
      inc_mat[, 1L]
    } else if (combine_method == "average") {
      combine_streams(AverageCombiner(), inc_mat, log = TRUE)
    } else {
      combine_streams(ProductCombiner(), inc_mat, log = TRUE)
    }
  }

  # ── run detector with resets after false alarms ──
  det <- ShiryaevRobertsDetector(alpha = alpha, criterion = criterion,
                                  multiple_alarms = TRUE)
  out <- run_detector(det, inc, log = TRUE)

  all_alarms  <- out$alarm_times
  first_alarm <- if (length(all_alarms) > 0L) all_alarms[1L] else Inf
  post_alarms <- all_alarms[all_alarms > nu]
  tau_detect  <- if (length(post_alarms) > 0L) post_alarms[1L] else Inf

  list(
    stopping_time = first_alarm,
    false_alarm   = any(all_alarms <= nu),
    delay         = if (is.finite(tau_detect)) tau_detect - nu else NA_real_
  )
}

# ── 4. Main simulation loop ──────────────────────────────────────────────────

results <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  row  <- param_grid[i, ]
  K    <- as.integer(row$K)
  p    <- as.integer(row$p)
  nu   <- as.integer(row$nu)

  phi_sc <- AR_COEFS[[as.character(p)]]
  mu0    <- rep(MU_PRE, K)
  eta    <- rep(MU_PRE, K)   # null mean = pre-change mean

  # Post-change mean vector
  mu1 <- if (isTRUE(row$sparse)) {
    m <- mu0; m[1L] <- mu0[1L] + row$delta_norm; m
  } else {
    mu0 + row$delta_norm / sqrt(K)
  }

  # Innovation covariance
  sig2 <- SIGMA_INN^2
  Sig  <- if (isTRUE(row$independent) || K == 1L) {
    diag(sig2, K)
  } else {
    S <- diag(sig2, K)
    S[upper.tri(S)] <- RHO_CORR * sig2
    S[lower.tri(S)] <- RHO_CORR * sig2
    S
  }

  # Detector variants to run for this condition
  det_variants <- list(
    list(type = "oracle",  combine = "joint"),
    list(type = "misspec", combine = "joint"),
    list(type = "bounded", combine = "average")
  )
  if (K > 1L && isTRUE(row$independent))
    det_variants <- c(det_variants, list(list(type = "bounded", combine = "product")))

  cond_rows <- vector("list", length(det_variants))

  for (d in seq_along(det_variants)) {
    dv   <- det_variants[[d]]
    reps <- vector("list", N_REP)

    for (r in seq_len(N_REP)) {
      set.seed(10000L * i + 100L * d + r)
      reps[[r]] <- tryCatch(
        .run_rep(
          N              = N_OBS,
          K              = K,
          nu             = nu,
          mu0            = mu0,
          mu1            = mu1,
          phi_scalars    = phi_sc,
          sigma_mat      = Sig,
          eta_vec        = eta,
          alpha          = row$alpha,
          criterion      = row$criterion,
          detector_type  = dv$type,
          combine_method = dv$combine
        ),
        error = function(e) {
          list(stopping_time = NA_real_, false_alarm = NA, delay = NA_real_)
        }
      )
    }

    stops  <- vapply(reps, `[[`, numeric(1),  "stopping_time")
    delays <- vapply(reps, `[[`, numeric(1),  "delay")
    fa     <- vapply(reps, `[[`, logical(1),  "false_alarm")

    cond_rows[[d]] <- data.frame(
      p           = p,
      K           = K,
      nu          = nu,
      delta       = row$delta_norm,
      sparse      = row$sparse,
      alpha       = row$alpha,
      criterion   = row$criterion,
      independent = row$independent,
      detector    = dv$type,
      combine     = dv$combine,
      ARL         = mean(stops,  na.rm = TRUE),
      CAD         = mean(delays, na.rm = TRUE),
      PFA         = mean(fa,     na.rm = TRUE),
      power       = mean(!is.na(delays)),
      stringsAsFactors = FALSE
    )
  }

  results[[i]] <- do.call(rbind, cond_rows)

  if (i %% 5L == 0L)
    message(sprintf("[%d / %d] p=%d  K=%d  nu=%d  delta=%.3f",
                    i, nrow(param_grid), p, K, nu, row$delta_norm))
}

results <- do.call(rbind, results)

# ── 5. Save ──────────────────────────────────────────────────────────────────

.pkg_root <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "DESCRIPTION")) && d != dirname(d)) d <- dirname(d)
  d
})
out_dir <- file.path(.pkg_root, "simulations", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(results, file.path(out_dir, "joint_sim_results.csv"), row.names = FALSE)
message("Saved to ", file.path(out_dir, "joint_sim_results.csv"))
