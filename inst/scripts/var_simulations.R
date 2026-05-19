# var_simulations.R
# Simulation study following Section 6 of Spertus et al.
#
# Three detector types are compared on K-stream VAR(p) data:
#   oracle    — GaussianVARModel with the true pre/post parameters
#   misspec   — GaussianVARModel with post-mean halved (mu1 replaced by (mu0+mu1)/2)
#   bounded   — BoundedModel per stream, combined via average / product / UP
#
# Innovation SD is set to sigma_inn = 0.06 so that the stationary distribution
# keeps data approximately in [0.06, 0.53] ⊂ [0,1], satisfying the bounded-TSM
# assumption.  (Worst-case AR(2) stationary SD ≈ 0.08; 3-sigma range from
# mu_pre = 0.3 is [0.06, 0.54].)
#
# ── Timing estimate (single core) ───────────────────────────────────────────
#   Quick mode (RUN_FULL = FALSE): ~5–10 min  [N=1000, 100 reps, 20 conditions]
#   Full  mode (RUN_FULL = TRUE) : ~2–4 days  [N=3000, 500 reps, ~17K conditions]
#
# compute_increments for both BoundedModel and GaussianVARModel is vectorised
# (O(N) and O(N*K^2) respectively), so per-replicate cost is negligible.
# The full grid is slow due to sheer condition count (~17K x 5 detector variants
# x 500 reps); use parallel::mclapply or future.apply on the outer loop.
# ────────────────────────────────────────────────────────────────────────────

library(multichangepoints)
library(mvtnorm)   # already imported by multichangepoints

# ── 0. Mode switch ──────────────────────────────────────────────────────────

RUN_FULL <- FALSE          # set TRUE for the full Section 6 grid

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

quick_params <- list(
  p           = 0L,
  K           = c(1L, 2L),
  nu          = c(50L, 100L),
  delta_norm  = seq(0.05, 0.20, length.out = 5),
  sparse      = FALSE,
  alpha       = 0.05,
  criterion   = "ARL",
  independent = TRUE
)

full_params <- list(
  p           = c(0L, 1L, 2L),
  K           = c(1L, 2L, 10L),
  nu          = c(10L, 50L, 100L, 1000L),
  delta_norm  = seq(0.01, 0.20, length.out = 20),
  sparse      = c(FALSE, TRUE),
  alpha       = c(0.001, 0.01, 0.1),
  criterion   = c("ARL", "PFA"),
  independent = c(TRUE, FALSE)
)

param_grid <- do.call(expand.grid,
  c(if (RUN_FULL) full_params else quick_params,
    list(stringsAsFactors = FALSE)))

N_OBS <- if (RUN_FULL) 3000L else 1000L
N_REP <- if (RUN_FULL) 500L  else 100L

# ── 3. Data generation ──────────────────────────────────────────────────────

# Generate an N_total x K VAR(p) path with a mean shift at time nu.
# The post-change segment is initialised from the end of the pre-change segment.
# phi_scalars: numeric vector of length p (applied equally to all K streams)
.generate_var_path <- function(N_total, K, mu_pre, mu_post,
                                nu, phi_scalars, sigma_mat) {
  p     <- length(phi_scalars)
  N_pad <- N_total + p                    # first p rows are burn-in
  x     <- matrix(0, nrow = N_pad, ncol = K)
  for (i in seq_len(p)) x[i, ] <- mu_pre # initialise at pre-change mean

  for (t in (p + 1L):N_pad) {
    obs_idx  <- t - p                     # observation index 1..N_total
    mu_now   <- if (obs_idx <= nu) mu_pre else mu_post
    cond_mu  <- mu_now
    for (j in seq_len(p))
      cond_mu <- cond_mu + phi_scalars[[j]] * (x[t - j, ] - mu_now)
    x[t, ] <- cond_mu + as.numeric(rmvnorm(1, sigma = sigma_mat))
  }
  x[(p + 1L):N_pad, , drop = FALSE]
}

# ── 4. Single-replicate runner ───────────────────────────────────────────────

# Returns a named list: stopping_time, false_alarm, delay
.run_rep <- function(N, K, nu, mu_pre_vec, mu_post_vec,
                     phi_scalars, sigma_mat, eta_vec,
                     alpha, criterion, detector_type, combine_method) {
  x <- .generate_var_path(N, K, mu_pre_vec, mu_post_vec,
                           nu, phi_scalars, sigma_mat)

  # ── compute per-stream or joint increments ──
  if (detector_type %in% c("oracle", "misspec")) {
    mu_use  <- if (detector_type == "oracle") mu_post_vec
               else (mu_pre_vec + mu_post_vec) / 2
    Phi_lst <- lapply(phi_scalars, function(ph) ph * diag(K))
    m       <- GaussianVARModel(
      Phi_pre   = Phi_lst, Sigma_pre = sigma_mat, mean_pre  = mu_pre_vec,
      Phi_post  = Phi_lst, Sigma_post = sigma_mat, mean_post = mu_use
    )
    inc <- compute_increments(TSM(m), x, log = TRUE)

  } else {
    # BoundedModel: K separate streams, then combine
    inc_mat <- matrix(NA_real_, nrow = N, ncol = K)
    for (k in seq_len(K))
      inc_mat[, k] <- compute_increments(TSM(BoundedModel(eta = eta_vec[k])),
                                         x[, k], log = TRUE)

    inc <- if (K == 1L) {
      inc_mat[, 1L]
    } else if (combine_method == "average") {
      combine_streams(AverageCombiner(), inc_mat, log = TRUE)
    } else if (combine_method == "product") {
      combine_streams(ProductCombiner(), inc_mat, log = TRUE)
    } else {
      combine_streams(UniversalPortfolioCombiner(mode = "sparse"), inc_mat, log = TRUE)
    }
  }

  # ── run detector ──
  det <- ShiryaevRobertsDetector(alpha = alpha, criterion = criterion)
  out <- run_detector(det, inc, log = TRUE)
  tau <- out$stopping_time

  list(
    stopping_time   = tau,
    false_alarm     = is.finite(tau) && tau <= nu,
    delay           = if (is.finite(tau) && tau > nu) tau - nu else NA_real_
  )
}

# ── 5. Main simulation loop ──────────────────────────────────────────────────

results <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  row  <- param_grid[i, ]
  K    <- as.integer(row$K)
  p    <- as.integer(row$p)
  nu   <- as.integer(row$nu)

  phi_sc    <- AR_COEFS[[as.character(p)]]
  mu0       <- rep(MU_PRE, K)
  eta       <- rep(MU_PRE, K)   # null mean = pre-change mean

  # Post-change mean vector
  if (isTRUE(row$sparse)) {
    mu1       <- mu0
    mu1[1L]   <- mu0[1L] + row$delta_norm
  } else {
    mu1 <- mu0 + row$delta_norm / sqrt(K)
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
    list(type = "misspec", combine = "joint")
  )
  # Bounded combiners: average always; product only under independence; UP when K > 1
  det_variants <- c(det_variants,
    list(list(type = "bounded", combine = "average")))
  if (K > 1L && (isTRUE(row$independent) || K == 1L))
    det_variants <- c(det_variants, list(list(type = "bounded", combine = "product")))
  if (K > 1L)
    det_variants <- c(det_variants, list(list(type = "bounded", combine = "up")))

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
          mu_pre_vec     = mu0,
          mu_post_vec    = mu1,
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

    stops   <- vapply(reps, `[[`, numeric(1), "stopping_time")
    delays  <- vapply(reps, `[[`, numeric(1), "delay")
    fa      <- vapply(reps, `[[`, logical(1), "false_alarm")

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
      ARL         = mean(stops, na.rm = TRUE),
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

# ── 6. Save ──────────────────────────────────────────────────────────────────

out_dir  <- "simulation_output"
if (!dir.exists(out_dir)) dir.create(out_dir)
saveRDS(results, file.path(out_dir, "var_sim_results.rds"))
write.csv(results, file.path(out_dir, "var_sim_results.csv"), row.names = FALSE)
message("Saved to ", out_dir)
