# localization_graph_simulations.R
# Simulation study for graph-structured localization (Sections 4.3, 5.3 of
# Spertus et al. 2026).
#
# K = 6 streams sit on a ProximityGraph. Under the propagation-on-a-graph
# configuration model of Section 4.3, a change originates at a source node at
# time nu0 and spreads outward: nu_k = nu0 + scale * d_G(source, k). We compare
# an adaptive LocalizedDetector (Section 5.3) that spends its allowance
# graph-informed (zeta > 0, correctly-specified kernel bandwidth) against the
# uniform e-d-Bonferroni baseline it nests (zeta = 0), across three graph
# topologies and fast/slow propagation, both under the global null (no change)
# and under the propagation-configuration alternative. All detectors control
# the global ARL to level alpha = 0.001; by Theorem 2 the same spending
# allowance also controls the error over patience of the associated monitoring
# rule, so no separate EOP-specific machinery is needed.
#
# Output columns:
#   Condition: graph_type, speed, zeta, has_change, K, alpha, N, nu0, scale, n_rep
#   Type-I diagnostics (meaningful mainly when has_change = FALSE):
#     ARL          = mean global stopping time (capped at N+1)
#     frac_alarmed = empirical P(any stream alarms within the horizon)
#   Localization / efficiency (meaningful only when has_change = TRUE):
#     CAD_total           = mean total conditional average delay (Section 2.4):
#                            sum over truly-changing streams of (D_k - nu_k),
#                            averaged over reps where that stream did not fire
#                            before its own change-point
#     mean_frac_detected  = average fraction of truly-changing streams (those
#                            with nu_k <= N) that were validly detected by N
#     detect_frac_early/mid/late = fraction of truly-changing streams already
#                            flagged (R_tk > 1/alpha) at three checkpoints after
#                            nu0 — the accumulating-detections curve
#
# ── Timing (measured, 8 cores via parallel::mclapply over conditions) ──────
#   Quick mode (RUN_FULL = FALSE): ~7 sec   [N=800,  50 reps,  24 conditions]
#   Full  mode (RUN_FULL = TRUE) : ~2.7 min [N=3000, 300 reps, 24 conditions]
#   Grid: graph_type(3) x speed(2) x zeta(2) x has_change(2) = 24 conditions.
# ─────────────────────────────────────────────────────────────────────────────

library(multichangepoints)

# ── 0. Mode switch ──────────────────────────────────────────────────────────

RUN_FULL <- TRUE          # set TRUE for the full grid

# ── 1. Fixed simulation constants ───────────────────────────────────────────

K         <- 6L
ALPHA     <- 0.001
MU_PRE    <- 0.3
SIGMA_INN <- 0.06          # keeps data in ~[0.06, 0.54] subset of [0,1]
DELTA     <- 0.06          # post-change mean shift applied at each stream's own nu_k

N_OBS <- if (RUN_FULL) 3000L else 800L
NU0   <- if (RUN_FULL) 500L  else 200L
N_REP <- if (RUN_FULL) 300L  else 50L

REMAIN     <- N_OBS - NU0
SCALE_FAST <- REMAIN / 50    # all reachable nodes change shortly after nu0 (dense)
SCALE_SLOW <- REMAIN * 0.7   # only the source and its immediate (1-hop) neighbors
                             # change within the horizon; 2+ hops away never change
                             # (sparse) — the regime where graph-informed spending
                             # helps most, since boosting a near neighbor of an
                             # active node necessarily comes at the expense of
                             # farther nodes that would otherwise share the budget

# Kernel bandwidth (Section 5.3) operates on raw graph distance d_G, which is
# in hop units here (unit edge weights) — NOT on the time-propagation `scale`
# above. A "correctly specified" kernel matches the graph's hop structure
# (how many neighbors are near/far), independent of how fast the change
# happens to propagate in time; using `scale` itself here would make the
# kernel nearly flat under slow propagation (bandwidth >> 1-5 hop distances).
KERNEL_XI <- 1.5

# ── 2. Graph builders (unit edge weights) ───────────────────────────────────

hub_spoke_graph <- function(K) {
  W <- matrix(0, K, K)
  for (k in 2:K) W[1, k] <- W[k, 1] <- 1
  ProximityGraph(W)
}

fully_connected_graph <- function(K) {
  W <- matrix(1, K, K)
  diag(W) <- 0
  ProximityGraph(W)
}

linear_graph <- function(K) {
  W <- matrix(0, K, K)
  for (k in seq_len(K - 1L)) W[k, k + 1L] <- W[k + 1L, k] <- 1
  ProximityGraph(W)
}

GRAPH_BUILDERS <- list(
  hub_spoke       = hub_spoke_graph,
  fully_connected = fully_connected_graph,
  linear          = linear_graph
)

# Function: pick_source
# purpose: source node for the propagation model — fixed (hub / chain end) for
#          hub-spoke and linear graphs; uniformly random each replicate for the
#          fully-connected graph (per spec: "a random node changes").
pick_source <- function(graph_type, K) {
  if (graph_type == "fully_connected") sample.int(K, 1L) else 1L
}

# ── 3. Parameter grid ────────────────────────────────────────────────────────

param_grid <- expand.grid(
  graph_type = names(GRAPH_BUILDERS),
  speed      = c("fast", "slow"),
  zeta       = c(0, 0.8),
  has_change = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)

# ── 4. Single-replicate runner ───────────────────────────────────────────────

# Simulates K independent per-stream Gaussian series, each shifting its own
# mean from MU_PRE to MU_PRE + DELTA at its own nu_k (Inf = never). The
# existing DGP/GaussianVARModel machinery only supports one scalar nu shared
# by every stream, so per-stream nu_k needs this small purpose-built generator
# rather than an extension of that core simulation class.
.simulate_streams_propagation <- function(N, K, nu_vec, mu0, mu1, sigma) {
  x <- matrix(NA_real_, nrow = N, ncol = K)
  for (k in seq_len(K)) {
    nuk <- nu_vec[k]
    if (!is.finite(nuk) || nuk >= N) {
      x[, k] <- rnorm(N, mean = mu0, sd = sigma)
    } else {
      n_pre  <- max(0L, floor(nuk))
      n_post <- N - n_pre
      x[, k] <- c(rnorm(n_pre,  mean = mu0, sd = sigma),
                  rnorm(n_post, mean = mu1, sd = sigma))
    }
  }
  x
}

# Runs one replicate: simulate data under the propagation configuration,
# compute per-stream BoundedModel/EWMABet log-increments, run one adaptive
# LocalizedDetector, and return the per-stream stopping times and statistic
# paths needed for the condition-level summary metrics.
.run_rep_localization <- function(N, K, graph, kernel, zeta, alpha,
                                  nu_vec, mu0, mu1, sigma) {
  x <- .simulate_streams_propagation(N, K, nu_vec, mu0, mu1, sigma)

  inc_mat <- matrix(NA_real_, nrow = N, ncol = K)
  for (k in seq_len(K)) {
    bm_k <- BoundedModel(eta = mu0, bets = EWMABet(rho = 0.1, mu_init = mu0))
    inc_mat[, k] <- compute_increments(TSM(bm_k), x[, k], log = TRUE)
  }

  ld  <- LocalizedDetector(K = K, alpha = alpha, criterion = "ARL",
                            proximity_graph = graph, kernel = kernel, zeta = zeta)
  out <- run_detector(ld, inc_mat, log = TRUE)

  stream_stops <- vapply(out$stream_results, function(r) r$stopping_time, numeric(1))
  stat_mat     <- vapply(out$stream_results, function(r) r$statistic, numeric(N))

  list(stream_stops = stream_stops, stat_mat = stat_mat)
}

# ── 5. Per-condition aggregation ─────────────────────────────────────────────

# Checkpoints (post-nu0) for the accumulating-detections curve, expressed as
# fractions of the remaining horizon so they scale automatically between
# quick and full mode.
LAG_FRACS <- c(early = 0.05, mid = 0.20, late = 0.50)

.run_condition <- function(i) {
  row        <- param_grid[i, ]
  graph_type <- row$graph_type
  speed      <- row$speed
  zeta       <- row$zeta
  has_change <- row$has_change

  graph <- GRAPH_BUILDERS[[graph_type]](K)
  scale <- if (speed == "fast") SCALE_FAST else SCALE_SLOW
  kernel <- exponential_kernel(xi = KERNEL_XI)
  lag_checkpoints <- NU0 + round(LAG_FRACS * REMAIN)
  log_thresh <- log(1 / ALPHA)   # stat_mat is on the log scale (.run_rep_localization
                                 # runs the detector with log = TRUE)

  g_stops   <- numeric(N_REP)
  cad_totals <- rep(NA_real_, N_REP)
  frac_detected <- rep(NA_real_, N_REP)
  detect_frac <- matrix(NA_real_, nrow = N_REP, ncol = length(LAG_FRACS),
                        dimnames = list(NULL, names(LAG_FRACS)))

  for (r in seq_len(N_REP)) {
    set.seed(30000L * i + r)
    source <- pick_source(graph_type, K)
    nu_vec <- if (has_change) propagation_nu(graph, source, NU0, scale) else rep(Inf, K)

    rep_out <- tryCatch(
      .run_rep_localization(N = N_OBS, K = K, graph = graph, kernel = kernel,
                            zeta = zeta, alpha = ALPHA, nu_vec = nu_vec,
                            mu0 = MU_PRE, mu1 = MU_PRE + DELTA, sigma = SIGMA_INN),
      error = function(e) NULL
    )
    if (is.null(rep_out)) {
      g_stops[r] <- NA_real_
      next
    }

    stops   <- rep_out$stream_stops
    stat_mat <- rep_out$stat_mat

    g_stops[r] <- min(stops)

    changing <- which(is.finite(nu_vec) & nu_vec <= N_OBS)
    if (length(changing) > 0L) {
      valid <- stops[changing] >= nu_vec[changing] & stops[changing] <= N_OBS
      cad_totals[r]    <- sum((stops[changing] - nu_vec[changing])[valid])
      frac_detected[r] <- mean(valid)

      for (j in seq_along(LAG_FRACS)) {
        t_check <- min(lag_checkpoints[j], N_OBS)
        flagged <- stat_mat[t_check, changing] > log_thresh
        detect_frac[r, j] <- mean(flagged)
      }
    }
  }

  cond_row <- data.frame(
    graph_type = graph_type,
    speed      = speed,
    zeta       = zeta,
    has_change = has_change,
    K          = K,
    alpha      = ALPHA,
    N          = N_OBS,
    nu0        = NU0,
    scale      = scale,
    n_rep      = N_REP,
    ARL              = mean(pmin(g_stops, N_OBS + 1), na.rm = TRUE),
    frac_alarmed     = mean(is.finite(g_stops)),
    CAD_total        = mean(cad_totals, na.rm = TRUE),
    mean_frac_detected = mean(frac_detected, na.rm = TRUE),
    detect_frac_early  = mean(detect_frac[, "early"], na.rm = TRUE),
    detect_frac_mid    = mean(detect_frac[, "mid"],   na.rm = TRUE),
    detect_frac_late   = mean(detect_frac[, "late"],  na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  message(sprintf("[%d / %d] graph=%s speed=%s zeta=%.1f has_change=%s",
                  i, nrow(param_grid), graph_type, speed, zeta, has_change))
  cond_row
}

# ── 6. Run ────────────────────────────────────────────────────────────────────

n_cores <- if (.Platform$OS.type == "windows") 1L else min(8L, nrow(param_grid))
results <- parallel::mclapply(seq_len(nrow(param_grid)), .run_condition,
                              mc.cores = n_cores)
results <- do.call(rbind, results)
rownames(results) <- NULL

# ── 7. Save ──────────────────────────────────────────────────────────────────

.pkg_root <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "DESCRIPTION")) && d != dirname(d)) d <- dirname(d)
  d
})
out_dir <- file.path(.pkg_root, "simulations", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(results, file.path(out_dir, "localization_graph_sim_results.csv"), row.names = FALSE)
message("Saved to ", file.path(out_dir, "localization_graph_sim_results.csv"))
