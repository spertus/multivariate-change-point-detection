# localization_graph_simulations.R
# Simulation study for graph-structured localization (Sections 4.3, 5.3 of
# Spertus et al. 2026).
#
# K in {6, 24} streams sit on a ProximityGraph (topologies defined in the
# sourced localization_graph_builders.R, shared with the figures so the
# graphs simulated and the graphs drawn can never drift apart). Under the
# propagation-on-a-graph configuration model of Section 4.3, a change
# originates at a source node at time nu0 and spreads outward:
# nu_k = nu0 + scale * d_G(source, k). We compare four spending strategies,
# both under the global null (no change) and under the propagation-
# configuration alternative:
#   oracle         — static allowance (Theorem 1) concentrating most of the
#                    budget, from t=1, equally across the streams known (only
#                    to the simulator) to truly change within the horizon.
#                    Not a real procedure — an efficiency ceiling.
#   uniform        — the uniform e-d-Bonferroni baseline (zeta = 0).
#   graph_correct  — adaptive spending (Section 5.3) on the TRUE proximity
#                    graph that generated the configuration.
#   graph_misspec  — adaptive spending on a graph whose edge weights are all
#                    scaled down by MISSPEC_FACTOR, so the detector perceives
#                    every pair of nodes as closer than truth and reacts as if
#                    the change propagates MISSPEC_FACTOR times faster than it
#                    really does (overestimates the speed of propagation).
# All detectors control the global ARL to level alpha = 0.001; by Theorem 2
# the same spending allowance also controls the error over patience of the
# associated monitoring rule, so no separate EOP-specific machinery is needed.
#
# Graph sizes: K = 6 (as before) and K = 24. hub_spoke scales structurally,
# not just numerically: at K = 6 it is a single hub + 5 spokes; at K = 24 it
# becomes 4 hubs x 5 spokes each (multi_hub_graph in the sourced builders
# file), with the hubs themselves either disconnected (4 separate star
# components -- a clustered-dependency-graph scenario where change can only
# ever reach one cluster) or chained in a line (hub_1--hub_2--hub_3--hub_4,
# through which change can eventually reach every cluster). This
# hub_connected factor only applies to (graph_type == "hub_spoke", K == 24);
# it is fixed at "n/a" everywhere else. linear and fully_connected scale to
# K = 24 with no structural change (just more nodes/edges).
#
# Output columns:
#   Condition: graph_type, speed, strategy, has_change, K, hub_connected,
#              alpha, N, nu0, scale, n_rep
#   Type-I diagnostics (meaningful mainly when has_change = FALSE):
#     ARL          = mean global stopping time (capped at N+1)
#     frac_alarmed = empirical P(any stream alarms within the horizon)
#   Localization / efficiency (meaningful only when has_change = TRUE):
#     CAD_total           = mean total conditional average delay (Section 2.4):
#                            sum over truly-changing streams of (D_k - nu_k),
#                            excluding streams that falsely fired before their
#                            own nu_k; streams never detected within the
#                            horizon are right-censored at N (not excluded --
#                            excluding them would reward giving up on hard
#                            cases), mirroring the ARL censoring convention
#     mean_frac_detected  = average fraction of truly-changing streams (those
#                            with nu_k <= N) that were validly detected by N
#     detect_frac_early/mid/late = fraction of truly-changing streams already
#                            flagged (R_tk > 1/alpha) at three checkpoints after
#                            nu0 — the accumulating-detections curve
#
# ── Timing (measured, 8 cores via parallel::mclapply over conditions) ──────
#   Grid: graph_type(3) x speed(2) x strategy(4) x has_change(2) x K(2), with
#   hub_spoke x K=24 additionally split by hub_connected(2) = 112 conditions.
#   Quick mode (RUN_FULL = FALSE): N=800,  50 reps  per condition
#   Full  mode (RUN_FULL = TRUE) : N=3000, 300 reps per condition
#   (see console output of the most recent run for measured wall-clock time)
# ─────────────────────────────────────────────────────────────────────────────

library(multichangepoints)

# ── -1. Shared graph-topology builders (also used by ssc-2026/figures.R) ───

.pkg_root <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "DESCRIPTION")) && d != dirname(d)) d <- dirname(d)
  d
})
source(file.path(.pkg_root, "simulations", "localization_graph_builders.R"))

# ── 0. Mode switch ──────────────────────────────────────────────────────────

RUN_FULL <- TRUE          # set TRUE for the full grid

# ── 1. Fixed simulation constants ───────────────────────────────────────────

ALPHA     <- 0.001
MU_PRE    <- 0.3
SIGMA_INN <- 0.06          # keeps data in ~[0.06, 0.54] subset of [0,1]
DELTA     <- 0.06          # post-change mean shift applied at each stream's own nu_k

N_OBS <- if (RUN_FULL) 3000L else 800L
NU0   <- if (RUN_FULL) 500L  else 200L
N_REP <- if (RUN_FULL) 300L  else 50L

REMAIN     <- N_OBS - NU0
SCALE_FAST <- REMAIN / 50    # all reachable nodes change shortly after nu0 (dense),
                             # for every topology

# "Slow" propagation is calibrated per topology, not globally: hub-and-spoke and
# fully-connected have NO graded hop structure among non-source nodes -- every
# non-source node is exactly 1 hop away, so a single scale can only put ALL of
# them on the same side of the horizon (all-changing or none-changing), never a
# genuine partial/sparse split. To get an actually sparse scenario for these two
# topologies (only the source changes within the horizon -- a real efficiency
# test for targeted spending), their slow scale is pushed past REMAIN so even
# the 1-hop nodes fall outside it. The linear chain has graded hop distances
# (1..K-1), so a single scale already yields a genuine partial split (source +
# its immediate neighbor change; farther nodes don't) without needing this.
SCALE_SLOW_STAR  <- REMAIN * 1.2   # hub-and-spoke, fully-connected: only source changes
SCALE_SLOW_CHAIN <- REMAIN * 0.7   # linear: source + immediate (1-hop) neighbor change

# Kernel bandwidth (Section 5.3) operates on raw graph distance d_G, which is
# in hop units here (unit edge weights) — NOT on the time-propagation `scale`
# above. A "correctly specified" kernel matches the graph's hop structure
# (how many neighbors are near/far), independent of how fast the change
# happens to propagate in time; using `scale` itself here would make the
# kernel nearly flat under slow propagation (bandwidth >> 1-5 hop distances).
KERNEL_XI <- 1.5

# ── 2. Graph builders ─────────────────────────────────────────────────────
# (build_graph, hub_spoke_graph, fully_connected_graph, linear_graph,
#  multi_hub_graph, pick_source -- sourced above from
#  localization_graph_builders.R, shared with ssc-2026/figures.R)

# Focus parameter shared by both adaptive-graph strategies (correct/misspec) —
# only the assumed graph differs between them, so zeta is held fixed for a
# clean, single-factor comparison.
ZETA <- 0.8

# Function: oracle_allowance
# purpose: static allowance (Theorem 1) that concentrates ORACLE_SHARE of the
#          budget, equally, across streams known -- only to the simulator, not
#          a real detector -- to truly change within the horizon; the rest is
#          split equally among the remaining streams. Falls back to uniform
#          when nothing (or everything) changes, matching the convention used
#          throughout the package for "no information to act on".
# inputs:
#   nu_vec = numeric length-K vector of true change-points (Inf = never)
#   N      = integer horizon
#   K      = integer number of streams
#   share  = numeric scalar in (0,1); total budget given to changing streams
# outputs:
#   numeric length-K allowance vector summing to 1
ORACLE_SHARE <- 0.95
oracle_allowance <- function(nu_vec, N, K, share = ORACLE_SHARE) {
  changing <- which(is.finite(nu_vec) & nu_vec <= N)
  if (length(changing) == 0L || length(changing) == K) return(rep(1 / K, K))
  alloc <- rep((1 - share) / (K - length(changing)), K)
  alloc[changing] <- share / length(changing)
  alloc
}

# Function: misspecify_graph
# purpose: scale down every edge weight of a ProximityGraph by `factor`, so
#          d_G computed on the returned graph is `factor` times shorter than
#          on the true graph -- the detector perceives nodes as closer than
#          they are and reacts as if the change propagates `factor` times
#          faster than it truly does ("overestimates the speed of propagation").
MISSPEC_FACTOR <- 3
misspecify_graph <- function(graph, factor = MISSPEC_FACTOR) {
  ProximityGraph(graph@W / factor)
}

# Function: build_detector
# purpose: construct the LocalizedDetector for one of the four spending
#          strategies compared in this simulation (see file header).
build_detector <- function(strategy, K, alpha, true_graph, misspec_graph, kernel,
                           nu_vec, N) {
  if (strategy == "oracle") {
    LocalizedDetector(K = K, alpha = alpha, criterion = "ARL",
                       allowance = oracle_allowance(nu_vec, N, K))
  } else if (strategy == "uniform") {
    LocalizedDetector(K = K, alpha = alpha, criterion = "ARL")
  } else if (strategy == "graph_correct") {
    LocalizedDetector(K = K, alpha = alpha, criterion = "ARL",
                       proximity_graph = true_graph, kernel = kernel, zeta = ZETA)
  } else if (strategy == "graph_misspec") {
    LocalizedDetector(K = K, alpha = alpha, criterion = "ARL",
                       proximity_graph = misspec_graph, kernel = kernel, zeta = ZETA)
  } else {
    stop("unknown strategy: ", strategy, call. = FALSE)
  }
}

# ── 3. Parameter grid ────────────────────────────────────────────────────────

# hub_connected only applies to (graph_type == "hub_spoke", K == 24); it is
# fixed at "n/a" everywhere else rather than needlessly duplicating every
# other condition across a factor that has no effect on it.
.base_grid <- expand.grid(
  graph_type = c("hub_spoke", "fully_connected", "linear"),
  speed      = c("fast", "slow"),
  strategy   = c("oracle", "uniform", "graph_correct", "graph_misspec"),
  has_change = c(TRUE, FALSE),
  K          = c(6L, 24L),
  stringsAsFactors = FALSE
)

.is_multi_hub <- .base_grid$graph_type == "hub_spoke" & .base_grid$K == 24L

.grid_plain <- .base_grid[!.is_multi_hub, ]
.grid_plain$hub_connected <- "n/a"

.grid_hub <- .base_grid[.is_multi_hub, ]
.grid_hub_disc <- .grid_hub; .grid_hub_disc$hub_connected <- "disconnected"
.grid_hub_line <- .grid_hub; .grid_hub_line$hub_connected <- "line"

param_grid <- rbind(.grid_plain, .grid_hub_disc, .grid_hub_line)
rownames(param_grid) <- NULL

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
# compute per-stream BoundedModel/EWMABet log-increments, run one
# LocalizedDetector under the requested strategy, and return the per-stream
# stopping times and statistic paths needed for the condition-level summary
# metrics.
.run_rep_localization <- function(N, K, strategy, true_graph, misspec_graph, kernel,
                                  alpha, nu_vec, mu0, mu1, sigma) {
  x <- .simulate_streams_propagation(N, K, nu_vec, mu0, mu1, sigma)

  inc_mat <- matrix(NA_real_, nrow = N, ncol = K)
  for (k in seq_len(K)) {
    bm_k <- BoundedModel(eta = mu0, bets = EWMABet(rho = 0.1, mu_init = mu0))
    inc_mat[, k] <- compute_increments(TSM(bm_k), x[, k], log = TRUE)
  }

  ld  <- build_detector(strategy, K, alpha, true_graph, misspec_graph, kernel, nu_vec, N)
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
  row           <- param_grid[i, ]
  graph_type    <- row$graph_type
  speed         <- row$speed
  strategy      <- row$strategy
  has_change    <- row$has_change
  K             <- row$K
  hub_connected <- row$hub_connected

  true_graph    <- build_graph(graph_type, K, hub_connected)
  misspec_graph <- misspecify_graph(true_graph)
  scale_slow <- if (graph_type == "linear") SCALE_SLOW_CHAIN else SCALE_SLOW_STAR
  scale  <- if (speed == "fast") SCALE_FAST else scale_slow
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
    nu_vec <- if (has_change) propagation_nu(true_graph, source, NU0, scale) else rep(Inf, K)

    rep_out <- tryCatch(
      .run_rep_localization(N = N_OBS, K = K, strategy = strategy,
                            true_graph = true_graph, misspec_graph = misspec_graph,
                            kernel = kernel, alpha = ALPHA, nu_vec = nu_vec,
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
      # A stream "counts" toward CAD unless it fired before its own change-point
      # (an invalid detection, per Section 2.4's CAD_k definition); streams that
      # never fire within the horizon are NOT excluded -- excluding them would
      # reward a strategy for silently giving up on hard-to-detect streams
      # (survivorship bias). Instead their delay is right-censored at N_OBS,
      # the same conservative lower-bound convention already used for ARL
      # (mean(pmin(stop_times, N+1))) elsewhere in this codebase.
      not_false_alarm <- stops[changing] >= nu_vec[changing]
      censored_delay  <- pmin(stops[changing], N_OBS) - nu_vec[changing]
      cad_totals[r]    <- sum(censored_delay[not_false_alarm])
      frac_detected[r] <- mean(stops[changing] <= N_OBS & not_false_alarm)

      for (j in seq_along(LAG_FRACS)) {
        t_check <- min(lag_checkpoints[j], N_OBS)
        flagged <- stat_mat[t_check, changing] > log_thresh
        detect_frac[r, j] <- mean(flagged)
      }
    }
  }

  cond_row <- data.frame(
    graph_type    = graph_type,
    speed         = speed,
    strategy      = strategy,
    has_change    = has_change,
    K             = K,
    hub_connected = hub_connected,
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

  message(sprintf("[%d / %d] graph=%s K=%d hub_conn=%s speed=%s strategy=%s has_change=%s",
                  i, nrow(param_grid), graph_type, K, hub_connected, speed, strategy, has_change))
  cond_row
}

# ── 6. Run ────────────────────────────────────────────────────────────────────

n_cores <- if (.Platform$OS.type == "windows") 1L else min(8L, nrow(param_grid))
results <- parallel::mclapply(seq_len(nrow(param_grid)), .run_condition,
                              mc.cores = n_cores)
results <- do.call(rbind, results)
rownames(results) <- NULL

# ── 7. Save ──────────────────────────────────────────────────────────────────
# (.pkg_root already computed above, before sourcing the shared graph builders)

out_dir <- file.path(.pkg_root, "simulations", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(results, file.path(out_dir, "localization_graph_sim_results.csv"), row.names = FALSE)
message("Saved to ", file.path(out_dir, "localization_graph_sim_results.csv"))
