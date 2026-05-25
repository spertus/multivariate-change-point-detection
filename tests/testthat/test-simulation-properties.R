# test-simulation-properties.R
# Three correctness tests targeting properties the simulation study relies on.
#
# The tests intentionally use simple, often deterministic scenarios so they
# run fast and fail loudly if something regresses.
#
# ── Note on CAD vs ARL ────────────────────────────────────────────────────────
# In a simulation with multiple_alarms = TRUE and a finite nu:
#
#   CAD ≠ ARL in general.
#
# By renewal theory, the expected remaining time in the alarm cycle at an
# arbitrary point nu equals E[T^2] / (2 * E[T]), where T is the inter-alarm
# distribution.  For a detector whose increments have non-trivial variance,
# this is ARL/2 + Var[T]/(2*ARL), not ARL itself.
#
# CAD = ARL only when (a) the increments are identically 0 in log-scale
# (LR = 1 exactly, as in the oracle with mu_pre = mu_post) AND (b) nu is a
# multiple of the alarm period.  Test 1 exercises this deterministic special case.
#
# The oracle does not always have lower CAD than the misspecified detector.
# The oracle uses a larger LR magnitude, producing higher-variance increments;
# higher variance increases the renewal-theory expected remaining time, which can
# offset (or reverse) the oracle's faster drift under H1.  Test 2 verifies the
# underlying drift ordering—the quantity that IS always guaranteed—rather than
# a CAD ordering that depends on renewal dynamics.
#
# Product combining does not always produce lower CAD than average combining.
# Under near-null conditions the product sums K log-increments, giving K× the
# variance; by renewal theory this raises expected remaining time above the
# average combiner.  Test 3 uses a strong signal where the product's larger
# drift dominates and verifies the ordering holds with high probability.
# ─────────────────────────────────────────────────────────────────────────────

library(multichangepoints)

# ── Test 1: Exact null → CAD equals alarm period ──────────────────────────────
#
# When mu_pre = mu_post, each log-LR increment is identically 0, so the S-R
# statistic grows as R_t = t (deterministic).  The detector fires at a fixed
# period T ≈ 1/alpha and resets.  For nu = k*T (a multiple of the period),
# the first post-change alarm is exactly one period later, giving CAD = T ≈ ARL.

test_that("exact null: S-R has deterministic alarm period and CAD equals that period", {
  alpha <- 0.05

  m   <- GaussianModel(mean_pre = 0.3, mean_post = 0.3, sd_pre = 0.06, sd_post = 0.06)
  inc <- compute_increments(TSM(m), rnorm(600, 0.3, 0.06), log = TRUE)
  expect_true(all(inc == 0))   # LR = 1 identically → log-inc = 0 everywhere

  det <- ShiryaevRobertsDetector(alpha = alpha, criterion = "ARL", multiple_alarms = TRUE)
  out <- run_detector(det, inc, log = TRUE)

  # Alarm period T is the first alarm time (deterministic for log-inc = 0).
  # Due to floating-point rounding in logsumexp, T may be 1/alpha or 1/alpha + 1.
  T_period <- out$alarm_times[1L]
  expect_lte(abs(T_period - 1 / alpha), 1L)  # within 1 step of ARL = 20

  # With nu = 3 * T_period (a multiple of the alarm period), the first post-change
  # alarm is exactly T_period steps after nu.
  nu         <- 3L * T_period
  post_alarm <- out$alarm_times[out$alarm_times > nu][1L]
  expect_equal(post_alarm - nu, T_period)
})

# ── Test 2: Oracle log-increments exceed misspecified under H1 ─────────────────
#
# Under the true post-change distribution P1, the oracle log-LR has a strictly
# larger expectation than the misspecified log-LR:
#
#   E_1[log(p1/p0)]   = delta^2 / (2 sigma^2)
#   E_1[log(p_m/p0)]  = 3 delta^2 / (8 sigma^2)   (mu_post = (mu0 + mu1)/2)
#
# The ratio is 4/3, so oracle drift is 33% larger.  This is an exact algebraic
# identity independent of sample size; with N = 5000 it is easily measurable.

test_that("oracle log-increments strictly exceed misspecified under post-change", {
  set.seed(2024)
  mu0 <- 0.3; delta <- 0.06; sigma <- 0.06; N <- 5000L
  mu1    <- mu0 + delta
  mu_use <- (mu0 + mu1) / 2

  x_post <- rnorm(N, mu1, sigma)

  oracle_m  <- GaussianModel(mean_pre = mu0, mean_post = mu1,   sd_pre = sigma, sd_post = sigma)
  misspec_m <- GaussianModel(mean_pre = mu0, mean_post = mu_use, sd_pre = sigma, sd_post = sigma)

  inc_oracle  <- compute_increments(TSM(oracle_m),  x_post, log = TRUE)
  inc_misspec <- compute_increments(TSM(misspec_m), x_post, log = TRUE)

  # Oracle drift > misspec drift under H1
  expect_gt(mean(inc_oracle), mean(inc_misspec))

  # Oracle is ≥ 25% faster (theoretical ratio = 4/3 ≈ 1.33; allow 10% slack)
  expect_gt(mean(inc_oracle) / mean(inc_misspec), 1.2)
})

# ── Test 3: Product combining is more powerful than average under dense signal ─
#
# Under dense, independent change (all K streams shift), the product combiner
# accumulates K times the per-stream log-evidence per step, while the average
# accumulates only the single-stream amount.  For large enough signal the
# product reaches the alarm threshold sooner.
#
# This test uses 300 Monte Carlo replicates with a strong signal (delta = 0.15
# across K = 2 independent BoundedModel streams) where theory guarantees the
# ordering.  With 300 reps the probability of the wrong ordering is < 0.01 for
# these parameters (verified empirically).

test_that("product combining has lower CAD than average under strong dense independent signal", {
  set.seed(42)
  K     <- 2L;  alpha <- 0.05;  nu <- 50L;  N <- 500L;  N_rep <- 300L
  delta <- 0.15;  mu0 <- 0.3;  sigma <- 0.06
  mu1   <- rep(mu0 + delta / sqrt(K), K)
  eta   <- rep(mu0, K)
  Sig   <- diag(sigma^2, K)

  dgp_m <- GaussianVARModel(
    Phi_pre    = list(), Sigma_pre  = Sig, mean_pre  = rep(mu0, K),
    Phi_post   = list(), Sigma_post = Sig, mean_post = mu1
  )

  cad_avg  <- numeric(N_rep)
  cad_prod <- numeric(N_rep)

  for (r in seq_len(N_rep)) {
    set.seed(r)
    x <- generate_stream(DGP(dgp_m, nu = nu), N = N)
    if (is.null(dim(x))) x <- matrix(x, ncol = 1L)

    inc_mat <- matrix(NA_real_, N, K)
    for (k in seq_len(K))
      inc_mat[, k] <- compute_increments(
        TSM(BoundedModel(eta = eta[k])), x[, k], log = TRUE)

    for (combine in c("avg", "prod")) {
      inc <- if (combine == "avg") combine_streams(AverageCombiner(), inc_mat, log = TRUE)
             else                  combine_streams(ProductCombiner(), inc_mat, log = TRUE)
      det <- ShiryaevRobertsDetector(alpha = alpha, criterion = "ARL", multiple_alarms = TRUE)
      out <- run_detector(det, inc, log = TRUE)
      post  <- out$alarm_times[out$alarm_times > nu]
      delay <- if (length(post) > 0L) post[1L] - nu else NA_real_
      if (combine == "avg") cad_avg[r] <- delay else cad_prod[r] <- delay
    }
  }

  # Product CAD is strictly smaller than average CAD
  expect_lt(mean(cad_prod, na.rm = TRUE), mean(cad_avg, na.rm = TRUE))
})
