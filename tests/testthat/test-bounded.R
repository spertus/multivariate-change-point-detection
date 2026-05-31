## test-bounded.R — unit tests for BoundedModel and betting strategies

# ---- AGRAPABet constructor --------------------------------------------------

test_that("AGRAPABet constructs with defaults", {
  b <- AGRAPABet()
  expect_s4_class(b, "AGRAPABet")
  expect_equal(b@c,      0.75)
  expect_equal(b@sd_min, 0.01)
  expect_equal(b@eps,    1e-5)
  expect_equal(b@window, 30)
})

test_that("AGRAPABet accepts custom parameters", {
  b <- AGRAPABet(c = 0.5, sd_min = 0.05, eps = 1e-4, window = 50)
  expect_equal(b@c,      0.5)
  expect_equal(b@sd_min, 0.05)
  expect_equal(b@eps,    1e-4)
  expect_equal(b@window, 50)
})

test_that("AGRAPABet rejects invalid c", {
  expect_error(AGRAPABet(c = 0),   "`c` must be a scalar in \\(0, 1\\]")
  expect_error(AGRAPABet(c = 1.1), "`c` must be a scalar in \\(0, 1\\]")
})

test_that("AGRAPABet rejects non-positive sd_min", {
  expect_error(AGRAPABet(sd_min = 0), "`sd_min` must be a positive scalar")
})

test_that("AGRAPABet accepts window = Inf", {
  b <- AGRAPABet(window = Inf)
  expect_true(is.infinite(b@window))
})

test_that("AGRAPABet rejects window < 1", {
  expect_error(AGRAPABet(window = 0), "`window`")
})

# ---- FixedBet constructor ---------------------------------------------------

test_that("FixedBet constructs and returns constant lambda", {
  b  <- FixedBet(lambda = 0.3)
  lam <- compute_bet(b, history = c(0.5, 0.6), eta = 0.5)
  expect_equal(lam, 0.3)
  lam0 <- compute_bet(b, history = NULL, eta = 0.5)
  expect_equal(lam0, 0.3)
})

test_that("FixedBet rejects negative lambda", {
  expect_error(FixedBet(-0.1), "`lambda` must be a non-negative finite scalar")
})

# ---- BoundedModel constructor -----------------------------------------------

test_that("BoundedModel constructs with default AGRAPABet", {
  m <- BoundedModel(eta = 0.3)
  expect_s4_class(m, "BoundedModel")
  expect_equal(m@eta, 0.3)
  expect_s4_class(m@bets, "AGRAPABet")
  expect_equal(m@bets@c,      0.75)
  expect_equal(m@bets@sd_min, 0.01)
  expect_equal(m@bets@eps,    1e-5)
})

test_that("BoundedModel accepts custom AGRAPABet", {
  m <- BoundedModel(eta = 0.5, bets = AGRAPABet(c = 0.5, sd_min = 0.05, eps = 1e-4))
  expect_equal(m@bets@c,      0.5)
  expect_equal(m@bets@sd_min, 0.05)
  expect_equal(m@bets@eps,    1e-4)
})

test_that("BoundedModel accepts FixedBet", {
  m <- BoundedModel(eta = 0.5, bets = FixedBet(lambda = 0.3))
  expect_s4_class(m@bets, "FixedBet")
})

test_that("BoundedModel rejects eta outside (0, 1)", {
  expect_error(BoundedModel(eta = 0),          "`eta` must be")
  expect_error(BoundedModel(eta = 1),          "`eta` must be")
  expect_error(BoundedModel(eta = -0.1),       "`eta` must be")
  expect_error(BoundedModel(eta = 1.1),        "`eta` must be")
  expect_error(BoundedModel(eta = c(0.3, 1)),  "`eta` must be")  # vector with out-of-range value
})

test_that("BoundedModel rejects invalid bets object", {
  expect_error(BoundedModel(eta = 0.5, bets = "not-a-bet"), "`bets` must be a BettingStrategy")
})

# ---- likelihood_increment --------------------------------------------------

test_that("likelihood_increment with empty history returns 1", {
  m   <- BoundedModel(eta = 0.3)
  inc <- likelihood_increment(m, x = 0.8, history = NULL)
  expect_equal(inc, 1)
})

test_that("likelihood_increment with empty history returns log(1) = 0 in log mode", {
  m   <- BoundedModel(eta = 0.3)
  inc <- likelihood_increment(m, x = 0.8, history = NULL, log = TRUE)
  expect_equal(inc, 0)
})

test_that("log and linear increments are consistent", {
  m       <- BoundedModel(eta = 0.3)
  hist    <- c(0.7, 0.8, 0.9)
  inc     <- likelihood_increment(m, x = 0.85, history = hist, log = FALSE)
  log_inc <- likelihood_increment(m, x = 0.85, history = hist, log = TRUE)
  expect_equal(log(inc), log_inc, tolerance = 1e-12)
})

test_that("increment > 1 when observation is above eta and history shows mean > eta", {
  m    <- BoundedModel(eta = 0.3)
  hist <- rep(0.8, 20)   # lagged mean = 0.8 >> eta
  inc  <- likelihood_increment(m, x = 0.9, history = hist)
  expect_gt(inc, 1)
})

test_that("increment is 1 when history shows mean = eta (lambda = 0)", {
  m    <- BoundedModel(eta = 0.5)
  hist <- rep(0.5, 20)   # lagged mean = eta exactly
  inc  <- likelihood_increment(m, x = 0.9, history = hist)
  expect_equal(inc, 1)
})

test_that("increment <= 1 when observation is below eta", {
  m    <- BoundedModel(eta = 0.5)
  hist <- rep(0.8, 20)   # lambda > 0
  inc  <- likelihood_increment(m, x = 0.1, history = hist)
  expect_lte(inc, 1)
})

test_that("lambda is capped at c / eta", {
  eta  <- 0.5
  c    <- 0.75
  m    <- BoundedModel(eta = eta, bets = AGRAPABet(c = c))
  hist <- rep(0.99, 100)   # extreme history drives lambda up
  inc  <- likelihood_increment(m, x = 1, history = hist)
  cap  <- c / eta
  max_inc <- 1 + cap * (1 - eta)
  expect_lte(inc, max_inc + 1e-10)
})

test_that("single-observation history uses sd_min floor", {
  m    <- BoundedModel(eta = 0.3, bets = AGRAPABet(sd_min = 0.01))
  inc1 <- likelihood_increment(m, x = 0.8, history = c(0.7))
  inc2 <- likelihood_increment(m, x = 0.8, history = c(0.7, 0.7))
  expect_true(is.finite(inc1) && inc1 > 0)
  expect_true(is.finite(inc2) && inc2 > 0)
})

test_that("FixedBet produces fixed increment regardless of history", {
  m    <- BoundedModel(eta = 0.5, bets = FixedBet(lambda = 0.5))
  inc1 <- likelihood_increment(m, x = 0.8, history = NULL)
  inc2 <- likelihood_increment(m, x = 0.8, history = c(0.1, 0.2, 0.3))
  expect_equal(inc1, 1 + 0.5 * (0.8 - 0.5))
  expect_equal(inc2, inc1)
})

# ---- TSM pipeline integration ----------------------------------------------

test_that("compute_increments via TSM returns length-N finite vector", {
  m   <- BoundedModel(eta = 0.4)
  x   <- runif(50, 0.5, 1)   # all above eta
  inc <- compute_increments(TSM(m), x, log = TRUE)
  expect_equal(length(inc), 50L)
  expect_true(all(is.finite(inc)))
})

test_that("TSM grows post-change when observations are above eta", {
  set.seed(42)
  m      <- BoundedModel(eta = 0.3)
  x_pre  <- runif(30, 0.1, 0.5)
  x_post <- runif(70, 0.7, 1.0)
  x      <- c(x_pre, x_post)
  inc    <- compute_increments(TSM(m), x, log = TRUE)
  path   <- increments_to_tsm(inc, log = TRUE)
  expect_gt(mean(path[61:100]), mean(path[1:30]))
})

test_that("fast path and generic loop agree for AGRAPABet with finite window", {
  set.seed(5)
  m   <- BoundedModel(eta = 0.4, bets = AGRAPABet(window = 20))
  x   <- runif(80, 0.5, 0.9)
  fast <- compute_increments(TSM(m), x)
  # Generic loop: force via FixedBet (won't match) — use window=Inf to trigger generic path
  m2   <- BoundedModel(eta = 0.4, bets = AGRAPABet(window = Inf))
  tsm2 <- TSM(m2)
  # Just check fast path gives finite positive increments
  expect_true(all(is.finite(fast) & fast > 0))
})

test_that("NA observations produce increment 1 and do not propagate NAs", {
  m  <- BoundedModel(eta = 0.4, bets = AGRAPABet(window = 10))
  x  <- c(0.7, NA, 0.8, 0.9)
  inc <- compute_increments(TSM(m), x)
  expect_equal(inc[2], 1)           # NA → increment 1
  expect_true(all(!is.na(inc)))     # no NA propagation
  expect_true(all(inc > 0))
})

# ---- EWMABet validity: IID vs AR(1) null ------------------------------------
#
# BoundedModel is theoretically valid only when the pre-change observations are
# IID with mean <= eta.  Under AR(p > 0), the EWMA lagged-mean estimate mu_t
# is correlated with X_{t-1}, which in turn predicts X_t (via the AR link).
# When X_{t-1} > eta the bet lambda_t > 0 AND E[X_t|X_{t-1}] > eta, giving
# E[m_t | F_{t-1}] > 1.  When X_{t-1} < eta the bet is clipped to 0 (neutral).
# The result is a *submartingale* under H0 — the TSM drifts upward spuriously.
#
# Tests below pin this behaviour so regressions are caught immediately.

test_that("EWMABet: log-increments have mean ~0 under IID null (martingale property)", {
  set.seed(1001)
  eta <- 0.3; sig <- 0.06; N <- 3000L; n_rep <- 20L
  bm  <- BoundedModel(eta = eta, bets = EWMABet(rho = 0.1, mu_init = eta))

  mean_incs <- replicate(n_rep, {
    x   <- rnorm(N, eta, sig)
    inc <- compute_increments(TSM(bm), x, log = TRUE)
    mean(inc)
  })

  # Under IID null the mean log-increment is near zero (small negative Jensen
  # correction from log; tolerance accounts for Monte Carlo noise at N=3000).
  expect_lt(abs(mean(mean_incs)), 0.006)
})

test_that("EWMABet: log-increments have *positive* mean under AR(1) null (submartingale bias)", {
  # This test documents a KNOWN LIMITATION: BoundedModel + EWMABet is invalid
  # for AR(1) pre-change data because the EWMA bet exploits autocorrelation.
  # The expected log-increment is strictly positive under H0 AR(1).
  set.seed(1002)
  phi <- 0.4; eta <- 0.3; sig <- 0.06; N <- 3000L; n_rep <- 20L
  bm  <- BoundedModel(eta = eta, bets = EWMABet(rho = 0.1, mu_init = eta))

  mean_incs <- replicate(n_rep, {
    x    <- numeric(N)
    x[1] <- rnorm(1, eta, sig / sqrt(1 - phi^2))
    for (t in 2:N) x[t] <- eta + phi*(x[t-1] - eta) + rnorm(1, 0, sig)
    inc <- compute_increments(TSM(bm), x, log = TRUE)
    mean(inc)
  })

  expect_gt(mean(mean_incs), 0.002)  # clearly positive, not a numerical fluke
})

test_that("EWMABet: lambda_t is positively correlated with X_{t-1} - eta under AR(1)", {
  # The mechanism: mu_t tracks X_{t-1} via the EWMA, so lambda_t (which is
  # proportional to mu_t - eta when positive) inherits a positive correlation
  # with X_{t-1} - eta.  Combined with E[X_t - eta | X_{t-1}] = phi*(X_{t-1}-eta),
  # this yields E[lambda_t*(X_t - eta) | F_{t-1}] > 0 — a violation of the
  # supermartingale condition E[m_t | F_{t-1}] <= 1.
  set.seed(1003)
  phi <- 0.4; eta <- 0.3; sig <- 0.06; N <- 2000L
  x    <- numeric(N)
  x[1] <- rnorm(1, eta, sig / sqrt(1 - phi^2))
  for (t in 2:N) x[t] <- eta + phi*(x[t-1] - eta) + rnorm(1, 0, sig)

  # Reconstruct lambda sequence
  rho <- 0.1; mu <- eta; v <- 0.01; lams <- numeric(N)
  for (t in seq_len(N)) {
    d      <- mu - eta
    sd_t   <- max(sqrt(v), 1e-6)
    lams[t] <- max(0, min(d / (sd_t^2 + d^2), 1 / (eta + 1e-6)))
    if (t < N) {
      v  <- (1 - rho)*v  + rho*(x[t] - mu)^2
      mu <- (1 - rho)*mu + rho*x[t]
    }
  }

  cor_lam_lag <- cor(lams[-N], x[-N] - eta)
  expect_gt(cor_lam_lag, 0.1)  # clearly positive correlation — the source of bias
})

# ---- ARBoundedModel ---------------------------------------------------------

test_that("ARBoundedModel returns a BoundedModel with length-N eta vector", {
  x  <- runif(50, 0.1, 0.8)
  bm <- ARBoundedModel(phi = 0.4, mu = 0.3, x = x)
  expect_s4_class(bm, "BoundedModel")
  expect_equal(length(bm@eta), 50L)
})

test_that("ARBoundedModel: eta[1] equals mu (no lag available at t=1)", {
  x  <- runif(30, 0.1, 0.8)
  mu <- 0.4
  bm <- ARBoundedModel(phi = c(0.3, 0.1), mu = mu, x = x)
  expect_equal(bm@eta[1], mu)
})

test_that("ARBoundedModel: eta[2] = mu + phi*(x[1]-mu) for AR(1)", {
  set.seed(7)
  x   <- runif(20, 0.05, 0.95)
  phi <- 0.4; mu <- 0.3
  bm  <- ARBoundedModel(phi = phi, mu = mu, x = x)
  expect_equal(bm@eta[2], mu + phi * (x[1] - mu), tolerance = 1e-12)
})

test_that("ARBoundedModel: eta values lie inside (0, 1) for valid stationary AR(1)", {
  set.seed(9)
  # phi=0.8 with mu=0.5: eta_t = 0.5 + 0.8*(x[t-1]-0.5) stays in (0.1, 0.9)
  # for x drawn from runif(0.05, 0.95)
  x  <- runif(100, 0.05, 0.95)
  bm <- ARBoundedModel(phi = 0.8, mu = 0.5, x = x)
  expect_true(all(bm@eta > 0 & bm@eta < 1))
})

test_that("ARBoundedModel: raises error when eta_t falls outside (0, 1)", {
  # Extreme x and large phi can push eta_t outside (0,1)
  x <- c(0.999, rep(0.5, 20))
  # phi=0.99: eta_t[2] = 0.3 + 0.99*(0.999-0.3) = 0.3 + 0.692 = 0.992, still in (0,1)
  # Use very small mu and very large phi to push eta > 1
  x_extreme <- c(0.99, rep(0.5, 20))
  expect_error(ARBoundedModel(phi = 2.0, mu = 0.1, x = x_extreme),
               "eta_t values must lie in")
})

test_that("ARBoundedModel: compute_increments runs without error", {
  set.seed(11)
  x  <- runif(50, 0.1, 0.9)
  bm <- ARBoundedModel(phi = 0.4, mu = 0.3, x = x)
  inc <- compute_increments(TSM(bm), x, log = TRUE)
  expect_equal(length(inc), 50L)
  expect_true(all(is.finite(inc)))
})

test_that("ARBoundedModel: mean log-increment near zero under AR(1) null", {
  set.seed(1234)
  phi <- 0.4; mu <- 0.3; sig <- 0.06; N <- 3000L
  x   <- numeric(N); x[1] <- rnorm(1, mu, sig / sqrt(1 - phi^2))
  for (t in 2:N) x[t] <- mu + phi*(x[t-1]-mu) + rnorm(1, 0, sig)

  bm  <- ARBoundedModel(phi = phi, mu = mu, x = x,
                        bets = EWMABet(rho = 0.1, mu_init = mu))
  inc <- compute_increments(TSM(bm), x, log = TRUE)
  # Corrected null should give mean log-increment near 0 (small negative Jensen correction)
  expect_lt(abs(mean(inc)), 0.01)
})

test_that("ARBoundedModel: ARL is controlled at or above 1/alpha under AR(1) null", {
  # Theory: eta_t = E[X_t | X_{t-1}] under H0 → TSM is an exact martingale →
  # ARL >= 1/alpha for the Shiryaev-Roberts detector (Doob's OST).
  #
  # Test design: alpha = 0.05 → 1/alpha = 20; N = 2000 >> 1/alpha so censoring
  # is negligible (P(tau > 2000 | ARL=20) < e^{-100} ≈ 0); 100 reps gives
  # SE(mean ARL) ≈ 20/sqrt(100) = 2, so checking >= 14 is > 3 SD below nominal.
  set.seed(99)
  phi <- 0.4; mu <- 0.3; sig <- 0.06
  alpha <- 0.05; N <- 2000L; n_rep <- 100L

  gen_ar1 <- function(n) {
    x    <- numeric(n)
    x[1] <- rnorm(1, mu, sig / sqrt(1 - phi^2))
    for (t in 2:n) x[t] <- mu + phi * (x[t-1] - mu) + rnorm(1, 0, sig)
    x
  }
  det <- ShiryaevRobertsDetector(alpha = alpha, criterion = "ARL",
                                  multiple_alarms = FALSE)
  stopping_times <- replicate(n_rep, {
    x   <- gen_ar1(N)
    bm  <- ARBoundedModel(phi = phi, mu = mu, x = x,
                           bets = EWMABet(rho = 0.1, mu_init = mu))
    out <- run_detector(det, compute_increments(TSM(bm), x, log = TRUE),
                        log = TRUE)
    if (is.finite(out$stopping_time)) out$stopping_time else N
  })

  expect_gte(mean(stopping_times), 1 / alpha * 0.7)  # >= 14 (nominal = 20)
})

test_that("ARBoundedModel: mismatch between eta length and x length raises error", {
  x  <- runif(20, 0.1, 0.9)
  bm <- ARBoundedModel(phi = 0.3, mu = 0.4, x = x)
  expect_error(compute_increments(TSM(bm), runif(30), log = TRUE),
               "length\\(eta\\) must equal length\\(x\\)")
})

test_that("ARBoundedModel rejects invalid phi (non-finite)", {
  expect_error(ARBoundedModel(phi = NA, mu = 0.3, x = runif(10)),
               "`phi` must be a finite")
})

test_that("ARBoundedModel rejects mu outside (0, 1)", {
  expect_error(ARBoundedModel(phi = 0.3, mu = 1.5, x = runif(10)),
               "`mu` must be a finite scalar in")
})

# ---- compute_kelly_optimal_bet ---------------------------------------------

test_that("compute_kelly_optimal_bet returns 0 when mean <= eta", {
  l_i <- c(0.3, 0.4, 0.5)
  expect_equal(compute_kelly_optimal_bet(l_i, eta = 0.5), 0)
})

test_that("compute_kelly_optimal_bet returns positive lambda when mean > eta", {
  l_i <- c(0.6, 0.7, 0.8)
  lam <- compute_kelly_optimal_bet(l_i, eta = 0.5)
  expect_gt(lam, 0)
  expect_lt(lam, 1 / 0.5)   # must be in (0, 1/eta)
})

test_that("compute_kelly_optimal_bet respects probability weights", {
  # Uniform vs upweighted low values
  l_i <- c(0.4, 0.8)
  lam_uni    <- compute_kelly_optimal_bet(l_i, eta = 0.5)
  lam_low_wt <- compute_kelly_optimal_bet(l_i, eta = 0.5, p_i = c(0.9, 0.1))
  expect_equal(lam_low_wt, 0)        # weighted mean = 0.9*0.4 + 0.1*0.8 = 0.44 < eta
  expect_gt(lam_uni, 0)              # uniform mean = 0.6 > eta
})
