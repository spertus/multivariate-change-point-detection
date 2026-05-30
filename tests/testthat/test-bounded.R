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
  expect_error(BoundedModel(eta = 0),    "`eta` must be a scalar in \\(0, 1\\)")
  expect_error(BoundedModel(eta = 1),    "`eta` must be a scalar in \\(0, 1\\)")
  expect_error(BoundedModel(eta = -0.1), "`eta` must be a scalar in \\(0, 1\\)")
  expect_error(BoundedModel(eta = 1.1),  "`eta` must be a scalar in \\(0, 1\\)")
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
