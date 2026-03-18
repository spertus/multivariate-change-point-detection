# ---- Bets constructors ----

test_that("FixedBets constructor validates inputs", {
  expect_error(FixedBets(c = 0), "`c`")
  expect_error(FixedBets(c = 1.1), "`c`")
  expect_s4_class(FixedBets(c = 0.5), "FixedBets")
})

test_that("AGRAPABets constructor validates inputs", {
  expect_error(AGRAPABets(c = 0), "`c`")
  expect_error(AGRAPABets(sd_min = -1), "`sd_min`")
  expect_error(AGRAPABets(eps = -1), "`eps`")
  expect_s4_class(AGRAPABets(), "AGRAPABets")
})

test_that("PredictablePluginBets constructor validates inputs", {
  expect_error(PredictablePluginBets(alpha = 0), "`alpha`")
  expect_error(PredictablePluginBets(alpha = 1), "`alpha`")
  expect_error(PredictablePluginBets(c = 1.5), "`c`")
  expect_s4_class(PredictablePluginBets(), "PredictablePluginBets")
})

# ---- compute_bets ----

test_that("FixedBets returns constant bets of correct length", {
  bets <- FixedBets(c = 0.6)
  x   <- runif(10)
  lam <- compute_bets(bets, x, eta = 0.5)
  expect_equal(lam, rep(0.6, 10))
})

test_that("AGRAPABets returns non-negative bets of correct length", {
  bets <- AGRAPABets()
  x   <- runif(30)
  lam <- compute_bets(bets, x, eta = 0.5)
  expect_equal(length(lam), 30L)
  expect_true(all(lam >= 0))
})

test_that("AGRAPABets bets are larger when eta is well below the data mean", {
  bets   <- AGRAPABets()
  x      <- rep(0.9, 30)
  lam_lo <- compute_bets(bets, x, eta = 0.3)  # big gap
  lam_hi <- compute_bets(bets, x, eta = 0.85) # small gap
  expect_true(mean(lam_lo) > mean(lam_hi))
})

test_that("AGRAPABets bets are zero when lagged mean equals eta", {
  # lag_mu_0 = 0.5 by default; so at t=1 with eta=0.5, bet is 0
  bets <- AGRAPABets()
  x    <- c(0.5, 0.5, 0.5)
  lam  <- compute_bets(bets, x, eta = 0.5)
  expect_equal(lam[1], 0)
})

test_that("PredictablePluginBets returns bets in [0, c]", {
  bets <- PredictablePluginBets(c = 0.8, alpha = 0.05)
  x   <- runif(50)
  lam <- compute_bets(bets, x, eta = 0.5)
  expect_equal(length(lam), 50L)
  expect_true(all(lam >= 0))
  expect_true(all(lam <= 0.8))
})

test_that("compute_bets returns length-0 vector for empty x", {
  for (bets in list(FixedBets(), AGRAPABets(), PredictablePluginBets())) {
    lam <- compute_bets(bets, numeric(0), eta = 0.5)
    expect_equal(length(lam), 0L)
  }
})

# ---- BettingTSM constructor ----

test_that("BettingTSM constructor validates inputs", {
  bets <- AGRAPABets()
  expect_error(BettingTSM(bets = "not_bets", eta = 0.5), "`bets`")
  expect_error(BettingTSM(bets = bets, eta = 0), "`eta`")
  expect_error(BettingTSM(bets = bets, eta = 1.1), "`eta`")
  expect_error(BettingTSM(bets = bets, eta = 0.5, initial = 0), "`initial`")
  expect_s4_class(BettingTSM(bets = bets, eta = 0.5), "BettingTSM")
})

# ---- compute_increments ----

test_that("compute_increments returns positive vector of correct length", {
  tsm <- BettingTSM(bets = AGRAPABets(), eta = 0.5)
  x   <- runif(40)
  inc <- compute_increments(tsm, x)
  expect_equal(length(inc), 40L)
  expect_true(all(inc > 0))
})

test_that("compute_increments log mode matches exp of log", {
  tsm     <- BettingTSM(bets = AGRAPABets(), eta = 0.5)
  x       <- runif(20)
  inc     <- compute_increments(tsm, x, log = FALSE)
  log_inc <- compute_increments(tsm, x, log = TRUE)
  expect_equal(exp(log_inc), inc, tolerance = 1e-12)
})

test_that("compute_increments with FixedBets matches manual formula", {
  bets <- FixedBets(c = 0.5)
  eta  <- 0.4
  tsm  <- BettingTSM(bets = bets, eta = eta)
  x    <- c(0.6, 0.2, 0.8)
  inc  <- compute_increments(tsm, x)
  expected <- pmax(1 + 0.5 * (x - eta), .Machine$double.eps)
  expect_equal(inc, expected)
})

# ---- compute_tsm ----

test_that("compute_tsm starts at initial and has correct length", {
  tsm  <- BettingTSM(bets = FixedBets(c = 0.5), eta = 0.5, initial = 2)
  x    <- runif(20)
  path <- compute_tsm(tsm, x)
  expect_equal(length(path), 20L)
  inc1 <- compute_increments(tsm, x)[1]
  expect_equal(path[1], 2 * inc1)
})

test_that("compute_tsm grows under the alternative", {
  set.seed(42)
  tsm  <- BettingTSM(bets = AGRAPABets(), eta = 0.3)
  x    <- runif(200, min = 0.6, max = 1)   # mean 0.8 >> eta = 0.3
  path <- compute_tsm(tsm, x)
  expect_gt(path[200], 10)
})

test_that("compute_tsm stays near 1 under the null on average", {
  set.seed(42)
  tsm <- BettingTSM(bets = AGRAPABets(), eta = 0.5)
  final_vals <- replicate(200, {
    x <- runif(100, 0, 1)   # mean = 0.5 = eta
    compute_tsm(tsm, x)[100]
  })
  # by Ville's inequality: P(max >= 1/alpha) <= alpha
  # so E[final] should not be dramatically inflated
  expect_lt(mean(final_vals > 20), 0.10)
})

# ---- integration with detector ----

test_that("BettingTSM integrates with ShiryaevRobertsDetector", {
  set.seed(42)
  tsm <- BettingTSM(bets = AGRAPABets(), eta = 0.5)
  x   <- c(runif(60, 0, 0.5), runif(100, 0.7, 1))
  inc <- compute_increments(tsm, x)
  det <- ShiryaevRobertsDetector(alpha = 0.05, criterion = "ARL")
  out <- run_detector(det, inc)
  expect_true(is.finite(out$stopping_time))
  expect_equal(length(out$statistic), length(x))
})

test_that("BettingTSM log mode gives same stopping time as standard mode", {
  set.seed(7)
  tsm <- BettingTSM(bets = FixedBets(c = 0.5), eta = 0.4)
  x   <- runif(100, 0.6, 1)
  inc     <- compute_increments(tsm, x, log = FALSE)
  log_inc <- compute_increments(tsm, x, log = TRUE)
  det <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")
  out_std <- run_detector(det, inc,     log = FALSE)
  out_log <- run_detector(det, log_inc, log = TRUE)
  expect_equal(out_std$stopping_time, out_log$stopping_time)
})
