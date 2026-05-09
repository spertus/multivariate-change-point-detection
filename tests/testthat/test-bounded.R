## test-bounded.R — unit tests for BoundedModel

# ---- constructor -----------------------------------------------------------

test_that("BoundedModel constructs with valid eta", {
  m <- BoundedModel(eta = 0.3)
  expect_s4_class(m, "BoundedModel")
  expect_equal(m@eta,    0.3)
  expect_equal(m@c,      0.75)
  expect_equal(m@sd_min, 0.01)
  expect_equal(m@eps,    1e-5)
})

test_that("BoundedModel accepts custom tuning parameters", {
  m <- BoundedModel(eta = 0.5, c = 0.5, sd_min = 0.05, eps = 1e-4)
  expect_equal(m@c,      0.5)
  expect_equal(m@sd_min, 0.05)
  expect_equal(m@eps,    1e-4)
})

test_that("BoundedModel rejects eta outside (0, 1)", {
  expect_error(BoundedModel(eta = 0),   "`eta` must be a scalar in \\(0, 1\\)")
  expect_error(BoundedModel(eta = 1),   "`eta` must be a scalar in \\(0, 1\\)")
  expect_error(BoundedModel(eta = -0.1),"`eta` must be a scalar in \\(0, 1\\)")
  expect_error(BoundedModel(eta = 1.1), "`eta` must be a scalar in \\(0, 1\\)")
})

test_that("BoundedModel rejects invalid c", {
  expect_error(BoundedModel(eta = 0.5, c = 0),   "`c` must be a scalar in \\(0, 1\\]")
  expect_error(BoundedModel(eta = 0.5, c = 1.1), "`c` must be a scalar in \\(0, 1\\]")
})

test_that("BoundedModel rejects non-positive sd_min", {
  expect_error(BoundedModel(eta = 0.5, sd_min = 0), "`sd_min` must be a positive scalar")
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
  m    <- BoundedModel(eta = 0.3)
  hist <- c(0.7, 0.8, 0.9)
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
  m    <- BoundedModel(eta = eta, c = c)
  hist <- rep(0.99, 100)   # extreme history drives lambda up
  inc  <- likelihood_increment(m, x = 1, history = hist)
  cap  <- c / eta
  max_inc <- 1 + cap * (1 - eta)
  expect_lte(inc, max_inc + 1e-10)
})

test_that("single-observation history uses sd_min floor", {
  m    <- BoundedModel(eta = 0.3, sd_min = 0.01)
  inc1 <- likelihood_increment(m, x = 0.8, history = c(0.7))
  inc2 <- likelihood_increment(m, x = 0.8, history = c(0.7, 0.7))
  # both should be finite positive
  expect_true(is.finite(inc1) && inc1 > 0)
  expect_true(is.finite(inc2) && inc2 > 0)
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
  x_pre  <- runif(30, 0.1, 0.5)   # straddles eta
  x_post <- runif(70, 0.7, 1.0)   # clearly above eta
  x      <- c(x_pre, x_post)
  inc    <- compute_increments(TSM(m), x, log = TRUE)
  path   <- increments_to_tsm(inc, log = TRUE)
  expect_gt(mean(path[61:100]), mean(path[1:30]))
})
