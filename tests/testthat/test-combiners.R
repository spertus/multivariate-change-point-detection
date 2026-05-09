test_that("product and average combiners produce length-N output", {
  i1 <- runif(30, 0.8, 1.2)
  i2 <- runif(30, 0.8, 1.2)
  streams <- cbind(i1, i2)

  p <- ProductCombiner()
  a <- AverageCombiner()

  zp <- combine_streams(p, streams)
  za <- combine_streams(a, streams, weights = c(0.3, 0.7))

  expect_equal(length(zp), 30)
  expect_equal(length(za), 30)
})

test_that("product/average combiners support log-increment mode", {
  i1 <- runif(20, 0.8, 1.2)
  i2 <- runif(20, 0.8, 1.2)
  streams <- cbind(i1, i2)
  log_streams <- log(streams)

  p <- ProductCombiner()
  a <- AverageCombiner()

  zp <- combine_streams(p, streams, log = FALSE)
  zap <- combine_streams(p, log_streams, log = TRUE)
  za <- combine_streams(a, streams, weights = c(0.4, 0.6), log = FALSE)
  zal <- combine_streams(a, log_streams, weights = c(0.4, 0.6), log = TRUE)

  expect_equal(zap, log(zp), tolerance = 1e-10)
  expect_equal(zal, log(za), tolerance = 1e-10)
})

test_that("sparse UP (default) is positive", {
  s1 <- runif(20, 0.9, 1.3)
  s2 <- runif(20, 0.9, 1.3)
  up <- UniversalPortfolioCombiner()          # sparse, K+1 portfolios
  z  <- combine_streams(up, cbind(s1, s2))
  expect_true(all(z > 0))
})

test_that("dense UP is positive", {
  s1 <- runif(20, 0.9, 1.3)
  s2 <- runif(20, 0.9, 1.3)
  up <- UniversalPortfolioCombiner(mode = "dense", resolution = 4)
  z  <- combine_streams(up, cbind(s1, s2))
  expect_true(all(z > 0))
})

# ---- sparse UP structural properties ----------------------------------------

test_that("sparse UP: first step equals simple average for K=2", {
  # With K=2 and equal initial wealth across 3 portfolios {1,0},{0,1},{0.5,0.5},
  # the first combined increment = ((s1 + s2 + (s1+s2)/2)) / 3 = (s1+s2)/2.
  s1 <- 1.5; s2 <- 0.8
  up <- UniversalPortfolioCombiner()
  z  <- combine_streams(up, matrix(c(s1, s2), nrow = 1))
  expect_equal(z, (s1 + s2) / 2, tolerance = 1e-12)
})

test_that("sparse UP: first step equals simple average for K=3", {
  # With K=3, portfolios {1,0,0},{0,1,0},{0,0,1},{1/3,1/3,1/3} (equal wealth):
  # combined = (s1 + s2 + s3 + (s1+s2+s3)/3) / 4 = (s1+s2+s3)/3 = simple avg.
  s <- c(2.0, 0.5, 1.2)
  up <- UniversalPortfolioCombiner()
  z  <- combine_streams(up, matrix(s, nrow = 1))
  expect_equal(z, mean(s), tolerance = 1e-12)
})

test_that("sparse UP: adapts toward better stream over many steps", {
  set.seed(5)
  N  <- 200
  # stream 2 consistently doubles; stream 1 stays near 1
  s1 <- runif(N, 0.9, 1.1)
  s2 <- runif(N, 1.8, 2.2)
  up <- UniversalPortfolioCombiner()
  z  <- combine_streams(up, cbind(s1, s2))
  # At t=1, z equals the average (equi-weighted mix).
  # After adaptation toward s2, late z/s2 should exceed early z/s2.
  early_ratio <- mean(z[1:5]    / s2[1:5])     # z starts at avg < s2
  late_ratio  <- mean(z[150:N]  / s2[150:N])   # z converges toward s2
  expect_gt(late_ratio, early_ratio)
})

test_that("sparse UP log/non-log modes agree", {
  set.seed(3)
  streams     <- matrix(runif(40, 0.7, 1.4), 20, 2)
  log_streams <- log(streams)
  up  <- UniversalPortfolioCombiner()
  z   <- combine_streams(up, streams,     log = FALSE)
  zl  <- combine_streams(up, log_streams, log = TRUE)
  expect_equal(zl, log(z), tolerance = 1e-10)
})

# ---- dense UP ---------------------------------------------------------------

test_that("dense UP: errors when grid exceeds 50 000 portfolios", {
  # K=50, resolution=4 => choose(53, 49) = choose(53,4) = 292 825 > 50 000
  up <- UniversalPortfolioCombiner(mode = "dense", resolution = 4)
  streams <- matrix(runif(10 * 50), 10, 50)
  expect_error(combine_streams(up, streams), "50 000")
})

test_that("dense UP log/non-log modes agree", {
  set.seed(8)
  streams     <- matrix(runif(30 * 3, 0.8, 1.3), 30, 3)
  log_streams <- log(streams)
  up  <- UniversalPortfolioCombiner(mode = "dense", resolution = 4)
  z   <- combine_streams(up, streams,     log = FALSE)
  zl  <- combine_streams(up, log_streams, log = TRUE)
  expect_equal(zl, log(z), tolerance = 1e-10)
})

test_that("dense UP tracks the best portfolio at least as well as simple avg", {
  set.seed(11)
  N  <- 50
  # stream 1 dominates; best single-stream portfolio = {1, 0}
  s1 <- runif(N, 1.4, 1.8)
  s2 <- runif(N, 0.7, 1.0)
  streams <- cbind(s1, s2)
  up_dense <- UniversalPortfolioCombiner(mode = "dense", resolution = 4)
  avg      <- AverageCombiner()
  zd <- combine_streams(up_dense, streams)
  za <- combine_streams(avg,      streams)
  # UP wealth must exceed average wealth (UP contains the average portfolio)
  expect_gte(prod(zd), prod(za) * 0.99)
})

# ---- offline stream (leading NA) tests ----------------------------------------

# Helper: build a 2-stream log-increment matrix where stream 2 comes online at t=onset+1
make_offline_streams <- function(n, onset, seed = 1) {
  set.seed(seed)
  s1 <- rnorm(n, mean = 0.1, sd = 0.3)        # always online; slightly positive
  s2 <- c(rep(NA_real_, onset), rnorm(n - onset, mean = 0.1, sd = 0.3))
  cbind(s1, s2)
}

test_that("compute_increments returns NA for NA inputs and non-NA elsewhere", {
  model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm   <- TSM(model)
  x     <- c(NA, NA, 0.5, 1.0, -0.3)
  inc   <- compute_increments(tsm, x, log = TRUE)

  expect_equal(length(inc), 5L)
  expect_true(all(is.na(inc[1:2])))
  expect_true(all(!is.na(inc[3:5])))
})

test_that("compute_increments history ignores NA: same result as fresh stream with same x values", {
  model   <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm     <- TSM(model)
  x_late  <- c(NA, NA, -0.3, 0.8)
  x_fresh <- c(-0.3, 0.8)           # same values, no leading NAs

  inc_late  <- compute_increments(tsm, x_late,  log = TRUE)
  inc_fresh <- compute_increments(tsm, x_fresh, log = TRUE)

  # NAs are skipped; history is empty when -0.3 is first observed in both cases
  expect_equal(inc_late[3], inc_fresh[1], tolerance = 1e-12)
  expect_equal(inc_late[4], inc_fresh[2], tolerance = 1e-12)
})

test_that("ProductCombiner: NA-free output; offline stream contributes nothing (log)", {
  streams <- make_offline_streams(20, onset = 5)
  prd <- ProductCombiner()
  out <- combine_streams(prd, streams, log = TRUE)

  expect_equal(length(out), 20L)
  expect_false(anyNA(out))
  # During offline period t=1..5, product = stream 1 only
  expect_equal(out[1:5], streams[1:5, 1], tolerance = 1e-12)
  # After both online, product = sum of both log-increments
  expect_equal(out[6:20], rowSums(streams[6:20, ]), tolerance = 1e-12)
})

test_that("ProductCombiner: NA-free output in natural scale", {
  streams <- exp(make_offline_streams(15, onset = 4))
  streams[1:4, 2] <- NA
  prd <- ProductCombiner()
  out <- combine_streams(prd, streams, log = FALSE)

  expect_false(anyNA(out))
  expect_equal(out[1:4], streams[1:4, 1], tolerance = 1e-12)
})

test_that("AverageCombiner: offline stream not diluted — matches stream 1 alone during offline period", {
  streams <- make_offline_streams(20, onset = 8)
  avg <- AverageCombiner()
  out_log <- combine_streams(avg, streams, log = TRUE)
  out_nat <- combine_streams(avg, exp(streams[, , drop = FALSE]), log = FALSE)
  streams_nat <- exp(streams)
  streams_nat[1:8, 2] <- NA

  expect_false(anyNA(out_log))
  expect_false(anyNA(out_nat))

  # During offline period: average over 1 stream = that stream exactly
  expect_equal(out_log[1:8], streams[1:8, 1], tolerance = 1e-12)
  # Natural scale: same logic — weight is entirely on stream 1
  expect_equal(out_nat[1:8], streams_nat[1:8, 1], tolerance = 1e-12)
})

test_that("AverageCombiner: log and non-log modes agree with offline streams", {
  streams_log <- make_offline_streams(20, onset = 6)
  streams_nat <- exp(streams_log)
  streams_nat[1:6, 2] <- NA

  avg <- AverageCombiner()
  out_log <- combine_streams(avg, streams_log, log = TRUE)
  out_nat <- combine_streams(avg, streams_nat, log = FALSE)

  expect_equal(out_log, log(out_nat), tolerance = 1e-10)
})

test_that("UniversalPortfolioCombiner: NA-free output with offline stream", {
  streams <- make_offline_streams(20, onset = 5)
  up <- UniversalPortfolioCombiner(mode = "dense", resolution = 4)
  out <- combine_streams(up, streams, log = TRUE)

  expect_equal(length(out), 20L)
  expect_false(anyNA(out))
})

test_that("UniversalPortfolioCombiner: log and non-log modes agree with offline streams", {
  streams_log <- make_offline_streams(15, onset = 5)
  streams_nat <- exp(streams_log)
  streams_nat[1:5, 2] <- NA

  up <- UniversalPortfolioCombiner(mode = "dense", resolution = 4)
  out_log <- combine_streams(up, streams_log, log = TRUE)
  out_nat <- combine_streams(up, streams_nat, log = FALSE)

  expect_equal(out_log, log(out_nat), tolerance = 1e-10)
})

test_that("combined TSM with offline stream is non-negative (Ville validity check)", {
  # A valid TSM must stay non-negative; under the null the product should
  # stay near 1 in expectation. We check that it stays finite and positive.
  set.seed(42)
  n <- 100; onset <- 30
  model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 0.5, sd_post = 1)
  tsm   <- TSM(model)

  x1 <- rnorm(n)                               # null data, always online
  x2 <- c(rep(NA, onset), rnorm(n - onset))   # null data, comes online late

  inc1 <- compute_increments(tsm, x1, log = TRUE)
  inc2 <- compute_increments(tsm, x2, log = TRUE)
  streams <- cbind(inc1, inc2)

  for (Combiner in list(ProductCombiner(), AverageCombiner(),
                        UniversalPortfolioCombiner(mode = "dense", resolution = 4))) {
    combined <- combine_streams(Combiner, streams, log = TRUE)
    expect_false(anyNA(combined))
    expect_true(all(is.finite(combined)))
  }
})

test_that("ShiryaevRoberts detector runs end-to-end with offline stream combined increments", {
  set.seed(7)
  n <- 80; onset <- 20
  model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm   <- TSM(model)

  x1 <- c(rnorm(40), rnorm(40, mean = 1))     # change at t=41
  x2 <- c(rep(NA, onset), rnorm(20), rnorm(n - onset - 20, mean = 1))

  inc1 <- compute_increments(tsm, x1, log = TRUE)
  inc2 <- compute_increments(tsm, x2, log = TRUE)
  combined <- combine_streams(ProductCombiner(), cbind(inc1, inc2), log = TRUE)

  sr  <- ShiryaevRobertsDetector(alpha = 0.001, criterion = "ARL")
  out <- run_detector(sr, combined, log = TRUE)

  expect_true(is.numeric(out$statistic))
  expect_equal(length(out$statistic), n)
  expect_true(is.finite(out$stopping_time) || out$stopping_time == Inf)
})

test_that("product beats universal portfolio beats average", {
  s1 <- runif(20, 0.9, 1.3)
  s2 <- runif(20, 0.3, 1.9)
  s3 <- runif(20, 0.9, 3.0)

  p <- ProductCombiner()
  a <- AverageCombiner()
  up <- UniversalPortfolioCombiner(mode = "dense", resolution = 4)
  z_p <- combine_streams(p, cbind(s1, s2, s3))
  z_a <- combine_streams(a, cbind(s1, s2, s3))
  z_up <- combine_streams(up, cbind(s1, s2, s3))
  expect_true(all(z_p > 0))
  expect_true(all(z_up > 0))
  expect_true(all(z_a > 0))
  expect_true(z_p[length(z_p)] > z_up[length(z_up)])
  expect_true(z_up[length(z_up)] > z_a[length(z_a)])
})
