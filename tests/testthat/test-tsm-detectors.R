test_that("increment sequence has requested length", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- TSM(m)
  x <- rnorm(25)
  inc <- compute_increments(tsm, x)
  log_inc <- compute_increments(tsm, x, log = TRUE)
  expect_equal(length(inc), 25)
  expect_true(all(inc > 0))
  expect_equal(log_inc, log(inc), tolerance = 1e-10)

  z <- increments_to_tsm(inc)
  expect_equal(length(z), 25)
  expect_true(all(z > 0))
})

test_that("SR detector returns finite stop under strong post-change data", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 2, sd_post = 1)
  tsm <- TSM(m)
  x <- c(rnorm(40, 0, 1), rnorm(40, 3, 1))
  inc <- compute_increments(tsm, x)
  d <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")
  out <- run_detector(d, inc)
  expect_true(is.finite(out$stopping_time))
})

test_that("CUSUM detector statistic stays nonnegative", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 0.3, sd_post = 1)
  tsm <- TSM(m)
  inc <- compute_increments(tsm, rnorm(50))
  d <- CUSUMDetector(alpha = 0.1)
  out <- run_detector(d, inc)
  expect_true(all(out$statistic >= 0))
})

test_that("detector statistic is computed for full horizon after first alarm", {
  d <- ShiryaevRobertsDetector(alpha = 0.2, criterion = "ARL")
  out <- run_detector(d, evidence = rep(1.3, 20))
  expect_true(is.finite(out$stopping_time))
  expect_equal(length(out$statistic), 20)
  expect_true(all(out$statistic > 0))
})

test_that("joint multivariate Gaussian pipeline works without combiners", {
  model <- MultivariateGaussianModel(
    mu_pre = c(0, 0),
    Sigma_pre = diag(2),
    mu_post = c(1, 1),
    Sigma_post = diag(2)
  )
  tsm <- TSM(model)

  x <- cbind(
    c(rnorm(40, 0, 1), rnorm(40, 2, 1)),
    c(rnorm(40, 0, 1), rnorm(40, 2, 1))
  )

  inc <- compute_increments(tsm, x)
  det <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")
  out <- run_detector(det, inc)

  expect_equal(length(inc), nrow(x))
  expect_true(is.finite(out$stopping_time))
})

test_that("SR detector log mode matches standard mode for stopping time", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- TSM(m)
  x <- c(rnorm(40, 0, 1), rnorm(40, 2, 1))

  inc <- compute_increments(tsm, x, log = FALSE)
  log_inc <- compute_increments(tsm, x, log = TRUE)
  d <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")

  out_std <- run_detector(d, inc, log = FALSE)
  out_log <- run_detector(d, log_inc, log = TRUE)

  expect_equal(out_log$stopping_time, out_std$stopping_time)
  expect_equal(exp(out_log$statistic), out_std$statistic, tolerance = 1e-8)
})

test_that("composite Gaussian post model works through increment pipeline", {
  model <- GaussianModel(
    mean_pre = 0,
    sd_pre = 1,
    mean_post = c(-1, 2),
    sd_post = 1,
    method = "predictable",
    update_window = 5
  )
  tsm <- SimpleVsSimpleTSM(model)
  x <- c(rnorm(50, 0, 1), rnorm(50, 1.5, 1))
  inc <- compute_increments(tsm, x)
  det <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")
  out <- run_detector(det, inc)

  expect_equal(length(inc), length(x))
  expect_true(all(is.finite(inc)))
  expect_true(is.finite(out$stopping_time))
})

test_that("multivariate composite Gaussian model works through joint increment pipeline", {
  K <- 2
  model <- MultivariateGaussianModel(
    mu_pre = rep(0, K),
    Sigma_pre = diag(K),
    mu_post = cbind(rep(-1, K), rep(2, K)),
    Sigma_post = NULL,
    method = "predictable",
    update_window = 5
  )
  tsm <- TSM(model)
  x <- cbind(c(rnorm(40, 0, 1), rnorm(40, 1.5, 1)),
             c(rnorm(40, 0, 1), rnorm(40, 0.5, 1)))
  inc <- compute_increments(tsm, x)
  det <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")
  out <- run_detector(det, inc)

  expect_equal(length(inc), nrow(x))
  expect_true(all(is.finite(inc)))
  expect_true(is.finite(out$stopping_time))
})
