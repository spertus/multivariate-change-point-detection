test_that("increment sequence has requested length", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- SimpleVsSimpleTSM(m)
  x <- rnorm(25)
  inc <- compute_increments(tsm, x)
  expect_equal(length(inc), 25)
  expect_true(all(inc > 0))

  z <- increments_to_tsm(inc)
  expect_equal(length(z), 25)
  expect_true(all(z > 0))
})

test_that("SR detector returns finite stop under strong post-change data", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 2, sd_post = 1)
  tsm <- SimpleVsSimpleTSM(m)
  x <- c(rnorm(40, 0, 1), rnorm(40, 3, 1))
  inc <- compute_increments(tsm, x)
  d <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")
  out <- run_detector(d, inc)
  expect_true(is.finite(out$stopping_time))
})

test_that("CUSUM detector statistic stays nonnegative", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 0.3, sd_post = 1)
  tsm <- SimpleVsSimpleTSM(m)
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
