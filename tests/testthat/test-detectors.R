test_that("SR detector returns finite stop under strong evidence", {
  d <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")
  out <- run_detector(d, evidence = rep(1.3, 60))
  expect_true(is.finite(out$stopping_time))
})

test_that("CUSUM detector statistic stays nonnegative", {
  d <- CUSUMDetector(alpha = 0.1)
  out <- run_detector(d, evidence = exp(rnorm(50, mean = 0.1, sd = 0.2)))
  expect_true(all(out$statistic >= 0))
})

test_that("detector statistic is computed for full horizon after first alarm", {
  d <- ShiryaevRobertsDetector(alpha = 0.2, criterion = "ARL")
  out <- run_detector(d, evidence = rep(1.3, 20))
  expect_true(is.finite(out$stopping_time))
  expect_equal(length(out$statistic), 20)
  expect_true(all(out$statistic > 0))
})

test_that("SR detector log mode matches standard mode for stopping time", {
  inc <- c(rep(0.95, 30), rep(1.15, 50))
  d <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")

  out_std <- run_detector(d, inc, log = FALSE)
  out_log <- run_detector(d, log(inc), log = TRUE)

  expect_equal(out_log$stopping_time, out_std$stopping_time)
  expect_equal(exp(out_log$statistic), out_std$statistic, tolerance = 1e-8)
})

test_that("S-R PFA recursion injects precomputed pi_t terms with fixed threshold", {
  det <- ShiryaevRobertsDetector(
    alpha = 0.1,
    criterion = "PFA",
    spending = c(0.2, 0.1, 0.7),
    threshold = 10
  )

  out <- run_detector(det, evidence = c(1, 1, 1))
  expect_equal(out$statistic[1], 0.2)
  expect_equal(out$statistic[2], 0.3)
  expect_equal(out$statistic[3], 1.0)
  expect_false(out$alarm)
})
