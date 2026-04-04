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

# ---- multiple_alarms tests -------------------------------------------------

test_that("multiple_alarms=TRUE: no alarm when evidence is clearly below 1 (S-R decreasing)", {
  # Λ_t = 0.5 < 1: S-R recursion R_t = (1 + R_{t-1}) * 0.5 has fixed point at 1,
  # well below threshold = 1/0.05 = 20, so no alarm should ever fire.
  d <- ShiryaevRobertsDetector(alpha = 0.05, criterion = "ARL", multiple_alarms = TRUE)
  out <- run_detector(d, evidence = rep(0.5, 500), log = FALSE)
  expect_equal(length(out$alarm_times), 0L)
  expect_false(out$alarm)
})

test_that("multiple_alarms=TRUE: no alarm when evidence is clearly below 1 (log scale)", {
  # log(0.5) = -0.693; same fixed-point argument as above.
  d <- ShiryaevRobertsDetector(alpha = 0.05, criterion = "ARL", multiple_alarms = TRUE)
  out <- run_detector(d, evidence = rep(log(0.5), 500), log = TRUE)
  expect_equal(length(out$alarm_times), 0L)
  expect_false(out$alarm)
})

test_that("multiple_alarms=TRUE: multiple alarms raised when TSM grows continuously", {
  # Strong evidence: increments of 1.5 throughout. With alpha=0.05 (threshold=20)
  # S-R will alarm repeatedly after resetting.
  d   <- ShiryaevRobertsDetector(alpha = 0.05, criterion = "ARL", multiple_alarms = TRUE)
  out <- run_detector(d, evidence = rep(1.5, 500), log = FALSE)
  expect_true(length(out$alarm_times) > 1L)
})

test_that("multiple_alarms=TRUE: statistic resets to initial value after each alarm (log scale)", {
  # Use a very low threshold so alarms occur frequently; check that the statistic
  # drops back to -Inf (log of 1) immediately after each alarm.
  d   <- ShiryaevRobertsDetector(alpha = 0.5, criterion = "ARL", multiple_alarms = TRUE)  # threshold = 2
  inc <- rep(log(2), 100)   # each step doubles R_t: threshold hit at t=1, reset, again at t=2, ...
  out <- run_detector(d, evidence = inc, log = TRUE)

  # Should alarm at almost every step
  expect_true(length(out$alarm_times) > 10L)

  # The step after each alarm (except the last) the statistic should equal
  # the single-step value from a fresh start: log(1 + 1) + log(2) = log(2)
  for (i in seq_along(out$alarm_times)) {
    next_t <- out$alarm_times[i] + 1L
    if (next_t > length(out$statistic)) next
    expect_equal(out$statistic[next_t], log(2), tolerance = 1e-10)
  }
})

test_that("multiple_alarms=TRUE: backward-compatible stopping_time is first alarm", {
  d   <- ShiryaevRobertsDetector(alpha = 0.05, criterion = "ARL", multiple_alarms = TRUE)
  inc <- c(rep(1, 50), rep(1.5, 200))
  out <- run_detector(d, evidence = inc, log = FALSE)
  expect_equal(out$stopping_time, out$alarm_times[1L])
})

test_that("multiple_alarms=FALSE (default): alarm_times has at most one entry", {
  d   <- ShiryaevRobertsDetector(alpha = 0.05, criterion = "ARL")
  out <- run_detector(d, evidence = rep(1.5, 100), log = FALSE)
  expect_lte(length(out$alarm_times), 1L)
})

test_that("multiple_alarms=TRUE: PFA spending recycles after reset", {
  # Spending: pi = (0.5, 0.5) recycled; threshold = 2 (alpha=0.5).
  # With evidence=1 (log-inc 0): R_t = (pi_t + R_{t-1}) * 1.
  # Cycle 1: t=1: R=(0.5+0)*1=0.5; t=2: R=(0.5+0.5)*1=1.0 — no alarm (threshold=2).
  # Increase evidence so alarm triggers during first cycle; then verify spending restarts.
  spend <- c(0.6, 0.4)
  d <- ShiryaevRobertsDetector(alpha = 0.5, criterion = "PFA",
                                spending = spend, threshold = 2,
                                multiple_alarms = TRUE)
  # evidence=3: R_1 = (0.6 + 0) * 3 = 1.8; R_2 = (0.4 + 1.8) * 3 = 6.6 >= 2 => alarm at t=2
  # After reset, t=3 uses spend[1]=0.6 again: R_3 = (0.6 + 0) * 3 = 1.8
  out <- run_detector(d, evidence = rep(3, 6), log = FALSE)
  expect_true(length(out$alarm_times) >= 2L)
  # statistic at first step of each new cycle should equal (spend[1]) * evidence = 0.6 * 3 = 1.8
  for (i in seq_along(out$alarm_times)) {
    next_t <- out$alarm_times[i] + 1L
    if (next_t > length(out$statistic)) next
    expect_equal(out$statistic[next_t], 0.6 * 3, tolerance = 1e-10)
  }
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
