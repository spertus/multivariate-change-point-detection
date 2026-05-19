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

# ---- LocalizedDetector tests -----------------------------------------------

test_that("LocalizedDetector constructor stores K, alpha, criterion, threshold, invest", {
  ld <- LocalizedDetector(K = 4, alpha = 0.08, criterion = "ARL")
  expect_equal(ld@K, 4L)
  expect_equal(ld@alpha, 0.08)
  expect_equal(ld@criterion, "ARL")
  expect_equal(ld@threshold, 1 / 0.08)
  # default uniform Bonferroni allowance: 1-row matrix with every entry = 1/K
  expect_equal(ld@invest, matrix(rep(1/4, 4), nrow = 1L, ncol = 4L))
})

test_that("LocalizedDetector constructor respects custom threshold", {
  ld <- LocalizedDetector(K = 2, alpha = 0.1, criterion = "ARL", threshold = 50)
  expect_equal(ld@threshold, 50)
})

test_that("LocalizedDetector ARL: custom allowance matrix is stored and validated", {
  # Non-uniform but valid: stream 1 gets 70%, stream 2 gets 30%
  alloc <- matrix(c(0.7, 0.3), nrow = 1L)
  ld <- LocalizedDetector(K = 2, alpha = 0.05, criterion = "ARL", allowance = alloc)
  expect_equal(ld@invest, alloc)
  # Row-sums-not-1 should error
  bad <- matrix(c(0.6, 0.3), nrow = 1L)
  expect_error(LocalizedDetector(K = 2, alpha = 0.05, criterion = "ARL", allowance = bad),
               "sum to 1")
})

test_that("LocalizedDetector PFA: joint spending matrix is validated", {
  N <- 50
  # Valid N x K joint schedule
  base   <- matrix(runif(N * 2), nrow = N, ncol = 2)
  base   <- base / sum(base)
  ld <- LocalizedDetector(K = 2, alpha = 0.05, criterion = "PFA", spending = base)
  expect_equal(ld@invest, base)
  # Sum != 1 should error
  bad <- base * 2
  expect_error(LocalizedDetector(K = 2, alpha = 0.05, criterion = "PFA", spending = bad),
               "sum to 1")
})

test_that("LocalizedDetector uniform ARL gives same stopping times as manual Bonferroni", {
  # Manual Bonferroni: R_t = Lambda * (R_{t-1} + 1), threshold K/alpha
  # New: R_t = Lambda * (R_{t-1} + 1/K), threshold 1/alpha  — mathematically identical
  K <- 2; n <- 200; alpha <- 0.05; Lambda <- 1.3
  inc <- matrix(rep(Lambda, n * K), nrow = n, ncol = K)
  ld  <- LocalizedDetector(K = K, alpha = alpha, criterion = "ARL")
  out <- run_detector(ld, evidence = inc)
  # Manual computation with allowance 1/K and threshold 1/alpha
  R <- 0
  manual_stop <- Inf
  for (t in seq_len(n)) {
    R <- Lambda * (R + 1/K)
    if (R >= 1/alpha) { manual_stop <- t; break }
  }
  expect_equal(out$stream_results$stream_1$stopping_time, manual_stop)
})

test_that("LocalizedDetector run_detector returns K stream results", {
  set.seed(1)
  K <- 3
  n <- 30
  ld  <- LocalizedDetector(K = K, alpha = 0.05, criterion = "ARL")
  ev  <- matrix(rep(1.3, n * K), nrow = n, ncol = K)
  out <- run_detector(ld, evidence = ev)
  expect_equal(length(out$stream_results), K)
  expect_equal(names(out$stream_results), paste0("stream_", seq_len(K)))
  for (k in seq_len(K)) {
    expect_equal(length(out$stream_results[[k]]$statistic), n)
  }
})

test_that("LocalizedDetector identifies correct first-alarming stream", {
  # Stream 1: strong evidence (should alarm early); Stream 2: flat evidence (no alarm)
  n   <- 100
  ev  <- matrix(c(rep(1.4, n), rep(1.0, n)), nrow = n, ncol = 2)
  ld  <- LocalizedDetector(K = 2, alpha = 0.1, criterion = "ARL")
  out <- run_detector(ld, evidence = ev)
  expect_true(out$alarm)
  expect_equal(out$first_alarm_stream, 1L)
  expect_equal(out$stopping_time, out$stream_results$stream_1$stopping_time)
})

test_that("LocalizedDetector global stopping time is min of marginal stopping times", {
  # Construct two detectors with known stopping times
  n   <- 50
  ev  <- matrix(rep(1.5, n * 3), nrow = n, ncol = 3)
  ld  <- LocalizedDetector(K = 3, alpha = 0.05, criterion = "ARL")
  out <- run_detector(ld, evidence = ev)
  marginal_stops <- vapply(out$stream_results, function(r) r$stopping_time, numeric(1L))
  expect_equal(out$stopping_time, min(marginal_stops))
})

test_that("LocalizedDetector no alarm when evidence is below 1", {
  # Evidence < 1 means R_t has a stable fixed point below threshold; no alarm ever fires.
  n   <- 200
  ev  <- matrix(rep(0.5, n * 3), nrow = n, ncol = 3)
  ld  <- LocalizedDetector(K = 3, alpha = 0.05, criterion = "ARL")
  out <- run_detector(ld, evidence = ev)
  expect_false(out$alarm)
  expect_equal(out$stopping_time, Inf)
  expect_true(is.na(out$first_alarm_stream))
})

test_that("LocalizedDetector log mode matches standard mode for stopping times", {
  K   <- 2
  n   <- 80
  # Stream 1: evidence 1.4 throughout; Stream 2: evidence 0.9 throughout.
  # Both streams have the same evidence pattern, so log mode should give same stopping times.
  inc <- matrix(c(rep(1.4, n), rep(0.9, n)), nrow = n, ncol = K)
  ld  <- LocalizedDetector(K = K, alpha = 0.1, criterion = "ARL")
  out_std <- run_detector(ld, evidence = inc,      log = FALSE)
  out_log <- run_detector(ld, evidence = log(inc), log = TRUE)
  expect_equal(out_std$stopping_time, out_log$stopping_time)
  for (k in seq_len(K)) {
    expect_equal(
      out_std$stream_results[[k]]$stopping_time,
      out_log$stream_results[[k]]$stopping_time
    )
  }
})

test_that("LocalizedDetector rejects non-matrix evidence", {
  ld <- LocalizedDetector(K = 2, alpha = 0.05, criterion = "ARL")
  expect_error(run_detector(ld, evidence = c(1.1, 1.2, 1.3)), "numeric N-by-K matrix")
})

test_that("LocalizedDetector rejects wrong number of columns", {
  ld <- LocalizedDetector(K = 3, alpha = 0.05, criterion = "ARL")
  ev <- matrix(rep(1.1, 50 * 2), nrow = 50)
  expect_error(run_detector(ld, evidence = ev), "3 column")
})
