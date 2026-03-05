make_null_setup <- function() {
  model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- SimpleVsSimpleTSM(model)
  dgp_null <- DGP(
    generator = default_gaussian_dgp,
    pre_params = list(mean = 0, sd = 1),
    post_params = list(mean = 1, sd = 1),
    nu = Inf,
    name = "iid-gaussian-null"
  )
  list(model = model, tsm = tsm, dgp_null = dgp_null)
}

test_that("S-R ARL detector meets nominal ARL level under no change", {
  setup <- make_null_setup()
  alpha <- 0.1

  det <- ShiryaevRobertsDetector(alpha = alpha, criterion = "ARL")
  mc <- run_simulation(
    detector = det,
    tsm = setup$tsm,
    dgp = setup$dgp_null,
    n_rep = 300,
    N = 800,
    seed = 1203
  )

  # Finite-horizon ARL estimator is conservative (downward-biased).
  expect_gte(mc$ARL, 1 / alpha)
})

test_that("S-R PFA detector controls false alarm probability at level alpha", {
  setup <- make_null_setup()
  alpha <- 0.1

  det <- ShiryaevRobertsDetector(alpha = alpha, criterion = "PFA")
  mc <- run_simulation(
    detector = det,
    tsm = setup$tsm,
    dgp = setup$dgp_null,
    n_rep = 1500,
    N = 600,
    seed = 1204
  )

  # Monte Carlo tolerance for finite replication.
  expect_lte(mc$false_alarm_prob, alpha + 0.02)
})

test_that("geometric spending terms are fixed, non-renormalized geometric masses", {
  p <- 0.01
  w5 <- multichangepoints:::.geometric_spending(5, p = p)
  expect_equal(w5[1], p)
  expect_equal(w5[2], p * (1 - p))
  expect_equal(w5[5], p * (1 - p)^4)
  expect_lt(sum(w5), 1)
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
