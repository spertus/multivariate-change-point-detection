skip_if_not_integration()

make_null_setup <- function() {
  model    <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm      <- TSM(model)
  dgp_null <- DGP(model, nu = Inf, name = "iid-gaussian-null")
  list(model = model, tsm = tsm, dgp_null = dgp_null)
}

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

  expect_lte(mc$false_alarm_prob, alpha + 0.02)
})
