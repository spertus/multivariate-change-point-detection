test_that("default multivariate Gaussian DGP returns N-by-K matrix", {
  x <- default_multivariate_gaussian_dgp(
    N = 60,
    K = 3,
    nu = 30,
    pre_params = list(mu = c(0, 0, 0), Sigma = diag(3)),
    post_params = list(mu = c(1, 1, 1), Sigma = diag(3))
  )
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(60, 3))
})

test_that("default multivariate Gaussian DGP accepts `mean` alias for `mu`", {
  x <- default_multivariate_gaussian_dgp(
    N = 20,
    K = 3,
    nu = 10,
    pre_params = list(mean = c(0, 0, 0), Sigma = diag(3)),
    post_params = list(mean = c(1, 0, 0), Sigma = diag(3))
  )
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(20, 3))
})

test_that("expand_dgp_grid returns one DGP per parameter row", {
  template <- DGP(
    generator = default_gaussian_dgp,
    pre_params = list(mean = 0, sd = 1),
    post_params = list(mean = 1, sd = 1),
    nu = 50,
    name = "template"
  )

  out <- expand_dgp_grid(template, param_grid = list(mean = c(-1, 0, 1), sd = c(1, 2)))
  expect_length(out, 6)
  expect_true(all(vapply(out, function(obj) is(obj, "DGP"), logical(1))))
})

test_that("run_simulation returns expected columns", {
  dgp <- DGP(
    default_gaussian_dgp,
    pre_params = list(mean = 0, sd = 1),
    post_params = list(mean = 1, sd = 1),
    nu = 80,
    name = "gauss"
  )

  model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- TSM(model)
  det <- ShiryaevRobertsDetector(alpha = 0.1)

  out <- run_simulation(detector = det, tsm = tsm, dgp = dgp, n_rep = 20, N = 120, seed = 123)
  expect_true(all(c("detector", "dgp", "false_alarm_prob", "ARL", "ADD") %in% names(out)))
  expect_equal(nrow(out), 1)
})
