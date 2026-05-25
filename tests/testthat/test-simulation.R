test_that("DGP from MultivariateGaussianModel generates N-by-K matrix", {
  model <- MultivariateGaussianModel(
    mu_pre = c(0, 0, 0), Sigma_pre = diag(3),
    mu_post = c(1, 1, 1), Sigma_post = diag(3)
  )
  x <- generate_stream(DGP(model, nu = 30), N = 60)
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(60, 3))
})

test_that("DGP from MultivariateGaussianModel: sparse post-change", {
  model <- MultivariateGaussianModel(
    mu_pre = c(0, 0, 0), Sigma_pre = diag(3),
    mu_post = c(1, 0, 0), Sigma_post = diag(3)
  )
  x <- generate_stream(DGP(model, nu = 10), N = 20)
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(20, 3))
})

test_that("expand_dgp_grid returns one DGP per parameter row", {
  gen <- function(N, K, nu, pre_params, post_params) {
    stats::rnorm(N, pre_params$mean, pre_params$sd)
  }
  template <- DGP(
    generator  = gen,
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
  model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  dgp   <- DGP(model, nu = 80, name = "gauss")
  tsm   <- TSM(model)
  det   <- ShiryaevRobertsDetector(alpha = 0.1)

  out <- run_simulation(detector = det, tsm = tsm, dgp = dgp, n_rep = 20, N = 120, seed = 123)
  expect_true(all(c("detector", "dgp", "false_alarm_prob", "ARL", "ADD") %in% names(out)))
  expect_equal(nrow(out), 1)
})
