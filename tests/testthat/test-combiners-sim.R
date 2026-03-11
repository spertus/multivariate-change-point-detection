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

test_that("universal portfolio combiner is positive", {
  s1 <- runif(20, 0.9, 1.3)
  s2 <- runif(20, 0.9, 1.3)
  up <- UniversalPortfolioCombiner(resolution = 4)
  z <- combine_streams(up, cbind(s1, s2))
  expect_true(all(z > 0))
})

test_that("simulation returns expected columns", {
  dgp <- DGP(default_gaussian_dgp,
             pre_params = list(mean = 0, sd = 1),
             post_params = list(mean = 1, sd = 1),
             nu = 80,
             name = "gauss")

  model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- SimpleVsSimpleTSM(model)
  det <- ShiryaevRobertsDetector(alpha = 0.1)

  out <- run_simulation(detector = det, tsm = tsm, dgp = dgp, n_rep = 20, N = 120, seed = 123)
  expect_true(all(c("detector", "dgp", "false_alarm_prob", "ARL", "ADD") %in% names(out)))
  expect_equal(nrow(out), 1)
})

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

test_that("2-D multivariate Gaussian DGP works in joint-model simulation path", {
  dgp_mv <- DGP(
    generator = default_multivariate_gaussian_dgp,
    pre_params = list(mu = c(0, 0), Sigma = diag(2)),
    post_params = list(mu = c(1, 1), Sigma = diag(2)),
    nu = 50,
    name = "mv-gauss"
  )

  model_mv <- MultivariateGaussianModel(
    mu_pre = c(0, 0),
    Sigma_pre = diag(2),
    mu_post = c(1, 1),
    Sigma_post = diag(2)
  )
  tsm_mv <- SimpleVsSimpleTSM(model_mv)
  det <- ShiryaevRobertsDetector(alpha = 0.1, criterion = "ARL")

  out <- run_simulation(detector = det, tsm = tsm_mv, dgp = dgp_mv, n_rep = 15, N = 120, K = 2, seed = 777)
  expect_true(all(c("detector", "dgp", "false_alarm_prob", "ARL", "ADD") %in% names(out)))
  expect_equal(nrow(out), 1)
})

test_that("3-D multivariate Gaussian DGP generates a stream", {
  dgp_mv <- DGP(
    generator = default_multivariate_gaussian_dgp,
    pre_params = list(mu = c(0, 0, 0), Sigma = diag(3)),
    post_params = list(mu = c(1, 1, 1), Sigma = diag(3)),
    nu = 50,
    name = "mv-gauss"
  )
  
  out <- generate_stream(dgp_mv, N = 10, K = 3)
  expect_equal(ncol(out), 3)
})
