test_that("Gaussian model likelihood increments are positive", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  inc <- likelihood_increment(m, x = c(0, 1, 2))
  expect_equal(length(inc), 3)
  expect_true(all(inc > 0))
})

test_that("Bernoulli model increment favors post-change when x=1 and p_post>p_pre", {
  m <- BernoulliModel(p_pre = 0.2, p_post = 0.8)
  inc1 <- likelihood_increment(m, x = 1)
  inc0 <- likelihood_increment(m, x = 0)
  expect_gt(inc1, 1)
  expect_lt(inc0, 1)
})

test_that("AR1 model uses history", {
  m <- AR1Model(phi_pre = 0.1, sigma_pre = 1, mu_0 = 0, phi_post = 0.8, sigma_post = 1, mu_1 = 0, x0 = 0)
  d_no_hist <- model_density(m, x = 1, regime = "post", history = NULL)
  d_hist <- model_density(m, x = 1, regime = "post", history = c(10))
  expect_false(isTRUE(all.equal(d_no_hist, d_hist)))
})

test_that("AR1 model respects specified pre/post process means", {
  m <- AR1Model(phi_pre = 0.7, sigma_pre = 1, mu_0 = 2, phi_post = 0.7, sigma_post = 1, mu_1 = -1, x0 = 2)
  d_pre <- model_density(m, x = 2, regime = "pre", history = c(2))
  d_post <- model_density(m, x = 2, regime = "post", history = c(2))
  expect_gt(d_pre, d_post)
})

test_that("AR1 conditional mean uses intercept parameterization", {
  m <- AR1Model(phi_pre = 0.5, sigma_pre = 1, mu_0 = 2, phi_post = 0.5, sigma_post = 1, mu_1 = 2, x0 = 0)
  d_at_mean <- model_density(m, x = 4, regime = "pre", history = c(4))
  d_off_mean <- model_density(m, x = 3, regime = "pre", history = c(4))
  expect_gt(d_at_mean, d_off_mean)
})

test_that("Multivariate Gaussian model returns scalar density and valid increment", {
  m <- MultivariateGaussianModel(
    mu_pre = c(0, 0),
    Sigma_pre = diag(2),
    mu_post = c(1, 1),
    Sigma_post = diag(2)
  )

  d <- model_density(m, x = c(0, 0), regime = "pre")
  inc <- likelihood_increment(m, x = c(0.2, -0.1))

  expect_true(is.numeric(d) && length(d) == 1)
  expect_true(is.numeric(inc) && length(inc) == 1)
  expect_gt(inc, 0)
})

test_that("Gaussian constructor supports mix_weight_adapt argument", {
  m <- GaussianModel(
    mean_pre = 0,
    sd_pre = 1,
    mean_post = c(-1, 1),
    sd_post = 1,
    method = "mixture",
    mix_weight_adapt = "none"
  )

  expect_s4_class(m, "GaussianModel")
  expect_identical(m@mix_weight_adapt, "none")
})

test_that(".predictable_mean_estimate returns mean when inside bounds", {
  x <- rnorm(100, 2, 1)
  lag_mean <- .predictable_mean_estimate(
    history = x,
    update_window = 5,
    lower = 0,
    upper = 3,
    init_mean = 0
  )
  lag_bound_mean <- .predictable_mean_estimate(
    history = x,
    update_window = 5,
    lower = 0,
    upper = 1,
    init_mean = 0
  )
  expect_equal(lag_mean, mean(x))
  expect_equal(lag_bound_mean, 1)
})

test_that("Gaussian composite post model (predictable) respects update window", {
  m <- GaussianModel(
    mean_pre = 0,
    sd_pre = 1,
    mean_post = c(-1, 2),
    sd_post = 1,
    method = "predictable",
    update_window = 2
  )

  h <- c(2, 2, -1)
  inc <- likelihood_increment(m, x = 2, history = h)
  ref <- dnorm(2, mean = 2, sd = 1) / dnorm(2, mean = 0, sd = 1)
  expect_equal(inc, ref, tolerance = 1e-10)
})

test_that("Gaussian composite post model (mixture) returns finite positive increment", {
  m <- GaussianModel(
    mean_pre = 0,
    sd_pre = 1,
    mean_post = c(-1, 1),
    sd_post = 1,
    method = "mixture",
    grid_size = 5
  )

  inc <- likelihood_increment(m, x = 0.3, history = c(-0.1, 0.2))
  log_inc <- likelihood_increment(m, x = 0.3, history = c(-0.1, 0.2), log = TRUE)
  expect_gt(inc, 0)
  expect_true(is.finite(log_inc))
  expect_equal(log(inc), log_inc, tolerance = 1e-10)
})

test_that("Gaussian composite post model can estimate post-change sd from lagged data", {
  m <- GaussianModel(
    mean_pre = 0,
    sd_pre = 1,
    mean_post = c(-1, 2),
    sd_post = numeric(0),
    method = "predictable",
    update_window = 3
  )

  inc <- likelihood_increment(m, x = 2, history = c(2, 2, 2, -1))
  expect_true(is.finite(inc))
  expect_gt(inc, 0)
})

test_that("Multivariate Gaussian composite post model supports box constraints and unknown Sigma", {
  K <- 2
  m <- MultivariateGaussianModel(
    mu_pre = c(0, 0),
    Sigma_pre = diag(K),
    mu_post = cbind(c(-1, -1), c(2, 2)),
    Sigma_post = NULL,
    method = "predictable",
    update_window = 2
  )

  h <- rbind(c(0.2, 0.1), c(0.3, 0.1), c(0.25, 0.2))
  inc <- likelihood_increment(m, x = c(1, 1), history = h)
  expect_true(is.finite(inc))
  expect_gt(inc, 0)
})
