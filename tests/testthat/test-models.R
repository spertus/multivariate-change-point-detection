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
  # Under intercept parameterization, E[X_t | X_{t-1}=4] = 2 + 0.5 * 4 = 4
  d_at_mean <- model_density(m, x = 4, regime = "pre", history = c(4))
  d_off_mean <- model_density(m, x = 3, regime = "pre", history = c(4))
  expect_gt(d_at_mean, d_off_mean)
})
