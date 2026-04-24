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
  m <- AR1Model(phi_pre = 0.1, sigma_pre = 1, mean_pre = 0, phi_post = 0.8, sigma_post = 1, mean_post = 0, x0 = 0)
  d_no_hist <- model_density(m, x = 1, regime = "post", history = NULL)
  d_hist    <- model_density(m, x = 1, regime = "post", history = c(10))
  expect_false(isTRUE(all.equal(d_no_hist, d_hist)))
})

test_that("AR1 model respects specified pre/post long-run means", {
  m <- AR1Model(phi_pre = 0.7, sigma_pre = 1, mean_pre = 2, phi_post = 0.7, sigma_post = 1, mean_post = -1, x0 = 2)
  d_pre  <- model_density(m, x = 2,  regime = "pre",  history = c(2))
  d_post <- model_density(m, x = 2,  regime = "post", history = c(2))
  expect_gt(d_pre, d_post)
})

test_that("AR1 long-run mean: density peaks at unconditional mean when process is at stationarity", {
  # When X_{t-1} = m, the conditional mean is mu + phi * m = m*(1-phi) + phi*m = m,
  # so the density should peak at x = m.
  m_val <- 3
  phi   <- 0.6
  m <- AR1Model(phi_pre = phi, sigma_pre = 1, mean_pre = m_val,
                phi_post = phi, sigma_post = 1, mean_post = m_val, x0 = m_val)
  d_at_mean  <- model_density(m, x = m_val,     regime = "pre", history = c(m_val))
  d_off_mean <- model_density(m, x = m_val + 1, regime = "pre", history = c(m_val))
  expect_gt(d_at_mean, d_off_mean)
})

# ---- ARpModel tests ----------------------------------------------------------

test_that("ARpModel constructor accepts valid stationary AR(2)", {
  m <- ARpModel(phi_pre  = c(0.5, 0.2), sigma_pre  = 1, mean_pre  = 0,
                phi_post = c(0.3, 0.1), sigma_post = 1, mean_post = 2)
  expect_s4_class(m, "ARpModel")
  expect_equal(m@mean_pre,  0)
  expect_equal(m@mean_post, 2)
})

test_that("ARpModel rejects non-stationary pre-change model", {
  # phi_1 = 1.0 gives unit root — not stationary
  expect_error(
    ARpModel(phi_pre  = c(1.0), sigma_pre  = 1, mean_pre  = 0,
             phi_post = c(0.5), sigma_post = 1, mean_post = 0),
    "Pre-change AR model is not stationary"
  )
})

test_that("ARpModel rejects non-stationary post-change model", {
  expect_error(
    ARpModel(phi_pre  = c(0.5), sigma_pre  = 1, mean_pre  = 0,
             phi_post = c(0.8, 0.6), sigma_post = 1, mean_post = 0),
    "Post-change AR model is not stationary"
  )
})

test_that("ARpModel model_density returns higher density at conditional mean", {
  phi <- c(0.5, 0.2)
  m <- ARpModel(phi_pre = phi, sigma_pre = 1, mean_pre = 0,
                phi_post = phi, sigma_post = 1, mean_post = 0)
  # With history = c(1, 2): lags = [2, 1]; intercept = 0; cond mean = 0.5*2 + 0.2*1 = 1.2
  cond_mean  <- 0.5 * 2 + 0.2 * 1
  d_at_cm    <- model_density(m, x = cond_mean,     regime = "pre", history = c(1, 2))
  d_off_cm   <- model_density(m, x = cond_mean + 2, regime = "pre", history = c(1, 2))
  expect_gt(d_at_cm, d_off_cm)
})

test_that("ARpModel model_density uses x0 when history is shorter than order", {
  phi <- c(0.4, 0.3)
  m <- ARpModel(phi_pre = phi, sigma_pre = 1, mean_pre = 0,
                phi_post = phi, sigma_post = 1, mean_post = 0, x0 = 5)
  # history length 1 < p=2: lags = [hist[1], x0] = [1, 5]; cond mean = 0.4*1 + 0.3*5 = 1.9
  cond_mean <- 0.4 * 1 + 0.3 * 5
  d_cm <- model_density(m, x = cond_mean, regime = "pre", history = c(1))
  d_off <- model_density(m, x = cond_mean + 3, regime = "pre", history = c(1))
  expect_gt(d_cm, d_off)
})

test_that("ARpModel long-run mean is recovered: when process is at m, cond mean = m", {
  phi   <- c(0.3, 0.2)
  m_val <- 4
  m <- ARpModel(phi_pre = phi, sigma_pre = 0.5, mean_pre = m_val,
                phi_post = phi, sigma_post = 0.5, mean_post = m_val)
  # Stationary: X_{t-1} = X_{t-2} = m  =>  cond mean = intercept + phi1*m + phi2*m
  # = m*(1-phi1-phi2) + (phi1+phi2)*m = m
  intercept  <- m_val * (1 - sum(phi))
  cond_mean  <- intercept + sum(phi) * m_val
  expect_equal(cond_mean, m_val, tolerance = 1e-12)
  d_cm  <- model_density(m, x = m_val,     regime = "pre", history = rep(m_val, 3))
  d_off <- model_density(m, x = m_val + 1, regime = "pre", history = rep(m_val, 3))
  expect_gt(d_cm, d_off)
})

test_that("ARpModel: AR(1) special case matches AR1Model output", {
  phi_v  <- 0.6; sigma_v <- 0.8; m_val <- 2
  m_arp  <- ARpModel(phi_pre  = phi_v, sigma_pre  = sigma_v, mean_pre  = m_val,
                     phi_post = phi_v, sigma_post = sigma_v, mean_post = m_val)
  m_ar1  <- AR1Model(phi_pre  = phi_v, sigma_pre  = sigma_v, mean_pre  = m_val,
                     phi_post = phi_v, sigma_post = sigma_v, mean_post = m_val, x0 = 0)
  hist <- c(1, 3)
  d_arp <- model_density(m_arp, x = 2.5, regime = "pre", history = hist)
  d_ar1 <- model_density(m_ar1, x = 2.5, regime = "pre", history = hist)
  expect_equal(d_arp, d_ar1, tolerance = 1e-12)
})

test_that("ARpModel: compute_increments runs end-to-end via generic fallback", {
  m <- ARpModel(phi_pre  = c(0.5, 0.1), sigma_pre  = 1, mean_pre  = 0,
                phi_post = c(0.5, 0.1), sigma_post = 0.5, mean_post = 3)
  tsm <- TSM(m)
  x   <- c(rnorm(20), rnorm(20, mean = 3, sd = 0.5))
  inc <- compute_increments(tsm, x, log = TRUE)
  expect_equal(length(inc), 40L)
  expect_true(all(is.finite(inc)))
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

test_that("Gaussian constructor supports fixed prior weights in mixture mode", {
  m <- GaussianModel(
    mean_pre = 0,
    sd_pre = 1,
    mean_post = c(-1, 1),
    sd_post = 1,
    method = "mixture",
    prior_weights = c(1, 2, 3, 4, 5),
    grid_size = 5
  )

  expect_s4_class(m, "GaussianModel")
  expect_equal(sum(m@prior_weights), 15)
})

test_that("Gaussian composite post model (predictable) uses full lagged history", {
  m <- GaussianModel(
    mean_pre = 0,
    sd_pre = 1,
    mean_post = c(-1, 2),
    sd_post = 1,
    method = "predictable"
  )

  h <- c(2, 2, -1)
  inc <- likelihood_increment(m, x = 2, history = h)
  ref <- dnorm(2, mean = mean(h), sd = 1) / dnorm(2, mean = 0, sd = 1)
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
    method = "predictable"
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
    method = "predictable"
  )

  h <- rbind(c(0.2, 0.1), c(0.3, 0.1), c(0.25, 0.2))
  inc <- likelihood_increment(m, x = c(1, 1), history = h)
  expect_true(is.finite(inc))
  expect_gt(inc, 0)
})
