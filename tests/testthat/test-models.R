# ---- GaussianVARModel constructor -----------------------------------------------

test_that("GaussianModel wrapper returns GaussianVARModel (K=1 IID)", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  expect_s4_class(m, "GaussianVARModel")
  expect_equal(length(m@Phi_pre),   0L)
  expect_equal(length(m@Phi_post),  0L)
  expect_equal(m@mean_pre,  0)
  expect_equal(m@mean_post, 1)
  expect_equal(m@Sigma_pre,  matrix(1))
  expect_equal(m@Sigma_post, matrix(1))
})

test_that("MultivariateGaussianModel wrapper returns GaussianVARModel (K>1 IID)", {
  m <- MultivariateGaussianModel(mu_pre = c(0, 0), Sigma_pre = diag(2),
                                  mu_post = c(1, 1), Sigma_post = diag(2))
  expect_s4_class(m, "GaussianVARModel")
  expect_equal(length(m@Phi_pre), 0L)
  expect_equal(m@mean_pre,  c(0, 0))
  expect_equal(m@mean_post, c(1, 1))
})

test_that("GaussianVARModel accepts numeric vector Phi_pre as univariate AR coefficients", {
  phi <- c(0.5, 0.2, 0.1)
  m   <- GaussianVARModel(Phi_pre = phi, Sigma_pre = matrix(1), mean_pre = 0, mean_post = 1)
  expect_equal(length(m@Phi_pre), 3L)
  expect_equal(m@Phi_pre[[1]], matrix(0.5))
  expect_equal(m@Phi_pre[[2]], matrix(0.2))
  expect_equal(m@Phi_pre[[3]], matrix(0.1))
})

test_that("GaussianVARModel accepts bare matrix Phi_pre as VAR(1)", {
  K   <- 2
  Phi <- 0.3 * diag(K)
  m   <- GaussianVARModel(Phi_pre = Phi, Sigma_pre = diag(K),
                     mean_pre = c(0, 0), mean_post = c(1, 1))
  expect_equal(length(m@Phi_pre), 1L)
  expect_equal(m@Phi_pre[[1]], Phi)
})

test_that("GaussianVARModel x0 defaults to zero vector", {
  m <- GaussianVARModel(Phi_pre = list(), Sigma_pre = matrix(1),
                   mean_pre = 0, mean_post = 1)
  expect_equal(m@x0, 0)
})

test_that("GaussianVARModel rejects non-stationary pre-change model", {
  expect_error(
    GaussianVARModel(Phi_pre = matrix(1.0), Sigma_pre = matrix(1),
                mean_pre = 0, mean_post = 1),
    "Pre-change VAR model is not stationary"
  )
})

test_that("GaussianVARModel rejects non-stationary post-change model", {
  expect_error(
    GaussianVARModel(Phi_pre = matrix(0.3), Sigma_pre = matrix(1),
                mean_pre = 0, Phi_post = matrix(1.1), Sigma_post = matrix(1),
                mean_post = 1),
    "Post-change VAR model is not stationary"
  )
})

test_that("GaussianVARModel rejects mismatched mean lengths", {
  expect_error(
    GaussianVARModel(Phi_pre = list(), Sigma_pre = matrix(1),
                mean_pre = 0, mean_post = c(1, 2)),
    "same length"
  )
})

# ---- model_density for GaussianVARModel ------------------------------------------

test_that("model_density returns positive scalar for IID GaussianVARModel (K=1)", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  d <- model_density(m, x = 0, regime = "pre")
  expect_length(d, 1L)
  expect_gt(d, 0)
})

test_that("model_density peaks at the long-run mean for IID model (K=1)", {
  m    <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 3, sd_post = 1)
  d_cm <- model_density(m, x = 0, regime = "pre")
  d_off <- model_density(m, x = 2, regime = "pre")
  expect_gt(d_cm, d_off)
})

test_that("model_density: post density exceeds pre density near post mean (K=2 IID)", {
  m     <- MultivariateGaussianModel(mu_pre = c(0, 0), mu_post = c(3, 3))
  d_pre  <- model_density(m, x = c(3, 3), regime = "pre")
  d_post <- model_density(m, x = c(3, 3), regime = "post")
  expect_gt(d_post, d_pre)
})

test_that("model_density for AR(1) GaussianVARModel peaks at conditional mean", {
  phi <- 0.6
  m   <- GaussianVARModel(Phi_pre = matrix(phi), Sigma_pre = matrix(1),
                     mean_pre = 0, mean_post = 2)
  # history = 1 → cond mean = (1-phi)*0 + phi*1 = 0.6
  cm   <- phi * 1
  d_cm  <- model_density(m, x = cm,     regime = "pre", history = 1)
  d_off <- model_density(m, x = cm + 2, regime = "pre", history = 1)
  expect_gt(d_cm, d_off)
})

test_that("model_density for AR(2) GaussianVARModel uses both lags", {
  phi1 <- 0.5; phi2 <- 0.2
  m    <- GaussianVARModel(Phi_pre = c(phi1, phi2), Sigma_pre = matrix(1),
                      mean_pre = 0, mean_post = 3)
  # history = c(1, 2) → lag-1 = 2, lag-2 = 1; cond mean = phi1*2 + phi2*1 = 1.2
  cm   <- phi1 * 2 + phi2 * 1
  hist <- c(1, 2)
  d_cm  <- model_density(m, x = cm,     regime = "pre", history = hist)
  d_off <- model_density(m, x = cm + 3, regime = "pre", history = hist)
  expect_gt(d_cm, d_off)
})

test_that("model_density for VAR(1) GaussianVARModel peaks at conditional mean (K=2)", {
  K    <- 2
  Phi1 <- matrix(c(0.4, 0.0, 0.0, 0.3), K, K)
  m    <- GaussianVARModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                      mean_pre = c(0, 0), mean_post = c(2, 2))
  hist <- matrix(c(1, 2), nrow = 1)   # one lag row
  cm   <- as.numeric(Phi1 %*% c(1, 2))
  d_cm  <- model_density(m, x = cm,     regime = "pre", history = hist)
  d_off <- model_density(m, x = cm + 2, regime = "pre", history = hist)
  expect_gt(d_cm, d_off)
})

test_that("model_density uses x0 when history is empty (K=1 AR(1))", {
  x0  <- 5
  phi <- 0.4
  m   <- GaussianVARModel(Phi_pre = matrix(phi), Sigma_pre = matrix(1),
                     mean_pre = 0, mean_post = 3, x0 = x0)
  cm_x0   <- phi * x0
  cm_zero <- phi * 0
  d_x0   <- model_density(m, x = cm_x0,   regime = "pre", history = NULL)
  d_zero <- model_density(m, x = cm_zero, regime = "pre", history = NULL)
  expect_gt(d_x0, d_zero)
})

# ---- BernoulliModel ---------------------------------------------------------

test_that("BernoulliModel increment favors post-change when x=1 and p_post > p_pre", {
  m    <- BernoulliModel(p_pre = 0.2, p_post = 0.8)
  inc1 <- likelihood_increment(m, x = 1)
  inc0 <- likelihood_increment(m, x = 0)
  expect_gt(inc1, 1)
  expect_lt(inc0, 1)
})

# ---- likelihood_increment for GaussianVARModel -----------------------------------

test_that("likelihood_increment returns positive scalar for IID K=1", {
  m   <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  inc <- sapply(c(0, 1, 2), function(xi) likelihood_increment(m, x = xi))
  expect_length(inc, 3L)
  expect_true(all(inc > 0))
})

test_that("likelihood_increment: log=TRUE and linear are consistent", {
  m       <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 2, sd_post = 1)
  inc     <- likelihood_increment(m, x = 1.5, log = FALSE)
  log_inc <- likelihood_increment(m, x = 1.5, log = TRUE)
  expect_equal(log(inc), log_inc, tolerance = 1e-12)
})

test_that("likelihood_increment for AR(1) GaussianVARModel uses history", {
  phi <- 0.7
  m   <- GaussianVARModel(Phi_pre = matrix(phi), Sigma_pre = matrix(1),
                     mean_pre = 0, mean_post = 4)
  inc_no_hist <- likelihood_increment(m, x = 1, history = NULL,   log = TRUE)
  inc_hist    <- likelihood_increment(m, x = 1, history = c(10),  log = TRUE)
  expect_false(isTRUE(all.equal(inc_no_hist, inc_hist)))
})

test_that("likelihood_increment: near post mean gives positive log increment (K=2)", {
  m   <- MultivariateGaussianModel(mu_pre = c(0, 0), mu_post = c(3, 3))
  inc <- likelihood_increment(m, x = c(3, 3), log = TRUE)
  expect_gt(inc, 0)
})

test_that("compute_increments returns length-N finite vector for IID GaussianVARModel", {
  m   <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  x   <- rnorm(30)
  inc <- compute_increments(TSM(m), x, log = TRUE)
  expect_equal(length(inc), 30L)
  expect_true(all(is.finite(inc)))
})

test_that("compute_increments via matrix path works for MultivariateGaussianModel", {
  set.seed(1)
  m   <- MultivariateGaussianModel(mu_pre = c(0, 0), mu_post = c(1, 1))
  x   <- matrix(rnorm(40 * 2), nrow = 40, ncol = 2)
  inc <- compute_increments(TSM(m), x, log = TRUE)
  expect_equal(length(inc), 40L)
  expect_true(all(is.finite(inc)))
})
