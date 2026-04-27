## test-varp.R  — unit tests for VARpModel

# ---- constructor -----------------------------------------------------------

test_that("VARpModel constructs with default shared Phi and Sigma", {
  K <- 2
  Phi1 <- matrix(c(0.4, 0.1, 0.1, 0.3), K, K)
  Sig  <- diag(K)
  m <- VARpModel(Phi_pre = Phi1, Sigma_pre = Sig,
                 mean_pre = c(0, 0), mean_post = c(1, 1))
  expect_s4_class(m, "VARpModel")
  expect_equal(m@Phi_post,   list(Phi1))
  expect_equal(m@Sigma_post, Sig)
  expect_equal(m@mean_pre,   c(0, 0))
  expect_equal(m@mean_post,  c(1, 1))
})

test_that("VARpModel accepts a bare matrix as order-1 Phi", {
  K   <- 3
  Phi <- 0.3 * diag(K)
  m   <- VARpModel(Phi_pre = Phi, Sigma_pre = diag(K),
                   mean_pre = rep(0, K), mean_post = rep(2, K))
  expect_length(m@Phi_pre, 1L)
  expect_equal(m@Phi_pre[[1]], Phi)
})

test_that("VARpModel accepts order-2 Phi as a list", {
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  Phi2 <- 0.1 * diag(K)
  m    <- VARpModel(Phi_pre = list(Phi1, Phi2), Sigma_pre = diag(K),
                    mean_pre = c(0, 0), mean_post = c(1, 0))
  expect_length(m@Phi_pre, 2L)
})

test_that("VARpModel rejects non-stationary pre-change model", {
  K    <- 2
  # Eigenvalue of companion = eigenvalue of Phi1 = 1 (unit root)
  Phi1 <- diag(c(1.0, 0.3), K)
  expect_error(
    VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
              mean_pre = c(0, 0), mean_post = c(1, 1)),
    "Pre-change VAR model is not stationary"
  )
})

test_that("VARpModel rejects non-stationary post-change model", {
  K      <- 2
  Phi_ok <- 0.3 * diag(K)
  Phi_bad <- diag(c(0.5, 1.1), K)
  expect_error(
    VARpModel(Phi_pre = Phi_ok, Sigma_pre = diag(K),
              Phi_post = Phi_bad, Sigma_post = diag(K),
              mean_pre = c(0, 0), mean_post = c(1, 1)),
    "Post-change VAR model is not stationary"
  )
})

test_that("VARpModel rejects mismatched K in Phi and mean", {
  K   <- 2
  Phi <- 0.3 * diag(K)
  expect_error(
    VARpModel(Phi_pre = Phi, Sigma_pre = diag(K),
              mean_pre = c(0, 0, 0), mean_post = c(1, 1, 1))   # K=3 ≠ K=2
  )
})

test_that("VARpModel x0 defaults to zero vector of length K", {
  K <- 3
  m <- VARpModel(Phi_pre = 0.2 * diag(K), Sigma_pre = diag(K),
                 mean_pre = rep(0, K), mean_post = rep(1, K))
  expect_equal(m@x0, rep(0, K))
})

# ---- model_density ---------------------------------------------------------

test_that("model_density returns positive scalar", {
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                    mean_pre = c(0, 0), mean_post = c(2, 2))
  d <- model_density(m, x = c(0, 0), regime = "pre", history = NULL)
  expect_length(d, 1L)
  expect_gt(d, 0)
})

test_that("model_density peaks at conditional mean", {
  K    <- 2
  Phi1 <- matrix(c(0.4, 0.0, 0.0, 0.3), K, K)
  Sig  <- diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = Sig,
                    mean_pre = c(0, 0), mean_post = c(3, 3))
  hist <- matrix(c(1, 2), nrow = 1)  # one row = one lag
  # conditional mean = nu + Phi1 %*% [1, 2]
  nu   <- as.numeric((diag(K) - Phi1) %*% c(0, 0))   # = 0
  cm   <- as.numeric(nu + Phi1 %*% c(1, 2))           # = [0.4, 0.6]
  d_cm  <- model_density(m, x = cm,       regime = "pre", history = hist)
  d_off <- model_density(m, x = cm + 2,   regime = "pre", history = hist)
  expect_gt(d_cm, d_off)
})

test_that("model_density uses x0 when history is empty", {
  K    <- 2
  x0   <- c(5, 5)
  Phi1 <- 0.4 * diag(K)
  Sig  <- diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = Sig,
                    mean_pre = c(0, 0), mean_post = c(1, 1), x0 = x0)
  nu     <- as.numeric((diag(K) - Phi1) %*% c(0, 0))   # 0
  cm_x0  <- as.numeric(nu + Phi1 %*% x0)               # 0.4 * [5,5] = [2, 2]
  cm_zero <- as.numeric(nu + Phi1 %*% c(0, 0))          # [0, 0]
  d_x0   <- model_density(m, x = cm_x0,  regime = "pre", history = NULL)
  d_zero <- model_density(m, x = cm_zero, regime = "pre", history = NULL)
  expect_gt(d_x0, d_zero)
})

test_that("model_density pre vs post differ when means differ", {
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                    mean_pre = c(0, 0), mean_post = c(4, 4))
  # At x = [4,4], post-change density (centred near 4) > pre-change (centred near 0)
  d_pre  <- model_density(m, x = c(4, 4), regime = "pre",  history = NULL)
  d_post <- model_density(m, x = c(4, 4), regime = "post", history = NULL)
  expect_gt(d_post, d_pre)
})

test_that("model_density: long-run mean recovered at stationarity", {
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  m_val <- c(2, 3)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                    mean_pre = m_val, mean_post = m_val)
  # When all lags = m, cond mean = nu + Phi*m = (I-Phi)m + Phi*m = m
  hist  <- matrix(rep(m_val, 5), nrow = 5, byrow = TRUE)
  d_cm  <- model_density(m, x = m_val,     regime = "pre", history = hist)
  d_off <- model_density(m, x = m_val + 1, regime = "pre", history = hist)
  expect_gt(d_cm, d_off)
})

# ---- VAR(2) order -----------------------------------------------------------

test_that("VAR(2): model_density uses both lags correctly", {
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  Phi2 <- 0.1 * diag(K)
  m    <- VARpModel(Phi_pre = list(Phi1, Phi2), Sigma_pre = diag(K),
                    mean_pre = c(0, 0), mean_post = c(1, 1))
  hist <- matrix(c(1, 1,   # lag-2 row (older)
                   2, 2),  # lag-1 row (more recent)
                 nrow = 2, byrow = TRUE)
  nu   <- rep(0, K)
  cm   <- nu + as.numeric(Phi1 %*% c(2, 2)) + as.numeric(Phi2 %*% c(1, 1))
  d_cm  <- model_density(m, x = cm,     regime = "pre", history = hist)
  d_off <- model_density(m, x = cm + 2, regime = "pre", history = hist)
  expect_gt(d_cm, d_off)
})

# ---- likelihood increments and TSM pipeline --------------------------------

test_that("likelihood_increment: post > pre when observation is near post mean", {
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                    mean_pre = c(0, 0), mean_post = c(3, 3))
  # x clearly in post-change territory
  inc <- likelihood_increment(m, x = c(3, 3), history = NULL, log = TRUE)
  expect_gt(inc, 0)
})

test_that("likelihood_increment: near pre-change mean gives negative log-increment", {
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                    mean_pre = c(0, 0), mean_post = c(5, 5))
  inc <- likelihood_increment(m, x = c(0, 0), history = NULL, log = TRUE)
  expect_lt(inc, 0)
})

test_that("compute_increments returns length-N finite vector for VARpModel", {
  set.seed(1)
  K    <- 3
  Phi1 <- 0.3 * diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                    mean_pre = rep(0, K), mean_post = rep(1, K))
  x <- matrix(rnorm(50 * K), nrow = 50, ncol = K)
  inc <- compute_increments(TSM(m), x, log = TRUE)
  expect_equal(length(inc), 50L)
  expect_true(all(is.finite(inc)))
})

test_that("compute_increments: TSM grows post-change under correct alternative", {
  set.seed(7)
  K    <- 2
  Phi1 <- 0.3 * diag(K)
  Sig  <- diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = Sig,
                    mean_pre = c(0, 0), mean_post = c(2, 2))
  # Pre-change data: mean 0
  x_pre  <- matrix(rnorm(30 * K), nrow = 30) %*% chol(Sig)
  # Post-change data: mean 2
  x_post <- matrix(rnorm(70 * K, mean = 2), nrow = 70) %*% chol(Sig)
  x <- rbind(x_pre, x_post)
  inc  <- compute_increments(TSM(m), x, log = TRUE)
  path <- increments_to_tsm(inc, log = TRUE)
  # TSM should grow after the change on average
  expect_gt(mean(path[71:100]), mean(path[1:30]))
})

test_that("run_detector raises alarm on VAR stream with clear mean shift", {
  set.seed(99)
  K    <- 2
  Phi1 <- 0.2 * diag(K)
  m    <- VARpModel(Phi_pre = Phi1, Sigma_pre = diag(K),
                    mean_pre = c(0, 0), mean_post = c(3, 3))
  x_pre  <- matrix(rnorm(30 * K), 30, K)
  x_post <- matrix(rnorm(120 * K, mean = 3), 120, K)
  x <- rbind(x_pre, x_post)
  inc <- compute_increments(TSM(m), x, log = TRUE)
  sr  <- ShiryaevRobertsDetector(alpha = 0.001, criterion = "ARL")
  out <- run_detector(sr, inc, log = TRUE)
  expect_true(is.finite(out$stopping_time))
  expect_gt(out$stopping_time, 30)   # should not alarm before the change
})
