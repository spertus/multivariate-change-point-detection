test_that("increment sequence has requested length", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- TSM(m)
  x <- rnorm(25)

  inc <- compute_increments(tsm, x)
  log_inc <- compute_increments(tsm, x, log = TRUE)

  expect_equal(length(inc), 25)
  expect_true(all(inc > 0))
  expect_equal(log_inc, log(inc), tolerance = 1e-10)
})

test_that("increments_to_tsm handles standard and log modes", {
  inc <- c(1.1, 0.9, 1.05, 1.2)
  log_inc <- log(inc)

  z <- increments_to_tsm(inc, initial = 1, log = FALSE)
  log_z <- increments_to_tsm(log_inc, initial = 1, log = TRUE)

  expect_equal(log(z), log_z, tolerance = 1e-10)
})

test_that("compute_tsm returns a path with one value per observation", {
  m <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
  tsm <- TSM(m)
  x <- rnorm(30)

  z <- compute_tsm(tsm, x)
  log_z <- compute_tsm(tsm, x, log = TRUE)

  expect_equal(length(z), length(x))
  expect_equal(log(z), log_z, tolerance = 1e-10)
})
