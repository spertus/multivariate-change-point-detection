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

test_that("product beats universal portfolio beats average", {
  s1 <- runif(20, 0.9, 1.3)
  s2 <- runif(20, 0.3, 1.9)
  s3 <- runif(20, 0.9, 1.9)

  p <- ProductCombiner()
  a <- AverageCombiner()
  up <- UniversalPortfolioCombiner(resolution = 4)
  z_p <- combine_streams(p, cbind(s1, s2, s3))
  z_a <- combine_streams(a, cbind(s1, s2, s3))
  z_up <- combine_streams(up, cbind(s1, s2, s3))
  expect_true(all(z_p > 0))
  expect_true(all(z_up > 0))
  expect_true(all(z_a > 0))
  expect_true(z_p[length(z_p)] > z_up[length(z_up)])
  expect_true(z_up[length(z_up)] > z_a[length(z_a)])
})
