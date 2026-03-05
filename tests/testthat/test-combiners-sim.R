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
