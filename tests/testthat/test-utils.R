test_that(".assert_numeric_vector accepts valid vectors and rejects invalid inputs", {
  expect_no_error(multichangepoints:::.assert_numeric_vector(c(1, 2, 3), "x"))
  expect_error(multichangepoints:::.assert_numeric_vector(numeric(0), "x"))
  expect_no_error(multichangepoints:::.assert_numeric_vector(c(1, NA_real_), "x"))   # NAs allowed (offline streams)
  expect_error(multichangepoints:::.assert_numeric_vector(c(1, Inf), "x"))            # non-finite non-NA still rejected
  expect_error(multichangepoints:::.assert_numeric_vector("abc", "x"))
})

test_that(".logsumexp is numerically stable", {
  x <- c(1000, 999, 998)
  out <- multichangepoints:::.logsumexp(x)
  ref <- 1000 + log(sum(exp(x - 1000)))
  expect_equal(out, ref, tolerance = 1e-12)
})

test_that("geometric spending terms are fixed, non-renormalized geometric masses", {
  p <- 0.01
  w5 <- multichangepoints:::.geometric_spending(5, p = p)
  expect_equal(w5[1], p)
  expect_equal(w5[2], p * (1 - p))
  expect_equal(w5[5], p * (1 - p)^4)
  expect_lt(sum(w5), 1)
})

test_that("composition grid rows sum to m", {
  g <- multichangepoints:::.composition_grid(K = 3, m = 4)
  expect_true(is.matrix(g))
  expect_true(all(rowSums(g) == 4))
})
