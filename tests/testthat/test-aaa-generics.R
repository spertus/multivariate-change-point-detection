test_that("all package generics are registered", {
  expected <- c(
    "model_density",
    "likelihood_increment",
    "compute_increments",
    "compute_tsm",
    "update_detector",
    "run_detector",
    "combine_streams",
    "generate_stream"
  )

  for (g in expected) {
    expect_true(methods::isGeneric(g), info = sprintf("Generic '%s' should exist", g))
  }
})
