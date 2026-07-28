# ---- DependencyGraph construction / validation ------------------------------

test_that("DependencyGraph: blocks constructor builds correct block-diagonal adjacency", {
  g <- DependencyGraph(blocks = list(c(1, 2), c(3, 4, 5)), K = 5)
  expected <- matrix(0, 5, 5)
  expected[1, 2] <- expected[2, 1] <- 1
  expected[3, 4] <- expected[4, 3] <- 1
  expected[3, 5] <- expected[5, 3] <- 1
  expected[4, 5] <- expected[5, 4] <- 1
  expect_equal(g@A, expected)
  expect_equal(g@V, as.character(1:5))
})

test_that("DependencyGraph: adjacency-matrix constructor stores A directly", {
  A <- matrix(0, 3, 3)
  A[1, 2] <- A[2, 1] <- 1
  g <- DependencyGraph(A = A)
  expect_equal(g@A, A)
})

test_that("DependencyGraph: requires exactly one of A / blocks", {
  expect_error(DependencyGraph(), "Exactly one")
  A <- matrix(0, 2, 2)
  expect_error(DependencyGraph(A = A, blocks = list(1, 2)), "Exactly one")
})

test_that("DependencyGraph: rejects non-symmetric, non-zero-diagonal, non-binary matrices", {
  A_asym <- matrix(c(0, 1, 0, 0), 2, 2)
  expect_error(DependencyGraph(A = A_asym), "symmetric")

  A_diag <- matrix(c(1, 0, 0, 0), 2, 2)
  expect_error(DependencyGraph(A = A_diag), "zero diagonal")

  A_nonbin <- matrix(c(0, 0.5, 0.5, 0), 2, 2)
  expect_error(DependencyGraph(A = A_nonbin), "0 or 1")
})

test_that("dependency_blocks: labels connected components in visit order", {
  g <- DependencyGraph(blocks = list(c(1, 2), c(3, 4, 5)), K = 5)
  expect_equal(dependency_blocks(g), c(1L, 1L, 2L, 2L, 2L))
})

test_that("dependency_blocks: all-singleton graph gives K distinct blocks", {
  g <- DependencyGraph(A = matrix(0, 4, 4))
  expect_equal(dependency_blocks(g), 1:4)
})

test_that("dependency_blocks: fully-connected graph gives one block", {
  A <- matrix(1, 3, 3); diag(A) <- 0
  g <- DependencyGraph(A = A)
  expect_equal(dependency_blocks(g), c(1L, 1L, 1L))
})

# ---- ProximityGraph construction / validation --------------------------------

test_that("ProximityGraph: stores W and validates symmetric non-negative zero-diagonal", {
  W <- matrix(0, 3, 3)
  W[1, 2] <- W[2, 1] <- 2
  g <- ProximityGraph(W)
  expect_equal(g@W, W)

  W_asym <- matrix(c(0, 1, 2, 0), 2, 2)
  expect_error(ProximityGraph(W_asym), "symmetric")

  W_diag <- matrix(c(1, 0, 0, 0), 2, 2)
  expect_error(ProximityGraph(W_diag), "zero diagonal")

  W_neg <- matrix(0, 2, 2); W_neg[1, 2] <- W_neg[2, 1] <- -1
  expect_error(ProximityGraph(W_neg), "non-negative")
})

# ---- graph_distance_matrix ----------------------------------------------------

test_that("graph_distance_matrix: linear chain 1-2-3-4 with unit weights", {
  W <- matrix(0, 4, 4)
  for (i in 1:3) W[i, i + 1] <- W[i + 1, i] <- 1
  g <- ProximityGraph(W)
  D <- graph_distance_matrix(g)

  expect_equal(unname(D[1, ]), c(0, 1, 2, 3))
  expect_equal(unname(D[2, 4]), 2)
  expect_equal(unname(graph_distance(g, 1, 4)), 3)
})

test_that("graph_distance_matrix: hub-and-spoke — spoke-to-spoke goes through the hub", {
  K <- 4  # node 1 = hub, 2:4 = spokes, weight 2 per spoke edge
  W <- matrix(0, K, K)
  for (k in 2:K) W[1, k] <- W[k, 1] <- 2
  g <- ProximityGraph(W)
  D <- graph_distance_matrix(g)

  expect_equal(unname(D[1, 2]), 2)
  expect_equal(unname(D[2, 3]), 4)   # spoke -> hub -> spoke
  expect_equal(unname(D[3, 4]), 4)
})

test_that("graph_distance_matrix: unreachable nodes have distance Inf", {
  W <- matrix(0, 3, 3)
  W[1, 2] <- W[2, 1] <- 1   # node 3 isolated
  g <- ProximityGraph(W)
  D <- graph_distance_matrix(g)
  expect_equal(unname(D[1, 3]), Inf)
  expect_equal(unname(D[3, 3]), 0)
})

# ---- proximity kernels ---------------------------------------------------------

test_that("exponential_kernel: matches exp(-d/xi), zero at Inf", {
  k <- exponential_kernel(xi = 2)
  expect_equal(k(0), 1)
  expect_equal(k(2), exp(-1))
  expect_equal(k(Inf), 0)
  expect_equal(k(c(0, 2, Inf)), c(1, exp(-1), 0))
})

test_that("hard_cutoff_kernel: 1{d <= r}, zero at Inf", {
  k <- hard_cutoff_kernel(r = 2)
  expect_equal(k(0), 1)
  expect_equal(k(2), 1)
  expect_equal(k(2.0001), 0)
  expect_equal(k(Inf), 0)
})

# ---- proximity_allowance (Section 5.3) -----------------------------------------

test_that("proximity_allowance: zeta = 0 always gives uniform 1/K", {
  W <- matrix(0, 3, 3)
  W[1, 2] <- W[2, 1] <- W[2, 3] <- W[3, 2] <- 1
  g <- ProximityGraph(W)
  k <- exponential_kernel(xi = 1)

  expect_equal(proximity_allowance(g, active = integer(0), kernel = k, zeta = 0), rep(1/3, 3))
  expect_equal(proximity_allowance(g, active = c(1), kernel = k, zeta = 0), rep(1/3, 3))
})

test_that("proximity_allowance: empty active set gives uniform 1/K regardless of zeta", {
  W <- matrix(0, 3, 3)
  W[1, 2] <- W[2, 1] <- W[2, 3] <- W[3, 2] <- 1
  g <- ProximityGraph(W)
  k <- exponential_kernel(xi = 1)
  expect_equal(proximity_allowance(g, active = integer(0), kernel = k, zeta = 0.9), rep(1/3, 3))
})

test_that("proximity_allowance: all-zero kernel weights fall back to uniform 1/K", {
  W <- matrix(0, 3, 3)
  W[1, 2] <- W[2, 1] <- W[2, 3] <- W[3, 2] <- 1
  g <- ProximityGraph(W)
  zero_kernel <- function(d) rep(0, length(d))
  expect_equal(proximity_allowance(g, active = c(1), kernel = zero_kernel, zeta = 0.9),
               rep(1/3, 3))
})

test_that("proximity_allowance: worked example on a 3-node chain matches hand computation", {
  # Chain 1-2-3, unit weights: D[1,]=c(0,1,2), active={1}, kernel=exp(-d).
  # Node 1 is already active, so it competes only for the uniform floor
  # (1-zeta)/3; the entire zeta-bonus goes to nodes 2 and 3 by distance to 1.
  W <- matrix(0, 3, 3)
  W[1, 2] <- W[2, 1] <- 1
  W[2, 3] <- W[3, 2] <- 1
  g <- ProximityGraph(W)
  k <- exponential_kernel(xi = 1)
  zeta <- 0.5

  d_to_active <- c(0, 1, 2)          # distance of each node to node 1
  u <- exp(-d_to_active)
  u[1] <- 0                          # node 1 is active: excluded from the bonus
  expected <- (1 - zeta) / 3 + zeta * u / sum(u)

  out <- proximity_allowance(g, active = c(1), kernel = k, zeta = zeta)
  expect_equal(out, expected)
  expect_equal(sum(out), 1)
  # node 2 (closest non-active neighbor) gets the largest share; node 1
  # (active, baseline only) gets less than node 3 (farther, but still competing)
  expect_true(out[2] > out[3])
  expect_true(out[3] > out[1])
})

# ---- propagation_nu (Section 4.3) ----------------------------------------------

test_that("propagation_nu: linear chain propagation matches nu0 + scale * d_G", {
  W <- matrix(0, 4, 4)
  for (i in 1:3) W[i, i + 1] <- W[i + 1, i] <- 1
  g <- ProximityGraph(W)
  nu <- propagation_nu(g, source = 1, nu0 = 100, scale = 10)
  expect_equal(unname(nu), c(100, 110, 120, 130))
})

test_that("propagation_nu: hub-and-spoke propagation from the hub", {
  K <- 4
  W <- matrix(0, K, K)
  for (k in 2:K) W[1, k] <- W[k, 1] <- 2
  g <- ProximityGraph(W)
  nu <- propagation_nu(g, source = 1, nu0 = 50, scale = 5)
  expect_equal(unname(nu), c(50, 60, 60, 60))
})

test_that("propagation_nu: unreachable nodes get Inf (never change)", {
  W <- matrix(0, 3, 3)
  W[1, 2] <- W[2, 1] <- 1   # node 3 isolated
  g <- ProximityGraph(W)
  nu <- propagation_nu(g, source = 1, nu0 = 10, scale = 1)
  expect_equal(unname(nu[3]), Inf)
})
