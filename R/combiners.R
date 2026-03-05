# Class: Combiner
# purpose: base class for methods that combine K marginal evidence streams
# slots:
#   name = character label for reporting
setClass("Combiner", slots = c(name = "character"))

# Class: ProductCombiner
# purpose: combine by coordinate-wise product (typically for independent streams)
setClass("ProductCombiner", contains = "Combiner")

# Class: AverageCombiner
# purpose: combine by weighted average (robust and fast default)
setClass("AverageCombiner", contains = "Combiner")

# Class: UniversalPortfolioCombiner
# purpose: dependence-robust combiner via discrete universal portfolio mixing
# slots:
#   resolution = numeric scalar grid resolution on simplex
#   max_points = numeric scalar cap on number of grid portfolios
setClass("UniversalPortfolioCombiner", contains = "Combiner", slots = c(resolution = "numeric", max_points = "numeric"))

# Constructor: ProductCombiner
# inputs:
#   name = character label
# outputs:
#   ProductCombiner object
ProductCombiner <- function(name = "product") {
  new("ProductCombiner", name = name)
}

# Constructor: AverageCombiner
# inputs:
#   name = character label
# outputs:
#   AverageCombiner object
AverageCombiner <- function(name = "average") {
  new("AverageCombiner", name = name)
}

# Constructor: UniversalPortfolioCombiner
# inputs:
#   resolution = positive numeric scalar controlling simplex discretization
#   max_points = positive numeric scalar limiting portfolio count
#   name       = character label
# outputs:
#   UniversalPortfolioCombiner object
UniversalPortfolioCombiner <- function(resolution = 6, max_points = 300, name = "universal-portfolio") {
  stopifnot(resolution >= 1, max_points >= 1)
  new("UniversalPortfolioCombiner", name = name, resolution = resolution, max_points = max_points)
}

# Method: combine_streams for ProductCombiner
# inputs:
#   object  = ProductCombiner object
#   streams = numeric N-by-K matrix of marginal increment sequences
#   weights = ignored for this combiner
# outputs:
#   numeric length-N vector (row-wise products)
setMethod("combine_streams", "ProductCombiner", function(object, streams, weights = NULL) {
  if (!is.matrix(streams) || !is.numeric(streams)) {
    stop("`streams` must be a numeric matrix of marginal increment sequences.", call. = FALSE)
  }
  apply(streams, 1, prod)
})

# Method: combine_streams for AverageCombiner
# inputs:
#   object  = AverageCombiner object
#   streams = numeric N-by-K matrix of marginal increment sequences
#   weights = optional nonnegative numeric length-K vector (defaults flat)
# outputs:
#   numeric length-N vector (weighted averages)
setMethod("combine_streams", "AverageCombiner", function(object, streams, weights = NULL) {
  if (!is.matrix(streams) || !is.numeric(streams)) {
    stop("`streams` must be a numeric matrix of marginal increment sequences.", call. = FALSE)
  }
  k <- ncol(streams)
  if (is.null(weights)) {
    weights <- rep(1 / k, k)
  }
  if (length(weights) != k || any(weights < 0)) {
    stop("`weights` must be nonnegative and have one entry per stream.", call. = FALSE)
  }
  w <- weights / sum(weights)
  as.vector(streams %*% w)
})

# Method: combine_streams for UniversalPortfolioCombiner
# inputs:
#   object  = UniversalPortfolioCombiner object
#   streams = numeric N-by-K matrix of marginal increment sequences
#   weights = currently unused (reserved for future extension)
# outputs:
#   numeric length-N vector, mixed-portfolio one-step increments
setMethod("combine_streams", "UniversalPortfolioCombiner", function(object, streams, weights = NULL) {
  if (!is.matrix(streams) || !is.numeric(streams)) {
    stop("`streams` must be a numeric matrix of marginal increment sequences.", call. = FALSE)
  }

  n <- nrow(streams)
  k <- ncol(streams)
  raw_grid <- .composition_grid(k, as.integer(object@resolution))
  grid <- raw_grid / object@resolution

  if (nrow(grid) > object@max_points) {
    idx <- unique(round(seq(1, nrow(grid), length.out = object@max_points)))
    grid <- grid[idx, , drop = FALSE]
  }

  if (k > 1L) {
    baseline <- rbind(diag(k), rep(1 / k, k))
    grid <- unique(rbind(grid, baseline))
  }

  m <- nrow(grid)
  wealth <- rep(1, m)
  combined_inc <- numeric(n)

  for (t in seq_len(n)) {
    gross <- pmax(as.vector(grid %*% streams[t, ]), .Machine$double.eps)
    combined_inc[t] <- sum(wealth * gross) / sum(wealth)
    wealth <- wealth * gross
  }

  combined_inc
})
