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
#   mode       = "sparse" (K+1 portfolios: K unit vectors + uniform average) or
#                "dense"  (full simplex discretization at given resolution)
#   resolution = positive integer; simplex step count used only in dense mode
setClass("UniversalPortfolioCombiner", contains = "Combiner",
         slots = c(mode = "character", resolution = "numeric"))

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
#   mode       = "sparse" (default) or "dense"
#                "sparse": K+1 portfolios — the K unit vectors plus the uniform
#                          average {1/K,...,1/K}; O(K) portfolios, always cheap.
#                "dense":  full simplex discretization with `resolution` equally-
#                          spaced steps; choose(K-1+resolution, K-1) portfolios.
#                          Errors if that count exceeds 50 000.
#   resolution = positive integer controlling step count in dense mode (ignored
#                in sparse mode).  Default 4 gives 5 values per axis
#                {0, 0.25, 0.5, 0.75, 1}.
#   name       = character label
# outputs:
#   UniversalPortfolioCombiner object
UniversalPortfolioCombiner <- function(mode = "sparse", resolution = 4,
                                       name = "universal-portfolio") {
  stopifnot(mode %in% c("sparse", "dense"))
  stopifnot(resolution >= 1)
  new("UniversalPortfolioCombiner", name = name, mode = mode,
      resolution = as.numeric(resolution))
}

# Method: combine_streams for ProductCombiner
# inputs:
#   object  = ProductCombiner object
#   streams = numeric N-by-K matrix of marginal increment sequences
#   weights = ignored for this combiner
#   log     = logical; if TRUE combine log-increments and return log-increments
# outputs:
#   numeric length-N vector (row-wise products or sums of log-increments)
setMethod("combine_streams", "ProductCombiner", function(object, streams, weights = NULL, log = FALSE) {
  if (!is.matrix(streams) || !is.numeric(streams)) {
    stop("`streams` must be a numeric matrix of marginal increment sequences.", call. = FALSE)
  }
  # NA entries (offline streams) contribute 0 in log-space (factor 1), so omit them.
  if (log) {
    return(rowSums(streams, na.rm = TRUE))
  }
  apply(streams, 1, function(row) prod(row[!is.na(row)]))
})

# Method: combine_streams for AverageCombiner
# inputs:
#   object  = AverageCombiner object
#   streams = numeric N-by-K matrix of marginal increment sequences
#   weights = optional nonnegative numeric length-K vector (defaults flat)
#   log     = logical; if TRUE combine log-increments and return log-increments
# outputs:
#   numeric length-N vector of combined increments (or log-increments)
setMethod("combine_streams", "AverageCombiner", function(object, streams, weights = NULL, log = FALSE) {
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

  # Fast path: no NAs — use vectorised operations.
  if (!anyNA(streams)) {
    if (!log) {
      return(as.vector(streams %*% w))
    }
    log_w <- ifelse(w > 0, log(w), -Inf)
    out <- numeric(nrow(streams))
    for (t in seq_len(nrow(streams))) {
      out[t] <- .logsumexp(log_w + streams[t, ])
    }
    return(out)
  }

  # NA-aware path: at each t, renormalise weights over online (non-NA) streams
  # so that offline streams do not dilute the average.
  out <- numeric(nrow(streams))
  for (t in seq_len(nrow(streams))) {
    online <- !is.na(streams[t, ])
    w_t    <- w[online] / sum(w[online])
    if (!log) {
      out[t] <- sum(streams[t, online] * w_t)
    } else {
      log_w_t <- ifelse(w_t > 0, log(w_t), -Inf)
      out[t]  <- .logsumexp(log_w_t + streams[t, online])
    }
  }
  out
})

# Method: combine_streams for UniversalPortfolioCombiner
# inputs:
#   object  = UniversalPortfolioCombiner object
#   streams = numeric N-by-K matrix of marginal increment sequences
#   weights = currently unused (reserved for future extension)
#   log     = logical; if TRUE combine log-increments, return log-increments
# outputs:
#   numeric length-N vector of mixed-portfolio increments (or log-increments)
#
# Portfolio grid:
#   sparse mode: K+1 rows — diag(K) (unit vectors) + rep(1/K, K) (uniform)
#   dense  mode: full simplex discretization via .composition_grid();
#                errors if choose(K-1+resolution, K-1) > 50 000
setMethod("combine_streams", "UniversalPortfolioCombiner",
          function(object, streams, weights = NULL, log = FALSE) {
  if (!is.matrix(streams) || !is.numeric(streams)) {
    stop("`streams` must be a numeric matrix.", call. = FALSE)
  }

  n <- nrow(streams)
  k <- ncol(streams)

  if (object@mode == "sparse") {
    # K unit-vector portfolios + uniform average: exactly K+1 rows
    grid <- rbind(diag(k), matrix(rep(1 / k, k), nrow = 1L))
  } else {
    # dense: full simplex discretization
    res <- as.integer(object@resolution)
    n_portfolios <- choose(k - 1L + res, k - 1L)
    if (n_portfolios > 50000L) {
      stop(sprintf(
        paste0("Dense UP grid has %d portfolios (> 50 000). ",
               "Reduce `resolution` or K, or use mode = \"sparse\"."),
        n_portfolios), call. = FALSE)
    }
    raw_grid <- .composition_grid(k, res)
    grid     <- raw_grid / res
    # ensure unit vectors and uniform average are always included
    if (k > 1L) {
      baseline <- rbind(diag(k), matrix(rep(1 / k, k), nrow = 1L))
      grid <- unique(rbind(grid, baseline))
    }
  }

  m <- nrow(grid)

  if (!log) {
    wealth       <- rep(1, m)
    combined_inc <- numeric(n)
    for (t in seq_len(n)) {
      row_t <- streams[t, ]
      row_t[is.na(row_t)] <- 1   # offline stream contributes factor 1
      gross           <- pmax(as.vector(grid %*% row_t), .Machine$double.eps)
      combined_inc[t] <- sum(wealth * gross) / sum(wealth)
      wealth          <- wealth * gross
    }
    return(combined_inc)
  }

  # log-scale path
  log_grid        <- log(grid)
  log_wealth      <- rep(0, m)
  combined_log_inc <- numeric(n)
  for (t in seq_len(n)) {
    row_t <- streams[t, ]
    row_t[is.na(row_t)] <- 0   # offline stream contributes 0 in log-space
    log_gross <- numeric(m)
    for (i in seq_len(m)) {
      log_gross[i] <- .logsumexp(log_grid[i, ] + row_t)
    }
    combined_log_inc[t] <-
      .logsumexp(log_wealth + log_gross) - .logsumexp(log_wealth)
    log_wealth <- log_wealth + log_gross
  }
  combined_log_inc
})
