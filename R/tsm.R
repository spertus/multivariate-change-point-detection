# Class: TestSupermartingale
# purpose: base class for objects that convert data streams into evidence processes
# slots:
#   model   = Model subclass object used for likelihood increments
#   initial = numeric scalar, initial wealth/evidence
setClass("TestSupermartingale", slots = c(model = "Model", initial = "numeric"))

# Class: TSM
# purpose: generic TSM class that can be used with simple or composite models
setClass("TSM", contains = "TestSupermartingale")

# Backward-compatible alias class
setClass("SimpleVsSimpleTSM", contains = "TSM")

# Constructor: TSM
# inputs:
#   model   = Model subclass instance
#   initial = numeric scalar > 0
# outputs:
#   TSM object
TSM <- function(model, initial = 1) {
  stopifnot(is(model, "Model"), length(initial) == 1L, initial > 0)
  new("TSM", model = model, initial = initial)
}

# Backward-compatible constructor alias
SimpleVsSimpleTSM <- function(model, initial = 1) {
  TSM(model = model, initial = initial)
}

# Method: compute_increments for TSM
# inputs:
#   object = TSM object
#   x      = numeric vector (univariate stream) or matrix n-by-K (multivariate stream)
#   log    = logical; if TRUE return log-increments
# outputs:
#   numeric length-N vector with one-step increments (or log-increments)
setMethod("compute_increments", "TSM", function(object, x, log = FALSE) {
  if (is.null(dim(x))) {
    .assert_numeric_vector(as.numeric(x), "x")
    x <- as.numeric(x)
    n <- length(x)
    out <- numeric(n)
    history <- numeric(0)

    for (t in seq_len(n)) {
      out[t] <- likelihood_increment(object@model, x = x[t], history = history, log = log)
      history <- c(history, x[t])
    }
    return(out)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric vector or matrix.", call. = FALSE)
  }
  if (!is(object@model, "MultivariateModel")) {
    stop("Matrix input requires a model inheriting from `MultivariateModel`.", call. = FALSE)
  }

  n <- nrow(x)
  out <- numeric(n)
  history <- matrix(numeric(0), nrow = 0, ncol = ncol(x))

  for (t in seq_len(n)) {
    out[t] <- likelihood_increment(object@model, x = x[t, ], history = history, log = log)
    history <- rbind(history, x[t, , drop = FALSE])
  }

  out
})

# Method: compute_tsm for TSM
# inputs:
#   object = TSM object
#   x      = numeric vector or matrix sequence
#   log    = logical; if TRUE return log-TSM path
# outputs:
#   numeric length-N vector with cumulative TSM path (or log-TSM path)
setMethod("compute_tsm", "TSM", function(object, x, log = FALSE) {
  inc <- compute_increments(object, x, log = log)
  increments_to_tsm(inc, initial = object@initial, running_max = TRUE, log = log)
})
