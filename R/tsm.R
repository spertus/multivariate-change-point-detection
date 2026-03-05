# Class: TestSupermartingale
# purpose: base class for objects that convert data streams into evidence processes
# slots:
#   model   = Model subclass object used for likelihood increments
#   initial = numeric scalar, initial wealth/evidence
setClass("TestSupermartingale", slots = c(model = "Model", initial = "numeric"))

# Class: SimpleVsSimpleTSM
# purpose: concrete TSM for simple-vs-simple pre/post models
# slots: inherited from TestSupermartingale
setClass("SimpleVsSimpleTSM", contains = "TestSupermartingale")

# Constructor: SimpleVsSimpleTSM
# inputs:
#   model   = Model subclass instance
#   initial = numeric scalar > 0
# outputs:
#   SimpleVsSimpleTSM object
SimpleVsSimpleTSM <- function(model, initial = 1) {
  stopifnot(is(model, "Model"), length(initial) == 1L, initial > 0)
  new("SimpleVsSimpleTSM", model = model, initial = initial)
}

# Method: compute_tsm for SimpleVsSimpleTSM
# Method: compute_increments for SimpleVsSimpleTSM
# inputs:
#   object = SimpleVsSimpleTSM object
#   x      = numeric vector, observed sequence
# outputs:
#   numeric length-N vector with one-step increments
setMethod("compute_increments", "SimpleVsSimpleTSM", function(object, x) {
  .assert_numeric_vector(as.numeric(x), "x")
  x <- as.numeric(x)
  n <- length(x)
  out <- numeric(n)
  history <- numeric(0)

  for (t in seq_len(n)) {
    out[t] <- likelihood_increment(object@model, x = x[t], history = history)
    history <- c(history, x[t])
  }

  out
})

# Method: compute_tsm for SimpleVsSimpleTSM
# inputs:
#   object = SimpleVsSimpleTSM object
#   x      = numeric vector, observed sequence
# outputs:
#   numeric length-N vector with cumulative TSM path (running max of cumulative product)
setMethod("compute_tsm", "SimpleVsSimpleTSM", function(object, x) {
  inc <- compute_increments(object, x)
  increments_to_tsm(inc, initial = object@initial, running_max = TRUE)
})
