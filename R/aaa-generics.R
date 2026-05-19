# Generic: model_density
# purpose: evaluate pre/post model density (or mass) at new observation(s)
# inputs:
#   object  = Model subclass instance
#   x       = numeric scalar/vector, new observation(s)
#   regime  = character, one of c("pre", "post")
#   history = optional numeric vector of past observations (for conditional models)
# outputs:
#   numeric scalar/vector of densities or probabilities
setGeneric("model_density", function(object, x, regime = c("pre", "post"), history = NULL) {
  standardGeneric("model_density")
})

# Generic: likelihood_increment
# purpose: compute per-time likelihood ratio increment between post and pre models
# inputs:
#   object  = Model subclass instance
#   x       = numeric scalar/vector, new observation(s)
#   history = optional numeric vector of past observations
#   log     = logical; if TRUE return log-likelihood increments
# outputs:
#   numeric scalar/vector of increments (or log-increments if log=TRUE)
setGeneric("likelihood_increment", function(object, x, history = NULL, log = FALSE) {
  standardGeneric("likelihood_increment")
})

# Generic: compute_increments
# purpose: build one-step likelihood/wealth increments from data
# inputs:
#   object = TestSupermartingale subclass instance
#   x      = numeric vector, observed stream
#   log    = logical; if TRUE return log-increments
# outputs:
#   numeric length-N vector of increments (or log-increments if log=TRUE)
setGeneric("compute_increments", function(object, x, log = FALSE) {
  standardGeneric("compute_increments")
})

# Generic: compute_tsm
# purpose: build a TSM path from data (helper for plotting/diagnostics)
# inputs:
#   object = TestSupermartingale subclass instance
#   x      = numeric vector, observed stream
#   log    = logical; if TRUE return log-TSM path
# outputs:
#   numeric length-N vector of cumulative TSM values (or log-TSM if log=TRUE)
setGeneric("compute_tsm", function(object, x, log = FALSE) {
  standardGeneric("compute_tsm")
})

# Generic: update_detector
# purpose: single-step update for detector state/statistic/alarm status
# inputs:
#   object   = Detector subclass instance
#   evidence = numeric scalar increment at time t
#   t        = integer time index
#   state    = numeric current detector state
#   log      = logical; whether evidence/state are on log scale
# outputs:
#   list with elements state, alarm (logical), statistic (numeric)
setGeneric("update_detector", function(object, evidence, t, state, log = FALSE) {
  standardGeneric("update_detector")
})

# Generic: run_detector
# purpose: run a detector across a full evidence sequence
# inputs:
#   object   = Detector subclass instance
#   evidence = numeric vector of increments
#   log      = logical; if TRUE interpret evidence as log-increments and return log-statistic
# outputs:
#   list with statistic path, stopping_time, alarm flag, and criterion name
setGeneric("run_detector", function(object, evidence, log = FALSE) {
  standardGeneric("run_detector")
})

# Generic: combine_streams
# purpose: combine K marginal increment sequences into one joint increment sequence
# inputs:
#   object  = Combiner subclass instance
#   streams = numeric N-by-K matrix of marginal increments
#   weights = optional numeric length-K vector
#   log     = logical; if TRUE combine log-increments and return log-increments
# outputs:
#   numeric length-N combined increment sequence
setGeneric("combine_streams", function(object, streams, weights = NULL, log = FALSE) {
  standardGeneric("combine_streams")
})

# Generic: generate_stream
# purpose: generate synthetic stream(s) from a DGP object
# inputs:
#   object = DGP instance
#   N      = integer horizon
#   K      = integer number of streams
# outputs:
#   numeric vector (K=1) or numeric N-by-K matrix (K>1)
setGeneric("generate_stream", function(object, N, K = 1L) {
  standardGeneric("generate_stream")
})
