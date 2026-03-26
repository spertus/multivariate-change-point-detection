# ---- Bets classes ----

# Class: Bets
# purpose: base class for predictable bet strategies used in bounded betting TSMs
# slots:
#   name = character label for reporting
setClass("Bets", slots = c(name = "character"))

# Class: FixedBets
# purpose: constant bet lambda_t = c for all t
# slots:
#   c = numeric scalar bet size in (0, 1]
setClass("FixedBets", contains = "Bets", slots = c(c = "numeric"))

# Class: AGRAPABets
# purpose: approximate GRAPA bets (Waudby-Smith and Ramdas, 2024, Section B.3)
#   lambda_t = clip((mu_{t-1} - eta) / (sd_{t-1}^2 + (mu_{t-1} - eta)^2), 0, c / eta)
# slots:
#   c      = numeric scalar truncation factor; bet is capped at c / eta
#   sd_min = numeric scalar floor for the lagged SD estimate
#   eps    = numeric scalar added to eta in the cap to avoid division by zero
setClass("AGRAPABets", contains = "Bets", slots = c(c = "numeric", sd_min = "numeric", eps = "numeric"))

# Class: PredictablePluginBets
# purpose: predictable plug-in estimator of Waudby-Smith and Ramdas (2024)
#   lambda_t = min(sqrt(2 * log(2 / alpha) / (sd_{t-1} * t * log(t + 1))), c)
# slots:
#   c      = numeric scalar maximum bet
#   alpha  = numeric scalar level parameter used in confidence-motivated tuning
#   sd_min = numeric scalar floor for the lagged SD estimate
setClass("PredictablePluginBets", contains = "Bets", slots = c(c = "numeric", alpha = "numeric", sd_min = "numeric"))

# ---- Bets constructors ----

# Constructor: FixedBets
# inputs:
#   c    = numeric scalar in (0, 1], the constant bet size
#   name = character label
# outputs:
#   FixedBets object
FixedBets <- function(c = 0.5, name = "fixed") {
  if (length(c) != 1L || !is.finite(c) || c <= 0 || c > 1) {
    stop("`c` must be a scalar in (0, 1].", call. = FALSE)
  }
  new("FixedBets", name = name, c = c)
}

# Constructor: AGRAPABets
# inputs:
#   c      = numeric scalar truncation factor in (0, 1], default 0.75 (Waudby-Smith and Ramdas)
#   sd_min = numeric scalar floor for lagged SD, default 0.01
#   eps    = numeric scalar floor added to eta in cap denominator, default 1e-5
#   name   = character label
# outputs:
#   AGRAPABets object
AGRAPABets <- function(c = 0.75, sd_min = 0.01, eps = 1e-5, name = "agrapa") {
  if (length(c) != 1L || !is.finite(c) || c <= 0 || c > 1) {
    stop("`c` must be a scalar in (0, 1].", call. = FALSE)
  }
  if (length(sd_min) != 1L || sd_min <= 0) {
    stop("`sd_min` must be a positive scalar.", call. = FALSE)
  }
  if (length(eps) != 1L || eps < 0) {
    stop("`eps` must be a non-negative scalar.", call. = FALSE)
  }
  new("AGRAPABets", name = name, c = c, sd_min = sd_min, eps = eps)
}

# Constructor: PredictablePluginBets
# inputs:
#   c      = numeric scalar maximum bet in (0, 1], default 0.9
#   alpha  = numeric scalar in (0, 1), level parameter, default 0.05
#   sd_min = numeric scalar floor for lagged SD, default 0.01
#   name   = character label
# outputs:
#   PredictablePluginBets object
PredictablePluginBets <- function(c = 0.9, alpha = 0.05, sd_min = 0.01, name = "predictable-plugin") {
  if (length(c) != 1L || !is.finite(c) || c <= 0 || c > 1) {
    stop("`c` must be a scalar in (0, 1].", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a scalar in (0, 1).", call. = FALSE)
  }
  if (length(sd_min) != 1L || sd_min <= 0) {
    stop("`sd_min` must be a positive scalar.", call. = FALSE)
  }
  new("PredictablePluginBets", name = name, c = c, alpha = alpha, sd_min = sd_min)
}

# ---- compute_bets methods ----

# Method: compute_bets for FixedBets
# inputs:
#   object = FixedBets object
#   x      = numeric vector of observations in [0, 1]
#   eta    = numeric scalar null mean (unused; kept for generic compatibility)
# outputs:
#   numeric length-N vector of constant bets equal to object@c
setMethod("compute_bets", "FixedBets", function(object, x, eta) {
  rep(object@c, length(x))
})

# Method: compute_bets for AGRAPABets
# purpose: approximate Kelly-optimal bets via lagged mean and SD
# inputs:
#   object = AGRAPABets object
#   x      = numeric vector of observations in [0, 1]
#   eta    = numeric scalar null mean in (0, 1]
# outputs:
#   numeric length-N vector of AGRAPA bets; predictable and non-negative
setMethod("compute_bets", "AGRAPABets", function(object, x, eta) {
  ls      <- .lag_mean_sd(x)
  lag_mu  <- ls$mean
  lag_sd  <- pmax(ls$sd, object@sd_min)
  lam_untrunc <- (lag_mu - eta) / (lag_sd^2 + (lag_mu - eta)^2)
  cap <- if (eta > 0) object@c / (eta + object@eps) else Inf
  pmax(0, pmin(lam_untrunc, cap))
})

# Method: compute_bets for PredictablePluginBets
# purpose: confidence-sequence-motivated predictable bet schedule
# inputs:
#   object = PredictablePluginBets object
#   x      = numeric vector of observations in [0, 1]
#   eta    = numeric scalar null mean (unused; kept for generic compatibility)
# outputs:
#   numeric length-N vector of predictable plug-in bets, clipped at object@c
setMethod("compute_bets", "PredictablePluginBets", function(object, x, eta) {
  ls      <- .lag_mean_sd(x)
  lag_sd  <- pmax(ls$sd, object@sd_min)
  t       <- seq_along(x)
  lam_untrunc <- sqrt((2 * log(2 / object@alpha)) / (lag_sd * t * log(t + 1L)))
  pmin(lam_untrunc, object@c)
})

# ---- BettingTSM class ----

# Class: BettingTSM
# purpose: test supermartingale for [0,1]-bounded observations using predictable bets
#   M_t = initial * prod_{i=1}^{t} [1 + lambda_i * (X_i - eta)]
#   where lambda_i is a predictable bet depending only on X_1, ..., X_{i-1}
# slots:
#   bets    = Bets subclass object defining the betting strategy
#   eta     = numeric scalar null mean; upper bound on E[X_t] under H_0
#   initial = numeric scalar initial wealth (default 1)
setClass(
  "BettingTSM",
  slots = c(bets = "Bets", eta = "numeric", initial = "numeric")
)

# Constructor: BettingTSM
# inputs:
#   bets    = Bets subclass object
#   eta     = numeric scalar in (0, 1], null mean upper bound
#   initial = numeric scalar > 0, initial wealth
# outputs:
#   BettingTSM object
BettingTSM <- function(bets, eta, initial = 1) {
  if (!is(bets, "Bets")) {
    stop("`bets` must be a Bets subclass object.", call. = FALSE)
  }
  if (length(eta) != 1L || !is.finite(eta) || eta <= 0 || eta > 1) {
    stop("`eta` must be a scalar in (0, 1].", call. = FALSE)
  }
  if (length(initial) != 1L || !is.finite(initial) || initial <= 0) {
    stop("`initial` must be a positive finite scalar.", call. = FALSE)
  }
  new("BettingTSM", bets = bets, eta = eta, initial = initial)
}

# Method: compute_increments for BettingTSM
# purpose: compute one-step wealth increments 1 + lambda_t * (x_t - eta)
# inputs:
#   object = BettingTSM object
#   x      = numeric vector of observations in [0, 1]
#   log    = logical; if TRUE return log-increments
# outputs:
#   numeric length-N vector of increments (or log-increments if log=TRUE)
setMethod("compute_increments", "BettingTSM", function(object, x, log = FALSE) {
  x <- as.numeric(x)
  .assert_numeric_vector(x, "x")
  lam <- compute_bets(object@bets, x, eta = object@eta)
  inc <- pmax(1 + lam * (x - object@eta), .Machine$double.eps)
  if (log) base::log(inc) else inc
})

# Method: compute_tsm for BettingTSM
# purpose: compute cumulative wealth path M_t = initial * prod_{i=1}^{t} increment_i
# inputs:
#   object = BettingTSM object
#   x      = numeric vector of observations in [0, 1]
#   log    = logical; if TRUE return log-TSM path
# outputs:
#   numeric length-N vector of cumulative TSM values (or log-TSM if log=TRUE)
setMethod("compute_tsm", "BettingTSM", function(object, x, log = FALSE) {
  inc <- compute_increments(object, x, log = log)
  increments_to_tsm(inc, initial = object@initial, log = log)
})
