# Class: Detector
# purpose: base class for sequential change-point detector engines
# slots:
#   alpha     = numeric scalar nominal level
#   criterion = character, one of c("ARL", "PFA")
#   threshold = numeric scalar fixed threshold
#   spending  = numeric vector spending schedule (used only by S-R in PFA mode)
#   name      = character label
setClass(
  "Detector",
  slots = c(alpha = "numeric", criterion = "character", threshold = "numeric", spending = "numeric", name = "character")
)

# Class: ShiryaevRobertsDetector
# purpose: S-R recursion-based detector
setClass("ShiryaevRobertsDetector", contains = "Detector")

# Class: CUSUMDetector
# purpose: CUSUM recursion-based detector
setClass("CUSUMDetector", contains = "Detector")

# Constructor: ShiryaevRobertsDetector
# inputs:
#   alpha           = numeric scalar nominal level
#   criterion       = character in c("ARL", "PFA")
#   spending        = optional numeric schedule for PFA recursion injection
#   threshold       = optional numeric threshold; defaults to 1/alpha
#   geometric_p     = geometric parameter used only if criterion="PFA" and spending missing
#   spending_length = number of spending terms precomputed once at construction time
#   name            = character label
# outputs:
#   ShiryaevRobertsDetector object
ShiryaevRobertsDetector <- function(alpha = 0.05,
                                    criterion = c("ARL", "PFA"),
                                    spending = numeric(0),
                                    threshold = NULL,
                                    geometric_p = 0.01,
                                    spending_length = 1000,
                                    name = "shiryaev-roberts") {
  criterion <- match.arg(criterion)
  if (is.null(threshold)) {
    threshold <- 1 / alpha
  }

  if (criterion == "PFA") {
    if (length(spending) == 0L) {
      spending <- .geometric_spending(spending_length, p = geometric_p)
      spending <- spending / sum(spending)
    }

    if (any(spending < 0)) {
      stop("`spending` must be nonnegative.", call. = FALSE)
    }

    if (abs(sum(spending) - 1) > 1e-10) {
      stop("`spending` must sum to 1 over its support.", call. = FALSE)
    }
  }

  new("ShiryaevRobertsDetector",
      alpha = alpha,
      criterion = criterion,
      threshold = threshold,
      spending = spending,
      name = name)
}

# Constructor: CUSUMDetector
# inputs:
#   alpha     = numeric scalar nominal level
#   criterion = character in c("ARL", "PFA")
#   threshold = optional numeric threshold; defaults to log(1/alpha)
#   name      = character label
# outputs:
#   CUSUMDetector object
CUSUMDetector <- function(alpha = 0.05,
                          criterion = c("ARL", "PFA"),
                          threshold = NULL,
                          name = "cusum") {
  criterion <- match.arg(criterion)
  if (is.null(threshold)) {
    threshold <- log(1 / alpha)
  }
  new("CUSUMDetector",
      alpha = alpha,
      criterion = criterion,
      threshold = threshold,
      spending = numeric(0),
      name = name)
}

# Method: update_detector for ShiryaevRobertsDetector
# inputs:
#   object   = ShiryaevRobertsDetector object
#   evidence = numeric scalar increment (or log-increment if log=TRUE)
#   t        = integer time index
#   state    = numeric current S-R state (or log-state if log=TRUE)
#   log      = logical; if TRUE use log-scale recursion and return log-statistic
# outputs:
#   list with updated state, alarm flag, and statistic value at time t
setMethod("update_detector", "ShiryaevRobertsDetector", function(object, evidence, t, state, log = FALSE) {
  if (!log) {
    if (object@criterion == "ARL") {
      Rt <- (1 + state) * evidence
    } else {
      if (length(object@spending) < t) {
        pi_t <- 0
      } else {
        pi_t <- object@spending[t]
      }
      Rt <- (pi_t + state) * evidence
    }
    alarm <- Rt >= object@threshold
    return(list(state = Rt, alarm = alarm, statistic = Rt))
  }

  log_inc <- evidence
  if (object@criterion == "ARL") {
    log_add <- 0
  } else {
    if (length(object@spending) < t) {
      log_add <- -Inf
    } else {
      pi_t <- object@spending[t]
      log_add <- if (pi_t <= 0) -Inf else log(pi_t)
    }
  }

  logRt <- log_inc + .logsumexp(c(log_add, state))
  alarm <- logRt >= log(object@threshold)
  list(state = logRt, alarm = alarm, statistic = logRt)
})

# Method: update_detector for CUSUMDetector
# inputs:
#   object   = CUSUMDetector object
#   evidence = numeric scalar increment (or log-increment if log=TRUE)
#   t        = integer time index
#   state    = numeric current CUSUM state
#   log      = logical; if TRUE evidence is interpreted as log-increment
# outputs:
#   list with updated state, alarm flag, and statistic value at time t
setMethod("update_detector", "CUSUMDetector", function(object, evidence, t, state, log = FALSE) {
  log_inc <- if (log) evidence else log(pmax(evidence, .Machine$double.eps))
  St <- max(0, state + log_inc)
  if (object@criterion == "ARL") {
    alarm <- St >= object@threshold
  } else {
    stop("CUSUM not implemented with PFA control", call. = FALSE)
  }
  list(state = St, alarm = alarm, statistic = St)
})

# Method: run_detector for Detector
# inputs:
#   object   = Detector subclass object
#   evidence = numeric vector of one-step increments (or log-increments if log=TRUE)
#   log      = logical; if TRUE run recursion in log domain and return log-statistics
# outputs:
#   list with elements:
#     statistic     = numeric vector detector path over full input horizon
#     stopping_time = integer or Inf
#     alarm         = logical
#     criterion     = character criterion used
setMethod("run_detector", "Detector", function(object, evidence, log = FALSE) {
  if(!log && any(evidence < 0)){stop("evidence is negative but log is false!")}
  .assert_numeric_vector(evidence, "evidence")
  n <- length(evidence)
  stat <- numeric(n)
  state <- if (log && is(object, "ShiryaevRobertsDetector")) -Inf else 0
  stop_time <- Inf

  for (t in seq_len(n)) {
    step <- update_detector(object, evidence = evidence[t], t = t, state = state, log = log)
    state <- step$state
    stat[t] <- step$statistic
    if (isTRUE(step$alarm) && !is.finite(stop_time)) {
      stop_time <- t
    }
  }

  list(statistic = stat, stopping_time = stop_time, alarm = is.finite(stop_time), criterion = object@criterion, log = log)
})
