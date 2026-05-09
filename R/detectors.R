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
  slots = c(alpha = "numeric", criterion = "character", threshold = "numeric",
            spending = "numeric", name = "character", multiple_alarms = "logical")
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
                                    multiple_alarms = FALSE,
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
      multiple_alarms = multiple_alarms,
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
                          multiple_alarms = FALSE,
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
      multiple_alarms = multiple_alarms,
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

# Class: LocalizedDetector
# purpose: K independent ShiryaevRoberts detectors with Bonferroni-corrected
#          levels (alpha/K each) for simultaneous per-stream localization.
#          Controls global ARL >= 1/alpha (via S-R martingale + Bonferroni) and
#          global PFA <= alpha (by union bound). Returns per-stream results so
#          that the specific stream(s) and time(s) of change can be identified.
# slots:
#   detectors = length-K list of ShiryaevRobertsDetector objects, each at level alpha/K
#   alpha     = numeric scalar global nominal level
#   criterion = character, one of c("ARL", "PFA")
setClass("LocalizedDetector",
         slots = c(detectors = "list", alpha = "numeric", criterion = "character"))

# Constructor: LocalizedDetector
# inputs:
#   K               = positive integer number of streams
#   alpha           = numeric scalar global nominal level
#   criterion       = character in c("ARL", "PFA")
#   spending        = optional numeric schedule for PFA recursion (passed to each marginal detector)
#   threshold       = optional numeric threshold; defaults to K/alpha (= 1/(alpha/K))
#   geometric_p     = geometric parameter for auto-generated PFA spending schedule
#   spending_length = length of auto-generated spending schedule
#   multiple_alarms = logical; whether each marginal detector resets and re-alarms
# outputs:
#   LocalizedDetector object
LocalizedDetector <- function(K,
                               alpha           = 0.05,
                               criterion       = c("ARL", "PFA"),
                               spending        = numeric(0),
                               threshold       = NULL,
                               geometric_p     = 0.01,
                               spending_length = 1000,
                               multiple_alarms = FALSE) {
  K <- as.integer(K)
  if (length(K) != 1L || K < 1L) {
    stop("`K` must be a positive integer (number of marginal streams).", call. = FALSE)
  }
  criterion <- match.arg(criterion)

  alpha_k <- alpha / K
  dets <- lapply(seq_len(K), function(k) {
    ShiryaevRobertsDetector(
      alpha           = alpha_k,
      criterion       = criterion,
      spending        = spending,
      threshold       = threshold,
      geometric_p     = geometric_p,
      spending_length = spending_length,
      multiple_alarms = multiple_alarms,
      name            = paste0("stream-", k)
    )
  })

  new("LocalizedDetector",
      detectors = dets,
      alpha     = as.numeric(alpha),
      criterion = criterion)
}

# Method: run_detector for LocalizedDetector
# inputs:
#   object   = LocalizedDetector object
#   evidence = numeric N-by-K matrix of per-stream increments (or log-increments if log=TRUE);
#              column k feeds the k-th marginal ShiryaevRobertsDetector
#   log      = logical; if TRUE interpret evidence as log-increments
# outputs:
#   list with elements:
#     stream_results     = named list of K individual run_detector outputs (statistic,
#                          stopping_time, alarm_times, alarm, criterion, log)
#     stopping_time      = global stopping time: min of marginal stopping times (Inf if none)
#     first_alarm_stream = integer index of the first-alarming stream (NA_integer_ if none)
#     alarm              = logical; TRUE if any stream alarmed
#     criterion          = character criterion used
#     log                = logical
setMethod("run_detector", "LocalizedDetector", function(object, evidence, log = FALSE) {
  if (!is.matrix(evidence) || !is.numeric(evidence)) {
    stop("`evidence` must be a numeric N-by-K matrix for a LocalizedDetector.", call. = FALSE)
  }
  K <- length(object@detectors)
  if (ncol(evidence) != K) {
    stop(sprintf("`evidence` must have %d column(s) (one per stream).", K), call. = FALSE)
  }
  if (!log && any(evidence < 0, na.rm = TRUE)) {
    stop("evidence has negative entries but log is FALSE.", call. = FALSE)
  }

  stream_results <- lapply(seq_len(K), function(k) {
    run_detector(object@detectors[[k]], evidence = evidence[, k], log = log)
  })
  names(stream_results) <- paste0("stream_", seq_len(K))

  stop_times   <- vapply(stream_results, function(r) r$stopping_time, numeric(1L))
  global_stop  <- min(stop_times)
  global_alarm <- is.finite(global_stop)
  first_stream <- if (global_alarm) unname(which(stop_times == global_stop)[1L]) else NA_integer_

  list(
    stream_results     = stream_results,
    stopping_time      = global_stop,
    first_alarm_stream = first_stream,
    alarm              = global_alarm,
    criterion          = object@criterion,
    log                = log
  )
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
  n          <- length(evidence)
  stat       <- numeric(n)
  state0     <- if (log && is(object, "ShiryaevRobertsDetector")) -Inf else 0
  state      <- state0
  alarm_times <- integer(0)
  spend_idx  <- 0L   # local spending offset; advances after each reset

  for (t in seq_len(n)) {
    # For PFA with multiple alarms: translate global t to recycled spending index.
    t_local <- if (object@criterion == "PFA") t - spend_idx else t
    step    <- update_detector(object, evidence = evidence[t], t = t_local, state = state, log = log)
    state   <- step$state
    stat[t] <- step$statistic

    if (isTRUE(step$alarm)) {
      if (object@multiple_alarms || length(alarm_times) == 0L) {
        alarm_times <- c(alarm_times, t)
      }
      if (object@multiple_alarms) {
        state     <- state0          # reset to initial value
        spend_idx <- t               # recycle spending from t+1 onward
      }
    }
  }

  stop_time <- if (length(alarm_times) > 0L) alarm_times[1L] else Inf
  list(
    statistic     = stat,
    stopping_time = stop_time,
    alarm_times   = alarm_times,
    alarm         = length(alarm_times) > 0L,
    criterion     = object@criterion,
    log           = log
  )
})
