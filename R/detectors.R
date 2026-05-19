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
# purpose: K independent S-R recursions with a joint spending process (Theorem 3,
#          Spertus et al. 2026). Each stream k runs the recursion
#            R_tk = Lambda_tk * (R_{t-1,k} + invest_tk)
#          with a shared threshold 1/alpha, where invest_tk is drawn from the N_sched x K
#          `invest` matrix:
#            ARL (spending allowance gamma_t): sum_k gamma_tk = 1 at each t;
#                recycles the last row beyond N_sched.
#            PFA (joint spending schedule pi_t):  sum_t sum_k pi_tk = 1;
#                uses zero investment beyond N_sched.
#          The default is the Bonferroni allocation (Proposition 3): gamma_tk = 1/K for ARL,
#          pi_tk = pi_t_tilde / K for PFA where pi_tilde is a geometric schedule.
# slots:
#   K         = positive integer number of streams
#   alpha     = numeric scalar global nominal level
#   criterion = character, one of c("ARL", "PFA")
#   threshold = numeric scalar (default 1/alpha)
#   invest    = numeric N_sched x K matrix of per-stream investments
setClass("LocalizedDetector",
         slots = c(K         = "integer",
                   alpha     = "numeric",
                   criterion = "character",
                   threshold = "numeric",
                   invest    = "matrix"))

# Constructor: LocalizedDetector
# inputs:
#   K               = positive integer number of streams
#   alpha           = numeric scalar global nominal level
#   criterion       = character in c("ARL", "PFA")
#   allowance       = ARL only: NULL (uniform 1/K; Bonferroni), a length-K numeric vector
#                     (constant allocation), or an N x K matrix with rows summing to 1.
#   spending        = PFA only: NULL (geometric schedule split K ways), a length-N numeric
#                     vector (same schedule per stream, normalised and split K ways), or an
#                     N x K matrix with all entries summing to 1 (joint spending schedule).
#   threshold       = optional override (default 1/alpha)
#   geometric_p     = geometric parameter for auto-generated PFA spending schedule
#   spending_length = length of auto-generated PFA spending schedule
# outputs:
#   LocalizedDetector object
LocalizedDetector <- function(K,
                               alpha           = 0.05,
                               criterion       = c("ARL", "PFA"),
                               allowance       = NULL,
                               spending        = NULL,
                               threshold       = NULL,
                               geometric_p     = 0.01,
                               spending_length = 1000) {
  K <- as.integer(K)
  if (length(K) != 1L || K < 1L)
    stop("`K` must be a positive integer.", call. = FALSE)
  criterion <- match.arg(criterion)
  if (is.null(threshold)) threshold <- 1 / alpha

  if (criterion == "ARL") {
    if (is.null(allowance)) {
      invest <- matrix(rep(1 / K, K), nrow = 1L, ncol = K)
    } else {
      if (is.numeric(allowance) && !is.matrix(allowance))
        allowance <- matrix(allowance, nrow = 1L)
      invest <- as.matrix(allowance)
      if (ncol(invest) != K)
        stop("`allowance` must have K columns.", call. = FALSE)
      if (any(invest < -1e-12))
        stop("`allowance` must be non-negative.", call. = FALSE)
      if (any(abs(rowSums(invest) - 1) > 1e-8))
        stop("`allowance` rows must each sum to 1.", call. = FALSE)
    }
  } else {
    if (is.null(spending)) {
      sched  <- .geometric_spending(spending_length, p = geometric_p)
      sched  <- sched / sum(sched)
      invest <- matrix(rep(sched / K, K), nrow = spending_length, ncol = K)
    } else if (is.numeric(spending) && !is.matrix(spending)) {
      spending <- as.numeric(spending)
      spending <- spending / sum(spending)
      invest   <- matrix(rep(spending / K, K), nrow = length(spending), ncol = K)
    } else {
      invest <- as.matrix(spending)
      if (ncol(invest) != K)
        stop("`spending` must have K columns.", call. = FALSE)
      if (any(invest < -1e-12))
        stop("`spending` must be non-negative.", call. = FALSE)
      if (abs(sum(invest) - 1) > 1e-8)
        stop("`spending` must sum to 1.", call. = FALSE)
    }
  }

  new("LocalizedDetector",
      K         = K,
      alpha     = as.numeric(alpha),
      criterion = criterion,
      threshold = as.numeric(threshold),
      invest    = invest)
}

# Method: run_detector for LocalizedDetector
# inputs:
#   object   = LocalizedDetector object
#   evidence = numeric N-by-K matrix of per-stream increments (or log-increments if log=TRUE)
#   log      = logical; if TRUE interpret evidence as log-increments
# outputs:
#   list with elements:
#     stream_results     = named list of K per-stream results (statistic, stopping_time, alarm)
#     stopping_time      = global stopping time: min of per-stream stopping times (Inf if none)
#     first_alarm_stream = integer index of first-alarming stream (NA_integer_ if none)
#     alarm              = logical; TRUE if any stream alarmed
#     criterion          = character criterion used
#     log                = logical
setMethod("run_detector", "LocalizedDetector", function(object, evidence, log = FALSE) {
  if (!is.matrix(evidence) || !is.numeric(evidence))
    stop("`evidence` must be a numeric N-by-K matrix for a LocalizedDetector.", call. = FALSE)
  N <- nrow(evidence)
  K <- object@K
  if (ncol(evidence) != K)
    stop(sprintf("`evidence` must have %d column(s) (one per stream).", K), call. = FALSE)
  if (!log && any(evidence < 0, na.rm = TRUE))
    stop("evidence has negative entries but log is FALSE.", call. = FALSE)

  is_pfa     <- object@criterion == "PFA"
  invest     <- object@invest
  n_sched    <- nrow(invest)
  log_thresh <- log(object@threshold)

  log_stat  <- matrix(-Inf, nrow = N, ncol = K)
  log_state <- rep(-Inf, K)
  stop_times <- rep(Inf, K)

  for (t in seq_len(N)) {
    inv_t <- if (is_pfa) {
      if (t <= n_sched) invest[t, ] else rep(0, K)
    } else {
      invest[min(t, n_sched), ]
    }

    log_inc <- if (log) evidence[t, ] else log(pmax(evidence[t, ], .Machine$double.eps))
    log_inv <- ifelse(inv_t <= 0, -Inf, log(inv_t))

    # vectorised two-argument logsumexp: logsumexp(log_state, log_inv)
    mx      <- pmax(log_state, log_inv)
    log_sum <- ifelse(is.finite(mx), mx + log1p(exp(pmin(log_state, log_inv) - mx)), mx)

    log_state    <- log_inc + log_sum
    log_stat[t, ] <- log_state

    new_alarms <- which(log_state >= log_thresh & is.infinite(stop_times))
    if (length(new_alarms) > 0L) stop_times[new_alarms] <- t
  }

  global_stop  <- min(stop_times)
  global_alarm <- is.finite(global_stop)
  first_stream <- if (global_alarm) which.min(stop_times) else NA_integer_

  stream_results <- lapply(seq_len(K), function(k) {
    stat <- log_stat[, k]
    if (!log) stat <- exp(stat)
    list(statistic     = stat,
         stopping_time = stop_times[k],
         alarm         = is.finite(stop_times[k]),
         criterion     = object@criterion,
         log           = log)
  })
  names(stream_results) <- paste0("stream_", seq_len(K))

  list(stream_results     = stream_results,
       stopping_time      = global_stop,
       first_alarm_stream = first_stream,
       alarm              = global_alarm,
       criterion          = object@criterion,
       log                = log)
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
