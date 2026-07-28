# Class: DependencyGraph
# purpose: unweighted 0-1 graph recording which streams may be dependent (Section 4).
#          akl = 0 indicates streams k, l are mutually independent; akl = 1 indicates
#          they may be dependent. A clustered (disconnected) graph partitions V into
#          blocks that are mutually independent of one another.
# slots:
#   V = character length-K vector of stream labels
#   A = symmetric K-by-K 0/1 adjacency matrix, zero diagonal
setClass("DependencyGraph", slots = c(V = "character", A = "matrix"))

# Class: ProximityGraph
# purpose: edge-weighted graph used to model the change-point configuration and to
#          drive graph-structured spending (Sections 4.3, 5.3). Larger edge weight
#          akl means node k and l are farther apart (akl = 0 means no edge).
# slots:
#   V = character length-K vector of stream labels
#   W = symmetric K-by-K non-negative weight matrix, zero diagonal (0 = no edge)
setClass("ProximityGraph", slots = c(V = "character", W = "matrix"))

# ---- internal validation helpers -------------------------------------------

# Internal helper: .validate_graph_matrix
# purpose: shared validation for adjacency/weight matrices used by both graph classes
# inputs:
#   M          = numeric matrix to validate
#   label      = character, used in error messages ("adjacency"/"weight")
#   binary     = logical; if TRUE, require entries in {0, 1}
# outputs:
#   NULL (invisibly) on success; throws error on invalid input
.validate_graph_matrix <- function(M, label, binary = FALSE) {
  if (!is.matrix(M) || !is.numeric(M) || nrow(M) != ncol(M)) {
    stop(sprintf("`%s` matrix must be a square numeric matrix.", label), call. = FALSE)
  }
  if (!isTRUE(all.equal(M, t(M), check.attributes = FALSE))) {
    stop(sprintf("`%s` matrix must be symmetric.", label), call. = FALSE)
  }
  if (any(diag(M) != 0)) {
    stop(sprintf("`%s` matrix must have a zero diagonal (no self-edges).", label), call. = FALSE)
  }
  if (any(M < 0)) {
    stop(sprintf("`%s` matrix entries must be non-negative.", label), call. = FALSE)
  }
  if (binary && !all(M %in% c(0, 1))) {
    stop(sprintf("`%s` matrix entries must be 0 or 1.", label), call. = FALSE)
  }
  invisible(NULL)
}

# ---- constructors ------------------------------------------------------------

# Constructor: DependencyGraph
# inputs:
#   A      = optional symmetric K-by-K 0/1 adjacency matrix (zero diagonal)
#   blocks = optional list of integer index vectors partitioning 1:K into mutually
#            independent clusters; within a block, all pairs are marked dependent.
#            Exactly one of `A` / `blocks` must be supplied.
#   K      = integer number of nodes; required (and only used) with `blocks` unless
#            inferable as max(unlist(blocks))
#   names  = optional character length-K vector of node labels
# outputs:
#   DependencyGraph object
DependencyGraph <- function(A = NULL, blocks = NULL, K = NULL, names = NULL) {
  if (is.null(A) == is.null(blocks)) {
    stop("Exactly one of `A` or `blocks` must be supplied.", call. = FALSE)
  }
  if (is.null(A)) {
    if (is.null(K)) K <- max(unlist(blocks))
    K <- as.integer(K)
    A <- matrix(0, nrow = K, ncol = K)
    for (b in blocks) {
      b <- as.integer(b)
      if (length(b) > 1L) A[b, b] <- 1
    }
    diag(A) <- 0
  } else {
    K <- nrow(A)
  }
  .validate_graph_matrix(A, "adjacency", binary = TRUE)
  nm <- if (is.null(names)) as.character(seq_len(K)) else as.character(names)
  if (length(nm) != K) stop("`names` must have length K.", call. = FALSE)
  new("DependencyGraph", V = nm, A = A)
}

# Constructor: ProximityGraph
# inputs:
#   W     = symmetric K-by-K non-negative weight matrix (zero diagonal; 0 = no edge)
#   names = optional character length-K vector of node labels
# outputs:
#   ProximityGraph object
ProximityGraph <- function(W, names = NULL) {
  .validate_graph_matrix(W, "weight", binary = FALSE)
  K <- nrow(W)
  nm <- if (is.null(names)) as.character(seq_len(K)) else as.character(names)
  if (length(nm) != K) stop("`names` must have length K.", call. = FALSE)
  new("ProximityGraph", V = nm, W = W)
}

# ---- DependencyGraph: blocks --------------------------------------------------

# Function: dependency_blocks
# purpose: label the connected components (clusters) of a DependencyGraph via BFS
# inputs:
#   graph = DependencyGraph object
# outputs:
#   integer length-K vector of block labels (1, 2, ..., M)
dependency_blocks <- function(graph) {
  if (!is(graph, "DependencyGraph")) stop("`graph` must be a DependencyGraph.", call. = FALSE)
  A <- graph@A
  K <- nrow(A)
  block <- rep(0L, K)
  cur <- 0L
  for (start in seq_len(K)) {
    if (block[start] != 0L) next
    cur <- cur + 1L
    queue <- start
    block[start] <- cur
    while (length(queue) > 0L) {
      node  <- queue[1L]
      queue <- queue[-1L]
      nbrs  <- which(A[node, ] > 0 & block == 0L)
      if (length(nbrs) > 0L) {
        block[nbrs] <- cur
        queue <- c(queue, nbrs)
      }
    }
  }
  block
}

# ---- ProximityGraph: distances ------------------------------------------------

# Function: graph_distance_matrix
# purpose: all-pairs shortest-path distance via Floyd-Warshall over raw edge weights.
#          d_G(i, j) is the smallest sum of edge weights along a path from i to j;
#          Inf if no path exists (disconnected graph); 0 on the diagonal.
# inputs:
#   graph = ProximityGraph object
# outputs:
#   K-by-K numeric distance matrix, dimnamed by graph@V
graph_distance_matrix <- function(graph) {
  if (!is(graph, "ProximityGraph")) stop("`graph` must be a ProximityGraph.", call. = FALSE)
  W <- graph@W
  K <- nrow(W)
  D <- matrix(Inf, nrow = K, ncol = K)
  diag(D) <- 0
  edge <- W > 0
  D[edge] <- W[edge]
  for (m in seq_len(K)) {
    D <- pmin(D, outer(D[, m], D[m, ], "+"))
  }
  dimnames(D) <- list(graph@V, graph@V)
  D
}

# Function: graph_distance
# purpose: shortest-path distance between a single pair of nodes
# inputs:
#   graph = ProximityGraph object
#   i, j  = integer indices (or character labels matching graph@V)
# outputs:
#   numeric scalar distance (Inf if unreachable)
graph_distance <- function(graph, i, j) {
  graph_distance_matrix(graph)[i, j]
}

# ---- proximity kernels (Section 5.3) ------------------------------------------

# Constructor: exponential_kernel
# purpose: kappa(d) = exp(-d / xi); kappa(Inf) = 0
# inputs:
#   xi = positive numeric scalar bandwidth
# outputs:
#   vectorised function(d) -> numeric
exponential_kernel <- function(xi) {
  if (length(xi) != 1L || !is.finite(xi) || xi <= 0)
    stop("`xi` must be a positive finite scalar.", call. = FALSE)
  function(d) ifelse(is.infinite(d), 0, exp(-d / xi))
}

# Constructor: hard_cutoff_kernel
# purpose: kappa(d) = 1{d <= r}
# inputs:
#   r = non-negative numeric scalar cutoff radius
# outputs:
#   vectorised function(d) -> numeric (0/1)
hard_cutoff_kernel <- function(r) {
  if (length(r) != 1L || !is.finite(r) || r < 0)
    stop("`r` must be a non-negative finite scalar.", call. = FALSE)
  function(d) as.numeric(is.finite(d) & d <= r)
}

# ---- graph-structured adaptive allowance (Section 5.3) ------------------------

# Internal helper: .proximity_allowance_from_dist
# purpose: compute gamma_t from a precomputed distance matrix (hot-loop path, used by
#          LocalizedDetector so graph_distance_matrix() is not recomputed every step)
# inputs:
#   dist_mat = K-by-K numeric distance matrix (graph_distance_matrix() output)
#   active   = integer vector of currently-active (alarmed) node indices (A_{t-1});
#              length 0 = no active nodes
#   kernel   = vectorised proximity kernel function(d) -> numeric >= 0
#   zeta     = numeric scalar focus parameter in [0, 1)
# outputs:
#   numeric length-K allowance vector summing to 1
.proximity_allowance_from_dist <- function(dist_mat, active, kernel, zeta = 0) {
  K <- nrow(dist_mat)
  if (length(active) == 0L) return(rep(1 / K, K))
  d_to_active <- unname(apply(dist_mat[, active, drop = FALSE], 1, min))
  u <- unname(kernel(d_to_active))
  # An already-active node's own distance to A_{t-1} is 0 (the kernel's maximum),
  # so without this it would keep absorbing the largest share of the zeta-bonus
  # forever after firing -- budget better spent on not-yet-detected neighbors.
  u[active] <- 0
  if (sum(u) <= 0) return(rep(1 / K, K))
  (1 - zeta) / K + zeta * u / sum(u)
}

# Function: proximity_allowance
# purpose: compute the graph-structured spending allowance gamma_t (Section 5.3):
#            gamma_tk = (1 - zeta)/K + zeta * u_tk / sum(u_t),
#            u_tk = kappa(d_G(k, active)) * 1{k not in active}
#          Already-active nodes compete only for the uniform floor (1-zeta)/K; the
#          entire zeta-controlled bonus is redirected to not-yet-detected nodes by
#          distance to the active set. Convention: gamma_tk = 1/K when `active` is
#          empty, every node is active, or every kernel weight is zero.
# inputs:
#   graph  = ProximityGraph object
#   active = integer vector of currently-active (alarmed) node indices
#   kernel = vectorised proximity kernel function(d) -> numeric >= 0
#   zeta   = numeric scalar focus parameter in [0, 1); zeta = 0 recovers the uniform
#            (e-d-Bonferroni) allowance regardless of `active`/`kernel`
# outputs:
#   numeric length-K allowance vector summing to 1
proximity_allowance <- function(graph, active, kernel, zeta = 0) {
  dist_mat <- graph_distance_matrix(graph)
  .proximity_allowance_from_dist(dist_mat, active, kernel, zeta)
}

# ---- proximity model on the configuration (Section 4.3) -----------------------

# Function: propagation_nu
# purpose: build a change-point configuration nu propagating outward from a source
#          node along a ProximityGraph: nu_k = nu0 + scale * d_G(source, k).
#          Scaling `scale` toward 0 makes every reachable node change almost at once
#          (dense change); scaling it up makes only the source change within any
#          finite horizon (sparse change).
# inputs:
#   graph  = ProximityGraph object
#   source = integer index (or character label matching graph@V) of the source node
#   nu0    = numeric scalar change-point of the source node
#   scale  = positive numeric scalar; time units per unit of graph distance
# outputs:
#   numeric length-K vector of per-node change-points (Inf for unreachable nodes)
propagation_nu <- function(graph, source, nu0, scale = 1) {
  if (!is(graph, "ProximityGraph")) stop("`graph` must be a ProximityGraph.", call. = FALSE)
  d <- graph_distance_matrix(graph)[source, ]
  nu0 + scale * d
}
