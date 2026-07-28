# localization_graph_builders.R
# Shared graph-topology builders for the localization simulation
# (localization_graph_simulations.R) and its topology figures
# (ssc-2026/figures.R). Pure function definitions only, no top-level
# execution -- safe to source() from either place without side effects.
#
# Four topologies, each usable at any K:
#   linear           -- a path 1--2--...--K
#   fully_connected  -- every pair of nodes adjacent
#   hub_spoke        -- one or more "star" clusters (see hub_spoke_graph below)
#   clustered_linear -- n_chains disconnected linear chains (see below)
# All edge weights are unit.

# Function: multi_hub_graph
# purpose: n_hubs star clusters of spokes_per_hub spokes each. Hubs are nodes
#          1..n_hubs; hub h's spokes are the next spokes_per_hub node indices
#          in sequence (hub 1's spokes, then hub 2's spokes, ...). Hubs can
#          optionally be chained together.
# inputs:
#   n_hubs         = positive integer number of hub nodes
#   spokes_per_hub = positive integer number of spokes attached to each hub
#   hub_connected  = "disconnected" (hubs share no edges -- the graph is
#                    n_hubs separate star components) or "line" (hubs
#                    chained hub_1 -- hub_2 -- ... -- hub_{n_hubs})
# outputs:
#   ProximityGraph with K = n_hubs * (1 + spokes_per_hub) nodes
multi_hub_graph <- function(n_hubs, spokes_per_hub,
                            hub_connected = c("disconnected", "line")) {
  hub_connected <- match.arg(hub_connected)
  K    <- n_hubs * (1L + spokes_per_hub)
  W    <- matrix(0, K, K)
  hubs <- seq_len(n_hubs)
  for (h in hubs) {
    spoke_start <- n_hubs + (h - 1L) * spokes_per_hub + 1L
    spoke_end   <- spoke_start + spokes_per_hub - 1L
    for (s in spoke_start:spoke_end) W[h, s] <- W[s, h] <- 1
  }
  if (hub_connected == "line" && n_hubs > 1L) {
    for (h in seq_len(n_hubs - 1L)) {
      W[hubs[h], hubs[h + 1L]] <- W[hubs[h + 1L], hubs[h]] <- 1
    }
  }
  ProximityGraph(W)
}

# Function: hub_spoke_graph
# purpose: hub-and-spoke topology at either supported size.
#          K = 6:  a single hub (node 1) with 5 spokes (nodes 2-6);
#                  `hub_connected` is irrelevant (only one hub, nothing to
#                  connect it to).
#          K = 24: 4 hubs (nodes 1-4) with 5 spokes each (nodes 5-24), hubs
#                  disconnected (4 separate star components) or chained in a
#                  line (hub_1--hub_2--hub_3--hub_4) per `hub_connected`.
# inputs:
#   K             = 6 or 24
#   hub_connected = "disconnected" or "line" (only consulted when K == 24)
# outputs:
#   ProximityGraph
hub_spoke_graph <- function(K, hub_connected = "disconnected") {
  if (K == 6L)  return(multi_hub_graph(n_hubs = 1L, spokes_per_hub = 5L))
  if (K == 24L) return(multi_hub_graph(n_hubs = 4L, spokes_per_hub = 5L,
                                       hub_connected = hub_connected))
  stop("hub_spoke_graph currently supports K in {6, 24}.", call. = FALSE)
}

# Function: fully_connected_graph
# purpose: complete graph on K nodes, unit edge weights.
fully_connected_graph <- function(K) {
  W <- matrix(1, K, K)
  diag(W) <- 0
  ProximityGraph(W)
}

# Function: linear_graph
# purpose: path graph 1--2--...--K, unit edge weights.
linear_graph <- function(K) {
  W <- matrix(0, K, K)
  for (k in seq_len(K - 1L)) W[k, k + 1L] <- W[k + 1L, k] <- 1
  ProximityGraph(W)
}

# Function: clustered_linear_graph
# purpose: n_chains disconnected linear chains of equal length (K / n_chains
#          each). Chain c occupies consecutive node indices: chain 1 is nodes
#          1..L, chain 2 is L+1..2L, etc. (L = K / n_chains); no edges between
#          chains. The chain analogue of hub_spoke's "disconnected" clusters:
#          a change starting in chain 1 can propagate along that chain but can
#          never reach the other n_chains - 1 chains.
# inputs:
#   K        = integer graph size, must be divisible by n_chains
#   n_chains = positive integer number of disconnected chains
# outputs:
#   ProximityGraph
clustered_linear_graph <- function(K, n_chains = 4L) {
  if (K %% n_chains != 0L) stop("K must be divisible by n_chains.", call. = FALSE)
  L <- K %/% n_chains
  W <- matrix(0, K, K)
  for (c in seq_len(n_chains)) {
    offset <- (c - 1L) * L
    for (k in seq_len(L - 1L)) {
      i <- offset + k; j <- offset + k + 1L
      W[i, j] <- W[j, i] <- 1
    }
  }
  ProximityGraph(W)
}

# Function: build_graph
# purpose: single dispatch point used by both the simulation and its figures,
#          so the graphs actually simulated and the graphs drawn can never
#          drift apart.
# inputs:
#   graph_type    = "hub_spoke" | "fully_connected" | "linear" | "clustered_linear"
#   K             = integer graph size (6 or 24)
#   hub_connected = "disconnected" | "line"; only consulted for
#                   graph_type == "hub_spoke" and K == 24 (ignored otherwise)
# outputs:
#   ProximityGraph
build_graph <- function(graph_type, K, hub_connected = "disconnected") {
  switch(graph_type,
        hub_spoke        = hub_spoke_graph(K, hub_connected),
        fully_connected  = fully_connected_graph(K),
        linear           = linear_graph(K),
        clustered_linear = clustered_linear_graph(K),
        stop("unknown graph_type: ", graph_type, call. = FALSE))
}

# Function: pick_source
# purpose: source node for the propagation model (Section 4.3). The source is
#          never fixed to a particular node -- at the moment it is chosen, no
#          node has yet changed, so it is drawn uniformly at random from all K
#          nodes for every topology (hub, spoke, chain end, or otherwise).
#          graph_type is accepted for interface symmetry with build_graph()
#          but does not affect the draw.
pick_source <- function(graph_type, K) {
  sample.int(K, 1L)
}
