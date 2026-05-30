# figures.R --- explainer figures for the SSC 2026 talk
#
# Uses the `multichangepoints` package for all TSM and Shiryaev-Roberts
# computations (no hand-rolled recursions). For the motivating wastewater
# panel it reads the real CSV that ships with the package.
#
# Run from the ssc-2026/ folder:
#   Rscript figures.R
#
# Output: figures/*.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ----------------------------------------------------------------------
# Locate the package and load it via devtools::load_all (matches vignettes)
# ----------------------------------------------------------------------
PKG_DIR <- normalizePath(
  file.path(".."),
  mustWork = FALSE
)
if (!dir.exists(PKG_DIR)) {
  PKG_DIR <- "~/Dropbox/CANSSI/multi-change-points"
}
if (requireNamespace("devtools", quietly = TRUE)) {
  suppressMessages(devtools::load_all(PKG_DIR))
} else {
  library(multichangepoints)
}

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

theme_talk <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# ======================================================================
# 0. Intro: the change-point problem
#    Single Gaussian stream, true change-point nu, detection D, delay highlighted
# ======================================================================
set.seed(5)
nu_intro <- 80
D_intro  <- 115
N_intro  <- 180

x_intro <- c(rnorm(nu_intro, mean = 0,   sd = 1),
             rnorm(N_intro - nu_intro, mean = 2, sd = 1))
df_intro <- data.frame(
  t      = seq_len(N_intro),
  x      = x_intro,
  period = factor(ifelse(seq_len(N_intro) <= nu_intro, "pre", "post"),
                  levels = c("pre", "post"))
)

y_top    <- max(x_intro) * 0.92
y_bracket <- min(x_intro) - 0.5

p0 <- ggplot(df_intro, aes(t, x)) +
  # shade detection delay
  annotate("rect", xmin = nu_intro, xmax = D_intro, ymin = -Inf, ymax = Inf,
           fill = "#f4a261", alpha = 0.25) +
  # data, coloured by pre/post change
  geom_point(aes(colour = period), size = 1.2, alpha = 0.8) +
  scale_colour_manual(values = c(pre = "#457b9d", post = "#e63946"),
                      labels = c(pre = "pre-change", post = "post-change"),
                      name   = NULL) +
  # true change-point
  geom_vline(xintercept = nu_intro, linetype = "dashed", colour = "grey20", linewidth = 0.8) +
  annotate("text", x = nu_intro - 2, y = y_top,
           label = "nu", parse = TRUE, hjust = 1, size = 6, colour = "grey20") +
  # detection time
  geom_vline(xintercept = D_intro, colour = "#2a9d8f", linewidth = 1) +
  annotate("text", x = D_intro + 2, y = y_top,
           label = "italic(D)", parse = TRUE, hjust = 0, size = 6, colour = "#2a9d8f") +
  # detection delay bracket (two-headed arrow)
  annotate("segment",
           x = nu_intro, xend = D_intro,
           y = y_bracket, yend = y_bracket,
           arrow = arrow(ends = "both", length = unit(0.10, "inches")),
           colour = "#e07020", linewidth = 0.9) +
  annotate("text", x = (nu_intro + D_intro) / 2, y = y_bracket,
           label = "detection delay", vjust = -0.55, size = 3.8, colour = "#e07020") +
  labs(x = "time", y = expression(X[t]),
       title = "The change-point problem in one dimension") +
  theme_talk +
  theme(legend.position = "top")

ggsave("figures/intro_changepoint.png", p0, width = 8, height = 4.5, dpi = 220)

# ======================================================================
# 0b. ARL intuition: 5 false alarms before nu = 100, run lengths near 20
# ======================================================================
set.seed(17)
nu_arl       <- 100
N_arl        <- 140

# Hand-chosen false alarm times: 5 alarms, run lengths 18, 19, 18, 21, 19
# (mean = 19 ≈ 20; they are clearly not all equal)
false_alarms <- c(18, 37, 55, 76, 95)
run_starts   <- c(1, head(false_alarms, -1) + 1)   # 1, 19, 38, 56, 77
rls          <- false_alarms - run_starts + 1       # 18, 19, 18, 21, 19

x_arl <- c(rnorm(nu_arl, mean = 0, sd = 1),
           rnorm(N_arl - nu_arl, mean = 2, sd = 1))
df_arl <- data.frame(
  t      = seq_len(N_arl),
  x      = x_arl,
  period = factor(ifelse(seq_len(N_arl) <= nu_arl, "pre", "post"),
                  levels = c("pre", "post"))
)

# Layout heights (computed after data are known)
y_hi    <- max(x_arl)
y_brack <- y_hi + 0.55    # horizontal bracket line
y_rl    <- y_brack + 0.28 # run-length label
y_plot  <- y_rl   + 0.25  # top of coord window
y_nu    <- y_hi * 0.78    # nu label inside data area

# Bracket data frame (one row per run-length interval)
brack_df <- data.frame(
  x    = run_starts + 0.4,
  xend = false_alarms - 0.4,
  y    = y_brack,
  yend = y_brack,
  xmid = (run_starts + false_alarms) / 2,
  rl   = as.character(rls)
)

p0b <- ggplot(df_arl, aes(t, x)) +
  geom_point(aes(colour = period), size = 1.2, alpha = 0.8) +
  scale_colour_manual(values = c(pre = "#457b9d", post = "#e63946"),
                      labels = c(pre = "pre-change", post = "post-change"),
                      name = NULL) +
  # false alarm vertical lines
  geom_vline(xintercept = false_alarms, colour = "#e07020",
             linewidth = 0.65, linetype = "solid", alpha = 0.9) +
  # run-length brackets: two-headed arrows
  geom_segment(data = brack_df,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(ends = "both", length = unit(0.07, "inches")),
               colour = "#e07020", linewidth = 0.65, inherit.aes = FALSE) +
  # run-length numbers above brackets
  geom_text(data = brack_df,
            aes(x = xmid, y = yend + 0.25, label = rl),
            colour = "#e07020", size = 3.7, fontface = "bold", inherit.aes = FALSE) +
  # label one of the orange lines so audiences know what they are
  annotate("text", x = false_alarms[3] + 1.5, y = -2.3,
           label = "false alarm", colour = "#e07020",
           size = 3.3, hjust = 0, fontface = "italic") +
  annotate("segment", x = false_alarms[3] + 1.4, xend = false_alarms[3] + 0.2,
           y = -2.3, yend = -1.6,
           arrow = arrow(length = unit(0.08, "inches")),
           colour = "#e07020", linewidth = 0.5) +
  # true change-point
  geom_vline(xintercept = nu_arl, linetype = "dashed",
             colour = "grey20", linewidth = 0.9) +
  annotate("text", x = nu_arl + 2, y = y_nu,
           label = "nu", parse = TRUE, hjust = 0, size = 6, colour = "grey20") +
  coord_cartesian(ylim = c(NA, y_plot), clip = "off") +
  labs(x = "time", y = expression(X[t]),
       title = "Average run length (ARL) = expected time between false alarms",
       subtitle = "ARL = 20  →  expect ν / ARL = 100 / 20 = 5 false alarms before ν") +
  theme_talk +
  theme(legend.position = "top")

ggsave("figures/arl_intuition.png", p0b, width = 8, height = 4.8, dpi = 220)

# ======================================================================
# 0c. Graph-adaptive spending allowance: before and after one node fires
#     K = 6, star-of-stars topology: node 1 connects to nodes 2 and 3,
#     which in turn connect to nodes 4, 5, and 6 respectively.
#     Before: γ_k = 1/6  for all k.
#     After node 1 fires: γ_1 = 0, γ_2 = γ_3 = 2/7 (neighbors), γ_4 = γ_5 = γ_6 = 1/7.
# ======================================================================
K_g <- 6L

# Node positions (hand-placed for a readable two-row layout)
node_pos <- data.frame(
  node = 1:K_g,
  x    = c(2, 1, 3, 0, 2, 4),
  y    = c(3, 1.5, 1.5, 0, 0, 0)
)

# Edge list (undirected)
edges_raw <- data.frame(
  from = c(1, 1, 2, 2, 3),
  to   = c(2, 3, 4, 5, 6)
)

# Expand each edge into from/to coordinates for geom_segment
edges_coords <- function(pos, edge_df) {
  merge(edge_df, pos, by.x = "from", by.y = "node") |>
    (\(d) { names(d)[names(d) == "x"] <- "x_from"; names(d)[names(d) == "y"] <- "y_from"; d })() |>
    merge(pos, by.x = "to", by.y = "node") |>
    (\(d) { names(d)[names(d) == "x"] <- "x_to"; names(d)[names(d) == "y"] <- "y_to"; d })()
}

# Panel label strings (keep short so factor ordering is controlled explicitly)
panel_before <- "Before: uniform"
panel_after  <- "After node 1 fires"

# ---- Before: all blue, uniform allowance ----
nodes_before <- node_pos
nodes_before$role    <- "other"
nodes_before$gamma   <- "1/6"
nodes_before$panel   <- panel_before

edges_before         <- edges_coords(node_pos, edges_raw)
edges_before$panel   <- panel_before
edges_before$fired   <- FALSE

# ---- After node 1 fires: red/yellow/blue ----
nodes_after <- node_pos
nodes_after$role  <- ifelse(nodes_after$node == 1, "fired",
                    ifelse(nodes_after$node %in% c(2, 3), "neighbor", "other"))
nodes_after$gamma <- ifelse(nodes_after$node == 1, "0",
                    ifelse(nodes_after$node %in% c(2, 3), "2/7", "1/7"))
nodes_after$panel <- panel_after

edges_after       <- edges_coords(node_pos, edges_raw)
edges_after$panel <- panel_after
edges_after$fired <- edges_after$from == 1 | edges_after$to == 1

# ---- Combine both panels (Before left, After right) ----
nodes_g <- rbind(nodes_before, nodes_after)
nodes_g$role  <- factor(nodes_g$role, levels = c("fired", "neighbor", "other"))
nodes_g$panel <- factor(nodes_g$panel, levels = c(panel_before, panel_after))

edges_g <- rbind(edges_before, edges_after)
edges_g$panel <- factor(edges_g$panel, levels = c(panel_before, panel_after))

role_colours <- c(fired = "#e63946", neighbor = "#f4a261", other = "#457b9d")

p0c <- ggplot() +
  # Edges: grey normally, red when they touch the fired node (after-panel only)
  geom_segment(data = edges_g,
               aes(x = x_from, xend = x_to, y = y_from, yend = y_to,
                   colour = fired, linewidth = fired),
               lineend = "round", show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#e63946"),
                      guide = "none") +
  scale_linewidth_manual(values = c("FALSE" = 0.8, "TRUE" = 1.8),
                         guide = "none") +
  # Nodes
  geom_point(data = nodes_g,
             aes(x, y, fill = role),
             shape = 21, size = 11, colour = "white", stroke = 0.5,
             show.legend = TRUE) +
  scale_fill_manual(values = role_colours,
                    labels = c(fired = "fired (k₀)", neighbor = "neighbor", other = "other"),
                    name = NULL) +
  # Node labels: node number
  geom_text(data = nodes_g,
            aes(x, y, label = node),
            colour = "white", size = 4.2, fontface = "bold") +
  # Allowance labels below each node (use plain ASCII gamma label)
  geom_text(data = nodes_g,
            aes(x, y - 0.38, label = paste0("γ=", gamma)),
            colour = "grey20", size = 3.2) +
  facet_wrap(~ panel) +
  coord_fixed(ratio = 1, xlim = c(-0.7, 5.0), ylim = c(-0.75, 3.75)) +
  theme_void(base_size = 13) +
  theme(
    strip.text       = element_text(face = "bold", size = 11, margin = margin(b = 6)),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    plot.margin      = margin(6, 6, 6, 6)
  )

ggsave("figures/graph_spending.png", p0c, width = 10, height = 4.8, dpi = 220)

# ======================================================================
# 1. TSM log-wealth paths under no change vs. after a change
#    -- uses GaussianModel + TSM + compute_increments + increments_to_tsm
# ======================================================================
make_tsm_paths <- function(n_paths = 6, N = 240, nu = 120,
                           mean_post = 0.4, seed = 42) {
  set.seed(seed)
  model <- GaussianModel(mean_pre = 0, sd_pre = 1,
                         mean_post = mean_post, sd_post = 1)
  tsm <- TSM(model)
  dgp <- DGP(model, nu = nu)
  paths <- vector("list", n_paths)
  for (i in seq_len(n_paths)) {
    x <- generate_stream(dgp, N = N)
    z <- compute_tsm(tsm, x, log = TRUE)
    paths[[i]] <- data.frame(t = seq_len(N), logM = z, path = factor(i))
  }
  list(df = do.call(rbind, paths), nu = nu)
}

tp <- make_tsm_paths()
p1 <- ggplot(tp$df, aes(t, logM, group = path, colour = path)) +
  geom_line(alpha = 0.85, size = 0.6) +
  geom_vline(xintercept = tp$nu, linetype = "dashed", colour = "grey30") +
  geom_hline(yintercept = 0, colour = "grey60") +
  annotate("text", x = tp$nu + 2, y = max(tp$df$logM) * 0.9,
           label = "true change-point", hjust = 0, size = 4) +
  labs(title = "Test supermartingale: log-wealth over time",
       subtitle = "Bounded under P; grows after the change",
       x = "t", y = expression(log~M[t])) +
  theme_talk +
  theme(legend.position = "none")
ggsave("figures/tsm_growth.png", p1, width = 8, height = 4.2, dpi = 220)

# ======================================================================
# 2. Shiryaev--Roberts detector with ARL control at two levels
#    -- uses ShiryaevRobertsDetector + run_detector (mirrors workflow vignette)
# ======================================================================
make_sr_panel <- function(N = 240, nu = 120, mean_post = 0.4,
                          alphas = c(0.05, 0.001), seed = 7) {
  set.seed(seed)
  model <- GaussianModel(mean_pre = 0, sd_pre = 1,
                         mean_post = mean_post, sd_post = 1)
  tsm <- TSM(model)
  dgp <- DGP(model, nu = nu)
  x   <- generate_stream(dgp, N = N)
  inc <- compute_increments(tsm, x, log = TRUE)

  rows <- lapply(alphas, function(a) {
    sr <- ShiryaevRobertsDetector(alpha = a, criterion = "ARL")
    out <- run_detector(sr, inc, log = TRUE)
    data.frame(
      t = seq_len(N),
      logR = out$statistic,
      threshold = log(sr@threshold),
      alarm = if (is.finite(out$stopping_time)) out$stopping_time else NA_integer_,
      label = sprintf("ARL = %g", 1 / a)
    )
  })
  list(df = do.call(rbind, rows), nu = nu)
}

sp <- make_sr_panel()
p2 <- ggplot(sp$df, aes(t, logR)) +
  geom_line(colour = "#1f77b4", size = 0.8) +
  geom_hline(aes(yintercept = threshold), linetype = "dotted", colour = "red") +
  geom_vline(xintercept = sp$nu, linetype = "dashed", colour = "grey30") +
  geom_vline(aes(xintercept = alarm), colour = "forestgreen", size = 0.7, na.rm = TRUE) +
  facet_wrap(~ label, ncol = 2) +
  labs(title = "Shiryaev-Roberts detector at two ARL levels",
       subtitle = "Dashed: true change-point. Dotted: log threshold. Green: alarm.",
       x = "t", y = expression(log~R[t])) +
  theme_talk
ggsave("figures/sr_procedure.png", p2, width = 9, height = 4.2, dpi = 220)

# ======================================================================
# 3. Spending schedules: shapes + the corresponding PFA-controlled detectors
#    -- uses ShiryaevRobertsDetector(criterion = "PFA", spending = ...)
# ======================================================================
make_spending_panel <- function(N = 240, nu = 120, mean_post = 0.4,
                                alpha = 0.05, seed = 99) {
  set.seed(seed)
  model <- GaussianModel(mean_pre = 0, sd_pre = 1,
                         mean_post = mean_post, sd_post = 1)
  tsm <- TSM(model)
  dgp <- DGP(model, nu = nu)
  x   <- generate_stream(dgp, N = N)
  inc <- compute_increments(tsm, x, log = TRUE)

  schedules <- list(
    "geometric (mean 100)" = {s <- dgeom(0:(N - 1), prob = 1/100); s/sum(s)},
    "Poisson (mean = nu)"   = {s <- dpois(1:N, lambda = nu);       s/sum(s)},
    "heavy-tailed (1/t^1.2)" = {s <- 1/(1:N)^1.2;                  s/sum(s)}
  )
  sched_df <- do.call(rbind, lapply(names(schedules), function(nm) {
    data.frame(t = seq_len(N), pi = schedules[[nm]], kind = nm)
  }))
  sched_df$kind <- factor(sched_df$kind, levels = names(schedules))

  paths <- do.call(rbind, lapply(names(schedules), function(nm) {
    sr <- ShiryaevRobertsDetector(alpha = alpha, criterion = "PFA",
                                  spending = schedules[[nm]])
    out <- run_detector(sr, inc, log = TRUE)
    data.frame(
      t = seq_len(N),
      logR = out$statistic,
      threshold = log(sr@threshold),
      alarm = if (is.finite(out$stopping_time)) out$stopping_time else NA_integer_,
      kind = nm
    )
  }))
  paths$kind <- factor(paths$kind, levels = names(schedules))
  list(sched_df = sched_df, paths = paths, nu = nu)
}

sp3 <- make_spending_panel()
p3_top <- ggplot(sp3$sched_df, aes(t, pi, colour = kind)) +
  geom_line(size = 0.7) +
  labs(x = "t", y = expression(pi[t]), colour = NULL,
       title = "Three spending schedules") +
  theme_talk

p3_bot <- ggplot(sp3$paths, aes(t, logR, colour = kind)) +
  geom_line(size = 0.6) +
  geom_hline(aes(yintercept = threshold), linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = sp3$nu, linetype = "dashed", colour = "grey30") +
  geom_vline(aes(xintercept = alarm, colour = kind),
             linetype = "solid", size = 0.6, na.rm = TRUE) +
  facet_wrap(~ kind, nrow = 1) +
  labs(x = "t", y = expression(log~R[t]^pi),
       title = "PFA-controlled S-R statistic under each schedule",
       colour = NULL) +
  theme_talk +
  theme(legend.position = "none")

# stacked into one PNG via gridExtra-style layout, but stay dependency-light
ggsave("figures/spending_schedule.png", p3_top, width = 8.4, height = 3.2, dpi = 220)
ggsave("figures/spending_schedule_detectors.png", p3_bot, width = 9, height = 3.4, dpi = 220)

# ======================================================================
# 4. Sparse vs. dense change -- multivariate streams + helper TSMs
#    -- uses default_multivariate_gaussian_dgp + MultivariateGaussianModel
# ======================================================================
make_sparse_vs_dense <- function(N = 200, nu = 100, K = 6, Delta = 1.5,
                                 seed = 11) {
  set.seed(seed)
  mu_pre <- rep(0, K)
  mu_sparse <- mu_pre; mu_sparse[1] <- Delta
  mu_dense  <- mu_pre + Delta / sqrt(K)

  mk <- function(mu_post, label) {
    dgp <- DGP(
      MultivariateGaussianModel(
        mu_pre = mu_pre, Sigma_pre = diag(K),
        mu_post = mu_post, Sigma_post = diag(K)
      ),
      nu = nu
    )
    x <- generate_stream(dgp, N = N)
    list(x = x, label = label)
  }
  list(sparse = mk(mu_sparse, "sparse change (one stream)"),
       dense  = mk(mu_dense,  "dense change (all streams)"),
       nu = nu, K = K)
}

sd_data <- make_sparse_vs_dense()
to_df <- function(xobj) {
  x <- xobj$x; K <- ncol(x); N <- nrow(x)
  data.frame(
    t = rep(seq_len(N), K),
    val = as.vector(x),
    stream = factor(rep(seq_len(K), each = N)),
    kind = xobj$label
  )
}
df_sd <- rbind(to_df(sd_data$sparse), to_df(sd_data$dense))
df_sd$kind <- factor(df_sd$kind,
                     levels = c("sparse change (one stream)",
                                "dense change (all streams)"))

p4 <- ggplot(df_sd, aes(t, val, colour = stream)) +
  geom_line(alpha = 0.7, size = 0.4) +
  geom_vline(xintercept = sd_data$nu, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ kind, ncol = 2) +
  labs(title = "Sparse vs. dense post-change mean shift",
       x = "t", y = expression(X[tk])) +
  theme_talk +
  theme(legend.position = "none")
ggsave("figures/sparse_vs_dense.png", p4, width = 9, height = 3.6, dpi = 220)

# ======================================================================
# 5. Detector log-paths for the three combiners (Product / Average / UP)
#    under the same sparse vs. dense streams
#    -- uses combine_streams + ShiryaevRobertsDetector
# ======================================================================
combiner_panel <- function(sd_obj, alpha = 0.001) {
  K <- sd_obj$K
  marg_models <- lapply(seq_len(K), function(k) {
    GaussianModel(mean_pre = 0, sd_pre = 1,
                  mean_post = if (k == 1) 1.5 else 0.6, sd_post = 1)
  })
  marg_tsms <- lapply(marg_models, TSM)

  build_paths <- function(xobj) {
    x <- xobj$x
    inc_mat <- sapply(seq_len(K), function(k)
      compute_increments(marg_tsms[[k]], x[, k], log = TRUE))

    combiners <- list(
      "product"             = ProductCombiner(),
      "average"             = AverageCombiner(),
      "universal portfolio" = UniversalPortfolioCombiner()
    )

    rows <- lapply(names(combiners), function(nm) {
      inc <- combine_streams(combiners[[nm]], inc_mat, log = TRUE)
      sr <- ShiryaevRobertsDetector(alpha = alpha, criterion = "ARL")
      out <- run_detector(sr, inc, log = TRUE)
      data.frame(
        t = seq_along(out$statistic),
        logR = out$statistic,
        threshold = log(sr@threshold),
        combiner = nm,
        kind = xobj$label
      )
    })
    do.call(rbind, rows)
  }

  df <- rbind(build_paths(sd_obj$sparse), build_paths(sd_obj$dense))
  df$combiner <- factor(df$combiner,
                        levels = c("product", "average", "universal portfolio"))
  df
}

cp_df <- combiner_panel(sd_data)
p5 <- ggplot(cp_df, aes(t, logR, colour = combiner)) +
  geom_line(size = 0.7) +
  geom_hline(aes(yintercept = threshold), linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = sd_data$nu, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ kind, ncol = 2) +
  labs(title = "S-R log-statistic under three combiners",
       subtitle = "K = 6 streams; marginal Gaussian TSMs",
       x = "t", y = expression(log~R[t]),
       colour = NULL) +
  theme_talk
ggsave("figures/combiners_paths.png", p5, width = 9.5, height = 4.0, dpi = 220)

# ======================================================================
# 6. Real wastewater motivation plot from the package's CSV
#    -- a few high-volume sites, COVID and RSV trajectories
# ======================================================================
ww_csv <- file.path(PKG_DIR, "Data", "wastewater_aggregate.csv")
if (file.exists(ww_csv)) {
  ww <- suppressMessages(read_csv(ww_csv, show_col_types = FALSE))
  ww <- ww %>%
    filter(site != "", measureid %in% c("covN2", "rsv")) %>%
    mutate(weekstart = as.Date(weekstart))

  # Pick a handful of sites with the most COVID observations
  top_sites <- ww %>%
    filter(measureid == "covN2") %>%
    count(site, sort = TRUE) %>%
    slice_head(n = 4) %>%
    pull(site)

  ww_sub <- ww %>%
    filter(site %in% top_sites) %>%
    mutate(measureid = recode(measureid, covN2 = "COVID (N2)", rsv = "RSV"))

  p6 <- ggplot(ww_sub, aes(weekstart, w_avg, colour = measureid)) +
    geom_line(size = 0.5, alpha = 0.9) +
    facet_wrap(~ site, scales = "free_y", ncol = 2) +
    labs(title = "Wastewater viral signals (selected NWMP sites)",
         x = NULL, y = "weekly average (copies / mL)",
         colour = NULL) +
    theme_talk
  ggsave("figures/wastewater_motivation.png", p6,
         width = 9, height = 5, dpi = 220)
} else {
  message("Could not find ", ww_csv, " -- leaving wastewater motivation as placeholder.")
}

# ======================================================================
# 7. SMALL PILOT simulation for the regret-by-combiner slide.
#    n_rep is intentionally small. Scale up before the talk by setting
#    BIG_SIM <- TRUE and increasing the magnitudes grid / n_rep.
# ======================================================================
BIG_SIM <- FALSE   # <-- flip to TRUE for the conference-quality version
n_rep   <- if (BIG_SIM) 500 else 60
mags    <- if (BIG_SIM) seq(0.2, 2.0, length.out = 10) else seq(0.4, 1.8, length.out = 5)

run_regret_pilot <- function(K = 6, N = 300, nu = 100, alpha = 0.001,
                             sparse = TRUE, n_rep = 60,
                             mags = seq(0.4, 1.8, length.out = 5),
                             seed = 1) {
  set.seed(seed)
  out <- list()
  combiners <- list(product = ProductCombiner(),
                    average = AverageCombiner(),
                    `universal portfolio` = UniversalPortfolioCombiner())
  marg_tsms <- lapply(seq_len(K), function(k)
    TSM(GaussianModel(mean_pre = 0, sd_pre = 1,
                      mean_post = 0.6, sd_post = 1)))

  for (m in mags) {
    mu_post <- rep(0, K)
    if (sparse) mu_post[1] <- m else mu_post[] <- m / sqrt(K)
    dgp <- DGP(
      MultivariateGaussianModel(
        mu_pre = rep(0, K), Sigma_pre = diag(K),
        mu_post = mu_post, Sigma_post = diag(K)
      ),
      nu = nu)
    for (nm in names(combiners)) {
      delays <- replicate(n_rep, {
        x <- generate_stream(dgp, N = N)
        inc_mat <- sapply(seq_len(K), function(k)
          compute_increments(marg_tsms[[k]], x[, k], log = TRUE))
        inc <- combine_streams(combiners[[nm]], inc_mat, log = TRUE)
        sr <- ShiryaevRobertsDetector(alpha = alpha, criterion = "ARL")
        res <- run_detector(sr, inc, log = TRUE)
        if (is.finite(res$stopping_time) && res$stopping_time > nu)
          res$stopping_time - nu else NA_real_
      })
      out[[length(out) + 1]] <- data.frame(
        magnitude = m, combiner = nm,
        cad = mean(delays, na.rm = TRUE),
        scenario = if (sparse) "sparse" else "dense"
      )
    }
  }
  do.call(rbind, out)
}

sim_df <- rbind(
  run_regret_pilot(sparse = TRUE,  n_rep = n_rep, mags = mags, seed = 1),
  run_regret_pilot(sparse = FALSE, n_rep = n_rep, mags = mags, seed = 2)
)
sim_df$combiner <- factor(sim_df$combiner,
                          levels = c("product", "average", "universal portfolio"))

p_sim <- ggplot(sim_df, aes(magnitude, cad, colour = combiner, shape = combiner)) +
  geom_line(size = 0.7) + geom_point(size = 2.2) +
  facet_wrap(~ scenario, ncol = 2) +
  labs(title = "Conditional average delay vs. change magnitude",
       subtitle = sprintf("Pilot: K=6, n_rep=%d, alpha=0.001 (scale up via BIG_SIM)", n_rep),
       x = expression("change magnitude  "*"||"*mu[1]*"||"[2]),
       y = "CAD",
       colour = NULL, shape = NULL) +
  theme_talk
ggsave("figures/sim_regret_placeholder.png", p_sim,
       width = 9, height = 4.2, dpi = 220)

# ======================================================================
# 8. Section 6 simulation results: CAD vs change magnitude by detector
#    -- reads simulation_output/var_sim_results.rds produced by
#       inst/scripts/var_simulations.R (quick or full mode)
# ======================================================================
var_sim_rds <- file.path(PKG_DIR, "simulation_output", "var_sim_results.rds")
if (file.exists(var_sim_rds)) {
  var_res <- readRDS(var_sim_rds)

  # ---- 8a. K=1: oracle vs misspec vs bounded (ARL criterion, nu=100) ----
  res1 <- subset(var_res, K == 1 & criterion == "ARL" & nu == 100)
  res1$Detector <- factor(
    paste0(res1$detector, ifelse(res1$detector == "bounded",
                                 paste0(" (", res1$combine, ")"), "")),
    levels = c("oracle", "misspec", "bounded (average)")
  )
  res1$label <- factor(
    res1$Detector,
    labels = c("oracle (true params)",
               "misspecified (mu1/2)",
               "bounded (AGRAPA)")
  )

  p8a <- ggplot(res1, aes(delta, CAD, colour = label, shape = label)) +
    geom_line(size = 0.8) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c("steelblue", "tomato", "forestgreen")) +
    labs(
      title    = "Conditional average delay vs. change magnitude  (K = 1)",
      subtitle = "ARL criterion, nu = 100, p = 0, sigma = 0.06",
      x        = expression(delta ~ "(change in mean)"),
      y        = "CAD",
      colour   = NULL, shape = NULL
    ) +
    theme_talk
  ggsave("figures/sim_cad_k1.png", p8a, width = 7.5, height = 4.2, dpi = 220)

  # ---- 8b. K=2: oracle vs all bounded combiners (ARL criterion, nu=100) ----
  res2 <- subset(var_res, K == 2 & criterion == "ARL" & nu == 100)
  res2$label <- with(res2, ifelse(
    detector %in% c("oracle", "misspec"),
    detector,
    paste0("bounded (", combine, ")")
  ))
  res2$label <- factor(res2$label,
    levels = c("oracle", "misspec",
               "bounded (average)", "bounded (product)", "bounded (up)"),
    labels = c("oracle", "misspecified",
               "bounded: average", "bounded: product", "bounded: UP")
  )

  p8b <- ggplot(res2, aes(delta, CAD, colour = label, shape = label)) +
    geom_line(size = 0.8) +
    geom_point(size = 2.2) +
    scale_colour_manual(
      values = c("steelblue", "tomato",
                 "forestgreen", "darkorange", "purple")
    ) +
    labs(
      title    = "Conditional average delay vs. change magnitude  (K = 2)",
      subtitle = "ARL criterion, nu = 100, p = 0, sigma = 0.06, independent streams",
      x        = expression(delta ~ "(change in mean, equal across streams)"),
      y        = "CAD",
      colour   = NULL, shape   = NULL
    ) +
    theme_talk
  ggsave("figures/sim_cad_k2.png", p8b, width = 7.5, height = 4.2, dpi = 220)

  # ---- 8c. Power curves (K=1 and K=2 side-by-side) ----
  res_pwr <- subset(var_res, criterion == "ARL" & nu == 100 &
                      detector %in% c("oracle", "bounded") &
                      (K == 1 | (K == 2 & combine == "average")))
  res_pwr$Detector <- paste0(res_pwr$detector, " (K=", res_pwr$K, ")")

  p8c <- ggplot(res_pwr, aes(delta, power, colour = Detector, shape = Detector)) +
    geom_line(size = 0.8) +
    geom_point(size = 2.2) +
    labs(
      title    = "Detection power vs. change magnitude",
      subtitle = "Fraction of runs that stopped after the change-point",
      x        = expression(delta),
      y        = "Power",
      colour   = NULL, shape = NULL
    ) +
    theme_talk
  ggsave("figures/sim_power.png", p8c, width = 7.5, height = 4.2, dpi = 220)

  message("Wrote simulation result figures (8a, 8b, 8c).")
} else {
  message("var_sim_results.rds not found — skipping Section 6 figures.",
          "\nRun inst/scripts/var_simulations.R first.")
}

# ======================================================================
# 9. Placeholders for results-only figures that still need real data
# ======================================================================
placeholder <- function(filename, title,
                        subtitle = "Replace with real result before the talk",
                        w = 8, h = 4.5) {
  p <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "grey92", colour = "grey60") +
    annotate("text", x = 0.5, y = 0.6, label = title,
             size = 7, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.45, label = subtitle,
             size = 4.5, colour = "grey40") +
    annotate("text", x = 0.5, y = 0.25, label = "PLACEHOLDER",
             size = 9, colour = "red", alpha = 0.35) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void()
  ggsave(filename, p, width = w, height = h, dpi = 220)
}

placeholder("figures/arl_vs_pfa_placeholder.png",
            "ARL vs. PFA: regret as a function of alpha",
            "Two CAD curves at matched alpha; ARL-SR vs. PFA-SR with concentrated schedule")

placeholder("figures/sim_arl_vs_pfa_localization_placeholder.png",
            "Cost of stronger control / localization",
            "Side-by-side regret panels (ARL vs PFA; pooled vs localized)",
            w = 9, h = 4.2)

placeholder("figures/ww_detector_placeholder.png",
            "Wastewater detector on 2025-26 season",
            "Residuals (top) and log Rt with threshold (bottom)",
            w = 9, h = 4.6)

placeholder("figures/ww_localization_placeholder.png",
            "Adjacency-aware spending allowance",
            "Cartoon of WW sites with reallocation after firing",
            w = 7, h = 4.0)

message("Wrote figures/*.png")

# ======================================================================
# 10. Joint + localized simulation results
#     Reads simulations/output/{joint,localized}_sim_results.csv
#     Saves to figures/simulations/
# ======================================================================
sim_dir_cad <- "figures/simulations/cad"
sim_dir_arl <- "figures/simulations/arl"
if (!dir.exists(sim_dir_cad)) dir.create(sim_dir_cad, recursive = TRUE)
if (!dir.exists(sim_dir_arl)) dir.create(sim_dir_arl, recursive = TRUE)

joint_csv <- file.path(PKG_DIR, "simulations", "output", "joint_sim_results.csv")
local_csv <- file.path(PKG_DIR, "simulations", "output", "localized_sim_results.csv")

if (!file.exists(joint_csv) || !file.exists(local_csv)) {
  message("Simulation CSVs not found — skipping section 10.")
} else {

joint_raw <- read.csv(joint_csv,  stringsAsFactors = FALSE)
local_raw <- read.csv(local_csv,  stringsAsFactors = FALSE)

# ── Shared colour / shape palettes ────────────────────────────────────────
det_pal <- c(
  "oracle"           = "#1f77b4",
  "misspecified"     = "#d62728",
  "bounded: average" = "#2ca02c",
  "bounded: product" = "#ff7f0e"
)
det_shapes <- c("oracle" = 16, "misspecified" = 17,
                "bounded: average" = 15, "bounded: product" = 18)

dep_pal <- c(
  "indep + average" = "#2ca02c",
  "indep + product" = "#ff7f0e",
  "corr  + average" = "#9467bd"
)

p_pal    <- c("AR(0)" = "#2c7bb6", "AR(1)" = "#d7191c")
p_ltypes <- c("AR(0)" = "solid",   "AR(1)" = "dashed")

# ── Tidy joint data ───────────────────────────────────────────────────────
joint <- joint_raw %>%
  filter(delta > 0) %>%
  mutate(
    detector_label = case_when(
      detector == "oracle"  ~ "oracle",
      detector == "misspec" ~ "misspecified",
      detector == "bounded" & combine == "average" ~ "bounded: average",
      detector == "bounded" & combine == "product" ~ "bounded: product"
    ),
    detector_label = factor(detector_label, levels = names(det_pal)),
    K_label     = factor(paste0("K = ", K),   levels = c("K = 2", "K = 10")),
    sparse_label = factor(ifelse(sparse, "sparse change", "dense change"),
                          levels = c("sparse change", "dense change")),
    p_label     = factor(paste0("AR(", p, ")"), levels = c("AR(0)", "AR(1)"))
  )

# ── Tidy localized data ───────────────────────────────────────────────────
local <- local_raw %>%
  filter(delta > 0) %>%
  mutate(
    K_label     = factor(paste0("K = ", K),   levels = c("K = 2", "K = 10")),
    sparse_label = factor(ifelse(sparse, "sparse change", "dense change"),
                          levels = c("sparse change", "dense change")),
    p_label     = factor(paste0("AR(", p, ")"), levels = c("AR(0)", "AR(1)")),
    dep_label   = factor(ifelse(independent, "independent", "correlated"),
                         levels = c("independent", "correlated")),
    condition   = factor(
      paste0("AR(", p, "), ", ifelse(independent, "independent", "correlated")),
      levels = c("AR(0), independent", "AR(0), correlated",
                 "AR(1), independent", "AR(1), correlated"))
  )

# ── Shared axis scales (log-log) ──────────────────────────────────────────
x_log <- list(
  scale_x_log10(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2),
                labels = c("0.01", "0.02", "0.05", "0.1", "0.2")),
  xlab(expression(delta ~ "(shift magnitude)"))
)
y_log <- list(
  scale_y_log10(),
  ylab("CAD (log scale)")
)

# ── Plot J1: Oracle vs misspecified vs bounded ─────────────────────────────
# Filter: IID (p=0), independent streams, late change (nu=1000)
j1 <- joint %>% filter(p == 0, independent == TRUE, nu == 1000)

pJ1 <- ggplot(j1, aes(delta, CAD,
                       colour = detector_label, shape = detector_label)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.0, na.rm = TRUE) +
  scale_colour_manual(values = det_pal,    name = NULL) +
  scale_shape_manual( values = det_shapes, name = NULL) +
  x_log + y_log +
  facet_grid(K_label ~ sparse_label) +
  labs(
    title    = "Oracle vs. misspecified vs. bounded detectors",
    subtitle = expression(
      "IID streams (" * italic(p) * " = 0), independent, " *
      nu * " = 1000, " * alpha * " = 10"^{-4})
  ) +
  theme_talk +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(sim_dir_cad, "joint_oracle_vs_bounded.png"),
       pJ1, width = 8, height = 6.5, dpi = 220)

# ── Plot J2: Effect of AR(1) (bounded + average combiner) ─────────────────
# Shows: does serial correlation within streams hurt detection?
j2 <- joint %>%
  filter(independent == TRUE, nu == 1000,
         detector == "bounded", combine == "average")

pJ2 <- ggplot(j2, aes(delta, CAD,
                       colour = p_label, linetype = p_label)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  scale_colour_manual(values = p_pal,    name = "AR order") +
  scale_linetype_manual(values = p_ltypes, name = "AR order") +
  x_log + y_log +
  facet_grid(K_label ~ sparse_label) +
  labs(
    title    = "Effect of within-stream serial correlation on detection delay",
    subtitle = expression(
      "Bounded + average combiner, independent streams, " *
      nu * " = 1000, " * alpha * " = 10"^{-4})
  ) +
  theme_talk +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(sim_dir_cad, "joint_effect_of_ar.png"),
       pJ2, width = 8, height = 6.5, dpi = 220)

# ── Plot J3: Effect of cross-stream dependence ─────────────────────────────
# Compare: independent+average, independent+product, correlated+average.
# Product is valid only under independence; average is always valid.
j3 <- joint %>%
  filter(p == 0, nu == 1000, detector == "bounded") %>%
  mutate(dep_combine = factor(
    paste0(ifelse(independent, "indep", "corr"), " + ", combine),
    levels = c("indep + average", "indep + product", "corr  + average")
  )) %>%
  filter(!is.na(dep_combine))   # drop corr+product (not run)

pJ3 <- ggplot(j3, aes(delta, CAD,
                       colour = dep_combine, shape = dep_combine)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.0, na.rm = TRUE) +
  scale_colour_manual(values = dep_pal, name = NULL) +
  scale_shape_manual( values = c(15, 18, 17), name = NULL) +
  x_log + y_log +
  facet_grid(K_label ~ sparse_label) +
  labs(
    title    = "Effect of cross-stream dependence on detection delay",
    subtitle = expression(
      "Bounded combiner, IID observations (" * italic(p) * " = 0), " *
      nu * " = 1000, " * alpha * " = 10"^{-4})
  ) +
  theme_talk +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(sim_dir_cad, "joint_effect_of_dependence.png"),
       pJ3, width = 8, height = 6.5, dpi = 220)

# ── Plot J4: Early vs. late change-point (nu = 10 vs nu = 1000) ───────────
# Fixed: bounded+average, IID, independent, sparse (most interesting)
j4 <- joint %>%
  filter(p == 0, independent == TRUE,
         detector == "bounded", combine == "average") %>%
  mutate(nu_label = factor(paste0("nu == ", nu),
                           levels = c("nu == 10", "nu == 1000"),
                           labels = c(expression(nu == 10),
                                      expression(nu == 1000))))

pJ4 <- ggplot(j4, aes(delta, CAD,
                       colour = factor(nu), linetype = factor(nu))) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  scale_colour_manual(values = c("10" = "#e08214", "1000" = "#542788"),
                      name = expression(nu)) +
  scale_linetype_manual(values = c("10" = "dashed", "1000" = "solid"),
                        name = expression(nu)) +
  x_log + y_log +
  facet_grid(K_label ~ sparse_label) +
  labs(
    title    = "Effect of change-point location on detection delay",
    subtitle = expression(
      "Bounded + average combiner, IID, independent, " * alpha * " = 10"^{-4})
  ) +
  theme_talk +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(sim_dir_cad, "joint_effect_of_nu.png"),
       pJ4, width = 8, height = 6.5, dpi = 220)

# ── Plot L1: Localized detector — global CAD ──────────────────────────────
# All four combinations of p x independence; nu = 100 (only value available)
cond_pal <- c(
  "AR(0), independent" = "#1a9641",
  "AR(0), correlated"  = "#a6d96a",
  "AR(1), independent" = "#d7191c",
  "AR(1), correlated"  = "#fdae61"
)
cond_shapes <- c("AR(0), independent" = 16, "AR(0), correlated"  = 17,
                 "AR(1), independent" = 15, "AR(1), correlated"  = 18)

pL1 <- ggplot(local, aes(delta, CAD,
                          colour = condition, shape = condition)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  scale_colour_manual(values = cond_pal,    name = NULL) +
  scale_shape_manual( values = cond_shapes, name = NULL) +
  x_log + y_log +
  facet_grid(K_label ~ sparse_label) +
  labs(
    title    = "Localized detector: global conditional average delay",
    subtitle = expression(
      "Bonferroni spending across K streams, " *
      nu * " = 100, " * alpha * " = 10"^{-4})
  ) +
  theme_talk +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(sim_dir_cad, "localized_global_cad.png"),
       pL1, width = 8, height = 6.5, dpi = 220)

# ── Plot L2: Localized — per-stream CAD for sparse change (K=10) ──────────
# Highlight stream 1 (the affected stream) vs. unaffected streams 2-10.
# Reveals whether the localized detector correctly pinpoints the changed stream.
stream_cols <- paste0("stream_", 1:10, "_CAD")
l2_wide <- local %>%
  filter(K == 10, sparse == TRUE, p == 0, independent == TRUE)

l2_long <- l2_wide %>%
  select(delta, all_of(stream_cols)) %>%
  pivot_longer(all_of(stream_cols),
               names_to  = "stream",
               values_to = "stream_CAD") %>%
  mutate(
    stream_id  = as.integer(sub("stream_(\\d+)_CAD", "\\1", stream)),
    stream_grp = factor(ifelse(stream_id == 1L, "stream 1 (changed)",
                               "streams 2-10 (unchanged)"),
                        levels = c("stream 1 (changed)", "streams 2-10 (unchanged)"))
  )

pL2 <- ggplot(l2_long, aes(delta, stream_CAD, group = stream,
                             colour = stream_grp, alpha = stream_grp)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  scale_colour_manual(values = c("stream 1 (changed)"    = "#d62728",
                                 "streams 2-10 (unchanged)" = "#7f7f7f"),
                      name = NULL) +
  scale_alpha_manual(values = c("stream 1 (changed)" = 1.0,
                                "streams 2-10 (unchanged)" = 0.45),
                     name = NULL) +
  x_log + y_log +
  labs(
    title    = "Per-stream CAD: sparse change in K = 10 streams",
    subtitle = expression(
      "Only stream 1 has a shift; p = 0, independent, " *
      nu * " = 100, " * alpha * " = 10"^{-4})
  ) +
  theme_talk +
  theme(legend.position = "bottom")

ggsave(file.path(sim_dir_cad, "localized_stream_cad_sparse.png"),
       pL2, width = 7, height = 4.5, dpi = 220)

message("Wrote figures/simulations/cad/  (", 6L, " joint+localized CAD figures)")
} # end if (file.exists(...))

# ======================================================================
# 11. Bounded AGRAPA simulation results
#     Reads simulations/output/bounded_sim_results.csv  (Section A)
#           simulations/output/bounded_w_results.csv    (Section B)
#     Saves to figures/simulations/
# ======================================================================
bounded_csv <- file.path(PKG_DIR, "simulations", "output", "bounded_sim_results.csv")
w_csv       <- file.path(PKG_DIR, "simulations", "output", "bounded_w_results.csv")

if (!file.exists(bounded_csv)) {
  message("bounded_sim_results.csv not found — skipping section 11.")
} else {

bounded <- read.csv(bounded_csv, stringsAsFactors = FALSE)

strat_levels <- c("Full AGRAPA", "Window AGRAPA (W=30)", "EWMA (rho=0.10)",
                   "Per-clock W=30", "Oracle Kelly")
strat_labels <- c("Full AGRAPA", "Window AGRAPA (W=30)", "EWMA (rho=0.10)",
                   "Per-clock (W=30)", "Oracle Kelly")
strat_pal    <- c(
  "Full AGRAPA"           = "#1f77b4",
  "Window AGRAPA (W=30)"  = "#ff7f0e",
  "EWMA (rho=0.10)"       = "#2ca02c",
  "Per-clock W=30"        = "#d62728",
  "Oracle Kelly"          = "#9467bd"
)
strat_shapes <- c(
  "Full AGRAPA"           = 16,
  "Window AGRAPA (W=30)"  = 17,
  "EWMA (rho=0.10)"       = 15,
  "Per-clock W=30"        = 18,
  "Oracle Kelly"          = 8
)

bounded$strategy <- factor(bounded$strategy, levels = strat_levels,
                            labels = strat_labels)

# relabel after factor creation
names(strat_pal)    <- strat_labels
names(strat_shapes) <- strat_labels

# ── Panel A: CAD vs delta ────────────────────────────────────────────────
pB_cad <- ggplot(bounded, aes(delta, CAD, colour = strategy, shape = strategy)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.2, na.rm = TRUE) +
  scale_colour_manual(values = strat_pal,    name = NULL) +
  scale_shape_manual( values = strat_shapes, name = NULL) +
  labs(
    title    = "Conditional average delay — bounded [0,1] change detection",
    subtitle = expression(
      "Pre: Uniform(0,1); Post: Beta(4(0.5 + "*delta*"), 4(0.5 - "*delta*"))"~
      "  "*alpha*" = 0.001, "*nu*" = 50, N = 500"),
    x = expression(delta ~ "(mean shift)"),
    y = "CAD"
  ) +
  theme_talk

# ── Panel B: computation time vs delta ───────────────────────────────────
pB_time <- ggplot(bounded, aes(delta, comp_time_ms,
                                colour = strategy, shape = strategy)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.2, na.rm = TRUE) +
  scale_colour_manual(values = strat_pal,    name = NULL) +
  scale_shape_manual( values = strat_shapes, name = NULL) +
  labs(
    title = "Computation time per replication",
    x     = expression(delta ~ "(mean shift)"),
    y     = "ms / replication"
  ) +
  theme_talk

ggsave(file.path(sim_dir_cad, "bounded_cad.png"),
       pB_cad,  width = 8.5, height = 4.5, dpi = 220)
ggsave(file.path(sim_dir_cad, "bounded_comp_time.png"),
       pB_time, width = 8.5, height = 4.5, dpi = 220)

message("Wrote bounded strategy comparison figures (11a, 11b).")
} # end bounded_csv block

# ── Section B: W-sensitivity ──────────────────────────────────────────────
if (!file.exists(w_csv)) {
  message("bounded_w_results.csv not found — skipping section 11b.")
} else {

w_res <- read.csv(w_csv, stringsAsFactors = FALSE)
w_res$W             <- factor(w_res$W, levels = c(1, 5, 30, 100))
w_res$strategy_type <- factor(w_res$strategy_type,
                               levels = c("Window", "Per-clock"),
                               labels = c("Window AGRAPA", "Per-clock"))

w_pal    <- c("1" = "#1f77b4", "5" = "#ff7f0e", "30" = "#2ca02c", "100" = "#d62728")
w_shapes <- c("1" = 16,        "5" = 17,         "30" = 15,        "100" = 18)

pW_cad <- ggplot(w_res, aes(delta, CAD, colour = W, shape = W, group = W)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.0, na.rm = TRUE) +
  scale_colour_manual(values = w_pal,    name = "Window (W)") +
  scale_shape_manual( values = w_shapes, name = "Window (W)") +
  facet_wrap(~ strategy_type) +
  labs(
    title    = "CAD vs. delta by window size",
    subtitle = expression(
      "Pre: Uniform(0,1); Post: Beta(4(0.5 + "*delta*"), 4(0.5 - "*delta*"))"~
      "  "*alpha*" = 0.001"),
    x = expression(delta ~ "(mean shift)"),
    y = "CAD"
  ) +
  theme_talk +
  theme(strip.text = element_text(face = "bold"))

pW_time <- ggplot(w_res, aes(delta, comp_time_ms, colour = W, shape = W, group = W)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.0, na.rm = TRUE) +
  scale_colour_manual(values = w_pal,    name = "Window (W)") +
  scale_shape_manual( values = w_shapes, name = "Window (W)") +
  facet_wrap(~ strategy_type) +
  labs(
    title = "Computation time vs. delta by window size",
    x     = expression(delta ~ "(mean shift)"),
    y     = "ms / replication"
  ) +
  theme_talk +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(sim_dir_cad, "bounded_w_cad.png"),
       pW_cad,  width = 9, height = 4.5, dpi = 220)
ggsave(file.path(sim_dir_cad, "bounded_w_comp_time.png"),
       pW_time, width = 9, height = 4.5, dpi = 220)

message("Wrote W-sensitivity figures (11c, 11d).")
} # end w_csv block

# ── Section 11C: change-point location sensitivity ────────────────────────
nu_csv <- file.path(PKG_DIR, "simulations", "output", "bounded_nu_results.csv")
if (!file.exists(nu_csv)) {
  message("bounded_nu_results.csv not found — skipping section 11c.")
} else {

nu_res <- read.csv(nu_csv, stringsAsFactors = FALSE)

strat_levels_c <- c("Full AGRAPA", "Window AGRAPA (W=30)", "EWMA (rho=0.10)",
                    "Per-clock W=30", "Oracle Kelly")
strat_labels_c <- c("Full AGRAPA", "Window AGRAPA (W=30)", "EWMA (rho=0.10)",
                    "Per-clock (W=30)", "Oracle Kelly")
strat_pal_c    <- setNames(
  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"), strat_labels_c)
strat_shapes_c <- setNames(c(16, 17, 15, 18, 8), strat_labels_c)

nu_res$strategy <- factor(nu_res$strategy,
                           levels = strat_levels_c, labels = strat_labels_c)
nu_res$nu_label <- factor(paste0("nu == ", nu_res$nu),
                           levels = paste0("nu == ", c(10, 50, 200)))

pC_cad <- ggplot(nu_res, aes(delta, CAD, colour = strategy, shape = strategy)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.2, na.rm = TRUE) +
  scale_colour_manual(values = strat_pal_c,    name = NULL) +
  scale_shape_manual( values = strat_shapes_c, name = NULL) +
  facet_wrap(~ nu_label, labeller = label_parsed) +
  labs(
    title    = "Effect of change-point location on detection delay",
    subtitle = expression(
      "Full AGRAPA burns in pre-change bet as "*nu*" grows; windowed strategies adapt faster"),
    x = expression(delta ~ "(mean shift)"),
    y = "CAD"
  ) +
  theme_talk +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(sim_dir_cad, "bounded_nu_cad.png"),
       pC_cad, width = 10, height = 4.5, dpi = 220)

message("Wrote nu-sensitivity figure (11c).")
} # end nu_csv block
