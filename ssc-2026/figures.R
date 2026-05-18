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
# 1. TSM log-wealth paths under no change vs. after a change
#    -- uses GaussianModel + TSM + compute_increments + increments_to_tsm
# ======================================================================
make_tsm_paths <- function(n_paths = 6, N = 240, nu = 120,
                           mean_post = 0.4, seed = 42) {
  set.seed(seed)
  model <- GaussianModel(mean_pre = 0, sd_pre = 1,
                         mean_post = mean_post, sd_post = 1)
  tsm <- TSM(model)
  dgp <- DGP(
    generator = default_gaussian_dgp,
    pre_params = list(mean = 0, sd = 1),
    post_params = list(mean = mean_post, sd = 1),
    nu = nu
  )
  paths <- vector("list", n_paths)
  for (i in seq_len(n_paths)) {
    x <- generate_stream(dgp, N = N, K = 1)
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
  dgp <- DGP(generator = default_gaussian_dgp,
             pre_params = list(mean = 0, sd = 1),
             post_params = list(mean = mean_post, sd = 1),
             nu = nu)
  x <- generate_stream(dgp, N = N, K = 1)
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
  dgp <- DGP(generator = default_gaussian_dgp,
             pre_params = list(mean = 0, sd = 1),
             post_params = list(mean = mean_post, sd = 1),
             nu = nu)
  x <- generate_stream(dgp, N = N, K = 1)
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
      generator = default_multivariate_gaussian_dgp,
      pre_params = list(mu = mu_pre, Sigma = diag(K)),
      post_params = list(mu = mu_post, Sigma = diag(K)),
      nu = nu
    )
    x <- generate_stream(dgp, N = N, K = K)
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
    dgp <- DGP(generator = default_multivariate_gaussian_dgp,
               pre_params  = list(mu = rep(0, K), Sigma = diag(K)),
               post_params = list(mu = mu_post, Sigma = diag(K)),
               nu = nu)
    for (nm in names(combiners)) {
      delays <- replicate(n_rep, {
        x <- generate_stream(dgp, N = N, K = K)
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
# 8. Placeholders for results-only figures that still need real data
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
