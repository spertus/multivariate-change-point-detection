# wastewater_final_figures.R
# Manuscript-ready versions of the two final figures from
# vignettes/wastewater.Rmd's graph-aware localized detection pipeline:
#   1. Transformed signal, 2023-2025 (per-site log1p(covN2) normalized by each
#      site's own 2023 average)
#   2. Alarm plot: 2024-2025 normalized signal with detection dates under the
#      uniform vs. proximity allowance
#
# This script transcribes the vignette's pipeline exactly (same site
# selection, same log1p/mean-2023 normalization, same eta_k = mu_init_k = 1,
# same rho = 0.7, same proximity graph / xi = 200 / zeta = 0.9, same
# alpha = 1/52) so the two figures match the vignette's numbers precisely;
# it only changes presentation (font size, y-axis label, and -- for the
# alarm plot -- a 3-colour scheme so overlapping same-date alarms under both
# allowances are visible instead of one line hiding the other).
#
# Run from the ssc-2026/ folder:
#   Rscript wastewater_final_figures.R
#
# Output: figures/wastewater_transformed_signal.png  -- standalone panel (a), with title
#         figures/wastewater_alarm_plot.png          -- standalone panel (b), with title
#         figures/wastewater_combined.png            -- both panels stacked, (a)/(b)
#                                                        tagged, titles dropped in favour
#                                                        of the external figure caption

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(canadianmaps)
  library(maps)
  library(mgcv)
  library(patchwork)
})

PKG_DIR <- normalizePath(file.path(".."), mustWork = FALSE)
if (requireNamespace("devtools", quietly = TRUE)) {
  suppressMessages(devtools::load_all(PKG_DIR))
} else {
  library(multichangepoints)
}

FONT_SIZE <- 16
Y_LABEL   <- "Proportion of 2023 Average COVID-19 Signal"

# ── 1. Data + site selection (identical to vignette) ────────────────────────

ww <- suppressMessages(read_csv(file.path(PKG_DIR, "Data", "wastewater_aggregate.csv"),
                                show_col_types = FALSE)) %>%
  mutate(weekstart = as.Date(weekstart))

window_start <- as.Date("2023-01-01")
train_end    <- as.Date("2024-01-01")   # training = 2023 only; monitor 2024-2025
window_end   <- as.Date("2026-01-01")
spine_pipe   <- seq(window_start, window_end, by = "7 days")
n_expected   <- length(spine_pipe)

site_enrollment <- ww %>%
  filter(measureid == "covN2", site != "") %>%
  group_by(site) %>%
  summarize(first_obs = min(weekstart), last_obs = max(weekstart), .groups = "drop")

site_completeness <- ww %>%
  filter(measureid %in% c("covN2", "rsv"), site != "",
         weekstart >= window_start, weekstart <= window_end) %>%
  select(site, weekstart, measureid, w_avg) %>%
  pivot_wider(names_from = measureid, values_from = w_avg) %>%
  group_by(site) %>%
  summarize(
    n_covN2   = sum(!is.na(covN2)),
    n_rsv     = sum(!is.na(rsv)),
    pct_covN2 = round(n_covN2 / n_expected, 3),
    pct_rsv   = round(n_rsv   / n_expected, 3),
    .groups   = "drop"
  ) %>%
  left_join(site_enrollment, by = "site") %>%
  mutate(
    enrolled_by_start = first_obs <= window_start,
    through_end       = last_obs  >= window_end,
    meets_criteria    = enrolled_by_start & through_end &
                        pct_covN2 >= 0.80 & pct_rsv >= 0.80
  )

selected_sites <- site_completeness %>%
  filter(meets_criteria) %>%
  arrange(site) %>%
  pull(site)
K_sites <- length(selected_sites)

# ── 2. log1p transform, per-site 2023-average normalization ────────────────

impute_row <- function(row) {
  obs <- which(!is.na(row))
  if (length(obs) < 2L) return(row)
  first <- obs[1L]; last <- obs[length(obs)]
  row[first:last] <- approx(obs, row[obs], xout = first:last)$y
  row
}

pipe_long <- ww %>%
  filter(measureid == "covN2", site %in% selected_sites,
         weekstart >= window_start, weekstart <= window_end) %>%
  select(site, weekstart, covN2 = w_avg) %>%
  complete(site, weekstart = spine_pipe) %>%
  arrange(site, weekstart) %>%
  group_by(site) %>%
  mutate(covN2 = impute_row(covN2),
         log_covN2 = log1p(covN2)) %>%
  ungroup()

mean2023_k <- pipe_long %>%
  filter(weekstart < train_end) %>%
  group_by(site) %>%
  summarize(mean2023 = mean(log_covN2, na.rm = TRUE), .groups = "drop")

pipe_long <- pipe_long %>%
  left_join(mean2023_k, by = "site") %>%
  mutate(norm = log_covN2 / mean2023)

norm_matrix <- pipe_long %>%
  select(site, weekstart, norm) %>%
  pivot_wider(names_from = weekstart, values_from = norm) %>%
  column_to_rownames("site") %>%
  as.matrix()

train_idx <- which(spine_pipe <  train_end)
test_idx  <- which(spine_pipe >= train_end)

eta_k     <- setNames(rep(1, K_sites), rownames(norm_matrix))
mu_init_k <- setNames(rep(1, K_sites), rownames(norm_matrix))

test_mat   <- t(norm_matrix[, test_idx, drop = FALSE])
test_dates <- spine_pipe[test_idx]
N_test     <- nrow(test_mat)

inc_mat <- matrix(NA_real_, nrow = N_test, ncol = K_sites,
                  dimnames = list(NULL, rownames(norm_matrix)))
for (k in seq_len(K_sites)) {
  bm_k <- BoundedModel(eta = eta_k[k], bets = EWMABet(rho = 0.7, mu_init = mu_init_k[k]))
  inc_mat[, k] <- compute_increments(TSM(bm_k), test_mat[, k], log = TRUE)
}

# ── 3. Proximity graph (identical geocoding to the vignette) ───────────────

ca_cities <- world.cities[world.cities$country.etc == "Canada", ]

manual_coords <- tribble(
  ~site,                                        ~lat,    ~long,
  "Battery Point",                              46.14,  -60.19,
  "Battleford",                                 52.73, -108.30,
  "Birch Hills",                                52.98, -105.44,
  "City of Charlottetown and Town of Stratford",46.24,  -63.13,
  "Dominion-Bridgeport",                        46.13,  -60.00,
  "Miramichi",                                  47.02,  -65.50,
  "Pasqua",                                     50.58, -103.98,
  "Peel G.E. Booth",                            43.73,  -79.72,
  "Southey",                                    50.95, -104.49,
  "St. John's",                                 47.56,  -52.71,
  "St. Stephen",                                45.20,  -67.28,
  "Trenton",                                    45.61,  -62.63,
  "Vancouver Annacis Island",                   49.17, -122.92,
  "Vancouver Iona Island",                      49.22, -123.21,
  "Vancouver Lions Gate",                       49.35, -123.13,
  "Vancouver Lulu Island",                      49.14, -123.19,
  "Vancouver Northwest Langley",                49.17, -122.60,
  "Île-à-la-Crosse",                  55.45, -107.89
)

site_city <- ww %>%
  filter(site %in% selected_sites) %>%
  distinct(site, city, province) %>%
  mutate(city_clean = trimws(city))

site_geo_k <- site_city %>%
  left_join(ca_cities %>% select(name, lat, long, pop),
            by = c("city_clean" = "name"), relationship = "many-to-many") %>%
  group_by(site) %>%
  slice_max(pop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rows_update(manual_coords, by = "site", unmatched = "ignore") %>%
  filter(!is.na(lat)) %>%
  arrange(match(site, selected_sites))

stopifnot(nrow(site_geo_k) == K_sites, all(site_geo_k$site == rownames(norm_matrix)))

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad; dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

MIN_DIST_KM <- 5
W_prox <- matrix(0, K_sites, K_sites)
for (i in seq_len(K_sites)) for (j in seq_len(K_sites)) {
  if (i != j) W_prox[i, j] <- max(MIN_DIST_KM,
    haversine_km(site_geo_k$lat[i], site_geo_k$long[i],
                site_geo_k$lat[j], site_geo_k$long[j]))
}
prox_graph <- ProximityGraph(W_prox, names = site_geo_k$site)
xi <- 200

# ── 4. Detection: uniform vs. proximity, alpha = 1/52 ───────────────────────

alpha_global <- 1 / 52

ld_unif <- LocalizedDetector(K = K_sites, alpha = alpha_global, criterion = "ARL")
ld_prox <- LocalizedDetector(
  K = K_sites, alpha = alpha_global, criterion = "ARL",
  proximity_graph = prox_graph, kernel = exponential_kernel(xi = xi), zeta = 0.9
)

det_unif <- run_detector(ld_unif, evidence = inc_mat, log = TRUE)
det_prox <- run_detector(ld_prox, evidence = inc_mat, log = TRUE)

stop_to_date <- function(stops) {
  out <- rep("no alarm", length(stops))
  fin <- is.finite(stops)
  out[fin] <- format(test_dates[pmin(as.integer(stops[fin]), N_test)])
  out
}
results_df <- data.frame(
  site      = rownames(norm_matrix),
  uniform   = stop_to_date(vapply(det_unif$stream_results, `[[`, numeric(1L), "stopping_time")),
  proximity = stop_to_date(vapply(det_prox$stream_results, `[[`, numeric(1L), "stopping_time")),
  row.names = NULL
)
message("Detection results:")
print(results_df)

# ── 5. Figure 1: Transformed signal, 2023-2025 ──────────────────────────────

# Facet-label-only shortening (does not touch the `site` column used for
# joins/matching elsewhere) so the one long official site name doesn't get
# clipped once panels narrow in the side-by-side combined figure.
short_site_label <- function(site) {
  ifelse(site == "City of Charlottetown and Town of Stratford",
        "Charlottetown/Stratford", site)
}

p_signal <- pipe_long %>%
  mutate(period = if_else(weekstart < train_end, "2023 (training)", "2024-2025 (monitoring)"),
         site_label = short_site_label(site)) %>%
  ggplot(aes(x = weekstart, y = norm, colour = period, group = site)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = train_end, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  facet_wrap(~ site_label, ncol = 3) +
  scale_colour_manual(values = c("2023 (training)"        = "gray45",
                                 "2024-2025 (monitoring)" = "#2E86AB"),
                      name = NULL) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b\n'%y") +
  labs(x = NULL, y = Y_LABEL,
       title = "Per-site-normalized log1p(covN2), 2023-2025, by selected site",
       subtitle = "Dotted = eta_k = 1 (each site's own 2023 average); dashed = train/monitor split") +
  theme_bw(base_size = FONT_SIZE) +
  theme(legend.position   = "top",
        strip.background  = element_rect(fill = "#D6EAF8"),
        strip.text        = element_text(face = "bold", size = rel(0.75)),
        axis.text.x       = element_text(size = rel(0.65)),
        panel.grid.minor  = element_blank())

ggsave(file.path("figures", "wastewater_transformed_signal.png"), p_signal,
      width = 12, height = 11, dpi = 300)
message("Wrote figures/wastewater_transformed_signal.png")

# ── 6. Figure 2: Alarm plot with a 3-colour uniform/proximity/both scheme ──

site_names <- rownames(norm_matrix)

alarm_category <- results_df %>%
  mutate(
    unif_alarm = uniform   != "no alarm",
    prox_alarm = proximity != "no alarm",
    category = case_when(
      unif_alarm & prox_alarm  ~ "Both allowances",
      unif_alarm & !prox_alarm ~ "Uniform only",
      !unif_alarm & prox_alarm ~ "Proximity only",
      TRUE ~ NA_character_
    )
  )

alarm_df <- alarm_category %>%
  filter(!is.na(category)) %>%
  pivot_longer(cols = c(uniform, proximity), names_to = "allowance", values_to = "stop_date") %>%
  filter(stop_date != "no alarm") %>%
  mutate(alarm_date = as.Date(stop_date))

alarm_colours <- c("Uniform only"    = "#457b9d",
                   "Proximity only"  = "#e07020",
                   "Both allowances" = "#6A3D9A")

norm_test_df <- as.data.frame(t(norm_matrix[, test_idx, drop = FALSE])) %>%
  tibble::rownames_to_column("date_chr") %>%
  tidyr::pivot_longer(-date_chr, names_to = "site", values_to = "norm") %>%
  dplyr::mutate(date = as.Date(date_chr), site_label = short_site_label(site))

alarm_df <- alarm_df %>% mutate(site_label = short_site_label(site))

p_alarm <- ggplot(norm_test_df, aes(x = date, y = norm)) +
  geom_line(colour = "gray30", linewidth = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40", linewidth = 0.6) +
  geom_vline(data = alarm_df, aes(xintercept = alarm_date, colour = category),
             linetype = "dashed", linewidth = 0.9) +
  scale_colour_manual(values = alarm_colours, name = "Alarms under:") +
  facet_wrap(~ site_label, ncol = 3) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n'%y") +
  labs(x = NULL, y = Y_LABEL,
       title = "2024-2025 normalized signal by site, with detection dates",
       subtitle = "Dotted = eta_k = 1 (each site's own 2023 average)") +
  theme_bw(base_size = FONT_SIZE) +
  theme(legend.position   = "top",
        strip.background  = element_rect(fill = "#D6EAF8"),
        strip.text        = element_text(face = "bold", size = rel(0.75)),
        axis.text.x       = element_text(size = rel(0.65)),
        panel.grid.minor  = element_blank())

ggsave(file.path("figures", "wastewater_alarm_plot.png"), p_alarm,
      width = 12, height = 11, dpi = 300)
message("Wrote figures/wastewater_alarm_plot.png")

# ── 7. Combined (a)/(b) figure for the manuscript ───────────────────────────
# Titles/subtitles are dropped here (kept on the standalone PNGs above) since
# the external LaTeX caption carries that information for a combined figure.

# NB: `|` binds looser than `+` in R, so the (A | B) combination must be
# fully parenthesized before `+ plot_annotation(...)` is appended -- otherwise
# the annotation silently attaches to p_alarm alone and no tags are drawn.
p_combined <- ((p_signal + labs(title = NULL, subtitle = NULL)) |
               (p_alarm  + labs(title = NULL, subtitle = NULL))) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = rel(1.3), face = "bold"))

ggsave(file.path("figures", "wastewater_combined.png"), p_combined,
      width = 22, height = 12, dpi = 300, limitsize = FALSE)
message("Wrote figures/wastewater_combined.png")
