# wastewater_proximity_graph.R
# Builds the proximity-graph figure for Section 8 (Application) of the
# manuscript: the K = 13 NWMP sites used by vignettes/wastewater.Rmd's
# detection pipeline (see "Site selection" there), taken -- for now -- as a
# fully connected ProximityGraph whose edge weights are the great-circle
# distance (km) between sites. Geocoding reuses the same
# world.cities-join-plus-manual-coordinates approach as the vignette's
# "Sampling locations" map, restricted to the 13 selected sites.
#
# Run from the ssc-2026/ folder:
#   Rscript wastewater_proximity_graph.R
#
# Output: figures/wastewater_proximity_graph.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(canadianmaps)
  library(maps)
  library(ggrepel)
})

PKG_DIR <- normalizePath(file.path(".."), mustWork = FALSE)
if (requireNamespace("devtools", quietly = TRUE)) {
  suppressMessages(devtools::load_all(PKG_DIR))
} else {
  library(multichangepoints)
}

# ── 1. Reproduce the vignette's site-selection criteria ─────────────────────
# (enrolled by 2023-01-01, observed through 2026-01-01, >= 80% complete data
#  for both covN2 and RSV in that window -- see wastewater.Rmd, "Site
#  selection") to get the same K = 13 sites used by the detection pipeline.

ww <- suppressMessages(read_csv(file.path(PKG_DIR, "Data", "wastewater_aggregate.csv"),
                                show_col_types = FALSE)) %>%
  mutate(weekstart = as.Date(weekstart))

window_start <- as.Date("2023-01-01")
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

cat("Selected sites (K =", length(selected_sites), "):\n")
cat(paste(selected_sites, collapse = "\n"), "\n\n")

site_city <- ww %>%
  filter(site %in% selected_sites) %>%
  distinct(site, city, province) %>%
  mutate(city_clean = trimws(city))

# ── 2. Geocode: world.cities join, then manual overrides for sites whose
#      `city` field doesn't match a world.cities name (metro-area facility
#      names, bilingual/compound names) -- same fallback list as the
#      vignette's "Sampling locations" map (subset relevant to these 13). ──

ca_cities <- world.cities[world.cities$country.etc == "Canada", ]

site_geo <- site_city %>%
  left_join(
    ca_cities %>% select(name, lat, long, pop),
    by = c("city_clean" = "name"),
    relationship = "many-to-many"
  ) %>%
  group_by(site) %>%
  slice_max(pop, n = 1, with_ties = FALSE) %>%
  ungroup()

manual_coords <- tribble(
  ~site,                                        ~lat,    ~long,
  "City of Charlottetown and Town of Stratford",46.24,  -63.13,
  "St. John's",                                 47.56,  -52.71,
  "Vancouver Annacis Island",                   49.17, -122.92,
  "Vancouver Iona Island",                      49.22, -123.21
)

site_geo <- site_geo %>%
  rows_update(manual_coords, by = "site", unmatched = "ignore") %>%
  filter(!is.na(lat)) %>%
  arrange(site)

stopifnot(nrow(site_geo) == length(selected_sites))

# ── 3. Fully connected ProximityGraph: edge weight = great-circle distance
#      (km, haversine formula, Earth radius 6371 km) between every pair of
#      sites. All pairwise distances between distinct Canadian sites are
#      strictly positive, so W has no structural zeros off the diagonal --
#      i.e. the graph is fully connected by construction, matching the
#      "take it to be fully connected for now" placeholder in Section 8. ──

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# City-centroid geocoding puts same-city facilities (e.g. the two Toronto or
# two Montreal sites) at identical coordinates, giving a raw distance of 0 --
# which would both collapse to a structural "no edge" under ProximityGraph's
# W[i,j] == 0 convention and blow up the 1/distance plotting aesthetic below.
# Floor same-city pairs at a nominal 5 km (a reasonable order of magnitude
# for within-city treatment-plant separation, not resolved by city-level
# geocoding) so the graph stays fully connected as intended.
MIN_DIST_KM <- 5

K <- nrow(site_geo)
W <- matrix(0, K, K)
for (i in seq_len(K)) for (j in seq_len(K)) {
  if (i != j) W[i, j] <- max(MIN_DIST_KM,
                             haversine_km(site_geo$lat[i], site_geo$long[i],
                                         site_geo$lat[j], site_geo$long[j]))
}
proximity_graph <- ProximityGraph(W, names = site_geo$site)

saveRDS(proximity_graph, file.path(PKG_DIR, "ssc-2026", "wastewater_proximity_graph.rds"))
cat("Edge weight (km) summary:\n")
print(summary(W[upper.tri(W)]))

# ── 4. Plot: sites on a Canada basemap, all K(K-1)/2 edges drawn (fully
#      connected), with edge alpha/width mapped to inverse distance so
#      nearby (informative) edges stay visible while the plot doesn't
#      collapse into an unreadable hairball. ─────────────────────────────

lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lon_0=-91.52 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

sites_sf <- site_geo %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = lcc_crs)
coords_proj <- st_coordinates(sites_sf)
site_geo$x <- coords_proj[, 1]
site_geo$y <- coords_proj[, 2]

# Distinguishing display label per site: city name alone for the 9 singleton
# cities, "City – Facility" for the 4 cities with two pipeline sites, so all
# 13 labels stay both unique and geographically legible.
label_map <- c(
  "City of Charlottetown and Town of Stratford" = "Charlottetown/Stratford",
  "Edmonton Goldbar"          = "Edmonton",
  "Halifax Halifax"           = "Halifax",
  "Montreal North"            = "Montreal – North",
  "Montreal South"            = "Montreal – South",
  "Regina"                    = "Regina",
  "St. John's"                = "St. John's",
  "Toronto Ashbridges Bay"    = "Toronto – Ashbridges Bay",
  "Toronto Humber"            = "Toronto – Humber",
  "Vancouver Annacis Island"  = "Vancouver – Annacis Island",
  "Vancouver Iona Island"     = "Vancouver – Iona Island",
  "Winnipeg North End"        = "Winnipeg – North End",
  "Winnipeg South End"        = "Winnipeg – South End"
)
site_geo$label <- unname(label_map[site_geo$site])

# City-centroid geocoding puts same-city facility pairs at (near-)identical
# coordinates (Toronto, Montreal, Winnipeg x2 exactly; Vancouver's two manual
# overrides are already a few km apart). Jitter duplicated plotting
# coordinates deterministically (display only -- does not touch the W matrix
# / ProximityGraph, which keep the true geocoded/floored distances) so both
# points and labels are visible.
JITTER_M <- 18000
dup_key <- paste(round(site_geo$x), round(site_geo$y))
for (k in unique(dup_key)) {
  idx <- which(dup_key == k)
  if (length(idx) > 1L) {
    ang <- seq(0, 2 * pi, length.out = length(idx) + 1L)[seq_along(idx)]
    site_geo$x[idx] <- site_geo$x[idx] + JITTER_M * cos(ang)
    site_geo$y[idx] <- site_geo$y[idx] + JITTER_M * sin(ang)
  }
}

# Edges built from the (jittered) plotting coordinates so segments terminate
# exactly at the displayed points; the underlying W/proximity_graph object
# above is unaffected.
edge_df <- do.call(rbind, lapply(seq_len(K - 1L), function(i) {
  do.call(rbind, lapply((i + 1L):K, function(j) {
    if (j > K) return(NULL)
    data.frame(x = site_geo$x[i], y = site_geo$y[i],
              xend = site_geo$x[j], yend = site_geo$y[j],
              dist_km = W[i, j])
  }))
}))

p <- ggplot() +
  geom_prov(fill = "gray95", colour = "gray65", size = 0.25) +
  geom_segment(data = edge_df,
              aes(x = x, y = y, xend = xend, yend = yend,
                  alpha = 1 / dist_km, linewidth = 1 / dist_km),
              colour = "#2E86AB") +
  scale_alpha_continuous(range = c(0.05, 0.85), guide = "none") +
  scale_linewidth_continuous(range = c(0.15, 1.3), guide = "none") +
  geom_point(data = site_geo, aes(x = x, y = y),
            colour = "#c0392b", size = 2.6, alpha = 0.9) +
  geom_text_repel(data = site_geo, aes(x = x, y = y, label = label),
                  size = 3.0, colour = "gray20", segment.size = 0.2,
                  max.overlaps = 20, seed = 1) +
  coord_sf(crs = lcc_crs) +
  labs(
    title = "Proximity graph on the K = 13 NWMP pipeline sites",
    subtitle = "Fully connected; edge weight = great-circle distance (km); darker/thicker = closer",
    caption = "Site coordinates: city centroid (world.cities) or manual override for metro-area facilities;\nsame-city facility pairs jittered apart for display only."
  ) +
  theme_map() +
  theme(
    plot.title      = element_text(face = "bold", margin = margin(b = 4)),
    plot.subtitle   = element_text(colour = "gray30", size = 9.5, margin = margin(b = 6)),
    plot.caption    = element_text(colour = "gray50", size = 7.5, margin = margin(t = 4)),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave(file.path("figures", "wastewater_proximity_graph.png"), p,
      width = 8, height = 6.4, dpi = 300, bg = "white")
message("Wrote figures/wastewater_proximity_graph.png  (K = ", K, " sites, ",
        K * (K - 1) / 2, " edges)")
