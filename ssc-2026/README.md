# SSC 2026 --- CANSSI Postdoc Spotlight talk

Files for the 30-minute talk on Tuesday 2 June 2026.

## Layout
- `slides.tex` --- Beamer source (metropolis theme, 16:9). Compile with
  `latexmk -pdf slides.tex` or upload the whole folder to Overleaf.
- `figures.R` --- generates the explainer PNGs into `figures/` using the
  `multichangepoints` package (loaded from the sibling
  `multivariate-change-point-detection/` directory via `devtools::load_all`).
  Run once before the first compile:
  ```
  Rscript figures.R
  ```
  Required R packages: `devtools` (for `load_all`), `ggplot2`, `dplyr`,
  `tidyr`, `readr`. The package itself imports `mvtnorm`.
- `figures/` --- populated by `figures.R`.

## What the script produces

Generated from the package (real didactic content; uses
`GaussianModel`, `TSM`, `ShiryaevRobertsDetector`, `compute_increments`,
`run_detector`, `combine_streams`, `default_multivariate_gaussian_dgp`):
- `tsm_growth.png` --- six TSM log-wealth paths
- `sr_procedure.png` --- ARL-controlled S-R at two thresholds, with alarms
- `spending_schedule.png` --- three spending schedules
- `spending_schedule_detectors.png` --- the PFA-controlled S-R statistic under each
- `sparse_vs_dense.png` --- multivariate streams (sparse vs. dense change)
- `combiners_paths.png` --- product / average / universal-portfolio S-R paths
- `wastewater_motivation.png` --- real NWMP COVID and RSV signals at the four
  best-observed sites, read from
  `multivariate-change-point-detection/Data/wastewater_aggregate.csv`
- `sim_regret_placeholder.png` --- a *small pilot* simulation (60 reps per
  cell, 5 magnitudes, K=6) of CAD by combiner under sparse vs. dense change.
  Scale up to conference quality by setting `BIG_SIM <- TRUE` near the top
  of `figures.R` (defaults to 500 reps, 10 magnitudes); the pilot runs in
  well under a minute, the full version will take longer.

Still placeholders (need new simulation runs or new modeling work):
- `arl_vs_pfa_placeholder.png` --- ARL vs PFA regret comparison
- `sim_arl_vs_pfa_localization_placeholder.png` --- localization regret
- `ww_detector_placeholder.png` --- detector trajectory on 2025-26
- `ww_localization_placeholder.png` --- adjacency-aware allowance cartoon

## Overleaf setup
- Compiler: pdfLaTeX (XeLaTeX/LuaLaTeX if you want Fira fonts in Metropolis)
- Main document: `slides.tex`
- Upload the whole `ssc-2026/` folder including `figures/` after running
  `figures.R` locally.
