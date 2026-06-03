# multichangepoints

An R package for sequential multivariate change-point detection via test supermartingales (TSMs).  Accompanies Spertus et al. (2026).

## Classes

| Layer | Classes |
|---|---|
| **Models** | `GaussianVARModel` (+ wrappers `GaussianModel`, `MultivariateGaussianModel`); `BernoulliModel`; `BoundedModel` (scalar or vector `eta` + wrapper `ARBoundedModel` for bounded AR(*p*) data) |
| **Betting strategies** | `AGRAPABet`, `EWMABet`, `FixedBet` |
| **TSM** | `TSM` (alias `SimpleVsSimpleTSM`) |
| **Detectors** | `ShiryaevRobertsDetector` (ARL and PFA criteria), `CUSUMDetector`, `LocalizedDetector` |
| **Combiners** | `AverageCombiner`, `ProductCombiner`, `UniversalPortfolioCombiner` |
| **Simulation** | `DGP`, `run_simulation`, `generate_stream` |

## Installation

```r
# Install from GitHub (requires the remotes package):
# install.packages("remotes")
remotes::install_github("spertus/multivariate-change-point-detection")

# Or install from a local clone:
# devtools::install("path/to/multi-change-points")
```

## Quick start

```r
library(multichangepoints)
set.seed(1)

# IID Gaussian change
model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
x     <- c(rnorm(100, 0, 1), rnorm(100, 1, 1))
inc   <- compute_increments(TSM(model), x, log = TRUE)

det <- ShiryaevRobertsDetector(alpha = 0.01, criterion = "ARL")
run_detector(det, inc, log = TRUE)$stopping_time
```

```r
# [0,1]-bounded data with AR(1) pre-change structure
phi <- 0.4;  mu <- 0.3;  sigma <- 0.06;  N <- 500
x   <- numeric(N);  x[1] <- rnorm(1, mu, sigma / sqrt(1 - phi^2))
for (t in 2:N) x[t] <- mu + phi * (x[t-1] - mu) + rnorm(1, 0, sigma)

# ARBoundedModel corrects the conditional null mean for each time step,
# making the TSM a valid martingale under the AR(1) pre-change model.
bm  <- ARBoundedModel(phi = phi, mu = mu, x = x,
                      bets = EWMABet(rho = 0.1, mu_init = mu))
inc <- compute_increments(TSM(bm), x, log = TRUE)
run_detector(ShiryaevRobertsDetector(alpha = 0.05), inc, log = TRUE)$stopping_time
```

## Validity notes

- **`BoundedModel` with scalar `eta`** is theoretically valid only for IID pre-change data.  Using a fixed `eta` on AR(*p*) data (φ > 0) produces a submartingale under H₀, inflating false alarm rates.
- **`ARBoundedModel(phi, mu, x, bets)`** computes the time-varying conditional null mean η_t = E[X_t | X_{t-1},…,X_{t-p}] and passes it as a vector to `BoundedModel`, restoring validity.

## Wastewater vignette

```r
vignette("wastewater", package = "multichangepoints")
# or, from source:
rmarkdown::render("vignettes/wastewater.Rmd")
```

## Testing

```r
# Unit tests (fast, ~5 s):
devtools::test()

# Include integration/calibration tests:
Sys.setenv(RUN_INTEGRATION_TESTS = "true")
devtools::test()
```
