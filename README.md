# multichangepoints

`multichangepoints` is an R package scaffold for multivariate sequential change-point detection via test supermartingales. It is modularized into: 

- model classes (`GaussianModel`, `MultivariateGaussianModel`, `BernoulliModel`, `AR1Model`)
- test supermartingales (`TSM`)
- detectors (`ShiryaevRobertsDetector`, `CUSUMDetector`)
- multistream combiners (`AverageCombiner`, `ProductCombiner`, `UniversalPortfolioCombiner`)
- simulation helpers (`DGP`, `run_simulation`)

## Quick start

```r
library(multichangepoints)
set.seed(20260305)

model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
tsm <- TSM(model)

dgp <- DGP(
  generator = default_gaussian_dgp,
  pre_params = list(mean = 0, sd = 1),
  post_params = list(mean = 1, sd = 1),
  nu = 120,
  name = "iid-gaussian-change"
)

x <- generate_stream(dgp, N = 250, K = 1)
inc <- compute_increments(tsm, x)

det <- ShiryaevRobertsDetector(alpha = 0.05, criterion = "ARL")
out <- run_detector(det, inc)
out$stopping_time
```

## Workflow vignette

```r
vignette("workflow", package = "multichangepoints")
```

If the package is not installed yet, render the source vignette directly:

```r
rmarkdown::render("vignettes/workflow.Rmd")
```

## Testing

Quick unit tests (default; skips integration tests):

```r
testthat::test_local()
```

Full suite including integration/calibration tests:

```r
Sys.setenv(RUN_INTEGRATION_TESTS = "true")
testthat::test_local()
```
