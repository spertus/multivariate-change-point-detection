# multichangepoints

`multichangepoints` is an S4-first R package scaffold for modular sequential change-point detection with:

- model classes (`GaussianModel`, `MultivariateGaussianModel`, `BernoulliModel`, `AR1Model`)
- test supermartingales (`SimpleVsSimpleTSM`)
- detectors (`ShiryaevRobertsDetector`, `CUSUMDetector`)
- multistream combiners (`AverageCombiner`, `ProductCombiner`, `UniversalPortfolioCombiner`)
- simulation helpers (`DGP`, `run_simulation`)

## Quick start

```r
library(multichangepoints)
set.seed(20260305)

model <- GaussianModel(mean_pre = 0, sd_pre = 1, mean_post = 1, sd_post = 1)
tsm <- SimpleVsSimpleTSM(model)

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
