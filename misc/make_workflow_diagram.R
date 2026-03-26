## make_workflow_diagram.R
## Generates a one-page landscape PDF workflow diagram for multichangepoints.
## Output: misc/workflow_diagram.pdf

outpath <- "~/Dropbox/CANSSI/multi-change-points/misc/workflow_diagram.pdf"

pdf(outpath, width = 13, height = 9)
par(mar = c(0.3, 0.3, 1.1, 0.3), bg = "white")
plot.new()
plot.window(xlim = c(0, 13), ylim = c(0, 9))
title("multichangepoints  --  Package Workflow",
      font.main = 1, cex.main = 1.55, col.main = "gray12")

# ---- colour palette -------------------------------------------------------
c_spec   <- "#D6EAF8";  b_spec   <- "#1A5276"   # blue   – specification
c_core   <- "#D5F5E3";  b_core   <- "#1D6A39"   # green  – core objects
c_data   <- "#FDEBD0";  b_data   <- "#784212"   # amber  – data
c_proc   <- "#FEF9E7";  b_proc   <- "#7D6608"   # yellow – processing
c_detect <- "#E8DAEF";  b_detect <- "#6C3483"   # purple – detection
c_opt    <- "#EAECEE";  b_opt    <- "#515A5A"   # gray   – optional
c_out    <- "#FADBD8";  b_out    <- "#78281F"   # red    – output

# ---- helpers --------------------------------------------------------------

draw_box <- function(cx, cy, w, h, lines,
                     fill = "white", border = "black",
                     cex1 = 0.88, cex2 = 0.73, lwd = 2.0,
                     lty_border = "solid") {
  rect(cx - w/2, cy - h/2, cx + w/2, cy + h/2,
       col = fill, border = border, lwd = lwd, lty = lty_border)
  n <- length(lines)
  if (n == 1L) {
    text(cx, cy, lines[1], cex = cex1, font = 2)
  } else {
    span <- h * 0.56
    ys   <- seq(cy + span / 2, cy - span / 2, length.out = n)
    text(cx, ys[1], lines[1], cex = cex1, font = 2, col = "gray10")
    for (i in seq_along(lines)[-1])
      text(cx, ys[i], lines[i], cex = cex2, font = 1, col = "gray30")
  }
}

arr <- function(x0, y0, x1, y1,
                col = "gray40", lwd = 1.8, lty = 1, label = NULL) {
  arrows(x0, y0, x1, y1, length = 0.13, angle = 20,
         col = col, lwd = lwd, lty = lty)
  if (!is.null(label)) {
    mx <- (x0 + x1) / 2
    my <- (y0 + y1) / 2
    text(mx, my + 0.20, label, cex = 0.63, col = "gray42", font = 3)
  }
}

# ---- layout geometry ------------------------------------------------------
xs <- 6.5          # spine x-centre
ws <- 4.7          # spine box width
hs <- 0.92         # spine box height (tall boxes)
hm <- 0.78         # mid-size boxes
hx <- 0.64         # small boxes

xr <- 11.05        # right column x-centre
wr <- 3.10         # right column box width

xl <- 1.75         # left column x-centre
wl <- 3.10         # left column box width

# Vertical centres (top → bottom)
y_model  <- 8.35
y_tsm    <- 6.72
y_data   <- 5.18
y_comp   <- 3.62
y_detect <- 2.08
y_out    <- 0.70

# Edges used for arrows
spine_L  <- xs - ws / 2
spine_R  <- xs + ws / 2
dgp_R    <- xl + wl / 2
comb_L   <- xr - wr / 2
comb_R   <- xr + wr / 2
dcfg_L   <- xr - wr / 2

# ---- BOXES ----------------------------------------------------------------

## (1) Model Definition — spine top
draw_box(xs, y_model, ws, hs, fill = c_spec, border = b_spec,
         lines = c("Model Definition",
                   "GaussianModel  |  BernoulliModel  |  AR1Model",
                   "MultivariateGaussianModel",
                   "pre-change and post-change parameters",
                   "simple vs simple  or  composite (predictable plug-in / mixture)"))

## (2) TSM — spine
draw_box(xs, y_tsm, ws, hm, fill = c_core, border = b_core,
         lines = c("TSM  (Test Supermartingale)",
                   "TSM(model)  or  BettingTSM(bets, eta)",
                   "encapsulates the likelihood-ratio increment structure"))

## (3) Data — spine
draw_box(xs, y_data, ws - 0.8, hx, fill = c_data, border = b_data,
         lines = c("Data  :  x",
                   "numeric vector (univariate)  or  N x K matrix (multivariate)"))

## (4) compute_increments — spine
draw_box(xs, y_comp, ws, hm, fill = c_proc, border = b_proc,
         lines = c("compute_increments(tsm, x)",
                   "log-likelihood ratio increments   Lambda_t := log p1(x_t) - log p0(x_t)"))

## (5) Detector — spine
draw_box(xs, y_detect, ws, hm, fill = c_detect, border = b_detect,
         lines = c("run_detector(detector, increments)",
                   "ShiryaevRobertsDetector  |  CUSUMDetector",
                   "update_detector() for online / streaming use"))

## (6) Output — spine bottom
draw_box(xs, y_out, ws, hx, fill = c_out, border = b_out,
         lines = c("Output",
                   "statistic path  |  stopping_time  |  alarm flag"))

## (7) DGP — left, optional, same y as Data
draw_box(xl, y_data, wl, hs, fill = c_opt, border = b_opt, lty_border = "dashed",
         lines = c("DGP  (optional)",
                   "default_gaussian_dgp",
                   "default_multivariate_gaussian_dgp",
                   "generate_stream(dgp, N, K)"))

## (8) Combiner — right, optional, same y as compute_increments
draw_box(xr, y_comp, wr, hs, fill = c_opt, border = b_opt, lty_border = "dashed",
         lines = c("Combiner  (multi-stream)",
                   "combine_streams(combiner, streams)",
                   "ProductCombiner",
                   "AverageCombiner  |  UniversalPortfolioCombiner"))

## (9) Detector Configuration — right, same y as Detector
draw_box(xr, y_detect, wr, hs, fill = c_detect, border = b_detect,
         lines = c("Detector Config",
                   "alpha  (control level)",
                   "criterion:  ARL  or  PFA",
                   "spending schedule  (if PFA)"))

# ---- ARROWS ---------------------------------------------------------------

## Spine: vertical chain
arr(xs, y_model  - hs / 2,  xs, y_tsm    + hm / 2)
arr(xs, y_tsm    - hm / 2,  xs, y_data   + hx / 2)
arr(xs, y_data   - hx / 2,  xs, y_comp   + hm / 2)
arr(xs, y_comp   - hm / 2,  xs, y_detect + hm / 2)
arr(xs, y_detect - hm / 2,  xs, y_out    + hx / 2)

## DGP → Data  (horizontal, dashed)
arr(dgp_R, y_data, spine_L + 0.4, y_data,
    col = b_opt, lty = 2, label = "generate_stream()")

## compute_increments → Combiner  (horizontal right)
arr(spine_R, y_comp, comb_L, y_comp,
    col = b_opt, lty = 2, label = "K streams")

## Combiner → Detector  (diagonal down-left)
arr(xr, y_comp - hs / 2,
    spine_R - 0.15, y_detect + hm / 2,
    col = b_opt, lty = 2)

## Detector Config → Detector  (horizontal left)
arr(dcfg_L, y_detect, spine_R, y_detect, col = b_detect)

# ---- LEGEND ---------------------------------------------------------------
legend(0.15, 2.5,
       legend = c("Specification", "Core objects", "Data",
                  "Processing", "Detection / config",
                  "Optional (multi-stream)", "Output"),
       fill   = c(c_spec, c_core, c_data, c_proc, c_detect, c_opt, c_out),
       border = c(b_spec, b_core, b_data, b_proc, b_detect, b_opt, b_out),
       bty = "n", cex = 0.78,
       title = "Component type", title.font = 2, title.cex = 0.82)

## "optional" italic labels above optional boxes
text(xl,  y_data  + hs / 2 + 0.17, "optional", cex = 0.67, col = b_opt, font = 3)
text(xr,  y_comp  + hs / 2 + 0.17, "optional", cex = 0.67, col = b_opt, font = 3)

## Dashed border note
text(0.15, 0.92, "- - -  dashed border = optional path",
     cex = 0.67, col = b_opt, font = 3, adj = 0)

dev.off()
message("Written: ", outpath)
