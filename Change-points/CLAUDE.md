# CLAUDE.md — operating instructions for the Change-points project

This file tells Claude how to work in this repo. Read it at the start of every session.

## Project at a glance

Jake (postdoc, supervisor: Bouchra) is writing a methods paper on **sequential multivariate change-point detection** built on test supermartingales (TSMs) plugged into the Shiryaev–Roberts / Shiryaev–Roberts–Pollak procedure, with three combining strategies (product, average, universal portfolio) and several model classes (Gaussian, Gaussian AR(1), AR(p), multivariate Gaussian, Bernoulli, betting TSMs).

The flagship empirical application is **Canadian wastewater surveillance** (SARS-CoV-2 viral concentrations from Health Canada's National Wastewater Monitoring Program). Each city/site is a stream; the goal is to flag sustained departures from baseline that may indicate variant emergence or other actionable signals. Simulations target **VAR(p)** processes.

**Near-term deadline:** end of May 2026 conference. Aim for a mix of simulation results and interim empirical findings on the wastewater signals. Final two weeks of May are for assembling a presentation summarizing methods, simulations, and empirical conclusions.

## Repository layout

The two repos live in different parents on Jake's machine:

- **Code:** `~/Dropbox/CANSSI/multi-change-points/Change-points/multivariate-change-point-detection/` — the R package `multichangepoints` (public: https://github.com/spertus/multivariate-change-point-detection)
  - `R/` — package source (models, tsm, detectors, combiners, betting, simulation, utils; ~2.1k LOC)
  - `tests/testthat/` — unit + integration tests (integration gated by `RUN_INTEGRATION_TESTS=true`)
  - `vignettes/workflow.Rmd` — package walk-through (Gaussian + S-R)
  - `vignettes/wastewater.Rmd` — empirical wastewater pipeline
  - `Data/wastewater_aggregate.csv` — Health Canada aggregate (~7.5 MB)
- **Manuscript:** `~/Dropbox/CANSSI/change_points_manuscript/` — LaTeX, linked to Overleaf (private: https://github.com/spertus/change_points_manuscript).

The default Cowork workspace folder is `~/Dropbox/CANSSI/multi-change-points/Change-points/`, which contains only the code repo and this CLAUDE.md. The manuscript folder is one directory level higher and **not** mounted by default — Claude can't read or edit it unless Jake selects `~/Dropbox/CANSSI/` (or the manuscript folder itself) as the workspace folder for the session.

Note: the existing vignettes call `devtools::load_all("~/Dropbox/CANSSI/multi-change-points/")`. After the reorg, the package now lives one level deeper at `~/Dropbox/CANSSI/multi-change-points/Change-points/multivariate-change-point-detection/`. If a vignette fails to load the package, update the path or run from the package directory.

There is also a leftover `_broken_clone/` directory from a sandbox-side git failure; it can be deleted from the host.

## Operating principles

### Simulation discipline (most important)

**Always run a small simulation first.** Confirm the pipeline works end-to-end on a tiny configuration (e.g., `n_rep = 20`, `N = 200`, `K = 2`, one DGP, one combiner) and verify the output looks sane before scaling up. Then **check in with Jake before launching anything large** — anything with `n_rep > 200`, `N > 2000`, multi-DGP grids, or anything expected to take more than ~2 minutes of compute. State the expected runtime and resource cost in the check-in.

Reasonable size ladder:
- **Smoke test:** seconds. One DGP, one combiner, `n_rep = 10–20`. Confirm shapes/types/no errors.
- **Pilot:** minutes. A few DGPs, full combiner set, `n_rep = 100–200`. Confirm trends look plausible.
- **Production run:** confirm with Jake first. Set seeds, log to disk, save intermediate results so a crash doesn't lose everything.

Use `set.seed(20260305)` (or another logged seed) for reproducibility. Save raw simulation output to disk (RDS) before plotting/summarizing.

### Methodology guardrails (decisions already made — do not re-litigate)

- **Feed raw log-transformed measurements** to the detector with per-site `mu_pre`/`sd_pre` from a held-out training window. **Do not** apply log-ratio differencing or ARIMA residualization; both destroy the sustained mean-shift signal.
- For staggered start times, use the **"frozen at 1" device** (`pad_offline_increments()`): set marginal TSMs to 1 during offline periods. This preserves detector validity under all three combiners. NA entries in the data matrix denote offline.
- For the **average combiner**, renormalize weights to active streams only (don't dilute with frozen streams).
- `GaussianModel` is the default; reach for `AR1Model` / `ARpModel` if baseline autocorrelation is meaningful.
- Site-level scale heterogeneity is handled by per-site marginal TSMs; if a *joint* multivariate model is used, standardize first.

### Code style and workflow

- Implementation language is **R**. The package is **S4-first** — extend existing classes rather than introducing parallel S3 code paths.
- Run unit tests after non-trivial changes:
  ```r
  testthat::test_local("multivariate-change-point-detection")
  ```
- For full integration tests:
  ```r
  Sys.setenv(RUN_INTEGRATION_TESTS = "true"); testthat::test_local("multivariate-change-point-detection")
  ```
- Don't commit or push to either repo without Jake's go-ahead. Always show the diff first.
- For the manuscript repo, prefer minimal, surgical edits — it's Overleaf-linked, so churn shows up there too.

### Communication style

Jake works at a high technical level and is comfortable with formal probability and measure-theoretic framing. Drop in proofs, lemmas, and paper-ready language when it helps. Pair theory with concrete code. Be terse; skip restating what the diff already shows.

Application choices are filtered through "genuine interest + career relevance" (potential public-health work in SF). External-collaboration scope is currently constrained by Bouchra; don't suggest reaching out to outside collaborators without checking.

## Roadmap to end of May 2026

Jake will hand off a TODO list separately. High-level shape:

- **Simulation track:** VAR(p) processes — pre/post regimes, vary K, signal strength, change-time, autocorrelation structure. Compare combiners (product / average / universal portfolio) and detectors (S-R / CUSUM).
- **Empirical track:** wastewater. Calibrate per-site pre-change parameters from training window, run the live pipeline on hold-out, characterize alarms in epidemiological context (variant emergence, wave onset).
- **Final two weeks of May:** consolidate into a conference presentation (methods → sim results → empirical findings).

When Jake gives the TODO list, capture each item as a task and revisit this CLAUDE.md if priorities shift.
