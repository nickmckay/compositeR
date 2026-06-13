# CLAUDE.md — compositeR

`compositeR` creates paleoclimate composite records from multiple LiPD timeseries. As of June
2026 it is part of a six-package family and builds on **ens** (binning, gaussianize, BAM age
simulation) and **lipdViz** (plotting) — it no longer depends on the monolithic geoChronR.

Repo: nickmckay/compositeR. Branch: `refactor`.

## Package family (dependency DAG: ens ← lipdViz ← geoChronR; actR, compositeR & fluxcapacitoR on top)

| Repo (`~/GitHub/...`) | GitHub | Branch | Role |
|---|---|---|---|
| ens | nickmckay/ens | main | Ensemble methods: `bin`, `gaussianize`, `simulateBam` |
| lipdViz | nickmckay/lipdViz | main | Plotting (`plotTimeseriesEnsRibbons`, `mapTs`) |
| geoChronR-chronOnly | nickmckay/geoChronR-chronOnly | main | Age modeling |
| actR | **LinkedEarth/actR** | refactor | Abrupt-change detection |
| **compositeR** (this repo) | nickmckay/compositeR | refactor | Record compositing |
| fluxcapacitoR | nickmckay/fluxcapacitoR | main | Flux-focused varve age modeling |

## Architecture

Compositing methods are **pluggable via function-parameter dispatch**:
`compositeEnsembles2(fTS, binvec, binFun = ..., stanFun = ..., ...)`. Custom binning
(`simpleBinTs`, `sampleEnsembleThenBinTs`) or standardization (`standardizeOverInterval`,
`standardizeOverRandomInterval`, `standardizeMeanIteratively`) plug in without modifying the
package. Result is an S3 `paleoComposite` (a classed list of matrices) with `print`/`summary`
methods — the matrix-shaped half of the family's two-layer result convention.

- `compositeEnsembles()` (older single-ensemble API) is **deprecated**, superseded by
  `compositeEnsembles2()`.
- Output exposes a correctly-spelled `parameters` element; the historical misspelling
  `paramaters` is kept as an identical copy for one deprecation cycle — remove it later.

## Gotchas

- **Golden tests** (`tests/testthat/test-golden.R`, ref via `tools/capture_reference.R`) —
  compositeR's first tests; lock binning/standardization/compositeEnsembles2 numerics.
  Mostly R-RNG-portable (no rEDM), so simpler than actR's.
- Uses `simulateAutoCorrelatedUncertainty` with sd-first arg order (historical); it's now a
  thin wrapper over `ens::simulateAutoCorrelatedUncertainty` (which is n-first).
- CI: R CMD check on Windows/Linux/macOS with vignettes.

## Dev

`devtools::load_all()` · `devtools::document()` · `devtools::test()` · `devtools::check()`.
Needs `ens` + `lipdViz` installed from source. Commit work when complete.
Co-author trailer: `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
