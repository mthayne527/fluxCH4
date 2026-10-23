# fluxCH4: chamber flux algorithm for CH₄ and CO₂

[![Checks and examples](https://github.com/mthayne527/fluxCH4/actions/workflows/run-example.yaml/badge.svg)](https://github.com/mthayne527/fluxCH4/actions/workflows/run-example.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18241611.svg)](https://doi.org/10.5281/zenodo.18241611)

`calculate_gas_flux()` estimates diffusive (and, for CH₄, ebullitive) gas fluxes from closed-chamber
concentration time series. Each deployment is fitted with a generalized additive model (GAM). When
ebullition screening is enabled, boosted regression trees (BRT) learn the diffusive background from
auxiliary sensor data. Residual spikes above a dynamic quantile threshold are then classified as
ebullition events. Diffusive and ebullitive sequences are fitted separately.

The method is described in Thayne et al. (2026), *Biogeosciences* 23, 477–495, https://doi.org/10.5194/bg-23-477-2026.

Version 1.1.0 adds CO₂ (`gas = "CO2"`), returns the results as a data frame, and runs non-interactively
under `Rscript`. `calculate_methane_flux()` still works and calls `calculate_gas_flux(gas = "CH4")`.

## Requirements

R ≥ 4.1 with `mgcv`, `ggplot2` and `scales`. `dismo` and `gbm` are needed only when ebullition screening is enabled.

```r
install.packages(c("mgcv", "ggplot2", "scales", "dismo", "gbm"))
```

## Quick start

```r
source("R/calculate_gas_flux.R")

co2 = read.csv("inst/extdata/sample_data_CO2.csv")
res = calculate_gas_flux(co2, gas = "CO2", unit = "nmol/m2/s", pressure_units = "Pa")

ch4 = read.csv("inst/extdata/sample_data_CH4.csv")
res = calculate_gas_flux(ch4, gas = "CH4", unit = "nmol/m2/s", pressure_units = "mbar",
                         predictors = c("gas_pressure", "gas_temp", "ph", "orp", "wt"),
                         window_size_concentration_smoother = 3, window_size_residuals_smoother = 3,
                         quantile_threshold = 0.90, min_ebullition_sequence = 20, seed = 1)
```

Complete scripts are in `inst/examples/` (`run_example_CH4.R`, `run_example_CO2.R`). Run them from the repository root.

## Input data

One row per measurement; deployments are stacked and identified by `id`.

| Column | Description |
|---|---|
| `id` | Deployment identifier |
| `time` | Elapsed time within the deployment (s) |
| `concentration` | Chamber headspace mole fraction (ppm) of the selected gas |
| `volume` | Chamber headspace volume (m³) |
| `area` | Chamber footprint area (m²) |
| `temperature_celsius` | Headspace air temperature (°C) |
| `air_pressure` | Air pressure, in the unit given by `pressure_units` |
| *predictors* | Optional numeric columns for BRT ebullition screening, given by name or column number in `predictors` |

Supported `pressure_units`: `Pa`, `kPa`, `atm`, `bar`, `mbar`, `mmHg`, `psi`. Supported `unit`: `mg/m2/s`, `umol/m2/s`, `nmol/m2/s`.

## Main arguments

| Argument | Default | Meaning |
|---|---|---|
| `gas` | `"CH4"` | `"CH4"` or `"CO2"`; sets molecular weight and the ebullition default |
| `enable_ebullition_detection` | `NULL` | `NULL` gives `TRUE` for CH₄ and `FALSE` for CO₂ |
| `r_squared_diffusive` | `0.99` | GAM fits at or above this value are accepted without screening |
| `gam_knots` | `5` | Basis dimension *k* of the GAM smooth |
| `quantile_threshold`, `dynamic_multiplier` | `0.95`, `2` | Residual threshold for ebullition events |
| `window_size_concentration_smoother`, `window_size_residuals_smoother` | `0`, `5` | Moving-average windows (0 = none) |
| `min_ebullition_sequence`, `buffer_size` | `8`, `5` | Minimum event length and buffer (records) |
| `learning_rate`, `tree_complexity`, `bag_fraction` | 0.001–0.005, `5`, `1` | BRT settings (`dismo::gbm.step`) |
| `start_id` | `NULL` | Resume from this deployment, in sorted `id` order |
| `results_directory`, `save_directory` | `NULL` | Folders for per-deployment CSVs and diagnostic plots; nothing is written when `NULL` |
| `require_manual_review` | `FALSE` | Ask for review of fits below `r_squared_diffusive` (interactive sessions) |
| `seed` | `NULL` | Random seed for the BRT cross-validation, set per deployment |

## Output

The function returns a data frame invisibly, with one row per flux sequence:

`ID`, `Deployment_ID`, `Gas`, `Flux_Type` (`Diffusive` or `Ebullitive`), `Flux_Estimate`, `Unit`,
`R_squared`, `RMSE`, `Fit_Type` (`GAM` or `LM`), `Decision` (`automatic`, `reviewed` or `manual`).

`R_squared` is the adjusted R² reported by `mgcv` and can be negative for near-flat series. Low R² by
itself does not mean a poor flux: a small but stable concentration change gives a low R² and a near-zero flux.

## Interactive and non-interactive use

In an interactive session (RStudio, R console), the function shows each BRT classification and asks for
acceptance. You can also enter fitting intervals manually. Under `Rscript` or continuous integration,
nothing can be answered: the automatic classification is accepted and the row is marked `automatic`.

## Citation

Please cite the paper and the software version you used (see `CITATION.cff`):

Thayne, M. W., Kemper, K., Wille, C., Kalhori, A., and Sachs, T. (2026). Spatiotemporal variability and environmental controls on aquatic methane emissions in an Arctic permafrost catchment. *Biogeosciences*, 23, 477–495. https://doi.org/10.5194/bg-23-477-2026

## License

MIT. See `LICENSE`.
