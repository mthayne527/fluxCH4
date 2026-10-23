# fluxCH4 1.1.0

## New
- `calculate_gas_flux()` supports CO₂ (`gas = "CO2"`, molecular weight 44.01) as well as CH₄.
- The function returns all results as a data frame (invisibly). New columns: `Deployment_ID`, `Gas`, `Unit`, `Fit_Type` and `Decision`.
- A CO₂ sample data set (`inst/extdata/sample_data_CO2.csv`) and separate CH₄ and CO₂ example scripts are included.
- `tests/run_checks.R` provides non-interactive checks; GitHub Actions runs them together with both examples.
- New arguments: `gas`, `require_manual_review`, `seed`, `verbose`.

## Changed behaviour
- When R is not interactive, prompts are skipped and the automatic classification is accepted (`Decision = "automatic"`). Interactive sessions behave as before.
- `bag_fraction` now defaults to `1`, which gives deterministic BRT fits. The previous default was 0.15.
- `enable_ebullition_detection = NULL` (the default) enables screening for CH₄ and disables it for CO₂.
- `start_id` is positional in sorted `id` order. Previously it was compared lexically, which skipped IDs such as `e10` when starting from `e2`.
- `results_directory` and `save_directory` are optional; nothing is written when they are `NULL`.
- Package checks use `requireNamespace()` at the start. The function no longer attaches packages.

## Fixed
- An ebullition event still open at the end of a series is extended through the trailing buffer again.
- BRT predictions use the optimal number of trees from cross-validation.
- Debug output was removed.

## Repository
- The sample data moved from `R/` to `inst/extdata/`. The CH₄ sample's index column was dropped, so it cannot be picked up as a predictor.
- `R/calculate_methane_flux.R` was replaced by `R/calculate_gas_flux.R`. `calculate_methane_flux()` remains as a wrapper.
- The `run_example.R` script that broke CI was replaced. The LICENSE text was completed, and `CITATION.cff` was added.

# fluxCH4 1.0.1

- Version archived on Zenodo (https://doi.org/10.5281/zenodo.18241611) and cited in Thayne et al. (2026).
