## Example: CH4 fluxes with ebullition screening for five chamber deployments.
## Run from the repository root:
##   Rscript inst/examples/run_example_CH4.R
## In RStudio the function pauses for review of each BRT classification; under Rscript the
## automatic classification is accepted and reported.

source("R/calculate_gas_flux.R")

ch4 = read.csv("inst/extdata/sample_data_CH4.csv")   # LI-7810 records; air pressure in mbar

out_dir = file.path(tempdir(), "fluxCH4_example_CH4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

res_ch4 = calculate_gas_flux(
  ch4,
  results_directory = out_dir,          # per-deployment CSV files
  save_directory = out_dir,             # diagnostic plots
  gas = "CH4",
  unit = "nmol/m2/s",
  pressure_units = "mbar",
  r_squared_diffusive = 0.99,
  window_size_concentration_smoother = 3,
  window_size_residuals_smoother = 3,
  quantile_threshold = 0.90,
  dynamic_multiplier = 2,
  predictors = c("gas_pressure", "gas_temp", "ph", "orp", "wt"),  # BRT predictors for ebullition screening
  learning_rate = c(0.001, 0.002, 0.003, 0.004, 0.005),
  tree_complexity = 5,
  bag_fraction = 1,
  gam_knots = 5,
  min_ebullition_sequence = 20,
  seed = 1
)

print(res_ch4)
message("Results and plots written to: ", out_dir)
