## Example: CO2 fluxes (diffusive only) for five chamber deployments with positive exchange.
## Run from the repository root:
##   Rscript inst/examples/run_example_CO2.R
## Ebullition screening is off by default for CO2; negative fluxes (apparent uptake) would be kept.

source("R/calculate_gas_flux.R")

co2 = read.csv("inst/extdata/sample_data_CO2.csv")   # ABB/LGR records; air pressure in Pa

out_dir = file.path(tempdir(), "fluxCH4_example_CO2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

res_co2 = calculate_gas_flux(
  co2,
  results_directory = out_dir,
  save_directory = out_dir,
  gas = "CO2",
  unit = "nmol/m2/s",
  pressure_units = "Pa",
  gam_knots = 5,
  require_manual_review = FALSE   # set TRUE to review low-R2 fits interactively
)

print(res_co2)
message("Results and plots written to: ", out_dir)
