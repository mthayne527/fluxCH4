## Non-interactive checks run by GitHub Actions:  Rscript tests/run_checks.R
source("R/calculate_gas_flux.R")

check = function(ok, what) {
  if (!isTRUE(ok)) stop("CHECK FAILED: ", what, call. = FALSE)
  message("ok - ", what)
}
cols = c("ID", "Deployment_ID", "Gas", "Flux_Type", "Flux_Estimate", "Unit",
         "R_squared", "RMSE", "Fit_Type", "Decision")

## CO2: diffusive-only default, no output directories
co2 = read.csv("inst/extdata/sample_data_CO2.csv")
res = calculate_gas_flux(co2, gas = "CO2", unit = "nmol/m2/s", pressure_units = "Pa", verbose = FALSE)
check(is.data.frame(res) && identical(names(res), cols), "CO2 returns a data.frame with the documented columns")
check(setequal(res$Deployment_ID, unique(co2$id)), "CO2: one result per deployment")
check(all(is.finite(res$Flux_Estimate)) && all(res$Flux_Estimate > 0), "CO2 sample fluxes are finite and positive")
check(all(res$Gas == "CO2") && all(res$Flux_Type == "Diffusive"), "CO2 rows are diffusive CO2 estimates")

## Independent sanity check: linear-slope flux of each deployment within a factor of 2
lin = sapply(split(co2, co2$id), function(d) {
  mol = d$concentration / 1e6 * d$air_pressure * d$volume[1] / (8.314 * (d$temperature_celsius + 273.15))
  unname(coef(lm(mol ~ d$time))[2]) / d$area[1] * 1e9
})
ratio = res$Flux_Estimate[match(names(lin), res$Deployment_ID)] / lin
check(all(ratio > 0.5 & ratio < 2), "CO2 GAM fluxes agree with linear slopes within a factor of 2")

## Units scale consistently
res_mg = calculate_gas_flux(co2, gas = "CO2", unit = "mg/m2/s", pressure_units = "Pa", verbose = FALSE)
check(isTRUE(all.equal(res_mg$Flux_Estimate, res$Flux_Estimate * 1e-9 * 44.01 * 1000, tolerance = 1e-8)),
      "mg/m2/s equals nmol/m2/s x molecular weight")

## start_id is positional in the sorted deployment order
ids = sort(unique(co2$id))
res_start = calculate_gas_flux(co2, gas = "CO2", unit = "nmol/m2/s", pressure_units = "Pa",
                               start_id = ids[3], verbose = FALSE)
check(setequal(res_start$Deployment_ID, ids[3:length(ids)]), "start_id processes the remaining deployments")

## CH4 with ebullition screening, non-interactive acceptance, files written
ch4 = read.csv("inst/extdata/sample_data_CH4.csv")
out_dir = file.path(tempdir(), "fluxCH4_checks"); dir.create(out_dir, showWarnings = FALSE)
res4 = calculate_gas_flux(ch4, results_directory = out_dir, save_directory = out_dir, gas = "CH4",
                          unit = "nmol/m2/s", pressure_units = "mbar",
                          predictors = c("gas_pressure", "gas_temp", "ph", "orp", "wt"),
                          window_size_concentration_smoother = 3, window_size_residuals_smoother = 3,
                          quantile_threshold = 0.90, min_ebullition_sequence = 20, seed = 1, verbose = FALSE)
check(identical(names(res4), cols) && nrow(res4) >= 1, "CH4 returns results")
check(all(is.finite(res4$Flux_Estimate)), "CH4 fluxes are finite")
check(all(unique(ch4$id) %in% res4$Deployment_ID), "CH4: every deployment has a result")
check(length(list.files(out_dir, pattern = "_CH4\\.csv$")) == length(unique(res4$Deployment_ID)),
      "CH4 per-deployment CSV files are written")

## Backward-compatible wrapper
one = ch4[ch4$id == unique(ch4$id)[1], ]
a = calculate_methane_flux(one, unit = "nmol/m2/s", pressure_units = "mbar", enable_ebullition_detection = FALSE, verbose = FALSE)
b = calculate_gas_flux(one, gas = "CH4", unit = "nmol/m2/s", pressure_units = "mbar", enable_ebullition_detection = FALSE, verbose = FALSE)
check(isTRUE(all.equal(a, b)), "calculate_methane_flux() matches calculate_gas_flux(gas = 'CH4')")

message("All checks passed.")
