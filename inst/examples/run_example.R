## Example usage of calculate_methane_flux using (inst\examples\sample_data.csv)
## This script is intended to be copy/pasted or sourced from the repository root
## and does not assume the code is installed as a package. Make sure the
## function calculate_methane_flux is defined in your R session (source the
## R/calculate_methane_flux.R file if needed).

#Data frame format for running algorithm
sample_data = data.frame(
  id = flux_df$site,#unique identifier for each chamber deployment
  time = as.numeric(flux_df$numbered_column),#time of measurement in seconds from start of deployment
  concentration = flux_df$CH4_ppm,#make sure CH4 is in ppm not ppb
  volume = flux_df$volume,#volume of chamber
  area = flux_df$area,#area of chamber base
  temperature_celsius = flux_df$air_temperature,#temperature in celsius
  air_pressure = flux_df$air_pressure,#air pressure in Pa (or other units, see argument)
  gas_pressure = flux_df$cavity_pres_kpa,# here and below are predictor variables for BRT isolating ebullition events (add as many as you like)
  gas_temp = flux_df$Cavity_temp,
  ph = flux_df$pH,
  orp = flux_df$Redox.Potential,
  wt = flux_df$Water.Temperature
)

predictors = c(8:12) # these are predictors you have included for BRT fitting (column numbers in sample_data) 

#Run function
calculate_methane_flux(sample_data,
                       results_directory = "put your path here for .csv files",
                       save_directory = "put your path here for reviewing diagnostic plots",
                       plot_diagnostics = T,
                       unit = "nmol/m2/s", 
                       # start_id = "your starting id here if you want to start in the middle of a large dataset",
                       r_squared_diffusive = 0.99, 
                       window_size_concentration_smoother = 3,
                       window_size_residuals_smoother = 3,
                       quantile_threshold = 0.90, 
                       dynamic_multiplier = 2,
                       predictors = predictors,
                       learning_rate = c(0.001, 0.002, 0.003, 0.004, 0.005),
                       tree_complexity = 5, 
                       bag_fraction = 1,
                       gam_knots = 5,
                       min_ebullition_sequence = 20,
                       pressure_units = "mbar")
 
