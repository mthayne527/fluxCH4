## Example usage of calculate_methane_flux using packaged example data
library(fluxCH4)

example_file <- system.file("extdata", "example.csv", package = "fluxCH4")
if (file.exists(example_file) && file.info(example_file)$size > 0) {
	df <- read.csv(example_file)
} else {
	# fallback synthetic data
	set.seed(1)
	n <- 200
	df <- data.frame(
		id = rep("1", n),
		time = seq(1, n),
		temperature_celsius = rep(20, n),
		volume = rep(0.01, n),
		area = rep(0.1, n),
		air_pressure = rep(101325, n),
		concentration = 1 + 0.001 * seq(1, n) + rnorm(n, 0, 0.01)
	)
}

out <- calculate_methane_flux(df, results_directory = tempdir(), save_directory = tempdir(), plot_diagnostics = FALSE)
print(out)
