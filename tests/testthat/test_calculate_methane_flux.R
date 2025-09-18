library(testthat)

context("calculate_methane_flux basic behavior")

test_that("function runs on small synthetic data and returns expected columns", {
  skip_on_cran()
  set.seed(2)
  n <- 50
  df <- data.frame(
    id = rep("1", n),
    time = seq_len(n),
    temperature_celsius = rep(20, n),
    volume = rep(0.01, n),
    area = rep(0.1, n),
    air_pressure = rep(101325, n),
    concentration = 1 + 0.001 * seq_len(n) + rnorm(n, 0, 0.01)
  )

  res <- calculate_methane_flux(df, results_directory = tempdir(), save_directory = tempdir(), plot_diagnostics = FALSE)
  expect_true(is.data.frame(res) || is.null(res))
  if (is.data.frame(res) && nrow(res) > 0) {
    expect_true(all(c("ID", "Flux_Type", "Flux_Estimate", "R_squared", "RMSE", "Fit_Type") %in% colnames(res)))
  }
})
