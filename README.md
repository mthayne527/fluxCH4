# fluxCH4

Utilities for calculating methane flux from enclosure concentration time-series.

Installation

Install required packages first:

```r
install.packages(c("mgcv", "dismo", "gbm", "ggplot2", "scales"))
# fluxCH4

Utilities to calculate methane flux from chamber/enclosure concentration time-series. The main exported function is `calculate_methane_flux()` which detects ebullition events, fits smooth models, and returns per-sequence flux estimates.

Installation
------------

Install the R packages required by the code:

```r
install.packages(c("mgcv", "dismo", "gbm", "ggplot2", "scales"))
```

To work with this repository locally from the package root:

```r
# install devtools if needed
install.packages("devtools")
devtools::load_all('.')
```

Quick Usage
-----------

1. Prepare a data.frame with columns: `id`, `time`, `temperature_celsius`, `volume`, `area`, `air_pressure`, `concentration` (CH4 in ppm).
2. Call the function (example uses packaged example data):

```r
library(fluxCH4)
example_file <- system.file("extdata", "example.csv", package = "fluxCH4")
df <- read.csv(example_file)

# Run with diagnostics disabled for non-interactive runs
res <- calculate_methane_flux(df, results_directory = tempdir(), save_directory = tempdir(), plot_diagnostics = FALSE)
print(res)
```

Notes
-----
- The function will check for required packages and stop with an informative message if any are missing.
- In interactive sessions the function may prompt you to review and optionally reclassify detected sequences; in non-interactive runs it defaults to proceeding with the automated classification.
- Units: default output unit is `mg/m2/s`. You may request `umol/m2/s` or `nmol/m2/s` via the `unit` argument.

Testing and Continuous Integration
---------------------------------

- Unit tests (basic) are included under `tests/testthat/`.
- A GitHub Actions workflow (`.github/workflows/R-CMD-check.yaml`) is included to run `R CMD check` on push and pull requests.

Contributing
------------

See `CONTRIBUTING.md` for contribution guidance.

License
-------
MIT
