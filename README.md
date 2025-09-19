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
# fluxCH4

Small collection of R scripts to calculate methane flux from enclosure concentration
time-series. This repository is organized as standalone code (not an R package). The
main analysis function is provided in `R/calculate_methane_flux.R` and a runnable
example is supplied at `inst/examples/run_example.R` which uses the example data in
`R/sample_data.csv`.

**Dependencies:**
- Install required packages before running the example:

```r
install.packages(c("mgcv", "dismo", "gbm", "ggplot2", "scales"))
```

**Quick Start (copy/paste friendly)**
- From the repository root, open an R session and source the main function:

```r
source("R/calculate_methane_flux.R")
source("inst/examples/run_example.R")  # example will use R/sample_data.csv
```

- Or run the example from PowerShell / command line:

```powershell
Rscript -e "source('inst/examples/run_example.R')"
```

**Data**
- Example data is provided at `R/sample_data.csv`. The example script
	expects columns including: `id`, `time`, `concentration` (CH4 in ppm), `volume`,
	`area`, `temperature_celsius`, and `air_pressure`. Additional predictor columns
	(e.g., `gas_pressure`, `gas_temp`, `ph`, `orp`, `wt`) are optional and will be
	used for event detection if present.

**Notes**
- The scripts are intended to be run as scripts (copy/paste or source) rather than
	installed as a package. The example uses `tempdir()` for outputs by default so it
	is safe to run without modifying the repository.
- The code prompts for interactive review when run in an R studio session.

**Contributing**
- See `CONTRIBUTING.md` for contribution guidance.

**License**
- MIT
