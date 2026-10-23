#' Calculate chamber-based gas fluxes (CH4 or CO2) with optional ebullition detection
#'
#' @description
#' Computes fluxes from closed-chamber concentration time series. Each chamber
#' deployment (`id`) is converted from ppm to moles with the ideal gas law and fitted
#' with a generalized additive model (GAM) of gas amount against time. The flux is the
#' mean rate of change between successive fitted values divided by the chamber basal
#' area. For CH4, deployments whose GAM fit does not reach `r_squared_diffusive` can be
#' screened for ebullition with a boosted-regression-tree (BRT) anomaly detector, and the
#' resulting diffusive/ebullition sequences are reviewed interactively. For CO2,
#' ebullition detection is off by default, so every deployment receives a diffusive flux;
#' negative values (apparent uptake) are retained.
#'
#' @details
#' **Interactive and non-interactive use.** In an interactive R session (e.g. RStudio)
#' the function prompts for review of BRT classifications, optional manual entry of
#' fitting intervals, and (if `require_manual_review = TRUE`) review of low-R2 fits.
#' When R is not interactive (`Rscript`, continuous integration), no prompt can be
#' answered; the automatic classification or GAM fit is then accepted and the
#' decision is reported with `message()`.
#'
#' **Fit statistic.** `R_squared` is the adjusted R2 reported by `mgcv` for the GAM
#' (or the R2 of the linear fallback). It can be close to zero or negative for chambers
#' without measurable accumulation, where analyzer noise dominates a flat record. Such
#' deployments are not discarded.
#'
#' **Processing order.** Rows are sorted by `id` and `time`; `start_id` starts processing
#' at that deployment in this sorted order.
#'
#' **Reproducibility.** `dismo::gbm.step` assigns cross-validation folds at random. Set
#' `seed` for reproducible ebullition detection; with `bag_fraction = 1` (default) the
#' trees themselves are deterministic.
#'
#' @param data A data.frame with columns `id`, `time` (s, strictly increasing within each
#'   `id`), `temperature_celsius`, `volume` (m3), `area` (m2), `air_pressure` and
#'   `concentration` (ppm). Additional columns may be used as BRT predictors.
#' @param results_directory Existing directory for per-deployment result CSV files, or
#'   `NULL` to skip writing CSV files.
#' @param save_directory Existing directory for diagnostic plots, or `NULL` to skip plots.
#' @param celsius_to_kelvin Celsius-to-Kelvin offset (default 273.15).
#' @param plot_diagnostics Logical; save diagnostic plots (requires `save_directory`).
#' @param unit Output flux unit: `"mg/m2/s"` (default), `"umol/m2/s"` or `"nmol/m2/s"`.
#' @param start_id Optional `id` at which to start processing (see Processing order).
#' @param r_squared_diffusive R2 above which a whole-deployment GAM fit is accepted as
#'   diffusive without further screening (default 0.99).
#' @param window_size_concentration_smoother Moving-average window applied to the
#'   concentration series before BRT screening (default 0, no smoothing).
#' @param window_size_residuals_smoother Moving-average window for BRT residuals (default 5).
#' @param quantile_threshold Quantile of absolute concentration change used to flag jumps.
#' @param dynamic_multiplier Multiplier on the residual SD for the event threshold.
#' @param predictors Optional BRT predictor columns (names or positions). If `NULL`, all
#'   non-core columns are used.
#' @param learning_rate Vector of BRT learning rates tried in order.
#' @param tree_complexity BRT tree complexity.
#' @param bag_fraction BRT bag fraction (default 1, deterministic trees).
#' @param gam_knots GAM basis dimension `k`; sequences with fewer than `3 * gam_knots`
#'   points are not fitted.
#' @param min_ebullition_sequence Minimum number of consecutive flagged points for an
#'   ebullition sequence.
#' @param R_value Gas constant (default 8.314 J mol-1 K-1).
#' @param pressure_units Units of `air_pressure`: `"Pa"`, `"kPa"`, `"atm"`, `"bar"`,
#'   `"mbar"`, `"mmHg"` or `"psi"`.
#' @param buffer_size Number of unflagged points that end an event sequence (default 5).
#' @param min_trees_stop Stop trying learning rates once a BRT has at least this many trees.
#' @param max_trees Maximum number of BRT trees.
#' @param enable_ebullition_detection Logical. `NULL` (default) means `TRUE` for CH4 and
#'   `FALSE` for CO2.
#' @param gas `"CH4"` (default) or `"CO2"`; sets the molecular weight for mass units and
#'   labels output files.
#' @param require_manual_review Logical; when ebullition detection is off, offer an
#'   interactive review of deployments whose GAM R2 is below `r_squared_diffusive`.
#' @param seed Optional integer passed to `set.seed()` before BRT fitting of each
#'   deployment, for reproducible ebullition detection.
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly, a data.frame with one row per flux estimate: `ID` (deployment, or
#'   deployment_sequence for split deployments), `Deployment_ID`, `Gas`, `Flux_Type`
#'   (`"Diffusive"` or `"Ebullition"`), `Flux_Estimate`, `Unit`, `R_squared`, `RMSE`,
#'   `Fit_Type` (`"GAM"` or `"LM"`) and `Decision` (`"automatic"`, `"reviewed"` or
#'   `"manual"`). The same rows are written to CSV when `results_directory` is set.
#'
#' @examples
#' \dontrun{
#' ch4 <- read.csv("inst/extdata/sample_data_CH4.csv")
#' res_ch4 <- calculate_gas_flux(ch4, gas = "CH4", unit = "nmol/m2/s",
#'                               pressure_units = "mbar", predictors = 8:12)
#' co2 <- read.csv("inst/extdata/sample_data_CO2.csv")
#' res_co2 <- calculate_gas_flux(co2, gas = "CO2", unit = "nmol/m2/s",
#'                               pressure_units = "Pa")
#' }
calculate_gas_flux = function(
    data,
    results_directory = NULL,
    save_directory = NULL,
    celsius_to_kelvin = 273.15,
    plot_diagnostics = TRUE,
    unit = "mg/m2/s",
    start_id = NULL,
    r_squared_diffusive = 0.99,
    window_size_concentration_smoother = 0,
    window_size_residuals_smoother = 5,
    quantile_threshold = 0.95,
    dynamic_multiplier = 2,
    predictors = NULL,
    learning_rate = c(0.001, 0.002, 0.003, 0.004, 0.005),
    tree_complexity = 5,
    bag_fraction = 1,
    gam_knots = 5,
    min_ebullition_sequence = 8,
    R_value = 8.314,
    pressure_units = "Pa",
    buffer_size = 5,
    min_trees_stop = 1000,
    max_trees = 10000,
    enable_ebullition_detection = NULL,
    gas = "CH4",
    require_manual_review = FALSE,
    seed = NULL,
    verbose = TRUE
) {

  # ---- Argument checks ----
  gas = toupper(gas)
  if (!gas %in% c("CH4", "CO2")) stop("Unsupported gas: ", gas, ". Use 'CH4' or 'CO2'.")
  unit = match.arg(unit, c("mg/m2/s", "umol/m2/s", "nmol/m2/s"))
  pressure_units = match.arg(pressure_units, c("Pa", "kPa", "atm", "bar", "mbar", "mmHg", "psi"))
  if (is.null(enable_ebullition_detection)) enable_ebullition_detection = (gas == "CH4")

  needed = c("mgcv", "ggplot2", "scales")
  if (enable_ebullition_detection) needed = c(needed, "dismo", "gbm")
  missing_pkgs = needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Install them with: install.packages(c(\"", paste(missing_pkgs, collapse = "\", \""), "\"))")
  }

  required_columns = c("id", "time", "temperature_celsius", "volume", "area", "air_pressure", "concentration")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input data must contain columns: ", paste(required_columns, collapse = ", "))
  }
  for (d in list(list(results_directory, "results_directory"), list(save_directory, "save_directory"))) {
    if (!is.null(d[[1]]) && !dir.exists(d[[1]])) stop("Specified ", d[[2]], " does not exist: ", d[[1]])
  }
  if (plot_diagnostics && is.null(save_directory)) plot_diagnostics = FALSE

  molecular_weight = c(CH4 = 16.04, CO2 = 44.01)[[gas]]
  say = function(...) if (verbose) message(...)

  # ---- Utilities ----
  ask = function(prompt, default) {
    # Interactive prompt; in non-interactive sessions the default is used and reported.
    if (interactive()) {
      cat(prompt)
      tolower(trimws(readline()))
    } else {
      message(prompt, "[non-interactive session: using '", default, "']")
      default
    }
  }

  convert_pressure_to_Pa = function(pressure_value) {
    pressure_value * c(Pa = 1, kPa = 1000, atm = 101325, bar = 1e5, mbar = 100, mmHg = 133.322, psi = 6894.76)[[pressure_units]]
  }

  convert_ppm_to_mol = function(gas_ppm, pressure_P, volume_V, temperature_kelvin) {
    (gas_ppm / 1e6) * (pressure_P * volume_V) / (R_value * temperature_kelvin)
  }

  calculate_flux = function(slope, area_A) {
    flux_mol_m2_s = slope / area_A
    switch(unit,
           "nmol/m2/s" = flux_mol_m2_s * 1e9,
           "umol/m2/s" = flux_mol_m2_s * 1e6,
           "mg/m2/s"   = flux_mol_m2_s * molecular_weight * 1000)
  }

  result_row = function(ID, deployment, flux_type, flux, r_squared, rmse, fit_type, decision) {
    data.frame(ID = as.character(ID), Deployment_ID = as.character(deployment), Gas = rep(gas, length(ID)),
               Flux_Type = flux_type, Flux_Estimate = flux, Unit = rep(unit, length(ID)), R_squared = r_squared,
               RMSE = rmse, Fit_Type = fit_type, Decision = decision, stringsAsFactors = FALSE)
  }

  write_results = function(results_df, id) {
    if (is.null(results_df) || nrow(results_df) == 0) {
      say("No results to write for ID ", id)
      return(NULL)
    }
    if (!is.null(results_directory)) {
      utils::write.csv(results_df, file.path(results_directory, paste0("results_for_id_", id, "_", gas, ".csv")),
                       row.names = FALSE)
    }
    results_df
  }

  plot_series = function(df, fitted, id) {
    if (!plot_diagnostics) return(invisible(NULL))
    if (!"event" %in% names(df)) df$event = FALSE
    df$fitted = as.numeric(fitted)
    p = ggplot2::ggplot(df, ggplot2::aes(x = time, y = concentration)) +
      ggplot2::geom_point(ggplot2::aes(color = factor(event))) +
      ggplot2::geom_line(ggplot2::aes(y = fitted), color = "blue") +
      ggplot2::scale_color_manual(values = c("FALSE" = "green", "TRUE" = "orange"), name = "Event") +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::labs(title = paste("ID", id, "-", gas, "flux diagnostics"), x = "Time (s)", y = paste0(gas, " (mol)"))
    ggplot2::ggsave(file.path(save_directory, paste0("ID_", id, "_", gas, "_Time_Series_Plot.png")),
                    plot = p, width = 10, height = 5)
  }

  plot_combined = function(combined_data, id) {
    if (is.null(save_directory) || nrow(combined_data) == 0) return(NULL)
    combined_data = combined_data[order(combined_data$time), ]
    p = ggplot2::ggplot(combined_data, ggplot2::aes(x = time, y = concentration, color = process_type)) +
      ggplot2::geom_point() +
      ggplot2::geom_text(ggplot2::aes(label = sequence_id), vjust = -1, check_overlap = TRUE) +
      ggplot2::scale_color_manual(values = c("Diffusive" = "blue", "Ebullition" = "red")) +
      ggplot2::labs(title = paste("ID", id, "-", gas, "classified sequences"), x = "Time (s)", y = paste0(gas, " (mol)")) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 30))
    plot_file = file.path(save_directory, paste0("ID_", id, "_", gas, "_Combined_Time_Series_Plot.png"))
    ggplot2::ggsave(plot_file, plot = p, width = 10, height = 5)
    plot_file
  }

  fit_gam_or_lm = function(seq_data) {
    tryCatch(mgcv::gam(concentration ~ s(time, k = gam_knots), data = seq_data),
             error = function(e) {
               warning("GAM fitting failed for a sequence (", e$message, ") - using a linear model.")
               stats::lm(concentration ~ time, data = seq_data)
             })
  }

  flux_by_sequence = function(combined_data, id, decision) {
    results = list(); plot_data = list()
    area_A = combined_data$area[1]
    for (seq_id in unique(combined_data$sequence_id)) {
      seq_data = combined_data[combined_data$sequence_id == seq_id, ]
      if (nrow(seq_data) <= 1 || sum(!is.na(seq_data$concentration)) <= 1) next
      if (nrow(seq_data) < gam_knots * 3) {
        warning("Not enough points to fit sequence ", seq_id, " of ID ", id, " - skipped.")
        next
      }
      time_diff = diff(seq_data$time)
      if (any(time_diff <= 0)) {
        warning("Time values are not strictly increasing in sequence ", seq_id, " of ID ", id, " - skipped.")
        next
      }
      model = fit_gam_or_lm(seq_data)
      predicted = as.numeric(stats::predict(model, newdata = seq_data))
      fit_type = if (inherits(model, "gam")) "GAM" else "LM"
      r_squared = if (fit_type == "GAM") summary(model)$r.sq else summary(model)$r.squared
      flux = calculate_flux(diff(predicted) / time_diff, area_A)
      results[[length(results) + 1]] = result_row(paste(id, seq_id, sep = "_"), id, seq_data$process_type[1],
                                                  mean(flux, na.rm = TRUE), r_squared,
                                                  sqrt(mean((seq_data$concentration - predicted)^2, na.rm = TRUE)),
                                                  fit_type, decision)
      plot_data[[as.character(seq_id)]] = data.frame(time = seq_data$time, concentration = seq_data$concentration,
                                                     predicted = predicted, sequence_id = seq_id, fit_type = fit_type)
    }
    if (plot_diagnostics && length(plot_data) > 0) {
      pd = do.call(rbind, plot_data)
      p = ggplot2::ggplot(pd, ggplot2::aes(x = time, y = concentration)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(ggplot2::aes(y = predicted, color = fit_type)) +
        ggplot2::facet_wrap(~ sequence_id, scales = "free") +
        ggplot2::scale_color_manual(values = c("GAM" = "red", "LM" = "blue")) +
        ggplot2::labs(title = "Flux data and model fit for each sequence", x = "Time (s)", y = paste0(gas, " (mol)"))
      ggplot2::ggsave(file.path(save_directory, paste0(id, "_", gas, "_Flux_Plot.png")), plot = p, width = 10, height = 5)
    }
    if (length(results) == 0) return(NULL)
    do.call(rbind, results)
  }

  manual_sequences = function(subset_data, id) {
    # Returns a results data.frame, or NULL if no manual sequences were entered.
    if (!interactive()) {
      message("ID ", id, ": manual sequence entry needs an interactive session - no manual sequences used.")
      return(NULL)
    }
    if (ask("Would you like to manually enter sequences? (yes/no): ", "no") != "yes") {
      say("Proceeding without changes.")
      return(NULL)
    }
    subset_data$process_type = NA_character_
    subset_data$sequence_id = NA_real_
    cat("Available time range:", min(subset_data$time), "to", max(subset_data$time), "\n")
    cat("Enter ranges like: 100-150,D,1;200-250,E,2 (type 'done' to finish)\n")
    repeat {
      input = readline()
      if (tolower(trimws(input)) == "done") break
      for (info in unlist(strsplit(gsub("\\s+", "", input), ";"))) {
        parts = unlist(strsplit(info, ","))
        if (length(parts) != 3) { cat("Invalid input format:", info, "\n"); next }
        range_vals = suppressWarnings(as.numeric(unlist(strsplit(parts[1], "-"))))
        seq_type = tolower(parts[2]); seq_id = suppressWarnings(as.numeric(parts[3]))
        if (length(range_vals) != 2 || anyNA(range_vals) || !seq_type %in% c("d", "e") || is.na(seq_id)) {
          cat("Invalid input in:", info, "\n"); next
        }
        idx = subset_data$time >= range_vals[1] & subset_data$time <= range_vals[2]
        if (!any(idx)) { cat("No data found in range:", parts[1], "\n"); next }
        subset_data$process_type[idx] = ifelse(seq_type == "d", "Diffusive", "Ebullition")
        subset_data$sequence_id[idx] = seq_id
      }
    }
    valid = subset_data[!is.na(subset_data$process_type) & !is.na(subset_data$sequence_id), ]
    if (nrow(valid) == 0) { cat("No valid sequences identified.\n"); return(NULL) }
    flux_by_sequence(valid, id, "manual")
  }

  whole_deployment_result = function(subset_data, predicted, r_squared, rmse, id, decision) {
    flux = calculate_flux(diff(predicted) / diff(subset_data$time), subset_data$area[1])
    result_row(id, id, "Diffusive", mean(flux, na.rm = TRUE), r_squared, rmse, "GAM", decision)
  }

  # ---- Deployments to process ----
  data = data[order(data$id, data$time), , drop = FALSE]
  ids = unique(data$id)
  if (!is.null(start_id)) {
    pos = match(as.character(start_id), as.character(ids))
    if (is.na(pos)) stop("Invalid start_id entered: ", start_id)
    ids = ids[pos:length(ids)]
  }

  all_results = list()
  keep = function(df) if (!is.null(df) && nrow(df) > 0) all_results[[length(all_results) + 1]] <<- df

  # ---- Main loop ----
  for (id in ids) {
    subset_data = data[data$id == id, , drop = FALSE]
    subset_data$air_pressure = convert_pressure_to_Pa(subset_data$air_pressure)

    if (all(is.na(subset_data$temperature_celsius))) {
      warning("All temperature values are missing in ID ", id, " - skipped.")
      next
    }
    if (anyNA(subset_data$temperature_celsius)) {
      warning("Missing temperature values in ID ", id, " - imputed with the deployment mean.")
      subset_data$temperature_celsius[is.na(subset_data$temperature_celsius)] = mean(subset_data$temperature_celsius, na.rm = TRUE)
    }
    if (anyNA(subset_data$air_pressure)) {
      warning("Missing air pressure values in ID ", id, " - imputed with the deployment mean.")
      subset_data$air_pressure[is.na(subset_data$air_pressure)] = mean(subset_data$air_pressure, na.rm = TRUE)
    }
    subset_data$temperature_kelvin = subset_data$temperature_celsius + celsius_to_kelvin
    subset_data$concentration = convert_ppm_to_mol(subset_data$concentration, subset_data$air_pressure,
                                                   subset_data$volume[1], subset_data$temperature_kelvin)

    if (any(diff(subset_data$time) <= 0)) {
      warning("Time values are not strictly increasing in ID ", id, " - skipped.")
      next
    }

    gam_model = tryCatch(mgcv::gam(concentration ~ s(time, k = gam_knots), data = subset_data),
                         error = function(e) { warning("GAM fitting failed for ID ", id, ": ", e$message); NULL })
    if (is.null(gam_model)) {
      keep(write_results(manual_sequences(subset_data, id), id))
      next
    }

    predicted = as.numeric(stats::predict(gam_model, newdata = subset_data))
    rmse_gam = sqrt(mean((subset_data$concentration - predicted)^2, na.rm = TRUE))
    r_squared = summary(gam_model)$r.sq
    if (!is.finite(r_squared)) {
      warning("R-squared is not finite for ID ", id, ".")
      r_squared = -Inf
    }

    # ---- Diffusive-only mode (default for CO2) ----
    if (!enable_ebullition_detection) {
      plot_series(subset_data, predicted, id)
      if (require_manual_review && r_squared < r_squared_diffusive && interactive()) {
        cat("GAM fitted for ID", id, "with R-squared:", signif(r_squared, 4), "RMSE:", signif(rmse_gam, 4), "\n")
        if (ask("Use the GAM diffusive flux (yes) or enter sequences manually (no)? (yes/no): ", "yes") == "no") {
          manual = manual_sequences(subset_data, id)
          if (!is.null(manual)) { keep(write_results(manual, id)); next }
        }
        keep(write_results(whole_deployment_result(subset_data, predicted, r_squared, rmse_gam, id, "reviewed"), id))
        next
      }
      keep(write_results(whole_deployment_result(subset_data, predicted, r_squared, rmse_gam, id, "automatic"), id))
      next
    }

    # ---- Diffusive acceptance before ebullition screening ----
    if (r_squared > r_squared_diffusive) {
      plot_series(subset_data, predicted, id)
      keep(write_results(whole_deployment_result(subset_data, predicted, r_squared, rmse_gam, id, "automatic"), id))
      next
    }

    # ---- Optional smoothing before BRT screening ----
    ws = window_size_concentration_smoother
    if (ws > 1) {
      if (ws <= nrow(subset_data)) {
        sm = stats::filter(subset_data$concentration, rep(1 / ws, ws), sides = 2)
        sm[is.na(sm)] = subset_data$concentration[is.na(sm)]
        subset_data$concentration = as.numeric(sm)
      } else {
        warning("Concentration smoothing window exceeds the number of observations in ID ", id, " - not smoothed.")
      }
    }

    # ---- BRT event detection ----
    if (is.null(predictors)) {
      predictors = setdiff(names(subset_data), c("id", "time", "concentration", "temperature_kelvin", "R"))
    }
    pred_cols = if (is.numeric(predictors)) names(subset_data)[predictors] else predictors
    valid_predictors = pred_cols[!vapply(subset_data[pred_cols], function(x) all(is.na(x)), logical(1))]
    if (length(valid_predictors) == 0) {
      warning("All predictors are missing for ID ", id, " - BRT screening skipped.")
      next
    }

    if (!is.null(seed)) set.seed(seed)
    best_brt = NULL
    for (lr in learning_rate) {
      say("Training BRT with learning rate ", lr, " for ID ", id)
      brt = tryCatch(
        dismo::gbm.step(data = subset_data, gbm.x = valid_predictors,
                        gbm.y = which(colnames(subset_data) == "concentration"),
                        family = "gaussian", tree.complexity = tree_complexity, learning.rate = lr,
                        max.trees = max_trees, bag.fraction = bag_fraction,
                        silent = TRUE, plot.folds = FALSE, keep.fold.fit = TRUE),
        error = function(e) { warning("BRT fitting failed for ID ", id, " (learning rate ", lr, "): ", e$message); NULL })
      if (!is.null(brt)) {
        best_brt = brt
        if (!is.null(brt$n.trees) && brt$n.trees >= min_trees_stop) break
      }
    }
    if (is.null(best_brt)) {
      say("No valid BRT model for ID ", id)
      next
    }

    gbm_predictions = as.numeric(stats::predict(best_brt, newdata = subset_data, n.trees = best_brt$n.trees))
    gbm_residuals = subset_data$concentration - gbm_predictions
    wr = window_size_residuals_smoother
    if (wr > 1) {
      sr = stats::filter(gbm_residuals, rep(1 / wr, wr), sides = 2)
      sr[is.na(sr)] = gbm_residuals[is.na(sr)]
      smoothed_residuals = as.numeric(sr)
    } else {
      smoothed_residuals = gbm_residuals
    }
    dynamic_threshold = dynamic_multiplier * stats::sd(smoothed_residuals, na.rm = TRUE)
    concentration_diff = c(0, abs(diff(subset_data$concentration)))
    jump_threshold = stats::quantile(concentration_diff, quantile_threshold, na.rm = TRUE)
    subset_data$event = (abs(smoothed_residuals) > dynamic_threshold) | (concentration_diff > jump_threshold)
    plot_series(subset_data, gbm_predictions, id)

    # ---- Group flagged points into event sequences ----
    in_event = FALSE; event_start = NA; event_end = NA; buffer_count = 0; last_buffer = NA
    for (i in seq_len(nrow(subset_data))) {
      if (!is.na(subset_data$event[i]) && subset_data$event[i]) {
        if (!in_event) { event_start = i; in_event = TRUE }
        event_end = i; buffer_count = 0
      } else if (in_event) {
        buffer_count = buffer_count + 1; last_buffer = i
        if (buffer_count >= buffer_size) {
          subset_data$event[event_start:event_end] = (event_end - event_start + 1) >= min_ebullition_sequence
          in_event = FALSE; buffer_count = 0
        }
      }
    }
    if (in_event) {
      long_enough = (event_end - event_start + 1) >= min_ebullition_sequence
      # An event still open at the end of the record extends through its trailing buffer (as in v1.0.1).
      end_i = if (long_enough && buffer_count > 0) last_buffer else event_end
      subset_data$event[event_start:end_i] = long_enough
    }
    subset_data$event[is.na(subset_data$event)] = FALSE

    # ---- Review of classified sequences ----
    subset_data$process_type = ifelse(subset_data$event, "Ebullition", "Diffusive")
    combined_data = subset_data[order(subset_data$time), ]
    combined_data$sequence_id = cumsum(c(1, diff(as.numeric(factor(combined_data$process_type))) != 0))

    plot_file = plot_combined(combined_data, id)
    where = if (is.null(plot_file)) "(plots disabled)" else plot_file
    choice = ask(paste0("Review the classified sequences for ID ", id, " ", where,
                        "\nAre the classified sequences correct? (yes/no/reclassify/skip): "), "yes")

    if (choice == "reclassify") {
      repeat {
        input = ask("Enter '<sequence_id> D/E' (e.g., 2 D); type 'done' when finished: ", "done")
        if (input == "done") break
        parts = strsplit(input, "\\s+")[[1]]
        sid = suppressWarnings(as.numeric(parts[1]))
        if (length(parts) != 2 || is.na(sid) || !parts[2] %in% c("d", "e")) {
          cat("Invalid input. Use: <sequence_id> <D/E>\n"); next
        }
        if (!any(combined_data$sequence_id == sid)) { cat("Sequence ID not found.\n"); next }
        combined_data$process_type[combined_data$sequence_id == sid] = ifelse(parts[2] == "d", "Diffusive", "Ebullition")
        cat("Sequence", sid, "reclassified.\n")
      }
      plot_file = plot_combined(combined_data, id)
      choice = if (ask("Are the re-classified sequences correct? (yes/no): ", "yes") == "yes") "reviewed" else "no"
    }

    if (choice %in% c("yes", "reviewed")) {
      keep(write_results(flux_by_sequence(combined_data, id, if (choice == "yes") "automatic" else "reviewed"), id))
    } else if (choice == "skip") {
      say("Skipping ID ", id)
    } else if (choice == "no") {
      keep(write_results(manual_sequences(subset_data, id), id))
    } else {
      say("Invalid input - skipping ID ", id)
    }
  }

  final = if (length(all_results) > 0) do.call(rbind, all_results) else
    result_row(character(0), character(0), character(0), numeric(0), numeric(0), numeric(0), character(0), character(0))
  rownames(final) = NULL
  invisible(final)
}

#' Calculate CH4 flux (backward-compatible wrapper)
#'
#' `calculate_methane_flux()` is the v1.0.x interface. It calls [calculate_gas_flux()]
#' with `gas = "CH4"`; all other arguments are passed through unchanged.
#'
#' @param data,... Passed to [calculate_gas_flux()].
#' @return As [calculate_gas_flux()].
calculate_methane_flux = function(data, ...) {
  calculate_gas_flux(data, ..., gas = "CH4")
}
