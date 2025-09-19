#' Calculate methane flux from chamber time series
#'
#' This function takes a data.frame of chamber measurements and estimates
#' diffusive and ebullitive methane fluxes per sequence. The function is a
#' non-CRAN example intended for laboratory or field processing pipelines.
#'
#' @param data A data.frame containing required columns: `id`, `time`,
#'   `temperature_celsius`, `volume`, `area`, `air_pressure`, `concentration` (ppm).
#' @param results_directory Directory where per-ID results CSV files are written.
#' @param save_directory Directory where diagnostic plots will be saved.
#' @param celsius_to_kelvin Numeric offset to convert Celsius to Kelvin (default 273.15).
#' @param plot_diagnostics Logical, save diagnostic plots when TRUE.
#' @param unit Character, output flux units: one of `"mg/m2/s"` (default),
#'   `"umol/m2/s"`, or `"nmol/m2/s"`.
#' @param start_id Optional ID (character or numeric) to start processing from.
#' @param r_squared_diffusive Threshold R^2 above which a GAM fit is considered
#'   a diffusive-only time series (default 0.99).
#' @param window_size_concentration_smoother Integer smoothing window for raw
#'   concentration before training GBM (default 0: no smoothing).
#' @param window_size_residuals_smoother Integer smoothing window for GBM residuals (default 5).
#' @param quantile_threshold Numeric (0-1) quantile for concentration jumps used in event detection.
#' @param dynamic_multiplier Multiplier applied to residual standard deviation to set dynamic threshold.
#' @param predictors Optional character vector of predictor column names for GBM.
#' @param learning_rate Numeric vector of learning rates to try for GBM training.
#' @param tree_complexity Integer tree depth for GBM.
#' @param bag_fraction Fraction of data used per GBM tree.
#' @param gam_knots Default number of knots to use for GAM smoothing.
#' @param min_ebullition_sequence Minimum sequence length considered an ebullition event.
#' @param R_value Universal gas constant in SI units (default 8.314 J/mol/K).
#' @param pressure_units Units for `air_pressure` column. Supported: `Pa`, `kPa`, `atm`, `bar`, `mbar`, `mmHg`, `psi`.
#'
#' @return Invisibly returns a data.frame of all per-ID results (also writes CSVs).
#' @export
#' @examples
#' 
#' # See inst/examples/run_example.R for a runnable example
#'
calculate_methane_flux = function(data, results_directory = NULL, save_directory = NULL, celsius_to_kelvin = 273.15,
                                  plot_diagnostics = TRUE, unit = "mg/m2/s", start_id = NULL, r_squared_diffusive = 0.99,
                                  window_size_concentration_smoother = 0, window_size_residuals_smoother = 5,
                                  quantile_threshold = 0.95, dynamic_multiplier = 2, predictors = NULL,
                                  learning_rate = c(0.001, 0.002, 0.003, 0.004, 0.005),
                                  tree_complexity = 5, bag_fraction = 0.15, gam_knots = 5, min_ebullition_sequence = 8,
                                  R_value = 8.314, pressure_units = "Pa") {
  
  # Packages required from CRAN/Bioconductor (only non-base packages listed)
  necessary_packages <- c("mgcv", "dismo", "gbm", "ggplot2", "scales")
  # Check which required packages are missing
  missing_pkgs <- necessary_packages[!vapply(necessary_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      ". Install them with: install.packages(c(\"", paste(missing_pkgs, collapse = "\", \""), "\"))"
    )
  }
  
  required_columns = c("id", "time", "temperature_celsius", "volume", "area", "air_pressure", "concentration")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input data must contain columns: ", paste(required_columns, collapse = ", "))
  }
  
  if (!is.null(save_directory) && !dir.exists(save_directory)) {
    stop("Specified save directory does not exist: ", save_directory)
  }
  
  data = data[order(data$id, data$time), , drop = FALSE]
  
  # Validate unit and pressure_units arguments to fail early with clear messages
  unit <- match.arg(unit, c("mg/m2/s", "umol/m2/s", "nmol/m2/s"))
  pressure_units <- match.arg(pressure_units, c("Pa", "kPa", "atm", "bar", "mbar", "mmHg", "psi"))
  
  convert_pressure_to_Pa = function(pressure_value, pressure_units) {
    unit_conversions = list(
      Pa = 1,
      kPa = 1000,
      atm = 101325,
      bar = 100000,
      mbar = 100,
      mmHg = 133.322,
      psi = 6894.76
    )
    if (!pressure_units %in% names(unit_conversions)) stop("Unsupported pressure units: ", pressure_units)
    return(pressure_value * unit_conversions[[pressure_units]])
  }
  
  convert_ppm_to_mol = function(methane_ppm, pressure_P, volume_V, R, temperature_kelvin) {
    methane_fraction = methane_ppm / 1e6
    total_moles_air = (pressure_P * volume_V) / (R * temperature_kelvin)
    concentration_methane_mol = methane_fraction * total_moles_air
    return(concentration_methane_mol)
  }
  
  calculate_flux = function(slope, area_A, molecular_weight_CH4, unit) {
    flux_mol_m2_s = slope / area_A
    if (unit == "nmol/m2/s") {
      flux = flux_mol_m2_s * 1e9
    } else if (unit == "umol/m2/s") {
      flux = flux_mol_m2_s * 1e6
    } else if (unit == "mg/m2/s") {
      flux = flux_mol_m2_s * molecular_weight_CH4 * 1000
    } else {
      stop("Unsupported unit for flux calculation: ", unit)
    }
    return(flux)
  }
  
  plot_diagnostics_func = function(data, model, id, save_directory) {
    if (!"event" %in% names(data)) data$event = FALSE
    p = ggplot2::ggplot(data, ggplot2::aes(x = time, y = concentration)) +
      ggplot2::geom_point(ggplot2::aes(color = factor(event))) +
      ggplot2::geom_line(ggplot2::aes(y = predict(model, newdata = data)), color = "blue") +
      ggplot2::scale_color_manual(values = c("FALSE" = "green", "TRUE" = "orange"), name = "Event Type") +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::labs(title = paste("ID", id, "- Methane Flux Diagnostics"), x = "Time", y = "CH4-mol")
    
    if (!is.null(save_directory)) {
      ggplot2::ggsave(file.path(save_directory, paste0("ID_", id, "_Time_Series_Plot.png")), plot = p, width = 10, height = 5)
    }
  }
  
  plot_combined_diagnostics_func = function(combined_data, id, save_directory) {
    if (nrow(combined_data) == 0) {
      warning("No data available for plotting.")
      return(NULL)
    }
    combined_data = combined_data[order(combined_data$time), ]
    p = ggplot2::ggplot(combined_data, ggplot2::aes(x = time, y = concentration, color = process_type)) +
      ggplot2::geom_point() +
      ggplot2::geom_text(ggplot2::aes(label = sequence_id), vjust = -1, check_overlap = TRUE) +
      ggplot2::scale_color_manual(values = c("Diffusive" = "blue", "Ebullition" = "red")) +
      ggplot2::labs(title = paste("ID", id, "- Combined Methane Flux Diagnostics"), x = "Time", y = "CH4-mol") +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 30))
    
    if (!is.null(save_directory)) {
      plot_file_name = file.path(save_directory, paste0("ID_", id, "_Combined_Time_Series_Plot.png"))
      ggplot2::ggsave(plot_file_name, plot = p, width = 10, height = 5)
      return(plot_file_name)
    }
    return(NULL)
  }
  
  calculate_and_plot_flux_for_each_sequence = function(combined_data, area_A, molecular_weight_CH4, unit, save_directory, id_prefix) {
    results = list()
    plot_data_list = list()
    
    for (seq_id in unique(combined_data$sequence_id)) {
      seq_data = combined_data[combined_data$sequence_id == seq_id, ]
      if (nrow(seq_data) > 1 && sum(!is.na(seq_data$concentration)) > 1) {
        unique_seq_id = paste(id_prefix, seq_id, sep = "_")
        if (nrow(seq_data) >= gam_knots * 3) {
          gam_model = tryCatch({
            mgcv::gam(concentration ~ s(time, k = gam_knots), data = seq_data)
          }, error = function(e) {
            warning(sprintf("GAM fitting failed for sequence %s: %s - switching to linear model.", seq_id, e$message))
            lm(concentration ~ time, data = seq_data)
          })
          
          model_type = ifelse(class(gam_model)[1] == "gam", "GAM", "LM")
          predicted_values = predict(gam_model, newdata = seq_data)
          residuals = seq_data$concentration - predicted_values
          r_squared = ifelse(model_type == "GAM", summary(gam_model)$r.sq, summary(gam_model)$r.squared)
          
          time_diff = diff(seq_data$time)
          if (any(time_diff <= 0)) {
            warning("Time values are not strictly increasing in sequence ID ", seq_id)
            next
          }
          
          slope = diff(predicted_values) / time_diff
          flux = calculate_flux(slope, area_A, molecular_weight_CH4, unit)
          
          results[[length(results) + 1]] = data.frame(
            ID = unique_seq_id,
            Flux_Type = seq_data$process_type[1],
            Flux_Estimate = mean(flux, na.rm = TRUE),
            R_squared = r_squared,
            RMSE = sqrt(mean(residuals^2)),
            Fit_Type = model_type,
            stringsAsFactors = FALSE
          )
          
          plot_df = data.frame(time = seq_data$time, concentration = seq_data$concentration,
                               predicted = predicted_values, sequence_id = seq_id, fit_type = model_type)
          plot_data_list[[as.character(seq_id)]] = plot_df
        } else {
          warning(sprintf("Not enough data points to fit GAM for sequence ID %s - skipping.", seq_id))
          next
        }
      }
    }
    
    if (length(plot_data_list) > 0) {
      combined_plot_data = do.call(rbind, plot_data_list)
      if (!is.null(save_directory)) {
        plot_file_name = file.path(save_directory, paste0(id_prefix, "_Flux_Plot.png"))
        p = ggplot2::ggplot(combined_plot_data, ggplot2::aes(x = time, y = concentration)) +
          ggplot2::geom_point() +
          ggplot2::geom_line(ggplot2::aes(y = predicted, color = fit_type)) +
          ggplot2::facet_wrap(~ sequence_id, scales = "free") +
          ggplot2::scale_color_manual(values = c("GAM" = "red", "LM" = "blue")) +
          ggplot2::ggtitle("Flux Data and Model Fit for Each Sequence") +
          ggplot2::xlab("Time") + ggplot2::ylab("CH4-mol")
        ggplot2::ggsave(plot_file_name, plot = p, width = 10, height = 5)
      }
    } else {
      combined_plot_data = NULL
    }
    
    if (length(results) == 0) return(data.frame())
    return(do.call(rbind, results))
  }
  
  if (!is.null(start_id) && !(start_id %in% unique(data$id))) stop("Invalid start_id entered.")
  if (is.null(start_id)) start_id = as.character(min(unique(data$id)))
  valid_ids = unique(data$id)[unique(data$id) >= start_id]
  
  all_results = list()
  
  for (id in valid_ids) {
    subset_data = data[data$id == id, , drop = FALSE]
    diffusive_results = list()
    ebullition_results = list()
    molecular_weight_CH4 = 16.04
    
    subset_data$air_pressure = convert_pressure_to_Pa(subset_data$air_pressure, pressure_units)
    
    if (all(is.na(subset_data$temperature_celsius))) {
      warning(sprintf("All temperature values are missing in ID %s - skipping to next ID.", id))
      next
    } else {
      missing_temps = is.na(subset_data$temperature_celsius)
      if (any(missing_temps)) {
        warning(sprintf("Missing temperature values in ID %s - imputing with mean value.", id))
        subset_data$temperature_celsius[missing_temps] = mean(subset_data$temperature_celsius, na.rm = TRUE)
      }
    }
    
    missing_airpressure = is.na(subset_data$air_pressure)
    if (any(missing_airpressure)) {
      warning(sprintf("Missing air pressure values in ID %s - imputing with mean value.", id))
      subset_data$air_pressure[missing_airpressure] = mean(subset_data$air_pressure, na.rm = TRUE)
    }
    
    subset_data$temperature_kelvin = subset_data$temperature_celsius + celsius_to_kelvin
    
    volume_V = subset_data$volume[1]
    area_A = subset_data$area[1]
    
    subset_data$concentration = convert_ppm_to_mol(
      subset_data$concentration,
      subset_data$air_pressure,
      volume_V,
      R_value,
      subset_data$temperature_kelvin
    )
    
    time_diff = diff(subset_data$time)
    if (any(time_diff <= 0)) {
      warning(sprintf("Time values are not strictly increasing in ID %s - skipping to next ID.", id))
      next
    }
    
    gam_model = tryCatch({
      mgcv::gam(concentration ~ s(time, k = gam_knots), data = subset_data)
    }, error = function(e) {
      warning(sprintf("GAM fitting failed for ID %s: %s - skipping this ID.", id, e$message))
      return(NULL)
    })
    
    if (is.null(gam_model)) next
    
    predicted_values = predict(gam_model, newdata = subset_data)
    residuals = subset_data$concentration - predicted_values
    r_squared = summary(gam_model)$r.sq
    rmse_gam = sqrt(mean(residuals^2))
    
    if (r_squared > r_squared_diffusive) {
      if (plot_diagnostics) plot_diagnostics_func(subset_data, gam_model, id, save_directory)
      slope = diff(predict(gam_model, newdata = subset_data)) / diff(subset_data$time)
      flux = calculate_flux(slope, area_A, molecular_weight_CH4, unit)
      
      diffusive_results[[length(diffusive_results) + 1]] = data.frame(
        ID = id,
        Flux_Type = "Diffusive",
        Flux_Estimate = mean(flux, na.rm = TRUE),
        R_squared = r_squared,
        RMSE = rmse_gam,
        stringsAsFactors = FALSE
      )
      
      results_df = do.call(rbind, diffusive_results)
      if (!is.null(results_directory) && nrow(results_df) > 0) {
        results_file_path = file.path(results_directory, paste0("results_for_id_", id, ".csv"))
        utils::write.csv(results_df, results_file_path, row.names = FALSE)
      }
      all_results[[length(all_results) + 1]] = results_df
      next
    }
    
    window_size = window_size_concentration_smoother
    if (window_size > 1) {
      if (window_size <= nrow(subset_data)) {
        smoothed_concentration = stats::filter(subset_data$concentration, rep(1/window_size, window_size), sides = 2)
        smoothed_concentration[is.na(smoothed_concentration)] = subset_data$concentration[is.na(smoothed_concentration)]
        subset_data$concentration = smoothed_concentration
      } else {
        warning(sprintf("Window size for concentration smoothing is larger than the number of observations in ID %s - skipping smoothing.", id))
      }
    }
    
    if (is.null(predictors)) {
      predictors = setdiff(names(subset_data), c("id", "time", "concentration", "temperature_kelvin", "R"))
    }
    
    valid_predictors = predictors[!sapply(subset_data[predictors], function(x) all(is.na(x)))]
    if (length(valid_predictors) == 0) {
      warning(sprintf("All predictors have missing values for ID: %s. Skipping GBM training.", id))
      next
    }
    
    BRT = NULL
    for (lr in learning_rate) {
      message(sprintf("Training with learning rate: %s for ID: %s", lr, id))
      BRT = tryCatch({
        dismo::gbm.step(data = subset_data,
                        gbm.x = valid_predictors,
                        gbm.y = which(colnames(subset_data) == "concentration"),
                        family = "gaussian",
                        tree.complexity = tree_complexity,
                        learning.rate = lr,
                        max.trees = 10000,
                        bag.fraction = bag_fraction,
                        silent = TRUE,
                        plot.folds = FALSE,
                        keep.fold.fit = TRUE)
      }, error = function(e) {
        warning(sprintf("GBM fitting failed for ID %s with learning rate %s: %s", id, lr, e$message))
        return(NULL)
      })
      
      if (!is.null(BRT)) {
        gbm_predictions = predict(BRT, subset_data, n.trees = BRT$n.trees)
        gbm_residuals = subset_data$concentration - gbm_predictions
        
        window_size_residuals = window_size_residuals_smoother
        if (window_size_residuals > 1) {
          smoothed_residuals = stats::filter(gbm_residuals, rep(1/window_size_residuals, window_size_residuals), sides = 2)
          smoothed_residuals[is.na(smoothed_residuals)] = gbm_residuals[is.na(smoothed_residuals)]
        } else {
          smoothed_residuals = gbm_residuals
        }
        
        std_dev_residuals = stats::sd(smoothed_residuals, na.rm = TRUE)
        dynamic_threshold = dynamic_multiplier * std_dev_residuals
        concentration_diff = c(0, abs(diff(subset_data$concentration)))
        threshold_concentration_diff = stats::quantile(concentration_diff, quantile_threshold, na.rm = TRUE)
        
        subset_data$event = (abs(smoothed_residuals) > dynamic_threshold) | (concentration_diff > threshold_concentration_diff)
        if (plot_diagnostics) plot_diagnostics_func(subset_data, BRT, id, save_directory)
        
        if (!is.null(BRT$n.trees) && BRT$n.trees >= 1000) {
          message(sprintf("Stopping criterion met for learning rate: %s in ID: %s", lr, id))
          break
        }
      }
    }
    
    if (is.null(BRT)) {
      message(sprintf("No valid model trained for ID: %s", id))
      next
    }
    
    in_event = FALSE
    event_start = NULL
    event_end = NULL
    diffusive_buffer = vector("list", 5)
    buffer_count = 0
    
    for (i in seq_len(nrow(subset_data))) {
      if (!is.na(subset_data$event[i]) && subset_data$event[i]) {
        if (!in_event) {
          event_start = i
          in_event = TRUE
        }
        event_end = i
        buffer_count = 0
      } else {
        if (in_event) {
          buffer_count = buffer_count + 1
          diffusive_buffer[[buffer_count]] = i
          if (buffer_count == 5) {
            if ((event_end - event_start + 1) >= min_ebullition_sequence) {
              subset_data$event[event_start:event_end] = TRUE
            } else {
              subset_data$event[event_start:event_end] = FALSE
            }
            in_event = FALSE
            buffer_count = 0
          }
        }
      }
    }
    
    if (in_event) {
      if (buffer_count > 0 && buffer_count <= 5 && (event_end - event_start + 1) >= min_ebullition_sequence) {
        subset_data$event[event_start:max(unlist(diffusive_buffer))] = TRUE
      } else if ((event_end - event_start + 1) >= min_ebullition_sequence) {
        subset_data$event[event_start:event_end] = TRUE
      } else {
        subset_data$event[event_start:event_end] = FALSE
      }
    }
    
    subset_data$process_type = ifelse(subset_data$event, "Ebullition", "Diffusive")
    combined_data = subset_data
    combined_data = combined_data[order(combined_data$time), ]
    combined_data$sequence_id = cumsum(c(1, diff(as.numeric(factor(combined_data$process_type))) != 0))
    
    plot_file_name = plot_combined_diagnostics_func(combined_data, id, save_directory)
    
    # Interactive behavior: prompt only when in an interactive session.
    # In non-interactive sessions (e.g., tests/CI) accept classifications by default.
    if (interactive()) {
      cat("Review the combined plot saved at", plot_file_name, "\n")
      cat("Are the classified sequences correct? If not, you can reclassify a sequence, skip this ID, or choose to reclassify (yes/no/reclassify/skip): ")
      initial_input = tolower(readline())
    } else {
      initial_input = "yes"
    }
    
    skip_current_id = FALSE
    
    if (initial_input == "yes") {
      area_A_to_pass = combined_data$area[1]
      sequence_flux_results = calculate_and_plot_flux_for_each_sequence(combined_data, area_A_to_pass, molecular_weight_CH4, unit, save_directory, id)
      diffusive_results = subset(sequence_flux_results, Flux_Type == "Diffusive")
      ebullition_results = subset(sequence_flux_results, Flux_Type == "Ebullition")
      results_df = if (nrow(sequence_flux_results) > 0) sequence_flux_results else data.frame()
      if (!is.null(results_directory) && nrow(results_df) > 0) {
        results_file_path = file.path(results_directory, paste0("results_for_id_", id, ".csv"))
        utils::write.csv(results_df, results_file_path, row.names = FALSE)
      }
      all_results[[length(all_results) + 1]] = results_df
      next
    } else if (initial_input == "reclassify") {
      reclassify_loop = TRUE
      confirmation_input = "no"
      while (reclassify_loop) {
        cat("Enter the sequence ID and the classification to reclassify as diffusive (D) or ebullition (E) (e.g., 2 D), type 'done' when finished: ")
        while (TRUE) {
          input = tolower(readline())
          if (input == "done") break
          sequence_details = strsplit(input, " ")[[1]]
          if (length(sequence_details) == 2) {
            sequence_id = as.numeric(sequence_details[1])
            seq_type = sequence_details[2]
            if (!is.na(sequence_id) && seq_type %in% c("d", "e")) {
              process_type = ifelse(seq_type == "d", "Diffusive", "Ebullition")
              if (nrow(combined_data[combined_data$sequence_id == sequence_id, ]) > 0) {
                combined_data[combined_data$sequence_id == sequence_id, "process_type"] = process_type
                cat("Sequence ID", sequence_id, "has been reclassified as", process_type, ".\n")
              } else {
                cat("Sequence ID", sequence_id, "not found.\n")
              }
            } else {
              cat("Invalid classification. Please enter 'D' for diffusive or 'E' for ebullition.\n")
            }
          } else {
            cat("Invalid input. Please enter a sequence ID followed by 'D' or 'E'.\n")
          }
        }
        plot_file_name = plot_combined_diagnostics_func(combined_data, id, save_directory)
        cat("Review the updated combined plot saved at", plot_file_name, "\nAre the re-classified sequences correct? (yes/no): ")
        confirmation_input = tolower(readline())
        if (confirmation_input == "yes") { reclassify_loop = FALSE; break }
        cat("Continuing reclassification.\n")
      }
      
      if (confirmation_input == "yes") {
        area_A_to_pass = combined_data$area[1]
        sequence_flux_results = calculate_and_plot_flux_for_each_sequence(combined_data, area_A_to_pass, molecular_weight_CH4, unit, save_directory, id)
        results_df = if (nrow(sequence_flux_results) > 0) sequence_flux_results else data.frame()
        if (!is.null(results_directory) && nrow(results_df) > 0) {
          results_file_path = file.path(results_directory, paste0("results_for_id_", id, ".csv"))
          utils::write.csv(results_df, results_file_path, row.names = FALSE)
        }
        all_results[[length(all_results) + 1]] = results_df
        next
      }
      
    } else if (initial_input == "skip") {
      message(sprintf("Skipping ID %s", id))
      next
    } else if (initial_input == "no") {
      cat("Would you like to manually enter sequences? (yes/no): ")
      manual_input = if (interactive()) tolower(readline()) else "no"
      if (manual_input == "yes") {
        if (!"process_type" %in% names(subset_data)) subset_data$process_type = NA
        if (!"sequence_id" %in% names(subset_data)) subset_data$sequence_id = NA
        cat("Available time range:", min(subset_data$time), "to", max(subset_data$time), "\n")
        cat("Enter the range, classification, and numerical label for any sequence (e.g., 100-150,D,1;200-250,E,2), type 'done' when finished: ")
        repeat {
          input = readline()
          if (tolower(input) == "done") break
          input = gsub("\\s+", "", input)
          sequence_info = unlist(strsplit(input, ";"))
          for (info in sequence_info) {
            sequence_details = unlist(strsplit(info, ","))
            if (length(sequence_details) != 3) { cat("Invalid input format. Please follow the example.\n"); next }
            sequence_range = as.numeric(unlist(strsplit(sequence_details[1], "-")))
            sequence_type = tolower(sequence_details[2])
            sequence_id = as.numeric(sequence_details[3])
            if (length(sequence_range) != 2 || !sequence_type %in% c("d", "e") || is.na(sequence_id)) {
              cat("Invalid input. Make sure to enter a valid range followed by 'D' or 'E' for classification, and a numerical label.\n")
            } else {
              subset_index = subset_data$time >= sequence_range[1] & subset_data$time <= sequence_range[2]
              if (!any(subset_index)) { cat("No data found in the specified range.\n") } else {
                subset_data$process_type[subset_index] = ifelse(sequence_type == "d", "Diffusive", "Ebullition")
                subset_data$sequence_id[subset_index] = sequence_id
              }
            }
          }
        }
        
        if (!"event" %in% names(subset_data)) subset_data$event = FALSE
        else subset_data$event[is.na(subset_data$event)] = FALSE
        valid_data = subset_data[!is.na(subset_data$process_type) & !is.na(subset_data$sequence_id), ]
        if (nrow(valid_data) == 0) { cat("No valid sequences identified. Skipping to next ID.\n"); next }
        area_A_to_pass = valid_data$area[1]
        sequence_flux_results = calculate_and_plot_flux_for_each_sequence(valid_data, area_A_to_pass, molecular_weight_CH4, unit, save_directory, id)
        results_df = if (nrow(sequence_flux_results) > 0) sequence_flux_results else data.frame()
        if (!is.null(results_directory) && nrow(results_df) > 0) {
          results_file_path = file.path(results_directory, paste0("results_for_id_", id, ".csv"))
          utils::write.csv(results_df, results_file_path, row.names = FALSE)
        }
        all_results[[length(all_results) + 1]] = results_df
        next
      } else {
        message("Proceeding without changes.")
        next
      }
      
    } else {
      message("Invalid input. Skipping ID by default.")
      next
    }
  }
  
  final_results = do.call(rbind, all_results)
  invisible(final_results)
}


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
                       # start_id = "str5",
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
 