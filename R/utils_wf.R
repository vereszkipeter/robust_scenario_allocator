library(lubridate) # For date manipulation
library(dplyr)    # For data manipulation
library(tibble)   # For tibble creation

#' @title Generate Walk-Forward Validation Windows
#' @description Generates a sequence of training windows (from_date, to_date)
#' and corresponding validation dates for walk-forward analysis.
#' @param app_config A list containing application configuration, including
#'   `data$from`, `data$to`, `initial_window_months`, `roll_forward_months`,
#'   and `validation_start_date`.
#' @return A tibble with columns: `window_id`, `from_date`, `to_date`, `val_date`.
generate_wf_windows <- function(app_config, panel_data) {
  
  overall_from_date <- start(panel_data)
  overall_to_date <- end(panel_data)
  
  initial_window_months <- app_config$default$initial_window_months
  roll_forward_months <- app_config$default$roll_forward_months
  validation_start_date <- as.Date(app_config$default$validation_start_date)
  min_dcc_obs <- app_config$default$models$min_dcc_obs # Read min_dcc_obs


  # Ensure validation_start_date is after the initial window can be formed
  message("DEBUG: generate_wf_windows - overall_from_date: ", overall_from_date, ", overall_to_date: ", overall_to_date)
  message("DEBUG: generate_wf_windows - initial_window_months: ", initial_window_months, ", roll_forward_months: ", roll_forward_months)
  message("DEBUG: generate_wf_windows - validation_start_date: ", validation_start_date, ", min_dcc_obs: ", min_dcc_obs)
  
  if (validation_start_date < (overall_from_date %m+% months(initial_window_months))) {
    stop("validation_start_date is too early. It must allow for the initial window.")
  }
  
  windows <- tibble()
  current_to_date <- lubridate::ceiling_date(validation_start_date, "month") - days(1)
  window_id_counter <- 1
  
  while (current_to_date <= overall_to_date) {
    message("DEBUG: Loop iteration ", window_id_counter, ", current_to_date: ", current_to_date)
    # Ensure training_window_end is end of month, based on current_to_date (which will also be month-end)
    training_window_end <- current_to_date # current_to_date is already month-end
    training_window_start <- training_window_end %m-% months(initial_window_months)
    
    # Adjust training_window_start to the end of the month if it falls mid-month due to initial_window_months subtraction
    training_window_start <- lubridate::ceiling_date(training_window_start, "month") - days(1)

    # Ensure overall_from_date is also month-end aligned before comparison and assignment
    # Ensure overall_from_date_aligned correctly represents the month-end of the *first* data month.
    # E.g., if overall_from_date is "2000-01-01", overall_from_date_aligned should be "2000-01-31".
    overall_from_date_aligned <- lubridate::floor_date(overall_from_date, "month") %m+% months(1) - days(1)
    if (training_window_start < overall_from_date_aligned) {
        training_window_start <- overall_from_date_aligned
    }

    oos_from_date <- training_window_end + days(1)
    oos_to_date_candidate <- training_window_end %m+% months(roll_forward_months)
    oos_to_date_candidate <- lubridate::ceiling_date(oos_to_date_candidate, "month") - days(1) # Ensure month-end
    oos_to_date <- min(oos_to_date_candidate, lubridate::ceiling_date(overall_to_date, "month") - days(1)) # Ensure overall_to_date is also month-end

    actual_training_months <- ((lubridate::year(training_window_end) - lubridate::year(training_window_start)) * 12) + 
                              lubridate::month(training_window_end) - lubridate::month(training_window_start) + 1 
    
    message("DEBUG:   training_window_start: ", training_window_start, ", training_window_end: ", training_window_end)
    message("DEBUG:   actual_training_months: ", actual_training_months, ", min_dcc_obs: ", min_dcc_obs)
    message("DEBUG:   oos_from_date: ", oos_from_date, ", oos_to_date: ", oos_to_date)
    message("DEBUG:   Condition check: oos_from_date (", oos_from_date, ") <= oos_to_date (", oos_to_date, ") is ", oos_from_date <= oos_to_date)
    message("DEBUG:   Condition check: actual_training_months (", actual_training_months, ") >= min_dcc_obs (", min_dcc_obs, ") is ", actual_training_months >= min_dcc_obs)
    
    if (oos_from_date <= oos_to_date && actual_training_months >= min_dcc_obs) { 
      message("DEBUG:     Adding window ", window_id_counter)
      windows <- bind_rows(windows, tibble(
        window_id = window_id_counter,
        from_date = training_window_start,
        to_date = training_window_end,
        val_date = training_window_end, 
        oos_from_date = oos_from_date,
        oos_to_date = oos_to_date
      ))
    } else {
      message("DEBUG:     Skipping window ", window_id_counter, " due to failed condition.")
    }
    
    # Advance current_to_date, ensuring it remains a month-end date
    current_to_date <- lubridate::ceiling_date(current_to_date %m+% months(roll_forward_months), "month") - days(1)
    window_id_counter <- window_id_counter + 1
  }
  
  message("DEBUG: Generated ", nrow(windows), " walk-forward validation windows. (Minimum ", min_dcc_obs, " months training data required).")
  return(windows)
}

#' @title Plot Time Series of Returns
#' @description Saves a plot of monthly returns time series.
#' @param returns_df A data.frame with Date, Asset, and Returns columns.
#' @param output_plots_dir Directory to save the plot.
plot_returns_time_series <- function(returns_df, output_plots_dir) {
  p <- ggplot(returns_df, aes(x = Date, y = Returns, color = Asset)) +
    geom_line() +
    labs(title = "Monthly Returns Time Series", y = "Returns") +
    theme_minimal()
  ggsave(file.path(output_plots_dir, "returns_time_series.png"), plot = p, width = 10, height = 6)
}

#' @title Plot Returns Distributions
#' @description Saves a plot of monthly returns distributions.
#' @param returns_df A data.frame with Date, Asset, and Returns columns.
#' @param output_plots_dir Directory to save the plot.
plot_returns_distributions <- function(returns_df, output_plots_dir) {
  p <- ggplot(returns_df, aes(x = Returns, fill = Asset)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    geom_density(alpha = 0.2, aes(y = after_stat(density))) +
    labs(title = "Monthly Returns Distributions", x = "Returns", y = "Density/Count") +
    theme_minimal()
  ggsave(file.path(output_plots_dir, "returns_distributions.png"), plot = p, width = 10, height = 6)
}

#' @title Plot Returns Boxplots
#' @description Saves a plot of monthly returns boxplots.
#' @param returns_df A data.frame with Date, Asset, and Returns columns.
#' @param output_plots_dir Directory to save the plot.
plot_returns_boxplots <- function(returns_df, output_plots_dir) {
  p <- ggplot(returns_df, aes(x = Asset, y = Returns, fill = Asset)) +
    geom_boxplot() +
    labs(title = "Monthly Returns Boxplots", x = "Asset", y = "Returns") +
    theme_minimal()
  ggsave(file.path(output_plots_dir, "returns_boxplots.png"), plot = p, width = 8, height = 5)
}

#' @title Plot ACF of Returns
#' @description Saves a combined plot of ACF for each asset's returns.
#' @param returns_df A data.frame with Date, Asset, and Returns columns.
#' @param output_plots_dir Directory to save the plot.
plot_returns_acf <- function(returns_df, output_plots_dir) {
  acf_plots <- list()
  for (asset in unique(returns_df$Asset)) {
    asset_returns <- returns_df %>% filter(Asset == asset) %>% pull(Returns)
    p_acf <- forecast::Acf(asset_returns, plot = FALSE) %>%
      autoplot() +
      labs(title = paste("ACF of", asset), y = "ACF") +
      theme_minimal()
    acf_plots[[asset]] <- p_acf
  }
  combined_acf_plots <- wrap_plots(acf_plots, ncol = 2) +
    plot_annotation(title = "Autocorrelation Functions of Monthly Returns")
  ggsave(file.path(output_plots_dir, "returns_acf_plots.png"), plot = combined_acf_plots, width = 12, height = 8)
}

#' @title Plot PACF of Returns
#' @description Saves a combined plot of PACF for each asset's returns.
#' @param returns_df A data.frame with Date, Asset, and Returns columns.
#' @param output_plots_dir Directory to save the plot.
plot_returns_pacf <- function(returns_df, output_plots_dir) {
  pacf_plots <- list()
  for (asset in unique(returns_df$Asset)) {
    asset_returns <- returns_df %>% filter(Asset == asset) %>% pull(Returns)
    p_pacf <- forecast::Pacf(asset_returns, plot = FALSE) %>%
      autoplot() +
      labs(title = paste("PACF of", asset), y = "PACF") +
      theme_minimal()
    pacf_plots[[asset]] <- p_pacf
  }
  combined_pacf_plots <- wrap_plots(pacf_plots, ncol = 2) +
    plot_annotation(title = "Partial Autocorrelation Functions of Monthly Returns")
  ggsave(file.path(output_plots_dir, "returns_pacf_plots.png"), plot = combined_pacf_plots, width = 12, height = 8)
}

#' @title Safely Process a Single Walk-Forward Window
#' @description Wraps the `process_window` function in a `tryCatch` block
#' to ensure that failures in one window do not halt the entire pipeline.
#' @param ... Arguments passed directly to `process_window`.
#' @return A list containing `result` (output of `process_window` on success, `NULL` on failure)
#'   and `error` (error message on failure, `NULL` on success).
safely_process_window <- function(...) {
  window_args <- list(...)
  window_id <- window_args$window_id
  
  message("Processing window_id: ", window_id)
  
  result <- tryCatch({
    do.call(process_window, window_args)
  }, error = function(e) {
    message("Error processing window_id ", window_id, ": ", e$message)
    list(result = NULL, error = e$message)
  })
  
  if (!is.null(result$error)) {
    return(result) # Already contains error information
  } else {
    return(list(result = result, error = NULL))
  }
}