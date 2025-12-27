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
generate_wf_windows <- function(app_config) {
  
  overall_from_date <- as.Date(app_config$default$data$from)
  overall_to_date <- as.Date(app_config$default$data$to)
  initial_window_months <- app_config$default$initial_window_months
  roll_forward_months <- app_config$default$roll_forward_months
  validation_start_date <- as.Date(app_config$default$validation_start_date)
  
  # Ensure validation_start_date is after the initial window can be formed
  if (validation_start_date < (overall_from_date %m+% months(initial_window_months))) {
    stop("validation_start_date is too early. It must allow for the initial window.")
  }
  
  windows <- tibble()
  current_to_date <- validation_start_date
  window_id_counter <- 1
  
  while (current_to_date <= overall_to_date) {
    # The training window ends at current_to_date.
    # The training window starts `initial_window_months` before current_to_date.
    # The strategy decision is made at current_to_date for the period *after* it.
    
    # Calculate the training window end date. This is the rebalancing date.
    # This date will also be used to filter the *historical data* for model fitting.
    training_window_end <- current_to_date
    
    # Calculate the training window start date
    training_window_start <- training_window_end %m-% months(initial_window_months)
    
    # Check if the training window starts before the overall data start.
    # If so, adjust training_window_start to overall_from_date.
    # This means the initial windows might be shorter than initial_window_months,
    # until it reaches full length. Or, we can enforce minimum length.
    # For now, let's enforce minimum length - window must be at least initial_window_months long.
    if (training_window_start < overall_from_date) {
        # This implies we can't form a full initial window, so we stop if we enforce strict length
        # Or, we can adjust the start date. Let's adjust to be robust.
        training_window_start <- overall_from_date
    }

    # Define the out-of-sample period for this window's weights
    oos_from_date <- training_window_end + days(1)
    oos_to_date_candidate <- training_window_end %m+% months(roll_forward_months)
    oos_to_date <- min(oos_to_date_candidate, overall_to_date) # Truncate if it goes beyond overall_to_date

    # Only include windows that have a valid out-of-sample period
    if (oos_from_date <= oos_to_date) {
      # Store the window
      windows <- bind_rows(windows, tibble(
        window_id = window_id_counter,
        from_date = training_window_start,
        to_date = training_window_end,
        val_date = training_window_end, # Decision date, also end of training period
        oos_from_date = oos_from_date,
        oos_to_date = oos_to_date
      ))
    }
    
    # Roll forward the window
    current_to_date <- current_to_date %m+% months(roll_forward_months)
    window_id_counter <- window_id_counter + 1
  }
  
  # Filter out windows where the training data range is too short if needed.
  # For now, assume a robust generation where first window might be shorter.
  
  message("Generated ", nrow(windows), " walk-forward validation windows.")
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

