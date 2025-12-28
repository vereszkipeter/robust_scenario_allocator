# Custom simple return function
simple_return <- function(x) {
  x / stats::lag(x) - 1
}

# Custom inverse simple return function
undo_simple_return <- function(x, last_level) {
  # The first value of x is the return from last_level to the first period.
  # So the first level is last_level * (1 + x[1])
  # The second level is first_level * (1 + x[2]), and so on.
  # This can be calculated with a cumulative product.
  
  # Ensure x is a numeric vector
  x_numeric <- as.numeric(x)
  
  # Replace NA with 0 for the cumulative product, assuming no return is 0 return
  x_numeric[is.na(x_numeric)] <- 0
  
  # Calculate cumulative returns
  cumulative_returns <- cumprod(1 + x_numeric)
  
  # Apply to the last level
  reconstructed_levels <- last_level * cumulative_returns
  
  return(reconstructed_levels)
}

# Identity function for transformations that do nothing
identity_transform <- function(x) {
  return(x)
}

identity_transform_inverse <- function(x, last_level) {
  return(x)
}


# Month-over-month change for inflation-like series
mom_change <- function(x) {
  x / stats::lag(x, k = 1) - 1
}

# Inverse of month-over-month change
undo_mom_change <- function(sim_series, last_level) {
  reconstructed <- numeric(length(sim_series))
  reconstructed[1] <- last_level * (1 + sim_series[1])
  if (length(sim_series) > 1) {
    for (t in 2:length(sim_series)) {
      reconstructed[t] <- reconstructed[t-1] * (1 + sim_series[t])
    }
  }
  return(reconstructed)
}

# Ensure app_config has all necessary defaults to avoid NULLs
ensure_app_config_defaults <- function(config) {
  
  # Set default for dcc_solver if not present
  if (is.null(config$default$models$dcc_solver)) {
    config$default$models$dcc_solver <- "solnp"
  }
  
  # Set default for dcc_solver_control if not present
  if (is.null(config$default$models$dcc_solver_control)) {
    config$default$models$dcc_solver_control <- list(trace = 0)
  }
  
  # Set a default for min_dcc_obs if not present
  if (is.null(config$default$models$min_dcc_obs)) {
    config$default$models$min_dcc_obs <- 100 # A reasonable default
  }
  
  return(config)
}

# A safer version of process_window that returns a structured error object
# This allows the overall `tar_make` process to continue even if one window fails
safely_process_window <- function(...) {
  tryCatch({
    # The `...` will capture all arguments passed by purrr::pmap
    process_window(...)
  }, error = function(e) {
    # Extract window_id from the arguments passed to the function
    args <- list(...)
    window_id_val <- if (!is.null(args$window_id)) args$window_id else "UNKNOWN"
    
    # Log the detailed error for debugging
    message(sprintf("Error processing window_id %s: %s", window_id_val, e$message))
    
    # Return a list with an error field, and NULL for other expected outputs
    list(
      window_id = window_id_val,
      error = e$message,
      val_date = NULL,
      oos_from_date = NULL,
      oos_to_date = NULL,
      optimal_weights = NULL,
      base_strategy_pnl_on_scenarios = NULL,
      entropy_pooled_probabilities = NULL,
      oos_performance = NULL,
      ew_benchmark_oos_performance = NULL
    )
  })
}