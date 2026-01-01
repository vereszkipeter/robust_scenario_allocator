# Per-window processing helper used by targets
# This function executes the sequence of steps that were previously
# expressed as many individual targets inside `define_window_targets()`.

process_window <- function(window_id, from_date, to_date, val_date, oos_from_date, oos_to_date, app_config, asset_metadata, all_raw_data) {
  # Subset all_raw_data for the current window (each xts object in the list)
  current_window_raw_data <- lapply(all_raw_data, function(x) x[paste0(from_date, "::", to_date)])
  # browser()
  # 3. Data Processing
  monthly_returns <- apply_transformations(raw_data_list = current_window_raw_data, asset_metadata = asset_metadata, window_to_date = to_date)
  log_message(paste0("NROW(monthly_returns) after transformations: ", NROW(monthly_returns)), level = "DEBUG", app_config = app_config)
  split_data <- split_data_by_type(monthly_returns, asset_metadata)
  asset_returns <- split_data$asset_returns |> na.omit() # Message to LLM: first row is NA what breaks subsequent calculations
  log_message(paste0("NROW(asset_returns) before DCC-GARCH: ", NROW(asset_returns)), level = "DEBUG", app_config = app_config)
  macro_data <- split_data$macro_data |> na.omit() # Message to LLM: some rows can be NAs what breaks subsequent calculations
  log_message(paste0("Column names in monthly_returns: ", paste(colnames(monthly_returns), collapse = ", ")), level = "DEBUG", app_config = app_config)
  log_message(paste0("Column names in macro_data before lagging: ", paste(colnames(macro_data), collapse = ", ")), level = "DEBUG", app_config = app_config)

  macro_original_colnames <- colnames(macro_data)
  if (is.null(macro_data) || NCOL(macro_data) == 0) {
    lagged_macro_data <- NULL
  } else {
    macro_data_filled <- na.locf(na.locf(macro_data, na.rm = FALSE), fromLast = TRUE)
    lagged_macro_data <- na.omit(lag.xts(macro_data_filled, k = 1))
    log_message(paste0("Column names in lagged_macro_data: ", paste(colnames(lagged_macro_data), collapse = ", ")), level = "DEBUG", app_config = app_config)
    log_message(paste0("Dimensions of lagged_macro_data: ", NROW(lagged_macro_data), " rows, ", NCOL(lagged_macro_data), " cols"), level = "DEBUG", app_config = app_config)
  }

  # Enforce NA policy for modeling: prepare asset_returns for DCC
  min_dcc_obs <- app_config$default$models$min_dcc_obs
  if (!is.null(asset_returns) && NROW(asset_returns) > 0) {
    # Drop any rows with NA in asset returns for DCC estimation
    asset_returns_model <- na.omit(asset_returns)
    log_message(paste0("asset_returns_model rows after na.omit: ", NROW(asset_returns_model)), level = "DEBUG", app_config = app_config)
    if (NROW(asset_returns_model) < min_dcc_obs) {
      message("WARNING: Insufficient asset_returns observations (", NROW(asset_returns_model), ") for reliable DCC estimation. DCC fit may fallback.")
    }
  } else {
    asset_returns_model <- asset_returns
  }
  # browser()
  # 4. Model fitting
  fitted_rsbvar_output <- tryCatch(
    fit_rsbvar_model(macro_data = lagged_macro_data, bvar_lags = app_config$default$models$bvar_lags, app_config = app_config),
    error = function(e) { message("fit_rsbvar_model failed: ", e$message); return(NULL) }
  )

  fitted_rsbvar_model <- if (!is.null(fitted_rsbvar_output)) fitted_rsbvar_output$fitted_model else NULL
  fitted_rsbvar_spec <- if (!is.null(fitted_rsbvar_output)) fitted_rsbvar_output$spec else NULL

  fitted_dcc_garch_model <- tryCatch(
    fit_dcc_t_garch_model(asset_returns_model, app_config = app_config),
    error = function(e) { 
      message("fit_dcc_t_garch_model failed with a detailed error:")
      print(e)
      return(NULL) 
    }
  )

  fitted_generative_model <- list(
    rsbvar_fitted_model = fitted_rsbvar_model,
    rsbvar_spec = fitted_rsbvar_spec,
    macro_original_colnames = macro_original_colnames, # Add original macro column names
    dcc_garch_model = fitted_dcc_garch_model,
    macro_data = lagged_macro_data
  )

  # 5. Scenario generation
  # Robust extraction of the last historical level for each series from the current_window_raw_data
  last_historical_data_levels_list <- list()
  for (ticker_name in names(current_window_raw_data)) { # Iterate over names
    series <- current_window_raw_data[[ticker_name]] # Access list element
    last_finite_val_raw <- tail(series[is.finite(series)], 1) # Get the last finite value
        
        if (NROW(last_finite_val_raw) > 0) {
      # Standardize the index of last_finite_val to the window's to_date
      standardized_xts <- xts(as.numeric(last_finite_val_raw), order.by = to_date)
      colnames(standardized_xts) <- ticker_name
      last_historical_data_levels_list[[ticker_name]] <- standardized_xts
    } else {
      # Fallback: if no finite value at the end, use the mean of finite values in the series for the window
      finite_series_values <- as.numeric(series[is.finite(series)])
      if (length(finite_series_values) > 0) {
        fallback_value <- mean(finite_series_values, na.rm = TRUE)
        # Create an xts object for the fallback value, indexed at the window's to_date
        fallback_xts <- xts(fallback_value, order.by = to_date) # Use to_date as the index
        colnames(fallback_xts) <- ticker_name
        last_historical_data_levels_list[[ticker_name]] <- fallback_xts
        warning(paste0("Macro ticker '", ticker_name, "' had no finite last historical level. Using mean (", round(fallback_value, 2), ") as fallback for last_level."))
      } else {
        # If the entire series is NA/Inf, still don't add it to last_historical_data_levels_list
        # This will correctly lead to NA_real_ for last_level in generate_full_scenarios,
        # which `reconstruct_macro_series` already handles by returning NA series.
        warning(paste0("Macro ticker '", ticker_name, "' is entirely non-finite in the current window. Skipping last_level creation."))
      }
    }
  }
  
  if (length(last_historical_data_levels_list) > 0) {
    # Merge these single-observation xts objects into one xts row for consistency
    last_historical_data_levels <- do.call(merge, last_historical_data_levels_list)
  } else {
    # If no data is available to determine last levels, we can either stop or use a default/warning.
    # For now, let's stop as it's critical for inverse transformations.
    stop("Could not determine last historical data levels for any series from raw_data_list.")
  }

  simulated_scenarios <- tryCatch(
    generate_full_scenarios(
      fitted_generative_model,
      n_sim = app_config$default$n_simulations,
      horizon = app_config$default$horizon_months,
      asset_metadata = asset_metadata,
      last_historical_levels = last_historical_data_levels,
      app_config = app_config
    ),
    error = function(e) {
      message("generate_full_scenarios failed: ", e$message)
      # Attempt a graceful fallback: build an empirical DCC fallback using sample moments
      emp_mean <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) colMeans(as.matrix(asset_returns), na.rm = TRUE) else numeric(0)
      emp_cov <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) cov(as.matrix(asset_returns), use = "pairwise.complete.obs") else matrix(NA_real_)
      if (length(emp_mean) > 0 && !any(is.na(emp_cov))) {
        fallback_dcc <- list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns))
        # Retry scenario generation with fallback DCC object (which generate_full_scenarios handles)
        try({
          fitted_generative_model$dcc_garch_model <- fallback_dcc
          return(generate_full_scenarios(
            fitted_generative_model,
            n_sim = app_config$default$n_simulations,
            horizon = app_config$default$horizon_months,
            asset_metadata = asset_metadata,
            last_historical_levels = last_historical_data_levels,
            app_config = app_config
          ))
        }, silent = TRUE)
      }
      # As a last resort return NA-filled scenarios so diagnostics are obvious (avoid flat zero paths)
      na_asset_array <- array(NA_real_, dim = c(app_config$default$n_simulations, app_config$default$horizon_months, ncol(asset_returns)), dimnames = list(NULL, NULL, colnames(asset_returns)))
      return(list(macro_scenarios = NULL, asset_scenarios = na_asset_array))
    }
  )
# browser()
  # If simulated_scenarios is NULL or malformed, ensure a valid fallback
  if (is.null(simulated_scenarios) || is.null(simulated_scenarios$asset_scenarios)) {
    message("Simulated scenarios missing or malformed; creating fallback zero-return scenarios.")
    simulated_scenarios <- list(macro_scenarios = NULL, asset_scenarios = array(0, dim = c(app_config$default$n_simulations, app_config$default$horizon_months, ncol(asset_returns)), dimnames = list(NULL, NULL, colnames(asset_returns))))
  }

  # Perform sanity check on the first window's scenarios
  if (window_id == 1) {
    sanity_check_dir <- "output/plots/scenario_sanity_check"
    perform_scenario_sanity_check(
      simulated_scenarios = simulated_scenarios,
      historical_returns = asset_returns,
      historical_macro_data = macro_data, # Add historical macro data
      output_dir = sanity_check_dir,
      app_config = app_config
    )
  }

  # 6. Base strategy weights and evaluation on scenarios
  base_strategy_weights <- calculate_base_strategy_weights(asset_returns = asset_returns, app_config = app_config)
  base_strategy_pnl_on_scenarios <- evaluate_strategies_on_scenarios(
    base_strategy_weights = base_strategy_weights,
    simulated_scenarios = simulated_scenarios,
    asset_metadata = asset_metadata,
    app_config = app_config
  )

  # Anchoring / equilibrium
  erc_weights <- calculate_erc_weights(asset_returns)
  implied_equilibrium_returns <- calculate_implied_equilibrium_returns(
    erc_weights = erc_weights,
    asset_returns = asset_returns,
    risk_aversion_param = app_config$default$models$risk_aversion_param,
    app_config = app_config
  )

  # Views & entropy pooling
  short_rate_scenarios <- extract_short_rate_from_rsbvar_scenarios(simulated_scenarios$macro_scenarios, app_config = app_config)
  term_premium_model_output <- model_yield_curve_and_term_premium_view(short_rate_scenarios, app_config = app_config)
  adjusted_scenario_probabilities <- apply_sequential_entropy_pooling(
    simulated_scenarios = simulated_scenarios,
    implied_equilibrium_returns = implied_equilibrium_returns,
    term_premium_model_output = term_premium_model_output,
    app_config = app_config
  )

  entropy_pooled_scenarios <- list(scenarios = simulated_scenarios, probabilities = adjusted_scenario_probabilities)

  # Optimization
  optimal_weights <- perform_distributionally_robust_optimization(
    base_strategy_pnl_on_scenarios = base_strategy_pnl_on_scenarios,
    scenario_probabilities = entropy_pooled_scenarios$probabilities,
    app_config = app_config
  )

  # Out-of-sample performance for RSA
  oos_performance <- calculate_oos_performance(
    optimal_weights = optimal_weights,
    oos_from_date = oos_from_date,
    oos_to_date = oos_to_date,
    asset_metadata = asset_metadata,
    app_config = app_config
  )

  # Benchmarks
  ew_benchmark_weights <- calculate_equal_weight_benchmark(asset_names = (asset_metadata %>% filter(asset_class == "Asset") %>% pull(ticker)))
  ew_benchmark_oos_performance <- calculate_oos_performance(
    optimal_weights = ew_benchmark_weights,
    oos_from_date = oos_from_date,
    oos_to_date = oos_to_date,
    asset_metadata = asset_metadata,
    app_config = app_config
  )

  list(
    window_id = window_id,
    val_date = val_date,
    oos_from_date = oos_from_date,
    oos_to_date = oos_to_date,
    optimal_weights = optimal_weights,
    base_strategy_pnl_on_scenarios = base_strategy_pnl_on_scenarios,
    entropy_pooled_probabilities = entropy_pooled_scenarios$probabilities,
    oos_performance = oos_performance,
    ew_benchmark_oos_performance = ew_benchmark_oos_performance
  )
}
