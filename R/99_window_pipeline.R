# Per-window processing helper used by targets
# This function executes the sequence of steps that were previously
# expressed as many individual targets inside `define_window_targets()`.

process_window <- function(window_id, from_date, to_date, val_date, oos_from_date, oos_to_date, app_config, asset_metadata) {
  # 2. Data Ingestion
  raw_data <- get_all_raw_data(
    asset_metadata = asset_metadata,
    cache_dir = app_config$default$cache_dir,
    from = from_date,
    to = to_date
  )

  # 3. Data Processing
  monthly_returns <- apply_transformations(raw_data = raw_data, asset_metadata = asset_metadata, window_to_date = to_date)
  split_data <- split_data_by_type(monthly_returns, asset_metadata)
  asset_returns <- split_data$asset_returns
  macro_data <- split_data$macro_data

  # 1-month lag for macro data
  macro_original_colnames <- colnames(macro_data)
  lagged_macro_data <- na.omit(xts::lag.xts(macro_data, k = 1))

  # 4. Model fitting
  fitted_rsbvar_output <- tryCatch(
    fit_rsbvar_model(macro_data = lagged_macro_data, bvar_lags = app_config$default$models$bvar_lags, app_config = app_config),
    error = function(e) { message("fit_rsbvar_model failed: ", e$message); return(NULL) }
  )

  fitted_rsbvar_model <- if (!is.null(fitted_rsbvar_output)) fitted_rsbvar_output$fitted_model else NULL
  fitted_rsbvar_spec <- if (!is.null(fitted_rsbvar_output)) fitted_rsbvar_output$spec else NULL

  fitted_dcc_garch_model <- tryCatch(
    fit_dcc_t_garch_model(asset_returns, app_config = app_config),
    error = function(e) { message("fit_dcc_t_garch_model failed: ", e$message); return(NULL) }
  )

  fitted_generative_model <- list(
    rsbvar_fitted_model = fitted_rsbvar_model,
    rsbvar_spec = fitted_rsbvar_spec,
    macro_original_colnames = macro_original_colnames, # Add original macro column names
    dcc_garch_model = fitted_dcc_garch_model
  )

  # 5. Scenario generation
  last_historical_data_levels <- tail(raw_data, 1)
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
      # Return a minimal placeholder so downstream steps can continue in a degraded mode
      return(list(macro_scenarios = NULL, asset_scenarios = array(0, dim = c(app_config$default$n_simulations, app_config$default$horizon_months, ncol(asset_returns)), dimnames = list(NULL, NULL, colnames(asset_returns)))))
    }
  )

  # If simulated_scenarios is NULL or malformed, ensure a valid placeholder
  if (is.null(simulated_scenarios) || is.null(simulated_scenarios$asset_scenarios)) {
    message("Simulated scenarios missing or malformed; creating placeholder zero-return scenarios.")
    simulated_scenarios <- list(macro_scenarios = NULL, asset_scenarios = array(0, dim = c(app_config$default$n_simulations, app_config$default$horizon_months, ncol(asset_returns)), dimnames = list(NULL, NULL, colnames(asset_returns))))
  }

  # Perform sanity check on the first window's scenarios
  if (window_id == 1) {
    sanity_check_dir <- "output/plots/scenario_sanity_check"
    perform_scenario_sanity_check(
      simulated_scenarios = simulated_scenarios,
      historical_returns = asset_returns,
      historical_macro_data = macro_data, # Add historical macro data
      output_dir = sanity_check_dir
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
    risk_aversion_param = app_config$default$models$risk_aversion_param
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
