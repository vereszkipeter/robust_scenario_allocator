# Assuming get_all_raw_data, apply_transformations are sourced from R/01_data.R
# Assuming data_cache_dir is passed from app_config

#' @title Evaluate Base Strategies on Simulated Scenarios
#' @description Calculates the cumulative P&L for each base strategy across all
#'   simulated asset scenarios. This effectively "backtests" each strategy
#'   on the generated future paths.
#' @param base_strategy_weights A list where each element is a named vector of
#'   weights for a strategy, with asset tickers as names.
#' @param simulated_scenarios A list containing:
#'   - `macro_scenarios`: 3D array (n_sim, horizon, n_macro_vars)
#'   - `asset_scenarios`: 3D array (n_sim, horizon, n_assets) (simple returns)
#' @param asset_metadata A tibble containing asset information, used to align
#'   asset tickers with the scenario data.
#' @param app_config A list containing application configuration.
#' @return A 3D array (horizon, n_strategies, n_sim) of cumulative P&L (equity curves)
#'   for each strategy and each simulation, starting from 1.
evaluate_strategies_on_scenarios <- function(base_strategy_weights, simulated_scenarios, asset_metadata, app_config) {

  asset_scenarios <- simulated_scenarios$asset_scenarios
  n_sim <- dim(asset_scenarios)[1]
  horizon <- dim(asset_scenarios)[2]
  n_assets_in_scenarios <- dim(asset_scenarios)[3]
  asset_names_in_scenarios <- dimnames(asset_scenarios)[[3]]

  strategy_names <- names(base_strategy_weights)
  n_strategies <- length(strategy_names)

  # Initialize a 3D array to store P&L: (horizon + 1) x n_strategies x n_sim
  # The +1 is for the initial capital of 1.
  strategy_pnl_on_scenarios <- array(1, dim = c(horizon + 1, n_strategies, n_sim),
                                     dimnames = list(
                                       c("Initial", paste0("Month_", 1:horizon)),
                                       strategy_names,
                                       paste0("Sim_", 1:n_sim)
                                     ))
  
  # Loop through each simulation
  for (s in 1:n_sim) {
    # Extract asset returns for the current simulation: horizon x n_assets_in_scenarios
    # ensure it's a matrix for matrix multiplication
    current_sim_asset_returns <- matrix(asset_scenarios[s, , ], ncol = n_assets_in_scenarios)
    colnames(current_sim_asset_returns) <- asset_names_in_scenarios

    # Loop through each strategy
    for (strat_idx in 1:n_strategies) {
      strategy_name <- strategy_names[strat_idx]
      weights <- base_strategy_weights[[strategy_name]] # Named vector of weights
      
      # DEBUGGING 1 START
      # message("DEBUG: Strategy Name: ", strategy_name)
      # message("DEBUG: Weights class: ", class(weights))
      # message("DEBUG: Weights length: ", length(weights))
      # DEBUGGING 1 END

      # Ensure weights are aligned with asset_names_in_scenarios
      # Handle cases where strategy might have weights for assets not in scenarios
      # or vice-versa.
      ordered_weights <- weights[asset_names_in_scenarios]
      ordered_weights[is.na(ordered_weights)] <- 0 # Assign 0 weight to missing assets
      
      # DEBUGGING 2 START
      # message("DEBUG: Ordered weights class (before as.numeric): ", class(ordered_weights))
      # message("DEBUG: Ordered weights length (before as.numeric): ", length(ordered_weights))
      # DEBUGGING 2 END

      # Ensure ordered_weights is a numeric vector
      ordered_weights <- as.numeric(ordered_weights)
      
      # DEBUGGING 3 START
      # message("DEBUG: Ordered weights class (after as.numeric): ", class(ordered_weights))
      # message("DEBUG: Ordered weights length (after as.numeric): ", length(ordered_weights))
      # message("DEBUG: Current sim asset returns dim: ", paste(dim(current_sim_asset_returns), collapse = "x"))
      # message("DEBUG: Current sim asset returns class: ", class(current_sim_asset_returns))
      # DEBUGGING 3 END

      portfolio_returns_current_strat <- current_sim_asset_returns %*% ordered_weights

      # Calculate cumulative P&L (equity curve) starting from 1
      # The first element of strategy_pnl_on_scenarios is 1 (initial capital)
      # cumprod(1 + returns) will give (1+R1), (1+R1)(1+R2), ...
      strategy_pnl_on_scenarios[2:(horizon + 1), strat_idx, s] <- cumprod(1 + portfolio_returns_current_strat)
    }
  }

  return(strategy_pnl_on_scenarios)
}


#' @title Calculate Out-of-Sample Performance
#' @description Fetches actual historical asset returns for a specified out-of-sample period,
#'   applies optimal portfolio weights, and computes key performance metrics.
#' @param optimal_weights A numeric vector of optimal portfolio weights (from the meta-learner).
#' @param oos_from_date A character string or Date object for the start of the out-of-sample period.
#' @param oos_to_date A character string or Date object for the end of the out-of-sample period.
#' @param asset_metadata A tibble containing asset information, used to align tickers.
#' @param app_config A list containing application configuration.
#' @return A list containing out-of-sample performance metrics and the OOS portfolio returns.
calculate_oos_performance <- function(optimal_weights, oos_from_date, oos_to_date, asset_metadata, app_config) {
  
  # Ensure necessary functions are available (get_all_raw_data, apply_transformations)
  # These should be sourced globally by _targets.R.

  # 1. Fetch historical data for the OOS period
  oos_raw_data <- get_all_raw_data(
    asset_metadata = asset_metadata,
    cache_dir = app_config$default$cache_dir,
    from = oos_from_date,
    to = oos_to_date
  )

  # 2. Transform raw data to monthly returns for the OOS period
  oos_monthly_returns <- apply_transformations(
    raw_data = oos_raw_data,
    asset_metadata = asset_metadata
  )
  
  # Filter only asset returns (excluding macro variables)
  oos_asset_returns <- oos_monthly_returns[, colnames(oos_monthly_returns) %in% 
                                             (asset_metadata %>% filter(asset_class == "Asset") %>% pull(ticker))]
  
  if (NROW(oos_asset_returns) == 0) {
    warning("No out-of-sample asset returns found for period: ", oos_from_date, " to ", oos_to_date)
    return(list(
      oos_returns = NULL,
      metrics = list(
        cumulative_return = NA,
        sharpe_ratio = NA,
        max_drawdown = NA,
        cdaR = NA
      )
    ))
  }

  # 3. Align weights with OOS asset returns
  asset_names_in_oos <- colnames(oos_asset_returns)
  ordered_weights <- optimal_weights[asset_names_in_oos]
  ordered_weights[is.na(ordered_weights)] <- 0 # Assign 0 weight to missing assets
  
  # Ensure weights sum to 1, or handle if there are missing assets
  if (sum(ordered_weights) == 0 && length(ordered_weights) > 0) {
      warning("All optimal weights became zero for OOS period due to missing assets. Cannot compute OOS performance.")
      return(list(
        oos_returns = NULL,
        metrics = list(
          cumulative_return = NA,
          sharpe_ratio = NA,
          max_drawdown = NA,
          cdaR = NA
        )
      ))
  } else if (sum(ordered_weights) != 1 && sum(ordered_weights) != 0) {
    # If not zero, re-normalize to 1 to ensure full investment if some assets were missing
    ordered_weights <- ordered_weights / sum(ordered_weights)
  }

  # 4. Calculate OOS portfolio returns
  oos_portfolio_returns <- oos_asset_returns %*% ordered_weights
  colnames(oos_portfolio_returns) <- "Portfolio"

  # 4.1 Apply transaction costs
  transaction_cost <- app_config$default$models$transaction_cost_flat_rate
  if (NROW(oos_portfolio_returns) > 0 && !is.null(transaction_cost) && transaction_cost > 0) {
    # Deduct transaction cost from the first return. This is a one-time cost for rebalancing.
    oos_portfolio_returns[1, 1] <- oos_portfolio_returns[1, 1] - transaction_cost
    message("Applied transaction cost of ", transaction_cost, " to OOS returns starting ", index(oos_portfolio_returns[1,]))
  }
  
  # 5. Compute performance metrics
  metrics <- list()
  
  # Cumulative Return
  metrics$cumulative_return <- PerformanceAnalytics::Return.cumulative(oos_portfolio_returns)
  
  # Sharpe Ratio (assuming 0 risk-free rate for simplicity for now)
  metrics$sharpe_ratio <- PerformanceAnalytics::SharpeRatio.annualized(oos_portfolio_returns, Rf = 0, scale = 12)
  
  # Maximum Drawdown
  metrics$max_drawdown <- PerformanceAnalytics::maxDrawdown(oos_portfolio_returns)
  
  # CDaR (Conditional Drawdown at Risk) - using a simple percentile-based approach on historical data
  # This is for realized historical returns, not simulated.
  # Use PerformanceAnalytics::ES (Expected Shortfall) on drawdowns.
  # First, calculate drawdowns from the OOS portfolio returns.
  oos_portfolio_pnl <- PerformanceAnalytics::Return.cumulative(oos_portfolio_returns, geometric = TRUE) + 1 # Equity curve
  oos_portfolio_drawdowns <- PerformanceAnalytics::Drawdowns(oos_portfolio_pnl) # Returns negative values
  
  # Calculate ES of these drawdowns (at a given alpha)
  cvar_alpha <- app_config$default$models$cvar_alpha # Use alpha from config
  metrics$cdaR <- PerformanceAnalytics::ES(oos_portfolio_drawdowns, p = cvar_alpha, method = "historical") # ES of drawdowns
  
  # Calculate PSR and DSR
  confidence_level_psr_dsr <- app_config$default$models$confidence_level
  # For DSR, 'n_strategies' (M) should be the number of strategies considered in the entire selection process.
  # If we are evaluating the single optimal portfolio selected by the meta-learner, M=1 for that portfolio's DSR.
  # If we were evaluating all base strategies and picking the best, M would be `n_strategies_base`.
  # For now, let's pass a constant M=1, assuming we're evaluating the single "final" portfolio.
  # Future enhancement: dynamically pass M (number of base strategies considered) to reflect multiple comparisons.
  psr_dsr_results <- calculate_psr_dsr(
    R = oos_portfolio_returns,
    confidence_level = confidence_level_psr_dsr,
    n_strategies = 1 # Assuming evaluation of a single selected portfolio
  )
  metrics$psr <- psr_dsr_results$psr
  metrics$dsr <- psr_dsr_results$dsr

      return(list(
        oos_returns = oos_portfolio_returns,
        metrics = metrics
      ))
  }
  
  #' @title Compare Simulated and Historical Distributions
  #' @description Calculates mean, standard deviation, and performs a Kolmogorov-Smirnov test
  #'   between simulated terminal values and historical data for a given variable.
  #' @param simulated_data A numeric vector of simulated terminal values.
  #' @param historical_data A numeric vector of historical values.
  #' @param var_name The name of the variable being compared.
  #' @return A character string summarizing the comparison.
  compare_distributions <- function(simulated_data, historical_data, var_name) {
    sim_mean <- mean(simulated_data, na.rm = TRUE)
    sim_sd <- sd(simulated_data, na.rm = TRUE)
    hist_mean <- mean(historical_data, na.rm = TRUE)
    hist_sd <- sd(historical_data, na.rm = TRUE)
    
    # Ensure historical_data has sufficient non-NA values for ks.test
    # ks.test requires at least 2 non-NA values
    if (length(na.omit(historical_data)) < 2 || length(na.omit(simulated_data)) < 2) {
        ks_test_result <- "KS test skipped (insufficient data)"
    } else {
        ks_test <- ks.test(simulated_data, historical_data)
        ks_test_result <- paste0("KS Test p-value: ", round(ks_test$p.value, 4))
    }
  
    return(
      paste0(
        "\n--- ", var_name, " Comparison ---\n",
        "Simulated Mean: ", round(sim_mean, 5), ", SD: ", round(sim_sd, 5), "\n",
        "Historical Mean: ", round(hist_mean, 5), ", SD: ", round(hist_sd, 5), "\n",
        ks_test_result, "\n"
      )
    )
  }
  
#' @title Perform Sanity Check on Simulated Scenarios
#' @description Generates diagnostic plots to visually inspect the quality of the
#'   simulated scenarios against historical data.
#' @param simulated_scenarios A list containing `asset_scenarios` and `macro_scenarios`.
#' @param historical_returns An xts object of historical monthly returns.
#' @param historical_macro_data An xts object of historical monthly macro data.
#' @param output_dir The directory to save the diagnostic plots.
#' @return The file path to a dummy indicator file.
perform_scenario_sanity_check <- function(simulated_scenarios, historical_returns, historical_macro_data, output_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  sanity_check_messages <- c("Scenario Sanity Check Report:\n")

  asset_scenarios <- simulated_scenarios$asset_scenarios
  macro_scenarios <- simulated_scenarios$macro_scenarios # Extract macro scenarios
  
  n_sim <- dim(asset_scenarios)[1]
  horizon <- dim(asset_scenarios)[2]
  
  # --- ASSET SCENARIOS SANITY CHECK ---
  if (is.null(asset_scenarios)) {
    warning("Asset scenarios are NULL. Skipping asset scenarios sanity check.")
    sanity_check_messages <- c(sanity_check_messages, "Asset scenarios are NULL. Skipping asset scenarios sanity check.\n")
  } else {
    asset_names <- dimnames(asset_scenarios)[[3]]
    
    # --- 1. Fan Chart of Simulated Asset Paths ---
    # Reshape data for ggplot
    sims_long <- as.data.frame.table(asset_scenarios, responseName = "Return") %>%
      rename(Sim = Var1, Horizon = Var2, Asset = Var3) %>%
      mutate(Horizon = as.numeric(Horizon))
    
    # Calculate quantiles for fan chart
    fan_chart_data <- sims_long %>%
      group_by(Asset, Horizon) %>%
      summarise(
        p10 = quantile(Return, 0.10),
        p25 = quantile(Return, 0.25),
        p50 = quantile(Return, 0.50),
        p75 = quantile(Return, 0.75),
        p90 = quantile(Return, 0.90)
      )
    
    # Plot fan chart for each asset
    for (asset_name in asset_names) {
      p <- ggplot(filter(fan_chart_data, Asset == asset_name), aes(x = Horizon)) +
        geom_ribbon(aes(ymin = p10, ymax = p90), fill = "skyblue", alpha = 0.4) +
        geom_ribbon(aes(ymin = p25, ymax = p75), fill = "dodgerblue", alpha = 0.6) +
        geom_line(aes(y = p50), color = "darkblue", linewidth = 1) +
        labs(
          title = paste("Fan Chart of Simulated Asset Returns for", asset_name),
          subtitle = paste(n_sim, "simulations over", horizon, "months"),
          x = "Horizon (Months)",
          y = "Simulated Monthly Return"
        ) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0("fanchart_asset_", asset_name, ".png")), plot = p, width = 10, height = 6)
    }
    
    # --- 2. Histograms of Terminal Asset Returns vs. Historical ---
    # Get terminal values (end of horizon)
    terminal_returns <- as.data.frame(asset_scenarios[, horizon, ])
    
    # Plot histograms and perform statistical comparison
    for (asset_name in asset_names) {
      # Ensure historical_returns has the correct column for the current asset
      if (!(asset_name %in% colnames(historical_returns))) {
        warning(paste("Historical returns for asset", asset_name, "not found. Skipping histogram and statistical comparison for this asset."))
        sanity_check_messages <- c(sanity_check_messages, paste("Historical returns for asset", asset_name, "not found. Skipping statistical comparison.\n"))
        next
      }
      p <- ggplot() +
        geom_histogram(data = terminal_returns, aes(x = .data[[asset_name]], y = after_stat(density)), bins = 50, fill = "lightblue", alpha = 0.7) +
        geom_density(data = as.data.frame(historical_returns), aes(x = .data[[asset_name]]), color = "red", size = 1) +
        labs(
          title = paste("Distribution of Terminal vs. Historical Asset Returns for", asset_name),
          subtitle = "Blue: Simulated (at horizon), Red: Historical",
          x = "Monthly Return",
          y = "Density"
        ) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0("histogram_asset_", asset_name, ".png")), plot = p, width = 10, height = 6)

      # Add statistical comparison for asset
      sim_data_asset <- terminal_returns[[asset_name]]
      hist_data_asset <- as.numeric(historical_returns[, asset_name]) # Convert xts column to numeric vector
      sanity_check_messages <- c(sanity_check_messages, compare_distributions(simulated_data = sim_data_asset, historical_data = hist_data_asset, var_name = paste0("Asset: ", asset_name)))
    }
  }
  
  # --- MACRO SCENARIOS SANITY CHECK ---
  if (is.null(macro_scenarios)) {
    warning("Macro scenarios are NULL. Skipping macro scenarios sanity check.")
    sanity_check_messages <- c(sanity_check_messages, "Macro scenarios are NULL. Skipping macro scenarios sanity check.\n")
  } else {
    macro_names <- dimnames(macro_scenarios)[[3]]
    
    # --- 3. Fan Chart of Simulated Macro Paths ---
    macro_sims_long <- as.data.frame.table(macro_scenarios, responseName = "Value") %>%
      rename(Sim = Var1, Horizon = Var2, MacroVar = Var3) %>%
      mutate(Horizon = as.numeric(Horizon))
    
    macro_fan_chart_data <- macro_sims_long %>%
      group_by(MacroVar, Horizon) %>%
      summarise(
        p10 = quantile(Value, 0.10),
        p25 = quantile(Value, 0.25),
        p50 = quantile(Value, 0.50),
        p75 = quantile(Value, 0.75),
        p90 = quantile(Value, 0.90)
      )
    
    for (macro_name in macro_names) {
      p <- ggplot(filter(macro_fan_chart_data, MacroVar == macro_name), aes(x = Horizon)) +
        geom_ribbon(aes(ymin = p10, ymax = p90), fill = "lightgreen", alpha = 0.4) +
        geom_ribbon(aes(ymin = p25, ymax = p75), fill = "forestgreen", alpha = 0.6) +
        geom_line(aes(y = p50), color = "darkgreen", linewidth = 1) +
        labs(
          title = paste("Fan Chart of Simulated Macro Variable for", macro_name),
          subtitle = paste(n_sim, "simulations over", horizon, "months"),
          x = "Horizon (Months)",
          y = "Simulated Macro Value"
        ) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0("fanchart_macro_", macro_name, ".png")), plot = p, width = 10, height = 6)
    }
    
    # --- 4. Histograms of Terminal Macro Values vs. Historical ---
    terminal_macro_values <- as.data.frame(macro_scenarios[, horizon, ])
    
    for (macro_name in macro_names) {
      if (!(macro_name %in% colnames(historical_macro_data))) {
        warning(paste("Historical data for macro variable", macro_name, "not found. Skipping histogram and statistical comparison for this variable."))
        sanity_check_messages <- c(sanity_check_messages, paste("Historical data for macro variable", macro_name, "not found. Skipping statistical comparison.\n"))
        next
      }
      p <- ggplot() +
        geom_histogram(data = terminal_macro_values, aes(x = .data[[macro_name]], y = after_stat(density)), bins = 50, fill = "lightcoral", alpha = 0.7) +
        geom_density(data = as.data.frame(historical_macro_data), aes(x = .data[[macro_name]]), color = "darkred", size = 1) +
        labs(
          title = paste("Distribution of Terminal vs. Historical Macro Variable for", macro_name),
          subtitle = "Red: Simulated (at horizon), Dark Red: Historical",
          x = "Macro Variable Value",
          y = "Density"
        ) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0("histogram_macro_", macro_name, ".png")), plot = p, width = 10, height = 6)

      # Add statistical comparison for macro variable
      sim_data_macro <- terminal_macro_values[[macro_name]]
      hist_data_macro <- as.numeric(historical_macro_data[, macro_name]) # Convert xts column to numeric vector
      sanity_check_messages <- c(sanity_check_messages, compare_distributions(simulated_data = sim_data_macro, historical_data = hist_data_macro, var_name = paste0("Macro: ", macro_name)))
    }
  }
  
  # Create a dummy file to satisfy targets
  output_file <- file.path(output_dir, "scenario_sanity_check_complete.txt")
  writeLines(sanity_check_messages, output_file) # Write all collected messages
  message("Scenario sanity check report generated:", output_file)
  return(output_file)
}