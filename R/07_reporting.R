library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer
library(purrr) # For map_dfr
library(PerformanceAnalytics) # For financial metrics like CDaR
library(Hmisc) # For weighted quantiles
library(yaml) # For writing reports in YAML

#' @title Plot Fan Chart of Final RSA Portfolio
#' @description Generates a fan chart showing the distribution of the final RSA
#' portfolio's cumulative P&L over time, based on simulated scenarios and adjusted probabilities.
#' The chart displays key quantiles (e.g., 5%, 25%, 50%, 75%, 95%) to visualize uncertainty.
#' @param base_strategy_pnl_on_scenarios A 3D array (time, strategy, simulation)
#'   of cumulative P&L for each base strategy across all scenarios.
#' @param optimal_weights A numeric vector of optimal meta-weights for combining base strategies.
#' @param scenario_probabilities A numeric vector of adjusted probabilities for each scenario.
#' @param output_plots_dir The directory path where the fan chart plot should be saved.
#' @return The file path to the saved fan chart plot.
plot_fan_chart <- function(base_strategy_pnl_on_scenarios, optimal_weights, scenario_probabilities, output_plots_dir) {
  
  if (!dir.exists(output_plots_dir)) {
    dir.create(output_plots_dir, recursive = TRUE)
  }

  n_time <- dim(base_strategy_pnl_on_scenarios)[1]
  n_strategies <- dim(base_strategy_pnl_on_scenarios)[2]
  n_sim <- dim(base_strategy_pnl_on_scenarios)[3]
  
  # 1. Calculate combined portfolio P&L for each scenario
  # Result: n_time x n_sim matrix
  combined_pnl_scenarios <- matrix(NA, nrow = n_time, ncol = n_sim)
  
  for (s in 1:n_sim) {
    # base_strategy_pnl_on_scenarios[, , s] is a (n_time x n_strategies) matrix
    # optimal_weights is a (n_strategies x 1) vector
    combined_pnl_scenarios[, s] <- base_strategy_pnl_on_scenarios[, , s] %*% optimal_weights
  }
  
  # 2. For each time step, calculate weighted quantiles
  quantiles_to_plot <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  quantile_names <- c("q05", "q25", "q50", "q75", "q95")
  
  fan_data <- map_dfr(1:n_time, function(t) {
    pnl_at_t <- combined_pnl_scenarios[t, ]
    
    # Calculate weighted quantiles
    weighted_quantiles <- Hmisc::wtd.quantile(
      x = pnl_at_t,
      weights = scenario_probabilities,
      probs = quantiles_to_plot,
      type = "i" # Ensure interpolation type is consistent
    )
    
    data.frame(
      Time = t,
      q05 = weighted_quantiles[1],
      q25 = weighted_quantiles[2],
      q50 = weighted_quantiles[3],
      q75 = weighted_quantiles[4],
      q95 = weighted_quantiles[5]
    )
  })
  
  # 3. Plot using ggplot2
  p <- ggplot(fan_data, aes(x = Time)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = "blue", alpha = 0.2) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "blue", alpha = 0.4) +
    geom_line(aes(y = q50), color = "darkblue", linewidth = 1) +
    labs(
      title = "Fan Chart of Final RSA Portfolio Cumulative P&L",
      y = "Cumulative P&L (Starting at 1)",
      x = "Months into Future"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  file_path <- file.path(output_plots_dir, "rsa_fan_chart.png")
  ggsave(file_path, plot = p, width = 10, height = 6)
  
  message("Fan Chart saved to: ", file_path)
  return(file_path)
}

#' @title Calculate Validation Metrics for a Single Window
#' @description Calculates key validation metrics (CDaR, Sharpe Ratio, Effective Number of Strategies, Leverage)
#'   for the RSA portfolio based on simulated scenarios for a single walk-forward window.
#'   This is an internal helper function.
#' @param base_strategy_pnl_on_scenarios A 3D array (time, strategy, simulation)
#'   of cumulative P&L for each base strategy across all scenarios for this window.
#' @param optimal_weights A numeric vector of optimal meta-weights for this window.
#' @param scenario_probabilities A numeric vector of adjusted probabilities for each scenario for this window.
#' @param app_config A list containing application configuration parameters.
#' @return A list containing calculated metrics for the window.
calculate_window_metrics <- function(base_strategy_pnl_on_scenarios, optimal_weights, scenario_probabilities, app_config) {

  n_time <- dim(base_strategy_pnl_on_scenarios)[1]
  n_sim <- dim(base_strategy_pnl_on_scenarios)[3]

  # 1. Calculate combined portfolio P&L for each scenario
  combined_pnl_scenarios <- matrix(NA, nrow = n_time, ncol = n_sim)
  for (s in 1:n_sim) {
    combined_pnl_scenarios[, s] <- base_strategy_pnl_on_scenarios[, , s] %*% optimal_weights
  }

  # 2. Derive combined portfolio returns for each scenario
  combined_returns_scenarios <- apply(combined_pnl_scenarios, 2, function(x) {
    c(NA, diff(x) / head(x, -1))
  })
  combined_returns_scenarios <- combined_returns_scenarios[-1, ] # Remove first NA row

  metrics <- list()

  # --- Metrics ---

  # CDaR (Conditional Drawdown at Risk)
  max_drawdowns_per_scenario <- apply(combined_pnl_scenarios, 2, function(x) {
    PerformanceAnalytics::maxDrawdown(xts(x, order.by = 1:length(x)))
  })
  
  sorted_drawdowns <- sort(max_drawdowns_per_scenario)
  sorted_probabilities <- scenario_probabilities[order(max_drawdowns_per_scenario)]
  
  alpha_cvar <- app_config$default$models$cvar_alpha # Use alpha from config
  
  cum_probs <- cumsum(sorted_probabilities)
  var_index <- which(cum_probs >= (1 - alpha_cvar))[1]
  
  if (!is.na(var_index) && var_index <= length(sorted_drawdowns)) {
    cvar_drawdowns <- sorted_drawdowns[var_index:length(sorted_drawdowns)]
    cvar_probabilities <- sorted_probabilities[var_index:length(sorted_probabilities)]
    cvar_probabilities <- cvar_probabilities / sum(cvar_probabilities)
    metrics$cdaR <- sum(cvar_drawdowns * cvar_probabilities)
  } else {
    metrics$cdaR <- NA
  }

  # Sharpe Ratio
  expected_returns_per_scenario_path <- apply(combined_returns_scenarios, 2, mean)
  std_dev_returns_per_scenario_path <- apply(combined_returns_scenarios, 2, sd)
  
  expected_portfolio_return_avg <- sum(expected_returns_per_scenario_path * scenario_probabilities)
  expected_portfolio_sd_avg <- sum(std_dev_returns_per_scenario_path * scenario_probabilities)
  
  risk_free_rate <- 0 # Assume 0 for now, could be passed in app_config
  
  if (expected_portfolio_sd_avg > 0) {
    metrics$sharpe_ratio <- (expected_portfolio_return_avg - risk_free_rate) / expected_portfolio_sd_avg
  } else {
    metrics$sharpe_ratio <- NA
  }
  
  # --- Sanity Checks ---

  # Effective Number of Strategies
  metrics$effective_num_strategies <- 1 / sum(optimal_weights^2)
  
  # Leverage (for meta-weights)
  metrics$meta_leverage <- sum(abs(optimal_weights))
  
  return(metrics)
}


#' @title Generate Global Validation Report for Walk-Forward Analysis
#' @description Aggregates validation metrics from all walk-forward windows
#'   and generates a comprehensive report.
#' @param all_window_results A list of lists, where each inner list contains
#'   results for a single walk-forward window (`window_id`, `val_date`,
#'   `oos_from_date`, `oos_to_date`, `optimal_weights`, `base_strategy_pnl_on_scenarios`,
#'   `entropy_pooled_probabilities`, `oos_performance`).
#' @param app_config A list containing application configuration parameters.
#' @param output_report_dir The directory path where the global report should be saved.
#' @param output_plots_dir The directory path where global plots should be saved.
#' @return The file path to the saved global validation report.
generate_global_validation_report <- function(all_window_results, app_config, output_report_dir, output_plots_dir) {
  
  if (!dir.exists(output_report_dir)) {
    dir.create(output_report_dir, recursive = TRUE)
  }
  if (!dir.exists(output_plots_dir)) {
    dir.create(output_plots_dir, recursive = TRUE)
  }

  all_in_sample_metrics_df <- tibble()
  all_oos_metrics_df <- tibble()
  all_ew_benchmark_oos_metrics_df <- tibble() # New: for benchmark metrics
  
  # List to hold OOS returns for combined equity curve
  all_oos_returns_list <- list()
  all_ew_benchmark_oos_returns_list <- list() # New: for benchmark OOS returns

  # Loop through each window's results
  successful_results <- list()
  failed_windows_info <- list()

  for (res_wrapped in all_window_results) {
    if (is.null(res_wrapped$error)) {
      successful_results[[length(successful_results) + 1]] <- res_wrapped$result
    } else {
      
      # Assuming window_id is always present in the result, even if partial
      window_id_failed <- if (!is.null(res_wrapped$result$window_id)) res_wrapped$result$window_id else "UNKNOWN"
      failed_windows_info[[length(failed_windows_info) + 1]] <- list(
        window_id = window_id_failed,
        error = res_wrapped$error
      )
      message("Skipping window due to error: Window ID ", window_id_failed, ", Error: ", res_wrapped$error)
    }
  }

  if (length(failed_windows_info) > 0) {
    warning(length(failed_windows_info), " walk-forward window(s) failed. Check logs for details.")
    # Optionally, save failed_windows_info to a separate log file
    # for example:
    # yaml::write_yaml(failed_windows_info, file = file.path(output_report_dir, "failed_windows_log.yml"))
  }
  
  if (length(successful_results) == 0) {
    stop("No successful walk-forward windows were processed. Cannot generate global report.")
  }

  # Loop through each *successful* window's results
  for (res in successful_results) {
    # Calculate In-Sample (Simulated Future) Metrics
    in_sample_metrics <- calculate_window_metrics(
      base_strategy_pnl_on_scenarios = res$base_strategy_pnl_on_scenarios,
      optimal_weights = res$optimal_weights,
      scenario_probabilities = res$entropy_pooled_probabilities,
      app_config = app_config
    )
    
    in_sample_metric_row <- as_tibble(in_sample_metrics) %>%
      mutate(window_id = res$window_id, val_date = as.Date(res$val_date))
    
    all_in_sample_metrics_df <- bind_rows(all_in_sample_metrics_df, in_sample_metric_row)

    # Extract Out-of-Sample (Realized Historical) Metrics for RSA
    oos_metrics_raw <- res$oos_performance$metrics
    oos_returns_xts <- res$oos_performance$oos_returns

    if (!is.null(oos_metrics_raw)) {
      oos_metric_row <- as_tibble(oos_metrics_raw) %>%
        rename_with(~ paste0("oos_", .x)) %>% # Prefix OOS metrics
        mutate(window_id = res$window_id, val_date = as.Date(res$val_date),
               oos_from_date = as.Date(res$oos_from_date), oos_to_date = as.Date(res$oos_to_date))
      all_oos_metrics_df <- bind_rows(all_oos_metrics_df, oos_metric_row)
    }

    if (!is.null(oos_returns_xts) && NROW(oos_returns_xts) > 0) {
      all_oos_returns_list[[as.character(res$val_date)]] <- oos_returns_xts
    }

    # New: Extract Out-of-Sample (Realized Historical) Metrics for EW Benchmark
    ew_oos_metrics_raw <- res$ew_benchmark_oos_performance$metrics
    ew_oos_returns_xts <- res$ew_benchmark_oos_performance$oos_returns

    if (!is.null(ew_oos_metrics_raw)) {
      ew_oos_metric_row <- as_tibble(ew_oos_metrics_raw) %>%
        rename_with(~ paste0("ew_oos_", .x)) %>% # Prefix EW OOS metrics
        mutate(window_id = res$window_id, val_date = as.Date(res$val_date),
               oos_from_date = as.Date(res$oos_from_date), oos_to_date = as.Date(res$oos_to_date))
      all_ew_benchmark_oos_metrics_df <- bind_rows(all_ew_benchmark_oos_metrics_df, ew_oos_metric_row)
    }

    if (!is.null(ew_oos_returns_xts) && NROW(ew_oos_returns_xts) > 0) {
      all_ew_benchmark_oos_returns_list[[as.character(res$val_date)]] <- ew_oos_returns_xts
    }
  }
  
  # --- Aggregation and Global Report ---
  
  # In-Sample Metrics Summary
  in_sample_summary_metrics <- all_in_sample_metrics_df %>%
    summarise(
      mean_cdaR = mean(cdaR, na.rm = TRUE),
      sd_cdaR = sd(cdaR, na.rm = TRUE),
      mean_sharpe_ratio = mean(sharpe_ratio, na.rm = TRUE),
      sd_sharpe_ratio = sd(sharpe_ratio, na.rm = TRUE),
      mean_effective_num_strategies = mean(effective_num_strategies, na.rm = TRUE),
      mean_meta_leverage = mean(meta_leverage, na.rm = TRUE)
    )

  # Out-of-Sample Metrics Summary for RSA
  oos_summary_metrics <- all_oos_metrics_df %>%
    summarise(
      mean_oos_cumulative_return = mean(oos_cumulative_return, na.rm = TRUE),
      mean_oos_sharpe_ratio = mean(oos_sharpe_ratio, na.rm = TRUE),
      mean_oos_max_drawdown = mean(oos_max_drawdown, na.rm = TRUE),
      mean_oos_cdaR = mean(oos_cdaR, na.rm = TRUE),
      mean_oos_psr = mean(oos_psr, na.rm = TRUE),
      mean_oos_dsr = mean(oos_dsr, na.rm = TRUE)
    )

  # New: Out-of-Sample Metrics Summary for EW Benchmark
  ew_benchmark_oos_summary_metrics <- all_ew_benchmark_oos_metrics_df %>%
    summarise(
      mean_ew_oos_cumulative_return = mean(ew_oos_cumulative_return, na.rm = TRUE),
      mean_ew_oos_sharpe_ratio = mean(ew_oos_sharpe_ratio, na.rm = TRUE),
      mean_ew_oos_max_drawdown = mean(ew_oos_max_drawdown, na.rm = TRUE),
      mean_ew_oos_cdaR = mean(ew_oos_cdaR, na.rm = TRUE),
      mean_ew_oos_psr = mean(ew_oos_psr, na.rm = TRUE),
      mean_ew_oos_dsr = mean(ew_oos_dsr, na.rm = TRUE)
    )
  
  # Plot In-Sample metrics over time
  if (nrow(all_in_sample_metrics_df) > 1) {
    p_cdaR_is <- ggplot(all_in_sample_metrics_df, aes(x = val_date, y = cdaR)) +
      geom_line() + geom_point() +
      labs(title = "In-Sample CDaR (95%) per Validation Window", x = "Validation Date", y = "CDaR") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_in_sample_cdaR_over_time.png"), plot = p_cdaR_is, width = 10, height = 6)

    p_sharpe_is <- ggplot(all_in_sample_metrics_df, aes(x = val_date, y = sharpe_ratio)) +
      geom_line() + geom_point() +
      labs(title = "In-Sample Sharpe Ratio per Validation Window", x = "Validation Date", y = "Sharpe Ratio") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_in_sample_sharpe_over_time.png"), plot = p_sharpe_is, width = 10, height = 6)
  }

  # Plot Out-of-Sample metrics over time for RSA
  if (nrow(all_oos_metrics_df) > 1) {
    p_cdaR_oos <- ggplot(all_oos_metrics_df, aes(x = val_date, y = oos_cdaR)) +
      geom_line() + geom_point() +
      labs(title = "RSA Out-of-Sample CDaR (95%) per Validation Window", x = "Validation Date", y = "CDaR") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_cdaR_over_time_rsa.png"), plot = p_cdaR_oos, width = 10, height = 6)

    p_sharpe_oos <- ggplot(all_oos_metrics_df, aes(x = val_date, y = oos_sharpe_ratio)) +
      geom_line() + geom_point() +
      labs(title = "RSA Out-of-Sample Sharpe Ratio per Validation Window", x = "Validation Date", y = "Sharpe Ratio") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_sharpe_over_time_rsa.png"), plot = p_sharpe_oos, width = 10, height = 6)

    p_psr_oos <- ggplot(all_oos_metrics_df, aes(x = val_date, y = oos_psr)) +
      geom_line() + geom_point() +
      labs(title = "RSA Out-of-Sample Probabilistic Sharpe Ratio per Validation Window", x = "Validation Date", y = "PSR") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_psr_over_time_rsa.png"), plot = p_psr_oos, width = 10, height = 6)

    p_dsr_oos <- ggplot(all_oos_metrics_df, aes(x = val_date, y = oos_dsr)) +
      geom_line() + geom_point() +
      labs(title = "RSA Out-of-Sample Deflated Sharpe Ratio per Validation Window", x = "Validation Date", y = "DSR") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_dsr_over_time_rsa.png"), plot = p_dsr_oos, width = 10, height = 6)
  }

  # New: Plot Out-of-Sample metrics over time for EW Benchmark
  if (nrow(all_ew_benchmark_oos_metrics_df) > 1) {
    p_ew_cdaR_oos <- ggplot(all_ew_benchmark_oos_metrics_df, aes(x = val_date, y = ew_oos_cdaR)) +
      geom_line() + geom_point() +
      labs(title = "EW Benchmark OOS CDaR (95%) per Validation Window", x = "Validation Date", y = "CDaR") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_cdaR_over_time_ew.png"), plot = p_ew_cdaR_oos, width = 10, height = 6)

    p_ew_sharpe_oos <- ggplot(all_ew_benchmark_oos_metrics_df, aes(x = val_date, y = ew_oos_sharpe_ratio)) +
      geom_line() + geom_point() +
      labs(title = "EW Benchmark OOS Sharpe Ratio per Validation Window", x = "Validation Date", y = "Sharpe Ratio") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_sharpe_over_time_ew.png"), plot = p_ew_sharpe_oos, width = 10, height = 6)

    p_ew_psr_oos <- ggplot(all_ew_benchmark_oos_metrics_df, aes(x = val_date, y = ew_oos_psr)) +
      geom_line() + geom_point() +
      labs(title = "EW Benchmark OOS Probabilistic Sharpe Ratio per Validation Window", x = "Validation Date", y = "PSR") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_psr_over_time_ew.png"), plot = p_ew_psr_oos, width = 10, height = 6)

    p_ew_dsr_oos <- ggplot(all_ew_benchmark_oos_metrics_df, aes(x = val_date, y = ew_oos_dsr)) +
      geom_line() + geom_point() +
      labs(title = "EW Benchmark OOS Deflated Sharpe Ratio per Validation Window", x = "Validation Date", y = "DSR") +
      theme_minimal()
    ggsave(file.path(output_plots_dir, "global_oos_dsr_over_time_ew.png"), plot = p_ew_dsr_oos, width = 10, height = 6)
  }

  # Combined Out-of-Sample Equity Curve (RSA vs. Benchmarks)
  if (length(all_oos_returns_list) > 0 || length(all_ew_benchmark_oos_returns_list) > 0) {
    combined_rsa_oos_returns <- do.call(rbind, all_oos_returns_list)
    combined_rsa_oos_returns <- combined_rsa_oos_returns[!duplicated(index(combined_rsa_oos_returns))]
    
    combined_ew_oos_returns <- do.call(rbind, all_ew_benchmark_oos_returns_list)
    combined_ew_oos_returns <- combined_ew_oos_returns[!duplicated(index(combined_ew_oos_returns))]

    # Only proceed if there is valid data for at least one series
    if (NROW(combined_rsa_oos_returns) > 0 || NROW(combined_ew_oos_returns) > 0) {
      
      # Merge returns and compute equity curves
      merged_returns <- merge(combined_rsa_oos_returns, combined_ew_oos_returns)
      colnames(merged_returns) <- c("RSA", "EW_Benchmark")

      # Remove NAs at the beginning due to different start dates (if any)
      merged_returns <- na.omit(merged_returns)

      if (NROW(merged_returns) > 0) {
        # Compute cumulative equity curve for all series
        cumulative_equity <- PerformanceAnalytics::Return.cumulative(merged_returns, geometric = TRUE) + 1
        
        # Add initial value (1) to the equity curve
        initial_date <- min(index(merged_returns)) - days(1)
        initial_value <- xts(matrix(1, nrow = 1, ncol = ncol(cumulative_equity)), order.by = initial_date)
        colnames(initial_value) <- colnames(cumulative_equity)
        
        full_equity_curve <- rbind(initial_value, cumulative_equity)

        equity_df <- data.frame(Date = index(full_equity_curve), coredata(full_equity_curve)) %>%
          pivot_longer(cols = -Date, names_to = "Strategy", values_to = "Value")
        
        p_equity <- ggplot(equity_df, aes(x = Date, y = Value, color = Strategy)) +
          geom_line() +
          labs(title = "Combined Out-of-Sample Equity Curve (RSA vs. Benchmarks)", y = "Cumulative Value", x = "Date") +
          theme_minimal()
        ggsave(file.path(output_plots_dir, "global_oos_equity_curve_comparison.png"), plot = p_equity, width = 12, height = 7)
      }
    }
  }
  
  # Prepare global report content
  report_content <- list(
    Title = "RSA Global Walk-Forward Validation Report",
    Timestamp = Sys.time(),
    In_Sample_Summary_Metrics = as.list(in_sample_summary_metrics),
    Out_of_Sample_Summary_Metrics = as.list(oos_summary_metrics),
    Out_of_Sample_Benchmark_Summary_Metrics = as.list(ew_benchmark_oos_summary_metrics), # New
    In_Sample_Metrics_Per_Window = all_in_sample_metrics_df %>% arrange(val_date) %>% as.list(),
    Out_of_Sample_Metrics_Per_Window = all_oos_metrics_df %>% arrange(val_date) %>% as.list(),
    Out_of_Sample_Benchmark_Metrics_Per_Window = all_ew_benchmark_oos_metrics_df %>% arrange(val_date) %>% as.list() # New
  )

  # Calculate and add Relative Performance Metrics (Alpha, Information Ratio)
  # Requires combined OOS returns for RSA and Benchmark
  if (length(all_oos_returns_list) > 0 && length(all_ew_benchmark_oos_returns_list) > 0) {
    combined_rsa_oos_returns <- do.call(rbind, all_oos_returns_list)
    combined_rsa_oos_returns <- combined_rsa_oos_returns[!duplicated(index(combined_rsa_oos_returns))]
    
    combined_ew_oos_returns <- do.call(rbind, all_ew_benchmark_oos_returns_list)
    combined_ew_oos_returns <- combined_ew_oos_returns[!duplicated(index(combined_ew_oos_returns))]

    # Merge returns to ensure alignment for relative performance calculations
    aligned_returns <- merge(combined_rsa_oos_returns, combined_ew_oos_returns)
    colnames(aligned_returns) <- c("RSA_Returns", "EW_Benchmark_Returns")
    aligned_returns <- na.omit(aligned_returns) # Remove NAs due to differing start/end

    if (NROW(aligned_returns) > 0) {
      # Assuming 0 risk-free rate
      # Alpha: Jensen's Alpha against the EW Benchmark
      # Needs a CAPM-like regression: Excess_RSA = alpha + beta * Excess_Benchmark
      # Or, a simpler alpha: mean(RSA_Returns - EW_Benchmark_Returns)
      
      # For now, let's calculate a simple active return (difference) and its information ratio
      active_returns <- aligned_returns$RSA_Returns - aligned_returns$EW_Benchmark_Returns
      
      relative_metrics <- list(
        mean_active_return = mean(active_returns, na.rm = TRUE),
        information_ratio = PerformanceAnalytics::InformationRatio(active_returns, benchmark_returns = aligned_returns$EW_Benchmark_Returns)
      )
      report_content$Relative_Performance_Metrics <- relative_metrics
    }
  }
  
  
  # Save to a text file
  report_file_path <- file.path(output_report_dir, "rsa_global_validation_report.txt")
  yaml::write_yaml(report_content, file = report_file_path)
  
  message("Global validation report saved to: ", report_file_path)
  return(report_file_path) # Return file path for targets
}
