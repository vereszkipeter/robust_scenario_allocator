options(xts.warn_dplyr_breaks_lag = FALSE)

library(ggplot2)
library(patchwork) # For combining ggplot2 plots
library(xts)       # For time series operations like index()
library(forecast)  # For Acf and Pacf plotting functions
library(PerformanceAnalytics) # For PSR and DSR calculations

#' @title Calculate Log Returns
#' @description Calculates the log returns of a financial time series using
#' the PerformanceAnalytics package for robustness.
#' @param x An xts object, typically of adjusted prices.
#' @return An xts object of log returns.
calculate_log_return <- function(x) {
  # NOTE: This function is not currently used in asset_metadata for transformations.
  # Instead, `base::log` is used for the VIXCLS series.
  # Wrapper around PerformanceAnalytics for consistent, robust calculation
  PerformanceAnalytics::Return.calculate(x, method = "log")
}

#' @title Inverse of Log Returns (Simple Returns)
#' @description Converts log returns back to simple returns (factor).
#' @param x An xts object of log returns.
#' @return An xts object of simple returns.
exp_return <- function(x) {
  exp(x)
}

#' @title Calculate Month-over-Month Change (Simple Return)
#' @description Calculates the month-over-month percentage change for a series.
#' @param x An xts object.
#' @return An xts object of month-over-month changes.
mom_change <- function(x) {
  # Wrapper around PerformanceAnalytics for consistent, robust calculation
  PerformanceAnalytics::Return.calculate(x, method = "simple")
}

#' @title Inverse of Month-over-Month Change
#' @description Reconstructs the level of a series from its month-over-month changes.
#' Note: This function provides the period-over-period multiplier. For full series
#' reconstruction, it needs to be applied cumulatively with a starting value.
#' @param x An xts object of mom changes.
#' @return An xts object of period-over-period multipliers (1 + mom_change).
undo_mom_change <- function(x) {
  1 + x
}

#' @title Calculate Simple Returns
#' @description Calculates the simple returns of a financial time series.
#' @param x An xts object, typically of adjusted prices.
#' @return An xts object of simple returns.
simple_return <- function(x) {
  PerformanceAnalytics::Return.calculate(x, method = "simple")
}

#' @title Inverse of Simple Returns (Factor)
#' @description Converts simple returns back to a multiplier (1 + return).
#' @param x An xts object of simple returns.
#' @return An xts object of multipliers (1 + simple_return).
undo_simple_return <- function(x) {
  1 + x
}


#' @title Identity Transformation
#' @description Returns the input series unchanged. Useful for data that doesn't
#' require transformation (e.g., interest rates already in percentage form).
#' @param x An xts object.
#' @return The input xts object.
identity_transform <- function(x) {
  x
}

#' @title Inverse Identity Transformation
#' @description Returns the input series unchanged. The inverse of identity_transform.
#' @param x An xts object.
#' @return The input xts object.
identity_transform_inverse <- function(x) {
  x
}


#' @title Perform Data Sanity Check and Visualization
#' @description This function performs a sanity check on the processed monthly returns
#'   data. It generates summary statistics, checks for NAs, and creates several
#'   plots to visualize the data quality and characteristics, saving them to a
#'   specified directory. It now includes ACF and PACF plots.
#' @param monthly_returns An xts object of multivariate monthly asset returns.
#' @param output_plots_dir The directory path where plots should be saved.
#' @return The file path to a dummy indicator file, signifying completion.
perform_data_sanity_check <- function(monthly_returns, output_plots_dir) {
  
  # Ensure output directory exists
  if (!dir.exists(output_plots_dir)) {
    dir.create(output_plots_dir, recursive = TRUE)
  }
  
  # 1. Summary Statistics
  print("Summary Statistics for Monthly Returns:")
  print(summary(monthly_returns))
  
  # 2. NA Check
  print("NA Count per Asset:")
  na_counts <- colSums(is.na(monthly_returns))
  print(na_counts)
  if (any(na_counts > 0)) {
    warning("NA values detected in monthly returns. `na.locf` may have been insufficient or new NAs introduced.")
  }
  
  # Convert xts to data.frame for ggplot2
  returns_df <- data.frame(
    Date = index(monthly_returns),
    coredata(monthly_returns)
  ) %>%
    pivot_longer(
      cols = -Date,
      names_to = "Asset",
      values_to = "Returns"
    )
  
  # 3. Generate all plots by calling the new dedicated functions
  plot_returns_time_series(returns_df, output_plots_dir)
  plot_returns_distributions(returns_df, output_plots_dir)
  plot_returns_boxplots(returns_df, output_plots_dir)
  plot_returns_acf(returns_df, output_plots_dir)
  plot_returns_pacf(returns_df, output_plots_dir)
  
  # Create a dummy file to satisfy targets' format = "file"
  output_file <- file.path(output_plots_dir, "data_sanity_check_report_complete.txt")
  writeLines(c("Data sanity check completed.", 
               paste("Plots saved to:", output_plots_dir)), 
             output_file)
  
  print(paste("Data sanity check report generated:", output_file))
  
  return(output_file)
}


#' @title Calculate Probabilistic Sharpe Ratio (PSR) and Deflated Sharpe Ratio (DSR)
#' @description Computes PSR and DSR using functions from the `PerformanceAnalytics` package.
#' @param R An xts object of portfolio returns.
#' @param confidence_level Numeric, the confidence level for PSR/DSR (e.g., 0.95).
#' @param n_strategies Numeric, the number of strategies or trials to adjust DSR for (M).
#'   Defaults to 1 for evaluating a single portfolio.
#' @return A list containing PSR (`psr`) and DSR (`dsr`).
calculate_psr_dsr <- function(R, confidence_level, n_strategies = 1) {
  
  # Ensure R is an xts object
  if (!xts::is.xts(R)) {
    stop("Input 'R' must be an xts object of returns.")
  }
  
  # Calculate Probabilistic Sharpe Ratio
  psr <- PerformanceAnalytics::SharpeRatio.Probabilistic(
    R,
    Rf = 0, # Assuming 0 risk-free rate as per project convention, or pass explicitly
    MAR = 0, # Minimum Acceptable Return, 0 if Rf is 0
    p = confidence_level,
    func = "StdDev" # Use StdDev for standard Sharpe Ratio input
  )
  
  # Calculate Deflated Sharpe Ratio
  # M is the number of trials or strategies run.
  # For evaluating a single portfolio, M=1 is appropriate.
  # If we were comparing multiple portfolios and picking the best, M would be higher.
  dsr <- PerformanceAnalytics::SharpeRatio.Deflated(
    R,
    Rf = 0, # Assuming 0 risk-free rate
    MAR = 0, # Minimum Acceptable Return
    p = confidence_level,
    M = n_strategies # Number of trials/strategies
  )
  
  return(list(psr = psr, dsr = dsr))
}


#' @title Ensure App Config Defaults
#' @description Populate commonly required nested defaults in `app_config` so
#' functions downstream can safely read values like `risk_aversion_param`,
#' `term_premium_min`, `term_premium_max`, `dro_epsilon`, and
#' `short_rate_variable_name` without failing. Returns the (possibly updated)
#' `app_config` list.
ensure_app_config_defaults <- function(app_config) {
  if (is.null(app_config)) app_config <- list()
  
  # Ensure top-level 'default' section exists
  if (is.null(app_config$default)) app_config$default <- list()

  # Set defaults for top-level keys within 'default' if not already present
  if (is.null(app_config$default$cache_dir)) app_config$default$cache_dir <- "data/cache"
  if (is.null(app_config$default$n_simulations)) app_config$default$n_simulations <- 100
  if (is.null(app_config$default$horizon_months)) app_config$default$horizon_months <- 60
  if (is.null(app_config$default$initial_window_months)) app_config$default$initial_window_months <- 120
  if (is.null(app_config$default$roll_forward_months)) app_config$default$roll_forward_months <- 12
  if (is.null(app_config$default$validation_start_date)) app_config$default$validation_start_date <- "2015-01-01"

  # Ensure 'models' sub-section exists within 'default'
  if (is.null(app_config$default$models)) app_config$default$models <- list()

  # Now populate defaults for keys within 'app_config$default$models'
  mods <- app_config$default$models
  if (is.null(mods$bvar_lags)) mods$bvar_lags <- 2 # Added this default based on config.yml
  if (is.null(mods$rsbvar_n_iter_mcmc)) mods$rsbvar_n_iter_mcmc <- 2000
  if (is.null(mods$rsbvar_n_burnin_mcmc)) mods$rsbvar_n_burnin_mcmc <- 1000
  if (is.null(mods$rsbvar_tau)) mods$rsbvar_tau <- 0.1
  if (is.null(mods$rsbvar_rho)) mods$rsbvar_rho <- 0.5
  if (is.null(mods$rsbvar_M)) mods$rsbvar_M <- 2
  if (is.null(mods$risk_aversion_param)) mods$risk_aversion_param <- 3.0
  if (is.null(mods$short_rate_variable_name)) mods$short_rate_variable_name <- "SHORT_RATE"
  if (is.null(mods$term_premium_min)) mods$term_premium_min <- -0.005
  if (is.null(mods$term_premium_max)) mods$term_premium_max <- 0.015
  if (is.null(mods$dro_epsilon)) mods$dro_epsilon <- 0.05
  if (is.null(mods$entropy_penalty_kappa)) mods$entropy_penalty_kappa <- 0.001
  if (is.null(mods$base_term_premium)) mods$base_term_premium <- 0.002
  if (is.null(mods$tp_sensitivity)) mods$tp_sensitivity <- 0.5
  if (is.null(mods$transaction_cost_flat_rate)) mods$transaction_cost_flat_rate <- 0
  if (is.null(mods$cvar_alpha)) mods$cvar_alpha <- 0.05
  if (is.null(mods$cdar_beta)) mods$cdar_beta <- 0.05 # Added this default based on config.yml
  if (is.null(mods$confidence_level)) mods$confidence_level <- 0.95
  if (is.null(mods$ep_shrinkage_alpha)) mods$ep_shrinkage_alpha <- 0.1
  if (is.null(mods$ep_return_tolerance)) mods$ep_return_tolerance <- 0.001
  if (is.null(mods$hrp_linkage_method)) mods$hrp_linkage_method <- "single" # Default as per HierPortfolios example
  if (is.null(mods$dcc_solver)) mods$dcc_solver <- "solnp" # Added this default
  if (is.null(mods$dcc_solver_control)) mods$dcc_solver_control <- list(trace = 0) # Added this default
  if (is.null(mods$min_dcc_obs)) mods$min_dcc_obs <- 100 # Minimum observations for DCC-GARCH

  app_config$default$models <- mods
  
  # Ensure app_config$default$data is properly handled
  if (is.null(app_config$default$data)) app_config$default$data <- list()
  if (is.null(app_config$default$data$tickers)) app_config$default$data$tickers <- c("SPY", "AGG", "GLD", "QQQ") # Example default
  if (is.null(app_config$default$data$from)) app_config$default$data$from <- "2010-01-01" # Example default
  if (is.null(app_config$default$data$to)) app_config$default$data$to <- Sys.Date() # Example default
  
  return(app_config)
}


#' @title Safe CVXR solver wrapper
#' @description Try solving a CVXR Problem using a list of solvers in order.
#' If all solvers fail and `allow_fallback` is FALSE, this function stops with a
#' clear error so callers can fail fast. If `allow_fallback` is TRUE, returns a
#' list with `result = NULL` and `error` information.
#' @param problem A CVXR Problem object.
#' @param solvers Character vector of solver names to try, in order.
#' @param allow_fallback Logical; if FALSE (default) stop on failure, else return structured error.
#' @return On success, the CVXR result object; on failure and `allow_fallback=TRUE`, a list with `result=NULL` and `error`.
safe_solve_cvxr <- function(problem, solvers = c("ECOS", "SCS"), allow_fallback = FALSE) {
  last_error <- NULL
  last_status <- "unknown_status"
  for (solver in solvers) {
    tryCatch({
      res <- CVXR::solve(problem, solver = solver)
      if (res$status %in% c("optimal", "optimal_inaccurate")) {
        return(res)
      } else {
        last_status <- res$status
        last_error <- paste0("Solver '", solver, "' returned status: ", res$status)
      }
    }, error = function(e) {
      last_error <- paste0("Solver '", solver, "' failed with error: ", e$message)
      last_status <- "solver_error"
    }, warning = function(w) {
      message("Solver '", solver, "' issued warning: ", w$message)
    })
  }

  # If we reach here, all solvers failed or returned non-optimal statuses
  msg <- if (!is.null(last_error)) {
    paste0("All CVXR solvers failed. Last attempt: ", last_error, ". Final status: ", last_status)
  } else {
    paste0("CVXR solve failed with unknown error. Final status: ", last_status)
  }
  
  if (allow_fallback) {
    return(list(result = NULL, error = msg))
  }
  stop(msg)
}


