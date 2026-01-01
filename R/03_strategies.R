library(PortfolioAnalytics) # For Min-CVaR, Max-Sharpe
library(CVXR)               # For Min-CDaR (if custom LP), Entropy Pooling (used later)
library(HierPortfolios)     # For HRP/HERC
library(riskParityPortfolio)# For ERC
library(PerformanceAnalytics)# For general financial functions, potentially CDaR
library(xts)                # For time series objects
library(dplyr)              # For data manipulation
library(MASS)               # For Ledoit-Wolf shrinkage (cov.shrink)
library(rlang)              # For the `%||%` operator

# --- Helper Functions for Strategy Implementations ---

#' @title Get Equal Weight Fallback Weights
#' @description Creates an equal-weighted portfolio as a fallback and logs a warning.
#' @param R An xts object of asset returns.
#' @param strategy_name The name of the strategy that failed.
#' @param app_config A list containing application configuration for logging.
#' @return A named numeric vector of equal weights.
get_equal_weight_fallback <- function(R, strategy_name, app_config) {
  log_message(paste("Error calculating", strategy_name, "weights. Returning equal weights as fallback."), level = "WARN", app_config = app_config)
  n_assets <- ncol(R)
  weights <- rep(1 / n_assets, n_assets)
  names(weights) <- colnames(R)
  return(weights)
}


#' @title Calculate Ledoit-Wolf Shrinkage Covariance Matrix
#' @description Computes the Ledoit-Wolf shrinkage covariance matrix.
#' @param R An xts object of asset returns.
#' @return A matrix representing the Ledoit-Wolf shrinkage covariance.
calculate_lw_covariance <- function(R) {
  # browser()
  cov_lw <- corpcor::cov.shrink(R, verbose = FALSE)
  return(cov_lw)
}

# --- Individual Base Strategy Implementations ---

#' @title Calculate Minimum CVaR Weights
#' @description Calculates portfolio weights that minimize Conditional Value-at-Risk (CVaR).
#' Uses `PortfolioAnalytics` and a linear programming solver.
#' @param asset_returns An xts object of asset returns.
#' @param app_config A list containing application configuration.
#' @return A vector of portfolio weights.
calculate_min_cvar_weights <- function(asset_returns, app_config) {
  log_message("Starting Min-CVaR (PortfolioAnalytics) weight calculation.", level = "DEBUG", app_config = app_config)
  alpha_cvar <- app_config$default$models$cvar_alpha
  portf <- portfolio.spec(assets = colnames(asset_returns))
  portf <- add.constraint(portf, type = "full_investment")
  portf <- add.constraint(portf, type = "long_only")
  portf <- add.objective(portf, type = "risk", name = "CVaR", arguments = list(p = alpha_cvar))

  result <- try(optimize.portfolio(R = asset_returns, portfolio = portf, optimize_method = "ROI", trace = FALSE), silent = TRUE)

  if (inherits(result, "try-error")) {
    error_msg <- paste("Min-CVaR (PortfolioAnalytics) optimization failed. Error:", as.character(result))
    log_message(error_msg, level = "ERROR", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "Min-CVaR", app_config))
  } else {
    weights <- PortfolioAnalytics::extractWeights(result)
    log_message("Successfully calculated Min-CVaR (PortfolioAnalytics) weights.", level = "DEBUG", app_config = app_config)
    return(as.numeric(weights))
  }
}


#' @title Calculate Minimum CVaR Weights using CVXR
#' @description Calculates portfolio weights that minimize Conditional Value-at-Risk (CVaR)
#' using a custom CVXR implementation based on the problem's canonical convex (LP) formulation.
#' This is a replacement for the previous non-DCP compliant CDaR implementation.
#' @param asset_returns An xts object of asset returns.
#' @param app_config A list containing application configuration.
#' @return A vector of portfolio weights.
calculate_min_cdar_weights <- function(asset_returns, app_config) {
  log_message("Starting Min-CVaR (CVXR) weight calculation.", level = "DEBUG", app_config = app_config)
  
  asset_returns_clean <- tryCatch({
    stats::na.omit(asset_returns)
  }, error = function(e) asset_returns)
  
  if (NROW(asset_returns_clean) < 2 || NCOL(asset_returns_clean) < 1) {
    log_message("Insufficient data after na.omit for Min-CVaR (CVXR); returning equal weights.", level = "WARN", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "Min-CVaR (CVXR)", app_config))
  }

  asset_returns_mat <- as.matrix(asset_returns_clean)
  T_obs <- nrow(asset_returns_mat)
  N_assets <- ncol(asset_returns_mat)
  
  alpha_cvar <- app_config$default$models$cvar_alpha %||% 0.95
  gamma_l2 <- app_config$default$models$cdar_l2_gamma %||% 1e-7 # Small regularization

  # --- Define CVXR Variables ---
  w <- CVXR::Variable(N_assets)    # Portfolio weights
  zeta <- CVXR::Variable(1)        # VaR, the variable for the optimization problem
  z <- CVXR::Variable(T_obs)       # Auxiliary variables for losses exceeding VaR

  # --- Define Loss ---
  # The loss is the negative of portfolio return for each observation.
  loss <- - (asset_returns_mat %*% w)

  # --- Objective Function ---
  # Minimize zeta + average of losses exceeding zeta. This is the canonical CVaR formulation.
  # A small L2 regularization on weights is added for stability.
  objective <- CVXR::Minimize(zeta + (1 / (T_obs * (1 - alpha_cvar))) * sum(z) + gamma_l2 * sum_squares(w))

  # --- Constraints ---
  constraints <- list(
    sum(w) == 1,      # Full investment
    w >= 0,           # Long only
    z >= 0,           # Losses z must be non-negative
    z >= loss - zeta  # If loss > zeta, then z must be at least loss - zeta
  )
  
  # Problem definition and solving
  problem <- CVXR::Problem(objective, constraints)
  
  # Use a safe solver wrapper if available, otherwise solve directly
  solve_fn <- if (exists("safe_solve_cvxr")) safe_solve_cvxr else function(p, ...) CVXR::solve(p, ...)
  
  solve_result <- tryCatch({
    solve_fn(problem, allow_fallback = TRUE, verbose = FALSE)
  }, error = function(e) {
    list(result = NULL, error = paste("CVXR solver failed:", e$message))
  })

  # Check solver outcome
  if (!is.null(solve_result$error) || is.null(solve_result$result) || !(tolower(solve_result$result$status) %in% c("optimal", "optimal_inaccurate"))) {
    
    error_msg <- if (!is.null(solve_result$error)) {
      paste("Min-CVaR (CVXR) optimization failed with solver error:", solve_result$error)
    } else if (is.null(solve_result$result)) {
      "Min-CVaR (CVXR) optimization failed: CVXR result is NULL."
    } else {
      paste("Min-CVaR (CVXR) did not find an optimal solution. Status:", solve_result$result$status)
    }
    log_message(error_msg, level = "ERROR", app_config = app_config)
    
    # Save diagnostics for debugging
    # (consider disabling this in production if it's too noisy)
    # diag_dir <- "output/cvar_cvxr_diagnostics"
    # if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
    # timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    # diag_file <- file.path(diag_dir, paste0("cvar_cvxr_fail_", timestamp, ".rds"))
    # saveRDS(list(asset_returns = asset_returns_mat, problem = problem, solve_result = solve_result), file = diag_file)
    # log_message(paste("Saved Min-CVaR (CVXR) diagnostic data to", diag_file), level = "INFO", app_config = app_config)
    
    return(get_equal_weight_fallback(asset_returns, "Min-CVaR (CVXR)", app_config))

  } else {
    # Success
    cvxr_result <- solve_result$result
    optimized_weights <- as.numeric(cvxr_result$getValue(w))
    names(optimized_weights) <- colnames(asset_returns)
    log_message("Successfully calculated Min-CVaR (CVXR) weights.", level = "DEBUG", app_config = app_config)
    return(optimized_weights)
  }
}


#' @title Calculate Hierarchical Risk Parity (HRP) Weights
#' @description Calculates portfolio weights using Hierarchical Risk Parity.
#' @param asset_returns An xts object of asset returns.
#' @param app_config A list containing application configuration.
#' @return A vector of portfolio weights.
calculate_hrp_weights <- function(asset_returns, app_config) {
  log_message("Starting HRP weight calculation.", level = "DEBUG", app_config = app_config)
  
  # HRP is sensitive to NA values in the covariance matrix
  asset_returns_clean <- na.omit(asset_returns)
  if (NROW(asset_returns_clean) < 2) {
    log_message("Insufficient data for HRP after na.omit; returning equal weights.", level = "WARN", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "HRP", app_config))
  }

  cov_mat <- cov(asset_returns_clean)
  hrp_linkage_method <- app_config$default$models$hrp_linkage_method %||% "single"

  hrp_result <- try(HierPortfolios::HRP_Portfolio(covar = cov_mat, linkage = hrp_linkage_method, graph = FALSE), silent = TRUE)

  if (inherits(hrp_result, "try-error") || !("weights" %in% names(hrp_result))) {
    error_msg <- if (inherits(hrp_result, "try-error")) as.character(hrp_result) else "HRP_Portfolio did not return 'weights'."
    log_message(paste("HRP calculation failed:", error_msg), level = "ERROR", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "HRP", app_config))
  } else {
    weights <- hrp_result$weights
    if (is.matrix(weights)) {
      weights <- as.vector(weights)
    }
    names(weights) <- colnames(asset_returns)
    log_message("Successfully calculated HRP weights.", level = "DEBUG", app_config = app_config)
    return(weights)
  }
}

#' @title Calculate Equal Risk Contribution (ERC) Weights
#' @description Calculates portfolio weights for Equal Risk Contribution.
#' @param asset_returns An xts object of asset returns.
#' @param app_config A list containing application configuration.
#' @return A vector of portfolio weights.
calculate_erc_weights <- function(asset_returns, app_config) {
  log_message("Starting ERC weight calculation.", level = "DEBUG", app_config = app_config)

  asset_returns_clean <- na.omit(asset_returns)
  if (NROW(asset_returns_clean) < 2) {
    log_message("Insufficient data for ERC after na.omit; returning equal weights.", level = "WARN", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "ERC", app_config))
  }

  cov_mat <- cov(asset_returns_clean)
  n <- ncol(asset_returns_clean)
  bvec <- rep(1 / n, n)
  result <- try(riskParityPortfolio::riskParityPortfolio(Sigma = cov_mat, b = bvec), silent = TRUE)

  if (inherits(result, "try-error")) {
    error_msg <- paste("ERC calculation failed. Error:", as.character(result))
    log_message(error_msg, level = "ERROR", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "ERC", app_config))
  } else {
    weights <- as.vector(result$w)
    names(weights) <- colnames(asset_returns_clean)
    log_message("Successfully calculated ERC weights.", level = "DEBUG", app_config = app_config)
    return(weights)
  }
}

#' @title Calculate Robust Mean-Variance (Max Sharpe) Weights
#' @description Calculates portfolio weights that maximize the Sharpe ratio,
#' using a Ledoit-Wolf shrinkage covariance matrix.
#' @param asset_returns An xts object of asset returns.
#' @param app_config A list containing application configuration.
#' @param portfolio_name A character string for the portfolio name.
#' @return A vector of portfolio weights.
calculate_max_sharpe_weights <- function(asset_returns, app_config, portfolio_name = "MaxSharpe_LW_Portfolio") {
  log_message("Starting Max-Sharpe (Ledoit-Wolf) weight calculation.", level = "DEBUG", app_config = app_config)

  asset_returns_clean <- na.omit(asset_returns)
  if (NROW(asset_returns_clean) < 2) {
    log_message("Insufficient data for Max-Sharpe after na.omit; returning equal weights.", level = "WARN", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "Max-Sharpe LW", app_config))
  }
  
  cov_lw <- calculate_lw_covariance(asset_returns_clean)
  mu <- colMeans(asset_returns_clean)

  portf <- portfolio.spec(assets = colnames(asset_returns_clean))
  portf <- add.constraint(portf, type = "full_investment")
  portf <- add.constraint(portf, type = "long_only")
  portf <- add.objective(portf, type = "return", name = "mean")
  portf <- add.objective(portf, type = "risk", name = "StdDev")

  result <- try(optimize.portfolio(
    R = asset_returns_clean,
    portfolio = portf,
    optimize_method = "random",
    search_size = 2000,
    trace = FALSE,
    mu = mu,
    sigma = cov_lw
  ), silent = TRUE)

  if (inherits(result, "try-error")) {
    error_msg <- paste("Max-Sharpe (Ledoit-Wolf) optimization failed. Error:", as.character(result))
    log_message(error_msg, level = "ERROR", app_config = app_config)
    return(get_equal_weight_fallback(asset_returns, "Max-Sharpe LW", app_config))
  } else {
    weights <- PortfolioAnalytics::extractWeights(result)
    log_message("Successfully calculated Max-Sharpe (Ledoit-Wolf) weights.", level = "DEBUG", app_config = app_config)
    return(as.numeric(weights))
  }
}


#' @title Calculate Weights for All Base Strategies
#' @description A wrapper function to calculate weights for all specified base strategies.
#' @param asset_returns An xts object of asset returns.
#' @param app_config A list containing application configuration.
#' @return A list where each element is a named vector of weights for a strategy.
calculate_base_strategy_weights <- function(asset_returns, app_config) {
  
  strategy_weights <- list()
  
  # Strategy 1: Min-CVaR
  strategy_weights[["MinCVaR"]] <- calculate_min_cvar_weights(asset_returns, app_config)
  
  # Strategy 2: Min-CDaR
  strategy_weights[["MinCDaR"]] <- calculate_min_cdar_weights(asset_returns, app_config)
  
  # Strategy 3: HRP
  strategy_weights[["HRP"]] <- calculate_hrp_weights(asset_returns, app_config)
  
  # Strategy 4: ERC
  strategy_weights[["ERC"]] <- calculate_erc_weights(asset_returns, app_config = app_config)
  
  # Strategy 5: Robust Max-Sharpe (Ledoit-Wolf)
  strategy_weights[["MaxSharpeLW"]] <- calculate_max_sharpe_weights(asset_returns, app_config)
  
  return(strategy_weights)
}

#' @title Calculate Equal Weight Benchmark Weights
#' @description Calculates equal weights for a given set of assets.
#' @param asset_names A character vector of asset tickers.
#' @return A named numeric vector of equal weights.
calculate_equal_weight_benchmark <- function(asset_names) {
  n_assets <- length(asset_names)
  weights <- rep(1 / n_assets, n_assets)
  names(weights) <- asset_names
  return(weights)
}
