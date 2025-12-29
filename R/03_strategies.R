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
#' @description Creates an equal-weighted portfolio as a fallback.
#' @param R An xts object of asset returns.
#' @param strategy_name The name of the strategy that failed, for a more informative warning.
#' @return A named numeric vector of equal weights.
get_equal_weight_fallback <- function(R, strategy_name) {
  warning("Error calculating ", strategy_name, " weights. Returning equal weights as fallback.")
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
#' @param portfolio_name A character string for the portfolio name.
#' @param alpha_cvar Numeric, the alpha level for CVaR (e.g., 0.05 for 95% CVaR).
#' @return A vector of portfolio weights.
calculate_min_cvar_weights <- function(asset_returns, app_config) {
  alpha_cvar <- app_config$default$models$cvar_alpha
  portf <- portfolio.spec(assets = colnames(asset_returns))
  portf <- add.constraint(portf, type = "full_investment")
  portf <- add.constraint(portf, type = "long_only")
  portf <- add.objective(portf, type = "risk", name = "CVaR", arguments = list(p = alpha_cvar))
  
  result <- try(optimize.portfolio(R = asset_returns, portfolio = portf, optimize_method = "ROI", trace = FALSE), silent = TRUE)

  if (inherits(result, "try-error")) {
    return(get_equal_weight_fallback(asset_returns, "Min-CVaR"))
  } else {
    weights <- PortfolioAnalytics::extractWeights(result)
    return(as.numeric(weights))
  }
}


#' @title Calculate Minimum CDaR (Conditional Drawdown at Risk) Weights
#' @description Calculates portfolio weights that minimize Conditional Drawdown at Risk (CDaR)
#' using a custom CVXR implementation based on the problem's convex formulation.
#' @param asset_returns An xts object of asset returns.
#' @param app_config A list containing application configuration.
#' @return A vector of portfolio weights.
calculate_min_cdar_weights <- function(asset_returns, app_config) {
  # Preprocess: remove rows with NA and ensure enough data
  asset_returns_clean <- tryCatch({
    stats::na.omit(asset_returns)
  }, error = function(e) asset_returns)
  if (NROW(asset_returns_clean) < 2 || NCOL(asset_returns_clean) < 1) {
    warning("Insufficient data after na.omit for Min-CDaR; returning equal weights.")
    return(get_equal_weight_fallback(asset_returns, "Min-CDaR"))
  }

  asset_returns_mat <- as.matrix(asset_returns_clean)
  T_obs <- nrow(asset_returns_mat)
  N_assets <- ncol(asset_returns_mat)
  
  # Retrieve CDaR alpha from config, with a fallback
  alpha_cdar <- app_config$default$models$cdar_beta %||% 0.95

  # Define CVXR variables
  w <- CVXR::Variable(N_assets)       # Portfolio weights
  u <- CVXR::Variable(T_obs)          # Drawdown at each time point
  zeta <- CVXR::Variable(1)           # CDaR threshold (like VaR for CVaR)
  z <- CVXR::Variable(T_obs)          # Positive deviation from the threshold

  # Portfolio returns at each time point (as an expression)
  portfolio_returns <- asset_returns_mat %*% w

  # Objective function
  objective <- CVXR::Minimize(zeta + (1 / (T_obs * (1 - alpha_cdar))) * sum(z))

  # Base constraints
  constraints <- list(
    sum(w) == 1,
    w >= 0,
    z >= 0,
    z >= u - zeta,
    u >= 0
  )
  
  # Add drawdown constraints: u_t >= u_{t-1} - r_t'w, with u_0 = 0
  # This is a standard LP formulation for drawdown.
  drawdown_constraints <- list()
  if (T_obs > 0) {
    drawdown_constraints[[1]] <- u[1] >= -portfolio_returns[1]
    if (T_obs > 1) {
      for (t in 2:T_obs) {
        # The drawdown at time t is at least the drawdown at t-1 minus the portfolio return at t
        drawdown_constraints[[t]] <- u[t] >= u[t-1] - portfolio_returns[t]
      }
    }
  }

  all_constraints <- c(constraints, drawdown_constraints)

  # Problem definition and solving
  problem <- CVXR::Problem(objective, all_constraints)
  
  # Use the safe solver utility which should be available
  solve_result <- tryCatch({
    safe_solve_cvxr(problem, allow_fallback = TRUE, verbose = TRUE)
  }, error = function(e) {
    list(result = NULL, error = paste("CVXR solver failed:", e$message))
  })

  # Check if solver failed or did not produce an optimal result
  if (!is.null(solve_result$error) || is.null(solve_result$result)) {
    warning(paste("Custom Min-CDaR optimization failed with CVXR:", solve_result$error))
    return(get_equal_weight_fallback(asset_returns, "Min-CDaR (CVXR)"))
  }

  cvxr_result <- solve_result$result
  if (tolower(cvxr_result$status) %in% c("optimal", "optimal_inaccurate")) {
    optimized_weights <- as.numeric(cvxr_result$getValue(w))
    names(optimized_weights) <- colnames(asset_returns)
    return(optimized_weights)
  } else {
    warning(paste("Custom Min-CDaR with CVXR did not find an optimal solution. Status:", cvxr_result$status))
    return(get_equal_weight_fallback(asset_returns, "Min-CDaR (CVXR)"))
  }
}


#' @title Calculate Hierarchical Risk Parity (HRP) Weights
#' @description Calculates portfolio weights using Hierarchical Risk Parity.
#' @param asset_returns An xts object of asset returns.
#' @return A vector of portfolio weights.
calculate_hrp_weights <- function(asset_returns, app_config) {
  cov_mat <- cov(asset_returns)
  hrp_linkage_method <- app_config$default$models$hrp_linkage_method %||% "single"

  hrp_result <- try(HierPortfolios::HRP_Portfolio(covar = cov_mat, linkage = hrp_linkage_method, graph = FALSE), silent = TRUE)

  if (inherits(hrp_result, "try-error") || !("weights" %in% names(hrp_result))) {
    return(get_equal_weight_fallback(asset_returns, "HRP"))
  } else {
    weights <- hrp_result$weights
    if (is.matrix(weights)) {
      weights <- as.vector(weights)
    }
    names(weights) <- colnames(asset_returns)
    return(weights)
  }
}

#' @title Calculate Equal Risk Contribution (ERC) Weights
#' @description Calculates portfolio weights for Equal Risk Contribution.
#' @param asset_returns An xts object of asset returns.
#' @return A vector of portfolio weights.
calculate_erc_weights <- function(asset_returns) {
  cov_mat <- cov(asset_returns)
  n <- ncol(asset_returns)
  bvec <- rep(1 / n, n)
  result <- try(riskParityPortfolio::riskParityPortfolio(Sigma = cov_mat, b = bvec), silent = TRUE)

  if (inherits(result, "try-error")) {
    return(get_equal_weight_fallback(asset_returns, "ERC"))
  } else {
    weights <- as.vector(result$w)
    names(weights) <- colnames(asset_returns)
    return(weights)
  }
}

#' @title Calculate Robust Mean-Variance (Max Sharpe) Weights
#' @description Calculates portfolio weights that maximize the Sharpe ratio,
#' using a Ledoit-Wolf shrinkage covariance matrix.
#' @param asset_returns An xts object of asset returns.
#' @param portfolio_name A character string for the portfolio name.
#' @return A vector of portfolio weights.
calculate_max_sharpe_weights <- function(asset_returns, portfolio_name = "MaxSharpe_LW_Portfolio") {
  cov_lw <- calculate_lw_covariance(asset_returns)
  mu <- colMeans(asset_returns)
  
  portf <- portfolio.spec(assets = colnames(asset_returns))
  portf <- add.constraint(portf, type = "full_investment")
  portf <- add.constraint(portf, type = "long_only")
  portf <- add.objective(portf, type = "return", name = "mean")
  portf <- add.objective(portf, type = "risk", name = "StdDev")
  
  result <- try(optimize.portfolio(
    R = asset_returns,
    portfolio = portf,
    optimize_method = "random",
    search_size = 2000,
    trace = FALSE,
    mu = mu,
    sigma = cov_lw
  ), silent = TRUE)

  if (inherits(result, "try-error")) {
    return(get_equal_weight_fallback(asset_returns, "Max-Sharpe LW"))
  } else {
    weights <- PortfolioAnalytics::extractWeights(result)
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
  strategy_weights[["ERC"]] <- calculate_erc_weights(asset_returns)
  
  # Strategy 5: Robust Max-Sharpe (Ledoit-Wolf)
  strategy_weights[["MaxSharpeLW"]] <- calculate_max_sharpe_weights(asset_returns)
  
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
