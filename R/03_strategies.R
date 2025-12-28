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
#' @description Calculates portfolio weights that minimize Conditional Drawdown at Risk (CDaR).
#' Uses `PortfolioAnalytics` and a linear programming solver to achieve this.
#' @param asset_returns An xts object of asset returns.
#' @param portfolio_name A character string for the portfolio name.
#' @param beta_cdar Numeric, the beta level for CDaR (e.g., 0.05 for 95% CDaR).
#' @return A vector of portfolio weights.
calculate_min_cdar_weights <- function(asset_returns, app_config) {
  beta_cdar <- app_config$default$models$cdar_beta
  
  portf <- portfolio.spec(assets = colnames(asset_returns))
  portf <- add.constraint(portf, type = "full_investment")
  portf <- add.constraint(portf, type = "long_only")
  portf <- add.objective(portf, type = "risk", name = "CDaR", arguments = list(p = beta_cdar))
  
  # Preprocess: remove rows with NA and ensure enough data
  asset_returns_clean <- tryCatch({
    stats::na.omit(asset_returns)
  }, error = function(e) asset_returns)
  if (NROW(asset_returns_clean) < 2 || NCOL(asset_returns_clean) < 1) {
    warning("Insufficient data after na.omit for Min-CDaR; returning equal weights.")
    return(get_equal_weight_fallback(asset_returns, "Min-CDaR"))
  }

  # Try methods in order: ROI, GLPK (if available), then randomized search as fallback
  result <- try(optimize.portfolio(R = asset_returns_clean, portfolio = portf, optimize_method = "ROI", trace = FALSE), silent = TRUE)

  if (inherits(result, "try-error")) {
    roi_msg <- as.character(result)
    warning("Min-CDaR ROI attempt failed: ", roi_msg)

    # Try GLPK if available (preferred LP solver)
    if (requireNamespace("Rglpk", quietly = TRUE) || requireNamespace("ROI.plugin.glpk", quietly = TRUE)) {
      warning("Attempting Min-CDaR with GLPK solver (optimize_method = 'glpk').")
      result_glpk <- try(optimize.portfolio(R = asset_returns_clean, portfolio = portf, optimize_method = "glpk", trace = FALSE), silent = TRUE)
      if (!inherits(result_glpk, "try-error")) {
        weights <- PortfolioAnalytics::extractWeights(result_glpk)
        return(as.numeric(weights))
      } else {
        warning("Min-CDaR GLPK attempt failed: ", as.character(result_glpk))
      }
    } else {
      warning("GLPK not available; skipping GLPK attempt.")
    }

    # Fallback: randomized search to get an approximate solution
    random_search_size <- tryCatch({
      app_config$default$models$mincdar_random_search_size
    }, error = function(e) 2000)
    warning("Falling back to randomized search for Min-CDaR (search_size=", random_search_size, ").")
    result_rand <- try(optimize.portfolio(R = asset_returns_clean, portfolio = portf, optimize_method = "random", search_size = as.integer(random_search_size), trace = FALSE), silent = TRUE)
    if (!inherits(result_rand, "try-error")) {
      weights <- PortfolioAnalytics::extractWeights(result_rand)
      return(as.numeric(weights))
    }

    warning("All Min-CDaR optimization attempts failed; returning equal weights.")
    return(get_equal_weight_fallback(asset_returns, "Min-CDaR"))
  } else {
    weights <- PortfolioAnalytics::extractWeights(result)
    return(as.numeric(weights))
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
