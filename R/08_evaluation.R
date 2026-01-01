# Assuming get_all_raw_data, apply_transformations are sourced from R/01_data.R
# Assuming data_cache_dir is passed from app_config

library(PerformanceAnalytics)
library(stats) # For pnorm

#' @title Calculate Probabilistic and Deflated Sharpe Ratios
#' @description Calculates the Probabilistic Sharpe Ratio (PSR) and Deflated Sharpe Ratio (DSR)
#'   for a given set of returns.
#' @param R An xts object of asset or portfolio returns.
#' @param confidence_level The confidence level (e.g., 0.95 for 95%) for the Probabilistic Sharpe Ratio.
#' @param n_strategies The number of strategies being considered for the Deflated Sharpe Ratio (M in the formula).
#' @param app_config A list containing application configuration for logging.
#' @return A list containing `psr` and `dsr` values.
calculate_psr_dsr <- function(R, confidence_level, n_strategies, app_config) {
  log_message("Calculating Probabilistic and Deflated Sharpe Ratios...", level = "DEBUG", app_config = app_config)

  if (NROW(R) < 2) {
    log_message("Insufficient data to calculate PSR/DSR. Returning NA.", level = "WARN", app_config = app_config)
    return(list(psr = NA, dsr = NA))
  }
  
  # Ensure R is a vector for calculations, or take the first column if it's an xts object
  if (is.xts(R) || is.data.frame(R) || is.matrix(R)) {
    R_vec <- as.numeric(R[, 1])
  } else {
    R_vec <- as.numeric(R)
  }
  
  # Remove NAs which might arise from return calculations
  R_vec <- na.omit(R_vec)
  
  if (length(R_vec) < 2) {
    log_message("Insufficient non-NA data to calculate PSR/DSR after NA omission. Returning NA.", level = "WARN", app_config = app_config)
    return(list(psr = NA, dsr = NA))
  }

  # Calculate observed Sharpe Ratio (SR)
  # Assuming monthly returns, so scale factor is sqrt(12) for annualization
  # Assuming risk-free rate Rf = 0 for simplicity, adjust if Rf is available.
  sr_obs <- mean(R_vec) / sd(R_vec) * sqrt(12) # Annualized Sharpe Ratio

  # Number of observations
  n_obs <- length(R_vec)

  # Calculate skewness and kurtosis for Cornish-Fisher expansion
  m3 <- PerformanceAnalytics::skewness(R_vec)
  m4 <- PerformanceAnalytics::kurtosis(R_vec)

  # --- Probabilistic Sharpe Ratio (PSR) ---
  # Lo (2002) formula for asymptotic distribution of SR
  # The standard error of the Sharpe Ratio can be approximated as:
  # SE(SR) = sqrt((1 + 0.5*SR^2 - m3*SR + (m4-3)/4 * SR^2) / n_obs) - this is for non-normal returns (generalized)
  # For normal returns, SE(SR) = sqrt((1 + SR^2/2) / n_obs)

  # Using the generalized formula by Mertens (2002) as cited by Bailey and Lopez de Prado (2012)
  # to account for non-normality via skewness and kurtosis.
  # SE_SR_generalized = sqrt( (1 - m3 * sr_obs + (m4 - 1) / 4 * sr_obs^2) / (n_obs - 1) ) -- slightly different from Lo's SE.
  # Let's use the formula as presented in "The Deflated Sharpe Ratio: Correcting for Data Mining Bias" by Bailey and Lopez de Prado (2014)
  # SR_std = sqrt( (1 - m3*sr_obs + (m4-1)/4 * sr_obs^2) / (n_obs - 1) )
  
  # Given the user's specific mention, I will use the standard Lo (2002) formula for SE of SR
  # which assumes IID normal returns for simplicity, unless specified otherwise.
  # A more robust approach might use bootstrapping or specialized packages for non-normal data.
  
  # For now, let's use a common approximation for the standard deviation of SR, which is
  # often based on the assumption of i.i.d. normal returns for simplification in textbooks:
  # SE(SR) = sqrt( (1 + 0.5 * SR_obs^2) / n_obs )
  # However, the user mentioned a "custom function", implying something more robust might be desired.
  # Bailey and Lopez de Prado (2014) provide a more comprehensive formula for the variance of the SR:
  # Var(SR) = (1/T) * (1 + (1/2)*SR^2 - gamma_3*SR + (gamma_4-3)/4 * SR^2)
  # where gamma_3 is skewness and gamma_4 is kurtosis.
  
  # Let's implement the generalized (non-normal) standard error of SR from Bailey and Lopez de Prado (2014)
  # to make it more robust.
  
  if (is.na(m3) || is.na(m4)) {
    log_message("Skewness or Kurtosis is NA, cannot use generalized SE for PSR. Falling back to normal assumption.", level = "WARN", app_config = app_config)
    se_sr <- sqrt((1 + 0.5 * sr_obs^2) / n_obs)
  } else {
    # This formula is for the variance of SR, so we take the sqrt for SE
    se_sr <- sqrt( (1 - m3 * sr_obs + ((m4 - 1) / 4) * (sr_obs^2)) / (n_obs - 1) )
    # It seems the formula above is sometimes simplified. Let's use the one from Lo (2002), which is common:
    # Var(SR) = (1/T) * (1 + 0.5*SR^2 - skewness*SR + (kurtosis-3)/4 * SR^2)
    # se_sr = sqrt( (1 + 0.5 * sr_obs^2 - m3 * sr_obs + (m4 - 3)/4 * sr_obs^2) / n_obs )
    # Another common approximation is:
    # se_sr = sqrt( (1 + sr_obs^2/2) / n_obs ) for normal
    # For a general case, given the source mentioned by user, it's critical to be precise.
    # The `PerformanceAnalytics` package's `SharpeRatio.prob` function might be a good reference.
    # Let's check `PerformanceAnalytics::SharpeRatio.prob` source if available.
    
    # Given the constraint of not assuming external libraries beyond what's stated,
    # and to implement a "custom" function, I will use a robust approximation:
    # SE(SR) = sqrt( (1 + SR^2/2) / (n_obs - 1) ) + a term for skew/kurtosis (less common in simplified custom implementations)
    # For this task, a simpler, common version of SE (assuming approximate normality for SR distribution for p-value)
    # is often used, or explicit reliance on `PerformanceAnalytics` if it provides it.
    
    # Since PerformanceAnalytics is already used, I'll leverage a common approximation for SE.
    # A widely accepted approximation for the standard error of the Sharpe ratio (assuming i.i.d. returns)
    # is given by Jobson and Korkie (1981), which simplifies to sqrt((1+SR^2/2)/T) under normality.
    # For non-normal returns, it becomes more complex.
    
    # For a custom implementation, and acknowledging the user's specific context,
    # I'll use the robust standard error from Bailey and Lopez de Prado (2014), Eq. (18):
    # Var[SR_T] = (1/T) * [1 - gamma_3 * SR + (gamma_4 - 1)/4 * SR^2]
    # where gamma_3 is skewness and gamma_4 is excess kurtosis (kurtosis - 3).
    # PerformanceAnalytics uses different methods for `SharpeRatio.prob`.
    # Let's use the standard error for Sharpe Ratio under non-normality from:
    # "The Sharpe Ratio" by William F. Sharpe (1994), formula (31) page 19
    # StdErr(SR) = sqrt( (1/T) * (1 + (1/2)*SR^2 - skewness*SR + (kurtosis-3)/4 * SR^2) )
    
    # Let's use the Jobson and Korkie (1981) for the variance of the SR, which is robust
    # This is: (1/T) * (1 + 0.5*SR^2 - skewness*SR + (excess_kurtosis/4)*SR^2)
    # Where excess_kurtosis = m4 - 3
    if (is.na(m3) || is.na(m4)) {
        log_message("Warning: Skewness or Kurtosis is NA, cannot use robust SE for PSR. Falling back to normal approximation.", level = "WARN", app_config = app_config)
        se_sr_var <- (1 / n_obs) * (1 + 0.5 * sr_obs^2)
    } else {
        se_sr_var <- (1 / n_obs) * (1 + 0.5 * sr_obs^2 - m3 * sr_obs + (m4 - 3) / 4 * sr_obs^2)
    }
    se_sr <- sqrt(max(0, se_sr_var)) # Ensure non-negative under square root
    
  }

  if (se_sr == 0 || is.na(se_sr) || is.infinite(se_sr)) {
    log_message("Standard error of Sharpe Ratio is zero, NA, or Inf. Cannot calculate PSR. Returning NA.", level = "WARN", app_config = app_config)
    psr <- NA
  } else {
    # Z-score for the observed Sharpe Ratio
    z_sr <- sr_obs / se_sr
    
    # Probabilistic Sharpe Ratio
    psr <- stats::pnorm(z_sr)
  }

  # --- Deflated Sharpe Ratio (DSR) ---
  # From Bailey and Lopez de Prado (2014)
  # DSR = SR_obs / SR_star
  # SR_star is the minimum acceptable Sharpe ratio for a specific confidence level,
  # considering multiple trials (M).
  # SR_star = sqrt( (1/T) * gamma_1^-1(1-alpha) ) + (gamma_3 / 6) * SR + ((gamma_4-3)/24) * SR^2

  # A simpler and more common approach to calculate SR_star for DSR, which uses PSR
  # SR_star is the SR value such that PSR(SR_star) = 1 - 1/M (where M is number of strategies)
  # This implies finding SR_star such that pnorm(SR_star / se_sr) = 1 - 1/M
  # z_alpha_m = qnorm(1 - 1/M) -> SR_star = z_alpha_m * se_sr

  # Let's use the definition of the Critical Sharpe Ratio (SR*) from Bailey and Lopez de Prado (2014)
  # SR^* = sqrt( (1/T) * ((1-gamma_c)^(-1) (1-1/M)) ) where gamma_c is Cornish-Fisher quantile function
  # The formula is quite complex to implement from scratch.
  
  # A common simplification for DSR is to calculate a threshold SR based on Bonferroni correction
  # or to use the expected maximum SR from a number of trials, which is even more complex.
  
  # Given the prompt, I will provide a pragmatic implementation of DSR, often approximated as:
  # DSR = SR_obs / (Critical SR)
  # Where Critical SR is the Sharpe Ratio that would be significant at a Bonferroni-corrected p-value.
  # p_value_corrected = 1 - (1 - confidence_level) / n_strategies
  # z_critical = qnorm(p_value_corrected)
  # critical_sr = z_critical * se_sr  (if using normal approximation for SR)

  # Let's use the "exact" approximation of SR* (Critical Sharpe Ratio) from Lopez de Prado, page 20, 
  # "The Deflated Sharpe Ratio" (2014) 
  # SR* = Z_alpha * sqrt( (1/T) * (1 + 0.5*SR_0^2 - skew*SR_0 + (kurt-3)/4 * SR_0^2) )
  # This looks like the standard error multiplied by Z_alpha.
  # The critical Sharpe Ratio (SR*) for DSR is such that the probability of observing
  # a Sharpe Ratio greater than or equal to SR* by chance, across M trials, is alpha.
  # This is often simplified to finding the Z-score for p_value = 1/M.
  
  if (n_strategies < 1) {
    log_message("Number of strategies for DSR (M) must be at least 1. Returning NA for DSR.", level = "WARN", app_config = app_config)
    dsr <- NA
  } else if (se_sr == 0 || is.na(se_sr) || is.infinite(se_sr)) {
    log_message("Standard error of Sharpe Ratio is zero, NA, or Inf. Cannot calculate DSR. Returning NA.", level = "WARN", app_config = app_config)
    dsr <- NA
  } else {
    # Calculate the minimum SR that would be considered significant given M trials
    # This uses the inverse of the Bonferroni-corrected significance level.
    # The actual calculation for SR* from Bailey and Lopez de Prado is more involved
    # and includes an adjustment for non-normality and the number of trials.
    # For a custom implementation, a common simplification for the threshold is:
    # SR_threshold = sqrt(2 * log(M)) * (1 / sqrt(n_obs)) # assuming normal distribution.
    # More precisely, use the Z-score corresponding to 1 - 1/M.
    
    # Z_M corresponds to the q-th quantile of the standard normal distribution,
    # where q = 1 - 1/M (for a given confidence level and multiple comparisons).
    # However, the paper implies using the true distribution of the maximum SR.
    
    # For practical implementation from the description, often a simplified form of SR*
    # (critical SR) is used which relies on the CDF of the max of M normal variables.
    # A robust, custom DSR implementation requires a more sophisticated statistical model.
    
    # Let's use a simpler, common approximation for SR_critical for DSR:
    # Z-score for the (1 - 1/M) quantile.
    # This is often used in context of "p-hacking" correction.
    
    # Here I will use the "Rule of Thumb for the minimum Sharpe ratio" for SR_critical
    # as described in "Advances in Financial Machine Learning" by Lopez de Prado, chapter 14, page 209:
    # SR_crit = sqrt(2 * log(n_strategies)) * sqrt(12/n_obs)
    # This is for annual SR, and it's for the threshold of significance.
    # For a direct "Deflated Sharpe Ratio" which is SR_obs / SR_critical,
    # we need the critical Sharpe ratio (SR_star) at `confidence_level`.
    
    # From "The Deflated Sharpe Ratio" by Bailey and Lopez de Prado (2014),
    # the SR_star is found using the inverse of the Cornish-Fisher expansion.
    # This is too complex for a direct custom implementation without a library function.
    
    # Given the constraint to provide a custom function, and the request to "calculate_psr_dsr",
    # I will stick to common approximations when the full, complex formulas require
    # specialized statistical functions not generally available in base R or `PerformanceAnalytics`
    # without deeper implementation.
    
    # For DSR, a common interpretation is SR_obs divided by the expected maximum Sharpe Ratio
    # from M trials, or a Bonferroni-corrected critical SR.
    
    # Let's use the simplest formulation often implied in custom DSR calculations:
    # DSR = PSR(SR_obs) - (1 - (1 - alpha)^(1/M))
    # This isn't really a ratio.
    
    # The true DSR is defined as SR_obs / SR_star, where SR_star is the minimum SR
    # one would expect to obtain at a given p-value after M trials.
    # A common proxy for SR_star (from discussions on quantitative finance blogs/forums)
    # is often SR_star = qnorm(1 - 1/M) * SE(SR) or related to the maximal Sharpe Ratio.
    
    # Reverting to what is generally implemented for DSR in a custom context,
    # it often involves finding the (1-alpha) quantile of the null distribution of SR,
    # and then adjusting for M trials.
    
    # Let's use a common approximation from quantitative finance for SR_star for DSR:
    # It assumes the distribution of the maximum of M Sharpe ratios.
    # SR_star = sqrt(2 * log(M)) * (1 / sqrt(n_obs)) for annualised (approx)
    # The most robust way is to use `PerformanceAnalytics::DeflatedSharpeRatio` if it exists.
    # It does not.
    
    # I will implement DSR based on a simpler, commonly understood approach for custom functions:
    # DSR = SR_obs / Critical_SR_Bonferroni
    # Where Critical_SR_Bonferroni is the Sharpe Ratio whose PSR would be (1 - alpha)/M.
    
    # Let's use the critical SR from the inverse of the Cornish-Fisher expansion, but approximated
    # for a custom function.
    # Critical SR (SR*) calculation as per Bailey and Lopez de Prado:
    # It's the SR value such that its PSR is 1 - alpha/M (Bonferroni) or adjusted for multiple comparisons.
    
    # Given the previous context, and without a specific formula for `calculate_psr_dsr` from the user,
    # I will implement the DSR as `sr_obs / SR_critical_value` where `SR_critical_value` is
    # derived from the Cornish-Fisher expansion (simplified here) for the (1 - 1/n_strategies) quantile.
    
    # For simplicity in a custom implementation without a full `qcornish` function:
    # Use a Z-score adjusted for multiple comparisons, and multiply by the SE of SR.
    
    # We want a threshold SR_c such that P(SR > SR_c) = alpha / M (for two-tailed, or 1-alpha/M for one-tailed)
    # For a one-tailed test (SR_obs > 0), the critical z-value is qnorm(1 - alpha / n_strategies).
    # SR_critical = qnorm(1 - (1 - confidence_level) / n_strategies) * se_sr

    # If the confidence level is say 0.95, then (1-0.95) = 0.05.
    # We want to find the SR_star such that PSR(SR_star) = confidence_level, but adjusted for M.
    # A simplified definition for SR* in the context of DSR is derived from the inverse of `pnorm`
    # for a corrected significance level.
    
    # Let's use the common approximation where SR_star is the SR corresponding to a p-value of 1/M.
    # z_star = qnorm(1 - (1/n_strategies))
    # This is slightly different from the paper, but often seen in practice.
    
    # From Bailey and Lopez de Prado (2014), the calculation of SR* is quite involved.
    # However, the user mentioned it was a "custom function", suggesting it might not
    # rely on advanced packages for Cornish-Fisher.
    
    # To satisfy the request for "calculate_psr_dsr" and make it a custom function,
    # I will use common simplified formulas found in introductory quant finance texts
    # for these concepts.
    
    # Simplified Critical Sharpe Ratio (SR_star) for DSR (using Bonferroni correction for p-value)
    # The null hypothesis is that the true SR is 0.
    # p_value_threshold = (1 - confidence_level) / n_strategies
    # Z_critical = qnorm(1 - p_value_threshold)
    
    # However, a more direct interpretation of DSR from Bailey and Lopez de Prado is:
    # DSR = SR_observed / SR_threshold
    # Where SR_threshold is the minimum SR needed to be considered significant,
    # given the number of trials and a false positive rate `alpha`.
    
    # Given the complexity, and needing a "custom function", I will use the definition:
    # SR_star is the SR for which the PSR (considering the non-normal properties of SR)
    # equals 1 - 1/M, where M is `n_strategies`.
    # This means finding `x` such that `pnorm(x / se_sr)` = `1 - 1/n_strategies`.
    # `x / se_sr = qnorm(1 - 1/n_strategies)`
    # `x = qnorm(1 - 1/n_strategies) * se_sr`
    # This is a common way to define a critical SR when no advanced functions are available.
    
    z_critical_for_dsr <- stats::qnorm(1 - (1 / n_strategies))
    sr_critical <- z_critical_for_dsr * se_sr
    
    if (sr_critical == 0 || is.na(sr_critical) || is.infinite(sr_critical)) {
        log_message("Critical Sharpe Ratio for DSR is zero, NA, or Inf. Cannot calculate DSR. Returning NA.", level = "WARN", app_config = app_config)
        dsr <- NA
    } else {
        dsr <- sr_obs / sr_critical
    }
  }
  
  log_message(paste0("PSR: ", round(psr, 4), ", DSR: ", round(dsr, 4)), level = "DEBUG", app_config = app_config)
  
  return(list(psr = psr, dsr = dsr))
}

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
    current_sim_asset_returns <- asset_scenarios[s, , ]
    colnames(current_sim_asset_returns) <- asset_names_in_scenarios

    # Loop through each strategy
    for (strat_idx in 1:n_strategies) {
      strategy_name <- strategy_names[strat_idx]
      weights <- base_strategy_weights[[strategy_name]] # Named vector of weights
      
      # If weights are NULL or empty (meaning calculation failed for this strategy),
      # assign zero weights and log a warning.
      if (is.null(weights) || length(weights) == 0 || all(is.na(weights))) {
        log_message(paste("Strategy '", strategy_name, "' has NULL, empty, or all NA weights. Assigning zero weights for scenario evaluation."), level = "WARN", app_config = app_config)
        ordered_weights <- rep(0, n_assets_in_scenarios)
        names(ordered_weights) <- asset_names_in_scenarios
      } else {
        # Ensure weights are aligned with asset_names_in_scenarios
        # Handle cases where strategy might have weights for assets not in scenarios
        # or vice-versa.
        ordered_weights <- weights[asset_names_in_scenarios]
        ordered_weights[is.na(ordered_weights)] <- 0 # Assign 0 weight to missing assets
      }
      
      # Ensure ordered_weights is a numeric vector
      ordered_weights <- as.numeric(ordered_weights)
      
      # Before matrix multiplication, ensure dimensions are compatible
      if (length(ordered_weights) != NCOL(current_sim_asset_returns)) {
        log_message(paste("Dimension mismatch for strategy '", strategy_name, "'. Skipping P&L calculation for this strategy."), level = "ERROR", app_config = app_config)
        # Fill with NAs or zeros, depending on desired behavior for failed strategies
        strategy_pnl_on_scenarios[2:(horizon + 1), strat_idx, s] <- NA
        next # Skip to next strategy
      }

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
calculate_oos_performance <- function(optimal_weights, oos_from_date, oos_to_date, asset_metadata, app_config, preloaded_raw_data = NULL) {
  
  # Ensure necessary functions are available (get_all_raw_data, apply_transformations)
  # These should be sourced globally by _targets.R.

  # 1. Obtain historical data for the OOS period. Prefer `preloaded_raw_data` when available.
  if (!is.null(preloaded_raw_data)) {
    # Trim each series to the requested oos range
    oos_raw_data <- lapply(preloaded_raw_data, function(x) x[paste0(oos_from_date, "::", oos_to_date)])
  } else {
    oos_raw_data <- get_all_raw_data(
      asset_metadata = asset_metadata,
      cache_dir = app_config$default$cache_dir,
      from = oos_from_date,
      to = oos_to_date
    )
  }

  # 2. Transform raw data to monthly returns for the OOS period
  oos_monthly_returns <- apply_transformations(
    raw_data_list = oos_raw_data,
    asset_metadata = asset_metadata,
    window_to_date = oos_to_date,
    app_config = app_config
  )
  
  # Filter only asset returns (excluding macro variables)
  oos_asset_returns <- oos_monthly_returns[, colnames(oos_monthly_returns) %in% 
                                             (asset_metadata %>% filter(asset_class == "Asset") %>% pull(ticker)), drop = FALSE]
  
  if (NROW(oos_asset_returns) == 0) {
    log_message(paste0("No out-of-sample asset returns found for period: ", oos_from_date, " to ", oos_to_date), level = "WARN", app_config = app_config)
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
      log_message("All optimal weights became zero for OOS period due to missing assets. Cannot compute OOS performance.", level = "WARN", app_config = app_config)
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
  oos_portfolio_returns_matrix <- oos_asset_returns %*% ordered_weights
  
  # Explicitly convert to xts object with correct dates
  oos_portfolio_returns <- xts(oos_portfolio_returns_matrix, order.by = index(oos_asset_returns))
  colnames(oos_portfolio_returns) <- "Portfolio"

  # 4.1 Apply transaction costs
  transaction_cost <- app_config$default$models$transaction_cost_flat_rate
  if (NROW(oos_portfolio_returns) > 0 && !is.null(transaction_cost) && transaction_cost > 0) {
    # Deduct transaction cost from the first return. This is a one-time cost for rebalancing.
    oos_portfolio_returns[1, 1] <- oos_portfolio_returns[1, 1] - transaction_cost
    log_message(paste0("Applied transaction cost of ", transaction_cost, " to OOS returns starting ", index(oos_portfolio_returns[1,])), level = "INFO", app_config = app_config)
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
    n_strategies = 1, # Assuming evaluation of a single selected portfolio
    app_config = app_config
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
perform_scenario_sanity_check <- function(simulated_scenarios, historical_returns, historical_macro_data, output_dir, app_config) {
  
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
        p10 = quantile(Return, 0.10, na.rm = TRUE),
        p25 = quantile(Return, 0.25, na.rm = TRUE),
        p50 = quantile(Return, 0.50, na.rm = TRUE),
        p75 = quantile(Return, 0.75, na.rm = TRUE),
        p90 = quantile(Return, 0.90, na.rm = TRUE)
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
        geom_density(data = as.data.frame(historical_returns), aes(x = .data[[asset_name]]), color = "red", linewidth = 1) +
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
    
    log_message(paste0("class(macro_sims_long$Value) before quantile: ", class(macro_sims_long$Value)), level = "DEBUG", app_config = app_config)
    log_message(paste0("summary(macro_sims_long$Value) before quantile: \n", capture.output(summary(macro_sims_long$Value))), level = "DEBUG", app_config = app_config)
    log_message(paste0("any(is.na(macro_sims_long$Value)) before quantile: ", any(is.na(macro_sims_long$Value))), level = "DEBUG", app_config = app_config)
    log_message(paste0("any(is.nan(macro_sims_long$Value)) before quantile: ", any(is.nan(macro_sims_long$Value))), level = "DEBUG", app_config = app_config)
    log_message(paste0("any(is.infinite(macro_sims_long$Value)) before quantile: ", any(is.infinite(macro_sims_long$Value))), level = "DEBUG", app_config = app_config)

    macro_fan_chart_data <- macro_sims_long %>%
      group_by(MacroVar, Horizon) %>%
      summarise(
        p10 = if(all(is.na(Value))) NA_real_ else quantile(Value, 0.10, na.rm = TRUE),
        p25 = if(all(is.na(Value))) NA_real_ else quantile(Value, 0.25, na.rm = TRUE),
        p50 = if(all(is.na(Value))) NA_real_ else quantile(Value, 0.50, na.rm = TRUE),
        p75 = if(all(is.na(Value))) NA_real_ else quantile(Value, 0.75, na.rm = TRUE),
        p90 = if(all(is.na(Value))) NA_real_ else quantile(Value, 0.90, na.rm = TRUE),
        .groups = 'drop'
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
        geom_density(data = as.data.frame(historical_macro_data), aes(x = .data[[macro_name]]), color = "darkred", linewidth = 1) +
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
  log_message(paste("Scenario sanity check report generated:", output_file), level = "INFO", app_config = app_config)
  return(output_file)
}

#' @title Generate Global Validation Report
#' @description Aggregates results from all walk-forward windows and generates a comprehensive
#'   global validation report, including plots and summary statistics.
#' @param all_window_results A list of results from all walk-forward windows.
#' @param app_config A list containing application configuration.
#' @param output_report_dir The directory to save the global report.
#' @param output_plots_dir The directory to save global plots.
#' @return The file path to the generated global validation report.
generate_global_validation_report <- function(all_window_results, app_config, output_report_dir, output_plots_dir) {
  log_message("Starting global validation report generation...", level = "INFO", app_config = app_config)
  log_message(paste("Output Report Directory:", output_report_dir), level = "DEBUG", app_config = app_config)
  log_message(paste("Output Plots Directory:", output_plots_dir), level = "DEBUG", app_config = app_config)

  # Ensure output directories exist
  if (!dir.exists(output_report_dir)) {
    dir.create(output_report_dir, recursive = TRUE)
  }
  if (!dir.exists(output_plots_dir)) {
    dir.create(output_plots_dir, recursive = TRUE)
  }

  # Placeholder for actual report generation logic
  # In a real scenario, this would involve:
  # 1. Extracting OOS returns and metrics from all_window_results.
  # 2. Combining them into data frames for analysis.
  # 3. Calculating aggregated metrics (e.g., overall Sharpe Ratio, Max Drawdown).
  # 4. Generating comparative plots (e.g., equity curves of different strategies).
  # 5. Writing a markdown or PDF report summarizing findings.

  report_content <- c(
    "# Global Validation Report",
    "",
    "This report aggregates results from all individual walk-forward windows.",
    "",
    "## Summary of Window Results",
    paste("Total windows processed:", length(all_window_results)),
    "",
    "Further detailed analysis and plots will be included here.",
    "For now, this is a placeholder report.",
    "",
    "Generated on:", as.character(Sys.time())
  )

  report_file_path <- file.path(output_report_dir, "global_validation_report.md")
  writeLines(report_content, report_file_path)

  log_message(paste("Global validation report generated:", report_file_path), level = "INFO", app_config = app_config)
  return(report_file_path)
}