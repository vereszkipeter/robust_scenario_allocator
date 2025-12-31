library(rugarch)
library(rmgarch)
library(bsvars)
library(xts)
library(MASS)
library(dplyr)
library(abind) # For combining arrays

#' @title Generate Full Scenarios from Fitted RS-BVAR and DCC-GARCH Models
#' @description Generates future scenarios for both macro variables and asset returns
#' by combining simulations from the fitted RS-BVAR and DCC-GARCH models.
#' @param fitted_generative_model A list containing:
#'   - `rsbvar_model`: A fitted bvars object for macro variables.
#'   - `dcc_garch_model`: A fitted DCCfit object for asset returns.
#' @param n_sim The number of simulations (scenarios) to generate.
#' @param horizon The time horizon in months for each scenario.
#' @param asset_metadata A tibble containing asset information, used to identify macro/asset tickers.
#' @param last_historical_levels An xts object containing the last historical *levels/prices*
#'   for all series (macro and asset) to be used as starting points for inverse transformations.
#' @return A list of 3D arrays with named elements `macro_scenarios` and `asset_scenarios`,
#'   each of dimensions `n_sim x horizon x n_variables`, with values in original levels/returns.
generate_full_scenarios <- function(fitted_generative_model, n_sim, horizon, asset_metadata, last_historical_levels, app_config, seed = NULL) {
  
  # Simplified and robust seed handling.
  final_seed <- get_valid_seed(app_config)
  # Global set.seed() removed from here. It will be applied just before dccsim.
  log_message(paste("Obtained validated seed for scenario generation:", final_seed), level = "DEBUG", app_config = app_config)
  
  rsbvar_fitted_model <- fitted_generative_model$rsbvar_fitted_model
  rsbvar_spec <- fitted_generative_model$rsbvar_spec
  macro_original_colnames <- fitted_generative_model$macro_original_colnames
  dcc_garch_model <- fitted_generative_model$dcc_garch_model

  log_message(paste0("class(rsbvar_fitted_model): ", paste(class(rsbvar_fitted_model), collapse = ", ")), level = "DEBUG", app_config = app_config)
  log_message(paste0("class(dcc_garch_model): ", paste(class(dcc_garch_model), collapse = ", ")), level = "DEBUG", app_config = app_config)

  # --- 1. Simulate Macro Variables if RS-BVAR model is provided ---
  macro_scenarios_array <- NULL 
  macro_var_names <- NULL
  if (!is.null(rsbvar_fitted_model) && inherits(rsbvar_fitted_model, "PosteriorBSVARMSH")) {
    log_message("Simulating macro scenarios from RS-BVAR model...", level = "INFO", app_config = app_config)
    macro_var_names <- macro_original_colnames
    n_macro_vars <- length(macro_var_names)
    
    rsbvar_forecasts <- bsvars::forecast(rsbvar_fitted_model, horizon = horizon)
    
    macro_forecasts_raw <- aperm(rsbvar_forecasts$forecasts, c(3, 2, 1))
    
    if (any(!is.finite(macro_forecasts_raw))) {
      warning("RS-BVAR forecast draws contain non-finite values. Skipping macro scenario generation.")
      macro_scenarios_array <- array(NA_real_, dim = c(n_sim, horizon, n_macro_vars), dimnames = list(NULL, NULL, macro_var_names))
    } else {
      n_draws_available <- dim(macro_forecasts_raw)[1]
      draw_indices <- sample(1:n_draws_available, n_sim, replace = n_sim > n_draws_available)
      if (n_sim > n_draws_available) {
        warning(paste0("Requested simulations (", n_sim, ") exceeds available RS-BVAR draws (", n_draws_available, "). Sampling with replacement."))
      }
      
      macro_sims_transformed <- macro_forecasts_raw[draw_indices, , , drop = FALSE]
      dimnames(macro_sims_transformed) <- list(NULL, NULL, macro_var_names)

      macro_scenarios_reconstructed <- array(NA, dim = dim(macro_sims_transformed), dimnames = dimnames(macro_sims_transformed))
      for (j in 1:n_macro_vars) {
        ticker <- macro_var_names[j]
        transform_name <- asset_metadata$transform_name[asset_metadata$ticker == ticker][1]
        for (s in 1:n_sim) {
          last_level <- as.numeric(tail(last_historical_levels[, ticker], 1))
          if (length(last_level) == 0) last_level <- NA
          reconstructed_series <- reconstruct_macro_series(macro_sims_transformed[s, , j], transform_name, last_level)
          macro_scenarios_reconstructed[s, , j] <- reconstructed_series
        }
      }
      macro_scenarios_array <- macro_scenarios_reconstructed
    }
  } else {
    log_message("RS-BVAR model not provided or invalid. Skipping macro scenario simulation.", level = "INFO", app_config = app_config)
  }

  # --- 2. Simulate Asset Returns from DCC-GARCH ---
  log_message("Simulating asset scenarios from DCC-GARCH model...", level = "INFO", app_config = app_config)
  if (is.null(dcc_garch_model) || !inherits(dcc_garch_model, "DCCfit")) {
    stop("A valid DCC-GARCH model ('DCCfit' object) must be provided.")
  }
    
  n_assets <- ncol(dcc_garch_model@model$modeldata$data)
  asset_var_names <- colnames(dcc_garch_model@model$modeldata$data)

  # Set seed directly before dccsim call for robustness
  set.seed(final_seed)
  log_message(paste("Applying seed", final_seed, "directly before dccsim call."), level = "DEBUG", app_config = app_config)

  sim_all <- tryCatch({
    rmgarch::dccsim(
      fitORspec = dcc_garch_model,
      n.sim = horizon,
      m.sim = n_sim,
      startMethod = "unconditional"
    )
  }, error = function(e) {
    stop(paste("generate_full_scenarios failed during dccsim execution:", e$message))
  })

  if (is.null(sim_all) || is.null(sim_all@msim$simX)) {
    stop("dccsim returned a NULL or incomplete object.")
  }

  # --- FIX: Robust reshaping of simulation output ---
  sim_list <- sim_all@msim$simX
  
  if (is.list(sim_list)) {
      # It's a list of matrices, bind them into a 3D array.
      # Expected dimensions of each matrix: horizon x n_assets
      # Resulting array dimension: horizon x n_assets x n_sim
      temp_array <- abind::abind(sim_list, along = 3)
  } else if (is.array(sim_list) && length(dim(sim_list)) == 3) {
      # It's already a 3D array, just assign it.
      # dccsim can sometimes return horizon x n_sim x n_assets
      temp_array <- sim_list
  } else {
      stop("Unexpected data structure for dccsim output (simX). Expected a list of matrices or a 3D array.")
  }

  # Permute the array to the standard dimension order: n_sim x horizon x n_assets
  # The original temp_array is likely [horizon, n_assets, n_sim] or a permutation thereof.
  # The target is [n_sim, horizon, n_assets]
  expected_dims <- c(n_sim, horizon, n_assets)
  current_dims <- dim(temp_array)
  
  if (setequal(current_dims, expected_dims)) {
      perm_map <- match(expected_dims, current_dims)
      asset_sims_transformed <- aperm(temp_array, perm_map)
  } else {
      stop(paste0("Dimension mismatch in dccsim output. Expected dims containing ", 
                  paste(expected_dims, collapse = ","), " but got ", 
                  paste(current_dims, collapse = ",")))
  }
  # --- End FIX ---
  
  dimnames(asset_sims_transformed) <- list(NULL, NULL, asset_var_names)
  
  # Scale simulated asset returns back down (they were scaled up by 100 for fitting)
  log_message("Scaling simulated asset returns back down by 100.", level = "DEBUG", app_config = app_config)
  asset_scenarios_array <- asset_sims_transformed / 100

  log_message(paste0("asset_scenarios_array final dim: ", paste(dim(asset_scenarios_array), collapse = ", ")), level = "DEBUG", app_config = app_config)
  log_message("Full scenarios generated successfully.", level = "INFO", app_config = app_config)
  
  return(list(macro_scenarios = macro_scenarios_array, asset_scenarios = asset_scenarios_array))
}



#' @title Extract Short Rate from RS-BVAR Macro Scenarios
#' @description Extracts the short rate time series from the simulated macro scenarios.
#' Assumes the short rate is identified by a specific name in the macro variables.
#' @param macro_scenarios A 3D array (n_sim, horizon, n_macro_vars) of simulated macro variables.
#' @param app_config A list containing application configuration parameters,
#'                   which may include the name of the short rate variable.
#' @return A 3D array (n_sim, horizon, 1) containing only the short rate scenarios.
extract_short_rate_from_rsbvar_scenarios <- function(macro_scenarios, app_config) {
  if (is.null(macro_scenarios)) {
    log_message("macro_scenarios is NULL. Returning NULL short rate scenarios.", level = "INFO", app_config = app_config)
    return(NULL) # Return NULL if macro_scenarios is NULL
  }

  # The name of the short rate variable. From asset_metadata, "BIL" is a good proxy.
  # If a specific name for the short rate variable is in the config, use it.
  # Otherwise, look for "BIL" or a similar proxy in the column names.
  
  # For now, let's assume one of the macro variables corresponds to the short rate.
  # In R/02_models.R, "BIL" is an asset, not a direct macro variable in BVAR.
  # Macro variables are CPIAUCSL, INDPRO, VIX.
  # We need to decide which macro variable from the BVAR output acts as the short rate proxy.
  # Or, the short rate is derived from these.
  
  # Let's assume for now that there is a variable explicitly named "SHORT_RATE"
  # or that we can infer it from the BVAR output names.
  # The `app_config` should ideally specify this.
  
  short_rate_var_name <- app_config$default$models$short_rate_variable_name

  if (is.null(short_rate_var_name) || short_rate_var_name == "") {
    stop("Configuration error: 'short_rate_variable_name' must be specified in app_config$default$models.")
  }
  
  macro_var_names <- dimnames(macro_scenarios)[[3]]

  if (!(short_rate_var_name %in% macro_var_names)) {
    stop(paste("Short rate variable '", short_rate_var_name, "' not found in macro scenarios. Available macro variables are: ", paste(macro_var_names, collapse = ", "), "."))
  }

  short_rate_scenarios <- macro_scenarios[, , short_rate_var_name, drop = FALSE]
  
  return(short_rate_scenarios)
}