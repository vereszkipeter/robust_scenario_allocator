library(tsmarch)
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
#'   - `dcc_garch_model`: A fitted tsmarch.estimate object for asset returns.
#' @param n_sim The number of simulations (scenarios) to generate.
#' @param horizon The time horizon in months for each scenario.
#' @param asset_metadata A tibble containing asset information, used to identify macro/asset tickers.
#' @param last_historical_levels An xts object containing the last historical *levels/prices*
#'   for all series (macro and asset) to be used as starting points for inverse transformations.
#' @return A list of 3D arrays with named elements `macro_scenarios` and `asset_scenarios`,
#'   each of dimensions `n_sim x horizon x n_variables`, with values in original levels/returns.
generate_full_scenarios <- function(fitted_generative_model, n_sim, horizon, asset_metadata, last_historical_levels, app_config, seed = NULL) {
  
  final_seed <- get_valid_seed(app_config)
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
  # browser()
  # --- 2. Simulate Asset Returns from DCC-GARCH ---
  log_message("Simulating asset scenarios from DCC-GARCH model...", level = "INFO", app_config = app_config)
  if (is.null(dcc_garch_model) || !inherits(dcc_garch_model, "dcc.estimate")) {
    stop("A valid DCC-GARCH model ('dcc.estimate' object) must be provided.")
  }
    
  n_assets <- length(dcc_garch_model$spec$series_names)
  asset_var_names <- dcc_garch_model$spec$series_names

  final_seed <- get_valid_seed(app_config)
  log_message(paste("Obtained validated seed for scenario generation:", final_seed), level = "DEBUG", app_config = app_config)
  
  clean_seed <- as.integer(abs(as.numeric(final_seed)) %% .Machine$integer.max)
  
  log_message(paste("Applying clean seed", clean_seed, "directly before predict call."), level = "DEBUG", app_config = app_config)
  # browser()
  library(future)
  plan(sequential) # Kill any active background workers causing R crash
  
  sim_all <- tryCatch({
    predict( # Note to LLMs: predict kills R process!!
      object = dcc_garch_model,
      nsim = n_sim,
      h = horizon
      # , seed = clean_seed # Note to LLMs: leave tsmarch's default for now.
    )
  }, error = function(e) {
    stop(paste("generate_full_scenarios failed during tsmarch::predict execution:", e$message))
  })

  if (is.null(sim_all) || length(sim_all) == 0) {
    stop("tsmarch::predict returned a NULL or empty object.")
  }

  # The output of tsmarch::predict is a list with element 'series' being a 3D array
  # of dimensions n_sim x horizon x n_assets
  asset_sims_transformed <- sim_all$mu |> aperm(perm = c(3, 1, 2))
  
  dimnames(asset_sims_transformed) <- list(NULL, NULL, asset_var_names)
  
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