library(rugarch)
library(rmgarch)
library(bsvars)
library(xts) # For time series manipulation
library(MASS) # For mvrnorm for placeholder, will remove once bvars simulation is proper
library(dplyr) # For lag, slice

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
generate_full_scenarios <- function(fitted_generative_model, n_sim, horizon, asset_metadata, last_historical_levels) {
  
  # Ensure reproducibility
  set.seed(123)
  
  rsbvar_model <- fitted_generative_model$rsbvar_model
  dcc_garch_model <- fitted_generative_model$dcc_garch_model

  # Basic input validation
  if (is.null(asset_metadata) || NROW(asset_metadata) == 0) stop("asset_metadata must be provided and non-empty")
  if (is.null(last_historical_levels) || NROW(last_historical_levels) == 0) warning("last_historical_levels is missing or empty; reconstruction may be inaccurate")

  # --- 1. Simulate Macro Variables if RS-BVAR model is provided ---
  macro_scenarios_array <- NULL # Initialize macro scenarios as NULL
  macro_var_names <- NULL
  if (!is.null(rsbvar_model) && inherits(rsbvar_model, "bsvars_rsbvar_obj")) { # Check for specific class
    message("Simulating macro scenarios from RS-BVAR model...")
    macro_var_names <- colnames(rsbvar_model$Y)
    n_macro_vars <- length(macro_var_names)
    
    # Use the forecast method for the bsvars object
    rsbvar_forecasts <- bsvars::forecast(
      object = rsbvar_model,
      horizon = horizon
    )
    
    # bsvars::forecast returns a list of class Forecasts.
    # The 'forecasts' component is an array of dimensions [N_vars, Horizon, N_draws].
    # We need to reorder this to [N_draws, Horizon, N_vars].
    macro_forecasts_raw <- aperm(rsbvar_forecasts$forecasts, c(3, 2, 1))
    
    n_draws_available <- dim(macro_forecasts_raw)[1]
    
    if (n_sim > n_draws_available) {
      warning(paste0("Number of requested simulations (", n_sim, ") exceeds available RS-BVAR forecast draws (", n_draws_available, "). Sampling with replacement."))
      draw_indices <- sample(1:n_draws_available, n_sim, replace = TRUE)
    } else {
      draw_indices <- sample(1:n_draws_available, n_sim, replace = FALSE)
    }
    
    macro_sims_transformed <- array(NA, dim = c(n_sim, horizon, n_macro_vars))
    for (idx in 1:n_sim) {
      macro_sims_transformed[idx, , ] <- macro_forecasts_raw[draw_indices[idx], , ]
    }
    dimnames(macro_sims_transformed)[[3]] <- macro_var_names

    # Apply inverse transformations for macro scenarios
    macro_scenarios_reconstructed <- array(NA, dim = dim(macro_sims_transformed))
    dimnames(macro_scenarios_reconstructed) <- dimnames(macro_sims_transformed)

    for (j in 1:n_macro_vars) {
      ticker <- macro_var_names[j]
      transform_name <- asset_metadata %>% dplyr::filter(ticker == !!ticker) %>% dplyr::pull(transform_name) %>% dplyr::first()

      for (s in 1:n_sim) {
        sim_series <- as.numeric(macro_sims_transformed[s, , j])
        last_level <- if (!is.null(last_historical_levels) && ticker %in% colnames(last_historical_levels)) {
          as.numeric(tail(last_historical_levels[, ticker], 1))
        } else {
          NA_real_
        }

        reconstructed_series <- reconstruct_macro_series(sim_series = sim_series, original_transform = transform_name, last_level = last_level)
        macro_scenarios_reconstructed[s, , j] <- reconstructed_series
      }
    }
    macro_scenarios_array <- macro_scenarios_reconstructed

  } else {
    message("RS-BVAR model not provided or not of expected type. Skipping macro scenario simulation.")
    # macro_scenarios_array remains NULL




  # --- 2. Simulate Asset Returns from DCC-GARCH (This part is mandatory) ---
  message("Simulating asset scenarios from DCC-GARCH model...")

  # Get number of assets and asset names from the fitted model's data
  # The actual data is in dcc_garch_model@model$modeldata@data
  n_assets <- ncol(dcc_garch_model@model$modeldata@data)
  asset_var_names <- colnames(dcc_garch_model@model$modeldata@data)
  
  asset_sims_transformed <- array(NA, dim = c(n_sim, horizon, n_assets))
  dimnames(asset_sims_transformed) <- list(NULL, NULL, asset_var_names) # Assign dimnames after creation

  for (i in 1:n_sim) {
    # Pass the DCCfit object directly and let dccsim handle starting values
    sim <- rmgarch::dccsim(
      dcc_garch_model, # Pass the DCCfit object directly
      n.sim = horizon, # Simulate for the given horizon
      m.sim = 1, # Generate 1 path for this simulation
      startMethod = "sample" # Let dccsim automatically extract starting values from the fit
      # Removed cluster = NULL
      # prereturns, presigma, preresiduals, preQ, preZ etc. will be automatically derived from the last observation of the fitted model
    )
    
    # Extract simulated series (returns)
    # sim@simulation$series[[1]] is a matrix `horizon x n_assets`.
    asset_sims_transformed[i, , ] <- sim@simulation$series[[1]]
  }
  
  # Apply inverse transformations for asset scenarios (log returns -> simple returns)
  asset_scenarios_reconstructed <- array(NA, dim = dim(asset_sims_transformed))
  dimnames(asset_scenarios_reconstructed) <- dimnames(asset_sims_transformed)

  for (j in 1:n_assets) {
    ticker <- asset_var_names[j]
    for (s in 1:n_sim) {
      sim_series <- as.numeric(asset_sims_transformed[s, , j])
      asset_scenarios_reconstructed[s, , j] <- reconstruct_asset_returns_from_log(sim_series)
    }
  }
  asset_scenarios_array <- asset_scenarios_reconstructed
  
  message("Full scenarios generated successfully.")
  
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
    message("macro_scenarios is NULL. Returning NULL short rate scenarios.")
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