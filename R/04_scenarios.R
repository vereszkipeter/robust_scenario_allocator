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
generate_full_scenarios <- function(fitted_generative_model, n_sim, horizon, asset_metadata, last_historical_levels, app_config, seed = NULL) {
  
  # Ensure reproducibility
  # set.seed(123) # Removed to avoid conflicts with targets' seed management and internal package seeding
  if (is.null(seed)) {
    seed <- get_valid_seed(app_config)
  }
  set.seed(seed)
  
  rsbvar_fitted_model <- fitted_generative_model$rsbvar_fitted_model
  rsbvar_spec <- fitted_generative_model$rsbvar_spec
  macro_original_colnames <- fitted_generative_model$macro_original_colnames
  dcc_garch_model <- fitted_generative_model$dcc_garch_model
# browser()
  # --- DEBUGGING START ---
  log_message(paste0("class(rsbvar_fitted_model): ", paste(class(rsbvar_fitted_model), collapse = ", ")), level = "DEBUG", app_config = app_config)
  log_message(paste0("class(rsbvar_spec): ", paste(class(rsbvar_spec), collapse = ", ")), level = "DEBUG", app_config = app_config)
  log_message(paste0("macro_original_colnames: ", paste(macro_original_colnames, collapse = ", ")), level = "DEBUG", app_config = app_config)
  log_message(paste0("class(dcc_garch_model): ", paste(class(dcc_garch_model), collapse = ", ")), level = "DEBUG", app_config = app_config)
  # --- DEBUGGING END ---

  # Basic input validation
  if (is.null(asset_metadata) || NROW(asset_metadata) == 0) stop("asset_metadata must be provided and non-empty")
  if (is.null(last_historical_levels) || NROW(last_historical_levels) == 0) warning("last_historical_levels is missing or empty; reconstruction may be inaccurate")

  # --- 1. Simulate Macro Variables if RS-BVAR model is provided ---
  macro_scenarios_array <- NULL # Initialize macro scenarios as NULL
  macro_var_names <- NULL
  if (!is.null(rsbvar_fitted_model) && inherits(rsbvar_fitted_model, "PosteriorBSVARMSH")) { # Check for specific class
    log_message("Simulating macro scenarios from RS-BVAR model...", level = "DEBUG", app_config = app_config)
    macro_var_names <- macro_original_colnames # Use the original column names
    log_message(paste0("macro_var_names (from original colnames): ", paste(macro_var_names, collapse = ", ")), level = "DEBUG", app_config = app_config)
    log_message(paste0("Class of rsbvar_spec$get_data_matrices()$Y: ", paste(class(rsbvar_spec$get_data_matrices()$Y), collapse = ", ")), level = "DEBUG", app_config = app_config)
    n_macro_vars <- length(macro_var_names)
    
    # Set a validated seed for macro simulation; fallback gracefully if invalid
    # set.seed(123) # Removed: targets handles global seeding


    
    # Use the forecast method for the bsvars object
    rsbvar_forecasts <- bsvars::forecast(
      rsbvar_fitted_model,
      horizon = horizon
    )
    
    # bsvars::forecast returns a list of class Forecasts.
    # The 'forecasts' component is an array of dimensions [N_vars, Horizon, N_draws].
    # We need to reorder this to [N_draws, Horizon, N_vars].
    macro_forecasts_raw <- aperm(rsbvar_forecasts$forecasts, c(3, 2, 1))
    
    n_draws_available <- dim(macro_forecasts_raw)[1]

    # Basic sanity checks for MCMC forecast draws: ensure no extreme or non-finite values
    any_nonfinite <- any(!is.finite(macro_forecasts_raw))
    any_large <- any(abs(macro_forecasts_raw) > 1e6)
    if (any_nonfinite || any_large) {
      warning("RS-BVAR forecast draws contain non-finite or extremely large values. This may indicate MCMC divergence or model misspecification. Skipping macro scenario reconstruction for this window.")
      # Create an NA-filled transformed array so later code handles missing macro scenarios gracefully
      macro_sims_transformed <- array(NA_real_, dim = c(n_sim, horizon, n_macro_vars))
      dimnames(macro_sims_transformed)[[3]] <- macro_var_names
      macro_scenarios_array <- array(NA_real_, dim = dim(macro_sims_transformed))
      dimnames(macro_scenarios_array) <- dimnames(macro_sims_transformed)
      # Skip further macro reconstruction
      goto_skip_macro <- TRUE
    } else {
      goto_skip_macro <- FALSE
    }
    
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

    if (!exists("goto_skip_macro") || !isTRUE(goto_skip_macro)) {
      for (j in 1:n_macro_vars) {
        ticker <- macro_var_names[j]
        log_message(paste0("Processing macro ticker for inverse transform: ", ticker), level = "DEBUG", app_config = app_config)

        filter_result <- asset_metadata %>% dplyr::filter(ticker == !!ticker)
        if (NROW(filter_result) == 0) {
          stop(paste0("Macro ticker '", ticker, "' not found in asset_metadata for inverse transformation."))
        }

        transform_name <- filter_result %>% dplyr::pull(transform_name) %>% dplyr::first()
        log_message(paste0("transform_name for ", ticker, ": ", transform_name), level = "DEBUG", app_config = app_config)

        if (is.null(transform_name) || length(transform_name) == 0) {
          stop(paste0("transform_name for macro ticker '", ticker, "' is NULL or empty."))
        }

        for (s in 1:n_sim) {
          sim_series <- as.numeric(macro_sims_transformed[s, , j])
          last_level <- if (!is.null(last_historical_levels) && ticker %in% colnames(last_historical_levels)) {
            as.numeric(tail(last_historical_levels[, ticker], 1))
          } else {
            NA_real_
          }

          log_message(paste0("Ticker: ", ticker, ", last_level for reconstruction: ", last_level), level = "DEBUG", app_config = app_config)
          log_message(paste0("Raw simulated series (sim_series) for ", ticker, ":"), level = "DEBUG", app_config = app_config)
          print(summary(sim_series))
          log_message(paste0("any(is.na(sim_series)) for ", ticker, ": ", any(is.na(sim_series))), level = "DEBUG", app_config = app_config)
          log_message(paste0("any(is.nan(sim_series)) for ", ticker, ": ", any(is.nan(sim_series))), level = "DEBUG", app_config = app_config)
          log_message(paste0("any(is.infinite(sim_series)) for ", ticker, ": ", any(is.infinite(sim_series))), level = "DEBUG", app_config = app_config)

          reconstructed_series <- reconstruct_macro_series(sim_series = sim_series, original_transform = transform_name, last_level = last_level)
          log_message(paste0("Ticker: ", ticker, ", reconstructed_series sum(is.na): ", sum(is.na(reconstructed_series))), level = "DEBUG", app_config = app_config)
          macro_scenarios_reconstructed[s, , j] <- reconstructed_series
        }
      }
      macro_scenarios_array <- macro_scenarios_reconstructed
    } else {
      log_message("Macro scenarios skipped due to RS-BVAR forecast diagnostic failure.", level = "INFO", app_config = app_config)
    }

  } else {
    log_message("RS-BVAR model not provided or not of expected type. Skipping macro scenario simulation.", level = "INFO", app_config = app_config)
    # macro_scenarios_array remains NULL
  } # Closing brace for the else block





  # --- 2. Simulate Asset Returns from DCC-GARCH (This part is mandatory) ---
  log_message("Simulating asset scenarios from DCC-GARCH model...", level = "INFO", app_config = app_config)
  # Set a validated seed for asset simulation; fallback gracefully if invalid

  # Add debug prints and checks for dcc_garch_model components
  log_message(paste0("Class of dcc_garch_model: ", paste(class(dcc_garch_model), collapse = ", ")), level = "DEBUG", app_config = app_config)
  if (is.null(dcc_garch_model)) {
    log_message("DCC-GARCH model is NULL. Will attempt fallback simulation if empirical moments provided.", level = "INFO", app_config = app_config)
    # If NULL, attempt to use fallback info if present
  }
  
  # If we have a fitted DCC object, inspect its internals; if we have a fallback list, skip these checks
  if (!is.null(dcc_garch_model) && inherits(dcc_garch_model, "DCCfit")) {
    log_message(paste0("Class of dcc_garch_model@model: ", paste(class(dcc_garch_model@model), collapse = ", ")), level = "DEBUG", app_config = app_config)
    if (is.null(dcc_garch_model@model)) {
      stop("dcc_garch_model@model is NULL before extracting asset names.")
    }

    log_message(paste0("Class of dcc_garch_model@model$modeldata: ", paste(class(dcc_garch_model@model$modeldata), collapse = ", ")), level = "DEBUG", app_config = app_config)
    if (is.null(dcc_garch_model@model$modeldata)) {
      stop("dcc_garch_model@model$modeldata is NULL before extracting asset names.")
    }

    log_message(paste0("Class of dcc_garch_model@model$modeldata$data: ", paste(class(dcc_garch_model@model$modeldata$data), collapse = ", ")), level = "DEBUG", app_config = app_config)
    if (is.null(dcc_garch_model@model$modeldata$data)) {
      stop("dcc_garch_model@model$modeldata$data is NULL or malformed before extracting asset names.")
    }
    log_message(paste0("dim(dcc_garch_model@model$modeldata$data): ", paste(dim(dcc_garch_model@model$modeldata$data), collapse = ", ")), level = "DEBUG", app_config = app_config)
  } else if (!is.null(dcc_garch_model) && is.list(dcc_garch_model) && !is.null(dcc_garch_model$fallback) && dcc_garch_model$fallback == TRUE) {
    log_message("DCC model unavailable; using empirical fallback for asset simulation.", level = "DEBUG", app_config = app_config)
  } else if (is.null(dcc_garch_model)) {
    stop("DCC-GARCH model is NULL before extracting asset names and no fallback provided.")
  } else {
    stop("Unsupported dcc_garch_model object type: ", paste(class(dcc_garch_model), collapse = ", "))
  }

  # If dcc_garch_model is a fitted DCCfit, use rmgarch dccsim; otherwise, fall back to empirical moments
  if (!is.null(dcc_garch_model) && inherits(dcc_garch_model, "DCCfit")) {
    # Get number of assets and asset names from the fitted model's data
    n_assets <- ncol(dcc_garch_model@model$modeldata$data)
    asset_var_names <- colnames(dcc_garch_model@model$modeldata$data)

    asset_sims_transformed <- array(NA, dim = c(n_sim, horizon, n_assets))
    dimnames(asset_sims_transformed) <- list(NULL, NULL, asset_var_names) # Assign dimnames after creation

        # --- New: Prepare arguments for dccsim with optional seed ---

        dccsim_args <- list(
          fitORspec = dcc_garch_model,
          n.sim = horizon,
          m.sim = n_sim,
          startMethod = "unconditional",
          rseed = seed
        )
        
        # Perform a single simulation call for all n_sim simulations using do.call
        sim_all <- do.call(rmgarch::dccsim, dccsim_args)
    print("DEBUG: str(sim_all) ----------------------------")
    str(sim_all)
    print("DEBUG: End dccsim inspection --------------------")
    log_message(paste0("class(sim_all): ", paste(class(sim_all), collapse = ", ")), level = "DEBUG", app_config = app_config)
    if (is.null(sim_all)) stop("dccsim returned NULL")
    if (is.null(sim_all@msim)) stop("sim_all@msim is NULL")
    if (is.null(sim_all@msim$simX)) stop("sim_all@msim$simX is NULL")
    log_message(paste0("After dccsim, NROW(sim_all@msim$simX) is ", NROW(sim_all@msim$simX), " and NCOL(sim_all@msim$simX) is ", NCOL(sim_all@msim$simX)), level = "DEBUG", app_config = app_config)

    # Extract simulated series (returns)
    asset_sims_transformed <- purrr::map_depth(sim_all@msim$simX, 1, ~ t(.x)) %>%
      simplify2array() %>%
      aperm(c(3, 1, 2))

    # Ensure dimnames are correctly set
    dimnames(asset_sims_transformed)[[3]] <- asset_var_names
  } else if (!is.null(dcc_garch_model) && is.list(dcc_garch_model) && !is.null(dcc_garch_model$fallback) && dcc_garch_model$fallback == TRUE) {
    # Use empirical mean and covariance fallback supplied by fit_dcc_t_garch_model
    emp_mean <- dcc_garch_model$emp_mean
    emp_cov <- dcc_garch_model$emp_cov
    asset_var_names <- dcc_garch_model$asset_names
    n_assets <- length(asset_var_names)

    asset_sims_transformed <- array(NA, dim = c(n_sim, horizon, n_assets))
    dimnames(asset_sims_transformed) <- list(NULL, NULL, asset_var_names)

    for (s in 1:n_sim) {
      draws <- MASS::mvrnorm(n = horizon, mu = emp_mean, Sigma = emp_cov)
      # draws is horizon x n_assets
      asset_sims_transformed[s, , ] <- draws
    }
  } else {
    stop("DCC-GARCH model unavailable and no fallback provided. Cannot simulate asset scenarios.")
  }
  
  # For asset returns, the simulated values from dccsim are already simple returns.
  # No inverse transformation is needed here; the array 'asset_sims_transformed'
  # already contains the simulated simple returns.
  asset_scenarios_array <- asset_sims_transformed
  # The original loop for inverse transformation of asset returns was:
  # for (j in 1:n_assets) {
  #   ticker <- asset_var_names[j]
  #   message("DEBUG: Processing asset ticker for inverse transform: ", ticker)
  #   for (s in 1:n_sim) {
  #     sim_series <- as.numeric(asset_sims_transformed[s, , j])
  #     message("DEBUG: sim_series class: ", class(sim_series), ", head: ", paste(head(sim_series), collapse = ", "))
  #     
  #     reconstructed_series <- reconstruct_asset_returns_from_log(sim_series)
  #     message("DEBUG: reconstructed_series class: ", class(reconstructed_series), ", head: ", paste(head(reconstructed_series), collapse = ", "))
  #     
  #     asset_scenarios_reconstructed[s, , j] <- reconstructed_series
  #   }
  # }
  # asset_scenarios_array <- asset_scenarios_reconstructed
  
  log_message(paste0("asset_scenarios_array class: ", class(asset_scenarios_array), ", dim: ", paste(dim(asset_scenarios_array), collapse = ", ")), level = "DEBUG", app_config = app_config)
  log_message(paste0("macro_scenarios_array class: ", class(macro_scenarios_array), ", dim: ", paste(dim(macro_scenarios_array), collapse = ", ")), level = "DEBUG", app_config = app_config)
  
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