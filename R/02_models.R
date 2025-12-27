library(bsvars) # Using the correct bsvars package
library(rmgarch) # Already loaded, but good to explicitly state for context

#' @title Fit a Regime-Switching Bayesian Vector Autoregression (RS-BVAR) Model
#' @description Fits an RS-BVAR model to macro variables using the `bsvars` package.
#' This captures macro dynamics (Growth, Inflation, Rates) under P-Measure.
#' @param macro_data An xts object containing the macro time series for the BVAR.
#' @param bvar_lags An integer specifying the number of lags for the BVAR model.
#' @return A fitted `bsvars` object.
#' @param app_config A list containing application configuration, including RS-BVAR MCMC parameters.
fit_rsbvar_model <- function(macro_data, bvar_lags, app_config) {
  library(bsvars) # Explicitly load bsvars inside the function
  
  message("DEBUG: Starting fit_rsbvar_model.")
  
  # bsvars expects a matrix or data.frame
  macro_matrix <- as.matrix(macro_data)

  # Check if macro_matrix has enough observations for the given lags
  if (nrow(macro_matrix) <= bvar_lags) {
    stop("Not enough observations in macro_data for the specified bvar_lags.")
  }
  
  # Retrieve RS-BVAR parameters from app_config
  n_iter_mcmc <- app_config$default$models$rsbvar_n_iter_mcmc
  n_burnin_mcmc <- app_config$default$models$rsbvar_n_burnin_mcmc
  M <- app_config$default$models$rsbvar_M
  
  # Define N (number of variables)
  N <- ncol(macro_matrix)
  
  # Define B matrix for identification (recursive identification for now)
  # This implies a specific ordering of variables, which should be documented or configured.
  # For now, assuming a standard recursive identification (lower triangular B matrix)
  B_matrix_initial <- matrix(TRUE, N, N)
  B_matrix_initial[upper.tri(B_matrix_initial)] <- FALSE # Lower triangular identification
  
  message("DEBUG: Attempting to specify bsvars model with p=", bvar_lags, ", M=", M, ".")
  spec <- tryCatch({
    bsvars::specify_bsvar_msh$new(
      data = macro_matrix,
      p = bvar_lags,
      M = M,
      B = B_matrix_initial # Pass the B matrix for identification
    )
  }, error = function(e) {
    message("ERROR: bsvars::specify_bsvar_msh$new failed: ", e$message)
    stop(paste("bsvars::specify_bsvar_msh$new failed:", e$message))
  })
  message("DEBUG: bsvars model specification complete. Class: ", paste(class(spec), collapse = ", "))
  
  message("DEBUG: Attempting to estimate bsvars model with S=", n_iter_mcmc, ".")
  fitted_model <- tryCatch({
    bsvars::estimate(
      spec,
      S = n_iter_mcmc,
      thin = 1,
      show_progress = FALSE # Hide progress for cleaner logs
    )
  }, error = function(e) {
    message("ERROR: bsvars::estimate failed: ", e$message)
    stop(paste("bsvars::estimate failed:", e$message))
  })
  message("DEBUG: bsvars model estimation complete. Class: ", paste(class(fitted_model), collapse = ", "))
  
  message("DEBUG: Attempting to normalize posterior draws.")
  fitted_model <- tryCatch({
    BB <- NULL
    # Try different paths to get the B matrix for normalization
    if (!is.null(fitted_model$last_draw) &&
        !is.null(fitted_model$last_draw$starting_values) &&
        !is.null(fitted_model$last_draw$starting_values$B)) {
      BB <- fitted_model$last_draw$starting_values$B
    } else if (!is.null(fitted_model$get_last_draw()) &&
               !is.null(fitted_model$get_last_draw()$B)) {
      BB <- fitted_model$get_last_draw()$B
    } else if (!is.null(fitted_model$posterior) &&
               !is.null(fitted_model$posterior$B)) {
      BB <- fitted_model$posterior$B
    }
    
    if (!is.null(BB) && !all(is.na(BB))) {
      B_hat <- diag(sign(diag(BB))) %*% BB
      bsvars::normalise_posterior(fitted_model, B_hat) # Call, modifies fitted_model in place
      result_model <- fitted_model # Assign modified object
    } else {
      warning(
        "Cannot normalise posterior: Could not find valid B matrix in fitted_model. Returning un-normalised model."
      )
      result_model <- fitted_model
    }
    result_model
  }, error = function(e) {
    message("ERROR: bsvars::normalise_posterior failed: ", e$message)
    stop(paste("bsvars::normalise_posterior failed:", e$message))
  })
  message("DEBUG: Posterior normalization attempt complete. Class: ", paste(class(fitted_model), collapse = ", "))
  
  message("DEBUG: fit_rsbvar_model finished successfully.")
  return(list(fitted_model = fitted_model, spec = spec))
}

#' @title Fit a Dynamic Conditional Correlation (DCC) GARCH Model with t-Copula
#' @description Fits a DCC-GARCH model with GJR-GARCH(1,1) marginals and Student-t
#' distribution to asset returns, as specified in GEMINI.md.
#' @param asset_returns An xts object containing the asset returns time series.
#' @return A fitted `DCCfit` object.
fit_dcc_t_garch_model <- function(asset_returns, app_config) {
  # Define univariate GARCH specification for each asset
  # Based on GEMINI.md: AR(1)-GJR-GARCH(1,1) with Skewed-t distribution
  # Assuming AR(1) in mean, check if app_config provides this detail
  # For now, let's use a simple AR(1) mean model.
  # The uspec must be a list of ugarchspec objects or a multispec object.
  
  num_assets <- ncol(asset_returns)
  
  # Create a list of ugarchspec objects, one for each asset
  uspec_list <- lapply(1:num_assets, function(i) {
    rugarch::ugarchspec(
      mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
      variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
      distribution.model = "std" # Student-t distribution
    )
  })
  
  # Combine univariate specifications into a multispec object
  multi_uspec <- rugarch::multispec(uspec_list)
  
  # Define DCC specification
  dcc_spec <- rmgarch::dccspec(
    uspec = multi_uspec,
    dccOrder = c(1, 1),
    # DCC(1,1) as per standard and snippets.md example
    distribution = "mvt" # Multivariate Student-t distribution
  )
  
  # Fit the DCC-GARCH model
  # dccfit can take an optional 'cluster' argument, which might be configured in app_config
  dcc_fit_model <- rmgarch::dccfit(
    dcc_spec,
    data = asset_returns,
    solver = app_config$default$models$dcc_solver,
    # e.g., "solnp"
    solver.control = app_config$default$models$dcc_solver_control
  )
  
  return(dcc_fit_model)
}