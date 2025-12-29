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

  # Check for minimum column and row requirements for bsvars
  if (NCOL(macro_matrix) < 2) {
    message("WARNING: macro_data has less than 2 columns (", NCOL(macro_matrix), "). Skipping RS-BVAR model fitting as bsvars requires at least 2 variables.")
    return(list(fitted_model = NULL, spec = NULL))
  }
  if (NROW(macro_matrix) < 3) { # bsvars::specify_bsvar_msh$new requires at least 3 rows
    message("WARNING: macro_data has less than 3 rows (", NROW(macro_matrix), "). Skipping RS-BVAR model fitting as bsvars requires at least 3 observations.")
    return(list(fitted_model = NULL, spec = NULL))
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

  # --- Convergence Diagnostics ---
  # Create a directory for diagnostics if it doesn't exist
  diag_dir <- "output/rsbvar_diagnostics"
  if (!dir.exists(diag_dir)) {
    dir.create(diag_dir, recursive = TRUE)
  }
  
  # Generate a unique window name from the data's start and end dates
  window_name <- "unidentified_window"
  if (requireNamespace("xts", quietly = TRUE) && xts::is.xts(macro_data) && nrow(macro_data) > 0) {
    start_date <- format(zoo::index(macro_data)[1], "%Y-%m-%d")
    end_date <- format(zoo::index(macro_data)[nrow(macro_data)], "%Y-%m-%d")
    window_name <- paste0("window_", start_date, "_to_", end_date)
  }
  
  # Save summary statistics, which include potential convergence diagnostics like R-hat
  summary_file <- file.path(diag_dir, paste0(window_name, "_summary.txt"))
  tryCatch({
    model_summary <- summary(fitted_model)
    utils::capture.output(print(model_summary), file = summary_file)
    message(paste("RS-BVAR model summary saved to", summary_file))
  }, error = function(e) {
    warning(paste("Could not save RS-BVAR model summary:", e$message))
  })
  
  # Save trace plots to a PDF for visual inspection
  trace_plot_file <- file.path(diag_dir, paste0(window_name, "_trace_plots.pdf"))
  tryCatch({
    grDevices::pdf(trace_plot_file, width = 8, height = 10)
    plot(fitted_model) # The default plot for bsvars objects is MCMC trace plots
    grDevices::dev.off()
    message(paste("RS-BVAR trace plots saved to", trace_plot_file))
  }, error = function(e) {
    warning(paste("Could not save RS-BVAR trace plots:", e$message))
  })
  # --- End Convergence Diagnostics ---
  
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
  # Remove NAs from the input data
  asset_returns <- na.omit(asset_returns)

  # Save the input data for inspection
  saveRDS(asset_returns, "dcc_input.rds")

  message("DEBUG: NROW(asset_returns) inside fit_dcc_t_garch_model: ", NROW(asset_returns))

  min_obs <- app_config$default$models$min_dcc_obs
  # Basic validation
  if (is.null(asset_returns) || NROW(asset_returns) < min_obs || NCOL(asset_returns) < 1) {
    message("Insufficient data for DCC-GARCH fitting. Returning fallback with empirical moments.")
    emp_mean <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) colMeans(asset_returns, na.rm = TRUE) else numeric(0)
    emp_cov <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) cov(as.matrix(asset_returns), use = "pairwise.complete.obs") else matrix(NA_real_)
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }

  num_assets <- ncol(asset_returns)

  # Compute empirical moments for fallback
  emp_mean <- colMeans(asset_returns, na.rm = TRUE)
  emp_cov <- cov(as.matrix(asset_returns), use = "pairwise.complete.obs")

  # --- Stage 1: Fit univariate GARCH models individually for diagnostics ---
  for (i in 1:num_assets) {
    asset_name <- colnames(asset_returns)[i]
    message("--- Fitting univariate GARCH for asset: ", asset_name, " ---")

    uspec <- rugarch::ugarchspec(
      mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
      variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
      distribution.model = "std"
    )

    # Fit the model for the single, NA-free series
    fit <- tryCatch({
      rugarch::ugarchfit(spec = uspec, data = asset_returns[, i, drop = TRUE], solver = app_config$default$models$dcc_solver)
    }, error = function(e) {
      message("ERROR: Univariate GARCH fit failed for asset: ", asset_name, "; will use empirical fallback.")
      NULL
    })

    if (is.null(fit)) {
      # Return a fallback object instead of stopping the whole pipeline
      return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
    }
  }
  message("--- All univariate GARCH models fitted successfully. ---")

  # --- Stage 2: If all univariate fits succeeded, proceed with DCC ---

  uspec_list <- lapply(1:num_assets, function(i) {
    rugarch::ugarchspec(
      mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
      variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
      distribution.model = "std"
    )
  })
  multi_uspec <- rugarch::multispec(uspec_list)
  message("DEBUG: multi_uspec created.")
  saveRDS(multi_uspec, "multi_uspec.rds") # Save multi_uspec for inspection

  message("DEBUG: Calling dccspec.")
  dcc_spec <- rmgarch::dccspec(
    uspec = multi_uspec,
    dccOrder = c(1, 1),
    distribution = "mvt" # Changed from "mvt" to "mvnorm"
  )
  saveRDS(dcc_spec, "dcc_spec.rds") # Save dcc_spec for inspection

  message("DEBUG: Calling dccfit with data dimensions: ", paste(dim(asset_returns), collapse = ", "))
  message("DEBUG: Summary of asset_returns before dccfit:")
  print(summary(asset_returns))

  # Additional data quality checks before dccfit
  if (anyNA(asset_returns)) {
    message("ERROR: NA values found in asset_returns just before dccfit. This should not happen after na.omit. Returning fallback.")
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }
  if (any(is.infinite(asset_returns))) {
    message("ERROR: Infinite values found in asset_returns just before dccfit. Returning fallback.")
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }
  if (NROW(asset_returns) < min_obs || NCOL(asset_returns) < 1) {
    message("ERROR: Insufficient data dimensions (", NROW(asset_returns), " rows, ", NCOL(asset_returns), " cols) in asset_returns just before dccfit. Returning fallback.")
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }

  message("DEBUG: Correlation matrix of asset_returns before dccfit:")
  tryCatch({
    print(cor(asset_returns, use = "pairwise.complete.obs"))
  }, error = function(e) {
    message("WARNING: Could not compute correlation matrix: ", e$message)
  })

  # Save the asset_returns data right before the dccfit call
  saveRDS(asset_returns, "dcc_input_for_fit.rds")

  dcc_fit_model <- tryCatch({
    rmgarch::dccfit(
      dcc_spec,
      data = asset_returns, # Pass the cleaned data
      solver = app_config$default$models$dcc_solver,
      solver.control = app_config$default$models$dcc_solver_control
    )
  }, error = function(e) {
    message("DCC fit failed with error: ", e$message)
    # Return fallback empirical moments to allow scenario generation to continue
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  })

  return(dcc_fit_model)
}