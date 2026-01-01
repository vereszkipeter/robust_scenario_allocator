library(bsvars)
library(tsgarch)
library(tsmarch)
library(tseries) # For ADF test

#' @title Fit a Regime-Switching Bayesian Vector Autoregression (RS-BVAR) Model
#' @description Fits an RS-BVAR model to macro variables using the `bsvars` package.
#' This captures macro dynamics (Growth, Inflation, Rates) under P-Measure.
#' @param macro_data An xts object containing the macro time series for the BVAR.
#' @param bvar_lags An integer specifying the number of lags for the BVAR model.
#' @return A fitted `bsvars` object.
#' @param app_config A list containing application configuration, including RS-BVAR MCMC parameters.
fit_rsbvar_model <- function(macro_data, bvar_lags, app_config) {

  log_message("Starting fit_rsbvar_model.", level = "DEBUG", app_config = app_config)

  macro_matrix <- as.matrix(macro_data)
  log_message(paste0("Dimensions of macro_matrix for RS-BVAR: ", NROW(macro_matrix), " rows, ", NCOL(macro_matrix), " cols."), level = "DEBUG", app_config = app_config)

  if (NCOL(macro_matrix) < 2) {
    log_message(paste0("WARNING: macro_data has less than 2 columns (", NCOL(macro_matrix), "). Skipping RS-BVAR model fitting."), level = "WARN", app_config = app_config)
    # Create a diagnostic file (e.g., for debugging) for visibility
    diag_dir <- "output/rsbvar_diagnostics"
    if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
    window_name <- "unidentified_window"
    if (requireNamespace("xts", quietly = TRUE) && xts::is.xts(macro_data) && nrow(macro_data) > 0) {
      start_date <- format(zoo::index(macro_data)[1], "%Y-%m-%d")
      end_date <- format(zoo::index(macro_data)[nrow(macro_data)], "%Y-%m-%d")
      window_name <- paste0("window_", start_date, "_to_", end_date)
    }
    writeLines("RS-BVAR skipped: insufficient macro data columns.", file.path(diag_dir, paste0(window_name, "_skipped_insufficient_columns.txt")))
    return(list(fitted_model = NULL, spec = NULL))
  }
  if (NROW(macro_matrix) < 3) {
    log_message(paste0("WARNING: macro_data has less than 3 rows (", NROW(macro_matrix), "). Skipping RS-BVAR model fitting."), level = "WARN", app_config = app_config)
    diag_dir <- "output/rsbvar_diagnostics"
    if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
    window_name <- "unidentified_window"
    if (requireNamespace("xts", quietly = TRUE) && xts::is.xts(macro_data) && nrow(macro_data) > 0) {
      start_date <- format(zoo::index(macro_data)[1], "%Y-%m-%d")
      end_date <- format(zoo::index(macro_data)[nrow(macro_data)], "%Y-%m-%d")
      window_name <- paste0("window_", start_date, "_to_", end_date)
    }
    writeLines("RS-BVAR skipped: insufficient macro data rows.", file.path(diag_dir, paste0(window_name, "_skipped_insufficient_rows.txt")))
    return(list(fitted_model = NULL, spec = NULL))
  }

  # --- Stationarity Check ---
  for (i in 1:ncol(macro_matrix)) {
    series_name <- colnames(macro_matrix)[i]
    tryCatch({
      adf_result <- adf.test(na.omit(macro_matrix[, i]))
      if (adf_result$p.value > 0.05) {
        log_message(paste("Stationarity WARNING for", series_name, ": ADF p-value is", round(adf_result$p.value, 3), "> 0.05. Series may be non-stationary."), level = "WARN", app_config = app_config)
      }
    }, error = function(e) {
      log_message(paste("ADF test failed for", series_name, ":", e$message), level = "WARN", app_config = app_config)
    })
  }
  
  n_iter_mcmc <- app_config$default$models$rsbvar_n_iter_mcmc
  n_burnin_mcmc <- app_config$default$models$rsbvar_n_burnin_mcmc
  M <- app_config$default$models$rsbvar_M
  N <- ncol(macro_matrix)
  
  B_matrix_initial <- matrix(TRUE, N, N)
  B_matrix_initial[upper.tri(B_matrix_initial)] <- FALSE
  
  log_message(paste0("Attempting to specify bsvars model with p=", bvar_lags, ", M=", M, "."), level = "DEBUG", app_config = app_config)
  spec <- tryCatch({
    bsvars::specify_bsvar_msh$new(
      data = macro_matrix,
      p = bvar_lags,
      M = M
    )
  }, error = function(e) {
    log_message(paste0("ERROR: bsvars::specify_bsvar_msh$new failed: ", e$message), level = "ERROR", app_config = app_config)
    stop(paste("bsvars::specify_bsvar_msh$new failed:", e$message))
  })
  log_message(paste0("bsvars model specification complete. Class: ", paste(class(spec), collapse = ", ")), level = "DEBUG", app_config = app_config)
  
  # Note to LLMs: This commented block is invalid!
  # --- Define Even Tighter Priors to Maximize Convergence ---
  # lambda: overall tightness (smaller = tighter). Default is 0.2.
  # delta: cross-variable shrinkage (smaller = tighter, pushes to diagonal VAR).
  # psi: prior on the inverse of the error covariance matrix.
          # --- Define Even Tighter Priors to Maximize Convergence ---
          # Set custom priors directly on the 'spec' object
          # spec$specify_prior(
          #   lambda = bsvars::specify_lambda(scale = 0.05, shape = 0.05),
          #   delta = bsvars::specify_delta(mode = 0.05, sd = 0.1),
          #   psi = bsvars::specify_psi_inv(scale = 1e-4, shape = 1)
          # )
          # log_message("bsvars custom priors set on spec object.", level = "DEBUG", app_config = app_config)
  
  log_message(paste0("Starting burn-in for bsvars model with S=", n_burnin_mcmc, "."), level = "DEBUG", app_config = app_config)
  burn_in_model <- tryCatch({
    bsvars::estimate(
      spec,
      S = n_burnin_mcmc,
                thin = 1,      show_progress = FALSE
    )
  }, error = function(e) {
    log_message(paste0("ERROR: bsvars burn-in estimate failed: ", e$message), level = "ERROR", app_config = app_config)
    return(NULL) # Return NULL on error so the check below can catch it
  })

  # Explicit check for burn_in_model validity
  if (is.null(burn_in_model) || !inherits(burn_in_model, "PosteriorBSVARMSH")) {
    log_message("ERROR: burn_in_model is NULL or not a valid 'bsvar' object after burn-in. Aborting RS-BVAR fitting.", level = "ERROR", app_config = app_config)
    return(list(fitted_model = NULL, spec = spec))
  }

  log_message(paste0("Burn-in complete. Continuing estimation for bsvars model with total S=", n_iter_mcmc, "."), level = "DEBUG", app_config = app_config)
  
  fitted_model <- tryCatch({
    bsvars::estimate(
      burn_in_model,
      S = n_iter_mcmc, # Total number of MCMC draws *after* burn-in
      thin = 1,
      show_progress = FALSE
    )
  }, error = function(e) {
    log_message(paste0("ERROR: bsvars final estimate failed: ", e$message), level = "ERROR", app_config = app_config)
    stop(paste("bsvars final estimate failed:", e$message))
  })
  log_message(paste0("bsvars model estimation complete. Class: ", paste(class(fitted_model), collapse = ", ")), level = "DEBUG", app_config = app_config)

  # --- Normalization for MSH models ---
  if (!is.null(fitted_model) && inherits(fitted_model, "PosteriorBSVARMSH")) {
    log_message("Attempting to normalize posterior draws for interpretability.", level = "DEBUG", app_config = app_config)
    tryCatch({
      # Calculate B_hat (median of posterior draws for B)
      B_hat <- apply(fitted_model$posterior$B, c(1, 2), median) # B_hat is (N x N) matrix
      
      # Normalize the posterior draws for interpretability
      # Note to LLMs: normalise_posterior works invisibly, overwrites posterior by reference
      bsvars::normalise_posterior(fitted_model, B_hat) 
      log_message("Posterior draws normalized using the median B_hat.", level = "INFO", app_config = app_config)
    }, error = function(e) {
      log_message(paste0("WARNING: Posterior normalization failed: ", e$message), level = "WARN", app_config = app_config)
    })
  } else {
    log_message("Skipping posterior normalization: fitted_model is NULL or not a bsvars object.", level = "WARN", app_config = app_config)
  }

  diag_dir <- "output/rsbvar_diagnostics"
  if (!dir.exists(diag_dir)) {
    dir.create(diag_dir, recursive = TRUE)
  }
  
  window_name <- "unidentified_window"
  if (requireNamespace("xts", quietly = TRUE) && xts::is.xts(macro_data) && nrow(macro_data) > 0) {
    start_date <- format(zoo::index(macro_data)[1], "%Y-%m-%d")
    end_date <- format(zoo::index(macro_data)[nrow(macro_data)], "%Y-%m-%d")
    window_name <- paste0("window_", start_date, "_to_", end_date)
  }
  
  summary_file <- file.path(diag_dir, paste0(window_name, "_summary.txt"))
  tryCatch({
    model_summary <- summary(fitted_model)
    utils::capture.output(print(model_summary), file = summary_file)
    message(paste("RS-BVAR model summary saved to", summary_file))
  }, error = function(e) {
    warning(paste("Could not save RS-BVAR model summary:", e$message))
  })
  
  trace_plot_file <- file.path(diag_dir, paste0(window_name, "_trace_plots.pdf"))
  
  posterior_ok <- TRUE
  if (!is.null(fitted_model) && !is.null(fitted_model$posterior)) {
    for (param in names(fitted_model$posterior)) {
      if (any(!is.finite(fitted_model$posterior[[param]]))) {
        posterior_ok <- FALSE
        warning(paste("Non-finite values found in posterior draws for parameter:", param, ". Skipping trace plot generation for window ", window_name, "."))
        break
      }
    }
  } else {
    posterior_ok <- FALSE
    warning(paste("Posterior draws not found in fitted_model. Skipping trace plot generation for window ", window_name, "."))
  }

  if (posterior_ok) {
    tryCatch({
      grDevices::pdf(trace_plot_file, width = 8, height = 10)
      
      # Now rely on the generic plot function for bsvars objects for simplicity and robustness
      # Note to LLMs: there's no plot method for PosteriorBSVARMSH class object: 
      # regime probability plot is here only for a placeholder
      plot(compute_regime_probabilities(fitted_model))
      
      grDevices::dev.off()
      log_message(paste("RS-BVAR trace plots saved to", trace_plot_file), level = "INFO", app_config = app_config)
      
    }, error = function(e) {
      if(names(grDevices::dev.cur()) == "pdf") {
         grDevices::dev.off()
      }
      warning(paste("Could not save RS-BVAR trace plots for window ", window_name, ":", e$message, ". Creating skipped file."))
      file.create(sub("\\.pdf$", "_skipped_plot_error.txt", trace_plot_file))
    })
  } else {
    file.create(sub("\\.pdf$", "_skipped_non_finite_posterior.txt", trace_plot_file))
  }
  
  log_message("Attempting to normalize posterior draws.", level = "DEBUG", app_config = app_config)
  # Normalization for MSH models is more complex and depends on identifying restrictions.
  # The default identification might not be sufficient. For now, we skip explicit normalization
  # as the main issue is convergence. A robust identification strategy is a separate, complex task.

  
  log_message("fit_rsbvar_model finished successfully.", level = "DEBUG", app_config = app_config)
  return(list(fitted_model = fitted_model, spec = spec))
}

#' @title Fit a Dynamic Conditional Correlation (DCC) GARCH Model with t-Copula
#' @description Fits a DCC-GARCH model with GJR-GARCH(1,1) marginals and Student-t
#' distribution to asset returns, as specified in GEMINI.md.
#' @param asset_returns An xts object containing the asset returns time series.
#' @return A fitted `tsmarch` object.
fit_dcc_t_garch_model <- function(asset_returns, app_config) {
  log_message("Starting DCC-GARCH model fitting with tsmarch.", level = "DEBUG", app_config = app_config)
  
  diag_dir <- "output/dcc_diagnostics"
  if (!dir.exists(diag_dir)) {
    dir.create(diag_dir, recursive = TRUE)
  }
  
  asset_returns <- na.omit(asset_returns)
  
  log_message("Scaling asset returns by 100 for numerical stability.", level = "INFO", app_config = app_config)
  asset_returns <- asset_returns * 100
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  saveRDS(asset_returns, file.path(diag_dir, paste0("dcc_input_scaled_", timestamp, ".rds")))

  log_message(paste0("NROW(asset_returns) inside fit_dcc_t_garch_model: ", NROW(asset_returns)), level = "DEBUG", app_config = app_config)

  min_obs <- app_config$default$models$min_dcc_obs
  if (is.null(asset_returns) || NROW(asset_returns) < min_obs || NCOL(asset_returns) < 1) {
    log_message(paste0("Insufficient data for DCC-GARCH fitting (", NROW(asset_returns), " obs, need ", min_obs, "). Returning fallback."), level = "WARN", app_config = app_config)
    emp_mean <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) colMeans(asset_returns, na.rm = TRUE) else numeric(0)
    emp_cov <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) cov(as.matrix(asset_returns), use = "pairwise.complete.obs") else matrix(NA_real_)
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }

  emp_mean <- colMeans(asset_returns, na.rm = TRUE)
  emp_cov <- cov(as.matrix(asset_returns), use = "pairwise.complete.obs")
  
  # browser()
  garch_model_list <- lapply(colnames(asset_returns), function(col_name) {
    spec <- tsgarch::garch_modelspec(
      y = asset_returns[, col_name],
      model = "gjrgarch", # Keep this consistent for all assets
      constant = TRUE,
      distribution = "norm"
    )
    
    # Try standard estimation
    fit <- tryCatch({
      tsmethods::estimate(spec, keep_tmb = TRUE, stationarity_constraint = 0.95)
    }, error = function(e) NULL)
    
    # Check for NaNs (using your persistence_summary helper logic)
    is_failed <- is.null(fit) || any(is.nan(fit$persistence_summary[2])) || any(is.nan(fit$vcoef))
    
    if (is_failed) {
      message(paste("Cleaning Hessian for:", col_name))
      # RE-FIT with vcov = FALSE. This keeps the GJR structure (5 params)
      # but prevents the NaNs that crash the DCC and predict()
      fit <- tsmethods::estimate(spec, keep_tmb = TRUE, 
                                 stationarity_constraint = 0.95, vcov = FALSE)
    }
    return(fit)
  })
  
  garch_model <- to_multi_estimate(garch_model_list)
  names(garch_model) <- colnames(asset_returns)
  
  dcc_spec <- tsmarch::dcc_modelspec(
    garch_model,
    dynamics = "dcc",
    dcc_order = c(1, 1),
    distribution = "mvt"
  )
  
  # # 1. Force a minimum alpha_1 (shocks must update correlation)
  # dcc_spec$parmatrix[parameter == "alpha_1", lower := 0.015]
  # 
  # # 2. Cap beta_1 (ensure correlations don't have a unit root)
  # dcc_spec$parmatrix[parameter == "beta_1", upper := 0.88]
  # 
  # # 3. Ensure the shape (degrees of freedom) stays in a stable range
  # # (mvt with shape near 2.0 causes the Cholesky/predict crash)
  # dcc_spec$parmatrix[parameter == "shape", lower := 5.0]
   
  log_message("dcc_modelspec created.", level = "DEBUG", app_config = app_config)
  saveRDS(dcc_spec, file.path(diag_dir, paste0("dcc_spec_", timestamp, ".rds")))

  log_message(paste0("Calling estimate with data dimensions: ", paste(dim(asset_returns), collapse = ", ")), level = "DEBUG", app_config = app_config)
  
  if (anyNA(asset_returns)) {
    log_message("NA values found in asset_returns just before estimate. This should not happen. Returning fallback.", level = "ERROR", app_config = app_config)
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }
  if (any(is.infinite(asset_returns))) {
    log_message("Infinite values found in asset_returns just before estimate. Returning fallback.", level = "ERROR", app_config = app_config)
    return(list(fallback = TRUE, emp_mean = emp_cov, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }

  saveRDS(asset_returns, file.path(diag_dir, paste0("dcc_input_for_fit_", timestamp, ".rds")))

  # browser()
  dcc_fit_model <- estimate(dcc_spec)

  if (inherits(dcc_fit_model, "dcc.estimate")) {
      log_message("DCC-GARCH model fitting with tsmarch successful.", level = "DEBUG", app_config = app_config)
  } else {
    log_message("WARNING: tsmarch::estimate did not return a valid tsmarch.estimate object. A fallback using empirical moments will be used.", level = "WARN", app_config = app_config)
    emp_mean <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) colMeans(asset_returns, na.rm = TRUE) else numeric(0)
    emp_cov <- if (!is.null(asset_returns) && NCOL(asset_returns) >= 1) cov(as.matrix(asset_returns), use = "pairwise.complete.obs") else matrix(NA_real_)
    return(list(fallback = TRUE, emp_mean = emp_mean, emp_cov = emp_cov, asset_names = colnames(asset_returns)))
  }
  
  return(dcc_fit_model)
}
