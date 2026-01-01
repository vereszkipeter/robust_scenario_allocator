library(CVXR) # For Entropy Pooling later
library(dplyr) # For data manipulation
library(xts)   # For time series objects
library(corpcor)  # For cov.shrink (Ledoit-Wolf shrinkage)

#' @title Calculate Implied Equilibrium Returns
#' @description Calculates implied equilibrium returns using a risk aversion parameter,
#' asset covariance, and equilibrium weights (e.g., ERC weights). This is based
#' on the Black-Litterman model's implied returns concept.
#' @param erc_weights A numeric vector of Equal Risk Contribution (ERC) weights.
#' @param asset_returns An xts object of historical asset returns.
#' @param risk_aversion_param A numeric value representing the risk aversion coefficient.
#' @return A named numeric vector of implied equilibrium excess returns.
calculate_implied_equilibrium_returns <- function(erc_weights, asset_returns, risk_aversion_param, app_config) {
  
  # Ensure weights are aligned with asset_returns columns
  # This function assumes erc_weights are already named correctly.
  
  # Calculate asset covariance matrix. Use Ledoit-Wolf shrinkage for robustness.
  # Ensure missing values are handled: drop rows with any NA before estimating covariance.
  if (any(is.na(asset_returns))) {
    asset_returns <- stats::na.omit(asset_returns)
    if (NROW(asset_returns) == 0) {
      log_message("All rows removed after na.omit on asset_returns; cannot compute covariance.", level = "ERROR", app_config = app_config)
      stop("All rows removed after na.omit on asset_returns; cannot compute covariance.")
    }
  }
  # cov.shrink from corpcor provides shrinkage estimation
  asset_cov_matrix <- corpcor::cov.shrink(as.matrix(asset_returns), verbose = FALSE)
  
  # Ensure asset_cov_matrix rows/cols align with erc_weights and asset_returns.
  # Assuming erc_weights are named.
  if (!all(names(erc_weights) %in% colnames(asset_cov_matrix))) {
    log_message("ERC weights names do not match covariance matrix column names.", level = "ERROR", app_config = app_config)
    stop("ERC weights names do not match covariance matrix column names.")
  }
  
  # Align erc_weights to the order of the covariance matrix
  aligned_erc_weights <- erc_weights[colnames(asset_cov_matrix)]
  
  # Calculate implied equilibrium excess returns: Pi = delta * Sigma * w_eq
  # where delta is risk_aversion_param, Sigma is asset_cov_matrix, w_eq is aligned_erc_weights
  implied_returns <- risk_aversion_param * asset_cov_matrix %*% aligned_erc_weights
  
  # The result `implied_returns` is a matrix (n_assets x 1). Convert to named vector.
  implied_returns <- as.numeric(implied_returns)
  names(implied_returns) <- colnames(asset_cov_matrix)
  
  log_message("Implied equilibrium returns calculated.", level = "DEBUG", app_config = app_config)
  return(implied_returns)
}


#' @title Model Yield Curve and Term Premium View
#' @description This function generates scenario-dependent term premium data.
#' It uses a simplified model where term premium varies with the short rate
#' from the provided `short_rate_scenarios`, along with configurable base
#' term premium and sensitivity parameters. This is an advanced placeholder
#' for a full yield curve model (e.g., Nelson-Siegel) but provides scenario-dependent
#' term premium for Entropy Pooling.
#' @param short_rate_scenarios A 3D array (n_sim, horizon, 1) of short rate scenarios.
#' @param app_config A list containing application configuration parameters.
#' @return A list containing `term_premium_scenarios` (3D array: n_sim, horizon, 1)
#'   and `term_premium_view_bounds` (list with min/max bounds).
model_yield_curve_and_term_premium_view <- function(short_rate_scenarios, app_config) {
  
  # Extract term premium bounds from app_config
  term_premium_min <- app_config$default$models$term_premium_min
  term_premium_max <- app_config$default$models$term_premium_max
  
  # Extract simplified term premium model parameters from app_config
  base_term_premium <- app_config$default$models$base_term_premium
  tp_sensitivity <- app_config$default$models$tp_sensitivity
  
  # Ensure short_rate_scenarios is not NULL
  if (is.null(short_rate_scenarios)) {
    log_message("short_rate_scenarios is NULL. Cannot generate scenario-dependent term premium. Returning NULL term_premium_scenarios.", level = "WARN", app_config = app_config)
    term_premium_view_bounds <- list(min = term_premium_min, max = term_premium_max)
    return(list(term_premium_scenarios = NULL, term_premium_view_bounds = term_premium_view_bounds))
  }

  n_sim <- dim(short_rate_scenarios)[1]
  horizon <- dim(short_rate_scenarios)[2]

  # Generate scenario-dependent term premium
  # Simplified model: TermPremium = base_term_premium + tp_sensitivity * ShortRate
  term_premium_scenarios <- array(NA, dim = c(n_sim, horizon, 1))
  
  for (s in 1:n_sim) {
    for (t in 1:horizon) {
      term_premium_scenarios[s, t, 1] <- base_term_premium + tp_sensitivity * short_rate_scenarios[s, t, 1]
    }
  }
  
  log_message("Generated scenario-dependent term premium (simplified model).", level = "DEBUG", app_config = app_config)
  
  term_premium_view_bounds <- list(
    min = term_premium_min,
    max = term_premium_max
  )
  
  # Return the scenario-dependent term premium and the static view bounds
  return(list(term_premium_scenarios = term_premium_scenarios, term_premium_view_bounds = term_premium_view_bounds))
}


#' @title Apply Sequential Entropy Pooling
#' @description Adjusts scenario probabilities using Entropy Pooling to incorporate views
#' and anchor them to implied equilibrium returns. This is an advanced placeholder
#' that performs a single-step Entropy Pooling using CVXR, approximating the concept
#' of sequential application by averaging views over the horizon.
#' @param simulated_scenarios A list containing `macro_scenarios` and `asset_scenarios`.
#' @param implied_equilibrium_returns A named numeric vector of implied equilibrium excess returns.
#' @param term_premium_model_output A list containing `term_premium_scenarios` and `term_premium_view_bounds`.
#' @param app_config A list containing application configuration parameters.
#' @return A numeric vector of adjusted scenario probabilities.
apply_sequential_entropy_pooling <- function(simulated_scenarios, implied_equilibrium_returns, term_premium_model_output, app_config) {
  
  # Extract relevant data
  asset_scenarios <- simulated_scenarios$asset_scenarios
  macro_scenarios <- simulated_scenarios$macro_scenarios # Not directly used for constraints here, but good to have
  n_sim <- dim(asset_scenarios)[1]
  horizon <- dim(asset_scenarios)[2]
  n_assets <- dim(asset_scenarios)[3]
  asset_names <- dimnames(asset_scenarios)[[3]]
  
  term_premium_scenarios <- term_premium_model_output$term_premium_scenarios
  term_premium_view_bounds <- term_premium_model_output$term_premium_view_bounds
  
  # --- CVXR-based Single-Step Entropy Pooling ---
  
  # 1. Prior Probabilities (Uniform)
  p_prior <- rep(1 / n_sim, n_sim)
  
  # 2. CVXR Variable for Adjusted Probabilities
  p_new <- Variable(n_sim, pos = TRUE, name = "adjusted_probabilities")
  
  # 3. Objective: Minimize Kullback-Leibler divergence (Entropy Pooling)
  # Use CVXR's kl_div atom for a DCP-safe formulation: kl_div(x, y) = x * log(x / y)
  objective <- Minimize(sum(kl_div(p_new, p_prior)))
  
  # 4. Constraints (Views)
  constraints <- list(
    sum(p_new) == 1, # Probabilities must sum to 1
    p_new >= 0 # Probabilities must be non-negative
  )
  
  # Prepare scenario data for constraints
  
  # Average Asset Returns per Scenario (across horizon)
  # Result: n_assets x n_sim matrix, where each column is an average return vector for one scenario
  avg_asset_returns_per_scenario <- matrix(NA, nrow = n_assets, ncol = n_sim)
  for (s in 1:n_sim) {
    # asset_scenarios[s, , ] is horizon x n_assets
    avg_asset_returns_per_scenario[, s] <- colMeans(asset_scenarios[s, , ])
  }
  
  # Implied Equilibrium Returns View (exact view: expected asset returns == implied_returns)
  # sum(p_new[s] * avg_asset_returns_per_scenario[asset_i, s]) == implied_equilibrium_returns[asset_i]
  # This can be written as matrix multiplication: avg_asset_returns_per_scenario %*% p_new == implied_equilibrium_returns
  
  # Ensure implied_equilibrium_returns are aligned with asset_names
  aligned_implied_returns <- implied_equilibrium_returns[asset_names]
  
  # Shrink the target returns to make the problem more feasible
  # Calculate the mean return of each asset from the simulations
  mean_sim_returns <- rowMeans(avg_asset_returns_per_scenario)
  names(mean_sim_returns) <- asset_names
  
  # Get shrinkage parameter from config
  alpha <- app_config$default$models$ep_shrinkage_alpha
  
  # Create a shrunk target vector
  shrunk_implied_returns <- alpha * aligned_implied_returns + (1 - alpha) * mean_sim_returns
  
  # Get tolerance for implied returns view from config
  ep_return_tolerance <- app_config$default$models$ep_return_tolerance

  # Replace the exact equality constraint with inequality constraints (range)
  # constraints <- c(constraints,
  #                  avg_asset_returns_per_scenario %*% p_new >= shrunk_implied_returns - ep_return_tolerance,
  #                  avg_asset_returns_per_scenario %*% p_new <= shrunk_implied_returns + ep_return_tolerance
  # )
  # message(paste0("Implied equilibrium returns view incorporated into Entropy Pooling with shrinkage and tolerance +/-", ep_return_tolerance, "."))
  
  # --- Handle Term Premium View if scenarios are available ---
  # if (!is.null(term_premium_scenarios)) {
  #   # Average Term Premium per Scenario (across horizon)
  #   avg_term_premium_per_scenario <- apply(term_premium_scenarios, 1, mean) # Result is vector of length n_sim
  #   
  #   # Term Premium View (range view: min <= expected term premium <= max)
  #   expected_term_premium <- sum(p_new * avg_term_premium_per_scenario)
  #   constraints <- c(constraints,
  #                    expected_term_premium >= term_premium_view_bounds$min,
  #                    expected_term_premium <= term_premium_view_bounds$max
  #   )
  #   message("Term Premium view incorporated into Entropy Pooling.")
  # } else {
  #   warning("Term Premium scenarios are NULL. Skipping Term Premium view in Entropy Pooling.")
  # }
  
  # Formulate and solve the problem
  problem <- Problem(objective, constraints)

  # Use safe CVXR wrapper; allow fallback to uniform probabilities on failure
  cvxr_res <- safe_solve_cvxr(problem, solvers = c("ECOS", "SCS"), allow_fallback = TRUE, app_config = app_config)

  # If CVXR returned an error message, include diagnostic details in the warning
  if (is.null(cvxr_res$result) && !is.null(cvxr_res$error)) {
    log_message(paste0("CVXR Entropy Pooling failed. Solver diagnostics: ", cvxr_res$error), level = "WARN", app_config = app_config)
  }

  # Check if CVXR solved successfully
  if (is.null(cvxr_res$result)) {
    log_message(paste0("CVXR Entropy Pooling failed to solve optimally or was infeasible. Returning uniform prior probabilities. Error: ", cvxr_res$error), level = "WARN", app_config = app_config)
    return(p_prior)
  }
  
  # If the solver returned a CVXR result object, extract probabilities
  adjusted_probabilities <- as.numeric(cvxr_res$result$getValue(p_new))
  log_message("Entropy Pooling solved successfully using CVXR.", level = "DEBUG", app_config = app_config)
  
  log_message("Applied (advanced placeholder) entropy pooling.", level = "INFO", app_config = app_config)
  return(adjusted_probabilities)
}
