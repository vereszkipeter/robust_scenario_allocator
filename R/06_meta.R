library(CVXR) # For convex optimization
library(dplyr) # For data manipulation

#' @title Perform Distributionally Robust Optimization (Meta-Learner)
#' @description Solves the meta-optimization problem to find optimal weights for
#'   the base strategies. This involves Distributionally Robust Optimization (DRO)
#'   over the strategies' P&L on scenarios, considering uncertainty in scenario
#'   probabilities within a Wasserstein ball. It also includes an entropy
#'   regularization on the meta-weights.
#' @param base_strategy_pnl_on_scenarios A 3D array (horizon, strategy, simulation)
#'   of cumulative P&L for each base strategy across all scenarios.
#' @param scenario_probabilities A numeric vector of adjusted probabilities for each scenario
#'   (e.g., from Entropy Pooling).
#' @param app_config A list containing application configuration, including DRO parameters.
#' @return A named numeric vector of optimal meta-weights for the base strategies.
perform_distributionally_robust_optimization <- function(base_strategy_pnl_on_scenarios, scenario_probabilities, app_config) {
  
  # Extract relevant dimensions and parameters
  horizon <- dim(base_strategy_pnl_on_scenarios)[1] - 1
  n_strategies <- dim(base_strategy_pnl_on_scenarios)[2]
  n_sim <- dim(base_strategy_pnl_on_scenarios)[3]
  strategy_names <- dimnames(base_strategy_pnl_on_scenarios)[[2]]
  
  dro_epsilon <- app_config$default$models$dro_epsilon
  entropy_penalty_kappa <- app_config$default$models$entropy_penalty_kappa
  alpha_cvar <- app_config$default$models$cvar_alpha

  terminal_pnl_strategies <- base_strategy_pnl_on_scenarios[horizon + 1, , ]
  
  if (is.null(dim(terminal_pnl_strategies))) {
    terminal_pnl_strategies <- matrix(terminal_pnl_strategies, nrow = n_strategies, ncol = n_sim)
  }

  # Define CVXR variables
  theta <- Variable(n_strategies, name = "meta_weights")
  zeta <- Variable(1, name = "cvar_value")
  alpha_dual <- Variable(1, name = "alpha_dual", pos = TRUE)
  beta_dual <- Variable(n_sim, name = "beta_dual", pos = TRUE)

  # Portfolio terminal P&L for each scenario.
  # We want an expression of dimension n_sim x 1.
  # terminal_pnl_strategies is n_strategies x n_sim. theta is n_strategies x 1.
  # So, t(terminal_pnl_strategies) %*% theta gives (n_sim x n_strategies) %*% (n_strategies x 1) -> n_sim x 1.
  portfolio_terminal_pnl_scenarios <- t(terminal_pnl_strategies) %*% theta
  
  # --- FIX: Vectorized DRO constraint ---
  dro_constraint <- beta_dual >= -portfolio_terminal_pnl_scenarios - zeta - alpha_dual
  
  # Constraints for meta-weights
  constraints <- list(
    sum(theta) == 1,
    theta >= 0,
    dro_constraint
  )
  
  # --- FIX: Correctly formulate entropy penalty using DCP-compliant sum_entr ---
  # The goal is to minimize the negative entropy: sum(theta * log(theta)).
  # CVXR's sum_entr(x) calculates sum(-x * log(x)).
  # Therefore, sum(theta * log(theta)) is equivalent to -sum_entr(theta).
  objective <- Minimize(
    zeta + (1 / (1 - alpha_cvar)) * (sum(scenario_probabilities * beta_dual) - alpha_dual * dro_epsilon) 
    - entropy_penalty_kappa * sum_entries(entr(theta))
  )
  
  problem <- Problem(objective, constraints)

  cvxr_res <- safe_solve_cvxr(problem, solvers = c("ECOS", "SCS"), allow_fallback = FALSE)

  optimal_weights <- as.numeric(cvxr_res$getValue(theta))
  # Guard: only assign names if lengths match
  if (!is.null(strategy_names) && length(optimal_weights) == length(strategy_names)) {
    names(optimal_weights) <- strategy_names
  } else {
    log_message(paste0("Strategy names length (", length(strategy_names), ") does not match optimal_weights length (", length(optimal_weights), "). Skipping names assignment."), level = "WARN", app_config = app_config)
  }
  message("DRO meta-optimization solved successfully using CVXR.")
  
  return(optimal_weights)
}