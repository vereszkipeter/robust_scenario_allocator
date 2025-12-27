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
  horizon <- dim(base_strategy_pnl_on_scenarios)[1] - 1 # P&L includes initial capital
  n_strategies <- dim(base_strategy_pnl_on_scenarios)[2]
  n_sim <- dim(base_strategy_pnl_on_scenarios)[3]
  strategy_names <- dimnames(base_strategy_pnl_on_scenarios)[[2]]
  
  # DRO parameters from app_config
  dro_epsilon <- app_config$default$models$dro_epsilon
  entropy_penalty_kappa <- app_config$default$models$entropy_penalty_kappa
  alpha_cvar <- app_config$default$models$cvar_alpha # Alpha for CVaR calculation

  # Extract terminal P&L for each strategy and each simulation
  # Dimensions: n_strategies x n_sim
  terminal_pnl_strategies <- base_strategy_pnl_on_scenarios[horizon + 1, , ]
  
  # Ensure terminal_pnl_strategies is a matrix for CVXR if only one strategy/sim
  if (is.null(dim(terminal_pnl_strategies))) {
    terminal_pnl_strategies <- matrix(terminal_pnl_strategies, nrow = n_strategies, ncol = n_sim)
  }

  # Define CVXR variables
  theta <- Variable(n_strategies, name = "meta_weights") # Weights for base strategies
  zeta <- Variable(1, name = "cvar_value")             # Auxiliary variable for CVaR
  alpha_dual <- Variable(1, name = "alpha_dual", pos = TRUE) # Dual variable for Wasserstein ball (lambda in formulation)
  beta_dual <- Variable(n_sim, name = "beta_dual", pos = TRUE) # Dual variables for each scenario (beta_s in formulation)

  # Portfolio terminal P&L for each scenario (vector of length n_sim)
  # Note: CVXR matrix multiplication ' %*% ' is automatically vectorized
  portfolio_terminal_pnl_scenarios <- t(theta) %*% terminal_pnl_strategies
  
  # Constraints for meta-weights
  constraints <- list(
    sum(theta) == 1,
    theta >= 0 # Long-only constraint
  )
  
  # Constraints for DRO dual problem
  # beta_dual[s] >= -(portfolio_terminal_pnl_scenarios[s] - zeta) - alpha_dual
  # This corresponds to the CVaR definition: Loss_s + zeta + alpha_dual <= beta_dual
  # Or: Loss_s - zeta - alpha_dual <= beta_dual
  # Where Loss_s = -portfolio_terminal_pnl_scenarios[s]
  # So: -portfolio_terminal_pnl_scenarios[s] - zeta - alpha_dual <= beta_dual
  
  for (s in 1:n_sim) {
    constraints <- c(constraints, 
                     beta_dual[s] >= -portfolio_terminal_pnl_scenarios[s] - zeta - alpha_dual)
  }
  
  # Objective: Minimize worst-case CVaR with entropy regularization
  # CVaR component: zeta + (1 / (1 - alpha_cvar)) * (sum(scenario_probabilities * beta_dual) - alpha_dual * dro_epsilon)
  # Entropy regularization: entropy_penalty_kappa * sum(theta * log(theta))
  # Note: sum(theta * log(theta)) is the negative entropy. CVXR's Entropy() function is -sum(x*log(x)).
  # So, for an entropy *penalty* (to encourage diversification), we want to add a positive term:
  # -entropy_penalty_kappa * sum(theta * log(theta)) or entropy_penalty_kappa * Entropy(theta).
  # We are minimizing, so adding Entropy(theta) encourages less entropy (more concentration).
  # To encourage diversification, we want to maximize entropy, so we MINIMIZE -Entropy(theta).
  # Hence: + entropy_penalty_kappa * sum(theta * log(theta)) (if sum(x*log(x)) is the definition of entropy)
  # CVXR's Entropy(x) is -sum(x*log(x)), so we want to add -entropy_penalty_kappa * Entropy(theta) to minimize.
  
  # The formulation should be: minimize worst_case_CVaR - entropy_penalty_kappa * Entropy(theta)
  # or minimize worst_case_CVaR + entropy_penalty_kappa * sum(theta * log(theta))
  
  # Let's use the explicit sum(theta*log(theta)) form for clarity of convexity in CVXR.
  # term sum(theta * log(theta)) is convex on theta > 0.
  
  # Correct objective function combining DRO-CVaR and Entropy Regularization
  objective <- Minimize(
    zeta + (1 / (1 - alpha_cvar)) * (sum(scenario_probabilities * beta_dual) - alpha_dual * dro_epsilon) + 
    entropy_penalty_kappa * sum(theta * log(theta))
  )
  
  # Formulate and solve the problem
  problem <- Problem(objective, constraints)

  # Use safe CVXR wrapper; fail fast on solver errors so issues surface during development
  cvxr_res <- safe_solve_cvxr(problem, solvers = c("ECOS", "SCS"), allow_fallback = FALSE)

  # Extract optimal meta-weights (will stop earlier if solver failed)
  optimal_weights <- as.numeric(cvxr_res$getValue(theta))
  names(optimal_weights) <- strategy_names
  message("DRO meta-optimization solved successfully using CVXR.")
  
  return(optimal_weights)
}