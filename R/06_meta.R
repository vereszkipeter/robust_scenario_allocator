library(CVXR)
library(dplyr)

#' @title Perform Distributionally Robust Optimization (Meta-Learner)
#' @description Solves the meta-optimization problem to find optimal weights for
#'   the base strategies. This function implements a Distributionally Robust 
#'   Optimization (DRO) framework using a Wasserstein uncertainty ball around 
#'   scenario probabilities. It minimizes the Conditional Value-at-Risk (CVaR) 
#'   of the portfolio loss, integrated with an entropy regularization term 
#'   to ensure diversification among strategy weights.
#' @param base_strategy_pnl_on_scenarios A 3D array (horizon + 1, strategy, simulation)
#'   of cumulative P&L for each base strategy.
#' @param scenario_probabilities A numeric vector of probabilities for each scenario
#'   (summing to 1).
#' @param app_config A list containing application configuration:
#'   - dro_epsilon: Radius of the Wasserstein ball.
#'   - entropy_penalty_kappa: Strength of the entropy regularization.
#'   - cvar_alpha: Confidence level for CVaR (e.g., 0.95).
#' @return A named numeric vector of optimal meta-weights for the base strategies.
perform_distributionally_robust_optimization <- function(base_strategy_pnl_on_scenarios, 
                                                         scenario_probabilities, 
                                                         app_config) {
  
  # --- Handle NULL or empty input scenarios ---
  if (is.null(base_strategy_pnl_on_scenarios) || 
      !is.array(base_strategy_pnl_on_scenarios) || 
      length(dim(base_strategy_pnl_on_scenarios)) < 3 || 
      dim(base_strategy_pnl_on_scenarios)[3] == 0) {
    return(NULL)
  }
  
  # --- Extract relevant dimensions and parameters ---
  # Assuming index horizon + 1 is the terminal P&L
  terminal_idx <- dim(base_strategy_pnl_on_scenarios)[1]
  n_strategies <- dim(base_strategy_pnl_on_scenarios)[2]
  n_sim        <- dim(base_strategy_pnl_on_scenarios)[3]
  strategy_names <- dimnames(base_strategy_pnl_on_scenarios)[[2]]
  
  if (n_strategies == 0 || n_sim == 0) {
    return(NULL)
  }
  
  # Configuration parameters
  dro_epsilon            <- app_config$default$models$dro_epsilon
  entropy_penalty_kappa  <- app_config$default$models$entropy_penalty_kappa
  alpha_cvar             <- app_config$default$models$cvar_alpha
  
  # Extract terminal P&L (strategies x simulations)
  terminal_pnl_strategies <- base_strategy_pnl_on_scenarios[terminal_idx, , ]
  
  # Ensure it's a matrix even if n_strategies = 1
  if (is.null(dim(terminal_pnl_strategies))) {
    terminal_pnl_strategies <- matrix(terminal_pnl_strategies, 
                                      nrow = n_strategies, 
                                      ncol = n_sim)
  }
  
  # --- CVXR Variables ---
  theta      <- Variable(n_strategies, name = "meta_weights")
  zeta       <- Variable(1, name = "cvar_threshold") # VaR level
  alpha_dual <- Variable(1, name = "wasserstein_dual", pos = TRUE) # Lambda in DRO
  beta_dual  <- Variable(n_sim, name = "scenario_slack", pos = TRUE) 
  
  # --- Portfolio Logic ---
  # We define loss as -P&L for CVaR minimization.
  # portfolio_loss is (n_sim x 1)
  portfolio_loss <- -t(terminal_pnl_strategies) %*% theta
  
  # --- Constraints ---
  # 1. Budget constraint: weights sum to 1
  # 2. No-shorting: weights are non-negative
  # 3. CVaR/DRO dual constraint: relates losses to the dual slack variables
  constraints <- list(
    sum(theta) == 1,
    theta >= 0,
    beta_dual >= portfolio_loss - zeta
  )
  
  # --- Objective Function ---
  # The objective combines three components:
  # 1. CVaR component: Threshold (zeta) + scaled expected exceedance
  # 2. DRO component: alpha_dual * epsilon (Penalty for uncertainty radius)
  #    Note: Must be POSITIVE to avoid unboundedness in minimization.
  # 3. Regularization: -kappa * Entropy (Maximizing entropy = minimizing negative entropy)
  #    entr(x) is concave, so -sum(entr(x)) is convex and DCP-compliant for Minimize().
  
  cvar_term    <- zeta + (1 / (1 - alpha_cvar)) * sum(scenario_probabilities * beta_dual)
  dro_term     <- alpha_dual * dro_epsilon
  entropy_term <- -entropy_penalty_kappa * sum(entr(theta))
  
  objective <- Minimize(cvar_term + dro_term + entropy_term)
  
  # --- Solve ---
  # ECOS is generally preferred for entropy/logarithmic problems
  problem  <- Problem(objective, constraints)
  cvxr_res <- safe_solve_cvxr(problem, solvers = c("ECOS", "SCS"), allow_fallback = TRUE)
  
  # --- Post-Processing ---
  if (cvxr_res$result$status != "optimal") {
    warning(paste("DRO solver did not reach optimal status. Status:", cvxr_res$status))
  }
  
  optimal_weights <- as.numeric(cvxr_res$result$getValue(theta))
  
  # Clean numerical noise (ensure exact non-negativity and sum to 1)
  optimal_weights[optimal_weights < 0] <- 0
  if (sum(optimal_weights) > 0) {
    optimal_weights <- optimal_weights / sum(optimal_weights)
  }
  
  if (!is.null(strategy_names) && length(optimal_weights) == length(strategy_names)) {
    names(optimal_weights) <- strategy_names
  }
  
  return(optimal_weights)
}