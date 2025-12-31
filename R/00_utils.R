# Custom simple return function
simple_return <- function(x) {
  x / stats::lag(x) - 1
}

# Custom inverse simple return function
undo_simple_return <- function(x, last_level) {
  # The first value of x is the return from last_level to the first period.
  # So the first level is last_level * (1 + x[1])
  # The second level is first_level * (1 + x[2]), and so on.
  # This can be calculated with a cumulative product.
  
  # Ensure x is a numeric vector
  x_numeric <- as.numeric(x)
  
  # Replace NA with 0 for the cumulative product, assuming no return is 0 return
  x_numeric[is.na(x_numeric)] <- 0
  
  # Calculate cumulative returns
  cumulative_returns <- cumprod(1 + x_numeric)
  
  # Apply to the last level
  reconstructed_levels <- last_level * cumulative_returns
  
  return(reconstructed_levels)
}

# Identity function for transformations that do nothing
identity_transform <- function(x) {
  return(x)
}

identity_transform_inverse <- function(x, last_level) {
  return(x)
}

# Diff function for transformations that do nothing
diff_transform <- function(x) {
  return(diff(x))
}

diff_transform_inverse <- function(x, last_level) {
  return(cumsum(c(last_level, x)))
}

# Month-over-month change for inflation-like series (log-difference)
mom_change <- function(x) {
  # For inflation and industrial production, values are typically positive.
  log_x <- log(x)
  diff_log_x <- diff(log_x, lag = 1)
  # Ensure the output is an xts object with the correct index
  return(diff_log_x)
}

# Inverse of month-over-month log-difference
undo_mom_change <- function(sim_series, last_level) {
  reconstructed_log_x <- numeric(length(sim_series))
  
  # Start with the log of the last_level
  log_last_level <- log(last_level)

  # Reconstruct the log series
  reconstructed_log_x[1] <- sim_series[1] + log_last_level
  if (length(sim_series) > 1) {
    for (t in 2:length(sim_series)) {
      reconstructed_log_x[t] <- sim_series[t] + reconstructed_log_x[t-1]
    }
  }
  # Convert back to levels
  reconstructed_x <- exp(reconstructed_log_x)
  # Ensure output is an xts object if input sim_series was xts
  if (inherits(sim_series, "xts")) {
    reconstructed_x <- xts(reconstructed_x, order.by = index(sim_series))
  }
  return(reconstructed_x)
}

# Ensure app_config has all necessary defaults to avoid NULLs
ensure_app_config_defaults <- function(config) {
  
  # Set default for dcc_solver if not present
  if (is.null(config$default$models$dcc_solver)) {
    config$default$models$dcc_solver <- "solnp"
  }
  
  # Set default for dcc_solver_control if not present
  if (is.null(config$default$models$dcc_solver_control)) {
    config$default$models$dcc_solver_control <- list(trace = 0)
  }
  
  # Set a default for min_dcc_obs if not present
  if (is.null(config$default$models$min_dcc_obs)) {
    config$default$models$min_dcc_obs <- 100 # A reasonable default
  }

  # Set a default for log_level if not present
  if (is.null(config$default$log_level)) {
    config$default$log_level <- "INFO"
  }
  
  return(config)
}

# Logging utility
log_message <- function(message, level = "INFO", app_config) {
  log_levels <- c("DEBUG", "INFO", "WARN", "ERROR")
  current_level <- match(app_config$default$log_level, log_levels)
  message_level <- match(level, log_levels)
  
  if (message_level >= current_level) {
    cat(sprintf("[%s] %s: %s\n", Sys.time(), level, message))
  }
}

#' Safe CVXR solver wrapper
#' Tries a sequence of CVXR solvers and returns either the CVXR result (on success)
#' or a structured list when `allow_fallback = TRUE`.
#' @param problem CVXR Problem object
#' @param solvers character vector of solver names to try in order
#' @param allow_fallback logical; if TRUE return list(result=NULL,error=msg) on failure,
#'        if FALSE then stop on failure
#' @param verbose logical; print solver attempts
#' @param solver_args list; additional named args passed to `CVXR::solve`
#' @return If `allow_fallback = FALSE` returns CVXR result object on success.
#'         If `allow_fallback = TRUE` returns list(result = CVXR_result_or_NULL, error = NULL_or_message)
safe_solve_cvxr <- function(problem, solvers = c("ECOS", "SCS", "OSQP"), allow_fallback = TRUE, verbose = FALSE, solver_args = list()) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    msg <- "Package 'CVXR' is required for safe_solve_cvxr"
    if (allow_fallback) return(list(result = NULL, error = msg)) else stop(msg)
  }

  errors <- list()
  for (s in solvers) {
    if (verbose) message("Trying CVXR solver: ", s)
    attempt <- tryCatch({
      # Construct positional args: first the problem object, then solver and any solver_args
      args <- c(list(problem), list(solver = s), solver_args)
      # Call CVXR::solve positionally to avoid named-argument mismatches
      res <- do.call(CVXR::solve, args)
      res
    }, error = function(e) e)

    if (inherits(attempt, "error")) {
      errors[[s]] <- attempt$message
      next
    }

    status <- NULL
    if (!is.null(attempt$status)) status <- as.character(attempt$status)
    # Accept optimal or optimal_inaccurate as success
    if (!is.null(status) && tolower(status) %in% c("optimal", "optimal_inaccurate", "optimal_inaccurate")) {
      if (allow_fallback) return(list(result = attempt, error = NULL)) else return(attempt)
    }

    # If solver returned something but not optimal, still accept when allow_fallback=FALSE
    if (!is.null(status) && !allow_fallback) {
      # return the result and let caller decide; this mirrors CVXR behavior
      return(attempt)
    }

    # Record non-error solver result as warning
    errors[[s]] <- paste0("status=", ifelse(is.null(status), "UNKNOWN", status))
  }

  # All solvers failed or none returned acceptable status
  msg <- paste(names(errors), unlist(errors), sep = ": ", collapse = "; ")
  if (allow_fallback) return(list(result = NULL, error = msg)) else stop(msg)
}
