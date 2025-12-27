# Helpers for reconstructing simulated series from transformed simulations
# Small, focused utilities to replace inline loops in generate_full_scenarios

#' Reconstruct a macro series from simulated transformed values
#' @param sim_series numeric vector of simulated transformed values (length = horizon)
#' @param original_transform character or function reference ("log_return", "mom_change", "identity")
#' @param last_level numeric scalar of the last historical level for the series
#' @return numeric vector of reconstructed levels (length = horizon)
reconstruct_macro_series <- function(sim_series, original_transform, last_level) {
  horizon <- length(sim_series)
  if (is.null(last_level) || !is.finite(last_level)) {
    warning("Missing or invalid last_level. Returning NA series.")
    return(rep(NA_real_, horizon))
  }

  if (is.function(original_transform)) {
    # Allow passing the actual function object (e.g., log_return)
    name <- deparse(substitute(original_transform))
  } else {
    name <- as.character(original_transform)
  }

  reconstructed <- numeric(horizon)

  if (grepl("log_return", name)) {
    multipliers <- exp(as.numeric(sim_series))
    reconstructed[1] <- last_level * multipliers[1]
    if (horizon > 1) {
      for (t in 2:horizon) reconstructed[t] <- reconstructed[t-1] * multipliers[t]
    }
  } else if (grepl("mom_change", name)) {
    multipliers <- 1 + as.numeric(sim_series)
    reconstructed[1] <- last_level * multipliers[1]
    if (horizon > 1) {
      for (t in 2:horizon) reconstructed[t] <- reconstructed[t-1] * multipliers[t]
    }
  } else if (grepl("identity", name) || grepl("level", name)) {
    # Simulated values are already levels/rates
    reconstructed <- as.numeric(sim_series)
  } else {
    # Fallback: try applying inverse transform if it's a function
    if (is.function(original_transform)) {
      # We expect an inverse_transform to be supplied elsewhere; best-effort fallback
      reconstructed <- as.numeric(sim_series)
      warning("Unknown original_transform; returning simulated values directly.")
    } else {
      reconstructed <- as.numeric(sim_series)
      warning("Unknown transform name; returning simulated values directly.")
    }
  }

  return(reconstructed)
}

#' Convert simulated log-returns (vector) to simple returns vector (e.g., 0.01 for 1%)
#' @param sim_series numeric vector of simulated log-returns
#' @return numeric vector of simple returns
reconstruct_asset_returns_from_log <- function(sim_series) {
  as.numeric(exp(as.numeric(sim_series)) - 1)
}
