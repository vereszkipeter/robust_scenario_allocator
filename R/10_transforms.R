# Helpers for reconstructing simulated series from transformed simulations
# Small, focused utilities to replace inline loops in generate_full_scenarios

#' Reconstruct a macro series from simulated transformed values
#' @param sim_series numeric vector of simulated transformed values (length = horizon)
#' @param original_transform character or function reference ("log_return", "mom_change", "identity")
#' @param last_level numeric scalar of the last historical level for the series
#' @return numeric vector of reconstructed levels (length = horizon)
reconstruct_macro_series <- function(sim_series, original_transform, last_level, app_config) {
  horizon <- length(sim_series)
  if (is.null(last_level) || !is.finite(last_level)) {
    log_message("Missing or invalid last_level. Returning NA series.", level = "WARN", app_config = app_config)
    return(rep(NA_real_, horizon))
  }

  if (is.function(original_transform)) {
    # Allow passing the actual function object (e.g., log_return)
    name <- deparse(substitute(original_transform))
  } else {
    name <- as.character(original_transform)
  }

  reconstructed <- numeric(horizon)
  lname <- tolower(name)

  # Log-level transform (e.g., log_transform -> inverse is exp)
  if (grepl("log", lname) && grepl("transform", lname) || grepl("log[_\\.-]?trans", lname)) {
    reconstructed <- as.numeric(exp(as.numeric(sim_series)))

  # Log-returns (many naming variants: log_return, log.return, log_returns, etc.)
  } else if ((grepl("log", lname) && grepl("return", lname)) || grepl("log[_\\.-]?ret", lname)) {
    multipliers <- exp(as.numeric(sim_series))
    reconstructed[1] <- last_level * multipliers[1]
    if (horizon > 1) {
      for (t in 2:horizon) reconstructed[t] <- reconstructed[t-1] * multipliers[t]
    }

  # Month-over-month or annual changes represented as relative changes
  } else if (grepl("mom", lname) || (grepl("month", lname) && grepl("change", lname)) || (grepl("annual", lname) && grepl("change", lname)) || grepl("_change", lname)) {
    multipliers <- 1 + as.numeric(sim_series)
    reconstructed[1] <- last_level * multipliers[1]
    if (horizon > 1) {
      for (t in 2:horizon) reconstructed[t] <- reconstructed[t-1] * multipliers[t]
    }

  # Identity / already-level series (rates or levels)
  } else if (grepl("identity", lname) || grepl("^level$", lname) || grepl("rate", lname)) {
    reconstructed <- as.numeric(sim_series)

  } else {
    # Fallback: if original_transform was a function and an inverse exists, try it
    if (is.function(original_transform)) {
      # Try common inverse naming conventions e.g., undo_<name>
      inv_name1 <- paste0("undo_", name)
      inv_name2 <- paste0("undo_", lname)
      inv_fun <- NULL
      if (exists(inv_name1, mode = "function", inherits = TRUE)) inv_fun <- get(inv_name1, mode = "function", inherits = TRUE)
      if (is.null(inv_fun) && exists(inv_name2, mode = "function", inherits = TRUE)) inv_fun <- get(inv_name2, mode = "function", inherits = TRUE)

      if (!is.null(inv_fun) && is.function(inv_fun)) {
        # Call inverse function with best-effort signature
        reconstructed <- tryCatch({
          inv_fun(sim_series, last_level)
        }, error = function(e) {
          log_message(paste0("Inverse transform function '", inv_name1, "'/'", inv_name2, "' failed; returning simulated values. Error: ", e$message), level = "WARN", app_config = app_config)
          as.numeric(sim_series)
        })
      } else {
        reconstructed <- as.numeric(sim_series)
        log_message(paste0("Unknown original_transform; returning simulated values directly. (tried names: ", name, ", ", lname, ")"), level = "WARN", app_config = app_config)
      }
    } else {
      reconstructed <- as.numeric(sim_series)
      log_message(paste0("Unknown transform name; returning simulated values directly. (", name, ")"), level = "WARN", app_config = app_config)
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
