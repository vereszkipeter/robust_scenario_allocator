library(xts)
library(lubridate)

#' Convert a series to monthly levels using a specified method and index label
#' @param series An xts series
#' @param method One of "last" or "first" to pick period aggregation
#' @param label One of "end" or "start" to label the monthly timestamp
to_monthly_levels <- function(series, method = c("last", "first"), label = c("end", "start")) {
  method <- match.arg(method)
  label <- match.arg(label)
  if (!xts::is.xts(series)) stop("`series` must be an xts object")

  FUN <- if (method == "last") last else first
  monthly <- xts::apply.monthly(series, FUN = FUN)
  monthly <- as.xts(monthly) # Force it to be xts

  if (NROW(monthly) == 0) return(monthly)

  if (label == "end") {
    index(monthly) <- lubridate::ceiling_date(as.Date(index(monthly)), "month") - lubridate::days(1)
  } else {
    index(monthly) <- as.Date(format(index(monthly), "%Y-%m-01"))
  }
  return(monthly)
}

#' Standardize a monthly-indexed series to month start or month end
standardize_monthly_index <- function(series, label = c("end", "start")) {
  label <- match.arg(label)
  if (!xts::is.xts(series)) stop("`series` must be an xts object")
  if (NROW(series) == 0) return(series)
  if (label == "end") {
    index(series) <- lubridate::ceiling_date(as.Date(index(series)), "month") - lubridate::days(1)
  } else {
    index(series) <- as.Date(format(index(series), "%Y-%m-01"))
  }
  return(series)
}

#' Save a cached series object (data + meta)
save_cached_series <- function(series, meta, cache_file) {
  obj <- list(data = series, meta = meta)
  saveRDS(obj, cache_file)
  invisible(cache_file)
}

#' Load a cached series object (returns NULL if missing)
load_cached_series <- function(cache_file) {
  if (!file.exists(cache_file)) return(NULL)
  readRDS(cache_file)
}

#' Check whether a cache metadata covers requested range and source
cache_covers_request <- function(meta, from, to, source) {
  if (is.null(meta)) return(FALSE)
  if (!is.null(meta$source) && meta$source != source) return(FALSE)
  if (!is.null(meta$from) && !is.null(meta$to)) {
    return(as.Date(meta$from) <= as.Date(from) && as.Date(meta$to) >= as.Date(to))
  }
  return(FALSE)
}

#' @title Get a valid integer seed
#' @description Returns a valid integer seed. Prefer `tar_seed()` when available,
#' but fall back to a derived integer from current time. Always returns an integer.
get_valid_seed <- function(app_config) {
  s <- NA_integer_
  log_message("get_valid_seed() called.", level = "DEBUG", app_config = app_config)
  # try tar_seed_get() if available
  if (exists("tar_seed_get", mode = "function")) {
    try({
      s_try <- tar_seed_get()
      log_message(paste0("tar_seed_get() returned: ", s_try), level = "DEBUG", app_config = app_config)
      # Accept scalar numeric or integer-like seeds
      if (!is.null(s_try) && length(s_try) == 1 && !is.na(s_try) && is.finite(s_try)) {
        s <- as.numeric(s_try)
      }
    }, silent = TRUE)
  }
  if (is.na(s) || !is.finite(s)) {
    log_message("Falling back to Sys.time() for seed generation.", level = "DEBUG", app_config = app_config)
    # fallback to seconds since epoch
    s <- as.numeric(Sys.time())
  }
  
  # Ensure the seed is a single, positive integer within the valid 32-bit range.
  # 1. Take absolute value to handle negatives from tar_seed_get()
  # 2. Use modulo to wrap it into the positive integer range.
  # 3. Subtract 1 from integer.max to avoid boundary issues.
  # 4. Final as.integer coercion.
  final_seed <- as.integer(abs(s) %% (.Machine$integer.max - 1))

  # It must not be NA or zero
  if(is.na(final_seed) || final_seed == 0) {
      final_seed <- as.integer(as.numeric(Sys.time()) %% 1e9)
  }

  log_message(paste0("Final seed value: ", final_seed), level = "DEBUG", app_config = app_config)
  return(final_seed)
}


#' @title Split Data by Asset Class
#' @description Splits an xts object into two xts objects based on asset class
#' specified in the asset_metadata.
#' @param data_xts An xts object containing merged data.
#' @param asset_metadata A tibble with asset information, including 'ticker' and 'asset_class'.
#' @return A list with two xts objects: 'asset_returns' and 'macro_data'.
split_data_by_type <- function(data_xts, asset_metadata) {
  
  asset_tickers <- asset_metadata |> 
    filter(asset_class == "Asset") |> 
    pull(ticker)
  
  macro_tickers <- asset_metadata |> 
    filter(asset_class != "Asset") |> 
    pull(ticker)
  
  # Find which tickers are actually present in the data
  available_asset_tickers <- intersect(asset_tickers, colnames(data_xts))
  available_macro_tickers <- intersect(macro_tickers, colnames(data_xts))
  
  asset_returns <- data_xts[, available_asset_tickers, drop = FALSE]
  macro_data <- data_xts[, available_macro_tickers, drop = FALSE]
  
  return(list(asset_returns = asset_returns, macro_data = macro_data))
}