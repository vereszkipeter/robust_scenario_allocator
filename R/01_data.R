library(quantmod) # For Yahoo data
library(fredr)    # For FRED data
library(xts)      # For time series objects
library(dplyr)    # For general data manipulation, including lag if needed for other ops

#' @title Download and Cache a Single Financial Time Series
#' @description Downloads a single time series from Yahoo Finance or FRED,
#' caches it locally, and returns an xts object of the relevant series.
#' @param ticker A character string representing the asset ticker or FRED series ID.
#' @param source A character string, either "Yahoo" or "FRED".
#' @param cache_dir A character string specifying the directory for caching data.
#' @param from A character string for the start date, "YYYY-MM-DD".
#' @param to A character string for the end date, "YYYY-MM-DD".
#' @return An xts object containing the downloaded and processed time series.
#' @examples
#' # download_and_cache_series("SPY", "Yahoo", "data/raw", "2000-01-01", "2023-12-31")
#' # download_and_cache_series("CPIAUCSL", "FRED", "data/raw", "2000-01-01", "2023-12-31")
download_and_cache_series <- function(ticker, source, cache_dir, from, to) {
  
  # Ensure cache directory exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  
  sanitized_ticker <- gsub("\\^", "", ticker) # Remove '^' for filename
  cache_file <- file.path(cache_dir, paste0(sanitized_ticker, ".rds"))
  
  if (file.exists(cache_file)) {
    message("Loading ", ticker, " from cache.")
    data <- readRDS(cache_file)
  } else {
    message("Downloading ", ticker, " from ", source, "...")
    if (source == "Yahoo") {
      data_env <- new.env()
      tryCatch({
        data <- quantmod::getSymbols(
          ticker,
          src = "yahoo",
          from = from,
          to = to,
          auto.assign = FALSE,
          warnings = FALSE
        )
        if (NCOL(data) > 1 && paste0(ticker, ".Adjusted") %in% colnames(data)) {
            data <- quantmod::Ad(data) # For OHLCV, use Adjusted Price
        } else if (NCOL(data) > 0) {
            data <- data[, NCOL(data)] # For indices, often the last column (e.g., Close)
        } else {
            warning("No suitable data column found for ", ticker)
            data <- NULL
        }
        colnames(data) <- ticker # Clean column name
      }, error = function(e) {
        warning("Failed to download Yahoo data for ", ticker, ": ", e$message)
        data <- NULL
      })
    } else if (source == "FRED") {
      tryCatch({
        fred_data <- tryCatch(
          fredr::fredr_series_observations(
            series_id = ticker,
            observation_start = as.Date(from),
            observation_end = as.Date(to)
          ),
          error = function(e) {
            warning("Failed to retrieve FRED series for ", ticker, ": ", e$message)
            return(NULL)
          }
        )

        # Proceed only if fred_data was retrieved and contains rows
        if (is.null(fred_data) || NROW(fred_data) == 0) {
          warning("No FRED data found for ", ticker)
          data <- NULL
        } else {
          # FRED data is typically daily; convert to xts and select the 'value' column
          data <- xts(fred_data$value, order.by = as.Date(fred_data$date))
          colnames(data) <- ticker
        }
      }, error = function(e) {
        warning("Failed to download FRED data for ", ticker, ": ", e$message)
        data <- NULL
      })
    } else {
      stop("Unsupported data source: ", source)
    }
    
    if (!is.null(data)) {
      saveRDS(data, cache_file)
    }
  }
  return(data)
}

#' @title Get All Raw Financial and Macro Data
#' @description Downloads and caches all specified financial and macroeconomic time series
#' based on the provided asset metadata. Merges them into a single xts object.
#' @param asset_metadata A tibble containing asset information, including ticker and source.
#' @param cache_dir A character string specifying the directory for caching data.
#' @param from A character string for the start date, "YYYY-MM-DD".
#' @param to A character string for the end date, "YYYY-MM-DD".
#' @return An xts object with all merged raw data series.
get_all_raw_data <- function(asset_metadata, cache_dir, from, to) {
  
  all_data_list <- list()
  
  for (i in seq_len(nrow(asset_metadata))) {
    row <- asset_metadata[i, ]
    ticker <- row$ticker
    source <- row$source
    
    series_data <- download_and_cache_series(
      ticker = ticker,
      source = source,
      cache_dir = cache_dir,
      from = from,
      to = to
    )
    
    if (!is.null(series_data) && xts::is.xts(series_data)) { # Add xts::is.xts check
      all_data_list[[ticker]] <- series_data
    } else {
      warning("Series data for ", ticker, " is not a valid xts object or is NULL. Skipping.")
    }
  }

  # Merge all individual xts objects into a single one
  if (length(all_data_list) > 0) {
    # Use do.call(merge, ...) to merge multiple xts objects
    merged_data <- do.call(merge, all_data_list)
    # Remove rows with any NA values that might result from different start dates
    merged_data <- na.omit(merged_data)
  } else {
    stop("No data could be downloaded or loaded.")
  }
  
  return(merged_data)
}

#' @title Apply Transformations and Standardize Data
#' @description Applies the specified transformation functions from asset_metadata
#' to the raw time series data and standardizes all series to monthly frequency.
#' Handles potential NAs and aligns data.
#' @param raw_data An xts object containing all merged raw data series.
#' @param asset_metadata A tibble containing asset information, including ticker and transformation functions.
#' @param window_to_date The 'to' date of the current window, used to check for incomplete last months.
#' @return An xts object with all transformed and standardized monthly data.
apply_transformations <- function(raw_data, asset_metadata, window_to_date) {

  transformed_list <- list()

  for (i in seq_len(nrow(asset_metadata))) {
    row <- asset_metadata[i, ]
    ticker <- row$ticker
    transform_func <- row$transform[[1]] # Correctly extract the function object

    if (!ticker %in% colnames(raw_data)) {
      warning("Ticker ", ticker, " not found in raw_data. Skipping transformation.")
      next
    }

    series <- raw_data[, ticker] # This is the original daily/raw data for the ticker

    # Ensure data is monthly before applying transformations like log_return or mom_change
    # For Yahoo (prices), get monthly *levels* (last day of month price) first
    if (row$source == "Yahoo") {
      # Get the last trading day's price for each month
      monthly_levels <- apply.monthly(series, FUN = last)
    } else { # For FRED data
             # If FRED data is daily/weekly, convert it to monthly by taking the last value of the month.
        if (periodicity(series)$scale == "daily" || periodicity(series)$scale == "weekly") {
            monthly_levels <- apply.monthly(series, FUN = last)
        } else { # Already monthly or other suitable low frequency
            monthly_levels <- series
        }
    }

    # IMPORTANT: Check for incomplete last month and drop if necessary
    # `window_to_date` is the actual end date of the window, so we compare it with the end of the month
    if (!is.null(window_to_date) && NROW(monthly_levels) > 0) {
      last_date_in_monthly_levels <- index(tail(monthly_levels, 1))
      # Check if the window_to_date is before the true end of the month for the last observation
      # Use `ceiling_date` to find the end of the month containing window_to_date
      end_of_month_for_window_to_date <- lubridate::ceiling_date(as.Date(window_to_date), "month") - lubridate::days(1)

      if (as.Date(window_to_date) < end_of_month_for_window_to_date) {
        # If the last monthly observation's date is within the same month as window_to_date,
        # but window_to_date is not the actual end of that month, it's an incomplete month.
        if (lubridate::month(last_date_in_monthly_levels) == lubridate::month(window_to_date) &&
            lubridate::year(last_date_in_monthly_levels) == lubridate::year(window_to_date)) {
          monthly_levels <- head(monthly_levels, -1)
          message("Dropped incomplete last month for ticker: ", ticker, ". Last full month is up to: ", index(tail(monthly_levels, 1)))
        }
      }
    }
    
    # Apply transformation function to the monthly levels
    transformed_series <- transform_func(monthly_levels)
    
    # Ensure column name is correct
    colnames(transformed_series) <- ticker 
    transformed_list[[ticker]] <- transformed_series
  }
  
  if (length(transformed_list) > 0) {
    # Merge all transformed xts objects. They should be monthly and aligned.
    merged_transformed_data <- do.call(merge, transformed_list)
    # Remove any NAs introduced by transformations (e.g., lag) or merging (due to different start dates)
    merged_transformed_data <- na.omit(merged_transformed_data)
  } else {
    stop("No data could be transformed.")
  }
  
  return(merged_transformed_data)
}

#' @title Split Data into Macro Variables and Asset Returns
#' @description Splits a combined xts object of monthly returns into two separate
#' xts objects: one for macro variables and one for asset returns, based on
#' the asset metadata.
#' @param monthly_returns An xts object containing combined monthly returns of
#' assets and macro variables.
#' @param asset_metadata A tibble containing asset information, including ticker
#' and asset_class.
#' @return A list containing two xts objects: `macro_data` and `asset_returns`.
split_data_by_type <- function(monthly_returns, asset_metadata) {
  
  macro_tickers <- asset_metadata %>%
    dplyr::filter(asset_class == "Makro Változó") %>%
    dplyr::pull(ticker)
  
  asset_tickers <- asset_metadata %>%
    dplyr::filter(asset_class != "Makro Változó") %>%
    dplyr::pull(ticker)
  
  # Ensure all tickers are present in monthly_returns
  if (!all(macro_tickers %in% colnames(monthly_returns))) {
    missing_macro <- setdiff(macro_tickers, colnames(monthly_returns))
    warning("Missing macro tickers in monthly_returns: ", paste(missing_macro, collapse = ", "))
    # Remove missing tickers from macro_tickers list
    macro_tickers <- macro_tickers[macro_tickers %in% colnames(monthly_returns)]
  }
  if (!all(asset_tickers %in% colnames(monthly_returns))) {
    missing_assets <- setdiff(asset_tickers, colnames(monthly_returns))
    warning("Missing asset tickers in monthly_returns: ", paste(missing_assets, collapse = ", "))
    # Remove missing tickers from asset_tickers list
    asset_tickers <- asset_tickers[asset_tickers %in% colnames(monthly_returns)]
  }
  
  # Select columns, handling cases where some tickers might be missing
  macro_data <- if(length(macro_tickers) > 0) monthly_returns[, macro_tickers, drop = FALSE] else NULL
  asset_returns <- if(length(asset_tickers) > 0) monthly_returns[, asset_tickers, drop = FALSE] else NULL
  
  # Check if data frames are empty (changed from stop to warning and return NULL)
  if (is.null(macro_data) || NCOL(macro_data) == 0) {
    warning("No macro data found after splitting.")
    macro_data <- NULL # Ensure it's explicitly NULL if empty
  }
  if (is.null(asset_returns) || NCOL(asset_returns) == 0) {
    stop("No asset returns found after splitting. This is critical for the pipeline.")
  }

  return(list(macro_data = macro_data, asset_returns = asset_returns))
}

#' @title Load Project Configuration
#' @description A simple wrapper around config::get() to help with {targets}
#' dependency resolution.
#' @return A list containing the project configuration.
load_config <- function() {
  config::get(file = "config.yml")
}
