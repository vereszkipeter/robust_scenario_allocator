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
download_and_cache_series <- function(ticker, source, cache_dir, from, to, force = FALSE) {
  
  # Ensure cache directory exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  
  sanitized_ticker <- gsub("\\^", "", ticker) # Remove '^' for filename
  cache_file <- file.path(cache_dir, paste0(sanitized_ticker, ".rds"))
  
  cached_obj <- NULL
  if (file.exists(cache_file)) {
    cached_obj <- load_cached_series(cache_file)
  }

  if (!is.null(cached_obj) && !force && cache_covers_request(cached_obj$meta, from, to, source)) {
    message("Loading ", ticker, " from cache (covers request).")
    data <- cached_obj$data
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
          data <- quantmod::Ad(data)
        } else if (NCOL(data) > 0) {
          data <- data[, NCOL(data)]
        } else {
          warning("No suitable data column found for ", ticker)
          data <- NULL
        }
        colnames(data) <- ticker
        data <- as.xts(data) # Explicitly ensure it's an xts object
      }, error = function(e) {
        warning("Failed to download Yahoo data for ", ticker, ": ", e$message)
        data <- NULL
      })
    } else if (source == "FRED") {
      if (Sys.getenv("FRED_API_KEY") == "") {
        warning("FRED_API_KEY environment variable is not set. FRED data download will likely fail.")
      }
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
      meta <- list(
        ticker = ticker,
        source = source,
        requested_from = as.character(from),
        requested_to = as.character(to),
        from = as.character(min(as.Date(index(data)))),
        to = as.character(max(as.Date(index(data)))),
        downloaded = as.character(Sys.time()),
        periodicity = tryCatch(periodicity(data)$scale, error = function(e) NA),
        n = NROW(data)
      )
      save_cached_series(data, meta, cache_file)
    }
  }
  return(data)
}

#' @title Get All Raw Financial and Macro Data
#' @description Downloads and caches all specified financial and macroeconomic time series
#' based on the provided asset metadata.
#' @param asset_metadata A tibble containing asset information, including ticker and source.
#' @param cache_dir A character string specifying the directory for caching data.
#' @param from A character string for the start date, "YYYY-MM-DD".
#' @param to A character string for the end date, "YYYY-MM-DD".
#' @return A list of xts objects, where each element is a raw data series for a ticker.
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

  if (length(all_data_list) == 0) {
    stop("No data could be downloaded or loaded.")
  }
  
  return(all_data_list) # Return the list of xts objects
}

#' @title Apply Transformations and Standardize Data
#' @description Applies the specified transformation functions from asset_metadata
#' to the raw time series data and standardizes all series to monthly frequency.
#' Handles potential NAs and aligns data.
#' @param merged_raw_data An xts object containing all raw data series merged by `get_all_raw_data`.
#' @param asset_metadata A tibble containing asset information, including ticker and transformation functions.
#' @param window_to_date The 'to' date of the current window, used to check for incomplete last months.
#' @return An xts object with all transformed and standardized monthly data.
apply_transformations <- function(raw_data_list, asset_metadata, window_to_date, monthly_label = "end", app_config) {

  # Convert each series in the list to monthly levels
  list_of_monthly_series <- lapply(names(raw_data_list), function(ticker_name) {
    x <- raw_data_list[[ticker_name]]
    
    if (NROW(x) == 0) {
      warning("Series ", ticker_name, " is empty. Skipping.")
      return(xts(order.by = as.Date(character(0)))) # Return empty xts
    }
    monthly_x <- to_monthly_levels(x, method = "last", label = monthly_label)
    # Ensure it's an xts object before returning from lapply
    if (!xts::is.xts(monthly_x)) {
      monthly_x <- as.xts(monthly_x)
    }
    return(monthly_x)
  })
  names(list_of_monthly_series) <- names(raw_data_list)

  merged_monthly <- tryCatch(
    {
      if (length(list_of_monthly_series) == 0) return(xts(order.by = as.Date(character(0)))) # Handle empty list case
      # Remove empty xts objects before merging to avoid issues
      non_empty_series <- list_of_monthly_series[sapply(list_of_monthly_series, NROW) > 0]
      if (length(non_empty_series) == 0) return(xts(order.by = as.Date(character(0))))
      
      Reduce(\(x, y) merge(x, y, all = TRUE), non_empty_series)
    },
    error = function(e) {
      message("ERROR: Failed to Reduce/merge monthly series. Error: ", e$message)
      stop(e) # Re-throw the error
    }
  )

  panel_monthly <- merged_monthly |>
    na.locf() |> # Fill NAs after merging
    apply.monthly(last) |> # Ensure strict monthly periodicity
    standardize_monthly_index(label = monthly_label)
  
  message("DEBUG: Summary of panel_monthly before individual transformations:")
  print(summary(panel_monthly))
  message("DEBUG: Structure of panel_monthly before individual transformations:")
  print(str(panel_monthly))

  # Ensure column names are set
  if (is.null(colnames(panel_monthly)) && length(colnames(merged_raw_data)) > 0) {
    colnames(panel_monthly) <- colnames(merged_raw_data)
  }

  transformed_list <- list()

  # Iterate through the tickers (columns of the panel_monthly)
  for (ticker in colnames(panel_monthly)) {
    row <- asset_metadata |> dplyr::filter(ticker == .env$ticker) # Find metadata for the current ticker
    if (nrow(row) == 0) {
      warning("Metadata not found for ticker: ", ticker, ". Skipping transformation.")
      next
    }
    transform_func <- row$transform[[1]] # Correctly extract the function object

    series <- panel_monthly[, ticker, drop = FALSE] # Extract individual series from the panel
    
    # Apply transformation function to the series from the monthly panel
    transformed_series <- tryCatch({
      # Pass app_config to the transform_func
      transform_func(series, app_config = app_config)
    }, error = function(e) {
      warning("Transformation failed for ticker: ", ticker, ". Error: ", e$message)
      return(xts(order.by = as.Date(character(0)))) # Return empty xts on error
    })

    # Ensure transformed_series is an xts object
    if (!xts::is.xts(transformed_series)) {
      warning("Transformation for ticker: ", ticker, " did not return an xts object. Converting to empty xts.")
      transformed_series <- xts(order.by = as.Date(character(0)))
    }

    if (NROW(transformed_series) > 0) {
      # Standardize monthly index to chosen convention (default month-end)
      transformed_series <- standardize_monthly_index(transformed_series, label = monthly_label)
      transformed_list[[ticker]] <- transformed_series
    } else {
      transformed_list[[ticker]] <- transformed_series # Keep empty if empty
    }
  }

    if (length(transformed_list) > 0) {
    # The panel_monthly was already merged and handled NAs.
    # Combine the transformed series robustly using pairwise Reduce to avoid
    # warnings from merge.xts when providing 'join' to do.call with many objects.
    if (length(transformed_list) == 1) {
      merged <- transformed_list[[1]]
    } else {
      merged <- Reduce(function(x, y) merge(x, y, join = "outer"), transformed_list)
    }

    # Ensure column names match tickers
    colnames(merged) <- names(transformed_list)
    # Remove columns that are entirely NA
    all_na <- sapply(as.list(merged), function(col) all(is.na(col)))
    if (any(all_na)) {
      warning("The following series are all NA after transformation and will be dropped: ", paste(names(merged)[all_na], collapse = ", "))
      merged <- merged[, !all_na, drop = FALSE]
    }
    return(merged)
  }
  # If no transformed series, return NULL
  return(NULL)
}
