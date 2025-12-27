source('R/00_utils.R')
source('R/01_data.R')
source('R/02_models.R')

# Ensure config defaults
app_config <- tryCatch(load_config(), error = function(e) {message("load_config failed: ", e$message); list()})
print("app_config after load_config():")
print(app_config)
app_config <- ensure_app_config_defaults(app_config)
print("app_config after ensure_app_config_defaults():")
print(app_config)

# Re-set fred key for this session
fred_key <- tryCatch(keyring::key_get('FRED_API_KEY','Peter'), error=function(e) NA)
if(is.na(fred_key)) fred_key <- Sys.getenv('FRED_API_KEY', unset=NA)
if(!is.na(fred_key)) fredr::fredr_set_key(fred_key)

# Define asset_metadata, similar to _targets.R
asset_metadata <- tibble::tribble(
  ~ticker, ~source,   ~asset_class,     ~transform_name,      ~transform,         ~inverse_transform,
  app_config$default$data$tickers[1], "Yahoo", "Asset", "log_return",       log_return,         exp_return, # SPY
  app_config$default$data$tickers[2], "Yahoo", "Asset", "log_return",       log_return,         exp_return, # AGG
  app_config$default$data$tickers[3], "Yahoo", "Asset", "log_return",       log_return,         exp_return, # GLD
  app_config$default$data$tickers[4], "Yahoo", "Asset", "log_return",       log_return,         exp_return, # QQQ
  "CPIAUCSL", "FRED", "Makro Változó", "mom_change",       mom_change,       undo_mom_change, # Inflation
  "FEDFUNDS", "FRED", "Makro Változó", "identity_transform", identity_transform, identity_transform_inverse, # Rates
  "INDPRO",   "FRED", "Makro Változó", "mom_change",       mom_change,       undo_mom_change  # Growth
)

# Download/Load raw data
print(paste("Cache directory:", app_config$default$cache_dir))
raw <- get_all_raw_data(asset_metadata, app_config$default$cache_dir, app_config$default$data$from, app_config$default$data$to)
cat('Raw data class:', class(raw), '\n')

# Apply transformations
monthly <- apply_transformations(raw, asset_metadata)
cat('Monthly data class:', class(monthly), 'dims:', dim(monthly), '\n')

# Split
split <- split_data_by_type(monthly, asset_metadata)
cat('Macro data present?', !is.null(split$macro_data), '\n')
if(!is.null(split$macro_data)){
  print(head(split$macro_data))
}

# Try to call fit_rsbvar_model
if(!is.null(split$macro_data)){
  cat(paste0('Attempting to fit RS-BVAR with bvar_lags=', app_config$default$models$bvar_lags, '...\n'))
  mod <- tryCatch(
    fit_rsbvar_model(
      macro_data = split$macro_data,
      bvar_lags = app_config$default$models$bvar_lags,
      app_config = app_config
    ),
    error = function(e) {cat('RS-BVAR fit error:', e$message, '\n'); NULL}
  )
  if(!is.null(mod)) {
    cat('RS-BVAR fit ok\n')
    cat('Class of returned model:', class(mod), '\n')
    str(mod)
  } else {
    cat('RS-BVAR model is NULL.\n')
  }
}