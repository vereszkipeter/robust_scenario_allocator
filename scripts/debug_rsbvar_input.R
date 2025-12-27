source('R/00_utils.R')
source('R/01_data.R')
source('R/02_models.R')

# Ensure config defaults
app_config <- tryCatch(load_config(), error = function(e) {list()})
app_config <- ensure_app_config_defaults(app_config)

# Re-set fred key for this session
fred_key <- tryCatch(keyring::key_get('FRED_API_KEY','Peter'), error=function(e) NA)
if(is.na(fred_key)) fred_key <- Sys.getenv('FRED_API_KEY', unset=NA)
if(!is.na(fred_key)) fredr::fredr_set_key(fred_key)

# Download/Load raw data
raw <- get_all_raw_data(asset_metadata, 'data/cache', '2000-01-01', '2023-12-31')
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
  cat('Attempting to fit RS-BVAR with sample lags=1 (dry-run)...\n')
  mod <- tryCatch(fit_rsbvar_model(split$macro_data, bvar_lags = 1), error = function(e) {cat('RS-BVAR fit error:', e$message, '\n'); NULL})
  if(!is.null(mod)) cat('RS-BVAR fit ok\n')
}
