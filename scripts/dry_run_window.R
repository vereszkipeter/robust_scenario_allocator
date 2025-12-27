# Dry run: ingest + transforms for first walk-forward window
library(yaml)

# Source required helpers
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

# Load and tweak config for a quick dry run
app_config <- ensure_app_config_defaults(yaml::read_yaml("config.yml"))
app_config$default$n_simulations <- 10
app_config$default$horizon_months <- 3
app_config$default$models$bvar_lags <- 1

# Recreate asset_metadata similar to _targets.R
asset_metadata <- tibble::tribble(
  ~ticker, ~source,   ~asset_class,     ~transform,         ~inverse_transform,
  app_config$default$data$tickers[1], "Yahoo", "Asset", log_return,         exp_return,
  app_config$default$data$tickers[2], "Yahoo", "Asset", log_return,         exp_return,
  app_config$default$data$tickers[3], "Yahoo", "Asset", log_return,         exp_return,
  app_config$default$data$tickers[4], "Yahoo", "Asset", log_return,         exp_return,
  "CPIAUCSL", "FRED", "Makro Változó", mom_change,       undo_mom_change,
  "FEDFUNDS", "FRED", "Makro Változó", identity_transform, identity_transform_inverse,
  "INDPRO",   "FRED", "Makro Változó", mom_change,       undo_mom_change
)

# Generate windows
wf_windows <- generate_wf_windows(app_config)
print(wf_windows)

# Pick the first window for a focused dry run
w <- wf_windows[1, ]
print(w)

# 1) Data ingestion (from cache if available)
raw <- tryCatch(
  get_all_raw_data(asset_metadata = asset_metadata, cache_dir = app_config$default$cache_dir, from = w$from_date, to = w$to_date),
  error = function(e) { message('get_all_raw_data failed: ', e$message); NULL }
)
if (!is.null(raw)) {
  message('Raw data dims: ', paste(dim(raw), collapse = ' x '))
} else {
  message('Raw data is NULL')
}

# 2) Transformations
monthly <- NULL
if (!is.null(raw)) {
  monthly <- tryCatch(apply_transformations(raw_data = raw, asset_metadata = asset_metadata), error = function(e) { message('apply_transformations failed: ', e$message); NULL })
  if (!is.null(monthly)) message('Monthly transformed dims: ', paste(dim(monthly), collapse = ' x '))
}

# 3) Split
split <- NULL
if (!is.null(monthly)) {
  split <- tryCatch(split_data_by_type(monthly, asset_metadata), error = function(e) { message('split_data_by_type failed: ', e$message); NULL })
  if (!is.null(split)) {
    message('Macro cols: ', paste(colnames(split$macro_data), collapse = ', '))
    message('Asset cols: ', paste(colnames(split$asset_returns), collapse = ', '))
  }
}

# 4) Lag macro
if (!is.null(split$macro_data)) {
  lagged_macro <- tryCatch({ na.omit(xts::lag.xts(split$macro_data, k = -1)) }, error = function(e) { message('lagging macro failed: ', e$message); NULL })
  if (!is.null(lagged_macro)) message('Lagged macro dims: ', paste(dim(lagged_macro), collapse = ' x '))
}

# Save a small result for inspection
saveRDS(list(raw = raw, monthly = monthly, split = split), file = 'output/dry_run_window_result.rds')
message('Dry run complete; result saved to output/dry_run_window_result.rds')
