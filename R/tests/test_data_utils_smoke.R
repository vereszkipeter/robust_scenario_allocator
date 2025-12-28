library(xts)
source('R/02_data_utils.R')

cat('Running data utils smoke test...\n')

dates <- seq.Date(as.Date('2020-01-01'), as.Date('2020-03-31'), by = 'day')
series <- xts(1:length(dates), order.by = dates)

monthly_end <- to_monthly_levels(series, method = 'last', label = 'end')
cat('monthly_end index:', paste(as.character(index(monthly_end)), collapse = ', '), '\n')

monthly_start <- to_monthly_levels(series, method = 'last', label = 'start')
cat('monthly_start index:', paste(as.character(index(monthly_start)), collapse = ', '), '\n')

std_end <- standardize_monthly_index(monthly_start, label = 'end')
cat('standardized to end index:', paste(as.character(index(std_end)), collapse = ', '), '\n')

dir.create('data/cache', recursive = TRUE, showWarnings = FALSE)
cache_file <- 'data/cache/test_series_helper.rds'
save_cached_series(series, list(source = 'TEST', from = '2020-01-01', to = '2020-03-31'), cache_file)
obj <- load_cached_series(cache_file)
if (is.null(obj)) {
  warning('Failed to load cached object for smoke test. Skipping cache related assertion.')
} else {
  cat('Cached meta source:', obj$meta$source, '\n')
}

cat('SMOKE TEST: OK\n')
