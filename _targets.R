library(targets)
library(tarchetypes)
library(lubridate) # For date manipulation in window generation

# Load functions from R directory
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

# Set global options for targets, including required packages
tar_option_set(
  packages = c(
    "quantmod",
    "xts",
    "zoo",
    "lubridate",
    "tidyverse",
    "yaml",
    "rugarch",
    "rmgarch",
    "bsvars",
    "MASS",
    "patchwork",
    "forecast",
    "fredr",
    "tibble",
    "PortfolioAnalytics",
    "HierPortfolios",
    "riskParityPortfolio",
    "PerformanceAnalytics",
    "CVXR",
    "Hmisc",
    "corpcor",
    "ROI.plugin.glpk",
    "ROI.plugin.quadprog"
  )
)

# Main pipeline definition for targets

# Main pipeline definition for targets
list(
  # 1. Configuration (Global)
  tar_target(
    app_config,
    ensure_app_config_defaults(yaml::read_yaml("config.yml")),
    cue = tar_cue("always")
  ),
  
  # Ensure cache directory is created from config (Global)
  tar_target(
    data_cache_dir,
    {
      if (!dir.exists(app_config$default$cache_dir)) {
        dir.create(app_config$default$cache_dir, recursive = TRUE)
      }
      app_config$default$cache_dir
    }
  ),

  tar_target(
    asset_metadata,
    tibble::tribble(
      ~ticker, ~source,   ~asset_class,     ~transform_name,      ~transform,         ~inverse_transform,
      app_config$default$data$tickers[1], "Yahoo", "Asset", "simple_return",      simple_return,      undo_simple_return, # SPY
      app_config$default$data$tickers[2], "Yahoo", "Asset", "simple_return",      simple_return,      undo_simple_return, # AGG
      app_config$default$data$tickers[3], "Yahoo", "Asset", "simple_return",      simple_return,      undo_simple_return, # GLD
      app_config$default$data$tickers[4], "Yahoo", "Asset", "simple_return",      simple_return,      undo_simple_return, # QQQ
      "VIXCLS", "FRED", "Makro Változó", "log_transform", base::log, base::exp, # VIX
      "CPIAUCSL", "FRED", "Makro Változó", "mom_change",       mom_change,       undo_mom_change, # Inflation
      "FEDFUNDS", "FRED", "Makro Változó", "identity_transform", identity_transform, identity_transform_inverse, # Rates
      "INDPRO",   "FRED", "Makro Változó", "mom_change",       mom_change,       undo_mom_change  # Growth
    )
  ),
  
  # Map the pipeline over each walk-forward window
  tar_target(
    wf_windows,
    generate_wf_windows(app_config)
  ),

  tar_target(
    all_window_results_combined,
    purrr::pmap(
      wf_windows,
      function(window_id, from_date, to_date, val_date, oos_from_date, oos_to_date) {
        process_window(
          window_id = window_id,
          from_date = from_date,
          to_date = to_date,
          val_date = val_date,
          oos_from_date = oos_from_date,
          oos_to_date = oos_to_date,
          app_config = app_config,
          asset_metadata = asset_metadata
        )
      }
    )
  ),

  tar_target(
    output_reports_dir,
    {dir.create("output/reports", recursive = TRUE, showWarnings = FALSE); "output/reports"},
    format = "file"
  ),
  tar_target(
    output_plots_dir_global, # A separate global plot dir for aggregated plots
    {dir.create("output/plots/global", recursive = TRUE, showWarnings = FALSE); "output/plots/global"},
    format = "file"
  ),
  
  # Global Validation Report (aggregates results from all windows)
  tar_target(
    global_validation_report,
    generate_global_validation_report(
      all_window_results = all_window_results_combined, # Pass the list of results
      app_config = app_config,
      output_report_dir = output_reports_dir,
      output_plots_dir = output_plots_dir_global
    ),
    format = "file"
  )
)