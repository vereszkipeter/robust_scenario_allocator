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
    "ROI.plugin.quadprog",
    "coda"
  ),
  seed = 123 # Set a global seed for reproducibility
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
  tar_target(
    all_raw_data,
    get_all_raw_data(
      asset_metadata = asset_metadata,
      cache_dir = data_cache_dir,
      from = app_config$default$data$from,
      to = app_config$default$data$to
    )
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
      "VIXCLS", "FRED", "Makro Változó", "log_transform", base::log, base::exp, # VIX
      "CPIAUCSL", "FRED", "Makro Változó", "mom_change",       mom_change,       undo_mom_change, # Inflation
      "FEDFUNDS", "FRED", "Makro Változó", "identity_transform", diff_transform, diff_transform_inverse, # Rates
      "INDPRO",   "FRED", "Makro Változó", "mom_change",       mom_change,       undo_mom_change  # Growth
    )
  ),
  
        tar_target(
          panel_monthly_full,
          apply_transformations(
            raw_data_list = all_raw_data,
            asset_metadata = asset_metadata,
            # For the full panel, use the overall 'to' date from config as the window_to_date
            window_to_date = as.Date(app_config$default$data$to), 
            monthly_label = "end",
            app_config = app_config
          )
        ),  
  # Map the pipeline over each walk-forward window
  tar_target(
    wf_windows,
    generate_wf_windows(app_config, panel_monthly_full)
  ),

  tar_target(
    all_window_results_combined,
    purrr::pmap(
      wf_windows,
      function(window_id, from_date, to_date, val_date, oos_from_date, oos_to_date) {
        safely_process_window(
          window_id = window_id,
          from_date = from_date,
          to_date = to_date,
          val_date = val_date,
          oos_from_date = oos_from_date,
          oos_to_date = oos_to_date,
          app_config = app_config,
          asset_metadata = asset_metadata,
          all_raw_data = all_raw_data # Pass the new all_raw_data target
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
  )
  # ,
  # 
  # # Global Validation Report (aggregates results from all windows)
  # tar_target(
  #   global_validation_report,
  #   generate_global_validation_report(
  #     all_window_results = all_window_results_combined, # Pass the list of results
  #     app_config = app_config,
  #     output_report_dir = output_reports_dir,
  #     output_plots_dir = output_plots_dir_global
  #   ),
  #   format = "file"
  # )
)
