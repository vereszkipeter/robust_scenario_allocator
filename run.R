#!/usr/bin/env Rscript

# This script runs the {targets} pipeline.
# It automatically sources the functions in the R/ directory
# and executes the pipeline defined in _targets.R.

# Load the targets library
library(targets)
library(fredr) # Load fredr package

# Set FRED API key: try keyring first, then fallback to environment variable `FRED_API_KEY`.
# Warn if no key is found so users know certain downloads may fail.
fred_key <- tryCatch(
	keyring::key_get("FRED_API_KEY", "Peter"),
	error = function(e) NA_character_
)
if (is.na(fred_key) || nchar(fred_key) == 0) {
	fred_key <- Sys.getenv("FRED_API_KEY", unset = NA_character_)
}
if (is.na(fred_key) || nchar(fred_key) == 0) {
	warning("FRED API key not found in keyring or FRED_API_KEY env var. Some FRED downloads may fail.")
} else {
	fredr::fredr_set_key(fred_key)
}

# Run the pipeline
cat("--- Running the Robust Scenario Allocator Pipeline ---\n")
targets::tar_make()
cat("--- Pipeline execution finished. ---\n")
