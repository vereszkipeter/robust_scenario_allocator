library(rmgarch)
library(xts)

# 1. Create sample asset returns data
set.seed(123)
n_obs <- 200
n_assets <- 4
returns_data <- xts(matrix(rnorm(n_obs * n_assets, sd = 0.01), ncol = n_assets),
                    order.by = seq(as.Date("2020-01-01"), by = "day", length.out = n_obs))
colnames(returns_data) <- paste0("Asset", 1:n_assets)

# 2. Define and fit a DCC-GARCH model (simplified for demonstration)
uspec_list <- lapply(1:n_assets, function(i) {
  rugarch::ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    distribution.model = "norm" # Normal for simplicity
  )
})
multi_uspec <- rugarch::multispec(uspec_list)
dcc_spec <- rmgarch::dccspec(uspec = multi_uspec, dccOrder = c(1, 1), distribution = "mvnorm")
dcc_fit_model <- rmgarch::dccfit(dcc_spec, data = returns_data)

# --- NEW INSPECTION CODE ---
message("Class of dcc_fit_model (output of rmgarch::dccfit): ", paste(class(dcc_fit_model), collapse = ", "))
message("Slot names of dcc_fit_model: ", paste(slotNames(dcc_fit_model), collapse = ", "))

if ("model" %in% slotNames(dcc_fit_model)) {
  internal_model_slot_content <- dcc_fit_model@model
  message("Class of dcc_fit_model@model: ", paste(class(internal_model_slot_content), collapse = ", "))

  if (inherits(internal_model_slot_content, "list") && "modeldata" %in% names(internal_model_slot_content)) {
    modeldata_list <- internal_model_slot_content$modeldata
    message("Class of dcc_fit_model@model$modeldata: ", paste(class(modeldata_list), collapse = ", "))

    if (inherits(modeldata_list, "list") && "data" %in% names(modeldata_list)) {
      dcc_fit_data <- modeldata_list$data
      message("Structure of dcc_fit_data (dcc_fit_model@model$modeldata$data):")
      str(dcc_fit_data)

      message("\nClass of dcc_fit_data:")
      print(class(dcc_fit_data))

      message("\nColumn names of dcc_fit_data:")
      print(colnames(dcc_fit_data))
    } else {
      message("dcc_fit_model@model$modeldata is not a list or does not contain 'data'.")
    }
  } else {
    message("dcc_fit_model@model is not a list or does not contain 'modeldata'.")
  }
} else {
  message("dcc_fit_model does not have a 'model' slot.")
}
# --- END NEW INSPECTION CODE ---


# 3. Call dccsim with similar parameters
horizon <- 60 # Same horizon as in generate_full_scenarios
n_sim_paths <- 1 # For a single path simulation as in generate_full_scenarios

sim <- rmgarch::dccsim(
  dcc_fit_model,
  n.sim = horizon,
  m.sim = n_sim_paths,
  startMethod = "sample"
)

# Also save to RDS for inspection if needed
saveRDS(sim, "dccsim_investigation_output.rds")
message("\n'sim' object saved to dccsim_investigation_output.rds")
