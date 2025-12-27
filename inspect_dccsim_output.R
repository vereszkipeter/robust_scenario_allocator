library(rmgarch)

# Load the saved sim object
sim <- readRDS("dccsim_investigation_output.rds")

message("Is 'sim' NULL? ", is.null(sim))
if (!is.null(sim) && inherits(sim, "DCCsim")) {
  message("Class of 'sim': ", paste(class(sim), collapse = ", "))
  message("Slot names of 'sim': ", paste(slotNames(sim), collapse = ", "))

  if (!is.null(sim@msim)) {
    message("Is 'sim@msim' NULL? ", is.null(sim@msim))
    message("Names of 'sim@msim': ", paste(names(sim@msim), collapse = ", "))

    # New section to specifically check simX
    if (!is.null(sim@msim$simX)) {
      message("\nIs 'sim@msim$simX' NULL? ", is.null(sim@msim$simX))
      message("Class of 'sim@msim$simX': ", paste(class(sim@msim$simX), collapse = ", "))
      message("Dimensions of 'sim@msim$simX': ", paste(dim(sim@msim$simX), collapse = ", "))
      
      str_output_simX <- capture.output(str(sim@msim$simX))
      message("\nStructure of sim@msim$simX:\n", paste(str_output_simX, collapse = "\n"))
    } else {
      message("\n'sim@msim$simX' is NULL.")
    }

    if (!is.null(sim@msim$series) && length(sim@msim$series) > 0) {
      message("Is 'sim@msim$series' NULL? ", is.null(sim@msim$series))
      message("Length of 'sim@msim$series': ", length(sim@msim$series))

      str_output_series <- capture.output(str(sim@msim$series))
      message("\nStructure of sim@msim$series:\n", paste(str_output_series, collapse = "\n"))

      if (!is.null(sim@msim$series[[1]])) {
        message("\nIs 'sim@msim$series[[1]]' NULL? ", is.null(sim@msim$series[[1]]))
        message("Length of 'sim@msim$series[[1]]': ", length(sim@msim$series[[1]]))
        message("Dimensions of 'sim@msim$series[[1]]': ", paste(dim(sim@msim$series[[1]]), collapse = ", "))

        str_output_series1 <- capture.output(str(sim@msim$series[[1]]))
        message("\nStructure of sim@msim$series[[1]]:\n", paste(str_output_series1, collapse = "\n"))
      } else {
        message("\n'sim@msim$series[[1]]' is NULL.")
      }
    } else {
      message("\n'sim@msim$series' is NULL or empty (checked again).")
    }
  } else {
    message("\n'sim@msim' is NULL.")
  }
} else {
  message("\n'sim' is NULL or not a DCCsim object.")
}
