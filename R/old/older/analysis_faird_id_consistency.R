
load("data_cleaned.RData")
etraps_daily <- data_list$etraps_daily  # Asignación

# Test with your data
faird_maize_original <- analyze_device_consistency_robust(
  etraps_daily, "FAIRD1", "FAIRD2", 
  metric = "biomass_mg", 
  log_transform = FALSE
)

faird_maize_log <- analyze_device_consistency_robust(
  etraps_daily, "FAIRD1", "FAIRD2", 
  metric = "biomass_mg", 
  log_transform = TRUE
)

# Test with your data
faird_meadow_original <- analyze_device_consistency_robust(
  etraps_daily, "FAIRD3", "FAIRD4", 
  metric = "biomass_mg", 
  log_transform = FALSE
)

faird_meadow_log <- analyze_device_consistency_robust(
  etraps_daily, "FAIRD3", "FAIRD4", 
  metric = "biomass_mg", 
  log_transform = TRUE
)

# Test with your data
id_maize_original <- analyze_device_consistency_robust(
  etraps_daily, "ID1", "ID2", 
  metric = "biomass_mg", 
  log_transform = FALSE
)

id_maize_log <- analyze_device_consistency_robust(
  etraps_daily, "ID1", "ID2", 
  metric = "biomass_mg", 
  log_transform = TRUE
)


# Compare results
cat("FAIR-D maize Original data - Residuals normal?", faird_maize_original$residuals_normal, "\n")
cat("FAIR-D maize Log-transformed - Residuals normal?", faird_maize_log$residuals_normal, "\n")
cat("FAIR-D meadow Original data - Residuals normal?", faird_meadow_original$residuals_normal, "\n")
cat("FAIR-D meadow Log-transformed - Residuals normal?", faird_meadow_log$residuals_normal, "\n")
cat("ID maize Original data - Residuals normal?", id_maize_original$residuals_normal, "\n")
cat("ID maize Log-transformed - Residuals normal?", id_maize_log$residuals_normal, "\n")

