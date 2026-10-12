# Test si los datos DIARIOS ya son normales sin transformar
load("data_cleaned.RData")
etraps_daily <- data_list$etraps_daily

# FAIRD1 vs FAIRD2 - datos diarios ORIGINALES
daily_biomass <- etraps_daily %>%
  filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  select(Date, Device, biomass_mg) %>%
  pivot_wider(names_from = Device, values_from = biomass_mg, values_fill = 0)

# Test normalidad en RESIDUOS (no en datos brutos)
model_original <- lm(FAIRD2 ~ FAIRD1, data = daily_biomass)
shapiro_original <- shapiro.test(residuals(model_original))

cat("Shapiro-Wilk test on RESIDUALS (original scale):\n")
cat("p-value =", shapiro_original$p.value, "\n")
cat("Normal?", shapiro_original$p.value > 0.05, "\n")

# Compare con log-transformed
daily_biomass_log <- daily_biomass %>%
  mutate(
    FAIRD1_log = log10(FAIRD1 + 1),
    FAIRD2_log = log10(FAIRD2 + 1)
  )

model_log <- lm(FAIRD2_log ~ FAIRD1_log, data = daily_biomass_log)
shapiro_log <- shapiro.test(residuals(model_log))

cat("\nShapiro-Wilk test on RESIDUALS (log scale):\n")
cat("p-value =", shapiro_log$p.value, "\n")
cat("Normal?", shapiro_log$p.value > 0.05, "\n")