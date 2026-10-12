# ==============================================================================
# RIGOROUS REGRESSION DIAGNOSTICS - CORRECTED
# ==============================================================================
# Fixed: Explicit namespace to avoid MASS::select() masking dplyr::select()
# ==============================================================================

library(tidyverse)
library(car)  # For influence diagnostics
library(MASS) # For robust regression - LOADED LAST to avoid masking issues

load("data_cleaned.RData")
etraps_daily <- data_list$etraps_daily

# ==============================================================================
# PREPARE DATA - Use explicit dplyr:: namespace
# ==============================================================================

cat("\n=== PREPARING DATA ===\n")

counts_data <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  dplyr::select(Date, Device, abundance) %>%  # EXPLICIT dplyr::
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

cat("Data prepared:", nrow(counts_data), "observations\n")
cat("Date range:", min(counts_data$Date), "to", max(counts_data$Date), "\n\n")

# ==============================================================================
# 1. STANDARD LINEAR REGRESSION
# ==============================================================================

cat("\n=== ORIGINAL SCALE DIAGNOSTICS ===\n")

model_original <- lm(FAIRD2 ~ FAIRD1, data = counts_data)

cat("R² =", round(summary(model_original)$r.squared, 3), "\n")
cat("Adjusted R² =", round(summary(model_original)$adj.r.squared, 3), "\n\n")

# Residual normality
shapiro_orig <- shapiro.test(residuals(model_original))
cat("Shapiro-Wilk (residuals):", round(shapiro_orig$p.value, 4), 
    ifelse(shapiro_orig$p.value > 0.05, "✓ Normal", "✗ NOT normal"), "\n")

# Homoscedasticity (Breusch-Pagan test)
bp_test <- car::ncvTest(model_original)
cat("Breusch-Pagan (homoscedasticity):", round(bp_test$p, 4),
    ifelse(bp_test$p > 0.05, "✓ Homoscedastic", "✗ Heteroscedastic"), "\n\n")

# ==============================================================================
# 2. INFLUENTIAL POINTS ANALYSIS
# ==============================================================================

cat("=== INFLUENTIAL POINTS (Cook's Distance) ===\n")

# Cook's Distance
cooks_d <- cooks.distance(model_original)
threshold_cooks <- 4/nrow(counts_data)  # Standard threshold
influential <- which(cooks_d > threshold_cooks)

counts_data$cooks_d <- cooks_d
counts_data$influential <- cooks_d > threshold_cooks

cat("Number of influential points:", length(influential), "/", nrow(counts_data), "\n")
cat("Cook's D threshold:", round(threshold_cooks, 4), "\n")

if(length(influential) > 0) {
  cat("\nInfluential observations:\n")
  print(counts_data[influential, c("Date", "FAIRD1", "FAIRD2", "cooks_d")])
}

# Leverage
leverage <- hatvalues(model_original)
threshold_leverage <- 2*2/nrow(counts_data)  # 2p/n threshold
high_leverage <- which(leverage > threshold_leverage)
cat("\nNumber of high-leverage points:", length(high_leverage), "\n")
cat("Leverage threshold:", round(threshold_leverage, 4), "\n\n")

# ==============================================================================
# 3. LOG-TRANSFORMED REGRESSION
# ==============================================================================

cat("=== LOG-TRANSFORMED DIAGNOSTICS ===\n")

counts_log <- counts_data %>%
  dplyr::mutate(
    FAIRD1_log = log10(FAIRD1 + 1),
    FAIRD2_log = log10(FAIRD2 + 1)
  )

model_log <- lm(FAIRD2_log ~ FAIRD1_log, data = counts_log)

cat("R² =", round(summary(model_log)$r.squared, 3), "\n")
cat("Adjusted R² =", round(summary(model_log)$adj.r.squared, 3), "\n\n")

# Residual normality
shapiro_log <- shapiro.test(residuals(model_log))
cat("Shapiro-Wilk (residuals):", round(shapiro_log$p.value, 4),
    ifelse(shapiro_log$p.value > 0.05, "✓ Normal", "✗ NOT normal"), "\n")

# Homoscedasticity
bp_test_log <- car::ncvTest(model_log)
cat("Breusch-Pagan (homoscedasticity):", round(bp_test_log$p, 4),
    ifelse(bp_test_log$p > 0.05, "✓ Homoscedastic", "✗ Heteroscedastic"), "\n\n")

# Cook's Distance for log model
cooks_d_log <- cooks.distance(model_log)
influential_log <- which(cooks_d_log > 4/nrow(counts_log))
cat("Influential points (log scale):", length(influential_log), "/", nrow(counts_log), "\n\n")

# ==============================================================================
# 4. ROBUST REGRESSION (for comparison)
# ==============================================================================

cat("=== ROBUST REGRESSION (M-estimator) ===\n")

# Robust regression - less sensitive to outliers
model_robust <- MASS::rlm(FAIRD2 ~ FAIRD1, data = counts_data)

# Calculate pseudo-R² for robust regression
residuals_robust <- residuals(model_robust)
ss_res_robust <- sum(residuals_robust^2)
ss_tot <- sum((counts_data$FAIRD2 - mean(counts_data$FAIRD2))^2)
pseudo_r2_robust <- 1 - ss_res_robust/ss_tot

cat("Pseudo-R² (robust):", round(pseudo_r2_robust, 3), "\n\n")

# ==============================================================================
# 5. COEFFICIENT COMPARISON
# ==============================================================================

cat("=== COEFFICIENT COMPARISON ===\n")
cat("\nOriginal scale:\n")
print(coef(model_original))
cat("\nLog scale:\n")
print(coef(model_log))
cat("\nRobust (original scale):\n")
print(coef(model_robust))

# ==============================================================================
# 6. COMPARISON TABLE
# ==============================================================================

cat("\n\n=== SUMMARY COMPARISON ===\n\n")

comparison <- tibble::tibble(
  Method = c("Original scale", "Log-transformed", "Robust regression"),
  R_squared = c(
    summary(model_original)$r.squared,
    summary(model_log)$r.squared,
    pseudo_r2_robust
  ),
  Shapiro_p = c(
    shapiro_orig$p.value,
    shapiro_log$p.value,
    NA
  ),
  BP_p = c(
    bp_test$p,
    bp_test_log$p,
    NA
  ),
  N_influential = c(
    length(influential),
    length(influential_log),
    NA
  ),
  Recommendation = c(
    ifelse(bp_test$p < 0.05 | length(influential) > 2, "⚠️ Problematic", "✓ OK"),
    ifelse(shapiro_log$p.value > 0.05 & bp_test_log$p > 0.05, "✓ OK", "⚠️ Check"),
    "✓ Robust to outliers"
  )
) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 3)))

print(comparison)

# ==============================================================================
# 7. DIAGNOSTIC PLOTS
# ==============================================================================

cat("\n\n=== CREATING DIAGNOSTIC PLOTS ===\n")

# Save plots as PDF
pdf("regression_diagnostics_plots.pdf", width = 12, height = 10)

# Layout for 6 plots
par(mfrow = c(3, 2))

# --- Original scale plots ---
plot(counts_data$FAIRD1, counts_data$FAIRD2, 
     main = "Original Scale", 
     xlab = "FAIRD1 Abundance", 
     ylab = "FAIRD2 Abundance",
     pch = 19, col = ifelse(counts_data$influential, "red", "black"))
abline(model_original, col = "blue", lwd = 2)
legend("topleft", legend = c("Normal", "Influential"), 
       pch = 19, col = c("black", "red"), bty = "n")

plot(model_original, which = 1, main = "Original: Residuals vs Fitted")
plot(model_original, which = 2, main = "Original: Q-Q Plot")

# --- Log scale plots ---
plot(counts_log$FAIRD1_log, counts_log$FAIRD2_log,
     main = "Log-Transformed Scale",
     xlab = "log10(FAIRD1 + 1)",
     ylab = "log10(FAIRD2 + 1)",
     pch = 19, col = "darkgreen")
abline(model_log, col = "blue", lwd = 2)

plot(model_log, which = 1, main = "Log: Residuals vs Fitted")
plot(model_log, which = 2, main = "Log: Q-Q Plot")

dev.off()

cat("Diagnostic plots saved to: regression_diagnostics_plots.pdf\n")

# ==============================================================================
# 8. FINAL RECOMMENDATION
# ==============================================================================

cat("\n\n=== FINAL RECOMMENDATION ===\n\n")

if(bp_test$p < 0.05 | length(influential) > 2) {
  cat("⚠️  ORIGINAL SCALE VIOLATES ASSUMPTIONS:\n")
  cat("   - Heteroscedasticity (p =", round(bp_test$p, 4), ")\n")
  cat("   - Influential points:", length(influential), "\n\n")
  
  if(bp_test_log$p > 0.05 & shapiro_log$p.value > 0.05) {
    cat("✅ RECOMMENDATION: Use LOG-TRANSFORMED model\n")
    cat("   - Report R² =", round(summary(model_log)$r.squared, 3), "\n")
    cat("   - Meets homoscedasticity assumption\n")
    cat("   - Residuals are normally distributed\n")
  } else {
    cat("⚠️  RECOMMENDATION: Use ROBUST REGRESSION\n")
    cat("   - Pseudo-R² =", round(pseudo_r2_robust, 3), "\n")
    cat("   - Less sensitive to outliers\n")
    cat("   - Alternative: Report both log-transformed and robust results\n")
  }
} else {
  cat("✅ ORIGINAL SCALE MEETS ASSUMPTIONS\n")
  cat("   - Use R² =", round(summary(model_original)$r.squared, 3), "\n")
  cat("   - No transformation needed\n")
}

cat("\n✅ Diagnostics complete!\n")
cat("Review 'regression_diagnostics_plots.pdf' for visual inspection\n")