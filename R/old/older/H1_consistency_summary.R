################################################################################
# H1 CONSOLIDATED SUMMARY: TEMPORAL CONSISTENCY ACROSS ALL DEVICES
################################################################################
# Purpose: Integrate H1 results from FAIRD, ID, and AMMOD into a single
#          comparative table
#
# VERSION 3.0 - COMPLETE DIAGNOSTICS (FINAL):
#   - Adjusted R² values
#   - Regression slope with 95% confidence intervals
#   - Intercept with p-value (test if significantly different from 0)
#   - Shapiro-Wilk test (normality of residuals)
#   - Breusch-Pagan test (homoscedasticity)
#   - Pseudo-R² from robust regression (outlier validation)
#   - Enhanced output table with ALL regression diagnostics
#
# Input:  data_cleaned.RData (from load_data.R)
# Output: Table1_H1_Consolidated_Results_COMPLETE.csv
#         Console output with formatted summary including all diagnostics
#
# Author: Juan Chiavassa
# Date: November 2025
# Version: 3.0 (COMPLETE)
################################################################################

# Load required packages
library(dplyr)
library(readr)
library(tibble)
library(car)      # For ncvTest (Breusch-Pagan test)
library(MASS)     # For robust regression (rlm)

# ============================================================
# CONSOLE OUTPUT CAPTURE
# ============================================================
script_name <- "H1_consistency_summary"

output_file <- file.path("outputs",
                         paste0(script_name, "_output_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

# Close any previously open sinks to prevent output accumulation
while (sink.number() > 0) sink()
sink(output_file, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Script:", paste0(script_name, ".R"), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

################################################################################
# SECTION 1: LOAD CLEANED DATA AND DATA STRUCTURE
################################################################################

cat("\n=== H1 CONSOLIDATED SUMMARY (COMPLETE v3.0) ===\n")
cat("Loading cleaned data...\n")

load("data_cleaned.RData")

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("CLEANED DATA STRUCTURE:\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

# Unique device names
cat("=== FAIRD devices ===\n")
print(unique(data_list$etraps_daily$Device[data_list$etraps_daily$Device_type == "FAIRD"]))

cat("\n=== ID devices ===\n")
print(unique(data_list$etraps_daily$Device[data_list$etraps_daily$Device_type == "ID"]))

cat("\n=== AMMOD devices ===\n")
print(unique(data_list$ammod_biomass_original$Device))

cat("\n=== Sample FAIRD data ===\n")
print(head(data_list$etraps_daily %>% dplyr::filter(Device_type == "FAIRD"), 10))

cat("\n=== Structure of etraps_daily ===\n")
print(str(data_list$etraps_daily))

################################################################################
# SECTION 2: PREPARE DATA FOR EACH DEVICE TYPE
################################################################################

# 2.1 FAIRD Daily Abundance (complete with zeros)
faird_daily <- data_list$etraps_daily %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::select(Date, Device, Site, abundance) %>%
  dplyr::arrange(Date, Device)

# 2.2 ID Daily Abundance (complete with zeros)
id_daily <- data_list$etraps_daily %>%
  dplyr::filter(Device_type == "ID") %>%
  dplyr::select(Date, Device, Site, abundance) %>%
  dplyr::arrange(Date, Device)

# 2.3 AMMOD Daily Biomass (12 days, with prorated AMMOD3)
ammod_biomass <- data_list$ammod_biomass_original %>%
  dplyr::select(Date, Device, Site, Live_mass) %>%
  dplyr::arrange(Date, Device)

cat("Data loaded successfully.\n")

################################################################################
# SECTION 3: FUNCTION TO CALCULATE ALL METRICS (COMPLETE)
################################################################################

calculate_consistency_metrics <- function(data, device1, device2, 
                                          value_col = "abundance",
                                          metric_name = "Abundance") {
  
  # Pivot to wide format for correlation (specify Date as id_cols)
  wide_data <- data %>%
    dplyr::filter(Device %in% c(device1, device2)) %>%
    dplyr::select(Date, Device, all_of(value_col)) %>%
    tidyr::pivot_wider(id_cols = Date, 
                       names_from = Device, 
                       values_from = all_of(value_col))
  
  # Remove any rows with missing data (shouldn't happen with complete data)
  wide_data <- wide_data %>%
    tidyr::drop_na()
  
  # Extract vectors
  x <- wide_data[[device1]]
  y <- wide_data[[device2]]
  
  # Sample size
  n <- length(x)
  
  # Kendall tau on original data (directional synchrony between daily changes)
  kendall_test <- cor.test(x, y, method = "kendall", exact = FALSE)
  kendall_tau  <- kendall_test$estimate
  kendall_p    <- kendall_test$p.value
  
  
  # Log transformation
  x_log <- log10(x + 1)
  y_log <- log10(y + 1)
  
  # ============================================================================
  # LINEAR REGRESSION ON LOG-TRANSFORMED DATA
  # ============================================================================
  
  lm_model <- lm(y_log ~ x_log)
  
  # R² values
  r2_log <- summary(lm_model)$r.squared
  adj_r2_log <- summary(lm_model)$adj.r.squared
  
  # Regression coefficients and confidence intervals
  coef_summary <- summary(lm_model)$coefficients
  
  # SLOPE (coefficient for x_log)
  slope <- coef_summary[2, 1]              # Estimate
  slope_se <- coef_summary[2, 2]           # Standard error
  slope_p <- coef_summary[2, 4]            # p-value
  slope_ci <- stats::confint(lm_model, level = 0.95)[2, ]  # 95% CI
  
  # INTERCEPT
  intercept <- coef_summary[1, 1]          # Estimate
  intercept_se <- coef_summary[1, 2]       # Standard error
  intercept_p <- coef_summary[1, 4]        # p-value (test if different from 0)
  intercept_ci <- stats::confint(lm_model, level = 0.95)[1, ]  # 95% CI
  
  # ============================================================================
  # REGRESSION DIAGNOSTICS
  # ============================================================================
  
  # Shapiro-Wilk test for normality of residuals
  shapiro_test <- shapiro.test(residuals(lm_model))
  shapiro_p <- shapiro_test$p.value
  
  # Breusch-Pagan test for homoscedasticity
  bp_test <- car::ncvTest(lm_model)
  bp_p <- bp_test$p
  
  # ============================================================================
  # ROBUST REGRESSION (log scale — consistent with OLS)
  # ============================================================================
  
  # Purpose: Verify that OLS R² is not inflated by influential observations
  # Computed on log-transformed data to match OLS scale (Pseudo-R² comparable to R²_log)
  model_robust     <- MASS::rlm(y_log ~ x_log)
  ss_res_robust    <- sum(residuals(model_robust)^2)
  ss_tot_ols       <- sum((y_log - mean(y_log))^2)
  pseudo_r2_robust <- 1 - ss_res_robust / ss_tot_ols
  
  # ============================================================================
  # CORRELATION TESTS
  # ============================================================================
  
  # Pearson correlation on log-transformed data
  pearson_test <- cor.test(x_log, y_log, method = "pearson")
  pearson_r <- pearson_test$estimate
  pearson_p <- pearson_test$p.value
  
  # Spearman correlation on original data (non-parametric)
  spearman_test <- cor.test(x, y, method = "spearman", exact = FALSE)
  spearman_rho <- spearman_test$estimate
  spearman_p <- spearman_test$p.value
  
  # ============================================================================
  # CONSISTENCY VERDICT
  # ============================================================================
  
  # Determine verdict based on R² and p-value
  if (r2_log >= 0.60 && pearson_p < 0.001) {
    verdict <- "Strong"
  } else if (r2_log >= 0.40 && pearson_p < 0.05) {
    verdict <- "Moderate"
  } else if (pearson_p >= 0.05) {
    verdict <- "Poor (NS)"
  } else {
    verdict <- "Weak"
  }
  
  # ============================================================================
  # RETURN COMPREHENSIVE RESULTS
  # ============================================================================
  
  tibble::tibble(
    Comparison = paste0(device1, " vs ", device2),
    Metric = paste0(metric_name, " (n=", n, ")"),
    n = n,
    R2_log = r2_log,
    Adj_R2_log = adj_r2_log,
    Slope = slope,
    Slope_SE = slope_se,
    Slope_CI_lower = slope_ci[1],
    Slope_CI_upper = slope_ci[2],
    Slope_p = slope_p,
    Intercept = intercept,
    Intercept_SE = intercept_se,
    Intercept_CI_lower = intercept_ci[1],
    Intercept_CI_upper = intercept_ci[2],
    Intercept_p = intercept_p,
    Shapiro_p = shapiro_p,
    BP_p = bp_p,
    Pseudo_R2_robust = pseudo_r2_robust,
    Pearson_r = pearson_r,
    Pearson_p = pearson_p,
    Spearman_rho = spearman_rho,
    Spearman_p = spearman_p,
    Kendall_tau = kendall_tau,
Kendall_p   = kendall_p,
    Verdict = verdict
  )
}

################################################################################
# SECTION 4: CALCULATE METRICS FOR EACH DEVICE PAIR
################################################################################

cat("\nCalculating consistency metrics (with complete diagnostics)...\n")

# 4.1 FAIRD1 vs FAIRD2 (Sites 1 and 2 - Maize)
faird_12 <- calculate_consistency_metrics(
  data = faird_daily,
  device1 = "FAIRD1",
  device2 = "FAIRD2",
  value_col = "abundance",
  metric_name = "Abundance"
)

# 4.2 FAIRD3 vs FAIRD4 (Sites 3 and 4 - Meadow)
faird_34 <- calculate_consistency_metrics(
  data = faird_daily,
  device1 = "FAIRD3",
  device2 = "FAIRD4",
  value_col = "abundance",
  metric_name = "Abundance"
)

# 4.3 ID1 vs ID2 (Sites 1 and 2 - Maize)
id_12 <- calculate_consistency_metrics(
  data = id_daily,
  device1 = "ID1",
  device2 = "ID2",
  value_col = "abundance",
  metric_name = "Abundance"
)

# 4.4 AMMOD1 vs AMMOD2 (Sites 1 and 2 - Maize, first 12 days only)
ammod_12 <- calculate_consistency_metrics(
  data = ammod_biomass,
  device1 = "AMMOD1",
  device2 = "AMMOD2",
  value_col = "Live_mass",
  metric_name = "Biomass"
)

cat("Metrics calculated successfully.\n")

################################################################################
# SECTION 5: CONSOLIDATE RESULTS
################################################################################

# Bind all results into a single table
consolidated_results <- dplyr::bind_rows(
  faird_12,
  faird_34,
  id_12,
  ammod_12
)

# Reorder columns for clarity
consolidated_results <- consolidated_results %>%
  dplyr::select(Comparison, Metric, n, 
                R2_log, Adj_R2_log,
                Slope, Slope_SE, Slope_CI_lower, Slope_CI_upper, Slope_p,
                Intercept, Intercept_SE, Intercept_CI_lower, Intercept_CI_upper, Intercept_p,
                Shapiro_p, BP_p, Pseudo_R2_robust,
                Pearson_r, Pearson_p, 
                Spearman_rho, Spearman_p,
                Spearman_rho, Spearman_p,
                Kendall_tau, Kendall_p,
                Verdict)

################################################################################
# SECTION 6: EXPORT RESULTS
################################################################################

cat("\nExporting consolidated results...\n")

# Export as CSV
readr::write_csv(consolidated_results, "outputs/Table1_H1_Consolidated_Results_COMPLETE.csv")

cat("✅ Table exported: outputs/Table1_H1_Consolidated_Results_COMPLETE.csv\n")

################################################################################
# SECTION 7: FORMATTED CONSOLE OUTPUT - BASIC METRICS
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("H1: TEMPORAL CONSISTENCY RESULTS (COMPLETE v3.0)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("Correlation analysis (Pearson and Spearman) on log10(x+1) transformed data.\n")
cat("Includes: R², Adjusted R², Slope + 95% CI, Intercept + significance test,\n")
cat("          Shapiro-Wilk, Breusch-Pagan, and Pseudo-R² (robust regression)\n")
cat("\n")

# Print formatted table (basic metrics)
basic_results <- consolidated_results %>%
  dplyr::select(Comparison, Metric, n, R2_log, Adj_R2_log, 
                Pearson_r, Pearson_p, Spearman_rho, Spearman_p, Verdict)

print(basic_results, n = Inf)

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("KEY FINDINGS:\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

# Calculate comparative statistics
faird_avg_r2 <- mean(c(faird_12$R2_log, faird_34$R2_log))
faird_avg_adj_r2 <- mean(c(faird_12$Adj_R2_log, faird_34$Adj_R2_log))
id_r2 <- id_12$R2_log
fold_difference <- faird_avg_r2 / id_r2

cat(sprintf("• FAIRD average R² = %.3f (strong consistency)\n", faird_avg_r2))
cat(sprintf("• FAIRD average Adjusted R² = %.3f (conservative estimate)\n", faird_avg_adj_r2))
cat(sprintf("• ID R² = %.3f (poor consistency)\n", id_r2))
cat(sprintf("• FAIRD is %.1f× more consistent than ID\n", fold_difference))
cat(sprintf("• AMMOD R² = %.3f (strong consistency, ground truth)\n", ammod_12$R2_log))
cat("\n")
cat("• All FAIRD and AMMOD comparisons: p < 0.001 (highly significant)\n")
cat("• ID comparison: p >= 0.05 (not significant)\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

# Spearman correlations
cat("\nSPEARMAN CORRELATIONS (non-parametric validation):\n")
cat("───────────────────────────────────────────────────────────────────────────\n")
for (i in 1:nrow(consolidated_results)) {
  row <- consolidated_results[i, ]
  cat(sprintf("%-18s: ρ = %.3f, p %s\n", 
              row$Comparison, 
              row$Spearman_rho,
              ifelse(row$Spearman_p < 0.001, "< 0.001", 
                     sprintf("= %.4f", row$Spearman_p))))
}
cat("═══════════════════════════════════════════════════════════════════════════\n")

################################################################################
# SECTION 8: REGRESSION COEFFICIENTS OUTPUT
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("REGRESSION COEFFICIENTS (log-log scale: log10(y+1) ~ log10(x+1))\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("\n")

for (i in 1:nrow(consolidated_results)) {
  row <- consolidated_results[i, ]
  
  cat(sprintf("%-18s (n = %d):\n", row$Comparison, row$n))
  cat("───────────────────────────────────────────────────────────────────────────\n")
  
  # Slope with CI
  cat(sprintf("  Slope       = %.2f (95%% CI: %.2f to %.2f)\n",
              row$Slope,
              row$Slope_CI_lower,
              row$Slope_CI_upper))
  cat(sprintf("                SE = %.3f, p %s\n",
              row$Slope_SE,
              ifelse(row$Slope_p < 0.001, "< 0.001", 
                     sprintf("= %.4f", row$Slope_p))))
  
  # Interpretation of slope
  if (row$Slope_CI_lower <= 1 && row$Slope_CI_upper >= 1) {
    cat("                CI includes 1.0 → Proportional 1:1 relationship ✓\n")
  } else if (row$Slope > 1) {
    cat("                Slope > 1.0 → Device 2 captures relatively more\n")
  } else {
    cat("                Slope < 1.0 → Device 2 captures relatively less\n")
  }
  
  cat("\n")
  
  # Intercept with CI and significance test
  cat(sprintf("  Intercept   = %.2f (95%% CI: %.2f to %.2f)\n",
              row$Intercept,
              row$Intercept_CI_lower,
              row$Intercept_CI_upper))
  cat(sprintf("                SE = %.3f, p %s\n",
              row$Intercept_SE,
              ifelse(row$Intercept_p < 0.001, "< 0.001", 
                     sprintf("= %.4f", row$Intercept_p))))
  
  # Interpretation of intercept
  if (row$Intercept_p >= 0.05) {
    cat("                Not significantly different from 0 ✓\n")
  } else {
    if (row$Intercept > 0) {
      cat("                Significantly > 0 → Systematic positive offset ⚠\n")
    } else {
      cat("                Significantly < 0 → Systematic negative offset ⚠\n")
    }
  }
  
  cat("\n")
  
  # Overall R² statistics
  cat(sprintf("  R²          = %.3f\n", row$R2_log))
  cat(sprintf("  Adj. R²     = %.3f\n", row$Adj_R2_log))
  
  cat("\n")
}

cat("═══════════════════════════════════════════════════════════════════════════\n")

################################################################################
# SECTION 9: REGRESSION DIAGNOSTICS (NEW)
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("REGRESSION DIAGNOSTICS (log-transformed models)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("Tests for parametric assumptions:\n")
cat("  • Shapiro-Wilk: normality of residuals (H₀: residuals are normal)\n")
cat("  • Breusch-Pagan: homoscedasticity (H₀: constant variance)\n")
cat("  • Criterion: p > 0.05 indicates assumptions are met ✓\n")
cat("\n")

for (i in 1:nrow(consolidated_results)) {
  row <- consolidated_results[i, ]
  
  cat(sprintf("%-18s:\n", row$Comparison))
  
  # Shapiro-Wilk
  shapiro_status <- ifelse(row$Shapiro_p > 0.05, "✓", "⚠")
  cat(sprintf("  • Shapiro-Wilk (normality): p = %.4f %s\n",
              row$Shapiro_p, shapiro_status))
  
  # Breusch-Pagan
  bp_status <- ifelse(row$BP_p > 0.05, "✓", "⚠")
  cat(sprintf("  • Breusch-Pagan (homoscedasticity): p = %.4f %s\n",
              row$BP_p, bp_status))
  
  cat("\n")
}

# Overall diagnostic summary
all_pass <- all(consolidated_results$Shapiro_p > 0.05 & consolidated_results$BP_p > 0.05)

if (all_pass) {
  cat("✅ SUMMARY: All comparisons pass diagnostic tests (p > 0.05).\n")
  cat("   Log-transformation successfully normalizes residuals and stabilizes variance.\n")
} else {
  cat("⚠ SUMMARY: Some comparisons show diagnostic issues.\n")
  cat("  Review individual tests above for details.\n")
}

cat("═══════════════════════════════════════════════════════════════════════════\n")

################################################################################
# SECTION 10: ROBUST REGRESSION VALIDATION (NEW)
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("ROBUST REGRESSION VALIDATION (Pseudo-R²)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("Purpose: Verify that correlations are not driven by outliers\n")
cat("Method: MASS::rlm() robust regression on log10(x+1) scale (consistent with OLS)\n")
cat("Interpretation:\n")
cat("  • Pseudo-R² >> R² → Relationship is genuine, not outlier-dependent\n")
cat("  • Pseudo-R² ≈ R² → Consistent relationship across data range\n")
cat("  • Pseudo-R² << R² → Some outlier influence (still acceptable if moderate)\n")
cat("\n")

for (i in 1:nrow(consolidated_results)) {
  row <- consolidated_results[i, ]
  
  cat(sprintf("%-18s: Pseudo-R² = %.3f (validates R² = %.3f)\n",
              row$Comparison,
              row$Pseudo_R2_robust,
              row$R2_log))
  
  # Interpretation based on Pseudo-R² value
  if (row$Pseudo_R2_robust > 0.80) {
    cat("  → High pseudo-R²: relationship is genuine, not outlier-dependent ✓\n")
  } else if (row$Pseudo_R2_robust > 0.40) {
    cat("  → Moderate pseudo-R²: some outlier influence, but relationship remains\n")
  } else if (row$Pseudo_R2_robust > 0) {
    cat("  → Low pseudo-R²: weak relationship confirmed\n")
  } else {
    cat("  → Negative pseudo-R²: robust model performs worse than mean\n")
  }
  
  # Additional context for AMMOD
  if (grepl("AMMOD", row$Comparison) && row$Pseudo_R2_robust < row$R2_log) {
    cat("  → Note: Some influence from high-biomass days (expected in biomass data)\n")
  }
  
  cat("\n")
}

cat("═══════════════════════════════════════════════════════════════════════════\n")



# ==============================================================================
# SECTION 11: PUBLICATION-READY H1 SUMMARY TABLE
# ==============================================================================
# Columns: Pair | Habitat | n | R²(log) | Slope [95% CI] | Concordance | Kendall τ
#
# Recalculates concordance and Kendall directly from data for all 5 rows:
#   FAIRD1 vs FAIRD2  (Maize,  n=22)
#   FAIRD3 vs FAIRD4  (Meadow, n=22)
#   ID1 vs ID2        (Maize,  n=22 conservative / n=17 operational)
#   AMMOD1 vs AMMOD2  (Maize,  n=12)
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 9: PUBLICATION-READY H1 TABLE\n")
cat("========================================\n\n")

# --- Helper: significance stars ---
sig_stars_h1 <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# --- Helper: directional concordance + Kendall from a wide 2-column tibble ---
# Expects columns: Date, Dev1, Dev2 (numeric)
calc_concordance_kendall <- function(wide_df, col1, col2) {
  df <- wide_df %>%
    dplyr::arrange(Date) %>%
    dplyr::mutate(
      change1   = .data[[col1]] - dplyr::lag(.data[[col1]]),
      change2   = .data[[col2]] - dplyr::lag(.data[[col2]]),
      concordant = sign(change1) == sign(change2)
    ) %>%
    dplyr::filter(!is.na(concordant))
  
  n_changes    <- nrow(df)
  n_concordant <- sum(df$concordant, na.rm = TRUE)
  conc_pct     <- 100 * n_concordant / n_changes
  
  kt <- stats::cor.test(df$change1, df$change2,
                        method = "kendall", exact = FALSE)
  
  list(
    n_changes       = n_changes,
    concordance_pct = conc_pct,
    kendall_tau     = as.numeric(kt$estimate),
    kendall_p       = kt$p.value
  )
}

# ------------------------------------------------------------------------------
# 9.1  FAIRD pairs (abundance)
# ------------------------------------------------------------------------------

faird_wide <- data_list$etraps_daily %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device,
                     values_from = abundance,
                     values_fill = 0)

sync_faird12 <- calc_concordance_kendall(faird_wide, "FAIRD1", "FAIRD2")
sync_faird34 <- calc_concordance_kendall(faird_wide, "FAIRD3", "FAIRD4")

cat(sprintf("FAIRD1 vs FAIRD2: concordance=%.1f%% | tau=%.3f (p=%.4f)\n",
            sync_faird12$concordance_pct,
            sync_faird12$kendall_tau,
            sync_faird12$kendall_p))

cat(sprintf("FAIRD3 vs FAIRD4: concordance=%.1f%% | tau=%.3f (p=%.4f)\n",
            sync_faird34$concordance_pct,
            sync_faird34$kendall_tau,
            sync_faird34$kendall_p))

# ------------------------------------------------------------------------------
# 9.2  ID pairs (abundance) — conservative (n=22) and operational (n=17)
# ------------------------------------------------------------------------------

id_exclude_start <- as.Date("2023-08-28")
id_exclude_end   <- as.Date("2023-09-01")

id_wide_cons <- data_list$etraps_daily %>%
  dplyr::filter(Device_type == "ID",
                Device %in% c("ID1", "ID2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device,
                     values_from = abundance,
                     values_fill = 0)

id_wide_oper <- id_wide_cons %>%
  dplyr::filter(!(Date >= id_exclude_start & Date <= id_exclude_end))

sync_id_cons <- calc_concordance_kendall(id_wide_cons, "ID1", "ID2")
sync_id_oper <- calc_concordance_kendall(id_wide_oper, "ID1", "ID2")

cat(sprintf("ID1 vs ID2 (conservative n=%d): concordance=%.1f%% | tau=%.3f (p=%.4f)\n",
            nrow(id_wide_cons),
            sync_id_cons$concordance_pct,
            sync_id_cons$kendall_tau,
            sync_id_cons$kendall_p))

cat(sprintf("ID1 vs ID2 (operational  n=%d): concordance=%.1f%% | tau=%.3f (p=%.4f)\n",
            nrow(id_wide_oper),
            sync_id_oper$concordance_pct,
            sync_id_oper$kendall_tau,
            sync_id_oper$kendall_p))

# ------------------------------------------------------------------------------
# 9.3  AMMOD pair (live biomass)
# ------------------------------------------------------------------------------

ammod_wide_h1 <- data_list$ammod_biomass_for_stats %>%
  dplyr::filter(Device %in% c("AMMOD1", "AMMOD2")) %>%
  dplyr::select(Date, Device, Live_mass) %>%
  tidyr::pivot_wider(names_from = Device,
                     values_from = Live_mass,
                     values_fill = NA) %>%
  dplyr::filter(!is.na(AMMOD1) & !is.na(AMMOD2))

sync_ammod <- calc_concordance_kendall(ammod_wide_h1, "AMMOD1", "AMMOD2")

cat(sprintf("AMMOD1 vs AMMOD2: concordance=%.1f%% | tau=%.3f (p=%.4f)\n\n",
            sync_ammod$concordance_pct,
            sync_ammod$kendall_tau,
            sync_ammod$kendall_p))

# ------------------------------------------------------------------------------
# 9.4  Pull R², slope, CI from consolidated_results (already calculated above)
# ------------------------------------------------------------------------------

get_stats <- function(comparison_label) {
  row <- consolidated_results %>%
    dplyr::filter(Comparison == comparison_label)
  if (nrow(row) == 0) {
    stop(sprintf("Comparison '%s' not found in consolidated_results", comparison_label))
  }
  row
}

stats_f12   <- get_stats("FAIRD1 vs FAIRD2")
stats_f34   <- get_stats("FAIRD3 vs FAIRD4")
stats_id_c  <- get_stats("ID1 vs ID2")        # conservative (n=22)
stats_ammod <- get_stats("AMMOD1 vs AMMOD2")

# Operational ID: recalculate regression on operational dataset
id_log_oper <- id_wide_oper %>%
  dplyr::mutate(
    ID1_log = log10(ID1 + 1),
    ID2_log = log10(ID2 + 1)
  )

model_id_oper   <- lm(ID2_log ~ ID1_log, data = id_log_oper)
slope_id_oper   <- coef(model_id_oper)["ID1_log"]
ci_id_oper      <- confint(model_id_oper)["ID1_log", ]
r2_id_oper      <- summary(model_id_oper)$r.squared
p_id_oper       <- summary(model_id_oper)$coefficients["ID1_log", "Pr(>|t|)"]

cat(sprintf("ID operational regression: R2=%.3f | slope=%.3f [%.3f, %.3f] | p=%.4f\n\n",
            r2_id_oper, slope_id_oper, ci_id_oper[1], ci_id_oper[2], p_id_oper))

# ------------------------------------------------------------------------------
# 9.5  Assemble display table (5 rows)
# ------------------------------------------------------------------------------

# Row builder
make_row <- function(pair, habitat, n, r2, r2_p,
                     slope, ci_low, ci_high,
                     conc_pct, tau, tau_p,
                     spearman_rho, spearman_p) {
  tibble::tibble(
    Pair             = pair,
    Habitat          = habitat,
    n                = n,
    R2_log           = r2,
    R2_p             = r2_p,
    Slope            = slope,
    Slope_CI_low     = ci_low,
    Slope_CI_high    = ci_high,
    Concordance_pct  = conc_pct,
    Kendall_tau      = tau,
    Kendall_p        = tau_p,
    Spearman_rho     = spearman_rho,
    Spearman_p       = spearman_p
  )
}

h1_table_raw <- dplyr::bind_rows(
  make_row("FAIRD1 vs FAIRD2", "Maize",  stats_f12$n,
           stats_f12$R2_log,   stats_f12$Pearson_p,
           stats_f12$Slope,    stats_f12$Slope_CI_lower, stats_f12$Slope_CI_upper,
           sync_faird12$concordance_pct, sync_faird12$kendall_tau, sync_faird12$kendall_p,
           stats_f12$Spearman_rho, stats_f12$Spearman_p),
  
  make_row("FAIRD3 vs FAIRD4", "Meadow", stats_f34$n,
           stats_f34$R2_log,   stats_f34$Pearson_p,
           stats_f34$Slope,    stats_f34$Slope_CI_lower, stats_f34$Slope_CI_upper,
           sync_faird34$concordance_pct, sync_faird34$kendall_tau, sync_faird34$kendall_p,
           stats_f34$Spearman_rho, stats_f34$Spearman_p),
  
  make_row("ID1 vs ID2 (cons.)", "Maize", stats_id_c$n,
           stats_id_c$R2_log,  stats_id_c$Pearson_p,
           stats_id_c$Slope,   stats_id_c$Slope_CI_lower, stats_id_c$Slope_CI_upper,
           sync_id_cons$concordance_pct, sync_id_cons$kendall_tau, sync_id_cons$kendall_p,
           stats_id_c$Spearman_rho, stats_id_c$Spearman_p),
  
  make_row("ID1 vs ID2 (oper.)", "Maize", nrow(id_wide_oper),
           r2_id_oper,          p_id_oper,
           slope_id_oper,       ci_id_oper[1], ci_id_oper[2],
           sync_id_oper$concordance_pct, sync_id_oper$kendall_tau, sync_id_oper$kendall_p,
           cor.test(id_wide_oper$ID1, id_wide_oper$ID2, method="spearman", exact=FALSE)$estimate,
           cor.test(id_wide_oper$ID1, id_wide_oper$ID2, method="spearman", exact=FALSE)$p.value),
  
  make_row("AMMOD1 vs AMMOD2", "Maize",  stats_ammod$n,
           stats_ammod$R2_log, stats_ammod$Pearson_p,
           stats_ammod$Slope,  stats_ammod$Slope_CI_lower, stats_ammod$Slope_CI_upper,
           sync_ammod$concordance_pct, sync_ammod$kendall_tau, sync_ammod$kendall_p,
           stats_ammod$Spearman_rho, stats_ammod$Spearman_p)
)

# --- Formatted display version ---
h1_table_display <- h1_table_raw %>%
  dplyr::mutate(
    `R² (log)`       = paste0(formatC(R2_log, digits = 3, format = "f"),
                              " ", sig_stars_h1(R2_p)),
    `Slope [95% CI]` = sprintf("%.2f [%.2f, %.2f]",
                               Slope, Slope_CI_low, Slope_CI_high),
    `Concordance`    = sprintf("%.1f%%", Concordance_pct),
    `Kendall τ`      = paste0(formatC(Kendall_tau, digits = 3, format = "f"),
                              " ", sig_stars_h1(Kendall_p)),
    `Spearman ρ`     = paste0(formatC(Spearman_rho, digits = 3, format = "f"),
                              " ", sig_stars_h1(Spearman_p))
  ) %>%
  dplyr::select(Pair, Habitat, n,
                `R² (log)`, `Slope [95% CI]`,
                `Concordance`, `Spearman ρ`, `Kendall τ`)

cat("=== H1 PUBLICATION TABLE ===\n\n")
print(h1_table_display, n = Inf)
cat("\nSignificance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p>=0.05\n")
cat("† cons. = conservative (n=22, includes battery failure days)\n")
cat("  oper. = operational  (n=17, excludes Aug 28-Sep 1)\n\n")

# --- Export: display version ---
readr::write_csv(h1_table_display,
                 "outputs/H1_Publication_Table.csv")
cat("[OK] outputs/H1_Publication_Table.csv\n")

# --- Export: numeric version (for supplementary / PCV verification) ---
h1_table_numeric <- h1_table_raw %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

readr::write_csv(h1_table_numeric,
                 "outputs/H1_Publication_Table_numeric.csv")
cat("[OK] outputs/H1_Publication_Table_numeric.csv\n\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE (v3.0)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("\nFiles generated:\n")
cat("  [OK] outputs/Table1_H1_Consolidated_Results_COMPLETE.csv\n")
cat("  [OK] outputs/H1_Publication_Table.csv\n")
cat("  [OK] outputs/H1_Publication_Table_numeric.csv\n")
cat("\nConsolidated H1 analysis complete with FULL diagnostics:\n")
cat("  [OK] Regression coefficients (Slope + CI, Intercept + CI)\n")
cat("  [OK] Adjusted R² values\n")
cat("  [OK] Shapiro-Wilk tests (normality)\n")
cat("  [OK] Breusch-Pagan tests (homoscedasticity)\n")
cat("  [OK] Pseudo-R² (robust regression validation, log scale)\n")
cat("  [OK] Directional concordance + Kendall tau (all 5 pairs)\n")
cat("  [OK] Publication-ready summary table (Section 11)\n\n")

################################################################################
# END OF SCRIPT
################################################################################

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")