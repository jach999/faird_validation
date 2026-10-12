# ==============================================================================
# NORMALITY CHECK: FAIRD-AMMOD VALIDATION (H2)
# ==============================================================================

library(tidyverse)
library(ggpubr)

# Load data
load("data_cleaned.RData")
etraps_daily <- data_list$etraps_daily
ammod_biomass <- data_list$ammod_biomass_for_stats

# Source functions
source("analysis_functions.R")

cat("\n", rep("=", 70), "\n", sep="")
cat("HYPOTHESIS 2: FAIRD-AMMOD VALIDATION - NORMALITY CHECKS\n")
cat(rep("=", 70), "\n\n", sep="")

# ==============================================================================
# SITE 1 (Maize-Edge)
# ==============================================================================
cat("--- SITE 1 (Maize-Edge) ---\n")

site1_result <- analyze_faird_ammod_correlation(
  etraps_daily = etraps_daily,
  ammod_biomass = ammod_biomass,
  site_name = "Site1"
)

if(!is.null(site1_result)) {
  cat(sprintf("Sample size: n = %d days\n", site1_result$n_days))
  cat(sprintf("Pearson r = %.3f, R² = %.3f, p = %.4f\n", 
              site1_result$pearson_r, 
              site1_result$r_squared,
              site1_result$pearson_p))
  cat(sprintf("Shapiro-Wilk p = %.4f → Residuals normal? %s\n", 
              site1_result$shapiro_p,
              ifelse(site1_result$shapiro_p > 0.05, "YES ✓", "NO ✗")))
  cat(sprintf("Spearman ρ = %.3f, p = %.4f\n\n", 
              site1_result$spearman_rho,
              site1_result$spearman_p))
}

# ==============================================================================
# SITE 2 (Maize-Core) - KEY VALIDATION SITE
# ==============================================================================
cat("--- SITE 2 (Maize-Core) ⭐ KEY SITE ---\n")

site2_result <- analyze_faird_ammod_correlation(
  etraps_daily = etraps_daily,
  ammod_biomass = ammod_biomass,
  site_name = "Site2"
)

if(!is.null(site2_result)) {
  cat(sprintf("Sample size: n = %d days\n", site2_result$n_days))
  cat(sprintf("Pearson r = %.3f, R² = %.3f, p = %.4f\n", 
              site2_result$pearson_r, 
              site2_result$r_squared,
              site2_result$pearson_p))
  cat(sprintf("Shapiro-Wilk p = %.4f → Residuals normal? %s\n", 
              site2_result$shapiro_p,
              ifelse(site2_result$shapiro_p > 0.05, "YES ✓", "NO ✗")))
  cat(sprintf("Spearman ρ = %.3f, p = %.4f\n\n", 
              site2_result$spearman_rho,
              site2_result$spearman_p))
}

# ==============================================================================
# SITE 3 (Meadow) - PROBLEMATIC SITE
# ==============================================================================
cat("--- SITE 3 (Meadow) ⚠️ PROBLEMATIC ---\n")

site3_result <- analyze_faird_ammod_correlation(
  etraps_daily = etraps_daily,
  ammod_biomass = ammod_biomass,
  site_name = "Site3"
)

if(!is.null(site3_result)) {
  cat(sprintf("Sample size: n = %d days\n", site3_result$n_days))
  cat(sprintf("Pearson r = %.3f, R² = %.3f, p = %.4f\n", 
              site3_result$pearson_r, 
              site3_result$r_squared,
              site3_result$pearson_p))
  cat(sprintf("Shapiro-Wilk p = %.4f → Residuals normal? %s\n", 
              site3_result$shapiro_p,
              ifelse(site3_result$shapiro_p > 0.05, "YES ✓", "NO ✗")))
  cat(sprintf("Spearman ρ = %.3f, p = %.4f\n\n", 
              site3_result$spearman_rho,
              site3_result$spearman_p))
}

# ==============================================================================
# SUMMARY TABLE: FAIRD-AMMOD VALIDATION
# ==============================================================================
cat("\n", rep("=", 70), "\n", sep="")
cat("SUMMARY TABLE: FAIRD-AMMOD NORMALITY & CORRELATION\n")
cat(rep("=", 70), "\n\n", sep="")

validation_table <- tibble(
  Site = c("Site1 (Maize-Edge)", "Site2 (Maize-Core)", "Site3 (Meadow)"),
  `n (days)` = c(site1_result$n_days, 
                 site2_result$n_days, 
                 site3_result$n_days),
  `Pearson r` = round(c(site1_result$pearson_r, 
                        site2_result$pearson_r, 
                        site3_result$pearson_r), 3),
  `R²` = round(c(site1_result$r_squared, 
                 site2_result$r_squared, 
                 site3_result$r_squared), 3),
  `Pearson p` = round(c(site1_result$pearson_p, 
                        site2_result$pearson_p, 
                        site3_result$pearson_p), 4),
  `Shapiro-Wilk p` = round(c(site1_result$shapiro_p, 
                             site2_result$shapiro_p, 
                             site3_result$shapiro_p), 4),
  `Residuals Normal?` = ifelse(c(site1_result$shapiro_p, 
                                 site2_result$shapiro_p, 
                                 site3_result$shapiro_p) > 0.05, 
                               "Yes", "No"),
  `Spearman ρ` = round(c(site1_result$spearman_rho, 
                         site2_result$spearman_rho, 
                         site3_result$spearman_rho), 3),
  `Spearman p` = round(c(site1_result$spearman_p, 
                         site2_result$spearman_p, 
                         site3_result$spearman_p), 4),
  Validation = case_when(
    c(site1_result$r_squared, 
      site2_result$r_squared, 
      site3_result$r_squared) >= 0.50 ~ "Strong ✓",
    c(site1_result$r_squared, 
      site2_result$r_squared, 
      site3_result$r_squared) >= 0.30 ~ "Moderate",
    c(site1_result$r_squared, 
      site2_result$r_squared, 
      site3_result$r_squared) >= 0.10 ~ "Weak",
    TRUE ~ "None ✗"
  )
)

print(validation_table)

# Save results
write_csv(validation_table, "Table_S2_FAIRD_AMMOD_Validation.csv")
cat("\n✅ Results saved to: Table_S2_FAIRD_AMMOD_Validation.csv\n")

# ==============================================================================
# COMBINED SUMMARY: H1 + H2
# ==============================================================================
cat("\n", rep("=", 70), "\n", sep="")
cat("COMPLETE SUMMARY: ALL NORMALITY TESTS (H1 + H2)\n")
cat(rep("=", 70), "\n\n", sep="")

# Load H1 results
h1_table <- read_csv("Table_S1_Consistency_Statistics.csv", show_col_types = FALSE)

# Create combined summary
combined_summary <- tibble(
  Hypothesis = c(rep("H1: Consistency", 3), rep("H2: Validation", 3)),
  Comparison = c(
    "FAIRD1 vs FAIRD2",
    "FAIRD3 vs FAIRD4", 
    "ID1 vs ID2",
    "FAIRD vs AMMOD (Site1)",
    "FAIRD vs AMMOD (Site2)",
    "FAIRD vs AMMOD (Site3)"
  ),
  `Shapiro p (Log)` = c(
    0.657, 0.695, 0.321,
    site1_result$shapiro_p,
    site2_result$shapiro_p,
    site3_result$shapiro_p
  ),
  `Normal?` = ifelse(`Shapiro p (Log)` > 0.05, "Yes ✓", "No ✗"),
  `Pearson R²` = c(
    0.524, 0.643, 0.196,
    site1_result$r_squared,
    site2_result$r_squared,
    site3_result$r_squared
  ),
  `Conclusion` = case_when(
    `Pearson R²` >= 0.60 ~ "Excellent",
    `Pearson R²` >= 0.50 ~ "Good",
    `Pearson R²` >= 0.30 ~ "Moderate",
    `Pearson R²` >= 0.10 ~ "Weak",
    TRUE ~ "Poor"
  )
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(combined_summary)

write_csv(combined_summary, "Summary_All_Normality_Tests.csv")
cat("\n✅ Complete summary saved to: Summary_All_Normality_Tests.csv\n")

