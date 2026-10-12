# Create comprehensive comparison table
normality_comparison <- tibble(
  Comparison = rep(c("FAIRD1 vs FAIRD2", "FAIRD3 vs FAIRD4", "ID1 vs ID2"), 2),
  Transformation = rep(c("Original", "Log₁₀(x+1)"), each = 3),
  Shapiro_p_original = c(
    faird_maize_original$residuals_normality_p,
    faird_meadow_original$residuals_normality_p,
    id_maize_original$residuals_normality_p,
    NA, NA, NA
  ),
  Shapiro_p_log = c(
    NA, NA, NA,
    faird_maize_log$residuals_normality_p,
    faird_meadow_log$residuals_normality_p,
    id_maize_log$residuals_normality_p
  ),
  Normal = c(
    rep(FALSE, 3),  # Original all FALSE
    rep(TRUE, 3)    # Log all TRUE
  ),
  Pearson_r_original = c(
    faird_maize_original$pearson_r,
    faird_meadow_original$pearson_r,
    id_maize_original$pearson_r,
    NA, NA, NA
  ),
  Pearson_r_log = c(
    NA, NA, NA,
    faird_maize_log$pearson_r,
    faird_meadow_log$pearson_r,
    id_maize_log$pearson_r
  ),
  Spearman_rho = c(
    faird_maize_original$spearman_rho,
    faird_meadow_original$spearman_rho,
    id_maize_original$spearman_rho,
    faird_maize_log$spearman_rho,
    faird_meadow_log$spearman_rho,
    id_maize_log$spearman_rho
  )
)

# Clean up for paper
normality_table <- tibble(
  Comparison = c("FAIRD1 vs FAIRD2", "FAIRD3 vs FAIRD4", "ID1 vs ID2"),
  `Shapiro-Wilk p (Original)` = round(c(
    faird_maize_original$residuals_normality_p,
    faird_meadow_original$residuals_normality_p,
    id_maize_original$residuals_normality_p
  ), 4),
  `Normal (Original)` = c("No", "No", "No"),
  `Shapiro-Wilk p (Log)` = round(c(
    faird_maize_log$residuals_normality_p,
    faird_meadow_log$residuals_normality_p,
    id_maize_log$residuals_normality_p
  ), 4),
  `Normal (Log)` = c("Yes", "Yes", "Yes"),
  `Pearson r (Log)` = round(c(
    faird_maize_log$pearson_r,
    faird_meadow_log$pearson_r,
    id_maize_log$pearson_r
  ), 3),
  `Spearman ρ (Original)` = round(c(
    faird_maize_original$spearman_rho,
    faird_meadow_original$spearman_rho,
    id_maize_original$spearman_rho
  ), 3)
)

print(normality_table)
write_csv(normality_table, "Table_S1_Normality_Tests.csv")

# Añade R² a tu tabla
normality_table_complete <- tibble(
  Comparison = c("FAIRD1 vs FAIRD2", "FAIRD3 vs FAIRD4", "ID1 vs ID2"),
  `Shapiro-Wilk p (Original)` = c(0.0449, 0.0006, 0.0055),
  `Normal (Original)` = c("No", "No", "No"),
  `Shapiro-Wilk p (Log)` = c(0.657, 0.695, 0.321),
  `Normal (Log)` = c("Yes", "Yes", "Yes"),
  `Pearson r (Log)` = c(0.724, 0.802, 0.443),
  `R² (Log)` = c(0.724^2, 0.802^2, 0.443^2),  # R² = r²
  `Spearman ρ (Original)` = c(
    round(faird_maize_original$spearman_rho, 3),
    round(faird_meadow_original$spearman_rho, 3),
    round(id_maize_original$spearman_rho, 3)
  )
) %>%
  mutate(
    `R² (Log)` = round(`R² (Log)`, 3),
    Interpretation = case_when(
      `R² (Log)` >= 0.70 ~ "Excellent consistency",
      `R² (Log)` >= 0.50 ~ "Good consistency",
      `R² (Log)` >= 0.30 ~ "Moderate consistency",
      TRUE ~ "Poor consistency"
    )
  )

print(normality_table_complete)

write_csv(normality_table_complete, "Table_S1_Consistency_Statistics.csv")
cat("\n✅ Complete summary saved to: Table_S1_Consistency_Statistics.csv\n")

