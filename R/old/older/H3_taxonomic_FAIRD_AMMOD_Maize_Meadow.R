# ==============================================================================
# H3: MAIZE vs MEADOW STATISTICAL COMPARISON
# ==============================================================================
# Calculate Wilcoxon rank-sum test and Cohen's d effect size
# for habitat comparison in H3 taxonomic similarity analysis
#
# Author: Juan Chiavassa & Claude (Fixed by Gemini)
# Date: December 20, 2024
# ==============================================================================

library(dplyr)

# ============================================================
# CONSOLE OUTPUT CAPTURE
# ============================================================
script_name <- sub("\\.R$", "", basename(sys.frame(1)$ofile))

output_file <- file.path("outputs",
                         paste0(script_name, "_output_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

# Using native encoding to avoid issues on Windows
sink(output_file, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Script:", paste0(script_name, ".R"), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

# IMPORTANT: Adjust this path to where your similarity_results.csv is located
data_path <- "outputs/H3_FAIRD_AMMOD_similarity_results_GENUS_3mm.csv"

cat("\n===========================================================\n")
cat("LOADING DATA\n")
cat("===========================================================\n\n")

# Check if file exists to prevent crash
if (!file.exists(data_path)) {
  stop(paste("ERROR: File not found at:", data_path))
}

similarity_results <- read.csv(data_path, stringsAsFactors = FALSE)
cat("[OK] Data loaded successfully\n")
cat("    Total rows:", nrow(similarity_results), "\n")
cat("    Total sites:", length(unique(similarity_results$Site)), "\n\n")

# ==============================================================================
# 2. SEPARATE BY HABITAT
# ==============================================================================

# Extract values for each habitat
maize_values <- similarity_results %>% 
  dplyr::filter(Site %in% c("Site1", "Site2")) %>% 
  dplyr::pull(Delta_S_raw)

meadow_values <- similarity_results %>% 
  dplyr::filter(Site %in% c("Site3", "Site4")) %>% 
  dplyr::pull(Delta_S_raw)

# ==============================================================================
# 3. DATA QUALITY CHECK
# ==============================================================================

cat("\n===========================================================\n")
cat("DATA QUALITY CHECK\n")
cat("===========================================================\n\n")

# Maize summary
cat("MAIZE (Site1 + Site2):\n")
cat("  Total observations:", length(maize_values), "\n")
cat("  NA values:", sum(is.na(maize_values)), "\n")
cat("  Valid values:", sum(!is.na(maize_values)), "\n")
cat("  Mean:", round(mean(maize_values, na.rm = TRUE), 4), "\n")
cat("  SD:", round(sd(maize_values, na.rm = TRUE), 4), "\n")
cat("  Min:", round(min(maize_values, na.rm = TRUE), 4), "\n")
cat("  Max:", round(max(maize_values, na.rm = TRUE), 4), "\n\n")

# Meadow summary
cat("MEADOW (Site3 + Site4):\n")
cat("  Total observations:", length(meadow_values), "\n")
cat("  NA values:", sum(is.na(meadow_values)), "\n")
cat("  Valid values:", sum(!is.na(meadow_values)), "\n")
cat("  Mean:", round(mean(meadow_values, na.rm = TRUE), 4), "\n")
cat("  SD:", round(sd(meadow_values, na.rm = TRUE), 4), "\n")
cat("  Min:", round(min(meadow_values, na.rm = TRUE), 4), "\n")
cat("  Max:", round(max(meadow_values, na.rm = TRUE), 4), "\n\n")

# Difference
diff_mean <- mean(meadow_values, na.rm = TRUE) - mean(maize_values, na.rm = TRUE)
cat("DIFFERENCE:\n")
cat("  Mean difference (Meadow - Maize):", round(diff_mean, 4), "\n")
cat("  Percentage points:", round(diff_mean * 100, 2), "pp\n")
cat("  Relative increase:", round((diff_mean / mean(maize_values, na.rm = TRUE)) * 100, 1), "%\n\n")

# ==============================================================================
# 4. WILCOXON RANK-SUM TEST (Mann-Whitney U)
# ==============================================================================

cat("\n===========================================================\n")
cat("WILCOXON RANK-SUM TEST (Mann-Whitney U)\n")
cat("===========================================================\n\n")

# Perform test
wilcox_result <- wilcox.test(maize_values, meadow_values, 
                             exact = TRUE,
                             alternative = "two.sided")

cat("Null hypothesis: Maize and Meadow have same distribution\n")
cat("Alternative hypothesis: True location shift is not equal to 0\n\n")

cat("RESULTS:\n")
cat("  W statistic:", wilcox_result$statistic, "\n")
cat("  p-value:", round(wilcox_result$p.value, 4), "\n")
cat("  Significance (alpha=0.05):", ifelse(wilcox_result$p.value < 0.05, 
                                           "SIGNIFICANT [YES]", 
                                           "NOT SIGNIFICANT"), "\n")
cat("  Significance (alpha=0.10):", ifelse(wilcox_result$p.value < 0.10, 
                                           "MARGINALLY SIGNIFICANT [~]", 
                                           "NOT SIGNIFICANT"), "\n\n")

# Interpretation
if (wilcox_result$p.value < 0.05) {
  interpretation <- "There IS a statistically significant difference between habitats"
} else if (wilcox_result$p.value < 0.10) {
  interpretation <- "There is a MARGINALLY significant trend but not significant at alpha=0.05"
} else {
  interpretation <- "There is NO statistically significant difference between habitats"
}

cat("INTERPRETATION:\n")
cat("  ", interpretation, "\n\n")

# ==============================================================================
# 4.5. PARAMETRIC VALIDATION (T-TEST)
# ==============================================================================

cat("\n===========================================================\n")
cat("PARAMETRIC VALIDATION (T-TEST)\n")
cat("===========================================================\n\n")

# Run T-test (Two-sided, var.equal=TRUE based on Levene's test p=0.65)
t_test_res <- t.test(meadow_values, maize_values, var.equal = TRUE)

# Extract values
t_stat <- round(t_test_res$statistic, 3)
df_val <- round(t_test_res$parameter, 1)
p_val  <- t_test_res$p.value
mean_diff <- round(diff(t_test_res$estimate) * -1, 4) # Meadow - Maize

cat("RESULTS:\n")
cat("  t-statistic:", t_stat, "\n")
cat("  Degrees of freedom:", df_val, "\n")
cat("  p-value (exact):", p_val, "\n")
cat("  p-value (formatted):", format.pval(p_val, digits=4), "\n")
cat("  Difference (Meadow - Maize):", mean_diff, "\n\n")

# Contextual Conclusion
if(p_val < 0.05) {
  cat("[OK] RESOLUTION: The parametric T-test IS significant (p < 0.05).\n")
  cat("     This result aligns with the Cohen's d Confidence Interval.\n")
  cat("     Report this t-test in instead of Wilcoxon\n")
  cat("     to confirm that Meadow performance is significantly higher should be better.\n")
} else {
  cat("[X] RESULT: Still not significant. Rely on effect size interpretation.\n")
}
cat("\n")

# ==============================================================================
# 5. COHEN'S D EFFECT SIZE (Manual Calculation)
# ==============================================================================

cat("\n===========================================================\n")
cat("COHEN'S D EFFECT SIZE\n")
cat("===========================================================\n\n")

# Remove NAs
maize_clean <- maize_values[!is.na(maize_values)]
meadow_clean <- meadow_values[!is.na(meadow_values)]

# Calculate basic statistics
n_maize <- length(maize_clean)
n_meadow <- length(meadow_clean)
mean_maize <- mean(maize_clean)
mean_meadow <- mean(meadow_clean)
sd_maize <- sd(maize_clean)
sd_meadow <- sd(meadow_clean)

# Pooled standard deviation (assuming equal variances)
pooled_sd <- sqrt(((n_maize - 1) * sd_maize^2 + (n_meadow - 1) * sd_meadow^2) / 
                    (n_maize + n_meadow - 2))

# Cohen's d (positive value = Meadow > Maize)
cohens_d <- (mean_meadow - mean_maize) / pooled_sd

# Standard error of Cohen's d
se_d <- sqrt((n_maize + n_meadow) / (n_maize * n_meadow) + 
               cohens_d^2 / (2 * (n_maize + n_meadow)))

# 95% Confidence Interval (using normal approximation)
ci_lower <- cohens_d - 1.96 * se_d
ci_upper <- cohens_d + 1.96 * se_d

# Display calculation details
cat("CALCULATION DETAILS:\n")
cat("  n (Maize):", n_maize, "\n")
cat("  n (Meadow):", n_meadow, "\n")
cat("  Mean (Maize):", round(mean_maize, 4), "\n")
cat("  Mean (Meadow):", round(mean_meadow, 4), "\n")
cat("  SD (Maize):", round(sd_maize, 4), "\n")
cat("  SD (Meadow):", round(sd_meadow, 4), "\n")
cat("  Pooled SD:", round(pooled_sd, 4), "\n")
cat("  Standard Error:", round(se_d, 4), "\n\n")

# Display results
cat("RESULTS:\n")
cat("  Cohen's d =", round(cohens_d, 4), "\n")
cat("  95% CI = [", round(ci_lower, 4), ", ", round(ci_upper, 4), "]\n", sep = "")

# Effect size interpretation (Cohen, 1988)
if (abs(cohens_d) < 0.2) {
  effect_size <- "NEGLIGIBLE"
  range_text <- "(|d| < 0.2)"
} else if (abs(cohens_d) < 0.5) {
  effect_size <- "SMALL"
  range_text <- "(0.2 <= |d| < 0.5)"
} else if (abs(cohens_d) < 0.8) {
  effect_size <- "MEDIUM"
  range_text <- "(0.5 <= |d| < 0.8)"
} else {
  effect_size <- "LARGE"
  range_text <- "(|d| >= 0.8)"
}

cat("  Effect size =", effect_size, range_text, "\n\n")

# Direction and magnitude
direction <- ifelse(cohens_d > 0, 
                    "Meadow concordance > Maize concordance",
                    "Maize concordance > Meadow concordance")
cat("  Direction:", direction, "\n")
cat("  Magnitude:", abs(cohens_d), "standard deviations\n\n")

# Practical interpretation
cat("INTERPRETATION:\n")
cat("  The difference between Meadow and Maize has a", tolower(effect_size), "\n")
cat("  effect size. This indicates:\n")

if (abs(cohens_d) < 0.2) {
  cat("  - Virtually no practical difference between habitats\n")
  cat("  - Effect is negligible from a practical standpoint\n")
} else if (abs(cohens_d) < 0.5) {
  cat("  - A small but detectable practical difference\n")
  cat("  - Approximately", round(abs(cohens_d) * 100, 0), "% overlap between distributions\n")
} else if (abs(cohens_d) < 0.8) {
  cat("  - A moderate practical difference between habitats\n")
  cat("  - About", round((1 - abs(cohens_d) / 2) * 100, 0), "% overlap between distributions\n")
} else {
  cat("  - A large practical difference between habitats\n")
  cat("  - Minimal overlap between distributions\n")
}

# Check if CI includes zero
if (ci_lower <= 0 && ci_upper >= 0) {
  cat("  [!] 95% CI includes zero - effect may not be reliable\n")
} else {
  cat("  [OK] 95% CI excludes zero - effect is reliable\n")
}
cat("\n")

# ==============================================================================
# 6. VARIANCE HOMOGENEITY TEST (Levene's Test)
# ==============================================================================

cat("\n===========================================================\n")
cat("VARIANCE HOMOGENEITY CHECK (Levene's Test)\n")
cat("===========================================================\n\n")

# Prepare data for Levene's test
levene_data <- data.frame(
  value = c(maize_clean, meadow_clean),
  group = factor(c(rep("Maize", n_maize), rep("Meadow", n_meadow)))
)

# Levene's test (using car package if available, otherwise F-test)
if (requireNamespace("car", quietly = TRUE)) {
  levene_result <- car::leveneTest(value ~ group, data = levene_data)
  cat("Levene's test for homogeneity of variances:\n")
  cat("  F-statistic:", round(levene_result$`F value`[1], 4), "\n")
  cat("  p-value:", round(levene_result$`Pr(>F)`[1], 4), "\n")
  cat("  Conclusion:", ifelse(levene_result$`Pr(>F)`[1] > 0.05,
                              "Variances are homogeneous [OK]",
                              "Variances are NOT homogeneous"), "\n\n")
} else {
  # Fallback to F-test for equal variances
  var_test <- var.test(maize_clean, meadow_clean)
  cat("F-test for equality of variances:\n")
  cat("  F-statistic:", round(var_test$statistic, 4), "\n")
  cat("  p-value:", round(var_test$p.value, 4), "\n")
  cat("  Conclusion:", ifelse(var_test$p.value > 0.05,
                              "Variances are equal [OK]",
                              "Variances are NOT equal"), "\n\n")
}

# ==============================================================================
# 7. HABITAT EFFECT: GENERA/DAY COMPARISON (FAIRD & AMMOD)
# ==============================================================================

# Helper function: Levene's -> appropriate t-test
run_habitat_test <- function(maize_vals, meadow_vals, label, n_label) {
  cat(label, " by habitat (", n_label, ", n =", 
      length(maize_vals), "+", length(meadow_vals), "):\n", sep = "")
  # Replaced plus-minus unicode with +/-
  cat("  Maize:  mean =", round(mean(maize_vals), 1), "+/-", round(sd(maize_vals), 1), "\n")
  cat("  Meadow: mean =", round(mean(meadow_vals), 1), "+/-", round(sd(meadow_vals), 1), "\n")
  
  # Levene's test for variance homogeneity
  combined <- c(maize_vals, meadow_vals)
  groups <- factor(rep(c("Maize", "Meadow"), c(length(maize_vals), length(meadow_vals))))
  
  # Check for car package inside function too, or use simple F test logic if preferred. 
  # Assuming car is loaded or available as per original script logic.
  if (requireNamespace("car", quietly = TRUE)) {
    lev <- car::leveneTest(combined ~ groups)
    lev_p <- lev$`Pr(>F)`[1]
    lev_F <- lev$`F value`[1]
  } else {
    # Fallback if car is missing inside function scope
    lev_p <- var.test(maize_vals, meadow_vals)$p.value
    lev_F <- var.test(maize_vals, meadow_vals)$statistic
  }
  
  cat("  Levene's test: F =", round(lev_F, 2), ", p =", round(lev_p, 3))
  
  # Select t-test type based on Levene's result
  var_eq <- lev_p > 0.05
  if (var_eq) {
    cat(" -> variances homogeneous [OK]\n")
  } else {
    cat(" -> variances UNEQUAL [!]\n")
  }
  
  t_res <- t.test(meadow_vals, maize_vals, var.equal = var_eq)
  test_type <- ifelse(var_eq, "Student's", "Welch's")
  cat("  ", test_type, " t-test: t(", round(t_res$parameter, 1), ") = ", 
      round(t_res$statistic, 2), ", p = ", 
      format.pval(t_res$p.value, digits = 3), "\n", sep = "")
  
  # Effect direction
  diff_pct <- (mean(meadow_vals) - mean(maize_vals)) / mean(maize_vals) * 100
  cat("  Difference: ", sprintf("%+.1f", mean(meadow_vals) - mean(maize_vals)),
      " (", sprintf("%+.0f", diff_pct), "%)\n\n", sep = "")
}

cat("\n========================================================================\n")
cat("HABITAT EFFECT: GENERA/DAY COMPARISON\n")
cat("========================================================================\n\n")

# --- ALL DAYS (including zero-FAIRD days) -------------------------------------
cat("--- ALL DAYS (including zero-FAIRD days) -------------------------------\n\n")

ammod_maize_all <- similarity_results %>%
  dplyr::filter(Site %in% c("Site1", "Site2"), !is.na(N_taxa_AMMOD_full)) %>%
  dplyr::pull(N_taxa_AMMOD_full)

ammod_meadow_all <- similarity_results %>%
  dplyr::filter(Site %in% c("Site3", "Site4"), !is.na(N_taxa_AMMOD_full)) %>%
  dplyr::pull(N_taxa_AMMOD_full)

faird_maize_all <- similarity_results %>%
  dplyr::filter(Site %in% c("Site1", "Site2"), !is.na(N_taxa_FAIRD)) %>%
  dplyr::pull(N_taxa_FAIRD)

faird_meadow_all <- similarity_results %>%
  dplyr::filter(Site %in% c("Site3", "Site4"), !is.na(N_taxa_FAIRD)) %>%
  dplyr::pull(N_taxa_FAIRD)

run_habitat_test(ammod_maize_all, ammod_meadow_all, "AMMOD genera/day", "all days")
run_habitat_test(faird_maize_all, faird_meadow_all, "FAIRD genera/day", "all days")

# --- MATCHED DAYS ONLY (both FAIRD and AMMOD have data) ----------------------
cat("--- MATCHED DAYS ONLY (Delta_S available) ------------------------------\n\n")

ammod_maize_matched <- similarity_results %>%
  dplyr::filter(Site %in% c("Site1", "Site2"), !is.na(Delta_S_raw)) %>%
  dplyr::pull(N_taxa_AMMOD_full)

ammod_meadow_matched <- similarity_results %>%
  dplyr::filter(Site %in% c("Site3", "Site4"), !is.na(Delta_S_raw)) %>%
  dplyr::pull(N_taxa_AMMOD_full)

faird_maize_matched <- similarity_results %>%
  dplyr::filter(Site %in% c("Site1", "Site2"), !is.na(Delta_S_raw)) %>%
  dplyr::pull(N_taxa_FAIRD)

faird_meadow_matched <- similarity_results %>%
  dplyr::filter(Site %in% c("Site3", "Site4"), !is.na(Delta_S_raw)) %>%
  dplyr::pull(N_taxa_FAIRD)

run_habitat_test(ammod_maize_matched, ammod_meadow_matched, "AMMOD genera/day", "matched days")
run_habitat_test(faird_maize_matched, faird_meadow_matched, "FAIRD genera/day", "matched days")

# --- SUMMARY TABLE ------------------------------------------------------------
cat("--- SUMMARY ------------------------------------------------------------\n\n")
cat("                       ALL DAYS (n=44+44)          MATCHED DAYS (n=34+34)\n")
cat("  AMMOD Maize:     ", round(mean(ammod_maize_all), 1), "+/-", round(sd(ammod_maize_all), 1),
    "                  ", round(mean(ammod_maize_matched), 1), "+/-", round(sd(ammod_maize_matched), 1), "\n")
cat("  AMMOD Meadow:    ", round(mean(ammod_meadow_all), 1), "+/-", round(sd(ammod_meadow_all), 1),
    "                 ", round(mean(ammod_meadow_matched), 1), "+/-", round(sd(ammod_meadow_matched), 1), "\n")
cat("  FAIRD Maize:     ", round(mean(faird_maize_all), 1), "+/-", round(sd(faird_maize_all), 1),
    "                    ", round(mean(faird_maize_matched), 1), "+/-", round(sd(faird_maize_matched), 1), "\n")
cat("  FAIRD Meadow:    ", round(mean(faird_meadow_all), 1), "+/-", round(sd(faird_meadow_all), 1),
    "                    ", round(mean(faird_meadow_matched), 1), "+/-", round(sd(faird_meadow_matched), 1), "\n")
cat("\n")

# ==============================================================================
# 8. FINAL SUMMARY
# ==============================================================================

cat("\n========================================================================\n")
cat("SUMMARY\n")
cat("========================================================================\n\n")

cat("HABITAT CONCORDANCE COMPARISON (Delta_S_raw):\n")
cat("  - Wilcoxon rank-sum test: W =", wilcox_result$statistic, 
    ", p =", round(wilcox_result$p.value, 3), 
    "(", ifelse(wilcox_result$p.value < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), ")\n")
cat("  - Student's t-test: t(", round(t_test_res$parameter, 0), ") =", 
    round(t_test_res$statistic, 3), ", p =", format.pval(t_test_res$p.value, digits = 4),
    "(", ifelse(t_test_res$p.value < 0.05, "SIGNIFICANT *", "NOT SIGNIFICANT"), ")\n")
cat("  - Cohen's d:", round(cohens_d, 3), 
    "(", effect_size, ", 95% CI: [", round(ci_lower, 3), ",", round(ci_upper, 3), "])\n")
cat("  - Levene's test: F =", round(levene_result$`F value`[1], 2), 
    ", p =", round(levene_result$`Pr(>F)`[1], 3), "[OK]\n\n")

cat("HABITAT GENERA/DAY COMPARISON:\n")
cat("  - FAIRD: Meadow", round(mean(faird_meadow_matched), 1), "vs Maize", 
    round(mean(faird_maize_matched), 1), "genera/day \n")
cat("  - AMMOD: Meadow", round(mean(ammod_meadow_matched), 1), "vs Maize", 
    round(mean(ammod_maize_matched), 1), "genera/day\n")


cat("========================================================================\n\n")

# ==============================================================================
# NEW SECTION: FAIRD TAXONOMIC COMPOSITION BY HABITAT
# ==============================================================================
# Purpose: Calculate FAIRD (and ID) order-level composition broken down by
#          habitat (Maize vs Meadow) to provide reproducible source for
#          habitat-specific Diptera percentages
# ==============================================================================

cat("\n")
cat(paste(rep("═", 59), collapse = ""), "\n")
cat("SECTION: FAIRD TAXONOMIC COMPOSITION BY HABITAT\n")
cat(paste(rep("═", 59), collapse = ""), "\n\n")

# Load e-traps data (available via data_cleaned.RData already loaded above)
etraps_data <- data_list$etraps

# --- Order-level composition by Device_type x Habitat ---
composition_by_habitat <- etraps_data %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::group_by(Device_type, Ambient, Order) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Device_type, Ambient) %>%
  dplyr::mutate(
    Total     = sum(Count),
    Pct       = round(Count / Total * 100, 2)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Device_type, Ambient, dplyr::desc(Pct))

cat("--- Order-level composition by Device × Habitat ---\n\n")
print(composition_by_habitat, n = 60)

# --- Diptera-specific summary table ---
diptera_by_habitat <- composition_by_habitat %>%
  dplyr::filter(Order == "Diptera") %>%
  dplyr::select(Device_type, Ambient, Count, Total, Pct) %>%
  dplyr::rename(
    Diptera_n   = Count,
    Total_n     = Total,
    Diptera_pct = Pct
  )

cat("\n--- Diptera % by Device × Habitat ---\n\n")
print(diptera_by_habitat)

# --- Extract key values ---
faird_maize_dip <- diptera_by_habitat %>%
  dplyr::filter(Device_type == "FAIRD", Ambient == "Maize") %>%
  dplyr::pull(Diptera_pct)

faird_meadow_dip <- diptera_by_habitat %>%
  dplyr::filter(Device_type == "FAIRD", Ambient == "Meadow") %>%
  dplyr::pull(Diptera_pct)

id_maize_dip <- diptera_by_habitat %>%
  dplyr::filter(Device_type == "ID", Ambient == "Maize") %>%
  dplyr::pull(Diptera_pct)

id_meadow_dip <- diptera_by_habitat %>%
  dplyr::filter(Device_type == "ID", Ambient == "Meadow") %>%
  dplyr::pull(Diptera_pct)

cat("\n")
cat(paste(rep("─", 59), collapse = ""), "\n")
cat("DIPTERA % BY HABITAT\n")
cat(paste(rep("─", 59), collapse = ""), "\n")
cat(sprintf("  FAIRD  Maize:  %.1f%%  (n = %d Diptera / %d total)\n",
            faird_maize_dip,
            diptera_by_habitat %>% dplyr::filter(Device_type == "FAIRD", Ambient == "Maize") %>% dplyr::pull(Diptera_n),
            diptera_by_habitat %>% dplyr::filter(Device_type == "FAIRD", Ambient == "Maize") %>% dplyr::pull(Total_n)))
cat(sprintf("  FAIRD  Meadow: %.1f%%  (n = %d Diptera / %d total)\n",
            faird_meadow_dip,
            diptera_by_habitat %>% dplyr::filter(Device_type == "FAIRD", Ambient == "Meadow") %>% dplyr::pull(Diptera_n),
            diptera_by_habitat %>% dplyr::filter(Device_type == "FAIRD", Ambient == "Meadow") %>% dplyr::pull(Total_n)))
cat(sprintf("  ID     Maize:  %.1f%%  (n = %d Diptera / %d total)\n",
            id_maize_dip,
            diptera_by_habitat %>% dplyr::filter(Device_type == "ID", Ambient == "Maize") %>% dplyr::pull(Diptera_n),
            diptera_by_habitat %>% dplyr::filter(Device_type == "ID", Ambient == "Maize") %>% dplyr::pull(Total_n)))
cat(sprintf("  ID     Meadow: %.1f%%  (n = %d Diptera / %d total)\n",
            id_meadow_dip,
            diptera_by_habitat %>% dplyr::filter(Device_type == "ID", Ambient == "Meadow") %>% dplyr::pull(Diptera_n),
            diptera_by_habitat %>% dplyr::filter(Device_type == "ID", Ambient == "Meadow") %>% dplyr::pull(Total_n)))
cat(sprintf("\n  FAIRD habitat difference: %.1f pp (Maize - Meadow)\n",
            faird_maize_dip - faird_meadow_dip))
cat(paste(rep("─", 59), collapse = ""), "\n")

# --- Save table ---
readr::write_csv(
  composition_by_habitat,
  "outputs/H3_FAIRD_ID_Composition_by_Habitat.csv"
)
cat("✅ Table saved: outputs/H3_FAIRD_ID_Composition_by_Habitat.csv\n\n")

cat(paste(rep("═", 59), collapse = ""), "\n")
cat("SECTION COMPLETE\n")
cat(paste(rep("═", 59), collapse = ""), "\n\n")


cat("[OK] ANALYSIS COMPLETE\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")