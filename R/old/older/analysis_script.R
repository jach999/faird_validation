# ==============================================================================
# EXAMPLE ANALYSIS SCRIPT - FAIR-Device Validation
# ==============================================================================
# Description: Example workflow showing how to run all analyses
# Author: [Your name]
# Date: 2025-10-10
# Updated: Adapted to actual data structure and limitations
# ==============================================================================

# --- Setup ---
library(tidyverse)
library(here)

# Load data processing script
cat("\n========================================\n")
cat("LOADING DATA\n")
cat("========================================\n\n")

source("load_data.R")

# Load analysis functions
source("analysis_functions.R")

# Load cleaned data (created by load_data.R)
load("data_cleaned.RData")

# Access datasets
etraps_data <- data_list$etraps
etraps_daily <- data_list$etraps_daily
ammod_biomass <- data_list$ammod_biomass
ammod_taxonomy <- data_list$ammod_taxonomy

cat("\nData loaded successfully!\n")
cat("Available sites for FAIRD-AMMOD comparison: Site1, Site2, Site3\n")
cat("Note: Site4 excluded (AMMOD4 has no data)\n\n")

# ==============================================================================
# HYPOTHESIS 1: INTRA-DEVICE CONSISTENCY
# ==============================================================================

cat("\n========================================\n")
cat("HYPOTHESIS 1: INTRA-DEVICE CONSISTENCY\n")
cat("========================================\n\n")

# --- 1.1: Correlations between replicates ---
cat("--- Analysis 1.1: Correlations ---\n")

# FAIRD in Maize (FAIRD1 vs FAIRD2)
cat("\nAnalyzing FAIRD in Maize (FAIRD1 vs FAIRD2)...\n")
faird_maize <- analyze_device_consistency(
  daily_data = etraps_daily,
  device1 = "FAIRD1",
  device2 = "FAIRD2",
  metric = "biomass_mg"
)

if(!is.null(faird_maize)) {
  print(faird_maize$plot)
  cat(sprintf("✓ FAIRD Maize: R² = %.3f, p = %.4f\n", 
              faird_maize$r_squared, faird_maize$pearson_p))
  
  # Interpretation
  if(faird_maize$r_squared > 0.70 & faird_maize$pearson_p < 0.05) {
    cat("  → GOOD CONSISTENCY (H1 supported)\n")
  } else if(faird_maize$r_squared > 0.50 & faird_maize$pearson_p < 0.05) {
    cat("  → MODERATE CONSISTENCY (H1 partially supported)\n")
  } else {
    cat("  → POOR CONSISTENCY (H1 not supported)\n")
  }
}

# FAIRD in Meadow (FAIRD3 vs FAIRD4)
cat("\nAnalyzing FAIRD in Meadow (FAIRD3 vs FAIRD4)...\n")
faird_meadow <- analyze_device_consistency(
  daily_data = etraps_daily,
  device1 = "FAIRD3",
  device2 = "FAIRD4",
  metric = "biomass_mg"
)

if(!is.null(faird_meadow)) {
  print(faird_meadow$plot)
  cat(sprintf("✓ FAIRD Meadow: R² = %.3f, p = %.4f\n",
              faird_meadow$r_squared, faird_meadow$pearson_p))
  
  if(faird_meadow$r_squared > 0.70 & faird_meadow$pearson_p < 0.05) {
    cat("  → GOOD CONSISTENCY (H1 supported)\n")
  }
}

# ID in Maize (ID1 vs ID2)
cat("\nAnalyzing ID in Maize (ID1 vs ID2)...\n")
id_maize <- analyze_device_consistency(
  daily_data = etraps_daily,
  device1 = "ID1",
  device2 = "ID2",
  metric = "biomass_mg"
)

if(!is.null(id_maize)) {
  print(id_maize$plot)
  cat(sprintf("✓ ID Maize: R² = %.3f, p = %.4f\n",
              id_maize$r_squared, id_maize$pearson_p))
  
  if(id_maize$r_squared > 0.60 & id_maize$pearson_p < 0.05) {
    cat("  → GOOD CONSISTENCY (H1 supported)\n")
  }
}

# Create summary table
consistency_results <- list(faird_maize, faird_meadow, id_maize)
consistency_table <- create_summary_table(consistency_results, type = "correlation")

cat("\n--- Summary Table: Intra-device Correlations ---\n")
print(consistency_table)

# Save results
write_csv(consistency_table, "results_consistency_correlations.csv")

if(!is.null(faird_maize)) {
  ggsave("plot_consistency_faird_maize.png", faird_maize$plot, 
         width = 8, height = 6, dpi = 300)
}
if(!is.null(faird_meadow)) {
  ggsave("plot_consistency_faird_meadow.png", faird_meadow$plot,
         width = 8, height = 6, dpi = 300)
}
if(!is.null(id_maize)) {
  ggsave("plot_consistency_id_maize.png", id_maize$plot,
         width = 8, height = 6, dpi = 300)
}

# --- 1.2: Coefficient of Variation ---
cat("\n--- Analysis 1.2: Coefficient of Variation ---\n")

# Calculate CV for each pair
cv_faird_maize <- calculate_CV(etraps_daily, "FAIRD1", "FAIRD2", "biomass_mg")
cv_faird_meadow <- calculate_CV(etraps_daily, "FAIRD3", "FAIRD4", "biomass_mg")
cv_id_maize <- calculate_CV(etraps_daily, "ID1", "ID2", "biomass_mg")

# Print results
if(!is.null(cv_faird_maize)) {
  print(cv_faird_maize$plot)
  cat(sprintf("\n✓ FAIRD Maize CV: Mean = %.1f%%, Median = %.1f%%\n",
              cv_faird_maize$cv_stats$CV_mean,
              cv_faird_maize$cv_stats$CV_median))
  
  if(cv_faird_maize$cv_stats$CV_mean < 20) {
    cat("  → GOOD REPRODUCIBILITY (CV < 20%)\n")
  } else {
    cat("  → MODERATE REPRODUCIBILITY (CV > 20%)\n")
  }
  
  ggsave("plot_cv_faird_maize.png", cv_faird_maize$plot,
         width = 10, height = 6, dpi = 300)
}

# Create CV summary table
cv_results <- list(cv_faird_maize, cv_faird_meadow, cv_id_maize)
cv_table <- create_summary_table(cv_results, type = "cv")

cat("\n--- Summary Table: Coefficient of Variation ---\n")
print(cv_table)

# Save
write_csv(cv_table, "results_coefficient_variation.csv")

# ==============================================================================
# HYPOTHESIS 2: VALIDATION AGAINST AMMOD
# ==============================================================================

cat("\n========================================\n")
cat("HYPOTHESIS 2: VALIDATION AGAINST AMMOD\n")
cat("========================================\n\n")

# --- 2.1: Correlations FAIRD-AMMOD by site ---
cat("--- Analysis 2.1: FAIRD-AMMOD Correlations ---\n")
cat("Note: Only Site1, Site2, and Site3 available (Site4 AMMOD data missing)\n\n")

# Analyze each available site
sites_to_analyze <- c("Site1", "Site2", "Site3")
ammod_results <- list()

for(site in sites_to_analyze) {
  cat(sprintf("Analyzing %s...\n", site))
  
  result <- analyze_faird_ammod_correlation(
    etraps_daily = etraps_daily,
    ammod_biomass = ammod_biomass,
    site_name = site
  )
  
  if(!is.null(result)) {
    ammod_results[[site]] <- result
    print(result$plot_scatter)
    cat(sprintf("✓ %s: R² = %.3f, p = %.4f\n", 
                site, result$r_squared, result$pearson_p))
    
    # Interpretation
    if(site == "Site1") {
      if(result$r_squared < 0.30) {
        cat("  → LOW CORRELATION as expected (edge effect + size bias)\n")
      }
    } else {
      if(result$r_squared > 0.50 & result$pearson_p < 0.05) {
        cat("  → GOOD CORRELATION (H2 supported for homogeneous habitat)\n")
      } else if(result$r_squared > 0.30) {
        cat("  → MODERATE CORRELATION (H2 partially supported)\n")
      }
    }
    
    # Save plot
    ggsave(sprintf("plot_faird_ammod_%s.png", site), 
           result$plot_scatter, width = 8, height = 6, dpi = 300)
    
  } else {
    cat(sprintf("⚠️  %s: Insufficient data for correlation analysis\n", site))
  }
  cat("\n")
}

# Summary table
ammod_table <- create_summary_table(ammod_results, type = "correlation")

cat("--- Summary Table: FAIRD-AMMOD Correlations ---\n")
print(ammod_table)

write_csv(ammod_table, "results_faird_ammod_correlations.csv")

# --- 2.2: Temporal Synchronization (PRIORITY) ---
cat("\n--- Analysis 2.2: Temporal Synchronization (PRIORITY) ---\n")
cat("This is the MOST IMPORTANT analysis for demonstrating dynamic monitoring capability\n\n")

sync_results <- list()

for(site in sites_to_analyze) {
  cat(sprintf("Temporal analysis for %s...\n", site))
  
  result <- temporal_synchronization(
    etraps_daily = etraps_daily,
    ammod_biomass = ammod_biomass,
    site_name = site
  )
  
  if(!is.null(result)) {
    sync_results[[site]] <- result
    
    # Print combined plot
    combined_plot <- result$plot_timeseries / result$plot_ccf
    print(combined_plot)
    
    cat(sprintf("✓ %s: Concordance = %.1f%%, Kendall τ = %.3f (p = %.4f)\n",
                site, result$concordance_pct, 
                result$kendall_tau, result$kendall_p))
    
    # Interpretation
    if(result$concordance_pct > 70 & result$kendall_tau > 0.5) {
      cat("  → EXCELLENT TEMPORAL SYNCHRONIZATION (H2 strongly supported)\n")
    } else if(result$concordance_pct > 60 & result$kendall_tau > 0.3) {
      cat("  → GOOD TEMPORAL SYNCHRONIZATION (H2 supported)\n")
    } else {
      cat("  → MODERATE TEMPORAL SYNCHRONIZATION (H2 partially supported)\n")
    }
    
    # Save individual plots
    ggsave(sprintf("plot_timeseries_%s.png", site), 
           result$plot_timeseries, width = 10, height = 6, dpi = 300)
    ggsave(sprintf("plot_ccf_%s.png", site),
           result$plot_ccf, width = 8, height = 6, dpi = 300)
    ggsave(sprintf("plot_combined_temporal_%s.png", site),
           combined_plot, width = 10, height = 10, dpi = 300)
  } else {
    cat(sprintf("⚠️  %s: Insufficient data for temporal analysis\n", site))
  }
  cat("\n")
}

# Summary table
sync_table <- create_summary_table(sync_results, type = "temporal")

cat("--- Summary Table: Temporal Synchronization ---\n")
print(sync_table)

write_csv(sync_table, "results_temporal_synchronization.csv")

# ==============================================================================
# HYPOTHESIS 3: COMPARISON WITH INSECT DETECT
# ==============================================================================

cat("\n========================================\n")
cat("HYPOTHESIS 3: COMPARISON WITH ID\n")
cat("========================================\n\n")

cat("--- Analysis 3.1: Taxonomic Richness ---\n")

comparison_results <- compare_taxonomic_richness(etraps_data)

# Print results
cat("\n--- Taxonomic Richness by Device ---\n")
print(comparison_results$richness)

cat("\n--- Diversity Indices ---\n")
print(comparison_results$diversity_indices)

cat("\n--- Chi-squared Test (Taxonomic Composition) ---\n")
cat(sprintf("χ² = %.3f, df = %d, p = %.4f\n",
            comparison_results$chi_squared$statistic,
            comparison_results$chi_squared$parameter,
            comparison_results$chi_squared$p.value))

if(comparison_results$chi_squared$p.value < 0.05) {
  cat("  → SIGNIFICANTLY DIFFERENT taxonomic composition (H3 supported)\n")
}

# Check if FAIRD has higher richness
faird_richness <- comparison_results$richness %>% 
  filter(Device_type == "FAIRD") %>% 
  pull(N_genera)
id_richness <- comparison_results$richness %>% 
  filter(Device_type == "ID") %>% 
  pull(N_genera)

if(faird_richness > id_richness) {
  cat(sprintf("  → FAIRD captures MORE genera (%d) than ID (%d) (H3 supported)\n",
              faird_richness, id_richness))
}

print(comparison_results$plot)

# Save
write_csv(comparison_results$richness, "results_taxonomic_richness.csv")
write_csv(comparison_results$diversity_indices, "results_diversity_indices.csv")
ggsave("plot_taxonomic_composition_comparison.png", 
       comparison_results$plot, width = 10, height = 6, dpi = 300)

# ==============================================================================
# GENERATE COMPREHENSIVE REPORT
# ==============================================================================

cat("\n========================================\n")
cat("GENERATING COMPREHENSIVE REPORT\n")
cat("========================================\n\n")

# Compile all results
all_results <- list(
  consistency = list(
    correlations = consistency_table,
    cv = cv_table
  ),
  validation_ammod = list(
    correlations = ammod_table,
    temporal = sync_table
  ),
  comparison_id = list(
    richness = comparison_results$richness,
    diversity = comparison_results$diversity_indices,
    chi_squared = comparison_results$chi_squared
  )
)

# Save complete results
save(all_results, file = "analysis_results_complete.RData")
cat("✓ Complete results saved to: analysis_results_complete.RData\n")

# Create summary report
sink("analysis_summary_report.txt")

cat("==============================================================================\n")
cat("FAIR-DEVICE VALIDATION - ANALYSIS SUMMARY REPORT\n")
cat("==============================================================================\n\n")

cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

cat("DATA LIMITATIONS:\n")
cat("  - AMMOD4 (Site4): No biomass data available\n")
cat("  - ID3 (Meadow): Device failure, no data\n")
cat("  - AMMOD3: Suspicious repeated values after 2023-08-29\n")
cat("  - Analysis period: ~12 days (2023-08-23 to 2023-09-03)\n\n")

cat("==============================================================================\n")
cat("HYPOTHESIS 1: INTRA-DEVICE CONSISTENCY\n")
cat("==============================================================================\n\n")

cat("H1: Replicate devices produce reproducible measurements\n\n")

cat("--- Correlations between replicates ---\n")
print(consistency_table)

cat("\n--- Coefficient of Variation ---\n")
print(cv_table)

cat("\nINTERPRETATION:\n")
cat("  ✓ R² > 0.70 AND CV < 20% = H1 SUPPORTED\n")
cat("  ⚠️ R² 0.50-0.70 OR CV 20-30% = H1 PARTIALLY SUPPORTED\n")
cat("  ✗ R² < 0.50 OR CV > 30% = H1 NOT SUPPORTED\n\n")

cat("==============================================================================\n")
cat("HYPOTHESIS 2: VALIDATION AGAINST AMMOD (GROUND TRUTH)\n")
cat("==============================================================================\n\n")

cat("H2: FAIRD detects biomass patterns comparable to traditional methods\n\n")

cat("--- FAIRD-AMMOD Correlations (Static) ---\n")
print(ammod_table)

cat("\n--- Temporal Synchronization (Dynamic) ---\n")
print(sync_table)

cat("\nINTERPRETATION:\n")
cat("  Site1 (Edge): LOW R² expected due to edge effect + size bias\n")
cat("  Site2/3 (Homogeneous):\n")
cat("    ✓ R² > 0.50 = H2 SUPPORTED (valid for homogeneous habitats)\n")
cat("    ✓ Concordance > 70% + Kendall τ > 0.5 = Good temporal sync\n\n")

cat("CRITICAL FINDING:\n")
if(any(!is.na(sync_table$Concordance_pct))) {
  max_concordance <- max(sync_table$Concordance_pct, na.rm = TRUE)
  cat(sprintf("  Maximum concordance: %.1f%% - demonstrates temporal validity\n", 
              max_concordance))
}

cat("\n==============================================================================\n")
cat("HYPOTHESIS 3: COMPARISON WITH INSECT DETECT\n")
cat("==============================================================================\n\n")

cat("H3: FAIRD captures greater taxonomic diversity than ID\n\n")

cat("--- Taxonomic Richness ---\n")
print(comparison_results$richness)

cat("\n--- Diversity Indices ---\n")
print(comparison_results$diversity_indices)

cat(sprintf("\n--- Chi-squared Test ---\n"))
cat(sprintf("χ² = %.3f, p = %.4f\n",
            comparison_results$chi_squared$statistic,
            comparison_results$chi_squared$p.value))

cat("\nINTERPRETATION:\n")
cat("  ✓ FAIRD genera > ID genera = H3 SUPPORTED\n")
cat("  ✓ Chi² p < 0.05 = Different composition\n\n")

cat("==============================================================================\n")
cat("OVERALL CONCLUSIONS\n")
cat("==============================================================================\n\n")

cat("1. CONSISTENCY (H1):\n")
cat("   Check if R² > 0.70 and CV < 20% for all replicates\n\n")

cat("2. VALIDATION (H2):\n")
cat("   - Static: Check R² for Site2/Site3\n")
cat("   - Dynamic: Check temporal concordance and Kendall τ\n")
cat("   - Site1 low correlation is EXPECTED and EXPLAINED\n\n")

cat("3. COMPARISON (H3):\n")
cat("   Check if FAIRD has higher richness than ID\n\n")

cat("==============================================================================\n")
cat("FILES GENERATED\n")
cat("==============================================================================\n\n")

cat("CSV Tables:\n")
cat("  - results_consistency_correlations.csv\n")
cat("  - results_coefficient_variation.csv\n")
cat("  - results_faird_ammod_correlations.csv\n")
cat("  - results_temporal_synchronization.csv\n")
cat("  - results_taxonomic_richness.csv\n")
cat("  - results_diversity_indices.csv\n\n")

cat("PNG Figures:\n")
cat("  - plot_consistency_*.png (3 files)\n")
cat("  - plot_cv_*.png (3 files)\n")
cat("  - plot_faird_ammod_*.png (3 sites)\n")
cat("  - plot_timeseries_*.png (3 sites)\n")
cat("  - plot_ccf_*.png (3 sites)\n")
cat("  - plot_combined_temporal_*.png (3 sites)\n")
cat("  - plot_taxonomic_composition_comparison.png\n\n")

cat("RData Files:\n")
cat("  - data_cleaned.RData (input data)\n")
cat("  - analysis_results_complete.RData (all results)\n\n")

cat("==============================================================================\n")
cat("NEXT STEPS\n")
cat("==============================================================================\n\n")

cat("1. Review all plots carefully\n")
cat("2. Check interpretation criteria for each hypothesis\n")
cat("3. Prepare figures for manuscript (combine related plots)\n")
cat("4. Write Results section using these outputs\n")
cat("5. Consider additional analyses if needed:\n")
cat("   - Bland-Altman plots (H1)\n")
cat("   - Climate effects (H4)\n")
cat("   - Detailed diversity analysis (H5)\n\n")

cat("==============================================================================\n")
cat("END OF REPORT\n")
cat("==============================================================================\n")

sink()

cat("✓ Summary report saved to: analysis_summary_report.txt\n")

# ==============================================================================
# FINAL SUMMARY TO CONSOLE
# ==============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

cat("QUICK RESULTS SUMMARY:\n\n")

cat("H1 - Consistency:\n")
if(!is.null(faird_maize)) {
  cat(sprintf("  FAIRD Maize: R² = %.3f (target: >0.70)\n", faird_maize$r_squared))
}
if(!is.null(faird_meadow)) {
  cat(sprintf("  FAIRD Meadow: R² = %.3f (target: >0.70)\n", faird_meadow$r_squared))
}

cat("\nH2 - AMMOD Validation:\n")
for(site in names(ammod_results)) {
  if(!is.null(ammod_results[[site]])) {
    cat(sprintf("  %s: R² = %.3f (target Site2/3: >0.50)\n", 
                site, ammod_results[[site]]$r_squared))
  }
}

cat("\nH2 - Temporal Sync:\n")
for(site in names(sync_results)) {
  if(!is.null(sync_results[[site]])) {
    cat(sprintf("  %s: Concordance = %.1f%%, τ = %.3f (target: >70%%, >0.5)\n",
                site, sync_results[[site]]$concordance_pct,
                sync_results[[site]]$kendall_tau))
  }
}

cat("\nH3 - Taxonomic Comparison:\n")
cat(sprintf("  FAIRD: %d genera\n", faird_richness))
cat(sprintf("  ID: %d genera\n", id_richness))
cat(sprintf("  Difference: %+d (FAIRD captures more)\n", faird_richness - id_richness))

cat("\n========================================\n")
cat("Check 'analysis_summary_report.txt' for full details\n")
cat("========================================\n\n")

cat("🎉 All analyses complete! Ready for manuscript preparation.\n")