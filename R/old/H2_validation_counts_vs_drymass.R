################################################################################
# H2 ALTERNATIVE ANALYSIS: FAIRD COUNTS vs AMMOD DRY MASS
################################################################################
# 
# Purpose: Explore correlation between raw measurements without inference
#   - FAIRD: abundance (counts) - direct detection
#   - AMMOD: Dry_mass_g - direct weighing
# 
# Rationale: Eliminates two sources of measurement error:
#   1. Biomass inference in FAIRD (from morphometry)
#   2. Dry-to-Live conversion in AMMOD (×4 factor)
# 
# Question: Is the poor correlation at Site1 due to measurement error
#           or genuine detection bias (size threshold)?
#
# Date: November 13, 2025
################################################################################

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)      # For robust regression

# ============================================================
# CONSOLE OUTPUT CAPTURE
# ============================================================
script_name <- sub("\\.R$", "", basename(sys.frame(1)$ofile))

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
# SECTION 1: LOAD DATA
################################################################################

cat("Loading source data...\n")

# Load e-traps daily data (includes counts/abundance)
etraps_daily <- readr::read_csv("etraps_daily_biomass_clean.csv", 
                                show_col_types = FALSE)

# Load AMMOD data (includes dry mass)
ammod_data <- readr::read_csv("ammod_biomass_for_stats.csv", 
                              show_col_types = FALSE)

cat("  E-traps daily data loaded:", nrow(etraps_daily), "rows\n")
cat("  AMMOD data loaded:", nrow(ammod_data), "rows\n")

################################################################################
# SECTION 2: PREPARE DATA FOR ANALYSIS
################################################################################

cat("\nPreparing data for counts vs dry_mass analysis...\n")

# Define overlap period (Aug 23 - Sep 3, 2023)
overlap_start <- as.Date("2023-08-23")
overlap_end <- as.Date("2023-09-03")

# Filter FAIRD data for overlap period and sites with AMMOD
faird_counts <- etraps_daily %>%
  dplyr::filter(
    Device_type == "FAIRD",
    Site %in% c("Site1", "Site2", "Site3"),
    Date >= overlap_start,
    Date <= overlap_end
  ) %>%
  dplyr::select(Date, Site, Ambient, Device,
                FAIRD_counts = abundance,
                FAIRD_biomass_inferred = biomass_mg) %>%
  dplyr::mutate(Date = as.Date(Date))

# Filter AMMOD data (already only has overlap period and 3 sites)
ammod_dry <- ammod_data %>%
  dplyr::filter(Site %in% c("Site1", "Site2", "Site3")) %>%
  dplyr::select(Date, Site, Device,
                AMMOD_dry_mass_g = Dry_mass_g,
                AMMOD_live_mass_mg = Live_mass) %>%
  dplyr::mutate(Date = as.Date(Date))

# Join FAIRD counts with AMMOD dry mass
h2_alternative <- dplyr::full_join(
  faird_counts, 
  ammod_dry,
  by = c("Date", "Site")
) %>%
  dplyr::arrange(Date, Site)

# Add log-transformed variables (for statistical analysis)
h2_alternative <- h2_alternative %>%
  dplyr::mutate(
    Log_counts = log10(FAIRD_counts + 1),
    Log_dry_mass = log10(AMMOD_dry_mass_g + 1),
    Log_faird_biomass = log10(FAIRD_biomass_inferred + 1),
    Log_ammod_live = log10(AMMOD_live_mass_mg + 1)
  )

cat("  Alternative H2 dataset created:", nrow(h2_alternative), "observations\n")
cat("  Sites included:", paste(unique(h2_alternative$Site), collapse = ", "), "\n")
cat("  Date range:", min(h2_alternative$Date), "to", max(h2_alternative$Date), "\n")

# Check for missing values
n_missing <- sum(is.na(h2_alternative$FAIRD_counts) | 
                   is.na(h2_alternative$AMMOD_dry_mass_g))
if (n_missing > 0) {
  cat("  WARNING:", n_missing, "observations with missing values\n")
}

################################################################################
# SECTION 3: ANALYSIS BY SITE - COUNTS vs DRY MASS
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("ANALYSIS 1: FAIRD COUNTS vs AMMOD DRY MASS (ALTERNATIVE)\n")
cat(rep("=", 80), "\n\n", sep = "")

# Initialize results storage
results_counts_dry <- list()

for (site in c("Site1", "Site2", "Site3")) {
  cat("Site:", site, "\n")
  cat(rep("-", 60), "\n", sep = "")
  
  # Filter data for this site
  data_site <- h2_alternative %>%
    dplyr::filter(Site == site, 
                  !is.na(FAIRD_counts),
                  !is.na(AMMOD_dry_mass_g))
  
  # Linear regression (log-transformed)
  model_log <- lm(Log_counts ~ Log_dry_mass, data = data_site)
  model_summary <- summary(model_log)
  
  # Robust regression (for comparison)
  model_robust <- MASS::rlm(Log_counts ~ Log_dry_mass, data = data_site)
  
  # Correlation tests
  cor_pearson <- cor.test(~ Log_counts + Log_dry_mass, 
                          data = data_site, 
                          method = "pearson")
  
  cor_spearman <- cor.test(~ Log_counts + Log_dry_mass,
                           data = data_site,
                           method = "spearman")
  
  # Store results
  results_counts_dry[[site]] <- list(
    site = site,
    n = nrow(data_site),
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared,
    slope = coef(model_log)[2],
    intercept = coef(model_log)[1],
    slope_se = summary(model_log)$coefficients[2, 2],
    p_value = model_summary$coefficients[2, 4],
    pearson_r = cor_pearson$estimate,
    pearson_p = cor_pearson$p.value,
    spearman_rho = cor_spearman$estimate,
    spearman_p = cor_spearman$p.value,
    robust_r2 = 1 - sum(model_robust$residuals^2) / 
      sum((data_site$Log_counts - mean(data_site$Log_counts))^2)
  )
  
  # Print results
  cat("  Sample size: n =", nrow(data_site), "\n")
  cat("  Linear regression (log-log):\n")
  cat("    R² =", round(model_summary$r.squared, 3), "\n")
  cat("    Slope =", round(coef(model_log)[2], 3),
      "± SE:", round(summary(model_log)$coefficients[2, 2], 3), "\n")
  cat("    p-value =", format.pval(model_summary$coefficients[2, 4], digits = 3), "\n")
  cat("  Pearson correlation:\n")
  cat("    r =", round(cor_pearson$estimate, 3),
      ", p =", format.pval(cor_pearson$p.value, digits = 3), "\n")
  cat("  Spearman correlation:\n")
  cat("    ρ =", round(cor_spearman$estimate, 3),
      ", p =", format.pval(cor_spearman$p.value, digits = 3), "\n")
  cat("\n")
}

################################################################################
# SECTION 4: ANALYSIS BY SITE - BIOMASS vs BIOMASS (ORIGINAL)
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("ANALYSIS 2: FAIRD BIOMASS vs AMMOD LIVE MASS (ORIGINAL)\n")
cat(rep("=", 80), "\n\n", sep = "")

# Initialize results storage
results_biomass_biomass <- list()

for (site in c("Site1", "Site2", "Site3")) {
  cat("Site:", site, "\n")
  cat(rep("-", 60), "\n", sep = "")
  
  # Filter data for this site
  data_site <- h2_alternative %>%
    dplyr::filter(Site == site,
                  !is.na(FAIRD_biomass_inferred),
                  !is.na(AMMOD_live_mass_mg))
  
  # Linear regression (log-transformed)
  model_log <- lm(Log_faird_biomass ~ Log_ammod_live, data = data_site)
  model_summary <- summary(model_log)
  
  # Robust regression
  model_robust <- MASS::rlm(Log_faird_biomass ~ Log_ammod_live, data = data_site)
  
  # Correlation tests
  cor_pearson <- cor.test(~ Log_faird_biomass + Log_ammod_live,
                          data = data_site,
                          method = "pearson")
  
  cor_spearman <- cor.test(~ Log_faird_biomass + Log_ammod_live,
                           data = data_site,
                           method = "spearman")
  
  # Store results
  results_biomass_biomass[[site]] <- list(
    site = site,
    n = nrow(data_site),
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared,
    slope = coef(model_log)[2],
    intercept = coef(model_log)[1],
    slope_se = summary(model_log)$coefficients[2, 2],
    p_value = model_summary$coefficients[2, 4],
    pearson_r = cor_pearson$estimate,
    pearson_p = cor_pearson$p.value,
    spearman_rho = cor_spearman$estimate,
    spearman_p = cor_spearman$p.value,
    robust_r2 = 1 - sum(model_robust$residuals^2) /
      sum((data_site$Log_faird_biomass - mean(data_site$Log_faird_biomass))^2)
  )
  
  # Print results
  cat("  Sample size: n =", nrow(data_site), "\n")
  cat("  Linear regression (log-log):\n")
  cat("    R² =", round(model_summary$r.squared, 3), "\n")
  cat("    Slope =", round(coef(model_log)[2], 3),
      "± SE:", round(summary(model_log)$coefficients[2, 2], 3), "\n")
  cat("    p-value =", format.pval(model_summary$coefficients[2, 4], digits = 3), "\n")
  cat("  Pearson correlation:\n")
  cat("    r =", round(cor_pearson$estimate, 3),
      ", p =", format.pval(cor_pearson$p.value, digits = 3), "\n")
  cat("  Spearman correlation:\n")
  cat("    ρ =", round(cor_spearman$estimate, 3),
      ", p =", format.pval(cor_spearman$p.value, digits = 3), "\n")
  cat("\n")
}

################################################################################
# SECTION 5: COMPARATIVE SUMMARY TABLE
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("COMPARATIVE SUMMARY: ALTERNATIVE vs ORIGINAL ANALYSIS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Create comparison table
comparison_table <- data.frame(
  Site = c("Site1", "Site2", "Site3"),
  
  # Alternative approach (counts vs dry mass)
  Alt_R2 = sapply(c("Site1", "Site2", "Site3"), 
                  function(s) results_counts_dry[[s]]$r_squared),
  Alt_Spearman = sapply(c("Site1", "Site2", "Site3"),
                        function(s) results_counts_dry[[s]]$spearman_rho),
  Alt_p_value = sapply(c("Site1", "Site2", "Site3"),
                       function(s) results_counts_dry[[s]]$p_value),
  
  # Original approach (biomass vs biomass)
  Orig_R2 = sapply(c("Site1", "Site2", "Site3"),
                   function(s) results_biomass_biomass[[s]]$r_squared),
  Orig_Spearman = sapply(c("Site1", "Site2", "Site3"),
                         function(s) results_biomass_biomass[[s]]$spearman_rho),
  Orig_p_value = sapply(c("Site1", "Site2", "Site3"),
                        function(s) results_biomass_biomass[[s]]$p_value),
  
  # Difference (Alternative - Original)
  Delta_R2 = sapply(c("Site1", "Site2", "Site3"),
                    function(s) results_counts_dry[[s]]$r_squared - 
                      results_biomass_biomass[[s]]$r_squared),
  Delta_Spearman = sapply(c("Site1", "Site2", "Site3"),
                          function(s) results_counts_dry[[s]]$spearman_rho -
                            results_biomass_biomass[[s]]$spearman_rho)
)

print(comparison_table, digits = 3)

cat("\nINTERPRETATION:\n")
for (i in 1:3) {
  site <- comparison_table$Site[i]
  delta_r2 <- comparison_table$Delta_R2[i]
  
  cat("\n", site, ":\n", sep = "")
  if (abs(delta_r2) < 0.05) {
    cat("  → Minimal difference (ΔR² =", round(delta_r2, 3), ")\n")
    cat("     Measurement error is NOT the primary issue\n")
    cat("     Detection bias likely dominates\n")
  } else if (delta_r2 > 0.05) {
    cat("  → IMPROVED correlation (ΔR² = +", round(delta_r2, 3), ")\n")
    cat("     Measurement error WAS contributing to poor correlation\n")
    cat("     Raw measurements reveal stronger relationship\n")
  } else {
    cat("  → WORSE correlation (ΔR² =", round(delta_r2, 3), ")\n")
    cat("     Unexpected result - investigate further\n")
  }
}

################################################################################
# SECTION 6: VISUALIZATION - COMPARATIVE PLOTS
################################################################################

cat("\n\nGenerating comparative figures...\n")

# Create output directory
if (!dir.exists("h2_alternative_figures")) {
  dir.create("h2_alternative_figures")
}

# Colors for sites
site_colors <- c("Site1" = "#E41A1C", "Site2" = "#377EB8", "Site3" = "#4DAF4A")

# FIGURE 1: Alternative approach (Counts vs Dry mass)
cat("  Creating Figure 1: Counts vs Dry mass...\n")

fig1 <- ggplot(h2_alternative, aes(x = Log_dry_mass, y = Log_counts, color = Site)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = site_colors) +
  labs(
    title = "H2 Alternative: FAIRD Counts vs AMMOD Dry Mass",
    subtitle = "Raw measurements (no biomass inference or conversion)",
    x = "log10(AMMOD Dry Mass + 1) [g]",
    y = "log10(FAIRD Counts + 1) [individuals]"
  ) +
  facet_wrap(~ Site, ncol = 3) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave("h2_alternative_figures/Fig1_Counts_vs_DryMass.png", 
       fig1, width = 12, height = 4, dpi = 300)

# FIGURE 2: Original approach (Biomass vs Biomass)
cat("  Creating Figure 2: Biomass vs Biomass...\n")

fig2 <- ggplot(h2_alternative, aes(x = Log_ammod_live, y = Log_faird_biomass, color = Site)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = site_colors) +
  labs(
    title = "H2 Original: FAIRD Biomass vs AMMOD Live Mass",
    subtitle = "With biomass inference (FAIRD) and dry-to-live conversion (AMMOD)",
    x = "log10(AMMOD Live Mass + 1) [mg]",
    y = "log10(FAIRD Biomass + 1) [mg]"
  ) +
  facet_wrap(~ Site, ncol = 3) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave("h2_alternative_figures/Fig2_Biomass_vs_Biomass.png",
       fig2, width = 12, height = 4, dpi = 300)

# FIGURE 3: R² Comparison barplot
cat("  Creating Figure 3: R² Comparison...\n")

comparison_long <- comparison_table %>%
  tidyr::pivot_longer(
    cols = c(Alt_R2, Orig_R2),
    names_to = "Analysis",
    values_to = "R_squared"
  ) %>%
  dplyr::mutate(
    Analysis = dplyr::recode(Analysis,
                             "Alt_R2" = "Counts vs Dry Mass",
                             "Orig_R2" = "Biomass vs Biomass")
  )

fig3 <- ggplot(comparison_long, aes(x = Site, y = R_squared, fill = Analysis)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(R_squared, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Counts vs Dry Mass" = "#FF7F00",
                               "Biomass vs Biomass" = "#984EA3")) +
  labs(
    title = "Comparison of R² Values: Alternative vs Original Analysis",
    subtitle = "Which approach yields stronger FAIRD-AMMOD correlation?",
    x = "Site",
    y = "R² (Coefficient of Determination)",
    fill = "Analysis Type"
  ) +
  ylim(0, max(comparison_long$R_squared) * 1.15) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

ggsave("h2_alternative_figures/Fig3_R2_Comparison.png",
       fig3, width = 8, height = 6, dpi = 300)

# FIGURE 4: Spearman ρ Comparison
cat("  Creating Figure 4: Spearman ρ Comparison...\n")

comparison_spearman <- comparison_table %>%
  tidyr::pivot_longer(
    cols = c(Alt_Spearman, Orig_Spearman),
    names_to = "Analysis",
    values_to = "Spearman_rho"
  ) %>%
  dplyr::mutate(
    Analysis = dplyr::recode(Analysis,
                             "Alt_Spearman" = "Counts vs Dry Mass",
                             "Orig_Spearman" = "Biomass vs Biomass")
  )

fig4 <- ggplot(comparison_spearman, aes(x = Site, y = Spearman_rho, fill = Analysis)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Spearman_rho, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Counts vs Dry Mass" = "#FF7F00",
                               "Biomass vs Biomass" = "#984EA3")) +
  labs(
    title = "Comparison of Spearman ρ: Alternative vs Original Analysis",
    subtitle = "Non-parametric correlation (robust to outliers)",
    x = "Site",
    y = "Spearman ρ",
    fill = "Analysis Type"
  ) +
  ylim(-0.2, 1) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

ggsave("h2_alternative_figures/Fig4_Spearman_Comparison.png",
       fig4, width = 8, height = 6, dpi = 300)

################################################################################
# SECTION 7: EXPORT RESULTS
################################################################################

cat("\nExporting results...\n")

# Export comparison table
readr::write_csv(comparison_table, 
                 "h2_alternative_figures/H2_Alternative_Comparison_Table.csv")
cat("  ✓ Exported: H2_Alternative_Comparison_Table.csv\n")

# Export full dataset with both approaches
readr::write_csv(h2_alternative,
                 "h2_alternative_figures/H2_Alternative_Full_Data.csv")
cat("  ✓ Exported: H2_Alternative_Full_Data.csv\n")

# Create detailed results summary
results_summary <- data.frame(
  Site = rep(c("Site1", "Site2", "Site3"), 2),
  Analysis = rep(c("Counts_vs_DryMass", "Biomass_vs_Biomass"), each = 3),
  N = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$n),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$n)
  ),
  R_squared = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$r_squared),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$r_squared)
  ),
  Slope = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$slope),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$slope)
  ),
  Slope_SE = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$slope_se),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$slope_se)
  ),
  P_value = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$p_value),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$p_value)
  ),
  Pearson_r = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$pearson_r),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$pearson_r)
  ),
  Spearman_rho = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$spearman_rho),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$spearman_rho)
  ),
  Robust_R2 = c(
    sapply(c("Site1", "Site2", "Site3"), function(s) results_counts_dry[[s]]$robust_r2),
    sapply(c("Site1", "Site2", "Site3"), function(s) results_biomass_biomass[[s]]$robust_r2)
  )
)

readr::write_csv(results_summary,
                 "h2_alternative_figures/H2_Alternative_Detailed_Results.csv")
cat("  ✓ Exported: H2_Alternative_Detailed_Results.csv\n")

################################################################################
# SECTION 8: FINAL SUMMARY AND RECOMMENDATIONS
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("FINAL SUMMARY AND CONCLUSIONS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("KEY FINDINGS:\n\n")

for (i in 1:3) {
  site <- comparison_table$Site[i]
  
  cat(site, ":\n")
  cat("  Alternative (Counts vs Dry Mass):\n")
  cat("    R² =", round(comparison_table$Alt_R2[i], 3),
      "| ρ =", round(comparison_table$Alt_Spearman[i], 3), "\n")
  cat("  Original (Biomass vs Biomass):\n")
  cat("    R² =", round(comparison_table$Orig_R2[i], 3),
      "| ρ =", round(comparison_table$Orig_Spearman[i], 3), "\n")
  cat("  Change: ΔR² =", round(comparison_table$Delta_R2[i], 3),
      "| Δρ =", round(comparison_table$Delta_Spearman[i], 3), "\n\n")
}

# Determine overall pattern
mean_delta_r2 <- mean(comparison_table$Delta_R2)
all_similar <- all(abs(comparison_table$Delta_R2) < 0.05)
any_improved <- any(comparison_table$Delta_R2 > 0.1)

if (all_similar) {
  cat("1. CONCLUSION: Measurement error is NOT the primary factor\n")
  cat("   → Both approaches yield similar correlations\n")
  cat("   → The poor Site1 correlation reflects genuine detection bias\n")
  cat("   → Your original biomass-biomass analysis was appropriate\n\n")
  
  cat("   'Supplementary analysis using raw counts and dry mass measurements\n")
  cat("    (eliminating biomass inference) yielded similar correlation patterns,\n")
  cat("    confirming that measurement error does not primarily explain the\n")
  cat("    observed site-specific variation in FAIRD-AMMOD agreement.'\n\n")
  
} else if (any_improved) {
  cat("1. CONCLUSION: Measurement error DOES contribute to variation\n")
  cat("   → Alternative approach shows improved correlations\n")
  cat("   → Biomass inference adds noise to the relationship\n")
  cat("   → Consider this in interpretation of Site1 results\n\n")
  
  cat("2. RECOMMENDATION: Include as supplementary analysis\n")
  cat("   'Alternative validation using detection counts versus dry mass\n")
  cat("    revealed stronger correlations (ΔR² = +", 
      round(mean_delta_r2, 2), " on average),\n")
  cat("    suggesting that biomass inference contributes to observed variation\n")
  cat("    in FAIRD-AMMOD agreement, though detection bias remains evident.'\n\n")
  
} else {
  cat("1. CONCLUSION: Mixed results - investigate further\n")
  cat("   → Pattern varies by site\n")
  cat("   → May indicate complex interaction of factors\n\n")
  
  cat("2. RECOMMENDATION: Report both approaches\n")
  cat("   Compare and discuss possible explanations\n\n")
}

cat(rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Output files created in: h2_alternative_figures/\n")
cat("  ✓ 4 PNG figures\n")
cat("  ✓ 3 CSV tables\n\n")

################################################################################
# END OF SCRIPT
################################################################################

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")