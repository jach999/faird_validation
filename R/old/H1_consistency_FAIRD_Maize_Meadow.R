# ==============================================================================
# QUICK ANALYSIS: FAIRD MAIZE vs FAIRD MEADOW CONSISTENCY
# ==============================================================================
# Purpose: Check if temporal patterns are consistent across habitats
# Comparison: Average(F1+F2) vs Average(F3+F4) over time
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

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

# Load data
source("load_data.R")
load("data_cleaned.RData")

etraps_daily <- data_list$etraps_daily

cat("\n=== FAIRD: MAIZE vs MEADOW CONSISTENCY ===\n\n")

# --- Calculate daily averages per habitat ---
faird_habitat <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2", "FAIRD3", "FAIRD4")) %>%
  dplyr::mutate(
    Habitat = case_when(
      Device %in% c("FAIRD1", "FAIRD2") ~ "Maize",
      Device %in% c("FAIRD3", "FAIRD4") ~ "Meadow"
    )
  ) %>%
  dplyr::group_by(Date, Habitat) %>%
  dplyr::summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Habitat,
    values_from = mean_abundance,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    Maize_log = log10(Maize + 1),
    Meadow_log = log10(Meadow + 1)
  )

cat("Data prepared:", nrow(faird_habitat), "days\n")
cat("Date range:", min(faird_habitat$Date), "to", max(faird_habitat$Date), "\n\n")

# --- Linear regression (log-transformed) ---
model_log <- lm(Meadow_log ~ Maize_log, data = faird_habitat)
model_summary <- summary(model_log)

cat("=== REGRESSION RESULTS (LOG-TRANSFORMED) ===\n")
cat("R² =", round(model_summary$r.squared, 3), "\n")
cat("Adjusted R² =", round(model_summary$adj.r.squared, 3), "\n")
cat("Slope =", round(coef(model_log)[2], 3), "\n")
cat("p-value =", format.pval(model_summary$coefficients[2, 4], digits = 3), "\n\n")

# --- Correlations ---
cor_pearson <- cor.test(~ Maize_log + Meadow_log, 
                        data = faird_habitat, 
                        method = "pearson")

cor_spearman <- cor.test(~ Maize + Meadow, 
                         data = faird_habitat, 
                         method = "spearman")

cat("Pearson r (log) =", round(cor_pearson$estimate, 3), 
    ", p =", format.pval(cor_pearson$p.value, digits = 3), "\n")
cat("Spearman ρ (original) =", round(cor_spearman$estimate, 3), 
    ", p =", format.pval(cor_spearman$p.value, digits = 3), "\n\n")

# --- Plot ---
# Load custom colors if available
if(file.exists("custom_colors.R")) {
  source("custom_colors.R")
  faird_color <- custom_colors["FAIRD"]
} else {
  faird_color <- "forestgreen"  # green
  }

plot_habitat <- ggplot(faird_habitat, aes(x = Maize_log, y = Meadow_log)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = faird_color,
    alpha = 0.3,
    linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  stat_cor(method = "pearson",aes(label = paste(after_stat(rr.label),after_stat(p.label),sep = "~`,`~")), 
           label.x = 0.1, label.y = 1.8, 
           size = 4) +
  labs(
    title = "FAIRD: Maize vs Meadow Temporal Consistency",
    subtitle = "Daily averages: (F1+F2)/2 vs (F3+F4)/2",
    x = "log₁₀(FAIRD Maize + 1)",
    y = "log₁₀(FAIRD Meadow + 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1)

print(plot_habitat)

# --- Save ---
ggsave("outputs/H1_FAIRD_Maize_vs_Meadow_Consistency.png", 
       plot = plot_habitat, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Maize_vs_Meadow_Consistency.pdf",
       plot = plot_habitat,
       width = 7, height = 6, device = cairo_pdf)

cat("✅ Plot saved:\n")
cat("   - outputs/H1_FAIRD_Maize_vs_Meadow_Consistency.png\n")
cat("   - outputs/H1_FAIRD_Maize_vs_Meadow_Consistency.pdf\n\n")

# --- Summary interpretation ---
if(model_summary$r.squared > 0.5) {
  cat("✅ INTERPRETATION: Strong cross-habitat consistency (R² > 0.5)\n")
  cat("   → FAIRD captures similar temporal patterns across landscapes\n")
} else if(model_summary$r.squared > 0.3) {
  cat("⚠️  INTERPRETATION: Moderate cross-habitat consistency (R² = 0.3-0.5)\n")
  cat("   → Habitat-specific dynamics present but some shared patterns\n")
} else {
  cat("❌ INTERPRETATION: Low cross-habitat consistency (R² < 0.3)\n")
  cat("   → Habitats show independent temporal dynamics\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")


# ==============================================================================
# TEMPORAL PATTERN ANALYSIS
# ==============================================================================

cat("\n=== TEMPORAL PATTERN ANALYSIS ===\n\n")

# Calculate day-by-day differences
temporal_analysis <- faird_habitat %>%
  dplyr::mutate(
    Difference_log = abs(Maize_log - Meadow_log),
    Direction = case_when(
      Maize_log > Meadow_log ~ "Maize > Meadow",
      Maize_log < Meadow_log ~ "Meadow > Maize",
      TRUE ~ "Equal"
    )
  )

# Summary statistics
mean_diff <- mean(temporal_analysis$Difference_log, na.rm = TRUE)
sd_diff <- sd(temporal_analysis$Difference_log, na.rm = TRUE)
max_diff <- max(temporal_analysis$Difference_log, na.rm = TRUE)

cat("Day-to-day divergence (log-scale):\n")
cat(sprintf("  Mean absolute difference: %.3f\n", mean_diff))
cat(sprintf("  SD of differences: %.3f\n", sd_diff))
cat(sprintf("  Maximum difference: %.3f\n\n", max_diff))

# Days where one habitat has much higher activity
threshold_diff <- 0.5  # log-scale difference threshold
divergent_days <- temporal_analysis %>%
  dplyr::filter(Difference_log > threshold_diff) %>%
  dplyr::arrange(desc(Difference_log))

if(nrow(divergent_days) > 0) {
  cat(sprintf("Days with large divergence (>%.1f log units): %d/%d days\n", 
              threshold_diff, nrow(divergent_days), nrow(temporal_analysis)))
  cat("\nTop divergent days:\n")
  print(divergent_days %>% 
          dplyr::select(Date, Maize, Meadow, Difference_log, Direction) %>%
          head(5))
} else {
  cat("✅ No days with large divergence - strong concordance throughout\n")
}

cat("\n=== INTERPRETATION ===\n")
if(mean_diff < 0.3) {
  cat("✅ Excellent temporal concordance: Mean difference < 0.3 log units\n")
  cat("   Both habitats show highly synchronized daily dynamics\n")
} else if(mean_diff < 0.5) {
  cat("✓ Good temporal concordance: Mean difference < 0.5 log units\n")
  cat("   General synchronization with some habitat-specific variation\n")
} else {
  cat("⚠️  Moderate temporal concordance: Mean difference > 0.5 log units\n")
  cat("   Substantial habitat-specific dynamics present\n")
}


# ==============================================================================
# TEMPORAL CONCORDANCE ANALYSIS (DIRECTIONAL SYNCHRONIZATION)
# ==============================================================================

cat("\n=== TEMPORAL CONCORDANCE (DIRECTIONAL SYNCHRONIZATION) ===\n\n")

# Calculate directional concordance (same as H2_validation.R)
temporal_concordance <- faird_habitat %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    Maize_change = Maize - lag(Maize),
    Meadow_change = Meadow - lag(Meadow),
    concordant = sign(Maize_change) == sign(Meadow_change)
  ) %>%
  dplyr::filter(!is.na(concordant))

# Calculate concordance percentage
n_changes <- nrow(temporal_concordance)
n_concordant <- sum(temporal_concordance$concordant, na.rm = TRUE)
concordance_pct <- (n_concordant / n_changes) * 100

cat(sprintf("Temporal Concordance: %.1f%%\n", concordance_pct))
cat(sprintf("  Days analyzed: %d (n-1 for lag comparison)\n", n_changes))
cat(sprintf("  Concordant days: %d\n", n_concordant))
cat(sprintf("  Discordant days: %d\n\n", n_changes - n_concordant))

# Show discordant days
if((n_changes - n_concordant) > 0) {
  cat("Days with opposite directional changes:\n")
  discordant_days <- temporal_concordance %>%
    dplyr::filter(!concordant) %>%
    dplyr::select(Date, Maize, Meadow, Maize_change, Meadow_change) %>%
    dplyr::mutate(
      Direction = case_when(
        Maize_change > 0 & Meadow_change < 0 ~ "Maize ↑, Meadow ↓",
        Maize_change < 0 & Meadow_change > 0 ~ "Maize ↓, Meadow ↑",
        TRUE ~ "Other"
      )
    )
  print(discordant_days)
  cat("\n")
}

# Kendall's tau test (non-parametric correlation of changes)
kendall_test <- cor.test(temporal_concordance$Maize_change, 
                         temporal_concordance$Meadow_change, 
                         method = "kendall", 
                         exact = FALSE)

cat("Kendall's tau (correlation of directional changes):\n")
cat(sprintf("  τ = %.3f\n", kendall_test$estimate))
cat(sprintf("  p-value = %.4f\n\n", kendall_test$p.value))

# Update the final concordance summary section
cat("\n=== CONCORDANCE SUMMARY ===\n")
cat(sprintf("Overall temporal concordance (R²): %.1f%%\n", model_summary$r.squared * 100))
cat(sprintf("Directional concordance: %.1f%% (%d/%d days)\n", 
            concordance_pct, n_concordant, n_changes))
cat(sprintf("Days with strong agreement (<0.5 log units): %d/%d (%.1f%%)\n", 
            nrow(temporal_analysis) - nrow(divergent_days),
            nrow(temporal_analysis),
            ((nrow(temporal_analysis) - nrow(divergent_days)) / nrow(temporal_analysis)) * 100))


cat("\n=== TEMPORAL ANALYSIS COMPLETE ===\n")


# ==============================================================================
# SECTION: TEMPORAL CONCORDANCE PLOT
# ==============================================================================

cat("\n=== GENERATING TEMPORAL CONCORDANCE PLOT ===\n\n")

# Prepare long format data for temporal plot
faird_temporal <- faird_habitat %>%
  dplyr::select(Date, Maize, Meadow) %>%
  tidyr::pivot_longer(
    cols = c(Maize, Meadow),
    names_to = "Habitat",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(
    Abundance_log = log10(Abundance + 1)
  )

# --- Temporal plot (log-transformed) ---
plot_temporal <- ggplot(faird_temporal, 
                        aes(x = Date, y = Abundance_log, 
                            color = Habitat, shape = Habitat)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c("Maize" = 16, "Meadow" = 17)) +  # Circle and triangle
  labs(
    title = "Temporal Concordance: FAIRD Maize vs Meadow",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       summary(model_log)$r.squared,
                       concordance_pct,
                       kendall_test$estimate,
                       kendall_test$p.value),
    x = "Date",
    y = "log₁₀(Daily Abundance + 1)",
    color = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(plot_temporal)

# Save temporal plot
ggsave("outputs/H1_FAIRD_Temporal_Concordance_Maize_vs_Meadow.png", 
       plot = plot_temporal, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Temporal_Concordance_Maize_vs_Meadow.pdf",
       plot = plot_temporal,
       width = 7, height = 6, device = cairo_pdf)

cat("✅ Temporal concordance plot saved:\n")
cat("   - outputs/H1_FAIRD_Temporal_Concordance_Maize_vs_Meadow.png\n")
cat("   - outputs/H1_FAIRD_Temporal_Concordance_Maize_vs_Meadow.pdf\n\n")

# ==============================================================================
# COMBINED FIGURE (A: scatter | B: temporal concordance)
# ==============================================================================

cat("\n=== GENERATING COMBINED FIGURE (A: scatter | B: temporal) ===\n\n")

library(patchwork)

plot_habitat_A <- plot_habitat +
  labs(title = "A) FAIRD: Maize vs Meadow Temporal Consistency")

plot_temporal_B <- plot_temporal +
  labs(title = "B) Temporal Concordance: FAIRD Maize vs Meadow")

combined_cross_habitat <- plot_habitat_A + plot_temporal_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_cross_habitat)

ggsave("outputs/H1_FAIRD_Combined_Maize_vs_Meadow.png",
       plot   = combined_cross_habitat,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Combined_Maize_vs_Meadow.pdf",
       plot   = combined_cross_habitat,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved:\n")
cat("   - outputs/H1_FAIRD_Combined_Maize_vs_Meadow.png\n")
cat("   - outputs/H1_FAIRD_Combined_Maize_vs_Meadow.pdf\n\n")

# ==============================================================================
# PUBLICATION TABLE: FAIRD CROSS-HABITAT CONSISTENCY (Table 7 - FAIRD row)
# Add this block BEFORE the final sink() call at the end of the script.
# ==============================================================================

cat("\n=== PUBLICATION TABLE: FAIRD CROSS-HABITAT CONSISTENCY ===\n\n")

# --- Slope 95% CI ---
ci_faird    <- confint(model_log)
slope_ci    <- sprintf("%.2f [%.2f, %.2f]",
                       coef(model_log)[2], ci_faird[2, 1], ci_faird[2, 2])

# --- Pseudo-R² via robust regression (Huber M-estimation, log scale) ---
robust_faird  <- MASS::rlm(Meadow_log ~ Maize_log, data = faird_habitat,
                           psi = MASS::psi.huber)
ss_res_r      <- sum(robust_faird$residuals^2)
ss_tot_ols    <- sum((faird_habitat$Meadow_log - mean(faird_habitat$Meadow_log))^2)
pseudo_r2     <- 1 - (ss_res_r / ss_tot_ols)
low_weight    <- sum(robust_faird$w < 1)

# --- Significance labels ---
sig_label <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

p_r2  <- model_summary$coefficients[2, 4]
p_tau <- kendall_test$p.value
p_rho <- cor_spearman$p.value

# --- Assemble table ---
table_faird_crosshabitat <- tibble::tibble(
  Device       = "FAIRD",
  Pair         = "Maize vs Meadow",
  Metric       = "Abundance",
  n            = nrow(faird_habitat),
  R2_log       = paste0(round(model_summary$r.squared, 3),
                        " ", sig_label(p_r2)),
  Slope_95CI   = slope_ci,
  Concordance  = paste0(round(concordance_pct, 1), "%"),
  Spearman_rho = paste0(round(cor_spearman$estimate, 3),
                        " ", sig_label(p_rho)),
  Kendall_tau  = paste0(round(kendall_test$estimate, 3),
                        " ", sig_label(p_tau)),
  Pseudo_R2_log    = round(pseudo_r2, 3),
  Low_weight_days  = low_weight
)

# --- Print to console ---
cat("Table 7 (FAIRD row) — ready for publication:\n\n")
print(table_faird_crosshabitat, width = Inf)

# --- Export CSV ---
csv_path <- file.path("outputs",
                      paste0("H1_FAIRD_crosshabitat_table7_",
                             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))

readr::write_csv(table_faird_crosshabitat, csv_path)
cat("\n✅ Table exported to:", csv_path, "\n")

# --- Robust regression summary to console ---
cat("\nRobust regression validation:\n")
cat(sprintf("  FAIRD Maize vs Meadow: Pseudo-R² = %.3f | Low-weight days: %d/%d\n",
            pseudo_r2, low_weight, nrow(faird_habitat)))

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")