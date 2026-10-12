# ==============================================================================
# QUICK ANALYSIS: AMMOD MAIZE vs MEADOW CONSISTENCY
# ==============================================================================
# Purpose: Compare AMMOD cross-habitat consistency with FAIRD results
# Comparisons: 
#   1. AMMOD1 (Maize Site1) vs AMMOD3 (Meadow Site3)
#   2. AMMOD2 (Maize Site2) vs AMMOD3 (Meadow Site3)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

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
load("data_cleaned.RData")

# Extract AMMOD biomass data
ammod_biomass <- data_list$ammod_biomass

cat("\n=== AMMOD: MAIZE vs MEADOW CONSISTENCY ===\n\n")

# --- Prepare data ---
ammod_wide <- ammod_biomass %>%
  dplyr::select(Date, Device, Live_mass) %>%
  tidyr::pivot_wider(
    names_from = Device,
    values_from = Live_mass,
    values_fill = NA  # Keep NAs explicit for AMMOD4
  ) %>%
  dplyr::filter(!is.na(AMMOD1) & !is.na(AMMOD2) & !is.na(AMMOD3)) %>%
  dplyr::mutate(
    AMMOD1_log = log10(AMMOD1 + 1),
    AMMOD2_log = log10(AMMOD2 + 1),
    AMMOD3_log = log10(AMMOD3 + 1)
  )

cat("Data prepared:", nrow(ammod_wide), "days with complete data\n")
cat("Date range:", min(ammod_wide$Date), "to", max(ammod_wide$Date), "\n\n")

# ==============================================================================
# ANALYSIS 1: AMMOD1 (Maize Site1) vs AMMOD3 (Meadow Site3)
# ==============================================================================

cat("=== COMPARISON 1: AMMOD1 (Maize Site1) vs AMMOD3 (Meadow) ===\n")

# Linear regression
model1 <- lm(AMMOD3_log ~ AMMOD1_log, data = ammod_wide)
summary1 <- summary(model1)

cat("R² =", round(summary1$r.squared, 3), "\n")
cat("Adjusted R² =", round(summary1$adj.r.squared, 3), "\n")
cat("Slope =", round(coef(model1)[2], 3), "\n")
cat("p-value =", format.pval(summary1$coefficients[2, 4], digits = 3), "\n\n")

# Correlations
cor_pearson1 <- cor.test(~ AMMOD1_log + AMMOD3_log, 
                         data = ammod_wide, 
                         method = "pearson")

cor_spearman1 <- cor.test(~ AMMOD1 + AMMOD3, 
                          data = ammod_wide, 
                          method = "spearman")

cat("Pearson r (log) =", round(cor_pearson1$estimate, 3), 
    ", p =", format.pval(cor_pearson1$p.value, digits = 3), "\n")
cat("Spearman ρ (original) =", round(cor_spearman1$estimate, 3), 
    ", p =", format.pval(cor_spearman1$p.value, digits = 3), "\n\n")

# ==============================================================================
# ANALYSIS 2: AMMOD2 (Maize Site2) vs AMMOD3 (Meadow Site3)
# ==============================================================================

cat("=== COMPARISON 2: AMMOD2 (Maize Site2) vs AMMOD3 (Meadow) ===\n")

# Linear regression
model2 <- lm(AMMOD3_log ~ AMMOD2_log, data = ammod_wide)
summary2 <- summary(model2)

cat("R² =", round(summary2$r.squared, 3), "\n")
cat("Adjusted R² =", round(summary2$adj.r.squared, 3), "\n")
cat("Slope =", round(coef(model2)[2], 3), "\n")
cat("p-value =", format.pval(summary2$coefficients[2, 4], digits = 3), "\n\n")

# Correlations
cor_pearson2 <- cor.test(~ AMMOD2_log + AMMOD3_log, 
                         data = ammod_wide, 
                         method = "pearson")

cor_spearman2 <- cor.test(~ AMMOD2 + AMMOD3, 
                          data = ammod_wide, 
                          method = "spearman")

cat("Pearson r (log) =", round(cor_pearson2$estimate, 3), 
    ", p =", format.pval(cor_pearson2$p.value, digits = 3), "\n")
cat("Spearman ρ (original) =", round(cor_spearman2$estimate, 3), 
    ", p =", format.pval(cor_spearman2$p.value, digits = 3), "\n\n")

# ==============================================================================
# PLOTS
# ==============================================================================

# Load custom colors if available
if(file.exists("custom_colors.R")) {
  source("custom_colors.R")
  ammod_color <- custom_colors["AMMOD"]
} else {
  ammod_color <- "#8B0000"  # Dark red
}

# --- Determine shared axis limits ---
all_ammod_values <- c(ammod_wide$AMMOD1_log, 
                      ammod_wide$AMMOD2_log, 
                      ammod_wide$AMMOD3_log)
max_val <- max(all_ammod_values, na.rm = TRUE)
min_val <- min(all_ammod_values, na.rm = TRUE)
buffer <- (max_val - min_val) * 0.05
shared_limits <- c(min_val - buffer, max_val + buffer)

cat("Shared axis limits:", round(shared_limits, 3), "\n\n")

# --- Plot 1: AMMOD1 vs AMMOD3 ---
plot_consistency_1_3 <- ggplot(ammod_wide, aes(x = AMMOD1_log, y = AMMOD3_log)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = ammod_color,
    alpha = 0.3,
    linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  stat_cor(method = "pearson",
           aes(label = paste(after_stat(rr.label),
                             after_stat(p.label),
                             sep = "~`,`~")), 
           label.x = shared_limits[1] + 0.1, 
           label.y = shared_limits[2] - 0.2, 
           size = 3.5) +
  labs(
    title = "A) AMMOD1 (Maize Site1) vs AMMOD3 (Meadow)",
    x = "log₁₀(AMMOD1 Live Mass + 1)",
    y = "log₁₀(AMMOD3 Live Mass + 1)"
  ) +
  coord_cartesian(xlim = shared_limits, ylim = shared_limits) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )

# --- Plot 2: AMMOD2 vs AMMOD3 ---
plot_consistency_2_3 <- ggplot(ammod_wide, aes(x = AMMOD2_log, y = AMMOD3_log)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = ammod_color,
    alpha = 0.3,
    linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  stat_cor(method = "pearson",
           aes(label = paste(after_stat(rr.label),
                             after_stat(p.label),
                             sep = "~`,`~")), 
           label.x = shared_limits[1] + 0.1, 
           label.y = shared_limits[2] - 0.2, 
           size = 3.5) +
  labs(
    title = "B) AMMOD2 (Maize Site2) vs AMMOD3 (Meadow)",
    x = "log₁₀(AMMOD2 Live Mass + 1)",
    y = "log₁₀(AMMOD3 Live Mass + 1)"
  ) +
  coord_cartesian(xlim = shared_limits, ylim = shared_limits) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )

# --- Combined plot ---
#combined <- plot_consistency_1_3 + plot_consistency_2_3 + plot_layout(ncol = 2)

print(plot_consistency_1_3)
print(plot_consistency_2_3)

# --- Save ---
ggsave("outputs/H1_ammod_Maize1_vs_Meadow_Consistency.png", 
       plot = plot_consistency_1_3, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ammod_Maize1_vs_Meadow_Consistency.pdf",
       plot = plot_consistency_1_3,
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ammod_Maize2_vs_Meadow_Consistency.png", 
       plot = plot_consistency_2_3, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ammod_Maize2_vs_Meadow_Consistency.pdf",
       plot = plot_consistency_2_3,
       width = 7, height = 6, dpi = 300, bg = "white")

cat("✅ Plots saved:\n")
cat("   - outputs/H1_ammod_Maize1_vs_Meadow_Consistency.png\n")
cat("   - outputs/H1_ammod_Maize1_vs_Meadow_Consistency.pdf\n\n")
cat("   - outputs/H1_ammod_Maize2_vs_Meadow_Consistency.png\n")
cat("   - outputs/H1_ammod_Maize2_vs_Meadow_Consistency.pdf\n\n")

# ==============================================================================
# TEMPORAL CONCORDANCE ANALYSIS — AMMOD MAIZE vs MEADOW
# Pairs: AMMOD1 vs AMMOD3 | AMMOD2 vs AMMOD3
# ==============================================================================

cat("\n=== TEMPORAL CONCORDANCE (AMMOD MAIZE vs MEADOW) ===\n\n")

# ------------------------------------------------------------------------------
# PART A: AMMOD1 (Maize Site1) vs AMMOD3 (Meadow Site3)
# ------------------------------------------------------------------------------
cat("--- AMMOD1 vs AMMOD3 ---\n\n")

temporal_concordance_mm1 <- ammod_wide %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    AMMOD1_change_mm = AMMOD1 - dplyr::lag(AMMOD1),
    AMMOD3_change_mm = AMMOD3 - dplyr::lag(AMMOD3),
    concordant_mm1   = sign(AMMOD1_change_mm) == sign(AMMOD3_change_mm)
  ) %>%
  dplyr::filter(!is.na(concordant_mm1))

n_changes_mm1      <- nrow(temporal_concordance_mm1)
n_concordant_mm1   <- sum(temporal_concordance_mm1$concordant_mm1, na.rm = TRUE)
concordance_pct_mm1 <- (n_concordant_mm1 / n_changes_mm1) * 100

cat(sprintf("Temporal Concordance (AMMOD1 vs AMMOD3): %.1f%%\n", concordance_pct_mm1))
cat(sprintf("  Days analyzed: %d (n-1 lag)\n", n_changes_mm1))
cat(sprintf("  Concordant days: %d | Discordant: %d\n\n",
            n_concordant_mm1, n_changes_mm1 - n_concordant_mm1))

if ((n_changes_mm1 - n_concordant_mm1) > 0) {
  cat("Days with opposite directional changes (AMMOD1 vs AMMOD3):\n")
  discordant_mm1 <- temporal_concordance_mm1 %>%
    dplyr::filter(!concordant_mm1) %>%
    dplyr::select(Date, AMMOD1, AMMOD3, AMMOD1_change_mm, AMMOD3_change_mm) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        AMMOD1_change_mm > 0 & AMMOD3_change_mm < 0 ~ "AMMOD1 up, AMMOD3 down",
        AMMOD1_change_mm < 0 & AMMOD3_change_mm > 0 ~ "AMMOD1 down, AMMOD3 up",
        TRUE ~ "Other"
      )
    )
  print(discordant_mm1)
  cat("\n")
}

kendall_mm1 <- cor.test(temporal_concordance_mm1$AMMOD1_change_mm,
                        temporal_concordance_mm1$AMMOD3_change_mm,
                        method = "kendall",
                        exact  = FALSE)

cat(sprintf("Kendall's tau (AMMOD1 vs AMMOD3): tau = %.3f | p = %.4f\n\n",
            kendall_mm1$estimate, kendall_mm1$p.value))

# ------------------------------------------------------------------------------
# PART B: AMMOD2 (Maize Site2) vs AMMOD3 (Meadow Site3)
# ------------------------------------------------------------------------------
cat("--- AMMOD2 vs AMMOD3 ---\n\n")

temporal_concordance_mm2 <- ammod_wide %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    AMMOD2_change_mm = AMMOD2 - dplyr::lag(AMMOD2),
    AMMOD3_change_mm = AMMOD3 - dplyr::lag(AMMOD3),
    concordant_mm2   = sign(AMMOD2_change_mm) == sign(AMMOD3_change_mm)
  ) %>%
  dplyr::filter(!is.na(concordant_mm2))

n_changes_mm2       <- nrow(temporal_concordance_mm2)
n_concordant_mm2    <- sum(temporal_concordance_mm2$concordant_mm2, na.rm = TRUE)
concordance_pct_mm2 <- (n_concordant_mm2 / n_changes_mm2) * 100

cat(sprintf("Temporal Concordance (AMMOD2 vs AMMOD3): %.1f%%\n", concordance_pct_mm2))
cat(sprintf("  Days analyzed: %d (n-1 lag)\n", n_changes_mm2))
cat(sprintf("  Concordant days: %d | Discordant: %d\n\n",
            n_concordant_mm2, n_changes_mm2 - n_concordant_mm2))

if ((n_changes_mm2 - n_concordant_mm2) > 0) {
  cat("Days with opposite directional changes (AMMOD2 vs AMMOD3):\n")
  discordant_mm2 <- temporal_concordance_mm2 %>%
    dplyr::filter(!concordant_mm2) %>%
    dplyr::select(Date, AMMOD2, AMMOD3, AMMOD2_change_mm, AMMOD3_change_mm) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        AMMOD2_change_mm > 0 & AMMOD3_change_mm < 0 ~ "AMMOD2 up, AMMOD3 down",
        AMMOD2_change_mm < 0 & AMMOD3_change_mm > 0 ~ "AMMOD2 down, AMMOD3 up",
        TRUE ~ "Other"
      )
    )
  print(discordant_mm2)
  cat("\n")
}

kendall_mm2 <- cor.test(temporal_concordance_mm2$AMMOD2_change_mm,
                        temporal_concordance_mm2$AMMOD3_change_mm,
                        method = "kendall",
                        exact  = FALSE)

cat(sprintf("Kendall's tau (AMMOD2 vs AMMOD3): tau = %.3f | p = %.4f\n\n",
            kendall_mm2$estimate, kendall_mm2$p.value))

# ------------------------------------------------------------------------------
# TEMPORAL CONCORDANCE PLOTS
# ------------------------------------------------------------------------------

# --- Plot AMMOD1 vs AMMOD3 ---
ammod_temporal_mm1 <- ammod_wide %>%
  dplyr::select(Date, AMMOD1, AMMOD3) %>%
  tidyr::pivot_longer(
    cols      = c(AMMOD1, AMMOD3),
    names_to  = "Device",
    values_to = "Biomass"
  ) %>%
  dplyr::mutate(Biomass_log = log10(Biomass + 1))

plot_temporal_mm1 <- ggplot(
  ammod_temporal_mm1,
  aes(x = Date, y = Biomass_log, color = Device, shape = Device)
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
    labels = c("AMMOD1" = "AMMOD1 (Site1 - Maize)",
               "AMMOD3" = "AMMOD3 (Site3 - Meadow)")
  ) +
  scale_shape_manual(
    values = c("AMMOD1" = 16, "AMMOD3" = 17),
    labels = c("AMMOD1" = "AMMOD1 (Site1 - Maize)",
               "AMMOD3" = "AMMOD3 (Site3 - Meadow)")
  ) +
  labs(
    title    = "Temporal Concordance: AMMOD1 vs AMMOD3",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       summary1$r.squared,
                       concordance_pct_mm1,
                       kendall_mm1$estimate,
                       kendall_mm1$p.value),
    x        = "Date",
    y        = "log\u2081\u2080(Live Biomass + 1)",
    color    = "Device",
    shape    = "Device"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(plot_temporal_mm1)

ggsave("outputs/H1_AMMOD_MM_Temporal_Concordance_S1_S3.png",
       plot   = plot_temporal_mm1,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_MM_Temporal_Concordance_S1_S3.pdf",
       plot   = plot_temporal_mm1,
       width  = 7, height = 6, device = cairo_pdf)

cat("Temporal concordance plot saved (AMMOD1 vs AMMOD3):\n")
cat("   - outputs/H1_AMMOD_MM_Temporal_Concordance_S1_S3.png\n")
cat("   - outputs/H1_AMMOD_MM_Temporal_Concordance_S1_S3.pdf\n\n")

# --- Plot AMMOD2 vs AMMOD3 ---
ammod_temporal_mm2 <- ammod_wide %>%
  dplyr::select(Date, AMMOD2, AMMOD3) %>%
  tidyr::pivot_longer(
    cols      = c(AMMOD2, AMMOD3),
    names_to  = "Device",
    values_to = "Biomass"
  ) %>%
  dplyr::mutate(Biomass_log = log10(Biomass + 1))

plot_temporal_mm2 <- ggplot(
  ammod_temporal_mm2,
  aes(x = Date, y = Biomass_log, color = Device, shape = Device)
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
    labels = c("AMMOD2" = "AMMOD2 (Site2 - Maize)",
               "AMMOD3" = "AMMOD3 (Site3 - Meadow)")
  ) +
  scale_shape_manual(
    values = c("AMMOD2" = 16, "AMMOD3" = 17),
    labels = c("AMMOD2" = "AMMOD2 (Site2 - Maize)",
               "AMMOD3" = "AMMOD3 (Site3 - Meadow)")
  ) +
  labs(
    title    = "Temporal Concordance: AMMOD2 vs AMMOD3",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       summary2$r.squared,
                       concordance_pct_mm2,
                       kendall_mm2$estimate,
                       kendall_mm2$p.value),
    x        = "Date",
    y        = "log\u2081\u2080(Live Biomass + 1)",
    color    = "Device",
    shape    = "Device"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(plot_temporal_mm2)

ggsave("outputs/H1_AMMOD_MM_Temporal_Concordance_S2_S3.png",
       plot   = plot_temporal_mm2,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_MM_Temporal_Concordance_S2_S3.pdf",
       plot   = plot_temporal_mm2,
       width  = 7, height = 6, device = cairo_pdf)

cat("Temporal concordance plot saved (AMMOD2 vs AMMOD3):\n")
cat("   - outputs/H1_AMMOD_MM_Temporal_Concordance_S2_S3.png\n")
cat("   - outputs/H1_AMMOD_MM_Temporal_Concordance_S2_S3.pdf\n\n")

cat("=== TEMPORAL CONCORDANCE COMPLETE ===\n\n")

# ==============================================================================
# COMBINED FIGURES (A: scatter | B: temporal concordance)
# Two versions: AMMOD1 vs AMMOD3 | AMMOD2 vs AMMOD3
# ==============================================================================

cat("\n=== GENERATING COMBINED FIGURES (A: scatter | B: temporal) ===\n\n")

# --- Combined: AMMOD1 vs AMMOD3 ---
plot_consistency_1_3_A <- plot_consistency_1_3 +
  labs(title = "A) AMMOD1 (Maize Site1) vs AMMOD3 (Meadow)")

plot_temporal_mm1_B <- plot_temporal_mm1 +
  labs(title = "B) Temporal Concordance: AMMOD1 vs AMMOD3")

combined_ammod_mm1 <- plot_consistency_1_3_A + plot_temporal_mm1_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_ammod_mm1)

ggsave("outputs/H1_AMMOD_Combined_S1_S3.png",
       plot   = combined_ammod_mm1,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_Combined_S1_S3.pdf",
       plot   = combined_ammod_mm1,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved (AMMOD1 vs AMMOD3):\n")
cat("   - outputs/H1_AMMOD_Combined_S1_S3.png\n")
cat("   - outputs/H1_AMMOD_Combined_S1_S3.pdf\n\n")

# --- Combined: AMMOD2 vs AMMOD3 ---
plot_consistency_2_3_A <- plot_consistency_2_3 +
  labs(title = "A) AMMOD2 (Maize Site2) vs AMMOD3 (Meadow)")

plot_temporal_mm2_B <- plot_temporal_mm2 +
  labs(title = "B) Temporal Concordance: AMMOD2 vs AMMOD3")

combined_ammod_mm2 <- plot_consistency_2_3_A + plot_temporal_mm2_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_ammod_mm2)

ggsave("outputs/H1_AMMOD_Combined_S2_S3.png",
       plot   = combined_ammod_mm2,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_Combined_S2_S3.pdf",
       plot   = combined_ammod_mm2,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved (AMMOD2 vs AMMOD3):\n")
cat("   - outputs/H1_AMMOD_Combined_S2_S3.png\n")
cat("   - outputs/H1_AMMOD_Combined_S2_S3.pdf\n\n")

# ==============================================================================
# SUMMARY COMPARISON: AMMOD vs FAIRD
# ==============================================================================

cat("\n=== CROSS-HABITAT CONSISTENCY COMPARISON ===\n\n")

comparison_table <- tibble::tibble(
  Device_System = c("FAIRD", "AMMOD (Site1-Site3)", "AMMOD (Site2-Site3)"),
  Comparison = c("Maize vs Meadow", "Site1 vs Site3", "Site2 vs Site3"),
  R2 = c(
    0.72,  # From previous FAIRD analysis (you'll need to update this)
    summary1$r.squared,
    summary2$r.squared
  ),
  p_value = c(
    6.3e-07,  # From previous FAIRD analysis (you'll need to update this)
    summary1$coefficients[2, 4],
    summary2$coefficients[2, 4]
  ),
  Interpretation = case_when(
    R2 > 0.6 ~ "Strong synchronization",
    R2 > 0.4 ~ "Moderate synchronization",
    R2 > 0.2 ~ "Weak synchronization",
    TRUE ~ "Poor synchronization"
  )
) %>%
  dplyr::mutate(
    R2 = round(R2, 3),
    p_value = format.pval(p_value, digits = 3)
  )

print(comparison_table)

cat("\n=== KEY FINDINGS ===\n")

if(summary1$r.squared > 0.5 & summary2$r.squared > 0.5) {
  cat("✅ AMMOD shows strong cross-habitat consistency (both R² > 0.5)\n")
  cat("   → Similar to FAIRD, regional drivers dominate temporal patterns\n")
} else if(summary1$r.squared < 0.3 | summary2$r.squared < 0.3) {
  cat("⚠️  AMMOD shows weaker cross-habitat consistency than FAIRD\n")
  cat("   → AMMOD may be more sensitive to local habitat composition\n")
  cat("   → Reflects broader taxonomic sampling of small Diptera\n")
} else {
  cat("📊 AMMOD shows moderate cross-habitat consistency\n")
  cat("   → Intermediate between FAIRD and complete independence\n")
}

# ==============================================================================
# PUBLICATION TABLE: AMMOD CROSS-HABITAT CONSISTENCY (Table 7 - AMMOD rows)
# Add this block BEFORE the final sink() call, replacing the
# "=== ANALYSIS COMPLETE ===" section at the end of the script.
# ==============================================================================

cat("\n=== PUBLICATION TABLE: AMMOD CROSS-HABITAT CONSISTENCY ===\n\n")

# --- Slope 95% CIs ---
ci1 <- confint(model1)
ci2 <- confint(model2)

slope_ci1 <- sprintf("%.2f [%.2f, %.2f]",
                     coef(model1)[2], ci1[2, 1], ci1[2, 2])
slope_ci2 <- sprintf("%.2f [%.2f, %.2f]",
                     coef(model2)[2], ci2[2, 1], ci2[2, 2])

# --- Pseudo-R² via robust regression (Huber M-estimation, log scale) ---
robust1     <- MASS::rlm(AMMOD3_log ~ AMMOD1_log, data = ammod_wide,
                         psi = MASS::psi.huber)
ss_res_r1   <- sum(robust1$residuals^2)
ss_tot_ols1 <- sum((ammod_wide$AMMOD3_log - mean(ammod_wide$AMMOD3_log))^2)
pseudo_r2_1 <- 1 - (ss_res_r1 / ss_tot_ols1)
low_weight1 <- sum(robust1$w < 1)

robust2     <- MASS::rlm(AMMOD3_log ~ AMMOD2_log, data = ammod_wide,
                         psi = MASS::psi.huber)
ss_res_r2   <- sum(robust2$residuals^2)
ss_tot_ols2 <- sum((ammod_wide$AMMOD3_log - mean(ammod_wide$AMMOD3_log))^2)
pseudo_r2_2 <- 1 - (ss_res_r2 / ss_tot_ols2)
low_weight2 <- sum(robust2$w < 1)

# --- Significance labels ---
sig_label <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

p1_r2  <- summary1$coefficients[2, 4]
p2_r2  <- summary2$coefficients[2, 4]
p1_tau <- kendall_mm1$p.value
p2_tau <- kendall_mm2$p.value
p1_rho <- cor_spearman1$p.value
p2_rho <- cor_spearman2$p.value

# --- Assemble table ---
table_ammod_crosshabitat <- tibble::tibble(
  Device      = "AMMOD",
  Pair        = c("AMMOD1 vs AMMOD3", "AMMOD2 vs AMMOD3"),
  Habitat     = c("Maize vs Meadow", "Maize vs Meadow"),
  n           = nrow(ammod_wide),
  R2_log      = paste0(round(c(summary1$r.squared,
                               summary2$r.squared), 3),
                       " ", c(sig_label(p1_r2), sig_label(p2_r2))),
  Slope_95CI  = c(slope_ci1, slope_ci2),
  Concordance = paste0(round(c(concordance_pct_mm1,
                               concordance_pct_mm2), 1), "%"),
  Spearman_rho = paste0(round(c(cor_spearman1$estimate,
                                cor_spearman2$estimate), 3),
                        " ", c(sig_label(p1_rho), sig_label(p2_rho))),
  Kendall_tau  = paste0(round(c(kendall_mm1$estimate,
                                kendall_mm2$estimate), 3),
                        " ", c(sig_label(p1_tau), sig_label(p2_tau))),
  Pseudo_R2_log     = round(c(pseudo_r2_1, pseudo_r2_2), 3),
  Low_weight_days   = c(low_weight1, low_weight2)
)

# --- Print to console ---
cat("Table 7 (AMMOD rows) — ready for publication:\n\n")
print(table_ammod_crosshabitat, width = Inf)

# --- Export CSV ---
csv_path <- file.path("outputs",
                      paste0("H1_AMMOD_crosshabitat_table7_",
                             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))

readr::write_csv(table_ammod_crosshabitat, csv_path)
cat("\n✅ Table exported to:", csv_path, "\n")

# --- Robust regression summary to console ---
cat("\nRobust regression validation:\n")
cat(sprintf("  AMMOD1 vs AMMOD3: Pseudo-R² = %.3f | Low-weight days: %d/%d\n",
            pseudo_r2_1, low_weight1, nrow(ammod_wide)))
cat(sprintf("  AMMOD2 vs AMMOD3: Pseudo-R² = %.3f | Low-weight days: %d/%d\n",
            pseudo_r2_2, low_weight2, nrow(ammod_wide)))

cat("\n=== ANALYSIS COMPLETE ===\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")