# ==============================================================================
# CORRELATION PLOTS: ORIGINAL vs LOG-TRANSFORMED
# ==============================================================================
# Publication-quality figures showing both scales for transparency
# ==============================================================================

library(tidyverse)
library(patchwork)  # For combining plots - install if needed: install.packages("patchwork")
library(scales)

load("data_cleaned.RData")
etraps_daily <- data_list$etraps_daily

cat("\n=== GENERATING CORRELATION PLOTS ===\n\n")

# ==============================================================================
# PREPARE DATA
# ==============================================================================

counts_data <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0) %>%
  dplyr::mutate(
    FAIRD1_log = log10(FAIRD1 + 1),
    FAIRD2_log = log10(FAIRD2 + 1)
  )

# Identify influential points (from diagnostics)
cooks_d <- cooks.distance(lm(FAIRD2 ~ FAIRD1, data = counts_data))
counts_data$influential <- cooks_d > 4/nrow(counts_data)

# Models
model_original <- lm(FAIRD2 ~ FAIRD1, data = counts_data)
model_log <- lm(FAIRD2_log ~ FAIRD1_log, data = counts_data)

# ==============================================================================
# PLOT A: ORIGINAL SCALE
# ==============================================================================

plot_original <- ggplot(counts_data, aes(x = FAIRD1, y = FAIRD2)) +
  # Points
  geom_point(aes(color = influential, size = influential), alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#34495E"),
                     labels = c("TRUE" = "Influential", "FALSE" = "Normal")) +
  scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3),
                    labels = c("TRUE" = "Influential", "FALSE" = "Normal")) +
  
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "#3498DB", fill = "#3498DB", 
              alpha = 0.2, linewidth = 1) +
  
  # 1:1 reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray50", linewidth = 0.8) +
  
  # Labels and theme
  labs(
    title = "A) Original Scale",
    x = "FAIRD1 Daily Abundance",
    y = "FAIRD2 Daily Abundance",
    color = "Point type",
    size = "Point type"
  ) +
  
  # Add statistics annotation
  annotate("text", x = Inf, y = -Inf, 
           label = sprintf("R² = 0.933***\nSlope = 0.77\nBP p = 0.043"),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20",
           fontface = "italic") +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  coord_fixed(ratio = 1)

# ==============================================================================
# PLOT B: LOG-TRANSFORMED SCALE
# ==============================================================================

plot_log <- ggplot(counts_data, aes(x = FAIRD1_log, y = FAIRD2_log)) +
  # Points (use same color scheme for consistency)
  geom_point(aes(color = influential, size = influential), alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#27AE60"),
                     labels = c("TRUE" = "Influential", "FALSE" = "Normal")) +
  scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3),
                    labels = c("TRUE" = "Influential", "FALSE" = "Normal")) +
  
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "#27AE60", fill = "#27AE60",
              alpha = 0.2, linewidth = 1) +
  
  # 1:1 reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.8) +
  
  # Labels and theme
  labs(
    title = "B) Log₁₀-Transformed Scale",
    x = "log₁₀(FAIRD1 + 1)",
    y = "log₁₀(FAIRD2 + 1)",
    color = "Point type",
    size = "Point type"
  ) +
  
  # Add statistics annotation
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("R² = 0.648***\nSlope = 0.84\nBP p = 0.235"),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20",
           fontface = "italic") +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  coord_fixed(ratio = 1)

# ==============================================================================
# COMBINE PLOTS
# ==============================================================================

combined_plot <- plot_original + plot_log +
  plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom")

# Print to screen
print(combined_plot)

# Save high-resolution version
ggsave("FAIRD_correlation_both_scales.png", 
       plot = combined_plot, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("FAIRD_correlation_both_scales.pdf",
       plot = combined_plot,
       width = 12, height = 6, device = cairo_pdf)

cat("\n✅ Plots saved:\n")
cat("   - FAIRD_correlation_both_scales.png (for presentations)\n")
cat("   - FAIRD_correlation_both_scales.pdf (for publication)\n\n")

# ==============================================================================
# ALTERNATIVE: SINGLE PLOT WITH DUAL AXES (more compact)
# ==============================================================================

# Create single plot with log-transformed data but showing original values as labels
plot_compact <- ggplot(counts_data, aes(x = FAIRD1_log, y = FAIRD2_log)) +
  geom_point(aes(color = influential), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#27AE60"),
                     labels = c("Influential point", "Normal observation"),
                     name = "") +
  geom_smooth(method = "lm", se = TRUE, color = "#27AE60", fill = "#27AE60",
              alpha = 0.2, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.8) +
  
  # Dual axis labels showing both log and original values
  scale_x_continuous(
    name = "FAIRD1",
    breaks = log10(c(1, 5, 10, 25, 50, 100) + 1),
    labels = c("0", "5", "10", "25", "50", "100"),
    sec.axis = sec_axis(~ ., name = "log₁₀(FAIRD1 + 1)")
  ) +
  scale_y_continuous(
    name = "FAIRD2",
    breaks = log10(c(1, 5, 10, 25, 50, 100) + 1),
    labels = c("0", "5", "10", "25", "50", "100"),
    sec.axis = sec_axis(~ ., name = "log₁₀(FAIRD2 + 1)")
  ) +
  
  labs(title = "FAIRD Replicate Consistency (Log-Transformed)") +
  
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("R² = 0.648, p < 0.001\nHomoscedastic (BP p = 0.235)\nn = 21 days"),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20") +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    axis.title.x.top = element_text(color = "gray50", size = 10),
    axis.title.y.right = element_text(color = "gray50", size = 10)
  ) +
  coord_fixed(ratio = 1)

print(plot_compact)

ggsave("FAIRD_correlation_log_compact.png",
       plot = plot_compact,
       width = 8, height = 7, dpi = 300, bg = "white")

ggsave("FAIRD_correlation_log_compact.pdf",
       plot = plot_compact,
       width = 8, height = 7, device = cairo_pdf)

cat("\n✅ Compact version also saved:\n")
cat("   - FAIRD_correlation_log_compact.png\n")
cat("   - FAIRD_correlation_log_compact.pdf\n")

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

cat("\n=== FIGURE RECOMMENDATIONS ===\n\n")
cat("Option 1: TWO-PANEL FIGURE (combined_plot)\n")
cat("  - Shows transparency: both original and transformed data\n")
cat("  - Clearly demonstrates why transformation was needed\n")
cat("  - Good for main text or supplement\n")
cat("  - Recommended for: Papers emphasizing methodological rigor\n\n")

cat("Option 2: COMPACT SINGLE PLOT (plot_compact)\n")
cat("  - Shows only log-transformed (recommended model)\n")
cat("  - Dual axes show both log and original values\n")
cat("  - More space-efficient\n")
cat("  - Recommended for: Space-limited journals\n\n")

cat("My recommendation: Use Option 1 (two-panel) in main text.\n")
cat("It demonstrates scientific rigor and transparency.\n\n")