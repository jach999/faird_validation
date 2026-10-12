# ==============================================================================
# H1: AMMOD INTRA-DEVICE CONSISTENCY ANALYSIS
# ==============================================================================
# Description: Complete analysis of AMMOD device consistency with regression diagnostics
# Date: October 30, 2025
# Decision: Use log₁₀(x+1) transformation (standardized, following H1_consistency.R)
# Note: Analysis limited to 12 days (first bottle batch only)
#       AMMOD3 excluded due to mechanical failure (Aug 30-Sep 3)
#       AMMOD4 excluded due to complete data loss
#       Only Maize (Site1 and Site2) available
# ==============================================================================

# --- Setup ---
library(tidyverse)
library(car)      # For ncvTest (Breusch-Pagan)
library(MASS)     # For robust regression
library(patchwork)  # For combining plots
library(ggpubr)     # For stat_cor

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

# Load custom colors
source("custom_colors.R")

cat("✅ Custom colors loaded\n")

# Load data
cat("\n========================================\n")
cat("H1: AMMOD CONSISTENCY ANALYSIS\n")
cat("========================================\n\n")

# Load pre-cleaned data
source("load_data.R")
load("data_cleaned.RData")

# CRITICAL: Use ammod_biomass_for_stats (excludes AMMOD3 problematic days)
ammod_biomass <- data_list$ammod_biomass_for_stats

cat("Note: Using ammod_biomass_for_stats (AMMOD3 Aug 30-Sep 3 excluded)\n")
cat("Note: Only first bottle batch available (12 days)\n")
cat("Note: AMMOD4 excluded (complete data loss)\n\n")

# ==============================================================================
# SECTION 1: REGRESSION DIAGNOSTICS
# ==============================================================================

cat("\n=== SECTION 1: REGRESSION DIAGNOSTICS ===\n\n")

# Prepare data for AMMOD Maize (AMMOD1 vs AMMOD2)
# CORRECTED: use values_fill = NA to avoid spurious zeros
counts_data_ammod <- ammod_biomass %>%
  dplyr::filter(Device %in% c("AMMOD1", "AMMOD2")) %>%
  dplyr::select(Date, Device, Live_mass) %>%
  tidyr::pivot_wider(names_from = Device, values_from = Live_mass,
                     values_fill = NA) %>%        # ← NA, not 0
  dplyr::filter(!is.na(AMMOD1) & !is.na(AMMOD2)) # ← drop incomplete rows

cat("Data prepared:", nrow(counts_data_ammod), "observations\n")
cat("Date range:", min(counts_data_ammod$Date), "to", max(counts_data_ammod$Date), "\n")
cat("Note: AMMOD3 excluded (mechanical failure), AMMOD4 excluded (data loss)\n\n")

# --- 1.1: Original Scale Diagnostics ---
cat("--- Original Scale Analysis ---\n")

model_original <- lm(AMMOD2 ~ AMMOD1, data = counts_data_ammod)

cat("R² =", round(summary(model_original)$r.squared, 3), "\n")
cat("Adjusted R² =", round(summary(model_original)$adj.r.squared, 3), "\n\n")

# Residual normality
shapiro_orig <- shapiro.test(residuals(model_original))
cat("Shapiro-Wilk (residuals): p =", round(shapiro_orig$p.value, 4), 
    ifelse(shapiro_orig$p.value > 0.05, "✓ Normal", "✗ NOT normal"), "\n")

# Homoscedasticity (Breusch-Pagan test)
bp_test <- car::ncvTest(model_original)
cat("Breusch-Pagan (homoscedasticity): p =", round(bp_test$p, 4),
    ifelse(bp_test$p > 0.05, "✓ Homoscedastic", "✗ Heteroscedastic"), "\n\n")

# Cook's Distance
cooks_d <- cooks.distance(model_original)
threshold_cooks <- 4/nrow(counts_data_ammod)
influential <- which(cooks_d > threshold_cooks)

cat("Cook's Distance Analysis:\n")
cat("  Threshold:", round(threshold_cooks, 4), "\n")
cat("  Influential points:", length(influential), "/", nrow(counts_data_ammod), "\n")

if(length(influential) > 0) {
  cat("\n  Influential observations:\n")
  influential_df <- counts_data_ammod[influential, ] %>%
    dplyr::mutate(cooks_d = cooks_d[influential]) %>%
    dplyr::arrange(desc(cooks_d))
  print(influential_df)
}

# --- 1.2: Log-Transformed Diagnostics ---
cat("\n\n--- Log-Transformed Analysis ---\n")

counts_log <- counts_data_ammod %>%
  dplyr::mutate(
    AMMOD1_log = log10(AMMOD1 + 1),
    AMMOD2_log = log10(AMMOD2 + 1)
  )

model_log <- lm(AMMOD2_log ~ AMMOD1_log, data = counts_log)

cat("R² =", round(summary(model_log)$r.squared, 3), "\n")
cat("Adjusted R² =", round(summary(model_log)$adj.r.squared, 3), "\n\n")

# Residual normality
shapiro_log <- shapiro.test(residuals(model_log))
cat("Shapiro-Wilk (residuals): p =", round(shapiro_log$p.value, 4),
    ifelse(shapiro_log$p.value > 0.05, "✓ Normal", "✗ NOT normal"), "\n")

# Homoscedasticity
bp_test_log <- car::ncvTest(model_log)
cat("Breusch-Pagan (homoscedasticity): p =", round(bp_test_log$p, 4),
    ifelse(bp_test_log$p > 0.05, "✓ Homoscedastic", "✗ Heteroscedastic"), "\n\n")

# Cook's Distance for log model
cooks_d_log <- cooks.distance(model_log)
influential_log <- which(cooks_d_log > threshold_cooks)
cat("Influential points (log scale):", length(influential_log), "/", nrow(counts_log), "\n\n")

# --- 1.3: Robust Regression Validation ---
cat("--- Robust Regression Validation ---\n")

# Original scale (validates OLS on raw biomass)
model_robust     <- MASS::rlm(AMMOD2 ~ AMMOD1, data = counts_data_ammod)
residuals_robust <- residuals(model_robust)
ss_res_robust    <- sum(residuals_robust^2)
ss_tot           <- sum((counts_data_ammod$AMMOD2 - mean(counts_data_ammod$AMMOD2))^2)
pseudo_r2_robust <- 1 - ss_res_robust / ss_tot
n_low_weight_orig <- sum(model_robust$w < 0.5)

cat(sprintf("Original scale (n=12): Pseudo-R2 = %.3f | Low-weight days (<0.5): %d/12\n\n",
            pseudo_r2_robust, n_low_weight_orig))

# --- 1.4: Summary Table ---
cat("=== DIAGNOSTIC SUMMARY ===\n\n")

diagnostic_summary <- tibble::tibble(
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
  Decision = c(
    ifelse(bp_test$p < 0.05 | length(influential) > 2, "⚠️ Issues", "✓ OK"),
    ifelse(shapiro_log$p.value > 0.05 & bp_test_log$p > 0.05, "✅ USE THIS", "⚠️ Check"),
    "✓ Validation"
  )
) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 3)))

print(diagnostic_summary)

cat("\n✅ DECISION: Use log-transformed data (R²=", round(summary(model_log)$r.squared, 3), ")\n", sep = "")
cat("    Rationale: Following standardized methodology from H1_consistency.R\n\n")

# ==============================================================================
# SECTION 2: DETAILED CONSISTENCY METRICS
# ==============================================================================

cat("\n=== SECTION 2: CONSISTENCY METRICS (LOG-TRANSFORMED) ===\n\n")

# Function for log-transformed consistency analysis
analyze_consistency_log <- function(daily_data, device1, device2, metric = "Live_mass") {
  
  df <- daily_data %>%
    dplyr::filter(Device %in% c(device1, device2)) %>%
    dplyr::select(Date, Device, !!sym(metric)) %>%
    tidyr::pivot_wider(names_from = Device, values_from = !!sym(metric), values_fill = NA) %>%
    dplyr::mutate(
      dev1_log = log10(!!sym(device1) + 1),
      dev2_log = log10(!!sym(device2) + 1)
    )
  
  if(nrow(df) < 3) {
    warning(sprintf("Insufficient data for %s vs %s", device1, device2))
    return(NULL)
  }
  
  # Pearson on log-transformed
  cor_pearson_log <- cor.test(df$dev1_log, df$dev2_log, method = "pearson")
  
  # Spearman on original (non-parametric verification)
  cor_spearman <- cor.test(df[[device1]], df[[device2]], method = "spearman")
  
  # OLS regression on log
  model_log     <- lm(dev2_log ~ dev1_log, data = df)
  r_squared_log <- summary(model_log)$r.squared
  adj_r_squared <- summary(model_log)$adj.r.squared
  
  # Slope and 95% CI from OLS (log scale)
  model_coef  <- coef(summary(model_log))
  slope       <- model_coef["dev1_log", "Estimate"]
  slope_se    <- model_coef["dev1_log", "Std. Error"]
  slope_ci    <- confint(model_log)["dev1_log", ]
  intercept   <- model_coef["(Intercept)", "Estimate"]
  intercept_p <- model_coef["(Intercept)", "Pr(>|t|)"]
  
  # Robust regression on log scale (validates primary OLS model)
  model_robust_log  <- MASS::rlm(dev2_log ~ dev1_log, data = df)
  ss_res_log        <- sum(residuals(model_robust_log)^2)
  ss_tot_log        <- sum((df$dev2_log - mean(df$dev2_log))^2)
  pseudo_r2_log     <- 1 - ss_res_log / ss_tot_log
  n_low_weight      <- sum(model_robust_log$w < 0.5)
  
  list(
    devices       = c(device1, device2),
    metric        = metric,
    n_days        = nrow(df),
    r_squared_log = r_squared_log,
    adj_r_squared = adj_r_squared,
    pearson_r_log = as.numeric(cor_pearson_log$estimate),
    pearson_p_log = cor_pearson_log$p.value,
    spearman_rho  = as.numeric(cor_spearman$estimate),
    spearman_p    = cor_spearman$p.value,
    slope         = slope,
    slope_se      = slope_se,
    slope_ci_low  = slope_ci[1],
    slope_ci_high = slope_ci[2],
    intercept     = intercept,
    intercept_p   = intercept_p,
    pseudo_r2_log = pseudo_r2_log,
    n_low_weight  = n_low_weight,
    data          = df,
    model         = model_log,
    model_robust  = model_robust_log
  )
}

# Helper function for significance stars
get_sig_stars <- function(p) {
  if(p < 0.001) return("***")
  if(p < 0.01) return("**")
  if(p < 0.05) return("*")
  return("")
}

# Analyze AMMOD1 vs AMMOD2
cat("--- AMMOD1 vs AMMOD2 (Maize) ---\n")
result_ammod <- analyze_consistency_log(ammod_biomass, "AMMOD1", "AMMOD2", "Live_mass")

if(!is.null(result_ammod)) {
  cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f%s\n",
              result_ammod$n_days, result_ammod$r_squared_log,
              result_ammod$adj_r_squared, result_ammod$pearson_p_log,
              get_sig_stars(result_ammod$pearson_p_log)))
  cat(sprintf("   Slope: %.3f +/- %.3f SE  [95%% CI: %.3f - %.3f]  (Intercept: %.3f, p=%.3f)\n",
              result_ammod$slope, result_ammod$slope_se,
              result_ammod$slope_ci_low, result_ammod$slope_ci_high,
              result_ammod$intercept, result_ammod$intercept_p))
  cat(sprintf("   Spearman rho=%.3f, p=%.4f\n",
              result_ammod$spearman_rho, result_ammod$spearman_p))
  cat(sprintf("   Pseudo-R2 (robust, log scale): %.3f  |  Low-weight days (<0.5): %d/%d\n\n",
              result_ammod$pseudo_r2_log, result_ammod$n_low_weight, result_ammod$n_days))
}

# ==============================================================================
# SECTION 3: EXPORT RESULTS TABLE
# ==============================================================================

cat("\n=== SECTION 3: EXPORT RESULTS ===\n\n")

if(!is.null(result_ammod)) {
  results_table <- tibble::tibble(
    Device_Pair = "AMMOD1 vs AMMOD2",
    Sites = "Site1 vs Site2 (Maize)",
    Metric = "Live biomass (mg)",
    N_days = result_ammod$n_days,
    R_squared_log = result_ammod$r_squared_log,
    Pearson_r_log = result_ammod$pearson_r_log,
    Pearson_p_log = result_ammod$pearson_p_log,
    Spearman_rho = result_ammod$spearman_rho,
    Spearman_p = result_ammod$spearman_p,
    Shapiro_p_log = shapiro_log$p.value,
    BP_p_log = bp_test_log$p,
    Pseudo_R2_robust = pseudo_r2_robust,
    Decision = "Log-transformed",
    Note = "Inter-site comparison, not true replicates"
  ) %>%
    dplyr::mutate(across(where(is.numeric), ~round(.x, 4)))
  
  readr::write_csv(results_table, "outputs/Table_H1_AMMOD_Consistency_Results.csv")
  
  cat("✅ Results exported to: outputs/Table_H1_AMMOD_Consistency_Results.csv\n\n")
  print(results_table)
}

# ==============================================================================
# SECTION 4: DIAGNOSTIC PLOTS (SUPPLEMENTARY)
# ==============================================================================

cat("\n=== SECTION 4: DIAGNOSTIC COMPARISON PLOTS ===\n\n")

# Prepare data with influential point classification
counts_data_ammod_plot <- counts_data_ammod %>%
  dplyr::mutate(
    AMMOD1_log = log10(AMMOD1 + 1),
    AMMOD2_log = log10(AMMOD2 + 1),
    influential = cooks_d > threshold_cooks  # Using original scale influential points
  )

# Create dynamic labels with actual values
label_original <- sprintf(
  "R² = %.3f\nShapiro p = %.3f\nBP p = %.3f",
  summary(model_original)$r.squared,
  shapiro_orig$p.value,
  bp_test$p
)

label_log <- sprintf(
  "R² = %.3f\nShapiro p = %.3f\nBP p = %.3f",
  summary(model_log)$r.squared,
  shapiro_log$p.value,
  bp_test_log$p
)

# --- Plot A: Original Scale ---
plot_original <- ggplot(counts_data_ammod_plot, aes(x = AMMOD1, y = AMMOD2)) +
  geom_point(aes(color = influential, size = influential), alpha = 0.7) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray30"),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  scale_size_manual(
    values = c("TRUE" = 4, "FALSE" = 3),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = custom_colors["AMMOD"],  # Color for AMMOD
    alpha = 0.3, linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.8) +
  annotate("text", x = Inf, y = -Inf,
           label = label_original,
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20",
           fontface = "italic") +
  labs(
    title = "A) Original Scale",
    x = "AMMOD1 Daily Biomass (mg)",
    y = "AMMOD2 Daily Biomass (mg)",
    color = "Point type",
    size = "Point type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  coord_fixed(ratio = 1)

# --- Plot B: Log-Transformed Scale ---
plot_log <- ggplot(counts_data_ammod_plot, aes(x = AMMOD1_log, y = AMMOD2_log)) +
  geom_point(aes(color = influential, size = influential), alpha = 0.7) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray30"),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  scale_size_manual(
    values = c("TRUE" = 4, "FALSE" = 3),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = custom_colors["AMMOD"],  # Color for AMMOD
    alpha = 0.3, linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.8) +
  annotate("text", x = Inf, y = -Inf,
           label = label_log,
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20",
           fontface = "italic") +
  labs(
    title = "B) Log₁₀-Transformed Scale",
    x = "log₁₀(AMMOD1 + 1)",
    y = "log₁₀(AMMOD2 + 1)",
    color = "Point type",
    size = "Point type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  coord_fixed(ratio = 1)

# --- Combined Plot ---
combined_plot <- plot_original + plot_log +
  plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom")

print(combined_plot)

# Save plots
ggsave("outputs/H1_AMMOD_Diagnostic_Comparison.png", 
       plot = combined_plot, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_Diagnostic_Comparison.pdf",
       plot = combined_plot,
       width = 12, height = 6, device = cairo_pdf)

cat("✅ Diagnostic plots saved:\n")
cat("   - outputs/H1_AMMOD_Diagnostic_Comparison.png\n")
cat("   - outputs/H1_AMMOD_Diagnostic_Comparison.pdf\n\n")

# ==============================================================================
# SECTION 5: MAIN TEXT FIGURE (LOG-TRANSFORMED ONLY)
# ==============================================================================

cat("\n=== SECTION 5: MAIN TEXT FIGURE ===\n\n")

# Single panel plot for AMMOD comparison (log-transformed)
plot_ammod_main <- ggplot(result_ammod$data, aes(x = dev1_log, y = dev2_log)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = custom_colors["AMMOD"],  # Color for AMMOD
    alpha = 0.3
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  stat_cor(method = "pearson",aes(label = paste(after_stat(rr.label),after_stat(p.label),sep = "~`,`~")), label.x = 2.5, label.y = 2.5, size = 3.5) +
  labs(
    title = "AMMOD Inter-Site Consistency (AMMOD1 vs AMMOD2)",
    subtitle = "Site1 vs Site2 (Maize) | n = 12 days",
    x = "log₁₀(AMMOD1 Biomass + 1)",
    y = "log₁₀(AMMOD2 Biomass + 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    panel.border = element_rect(color = "gray80", fill = NA)
  ) #+
  #coord_fixed(ratio = 1)

print(plot_ammod_main)

# Save main figure
ggsave("outputs/H1_AMMOD_Consistency.png", 
       plot = plot_ammod_main, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_Consistency.pdf",
       plot = plot_ammod_main,
       width = 7, height = 6, device = cairo_pdf)

cat("✅ Main text figure saved:\n")
cat("   - outputs/H1_AMMOD_Consistency.png\n")
cat("   - outputs/H1_AMMOD_Consistency.pdf\n\n")


# ==============================================================================
# TEMPORAL CONCORDANCE ANALYSIS (DIRECTIONAL SYNCHRONIZATION)
# ==============================================================================

cat("\n=== TEMPORAL CONCORDANCE (DIRECTIONAL SYNCHRONIZATION) ===\n\n")

# Calculate directional concordance (same methodology as H2_validation.R)
temporal_concordance_ammod <- counts_data_ammod %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    AMMOD1_change = AMMOD1 - lag(AMMOD1),
    AMMOD2_change = AMMOD2 - lag(AMMOD2),
    concordant = sign(AMMOD1_change) == sign(AMMOD2_change)
  ) %>%
  dplyr::filter(!is.na(concordant))

# Calculate concordance percentage
n_changes_ammod <- nrow(temporal_concordance_ammod)
n_concordant_ammod <- sum(temporal_concordance_ammod$concordant, na.rm = TRUE)
concordance_pct_ammod <- (n_concordant_ammod / n_changes_ammod) * 100

cat(sprintf("Temporal Concordance: %.1f%%\n", concordance_pct_ammod))
cat(sprintf("  Days analyzed: %d (n-1 for lag comparison)\n", n_changes_ammod))
cat(sprintf("  Concordant days: %d\n", n_concordant_ammod))
cat(sprintf("  Discordant days: %d\n\n", n_changes_ammod - n_concordant_ammod))

# Show discordant days if any
if((n_changes_ammod - n_concordant_ammod) > 0) {
  cat("Days with opposite directional changes:\n")
  discordant_days_ammod <- temporal_concordance_ammod %>%
    dplyr::filter(!concordant) %>%
    dplyr::select(Date, AMMOD1, AMMOD2, AMMOD1_change, AMMOD2_change) %>%
    dplyr::mutate(
      Direction = case_when(
        AMMOD1_change > 0 & AMMOD2_change < 0 ~ "AMMOD1 ↑, AMMOD2 ↓",
        AMMOD1_change < 0 & AMMOD2_change > 0 ~ "AMMOD1 ↓, AMMOD2 ↑",
        TRUE ~ "Other"
      )
    )
  print(discordant_days_ammod)
  cat("\n")
}

# Kendall's tau test (non-parametric correlation of changes)
kendall_test_ammod <- cor.test(temporal_concordance_ammod$AMMOD1_change, 
                               temporal_concordance_ammod$AMMOD2_change, 
                               method = "kendall", 
                               exact = FALSE)

cat("Kendall's tau (correlation of directional changes):\n")
cat(sprintf("  τ = %.3f\n", kendall_test_ammod$estimate))
cat(sprintf("  p-value = %.4f\n\n", kendall_test_ammod$p.value))

# ==============================================================================
# SECTION 6: TEMPORAL CONCORDANCE PLOT AND ANALYSIS
# ==============================================================================

cat("\n=== SECTION 6: TEMPORAL CONCORDANCE ===\n\n")

# --- Prepare long format data for temporal plot ---
ammod_temporal <- counts_data_ammod %>%
  dplyr::select(Date, AMMOD1, AMMOD2) %>%
  tidyr::pivot_longer(
    cols = c(AMMOD1, AMMOD2),
    names_to = "Device",
    values_to = "Biomass"
  ) %>%
  dplyr::mutate(
    Biomass_log = log10(Biomass + 1)
  )

# --- Create temporal plot (log-transformed) ---
plot_temporal_ammod <- ggplot(ammod_temporal, 
                              aes(x = Date, y = Biomass_log, 
                                  color = Device, shape = Device)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
    labels = c("AMMOD1" = "AMMOD1", 
               "AMMOD2" = "AMMOD2")
  ) +
  scale_shape_manual(
    values = c("AMMOD1" = 16, "AMMOD2" = 17),  # Circle and triangle
    labels = c("AMMOD1" = "AMMOD1", 
               "AMMOD2" = "AMMOD2")
  ) +
  labs(
    title = "Temporal Concordance: AMMOD1 vs AMMOD2",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       result_ammod$r_squared_log,
                       concordance_pct_ammod,
                       kendall_test_ammod$estimate,
                       kendall_test_ammod$p.value),
    x = "Date",
    y = "log₁₀(Live Biomass + 1)",
    color = "Device",
    shape = "Device"
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

print(plot_temporal_ammod)

# Save temporal plot
ggsave("outputs/H1_AMMOD_Temporal_Concordance.png", 
       plot = plot_temporal_ammod, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_Temporal_Concordance.pdf",
       plot = plot_temporal_ammod,
       width = 7, height = 6, device = cairo_pdf)

cat("✅ Temporal concordance plot saved:\n")
cat("   - outputs/H1_AMMOD_Temporal_Concordance.png\n")
cat("   - outputs/H1_AMMOD_Temporal_Concordance.pdf\n\n")

# ------------------------------------------------------------------------------
# COMBINED FIGURE (A: scatter | B: temporal concordance)
# ------------------------------------------------------------------------------

cat("--- Generating combined figure (A: scatter | B: temporal) ---\n\n")

plot_ammod_main_A <- plot_ammod_main +
  labs(title = "A) AMMOD Inter-Site Consistency (AMMOD1 vs AMMOD2)")

plot_temporal_ammod_B <- plot_temporal_ammod +
  labs(title = "B) Temporal Concordance: AMMOD1 vs AMMOD2")

combined_ammod <- plot_ammod_main_A + plot_temporal_ammod_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_ammod)

ggsave("outputs/H1_AMMOD_Combined.png",
       plot   = combined_ammod,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_AMMOD_Combined.pdf",
       plot   = combined_ammod,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved:\n")
cat("   - outputs/H1_AMMOD_Combined.png\n")
cat("   - outputs/H1_AMMOD_Combined.pdf\n\n")

# ==============================================================================
# TEMPORAL PATTERN ANALYSIS
# ==============================================================================

cat("\n=== TEMPORAL PATTERN ANALYSIS ===\n\n")

# Calculate day-by-day differences in log scale
temporal_analysis_ammod <- counts_data_ammod %>%
  dplyr::mutate(
    AMMOD1_log = log10(AMMOD1 + 1),
    AMMOD2_log = log10(AMMOD2 + 1),
    Difference_log = abs(AMMOD1_log - AMMOD2_log),
    Direction = case_when(
      AMMOD1_log > AMMOD2_log ~ "AMMOD1 > AMMOD2",
      AMMOD1_log < AMMOD2_log ~ "AMMOD2 > AMMOD1",
      TRUE ~ "Equal"
    )
  )

# Summary statistics
mean_diff_ammod <- mean(temporal_analysis_ammod$Difference_log, na.rm = TRUE)
sd_diff_ammod <- sd(temporal_analysis_ammod$Difference_log, na.rm = TRUE)
max_diff_ammod <- max(temporal_analysis_ammod$Difference_log, na.rm = TRUE)

cat("Day-to-day divergence (log-scale):\n")
cat(sprintf("  Mean absolute difference: %.3f\n", mean_diff_ammod))
cat(sprintf("  SD of differences: %.3f\n", sd_diff_ammod))
cat(sprintf("  Maximum difference: %.3f\n\n", max_diff_ammod))

# Days where one device has much higher activity
threshold_diff_ammod <- 0.5  # log-scale difference threshold
divergent_days_ammod <- temporal_analysis_ammod %>%
  dplyr::filter(Difference_log > threshold_diff_ammod) %>%
  dplyr::arrange(desc(Difference_log))

if(nrow(divergent_days_ammod) > 0) {
  cat(sprintf("Days with large divergence (>%.1f log units): %d/%d days\n", 
              threshold_diff_ammod, 
              nrow(divergent_days_ammod), 
              nrow(temporal_analysis_ammod)))
  cat("\nTop divergent days:\n")
  print(divergent_days_ammod %>% 
          dplyr::select(Date, AMMOD1, AMMOD2, Difference_log, Direction) %>%
          head(5))
} else {
  cat("✅ No days with large divergence - strong concordance throughout\n")
}

cat("\n=== INTERPRETATION ===\n")
if(mean_diff_ammod < 0.3) {
  cat("✅ Excellent temporal concordance: Mean difference < 0.3 log units\n")
  cat("   Both devices show highly synchronized daily dynamics\n")
} else if(mean_diff_ammod < 0.5) {
  cat("✓ Good temporal concordance: Mean difference < 0.5 log units\n")
  cat("   General synchronization with some site-specific variation\n")
} else {
  cat("⚠️  Moderate temporal concordance: Mean difference > 0.5 log units\n")
  cat("   Substantial site-specific dynamics present\n")
}

cat("\n=== TEMPORAL CONCORDANCE COMPLETE ===\n")


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("H1 AMMOD ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n")
cat("  ✅ outputs/Table_H1_AMMOD_Consistency_Results.csv\n")
cat("  ✅ outputs/H1_AMMOD_Diagnostic_Comparison.png (Supplementary)\n")
cat("  ✅ outputs/H1_AMMOD_Diagnostic_Comparison.pdf (Supplementary)\n")
cat("  ✅ outputs/H1_AMMOD_Consistency.png (Main text)\n")
cat("  ✅ outputs/H1_AMMOD_Consistency.pdf (Main text)\n\n")

# Dynamic final summary
cat("KEY FINDINGS:\n")
if (!is.null(result_ammod)) {
  cat(sprintf("  • AMMOD Inter-Site (AMMOD1 vs AMMOD2): R² = %.3f%s\n",
              result_ammod$r_squared_log,
              get_sig_stars(result_ammod$pearson_p_log)))
  cat(sprintf("  • Directional concordance: %.1f%% (%d/%d days)\n",
              concordance_pct_ammod,
              n_concordant_ammod,
              n_changes_ammod))
  
  # Interpretation
  if (result_ammod$r_squared_log < 0.5) {
    cat("  • INTERPRETATION: Poor to moderate inter-site consistency\n")
  } else if (result_ammod$r_squared_log < 0.7) {
    cat("  • INTERPRETATION: Moderate inter-site consistency\n")
  } else {
    cat("  • INTERPRETATION: Strong inter-site consistency\n")
  }
}

cat("\nSTATISTICAL DECISION:\n")
cat("  • Use log10(x+1) transformation (standardized)\n")
cat("  • Rationale: Consistent methodology with FAIRD and ID analyses\n")
cat(sprintf("  • Robust regression (log scale,      n=12): Pseudo-R2 = %.3f | Low-weight: %d/12\n",
            result_ammod$pseudo_r2_log, result_ammod$n_low_weight))
cat(sprintf("  • Robust regression (original scale, n=12): Pseudo-R2 = %.3f | Low-weight: %d/12\n\n",
            pseudo_r2_robust, n_low_weight_orig))

cat("COMPARISON WITH FAIRD AND ID:\n")
if (!is.null(result_ammod)) {
  cat(sprintf("  • AMMOD (different sites): R² = %.3f\n", result_ammod$r_squared_log))
 }

# ==============================================================================
cat("========================================\n\n")

cat("\n========================================\n")
cat("H1 AMMOD COMPLETE\n")
cat("========================================\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")
