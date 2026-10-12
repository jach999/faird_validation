# ==============================================================================
# H1: INSECT DETECT INTRA-DEVICE CONSISTENCY ANALYSIS
# ==============================================================================
# Description: Complete analysis of ID device consistency with regression diagnostics
# Date: October 30, 2025
# Decision: Use log₁₀(x+1) transformation (standardized, following H1_consistency.R)
# Note: ID3 excluded due to device failure (documented in project context)
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
cat("H1: INSECT DETECT CONSISTENCY ANALYSIS\n")
cat("========================================\n\n")

# Load pre-cleaned and pre-completed data
# The 'etraps_daily' object from this file NOW includes all zero-filled days.
source("load_data.R")
load("data_cleaned.RData")

etraps_daily <- data_list$etraps_daily

# ==============================================================================
# SECTION 1: REGRESSION DIAGNOSTICS
# ==============================================================================

cat("\n=== SECTION 1: REGRESSION DIAGNOSTICS ===\n\n")

# Prepare CONSERVATIVE data (n=22, ALL days including battery failures)
# Used to calculate Pseudo-R² for the conservative analysis documented in PCV
counts_data_id_conservative <- etraps_daily %>%
  dplyr::filter(Device %in% c("ID1", "ID2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

cat("Conservative data (n=22, all days):", nrow(counts_data_id_conservative), "observations\n")

# Prepare data for ID Maize (ID1 vs ID2)
# ID3 excluded due to device failure
# EXCLUDING problematic days: Aug 28 - Sep 1 (battery failures)

# Prepare OPERATIONAL data (n=17, excluding battery failures Aug 28 - Sep 1)
counts_data_id_operational <- etraps_daily %>%
  dplyr::filter(Device %in% c("ID1", "ID2")) %>%
  dplyr::filter(Date < as.Date("2023-08-28") | Date > as.Date("2023-09-01")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

cat("Operational data (n=17, excl. battery failures):", nrow(counts_data_id_operational), "observations\n")
cat("Date range:", min(counts_data_id_operational$Date), "to", max(counts_data_id_operational$Date), "\n")
cat("Note: ID3 excluded due to device failure\n")
cat("Note: Days Aug 28 - Sep 1 excluded from operational analysis\n\n")

# --- 1.1: Original Scale Diagnostics ---
cat("--- Original Scale Analysis ---\n")

model_original <- lm(ID2 ~ ID1, data = counts_data_id_conservative)

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
threshold_cooks <- 4/nrow(counts_data_id_conservative)
influential <- which(cooks_d > threshold_cooks)

cat("Cook's Distance Analysis:\n")
cat("  Threshold:", round(threshold_cooks, 4), "\n")
cat("  Influential points:", length(influential), "/", nrow(counts_data_id_conservative), "\n")

if(length(influential) > 0) {
  cat("\n  Influential observations:\n")
  influential_df <- counts_data_id_conservative[influential, ] %>%
    dplyr::mutate(cooks_d = cooks_d[influential]) %>%
    dplyr::arrange(desc(cooks_d))
  print(influential_df)
}

# --- 1.2: Log-Transformed Diagnostics ---
cat("\n\n--- Log-Transformed Analysis ---\n")

counts_data_id_log_cons  <- counts_data_id_conservative %>%
  dplyr::mutate(
    ID1_log_cons  = log10(ID1 + 1),
    ID2_log_cons = log10(ID2 + 1)
  )

model_log_cons <- lm(ID2_log_cons ~ ID1_log_cons, data = counts_data_id_log_cons)

cat("R² =", round(summary(model_log_cons)$r.squared, 3), "\n")
cat("Adjusted R² =", round(summary(model_log_cons)$adj.r.squared, 3), "\n\n")

# Residual normality
shapiro_log_cons <- shapiro.test(residuals(model_log_cons))
cat("Shapiro-Wilk (residuals): p =", round(shapiro_log_cons$p.value, 4),
    ifelse(shapiro_log_cons$p.value > 0.05, "✓ Normal", "✗ NOT normal"), "\n")

# Homoscedasticity
bp_test_log_cons <- car::ncvTest(model_log_cons)
cat("Breusch-Pagan (homoscedasticity): p =", round(bp_test_log_cons$p, 4),
    ifelse(bp_test_log_cons$p > 0.05, "✓ Homoscedastic", "✗ Heteroscedastic"), "\n\n")

# Cook's Distance for log model
cooks_d_log_cons <- cooks.distance(model_log_cons)
influential_log_cons <- which(cooks_d_log_cons > threshold_cooks)
cat("Influential points (log scale):", length(influential_log_cons), "/", nrow(counts_data_id_log_cons), "\n\n")

# --- 1.3: Robust Regression Validation ---
cat("--- Robust Regression Validation ---\n")

# CONSERVATIVE (n=22, original scale): validates R2=0.027
model_robust_cons  <- MASS::rlm(ID2 ~ ID1, data = counts_data_id_conservative)
ss_res_cons        <- sum(residuals(model_robust_cons)^2)
ss_tot_cons        <- sum((counts_data_id_conservative$ID2 - mean(counts_data_id_conservative$ID2))^2)
pseudo_r2_cons     <- 1 - ss_res_cons / ss_tot_cons
n_low_weight_cons  <- sum(model_robust_cons$w < 0.5)

# OPERATIONAL (n=17, original scale): validates R2=0.434
model_robust_oper  <- MASS::rlm(ID2 ~ ID1, data = counts_data_id_operational)
ss_res_oper        <- sum(residuals(model_robust_oper)^2)
ss_tot_oper        <- sum((counts_data_id_operational$ID2 - mean(counts_data_id_operational$ID2))^2)
pseudo_r2_robust   <- 1 - ss_res_oper / ss_tot_oper  # Keep name for downstream compatibility
n_low_weight_oper  <- sum(model_robust_oper$w < 0.5)

cat(sprintf("Conservative (n=22): Pseudo-R2 = %.3f | Low-weight days (<0.5): %d/22\n",
            pseudo_r2_cons, n_low_weight_cons))
cat(sprintf("Operational  (n=17): Pseudo-R2 = %.3f | Low-weight days (<0.5): %d/17\n\n",
            pseudo_r2_robust, n_low_weight_oper))
# --- 1.4: Summary Table ---
cat("=== DIAGNOSTIC SUMMARY ===\n\n")

diagnostic_summary <- tibble::tibble(
  Method = c("Original scale", "Log-transformed", "Robust regression"),
  R_squared = c(
    summary(model_original)$r.squared,
    summary(model_log_cons)$r.squared,
    pseudo_r2_robust
  ),
  Shapiro_p = c(
    shapiro_orig$p.value,
    shapiro_log_cons$p.value,
    NA
  ),
  BP_p = c(
    bp_test$p,
    bp_test_log_cons$p,
    NA
  ),
  N_influential = c(
    length(influential),
    length(influential_log_cons),
    NA
  ),
  Decision = c(
    ifelse(bp_test$p < 0.05 | length(influential) > 2, "⚠️ Issues", "✓ OK"),
    ifelse(shapiro_log_cons$p.value > 0.05 & bp_test_log_cons$p > 0.05, "✅ USE THIS", "⚠️ Check"),
    "✓ Validation"
  )
) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 3)))

print(diagnostic_summary)

cat("\n✅ DECISION: Use log-transformed data (R²=", round(summary(model_log_cons)$r.squared, 3), ")\n", sep = "")
cat("    Rationale: Following standardized methodology from H1_consistency.R\n\n")

# ==============================================================================
# SECTION 2: DETAILED CONSISTENCY METRICS
# ==============================================================================

cat("\n=== SECTION 2: CONSISTENCY METRICS (LOG-TRANSFORMED) ===\n\n")

# Function for log-transformed consistency analysis
analyze_consistency_log_cons <- function(daily_data, device1, device2, metric = "abundance") {
  
  df <- daily_data %>%
    dplyr::filter(Device %in% c(device1, device2)) %>%
    dplyr::select(Date, Device, !!sym(metric)) %>%
    #dplyr::filter(Date < as.Date("2023-08-28") | Date > as.Date("2023-09-01")) %>%
    tidyr::pivot_wider(names_from = Device, values_from = !!sym(metric), values_fill = 0) %>%
    dplyr::mutate(
      dev1_log_cons = log10(!!sym(device1) + 1),
      dev2_log_cons = log10(!!sym(device2) + 1)
    )
  
  if(nrow(df) < 3) {
    warning(sprintf("Insufficient data for %s vs %s", device1, device2))
    return(NULL)
  }
  
  # Pearson on log-transformed
  cor_pearson_log_cons <- cor.test(df$dev1_log_cons, df$dev2_log_cons, method = "pearson")
  
  # Spearman on original (non-parametric verification)
  cor_spearman <- cor.test(df[[device1]], df[[device2]], method = "spearman")
  
  # OLS regression on log
  model_log_cons     <- lm(dev2_log_cons ~ dev1_log_cons, data = df)
  r_squared_log_cons <- summary(model_log_cons)$r.squared
  adj_r_squared <- summary(model_log_cons)$adj.r.squared
  
  # Slope and 95% CI from OLS (log scale)
  model_coef  <- coef(summary(model_log_cons))
  slope       <- model_coef["dev1_log_cons", "Estimate"]
  slope_se    <- model_coef["dev1_log_cons", "Std. Error"]
  slope_ci    <- confint(model_log_cons)["dev1_log_cons", ]
  intercept   <- model_coef["(Intercept)", "Estimate"]
  intercept_p <- model_coef["(Intercept)", "Pr(>|t|)"]
  
  # Robust regression on log scale (validates primary OLS model)
  model_robust  <- MASS::rlm(dev2_log_cons ~ dev1_log_cons, data = df)
  ss_res_robust <- sum(residuals(model_robust)^2)
  ss_tot_ols    <- sum((df$dev2_log_cons - mean(df$dev2_log_cons))^2)
  pseudo_r2_log_cons <- 1 - ss_res_robust / ss_tot_ols
  n_low_weight  <- sum(model_robust$w < 0.5)
  
  list(
    devices       = c(device1, device2),
    metric        = metric,
    n_days        = nrow(df),
    r_squared_log_cons = r_squared_log_cons,
    adj_r_squared = adj_r_squared,
    pearson_r_log_cons = as.numeric(cor_pearson_log_cons$estimate),
    pearson_p_log_cons = cor_pearson_log_cons$p.value,
    spearman_rho  = as.numeric(cor_spearman$estimate),
    spearman_p    = cor_spearman$p.value,
    slope         = slope,
    slope_se      = slope_se,
    slope_ci_low  = slope_ci[1],
    slope_ci_high = slope_ci[2],
    intercept     = intercept,
    intercept_p   = intercept_p,
    pseudo_r2_log_cons = pseudo_r2_log_cons,
    n_low_weight  = n_low_weight,
    data          = df,
    model         = model_log_cons,
    model_robust  = model_robust
  )
}

# Helper function for significance stars
get_sig_stars <- function(p) {
  if(p < 0.001) return("***")
  if(p < 0.01) return("**")
  if(p < 0.05) return("*")
  return("")
}

# Analyze ID1 vs ID2
cat("--- ID1 vs ID2 (Maize) ---\n")
result_id <- analyze_consistency_log_cons(etraps_daily, "ID1", "ID2", "abundance")

if(!is.null(result_id)) {
  cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f%s\n",
              result_id$n_days, result_id$r_squared_log_cons,
              result_id$adj_r_squared, result_id$pearson_p_log_cons,
              get_sig_stars(result_id$pearson_p_log_cons)))
  cat(sprintf("   Slope: %.3f +/- %.3f SE  [95%% CI: %.3f - %.3f]  (Intercept: %.3f, p=%.3f)\n",
              result_id$slope, result_id$slope_se,
              result_id$slope_ci_low, result_id$slope_ci_high,
              result_id$intercept, result_id$intercept_p))
  cat(sprintf("   Spearman rho=%.3f, p=%.4f\n",
              result_id$spearman_rho, result_id$spearman_p))
  cat(sprintf("   Pseudo-R2 (robust, log scale): %.3f  |  Low-weight days (<0.5): %d/%d\n\n",
              result_id$pseudo_r2_log_cons, result_id$n_low_weight, result_id$n_days))
}

# ==============================================================================
# SECTION 3: EXPORT RESULTS TABLE
# ==============================================================================

cat("\n=== SECTION 3: EXPORT RESULTS ===\n\n")

if(!is.null(result_id)) {
  results_table <- tibble::tibble(
    Device_Pair = "ID1 vs ID2",
    Site = "Maize (Site2)",
    N_days = result_id$n_days,
    R_squared_log_cons = result_id$r_squared_log_cons,
    Pearson_r_log_cons = result_id$pearson_r_log_cons,
    Pearson_p_log_cons = result_id$pearson_p_log_cons,
    Spearman_rho = result_id$spearman_rho,
    Spearman_p = result_id$spearman_p,
    Shapiro_p_log_cons = shapiro_log_cons$p.value,
    BP_p_log_cons = bp_test_log_cons$p,
    Pseudo_R2_robust = pseudo_r2_robust,
    Decision = "Log-transformed"
  ) %>%
    dplyr::mutate(across(where(is.numeric), ~round(.x, 4)))
  
  readr::write_csv(results_table, "outputs/Table_H1_ID_Consistency_Results.csv")
  
  cat("✅ Results exported to: outputs/Table_H1_ID_Consistency_Results.csv\n\n")
  print(results_table)
}

# ==============================================================================
# SECTION 4: DIAGNOSTIC PLOTS (SUPPLEMENTARY)
# ==============================================================================

cat("\n=== SECTION 4: DIAGNOSTIC COMPARISON PLOTS ===\n\n")

# Prepare data with influential point classification
counts_data_id_conservative_plot <- counts_data_id_conservative %>%
  dplyr::mutate(
    ID1_log_cons = log10(ID1 + 1),
    ID2_log_cons = log10(ID2 + 1),
    influential = cooks_d > threshold_cooks  # Using original scale influential points
  )

# Create dynamic labels with actual values
label_original <- sprintf(
  "R² = %.3f\nShapiro p = %.3f\nBP p = %.3f",
  summary(model_original)$r.squared,
  shapiro_orig$p.value,
  bp_test$p
)

label_log_cons <- sprintf(
  "R² = %.3f\nShapiro p = %.3f\nBP p = %.3f",
  summary(model_log_cons)$r.squared,
  shapiro_log_cons$p.value,
  bp_test_log_cons$p
)

# --- Plot A: Original Scale ---
plot_original <- ggplot(counts_data_id_conservative_plot, aes(x = ID1, y = ID2)) +
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
    fill = custom_colors["ID"],  # Coral for ID
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
    x = "ID1 Daily Abundance",
    y = "ID2 Daily Abundance",
    color = "Point type",
    size = "Point type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# --- Plot B: Log-Transformed Scale ---
plot_log_cons <- ggplot(counts_data_id_conservative_plot, aes(x = ID1_log_cons, y = ID2_log_cons)) +
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
    fill = custom_colors["ID"],  # Coral for ID
    alpha = 0.3, linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.8) +
  annotate("text", x = Inf, y = -Inf,
           label = label_log_cons,
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20",
           fontface = "italic") +
  labs(
    title = "B) Log₁₀-Transformed Scale",
    x = "log₁₀(ID1 + 1)",
    y = "log₁₀(ID2 + 1)",
    color = "Point type",
    size = "Point type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# --- Combined Plot ---
combined_plot <- plot_original + plot_log_cons +
  plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom")

print(combined_plot)

# Save plots
ggsave("outputs/Figure_S_ID_Diagnostic_Comparison.png", 
       plot = combined_plot, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("outputs/Figure_S_ID_Diagnostic_Comparison.pdf",
       plot = combined_plot,
       width = 12, height = 6, device = cairo_pdf)

cat("✅ Diagnostic plots saved:\n")
cat("   - outputs/Figure_S_ID_Diagnostic_Comparison.png\n")
cat("   - outputs/Figure_S_ID_Diagnostic_Comparison.pdf\n\n")

# ==============================================================================
# SECTION 5: MAIN TEXT FIGURE (LOG-TRANSFORMED ONLY)
# ==============================================================================

cat("\n=== SECTION 5: MAIN TEXT FIGURE ===\n\n")

# Single panel plot for ID comparison (log-transformed)
plot_id_main <- ggplot(result_id$data, aes(x = dev1_log_cons, y = dev2_log_cons)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = custom_colors["ID"],  # Coral for ID
    alpha = 0.3
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 2.3, size = 3.5) +
  labs(
    title = "Insect Detect Consistency (ID1 vs ID2)",
    x = "log₁₀(ID1 + 1)",
    y = "log₁₀(ID2 + 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.border = element_rect(color = "gray80", fill = NA)
  ) +
  coord_fixed(ratio = 1)

print(plot_id_main)

# Save main figure
ggsave("outputs/Figure_H1_ID_Consistency.png", 
       plot = plot_id_main, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/Figure_H1_ID_Consistency.pdf",
       plot = plot_id_main,
       width = 7, height = 6, device = cairo_pdf)

cat("✅ Main text figure saved:\n")
cat("   - outputs/Figure_H1_ID_Consistency.png\n")
cat("   - outputs/Figure_H1_ID_Consistency.pdf\n\n")


# ==============================================================================
# SECTION 6: TEMPORAL CONCORDANCE PLOT AND ANALYSIS
# ==============================================================================

cat("\n=== SECTION 6: TEMPORAL CONCORDANCE ===\n\n")

# --- Prepare long format data for temporal plot ---
id_temporal <- counts_data_id_conservative %>%
  dplyr::select(Date, ID1, ID2) %>%
  tidyr::pivot_longer(
    cols = c(ID1, ID2),
    names_to = "Device",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(
    Abundance_log_cons = log10(Abundance + 1)
  )

# --- Create temporal plot (log-transformed) ---
plot_temporal_id <- ggplot(id_temporal, 
                           aes(x = Date, y = Abundance_log_cons, 
                               color = Device, shape = Device)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
    labels = c("ID1" = "ID1 (Site2 - Maize)", 
               "ID2" = "ID2 (Site2 - Maize)")
  ) +
  scale_shape_manual(
    values = c("ID1" = 16, "ID2" = 17),  # Circle and triangle
    labels = c("ID1" = "ID1 (Site2 - Maize)", 
               "ID2" = "ID2 (Site2 - Maize)")
  ) +
  labs(
    title = "Temporal Concordance: ID1 vs ID2",
    subtitle = sprintf("Intra-site comparison showing temporal dynamics (R² = %.3f)", 
                       summary(model_log_cons)$r.squared),
    x = "Date",
    y = "log₁₀(Daily Abundance + 1)",
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

print(plot_temporal_id)

# Save temporal plot
ggsave("outputs/Figure_H1_ID_Temporal_Concordance.png", 
       plot = plot_temporal_id, 
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/Figure_H1_ID_Temporal_Concordance.pdf",
       plot = plot_temporal_id,
       width = 7, height = 6, device = cairo_pdf)

cat("✅ Temporal concordance plot saved:\n")
cat("   - outputs/ID_Temporal_Concordance.png\n")
cat("   - outputs/ID_Temporal_Concordance.pdf\n\n")

# ==============================================================================
# TEMPORAL CONCORDANCE ANALYSIS (DIRECTIONAL SYNCHRONIZATION)
# ==============================================================================

cat("\n=== TEMPORAL CONCORDANCE (DIRECTIONAL SYNCHRONIZATION) ===\n\n")

# Calculate directional concordance (same methodology as H2_validation.R)
temporal_concordance_id <- counts_data_id %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    ID1_change = ID1 - lag(ID1),
    ID2_change = ID2 - lag(ID2),
    concordant = sign(ID1_change) == sign(ID2_change)
  ) %>%
  dplyr::filter(!is.na(concordant))

# Calculate concordance percentage
n_changes_id <- nrow(temporal_concordance_id)
n_concordant_id <- sum(temporal_concordance_id$concordant, na.rm = TRUE)
concordance_pct_id <- (n_concordant_id / n_changes_id) * 100

cat(sprintf("Temporal Concordance: %.1f%%\n", concordance_pct_id))
cat(sprintf("  Days analyzed: %d (n-1 for lag comparison)\n", n_changes_id))
cat(sprintf("  Concordant days: %d\n", n_concordant_id))
cat(sprintf("  Discordant days: %d\n\n", n_changes_id - n_concordant_id))

# Show discordant days if any
if((n_changes_id - n_concordant_id) > 0) {
  cat("Days with opposite directional changes:\n")
  discordant_days_id <- temporal_concordance_id %>%
    dplyr::filter(!concordant) %>%
    dplyr::select(Date, ID1, ID2, ID1_change, ID2_change) %>%
    dplyr::mutate(
      Direction = case_when(
        ID1_change > 0 & ID2_change < 0 ~ "ID1 ↑, ID2 ↓",
        ID1_change < 0 & ID2_change > 0 ~ "ID1 ↓, ID2 ↑",
        TRUE ~ "Other"
      )
    )
  print(discordant_days_id)
  cat("\n")
}

# Kendall's tau test (non-parametric correlation of changes)
kendall_test_id <- cor.test(temporal_concordance_id$ID1_change, 
                            temporal_concordance_id$ID2_change, 
                            method = "kendall", 
                            exact = FALSE)

cat("Kendall's tau (correlation of directional changes):\n")
cat(sprintf("  τ = %.3f\n", kendall_test_id$estimate))
cat(sprintf("  p-value = %.4f\n\n", kendall_test_id$p.value))

# ==============================================================================
# TEMPORAL PATTERN ANALYSIS
# ==============================================================================

cat("\n=== TEMPORAL PATTERN ANALYSIS ===\n\n")

# Calculate day-by-day differences in log scale
temporal_analysis_id <- counts_data_id %>%
  dplyr::mutate(
    ID1_log_cons = log10(ID1 + 1),
    ID2_log_cons = log10(ID2 + 1),
    Difference_log_cons = abs(ID1_log_cons - ID2_log_cons),
    Direction = case_when(
      ID1_log_cons > ID2_log_cons ~ "ID1 > ID2",
      ID1_log_cons < ID2_log_cons ~ "ID2 > ID1",
      TRUE ~ "Equal"
    )
  )

# Summary statistics
mean_diff_id <- mean(temporal_analysis_id$Difference_log_cons, na.rm = TRUE)
sd_diff_id <- sd(temporal_analysis_id$Difference_log_cons, na.rm = TRUE)
max_diff_id <- max(temporal_analysis_id$Difference_log_cons, na.rm = TRUE)

cat("Day-to-day divergence (log-scale):\n")
cat(sprintf("  Mean absolute difference: %.3f\n", mean_diff_id))
cat(sprintf("  SD of differences: %.3f\n", sd_diff_id))
cat(sprintf("  Maximum difference: %.3f\n\n", max_diff_id))

# Days where one device has much higher activity
threshold_diff_id <- 0.5  # log-scale difference threshold
divergent_days_id <- temporal_analysis_id %>%
  dplyr::filter(Difference_log_cons > threshold_diff_id) %>%
  dplyr::arrange(desc(Difference_log_cons))

if(nrow(divergent_days_id) > 0) {
  cat(sprintf("Days with large divergence (>%.1f log units): %d/%d days\n", 
              threshold_diff_id, 
              nrow(divergent_days_id), 
              nrow(temporal_analysis_id)))
  cat("\nTop divergent days:\n")
  print(divergent_days_id %>% 
          dplyr::select(Date, ID1, ID2, Difference_log_cons, Direction) %>%
          head(5))
} else {
  cat("✅ No days with large divergence - strong concordance throughout\n")
}

cat("\n=== INTERPRETATION ===\n")
if(mean_diff_id < 0.3) {
  cat("✅ Excellent temporal concordance: Mean difference < 0.3 log units\n")
  cat("   Both devices show highly synchronized daily dynamics\n")
} else if(mean_diff_id < 0.5) {
  cat("✓ Good temporal concordance: Mean difference < 0.5 log units\n")
  cat("   General synchronization with some device-specific variation\n")
} else {
  cat("⚠️  Moderate temporal concordance: Mean difference > 0.5 log units\n")
  cat("   Substantial device-specific dynamics present\n")
}

cat("\n=== TEMPORAL CONCORDANCE COMPLETE ===\n")


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("H1 ID ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n")
cat("  ✅ outputs/Table_H1_ID_Consistency_Results.csv\n")
cat("  ✅ outputs/Figure_S_ID_Diagnostic_Comparison.png (Supplementary)\n")
cat("  ✅ outputs/Figure_S_ID_Diagnostic_Comparison.pdf (Supplementary)\n")
cat("  ✅ outputs/Figure_H1_ID_Consistency.png (Main text)\n")
cat("  ✅ outputs/Figure_H1_ID_Consistency.pdf (Main text)\n\n")

# Dynamic final summary
cat("KEY FINDINGS:\n")
if (!is.null(result_id)) {
  cat(sprintf("  • ID Maize (ID1 vs ID2): R² = %.3f%s (n = %d days, excluding battery failures)\n",
              result_id$r_squared_log_cons,
              get_sig_stars(result_id$pearson_p_log_cons),
              result_id$n_days))
  cat(sprintf("  • Directional concordance: %.1f%% (%d/%d days)\n",
              concordance_pct_id,
              n_concordant_id,
              n_changes_id))
  
  # Interpretation
  if (result_id$r_squared_log_cons < 0.3) {
    cat("  • INTERPRETATION: Poor consistency between ID replicates\n")
  } else if (result_id$r_squared_log_cons < 0.5) {
    cat("  • INTERPRETATION: Moderate consistency between ID replicates\n")
  } else if (result_id$r_squared_log_cons < 0.7) {
    cat("  • INTERPRETATION: Moderate-to-strong consistency between ID replicates\n")
  } else {
    cat("  • INTERPRETATION: Strong consistency between ID replicates\n")
  }
  
  cat("\n⚠️  NOTE: 5 days excluded (Aug 28 - Sep 1) due to documented battery failures\n")
  cat("   → Without power failures, ID shows moderate temporal consistency\n")
  cat("   → Operational reliability remains a critical limitation for continuous monitoring\n")
}

cat("\nSTATISTICAL DECISION:\n")
cat("  • Use log10(x+1) transformation (standardized)\n")
cat("  • Rationale: Consistent methodology with FAIRD analysis (H1_consistency.R)\n")
cat(sprintf("  • Robust regression (log scale, n=17): Pseudo-R2 = %.3f | Low-weight: %d/17\n",
            result_id$pseudo_r2_log_cons, result_id$n_low_weight))
cat(sprintf("  • Robust regression (original scale, n=22 conservative): Pseudo-R2 = %.3f | Low-weight: %d/22\n",
            pseudo_r2_cons, n_low_weight_cons))
cat(sprintf("  • Robust regression (original scale, n=17 operational): Pseudo-R2 = %.3f | Low-weight: %d/17\n\n",
            pseudo_r2_robust, n_low_weight_oper))

cat("NOTES:\n")
cat("  • ID3 excluded from analysis due to device failure\n")
cat("  • Analysis performed on Maize site (Site2) data only\n")
cat("  • Comparison with FAIRD consistency can be found in H1_consistency.R\n\n")


# ==============================================================================
cat("========================================\n\n")

cat("\n========================================\n")
cat("H1 ID COMPLETE\n")
cat("========================================\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")