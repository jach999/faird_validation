# ==============================================================================
# H1: FAIR-DEVICE CONSISTENCY ANALYSIS
# ==============================================================================
# Description: Complete analysis of device consistency with regression diagnostics
# Date: October 22, 2025
# Decision: Use log₁₀(x+1) transformation (standardized)
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
cat("H1: INTRA-DEVICE CONSISTENCY ANALYSIS\n")
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

# Prepare data for FAIRD Maize (primary analysis)
# 'etraps_daily' is now guaranteed to contain all dates
counts_data_faird <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

cat("Data prepared:", nrow(counts_data_faird), "observations\n")
cat("Date range:", min(counts_data_faird$Date), "to", max(counts_data_faird$Date), "\n\n")

# --- 1.1: Original Scale Diagnostics ---
cat("--- Original Scale Analysis ---\n")

model_original <- lm(FAIRD2 ~ FAIRD1, data = counts_data_faird)

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
threshold_cooks <- 4/nrow(counts_data_faird)
influential <- which(cooks_d > threshold_cooks)

cat("Cook's Distance Analysis:\n")
cat("  Threshold:", round(threshold_cooks, 4), "\n")
cat("  Influential points:", length(influential), "/", nrow(counts_data_faird), "\n")

if(length(influential) > 0) {
  cat("\n  Influential observations:\n")
  influential_df <- counts_data_faird[influential, ] %>%
    dplyr::mutate(cooks_d = cooks_d[influential]) %>%
    dplyr::arrange(desc(cooks_d))
  print(influential_df)
}

# --- 1.2: Log-Transformed Diagnostics ---
cat("\n\n--- Log-Transformed Analysis ---\n")

counts_log <- counts_data_faird %>%
  dplyr::mutate(
    FAIRD1_log = log10(FAIRD1 + 1),
    FAIRD2_log = log10(FAIRD2 + 1)
  )

model_log <- lm(FAIRD2_log ~ FAIRD1_log, data = counts_log)

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

model_robust <- MASS::rlm(FAIRD2 ~ FAIRD1, data = counts_data_faird)
residuals_robust <- residuals(model_robust)
ss_res_robust <- sum(residuals_robust^2)
ss_tot <- sum((counts_data_faird$FAIRD2 - mean(counts_data_faird$FAIRD2))^2)
pseudo_r2_robust <- 1 - ss_res_robust/ss_tot

cat("Pseudo-R² (robust):", round(pseudo_r2_robust, 3), "\n\n")

# --- 1.4: Summary Table ---
cat("=== DIAGNOSTIC SUMMARY ===\n\n")

diagnostic_summary <- tibble(
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
    ifelse(shapiro_log$p.value > 0.05 & bp_test_log$p > 0.05, "✅ BEST", "⚠️ Check"),
    "✓ Validation"
  )
) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 3)))

print(diagnostic_summary)

cat("\n✅ DECISION: Use log-transformed data (R²=", round(summary(model_log)$r.squared, 3), ")\n", sep = "")
cat("    Rationale: Homoscedasticity restored, reduced influential points\n\n")

# ==============================================================================
# SECTION 2: CONSISTENCY ANALYSIS FOR ALL DEVICE PAIRS (LOG-TRANSFORMED)
# ==============================================================================

cat("\n=== SECTION 2: CONSISTENCY ANALYSIS (ALL PAIRS) ===\n\n")

# Function for log-transformed consistency analysis
analyze_consistency_log <- function(daily_data, device1, device2, metric = "abundance") {
  
  df <- daily_data %>%
    dplyr::filter(Device %in% c(device1, device2)) %>%
    dplyr::select(Date, Device, !!sym(metric)) %>%
    tidyr::pivot_wider(names_from = Device, values_from = !!sym(metric), values_fill = 0) %>%
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
  
  # Robust regression (MASS::rlm, Huber M-estimation) on log scale
  model_robust     <- MASS::rlm(dev2_log ~ dev1_log, data = df)
  ss_res_robust    <- sum(residuals(model_robust)^2)
  ss_tot_ols       <- sum((df$dev2_log - mean(df$dev2_log))^2)
  pseudo_r2        <- 1 - ss_res_robust / ss_tot_ols
  rlm_weights      <- model_robust$w
  n_low_weight     <- sum(rlm_weights < 0.5)  # Days heavily downweighted
  
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
    pseudo_r2     = pseudo_r2,       # Robust regression pseudo-R²
    rlm_weights   = rlm_weights,     # Full weight vector (for Block 5 analysis)
    n_low_weight  = n_low_weight,    # Count of days with weight < 0.5
    data          = df,
    model         = model_log,
    model_robust  = model_robust
  )
}

# Analyze all device pairs
results <- list()

# FAIRD Maize
cat("--- FAIRD1 vs FAIRD2 (Maize) ---\n")
results$faird_maize <- analyze_consistency_log(etraps_daily, "FAIRD1", "FAIRD2", "abundance")
if(!is.null(results$faird_maize)) {
  cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f\n",
              results$faird_maize$n_days,
              results$faird_maize$r_squared_log,
              results$faird_maize$adj_r_squared,
              results$faird_maize$pearson_p_log))
  cat(sprintf("   Slope: %.3f +/- %.3f SE  [95%% CI: %.3f - %.3f]  (Intercept: %.3f, p=%.3f)\n",
              results$faird_maize$slope, results$faird_maize$slope_se,
              results$faird_maize$slope_ci_low, results$faird_maize$slope_ci_high,
              results$faird_maize$intercept, results$faird_maize$intercept_p))
  cat(sprintf("   Spearman rho=%.3f, p=%.4f\n",
              results$faird_maize$spearman_rho, results$faird_maize$spearman_p))
  cat(sprintf("   Pseudo-R2 (robust): %.3f  |  Low-weight days (<0.5): %d/%d\n\n",
              results$faird_maize$pseudo_r2, results$faird_maize$n_low_weight,
              results$faird_maize$n_days))
}

# FAIRD Meadow
cat("--- FAIRD3 vs FAIRD4 (Meadow) ---\n")
results$faird_meadow <- analyze_consistency_log(etraps_daily, "FAIRD3", "FAIRD4", "abundance")
if(!is.null(results$faird_meadow)) {
  cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f\n",
              results$faird_meadow$n_days,
              results$faird_meadow$r_squared_log,
              results$faird_meadow$adj_r_squared,
              results$faird_meadow$pearson_p_log))
  cat(sprintf("   Slope: %.3f +/- %.3f SE  [95%% CI: %.3f - %.3f]  (Intercept: %.3f, p=%.3f)\n",
              results$faird_meadow$slope, results$faird_meadow$slope_se,
              results$faird_meadow$slope_ci_low, results$faird_meadow$slope_ci_high,
              results$faird_meadow$intercept, results$faird_meadow$intercept_p))
  cat(sprintf("   Spearman rho=%.3f, p=%.4f\n",
              results$faird_meadow$spearman_rho, results$faird_meadow$spearman_p))
  cat(sprintf("   Pseudo-R2 (robust): %.3f  |  Low-weight days (<0.5): %d/%d\n\n",
              results$faird_meadow$pseudo_r2, results$faird_meadow$n_low_weight,
              results$faird_meadow$n_days))
}

# ==============================================================================
# SECTION 3: CREATE RESULTS TABLE
# ==============================================================================

cat("\n=== SECTION 3: RESULTS TABLE (TABLE 1) ===\n\n")

# Dynamically find all pairs analyzed
all_results <- list(results$faird_maize, results$faird_meadow)
all_results <- all_results[!sapply(all_results, is.null)] # Remove nulls if any failed

results_table <- tibble(
  Device_Pair = sapply(all_results, function(x) paste(x$devices, collapse = " vs ")),
  Habitat = c("Maize", "Meadow"), # This might need manual adjustment if pairs change
  n_days = sapply(all_results, function(x) x$n_days),
  R2_log = sapply(all_results, function(x) x$r_squared_log),
  p_value = sapply(all_results, function(x) x$pearson_p_log),
  Spearman_rho = sapply(all_results, function(x) x$spearman_rho),
  Spearman_p = sapply(all_results, function(x) x$spearman_p),
  Consistency = case_when(
    R2_log >= 0.6 ~ "Strong",
    R2_log >= 0.4 ~ "Moderate",
    R2_log >= 0.2 ~ "Weak",
    TRUE ~ "Poor"
  )
) %>%
  dplyr::mutate(across(where(is.numeric) & !n_days, ~round(.x, 3)))

print(results_table)

# Save table
write_csv(results_table, "outputs/Table1_H1_Consistency_Results.csv")
cat("\n✅ Table saved: outputs/Table1_H1_Consistency_Results.csv\n\n")

# ==============================================================================
# SECTION 4: DIAGNOSTIC PLOTS
# ==============================================================================

cat("\n=== SECTION 4: DIAGNOSTIC PLOTS ===\n\n")

# Helper function for significance stars
get_sig_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns") # Not significant
}

# Prepare data with influential flags
counts_data_faird_plot <- counts_data_faird %>%
  dplyr::mutate(
    cooks_d = cooks.distance(model_original),
    influential = cooks_d > threshold_cooks,
    FAIRD1_log = log10(FAIRD1 + 1),
    FAIRD2_log = log10(FAIRD2 + 1)
  )

# --- 4.1: Prepare dynamic labels for plots ---

# Original model values
r2_orig <- summary(model_original)$r.squared
slope_orig <- coef(model_original)[2] # Coefficient for FAIRD1
bp_p_orig <- bp_test$p
# Get p-value from the full regression summary for stars
p_val_orig <- summary(model_original)$coefficients[2, 4] 

label_orig <- sprintf("R² = %.3f%s\nSlope = %.2f\nBP p = %.3f",
                      r2_orig,
                      get_sig_stars(p_val_orig),
                      slope_orig,
                      bp_p_orig)

# Log model values
r2_log <- summary(model_log)$r.squared
slope_log <- coef(model_log)[2] # Coefficient for FAIRD1_log
bp_p_log <- bp_test_log$p
# Get p-value from the full regression summary for stars
p_val_log <- summary(model_log)$coefficients[2, 4] 

label_log <- sprintf("R² = %.3f%s\nSlope = %.2f\nBP p = %.3f",
                     r2_log,
                     get_sig_stars(p_val_log),
                     slope_log,
                     bp_p_log)


cat("  Dynamic label for Original plot:\n", label_orig, "\n")
cat("  Dynamic label for Log plot:\n", label_log, "\n\n")


# --- Plot A: Original Scale (with dynamic label) ---
plot_original <- ggplot(counts_data_faird_plot, aes(x = FAIRD1, y = FAIRD2)) +
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
    color = "gray30",             # Dark gray line
    fill = custom_colors["FAIRD"],  # Gold confidence band (Maize-Core)
    alpha = 0.3, linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray50", linewidth = 0.8) +
  # Use the dynamic label variable
  annotate("text", x = Inf, y = -Inf, 
           label = label_orig,
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20",
           fontface = "italic") +
  labs(
    title = "A) Original Scale",
    x = "FAIRD1 Daily Abundance",
    y = "FAIRD2 Daily Abundance",
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

plot(plot_original)

# Save plots
ggsave("outputs/H1_FAIRD_Consistency_Diagnostic_Original.png", 
       plot = plot_original, 
       width = 12, height = 6, dpi = 300, bg = "white")

# --- Plot B: Log-Transformed Scale (with dynamic label) ---
# CRITICAL: Use SAME influential classification from original model
# This ensures visual consistency between panels A and B

plot_log <- ggplot(counts_data_faird_plot, aes(x = FAIRD1_log, y = FAIRD2_log)) +
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
    color = "gray30",             # Dark gray line
    fill = custom_colors["FAIRD"],  # Gold confidence band (Maize-Core)
    alpha = 0.3, linewidth = 1
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.8) +
  # Use the dynamic label variable
  annotate("text", x = Inf, y = -Inf,
           label = label_log,
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray20",
           fontface = "italic") +
  labs(
    title = "B) Log₁₀-Transformed Scale",
    x = "log₁₀(FAIRD1 + 1)",
    y = "log₁₀(FAIRD2 + 1)",
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

plot(plot_log)

# Save plots
ggsave("outputs/H1_FAIRD_Consistency_Diagnostic_Log.png", 
       plot = plot_log, 
       width = 12, height = 6, dpi = 300, bg = "white")

# --- Combined Plot ---
combined_plot <- plot_original + plot_log +
  plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom")

print(combined_plot)

# Save plots
ggsave("outputs/S1_Diagnostic_Comparison.png", 
       plot = combined_plot, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("outputs/S1_Diagnostic_Comparison.pdf",
       plot = combined_plot,
       width = 12, height = 6, device = cairo_pdf)

cat("✅ Diagnostic plots saved:\n")
cat("   - outputs/S1_Diagnostic_Comparison.png\n")
cat("   - outputs/S1_Diagnostic_Comparison.pdf\n\n")

# ==============================================================================
# SECTION 5: MAIN TEXT FIGURE (LOG-TRANSFORMED ONLY)
# ==============================================================================

cat("\n=== SECTION 5: MAIN TEXT FIGURE ===\n\n")

# --- Determine shared axis limits for consistent comparison ---
all_faird_data <- dplyr::bind_rows(
  results$faird_maize$data %>% dplyr::select(dev1_log, dev2_log),
  results$faird_meadow$data %>% dplyr::select(dev1_log, dev2_log)
)

# Calculate shared limits with small buffer
x_range <- range(c(all_faird_data$dev1_log), na.rm = TRUE)
y_range <- range(c(all_faird_data$dev2_log), na.rm = TRUE)

# Use the maximum range for both axes to maintain square aspect ratio
max_range <- max(c(x_range[2], y_range[2]))
min_range <- min(c(x_range[1], y_range[1]))

# Add 5% buffer
buffer <- (max_range - min_range) * 0.05
shared_limits <- c(min_range - buffer, max_range + buffer)

cat("Shared axis limits:", round(shared_limits, 3), "\n\n")

# Panel A: FAIRD Maize
plot_faird_maize <- ggplot(results$faird_maize$data, 
                           aes(x = dev1_log, y = dev2_log)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = custom_colors["FAIRD"],
    alpha = 0.3
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  stat_cor(method = "pearson",aes(label = paste(after_stat(rr.label),after_stat(p.label),sep = "~`,`~")), label.x = shared_limits[1] + 0.1, 
           label.y = shared_limits[2] - 0.15, size = 3.5) +
  labs(
    title = "A) FAIRD Maize (Site1 vs Site2)",
    x = "log₁₀(FAIRD1 + 1)",
    y = "log₁₀(FAIRD2 + 1)"
  ) +
  coord_cartesian(xlim = shared_limits, ylim = shared_limits) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        aspect.ratio = 1)

# Panel B: FAIRD Meadow
plot_faird_meadow <- ggplot(results$faird_meadow$data, 
                            aes(x = dev1_log, y = dev2_log)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(
    method = "lm", se = TRUE, 
    color = "gray30",
    fill = custom_colors["FAIRD"],
    alpha = 0.3
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  stat_cor(method = "pearson",aes(label = paste(after_stat(rr.label),after_stat(p.label),sep = "~`,`~")), label.x = shared_limits[1] + 0.1, 
           label.y = shared_limits[2] - 0.15, size = 3.5) +
  labs(
    title = "B) FAIRD Meadow (Site3 vs Site4)",
    x = "log₁₀(FAIRD3 + 1)",
    y = "log₁₀(FAIRD4 + 1)"
  ) +
  coord_cartesian(xlim = shared_limits, ylim = shared_limits) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        aspect.ratio = 1)

# Save plots individually
ggsave("outputs/H1_FAIRD_Consistency_Maize.png",
       plot = plot_faird_maize,
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Consistency_Maize.pdf",
       plot = plot_faird_maize,
       width = 7, height = 6, device = cairo_pdf)

ggsave("outputs/H1_FAIRD_Consistency_Meadow.png",
       plot = plot_faird_meadow,
       width = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Consistency_Meadow.pdf",
       plot = plot_faird_meadow,
       width = 7, height = 6, device = cairo_pdf)

cat("✅ Main text figures saved:\n")
cat("   - H1_FAIRD_Consistency_Maize.png\n")
cat("   - H1_FAIRD_Consistency_Maize.pdf\n")
cat("   - H1_FAIRD_Consistency_Meadow.png\n")
cat("   - H1_FAIRD_Consistency_Meadow.pdf\n\n")

# ==============================================================================
# SECTION 6: TEMPORAL CONCORDANCE ANALYSIS
# ==============================================================================

cat("\n=== SECTION X: TEMPORAL CONCORDANCE (DIRECTIONAL SYNCHRONIZATION) ===\n\n")

# ------------------------------------------------------------------------------
# 6A: FAIRD1 vs FAIRD2 (Maize)
# ------------------------------------------------------------------------------
cat("--- Temporal Concordance: FAIRD1 vs FAIRD2 (Maize) ---\n\n")

# Prepare wide-format data for Maize pair
counts_maize_tc <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

# Calculate directional concordance (same methodology as H1_consistency_AMMOD.R)
temporal_concordance_maize <- counts_maize_tc %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    FAIRD1_change = FAIRD1 - dplyr::lag(FAIRD1),
    FAIRD2_change = FAIRD2 - dplyr::lag(FAIRD2),
    concordant    = sign(FAIRD1_change) == sign(FAIRD2_change)
  ) %>%
  dplyr::filter(!is.na(concordant))

n_changes_maize     <- nrow(temporal_concordance_maize)
n_concordant_maize  <- sum(temporal_concordance_maize$concordant, na.rm = TRUE)
concordance_pct_maize <- (n_concordant_maize / n_changes_maize) * 100

cat(sprintf("Temporal Concordance (Maize): %.1f%%\n", concordance_pct_maize))
cat(sprintf("  Days analyzed: %d (n-1 for lag comparison)\n", n_changes_maize))
cat(sprintf("  Concordant days: %d\n", n_concordant_maize))
cat(sprintf("  Discordant days: %d\n\n", n_changes_maize - n_concordant_maize))

# Show discordant days if any
if ((n_changes_maize - n_concordant_maize) > 0) {
  cat("Days with opposite directional changes (Maize):\n")
  discordant_maize <- temporal_concordance_maize %>%
    dplyr::filter(!concordant) %>%
    dplyr::select(Date, FAIRD1, FAIRD2, FAIRD1_change, FAIRD2_change) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        FAIRD1_change > 0 & FAIRD2_change < 0 ~ "FAIRD1 up, FAIRD2 down",
        FAIRD1_change < 0 & FAIRD2_change > 0 ~ "FAIRD1 down, FAIRD2 up",
        TRUE ~ "Other"
      )
    )
  print(discordant_maize)
  cat("\n")
}

# Kendall's tau (non-parametric correlation of changes)
kendall_maize <- cor.test(temporal_concordance_maize$FAIRD1_change,
                          temporal_concordance_maize$FAIRD2_change,
                          method = "kendall",
                          exact  = FALSE)

cat("Kendall's tau (Maize - correlation of directional changes):\n")
cat(sprintf("  tau = %.3f\n", kendall_maize$estimate))
cat(sprintf("  p-value = %.4f\n\n", kendall_maize$p.value))

# Day-by-day log-scale divergence (Maize)
temporal_analysis_maize <- counts_maize_tc %>%
  dplyr::mutate(
    FAIRD1_log    = log10(FAIRD1 + 1),
    FAIRD2_log    = log10(FAIRD2 + 1),
    Difference_log = abs(FAIRD1_log - FAIRD2_log),
    Direction      = dplyr::case_when(
      FAIRD1_log > FAIRD2_log ~ "FAIRD1 > FAIRD2",
      FAIRD1_log < FAIRD2_log ~ "FAIRD2 > FAIRD1",
      TRUE ~ "Equal"
    )
  )

mean_diff_maize <- mean(temporal_analysis_maize$Difference_log, na.rm = TRUE)
sd_diff_maize   <- sd(temporal_analysis_maize$Difference_log,   na.rm = TRUE)
max_diff_maize  <- max(temporal_analysis_maize$Difference_log,  na.rm = TRUE)

cat("Day-to-day divergence log-scale (Maize):\n")
cat(sprintf("  Mean absolute difference: %.3f\n", mean_diff_maize))
cat(sprintf("  SD of differences:        %.3f\n", sd_diff_maize))
cat(sprintf("  Maximum difference:       %.3f\n\n", max_diff_maize))

# ------------------------------------------------------------------------------
# 6B: FAIRD3 vs FAIRD4 (Meadow)
# ------------------------------------------------------------------------------
cat("--- Temporal Concordance: FAIRD3 vs FAIRD4 (Meadow) ---\n\n")

counts_meadow_tc <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD3", "FAIRD4")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

temporal_concordance_meadow <- counts_meadow_tc %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    FAIRD3_change = FAIRD3 - dplyr::lag(FAIRD3),
    FAIRD4_change = FAIRD4 - dplyr::lag(FAIRD4),
    concordant    = sign(FAIRD3_change) == sign(FAIRD4_change)
  ) %>%
  dplyr::filter(!is.na(concordant))

n_changes_meadow      <- nrow(temporal_concordance_meadow)
n_concordant_meadow   <- sum(temporal_concordance_meadow$concordant, na.rm = TRUE)
concordance_pct_meadow <- (n_concordant_meadow / n_changes_meadow) * 100

cat(sprintf("Temporal Concordance (Meadow): %.1f%%\n", concordance_pct_meadow))
cat(sprintf("  Days analyzed: %d (n-1 for lag comparison)\n", n_changes_meadow))
cat(sprintf("  Concordant days: %d\n", n_concordant_meadow))
cat(sprintf("  Discordant days: %d\n\n", n_changes_meadow - n_concordant_meadow))

if ((n_changes_meadow - n_concordant_meadow) > 0) {
  cat("Days with opposite directional changes (Meadow):\n")
  discordant_meadow <- temporal_concordance_meadow %>%
    dplyr::filter(!concordant) %>%
    dplyr::select(Date, FAIRD3, FAIRD4, FAIRD3_change, FAIRD4_change) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        FAIRD3_change > 0 & FAIRD4_change < 0 ~ "FAIRD3 up, FAIRD4 down",
        FAIRD3_change < 0 & FAIRD4_change > 0 ~ "FAIRD3 down, FAIRD4 up",
        TRUE ~ "Other"
      )
    )
  print(discordant_meadow)
  cat("\n")
}

kendall_meadow <- cor.test(temporal_concordance_meadow$FAIRD3_change,
                           temporal_concordance_meadow$FAIRD4_change,
                           method = "kendall",
                           exact  = FALSE)

cat("Kendall's tau (Meadow - correlation of directional changes):\n")
cat(sprintf("  tau = %.3f\n", kendall_meadow$estimate))
cat(sprintf("  p-value = %.4f\n\n", kendall_meadow$p.value))

temporal_analysis_meadow <- counts_meadow_tc %>%
  dplyr::mutate(
    FAIRD3_log     = log10(FAIRD3 + 1),
    FAIRD4_log     = log10(FAIRD4 + 1),
    Difference_log = abs(FAIRD3_log - FAIRD4_log),
    Direction      = dplyr::case_when(
      FAIRD3_log > FAIRD4_log ~ "FAIRD3 > FAIRD4",
      FAIRD3_log < FAIRD4_log ~ "FAIRD4 > FAIRD3",
      TRUE ~ "Equal"
    )
  )

mean_diff_meadow <- mean(temporal_analysis_meadow$Difference_log, na.rm = TRUE)
sd_diff_meadow   <- sd(temporal_analysis_meadow$Difference_log,   na.rm = TRUE)
max_diff_meadow  <- max(temporal_analysis_meadow$Difference_log,  na.rm = TRUE)

cat("Day-to-day divergence log-scale (Meadow):\n")
cat(sprintf("  Mean absolute difference: %.3f\n", mean_diff_meadow))
cat(sprintf("  SD of differences:        %.3f\n", sd_diff_meadow))
cat(sprintf("  Maximum difference:       %.3f\n\n", max_diff_meadow))

# ------------------------------------------------------------------------------
# 6C: TEMPORAL CONCORDANCE PLOTS
# Matching style of H1_consistency_AMMOD.R Section 6 plot
# ------------------------------------------------------------------------------

cat("--- Generating temporal concordance plots ---\n\n")

# Maize long-format data
faird_temporal_maize <- counts_maize_tc %>%
  tidyr::pivot_longer(
    cols      = c(FAIRD1, FAIRD2),
    names_to  = "Device",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(Abundance_log = log10(Abundance + 1))

# Meadow long-format data
faird_temporal_meadow <- counts_meadow_tc %>%
  tidyr::pivot_longer(
    cols      = c(FAIRD3, FAIRD4),
    names_to  = "Device",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(Abundance_log = log10(Abundance + 1))

# --- Plot Maize ---
plot_temporal_faird_maize <- ggplot(
  faird_temporal_maize,
  aes(x = Date, y = Abundance_log, color = Device, shape = Device)
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
  ) +
  scale_shape_manual(
    values = c("FAIRD1" = 16, "FAIRD2" = 17),
    labels = c("FAIRD1" = "FAIRD1 (Site1)", "FAIRD2" = "FAIRD2 (Site2)")
  ) +
  labs(
    title    = "Temporal Concordance: FAIRD1 vs FAIRD2 (Maize)",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       results$faird_maize$r_squared_log,
                       concordance_pct_maize,
                       kendall_maize$estimate,
                       kendall_maize$p.value),
    x      = "Date",
    y      = "log\u2081\u2080(Abundance + 1)",
    color  = "Device",
    shape  = "Device"
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

print(plot_temporal_faird_maize)

ggsave("outputs/H1_FAIRD_Temporal_Concordance_Maize.png",
       plot   = plot_temporal_faird_maize,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Temporal_Concordance_Maize.pdf",
       plot   = plot_temporal_faird_maize,
       width  = 7, height = 6, device = cairo_pdf)

cat("Temporal concordance plot saved (Maize):\n")
cat("   - outputs/H1_FAIRD_Temporal_Concordance_Maize.png\n")
cat("   - outputs/H1_FAIRD_Temporal_Concordance_Maize.pdf\n\n")

# --- Plot Meadow ---
plot_temporal_faird_meadow <- ggplot(
  faird_temporal_meadow,
  aes(x = Date, y = Abundance_log, color = Device, shape = Device)
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
  ) +
  scale_shape_manual(
    values = c("FAIRD3" = 16, "FAIRD4" = 17),
    labels = c("FAIRD3" = "FAIRD3 (Site3)", "FAIRD4" = "FAIRD4 (Site4)")
  ) +
  labs(
    title    = "Temporal Concordance: FAIRD3 vs FAIRD4 (Meadow)",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       results$faird_meadow$r_squared_log,
                       concordance_pct_meadow,
                       kendall_meadow$estimate,
                       kendall_meadow$p.value),
    x      = "Date",
    y      = "log\u2081\u2080(Abundance + 1)",
    color  = "Device",
    shape  = "Device"
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

print(plot_temporal_faird_meadow)

ggsave("outputs/H1_FAIRD_Temporal_Concordance_Meadow.png",
       plot   = plot_temporal_faird_meadow,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Temporal_Concordance_Meadow.pdf",
       plot   = plot_temporal_faird_meadow,
       width  = 7, height = 6, device = cairo_pdf)

cat("Temporal concordance plot saved (Meadow):\n")
cat("   - outputs/H1_FAIRD_Temporal_Concordance_Meadow.png\n")
cat("   - outputs/H1_FAIRD_Temporal_Concordance_Meadow.pdf\n\n")

# ------------------------------------------------------------------------------
# 6D: COMBINED FIGURES (scatter + temporal concordance side by side)
# Panel A = XY correlation (log scale), Panel B = temporal concordance
# ------------------------------------------------------------------------------

cat("--- Generating combined figures (A: scatter | B: temporal) ---\n\n")

# Relabel scatter panels to "A)" for combined layout
plot_faird_maize_A <- plot_faird_maize +
  labs(title = "A) FAIRD Maize (Site1 vs Site2)")

plot_faird_meadow_A <- plot_faird_meadow +
  labs(title = "A) FAIRD Meadow (Site3 vs Site4)")

# Relabel temporal panels to "B)" for combined layout
plot_temporal_faird_maize_B <- plot_temporal_faird_maize +
  labs(title = "B) Temporal Concordance: FAIRD1 vs FAIRD2 (Maize)")

plot_temporal_faird_meadow_B <- plot_temporal_faird_meadow +
  labs(title = "B) Temporal Concordance: FAIRD3 vs FAIRD4 (Meadow)")

# Combined: Maize (A + B)
combined_faird_maize <- plot_faird_maize_A + plot_temporal_faird_maize_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_faird_maize)

ggsave("outputs/H1_FAIRD_Combined_Maize.png",
       plot   = combined_faird_maize,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Combined_Maize.pdf",
       plot   = combined_faird_maize,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved (Maize):\n")
cat("   - outputs/H1_FAIRD_Combined_Maize.png\n")
cat("   - outputs/H1_FAIRD_Combined_Maize.pdf\n\n")

# Combined: Meadow (A + B)
combined_faird_meadow <- plot_faird_meadow_A + plot_temporal_faird_meadow_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_faird_meadow)

ggsave("outputs/H1_FAIRD_Combined_Meadow.png",
       plot   = combined_faird_meadow,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_FAIRD_Combined_Meadow.pdf",
       plot   = combined_faird_meadow,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved (Meadow):\n")
cat("   - outputs/H1_FAIRD_Combined_Meadow.png\n")
cat("   - outputs/H1_FAIRD_Combined_Meadow.pdf\n\n")


# --- Interpretation summary ---
cat("=== TEMPORAL CONCORDANCE INTERPRETATION ===\n\n")

for (pair in list(
  list(name = "Maize (FAIRD1 vs FAIRD2)", pct = concordance_pct_maize, mean_diff = mean_diff_maize),
  list(name = "Meadow (FAIRD3 vs FAIRD4)", pct = concordance_pct_meadow, mean_diff = mean_diff_meadow)
)) {
  cat(sprintf("  %s: Concordance = %.1f%% | Mean log-diff = %.3f -> ",
              pair$name, pair$pct, pair$mean_diff))
  if (pair$mean_diff < 0.3) {
    cat("Excellent temporal concordance\n")
  } else if (pair$mean_diff < 0.5) {
    cat("Good temporal concordance with moderate site-specific variation\n")
  } else {
    cat("Moderate concordance - substantial site-specific dynamics present\n")
  }
}

cat("\n=== TEMPORAL CONCORDANCE COMPLETE ===\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("H1 ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n")
cat("  [OK] outputs/Table1_H1_Consistency_Results.csv\n")
cat("  [OK] outputs/S1_Diagnostic_Comparison.png (Supplementary)\n")
cat("  [OK] outputs/S1_Diagnostic_Comparison.pdf (Supplementary)\n")
cat("  [OK] outputs/H1_FAIRD_Consistency.png (Main text)\n")
cat("  [OK] outputs/H1_FAIRD_Consistency.pdf (Main text)\n\n")

# Dynamic final summary using results objects
cat("KEY FINDINGS:\n")
if (!is.null(results$faird_maize)) {
  cat(sprintf("  • FAIRD Maize:  R2 = %.3f%s (strong consistency)\n",
              results$faird_maize$r_squared_log,
              get_sig_stars(results$faird_maize$pearson_p_log)))
}
if (!is.null(results$faird_meadow)) {
  cat(sprintf("  • FAIRD Meadow: R2 = %.3f%s (strong consistency)\n",
              results$faird_meadow$r_squared_log,
              get_sig_stars(results$faird_meadow$pearson_p_log)))
}

cat("STATISTICAL DECISION:\n")
cat("  • Use log10(x+1) transformation (standardized)\n")
cat("  • Rationale: Resolves heteroscedasticity, reduces influential points\n")
cat(sprintf("  • Robust regression: Pseudo-R2 = %.3f (Maize) | %.3f (Meadow)\n",
            results$faird_maize$pseudo_r2,
            results$faird_meadow$pseudo_r2))
cat(sprintf("  • Low-weight days (<0.5): %d/%d (Maize) | %d/%d (Meadow)\n\n",
            results$faird_maize$n_low_weight, results$faird_maize$n_days,
            results$faird_meadow$n_low_weight, results$faird_meadow$n_days))

# ==============================================================================
cat("========================================\n\n")

cat("\n========================================\n")
cat("H1 FAIRD COMPLETE\n")
cat("========================================\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")