# ==============================================================================
# H1: INSECT DETECT INTRA-DEVICE CONSISTENCY ANALYSIS
# ==============================================================================
# Description: Complete analysis of ID device consistency with regression
#              diagnostics. All variables use _id suffix for traceability.
#              Two datasets compared:
#                - Conservative (n=22): ALL days including battery failures
#                - Operational  (n=17): Excluding battery failures Aug 28-Sep 1
# Date: February 2026
# Decision: Use log10(x+1) transformation (standardized)
# Note: ID3 excluded due to device failure (documented in project context)
# ==============================================================================

# --- Setup ---
library(tidyverse)
library(car)        # For ncvTest (Breusch-Pagan)
library(MASS)       # For robust regression
library(patchwork)  # For combining plots
library(ggpubr)     # For stat_cor

# ============================================================
# CONSOLE OUTPUT CAPTURE
# ============================================================
script_name_id <- sub("\\.R$", "", basename(sys.frame(1)$ofile))

output_file_id <- file.path("outputs",
                            paste0(script_name_id, "_output_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

# Close any previously open sinks to prevent output accumulation
while (sink.number() > 0) sink()
sink(output_file_id, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Script:", paste0(script_name_id, ".R"), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Load custom colors
source("custom_colors.R")
cat("Custom colors loaded\n")

# Load data
cat("\n========================================\n")
cat("H1: INSECT DETECT CONSISTENCY ANALYSIS\n")
cat("========================================\n\n")

source("load_data.R")
load("data_cleaned.RData")

etraps_daily_id <- data_list$etraps_daily

# ==============================================================================
# DATA PREPARATION: TWO DATASETS
# ==============================================================================

cat("\n=== DATA PREPARATION ===\n\n")

# --- CONSERVATIVE dataset (n=22): all days, including battery-failure days ---
counts_cons_id <- etraps_daily_id %>%
  dplyr::filter(Device %in% c("ID1", "ID2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device,
                     values_from = abundance,
                     values_fill = 0) %>%
  dplyr::mutate(
    ID1_log_cons_id = log10(ID1 + 1),
    ID2_log_cons_id = log10(ID2 + 1)
  )

# --- OPERATIONAL dataset (n=17): excluding battery failures Aug 28 - Sep 1 ---
counts_oper_id <- etraps_daily_id %>%
  dplyr::filter(Device %in% c("ID1", "ID2")) %>%
  dplyr::filter(Date < as.Date("2023-08-28") | Date > as.Date("2023-09-01")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device,
                     values_from = abundance,
                     values_fill = 0) %>%
  dplyr::mutate(
    ID1_log_oper_id = log10(ID1 + 1),
    ID2_log_oper_id = log10(ID2 + 1)
  )

cat("Conservative dataset (n=22, all days):", nrow(counts_cons_id), "observations\n")
cat("  Date range:", format(min(counts_cons_id$Date)), "to",
    format(max(counts_cons_id$Date)), "\n\n")

cat("Operational dataset (n=17, excl. battery failures):",
    nrow(counts_oper_id), "observations\n")
cat("  Date range:", format(min(counts_oper_id$Date)), "to",
    format(max(counts_oper_id$Date)), "\n")
cat("  Excluded: 2023-08-28 to 2023-09-01 (documented battery failures)\n")
cat("  ID3 excluded from all analyses (device failure)\n\n")

# ==============================================================================
# SECTION 1: FOUR LINEAR MODELS WITH REGRESSION DIAGNOSTICS
# ==============================================================================

cat("\n=== SECTION 1: REGRESSION DIAGNOSTICS (ALL FOUR MODELS) ===\n\n")

# --- Cook's Distance threshold (based on operational n=17 for comparability) ---
threshold_cooks_id <- 4 / nrow(counts_oper_id)
cat("Cook's Distance threshold (4/n, n=17):", round(threshold_cooks_id, 4), "\n\n")

# ------------------------------------------------------------------------------
# MODEL 1: Conservative — Original scale (n=22)
# ------------------------------------------------------------------------------
cat("--- Model 1: Conservative Original Scale (n=22) ---\n")

model_original_cons_id <- lm(ID2 ~ ID1, data = counts_cons_id)
shapiro_orig_cons_id   <- shapiro.test(residuals(model_original_cons_id))
bp_orig_cons_id        <- car::ncvTest(model_original_cons_id)
cooks_d_orig_cons_id   <- cooks.distance(model_original_cons_id)
threshold_cooks_cons_id <- 4 / nrow(counts_cons_id)
n_influential_orig_cons_id <- sum(cooks_d_orig_cons_id > threshold_cooks_cons_id)

cat(sprintf("R2 = %.3f | Adj.R2 = %.3f\n",
            summary(model_original_cons_id)$r.squared,
            summary(model_original_cons_id)$adj.r.squared))
cat(sprintf("Shapiro-Wilk p = %.4f %s\n",
            shapiro_orig_cons_id$p.value,
            ifelse(shapiro_orig_cons_id$p.value > 0.05, "[Normal]", "[NOT normal]")))
cat(sprintf("Breusch-Pagan p = %.4f %s\n",
            bp_orig_cons_id$p,
            ifelse(bp_orig_cons_id$p > 0.05, "[Homoscedastic]", "[Heteroscedastic]")))
cat(sprintf("Influential points (Cook's D > %.4f): %d/%d\n\n",
            threshold_cooks_cons_id, n_influential_orig_cons_id, nrow(counts_cons_id)))

# Robust regression - conservative original
model_robust_orig_cons_id <- MASS::rlm(ID2 ~ ID1, data = counts_cons_id)
ss_res_rcons_id           <- sum(residuals(model_robust_orig_cons_id)^2)
ss_tot_rcons_id           <- sum((counts_cons_id$ID2 - mean(counts_cons_id$ID2))^2)
pseudo_r2_orig_cons_id    <- 1 - ss_res_rcons_id / ss_tot_rcons_id
n_low_weight_orig_cons_id <- sum(model_robust_orig_cons_id$w < 0.5)

# ------------------------------------------------------------------------------
# MODEL 2: Conservative — Log scale (n=22)
# ------------------------------------------------------------------------------
cat("--- Model 2: Conservative Log Scale (n=22) ---\n")

model_log_cons_id      <- lm(ID2_log_cons_id ~ ID1_log_cons_id, data = counts_cons_id)
shapiro_log_cons_id    <- shapiro.test(residuals(model_log_cons_id))
bp_log_cons_id         <- car::ncvTest(model_log_cons_id)
cooks_d_log_cons_id    <- cooks.distance(model_log_cons_id)
n_influential_log_cons_id <- sum(cooks_d_log_cons_id > threshold_cooks_cons_id)

cat(sprintf("R2 = %.3f | Adj.R2 = %.3f\n",
            summary(model_log_cons_id)$r.squared,
            summary(model_log_cons_id)$adj.r.squared))
cat(sprintf("Shapiro-Wilk p = %.4f %s\n",
            shapiro_log_cons_id$p.value,
            ifelse(shapiro_log_cons_id$p.value > 0.05, "[Normal]", "[NOT normal]")))
cat(sprintf("Breusch-Pagan p = %.4f %s\n",
            bp_log_cons_id$p,
            ifelse(bp_log_cons_id$p > 0.05, "[Homoscedastic]", "[Heteroscedastic]")))
cat(sprintf("Influential points (Cook's D > %.4f): %d/%d\n\n",
            threshold_cooks_cons_id, n_influential_log_cons_id, nrow(counts_cons_id)))

# Robust regression - conservative log
model_robust_log_cons_id  <- MASS::rlm(ID2_log_cons_id ~ ID1_log_cons_id,
                                       data = counts_cons_id)
ss_res_rlcons_id           <- sum(residuals(model_robust_log_cons_id)^2)
ss_tot_rlcons_id           <- sum((counts_cons_id$ID2_log_cons_id -
                                     mean(counts_cons_id$ID2_log_cons_id))^2)
pseudo_r2_log_cons_id      <- 1 - ss_res_rlcons_id / ss_tot_rlcons_id
n_low_weight_log_cons_id   <- sum(model_robust_log_cons_id$w < 0.5)

# ------------------------------------------------------------------------------
# MODEL 3: Operational — Original scale (n=17)
# ------------------------------------------------------------------------------
cat("--- Model 3: Operational Original Scale (n=17) ---\n")

model_original_oper_id    <- lm(ID2 ~ ID1, data = counts_oper_id)
shapiro_orig_oper_id      <- shapiro.test(residuals(model_original_oper_id))
bp_orig_oper_id           <- car::ncvTest(model_original_oper_id)
cooks_d_orig_oper_id      <- cooks.distance(model_original_oper_id)
n_influential_orig_oper_id <- sum(cooks_d_orig_oper_id > threshold_cooks_id)

cat(sprintf("R2 = %.3f | Adj.R2 = %.3f\n",
            summary(model_original_oper_id)$r.squared,
            summary(model_original_oper_id)$adj.r.squared))
cat(sprintf("Shapiro-Wilk p = %.4f %s\n",
            shapiro_orig_oper_id$p.value,
            ifelse(shapiro_orig_oper_id$p.value > 0.05, "[Normal]", "[NOT normal]")))
cat(sprintf("Breusch-Pagan p = %.4f %s\n",
            bp_orig_oper_id$p,
            ifelse(bp_orig_oper_id$p > 0.05, "[Homoscedastic]", "[Heteroscedastic]")))
cat(sprintf("Influential points (Cook's D > %.4f): %d/%d\n\n",
            threshold_cooks_id, n_influential_orig_oper_id, nrow(counts_oper_id)))

# Robust regression - operational original
model_robust_orig_oper_id  <- MASS::rlm(ID2 ~ ID1, data = counts_oper_id)
ss_res_roper_id            <- sum(residuals(model_robust_orig_oper_id)^2)
ss_tot_roper_id            <- sum((counts_oper_id$ID2 - mean(counts_oper_id$ID2))^2)
pseudo_r2_orig_oper_id     <- 1 - ss_res_roper_id / ss_tot_roper_id
n_low_weight_orig_oper_id  <- sum(model_robust_orig_oper_id$w < 0.5)

# Flag influential points for plotting (operational original scale)
influential_orig_oper_id <- which(cooks_d_orig_oper_id > threshold_cooks_id)
cat("Influential observations (operational, original scale):\n")
if (length(influential_orig_oper_id) > 0) {
  print(counts_oper_id[influential_orig_oper_id, ] %>%
          dplyr::mutate(cooks_d = cooks_d_orig_oper_id[influential_orig_oper_id]) %>%
          dplyr::arrange(desc(cooks_d)))
} else {
  cat("  None\n")
}
cat("\n")

# ------------------------------------------------------------------------------
# MODEL 4: Operational — Log scale (n=17)
# ------------------------------------------------------------------------------
cat("--- Model 4: Operational Log Scale (n=17) ---\n")

model_log_oper_id         <- lm(ID2_log_oper_id ~ ID1_log_oper_id, data = counts_oper_id)
shapiro_log_oper_id       <- shapiro.test(residuals(model_log_oper_id))
bp_log_oper_id            <- car::ncvTest(model_log_oper_id)
cooks_d_log_oper_id       <- cooks.distance(model_log_oper_id)
n_influential_log_oper_id <- sum(cooks_d_log_oper_id > threshold_cooks_id)

cat(sprintf("R2 = %.3f | Adj.R2 = %.3f\n",
            summary(model_log_oper_id)$r.squared,
            summary(model_log_oper_id)$adj.r.squared))
cat(sprintf("Shapiro-Wilk p = %.4f %s\n",
            shapiro_log_oper_id$p.value,
            ifelse(shapiro_log_oper_id$p.value > 0.05, "[Normal]", "[NOT normal]")))
cat(sprintf("Breusch-Pagan p = %.4f %s\n",
            bp_log_oper_id$p,
            ifelse(bp_log_oper_id$p > 0.05, "[Homoscedastic]", "[Heteroscedastic]")))
cat(sprintf("Influential points (Cook's D > %.4f): %d/%d\n\n",
            threshold_cooks_id, n_influential_log_oper_id, nrow(counts_oper_id)))

# Robust regression - operational log
model_robust_log_oper_id  <- MASS::rlm(ID2_log_oper_id ~ ID1_log_oper_id,
                                       data = counts_oper_id)
ss_res_rloper_id           <- sum(residuals(model_robust_log_oper_id)^2)
ss_tot_rloper_id           <- sum((counts_oper_id$ID2_log_oper_id -
                                     mean(counts_oper_id$ID2_log_oper_id))^2)
pseudo_r2_log_oper_id      <- 1 - ss_res_rloper_id / ss_tot_rloper_id
n_low_weight_log_oper_id   <- sum(model_robust_log_oper_id$w < 0.5)

# Flag influential points for plotting (operational log scale)
influential_log_oper_id <- which(cooks_d_log_oper_id > threshold_cooks_id)

# ==============================================================================
# SECTION 2: DETAILED CONSISTENCY METRICS (FUNCTION + RESULTS)
# ==============================================================================

cat("\n=== SECTION 2: CONSISTENCY METRICS (LOG-TRANSFORMED) ===\n\n")

# Function: takes a pre-filtered WIDE dataset (device1 and device2 as columns)
# Returns: full list of metrics
analyze_consistency_log_id <- function(wide_data, device1, device2) {
  
  df_id <- wide_data %>%
    dplyr::mutate(
      dev1_log = log10(!!sym(device1) + 1),
      dev2_log = log10(!!sym(device2) + 1)
    )
  
  if (nrow(df_id) < 3) {
    warning(sprintf("Insufficient data for %s vs %s", device1, device2))
    return(NULL)
  }
  
  # Pearson on log-transformed
  cor_pearson_log_id <- cor.test(df_id$dev1_log, df_id$dev2_log, method = "pearson")
  
  # Spearman on original (non-parametric verification)
  cor_spearman_id <- cor.test(df_id[[device1]], df_id[[device2]], method = "spearman")
  
  # OLS regression on log
  model_ols_id  <- lm(dev2_log ~ dev1_log, data = df_id)
  r2_log_id     <- summary(model_ols_id)$r.squared
  adj_r2_id     <- summary(model_ols_id)$adj.r.squared
  
  # Slope and 95% CI
  coef_id       <- coef(summary(model_ols_id))
  slope_id      <- coef_id["dev1_log", "Estimate"]
  slope_se_id   <- coef_id["dev1_log", "Std. Error"]
  slope_ci_id   <- confint(model_ols_id)["dev1_log", ]
  intercept_id  <- coef_id["(Intercept)", "Estimate"]
  intercept_p_id <- coef_id["(Intercept)", "Pr(>|t|)"]
  
  # Robust regression on log scale (Huber M-estimation)
  model_rlm_id   <- MASS::rlm(dev2_log ~ dev1_log, data = df_id)
  ss_res_id      <- sum(residuals(model_rlm_id)^2)
  ss_tot_id      <- sum((df_id$dev2_log - mean(df_id$dev2_log))^2)
  pseudo_r2_id   <- 1 - ss_res_id / ss_tot_id
  n_low_wt_id    <- sum(model_rlm_id$w < 0.5)
  
  list(
    devices        = c(device1, device2),
    n_days         = nrow(df_id),
    r2_log         = r2_log_id,
    adj_r2         = adj_r2_id,
    pearson_r      = as.numeric(cor_pearson_log_id$estimate),
    pearson_p      = cor_pearson_log_id$p.value,
    spearman_rho   = as.numeric(cor_spearman_id$estimate),
    spearman_p     = cor_spearman_id$p.value,
    slope          = slope_id,
    slope_se       = slope_se_id,
    slope_ci_low   = slope_ci_id[1],
    slope_ci_high  = slope_ci_id[2],
    intercept      = intercept_id,
    intercept_p    = intercept_p_id,
    pseudo_r2      = pseudo_r2_id,
    n_low_weight   = n_low_wt_id,
    data           = df_id,
    model          = model_ols_id,
    model_robust   = model_rlm_id
  )
}

# Helper function for significance stars
get_sig_stars_id <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("")
}

# --- Run analysis for both datasets ---
cat("--- Consistency metrics: Conservative (n=22) ---\n")
result_cons_id <- analyze_consistency_log_id(counts_cons_id, "ID1", "ID2")

if (!is.null(result_cons_id)) {
  cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f%s\n",
              result_cons_id$n_days, result_cons_id$r2_log,
              result_cons_id$adj_r2, result_cons_id$pearson_p,
              get_sig_stars_id(result_cons_id$pearson_p)))
  cat(sprintf("  Slope: %.3f +/- %.3f SE [95%% CI: %.3f - %.3f] | Intercept: %.3f (p=%.3f)\n",
              result_cons_id$slope, result_cons_id$slope_se,
              result_cons_id$slope_ci_low, result_cons_id$slope_ci_high,
              result_cons_id$intercept, result_cons_id$intercept_p))
  cat(sprintf("  Spearman rho=%.3f, p=%.4f\n",
              result_cons_id$spearman_rho, result_cons_id$spearman_p))
  cat(sprintf("  Pseudo-R2 (robust, log): %.3f | Low-weight days (<0.5): %d/%d\n\n",
              result_cons_id$pseudo_r2, result_cons_id$n_low_weight,
              result_cons_id$n_days))
}

cat("--- Consistency metrics: Operational (n=17) ---\n")
result_oper_id <- analyze_consistency_log_id(counts_oper_id, "ID1", "ID2")

if (!is.null(result_oper_id)) {
  cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f%s\n",
              result_oper_id$n_days, result_oper_id$r2_log,
              result_oper_id$adj_r2, result_oper_id$pearson_p,
              get_sig_stars_id(result_oper_id$pearson_p)))
  cat(sprintf("  Slope: %.3f +/- %.3f SE [95%% CI: %.3f - %.3f] | Intercept: %.3f (p=%.3f)\n",
              result_oper_id$slope, result_oper_id$slope_se,
              result_oper_id$slope_ci_low, result_oper_id$slope_ci_high,
              result_oper_id$intercept, result_oper_id$intercept_p))
  cat(sprintf("  Spearman rho=%.3f, p=%.4f\n",
              result_oper_id$spearman_rho, result_oper_id$spearman_p))
  cat(sprintf("  Pseudo-R2 (robust, log): %.3f | Low-weight days (<0.5): %d/%d\n\n",
              result_oper_id$pseudo_r2, result_oper_id$n_low_weight,
              result_oper_id$n_days))
}

# ==============================================================================
# SECTION 3: COMPREHENSIVE RESULTS TABLE (ALL FOUR MODELS)
# ==============================================================================

cat("\n=== SECTION 3: COMPREHENSIVE RESULTS TABLE ===\n\n")

results_table_id <- tibble::tibble(
  Dataset     = c("Conservative", "Conservative", "Operational", "Operational"),
  Scale       = c("Original", "Log10(x+1)", "Original", "Log10(x+1)"),
  N_days      = c(nrow(counts_cons_id), nrow(counts_cons_id),
                  nrow(counts_oper_id),  nrow(counts_oper_id)),
  R2          = c(
    summary(model_original_cons_id)$r.squared,
    summary(model_log_cons_id)$r.squared,
    summary(model_original_oper_id)$r.squared,
    summary(model_log_oper_id)$r.squared
  ),
  Adj_R2      = c(
    summary(model_original_cons_id)$adj.r.squared,
    summary(model_log_cons_id)$adj.r.squared,
    summary(model_original_oper_id)$adj.r.squared,
    summary(model_log_oper_id)$adj.r.squared
  ),
  Pearson_p   = c(
    summary(model_original_cons_id)$coefficients[2, 4],
    summary(model_log_cons_id)$coefficients[2, 4],
    summary(model_original_oper_id)$coefficients[2, 4],
    summary(model_log_oper_id)$coefficients[2, 4]
  ),
  Shapiro_p   = c(
    shapiro_orig_cons_id$p.value,
    shapiro_log_cons_id$p.value,
    shapiro_orig_oper_id$p.value,
    shapiro_log_oper_id$p.value
  ),
  BP_p        = c(
    bp_orig_cons_id$p,
    bp_log_cons_id$p,
    bp_orig_oper_id$p,
    bp_log_oper_id$p
  ),
  N_influential = c(
    n_influential_orig_cons_id,
    n_influential_log_cons_id,
    n_influential_orig_oper_id,
    n_influential_log_oper_id
  ),
  Pseudo_R2_robust = c(
    pseudo_r2_orig_cons_id,
    pseudo_r2_log_cons_id,
    pseudo_r2_orig_oper_id,
    pseudo_r2_log_oper_id
  ),
  N_low_weight = c(
    n_low_weight_orig_cons_id,
    n_low_weight_log_cons_id,
    n_low_weight_orig_oper_id,
    n_low_weight_log_oper_id
  ),
  Decision = c(
    ifelse(bp_orig_cons_id$p < 0.05 | n_influential_orig_cons_id > 2,
           "Issues", "OK"),
    ifelse(shapiro_log_cons_id$p.value > 0.05 & bp_log_cons_id$p > 0.05,
           "USE THIS", "Check"),
    ifelse(bp_orig_oper_id$p < 0.05 | n_influential_orig_oper_id > 2,
           "Issues", "OK"),
    ifelse(shapiro_log_oper_id$p.value > 0.05 & bp_log_oper_id$p > 0.05,
           "USE THIS", "Check")
  )
) %>%
  dplyr::mutate(across(where(is.numeric) & !N_days & !N_influential & !N_low_weight,
                       ~round(.x, 4)))

cat("Four-model comparison table (ID1 vs ID2):\n\n")
print(results_table_id)

readr::write_csv(results_table_id,
                 "outputs/Table_H1_ID_AllModels_Results.csv")
cat("\nResults exported to: outputs/Table_H1_ID_AllModels_Results.csv\n\n")

# ==============================================================================
# SECTION 4: DIAGNOSTIC PLOTS — OPERATIONAL DATA (SUPPLEMENTARY)
# Two panels: original scale + log scale, points colored by Cook's D influence
# ==============================================================================

cat("\n=== SECTION 4: DIAGNOSTIC COMPARISON PLOTS (OPERATIONAL, n=17) ===\n\n")

# Prepare plot data with influential flags (operational, original scale Cook's D)
counts_oper_plot_id <- counts_oper_id %>%
  dplyr::mutate(
    influential_orig_id = cooks_d_orig_oper_id > threshold_cooks_id,
    influential_log_id  = cooks_d_log_oper_id  > threshold_cooks_id
  )

label_orig_oper_id <- sprintf(
  "R2 = %.3f\nShapiro p = %.3f\nBP p = %.3f",
  summary(model_original_oper_id)$r.squared,
  shapiro_orig_oper_id$p.value,
  bp_orig_oper_id$p
)

label_log_oper_id <- sprintf(
  "R2 = %.3f\nShapiro p = %.3f\nBP p = %.3f",
  summary(model_log_oper_id)$r.squared,
  shapiro_log_oper_id$p.value,
  bp_log_oper_id$p
)

# --- Plot A: Original Scale (operational) ---
plot_diag_orig_oper_id <- ggplot(
  counts_oper_plot_id,
  aes(x = ID1, y = ID2)
) +
  geom_point(aes(color = influential_orig_id, size = influential_orig_id),
             alpha = 0.7) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray30"),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  scale_size_manual(
    values = c("TRUE" = 4, "FALSE" = 3),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  geom_smooth(method = "lm", se = TRUE,
              color = "gray30",
              fill  = custom_colors["ID"],
              alpha = 0.3, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  annotate("text", x = Inf, y = -Inf,
           label = label_orig_oper_id,
           hjust = 1.1, vjust = -0.5, size = 3.5,
           color = "gray20", fontface = "italic") +
  labs(
    title    = "A) Original Scale — Operational (n=17)",
    subtitle = "Excluding battery-failure days (Aug 28 - Sep 1)",
    x        = "ID1 Daily Abundance",
    y        = "ID2 Daily Abundance",
    color    = "Point type",
    size     = "Point type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# --- Plot B: Log Scale (operational) ---
plot_diag_log_oper_id <- ggplot(
  counts_oper_plot_id,
  aes(x = ID1_log_oper_id, y = ID2_log_oper_id)
) +
  geom_point(aes(color = influential_log_id, size = influential_log_id),
             alpha = 0.7) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray30"),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  scale_size_manual(
    values = c("TRUE" = 4, "FALSE" = 3),
    labels = c("TRUE" = "Influential", "FALSE" = "Normal")
  ) +
  geom_smooth(method = "lm", se = TRUE,
              color = "gray30",
              fill  = custom_colors["ID"],
              alpha = 0.3, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  annotate("text", x = Inf, y = -Inf,
           label = label_log_oper_id,
           hjust = 1.1, vjust = -0.5, size = 3.5,
           color = "gray20", fontface = "italic") +
  labs(
    title    = "B) Log10 Scale — Operational (n=17)",
    subtitle = "Excluding battery-failure days (Aug 28 - Sep 1)",
    x        = "log10(ID1 + 1)",
    y        = "log10(ID2 + 1)",
    color    = "Point type",
    size     = "Point type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# --- Combined diagnostic plot ---
combined_diag_id <- plot_diag_orig_oper_id + plot_diag_log_oper_id +
  patchwork::plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom")

print(combined_diag_id)

ggsave("outputs/H1_ID_Diagnostic_Comparison.png",
       plot   = combined_diag_id,
       width  = 12, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ID_Diagnostic_Comparison.pdf",
       plot   = combined_diag_id,
       width  = 12, height = 6, device = cairo_pdf)

cat("Supplementary diagnostic plots saved:\n")
cat("   - outputs/H1_ID_Diagnostic_Comparison.png\n")
cat("   - outputs/H1_ID_Diagnostic_Comparison.pdf\n\n")

# ==============================================================================
# SECTION 5: MAIN TEXT FIGURES — SCATTER PLOTS (CLEAN, WITHOUT OUTLIER COLORING)
# A) Conservative log (n=22) | B) Operational log (n=17)
# ==============================================================================

cat("\n=== SECTION 5: MAIN TEXT SCATTER PLOTS ===\n\n")

# Shared axis limits across both datasets for visual comparability
all_log_vals_id <- c(
  result_cons_id$data$dev1_log,
  result_cons_id$data$dev2_log,
  result_oper_id$data$dev1_log,
  result_oper_id$data$dev2_log
)
axis_min_id <- min(all_log_vals_id, na.rm = TRUE)
axis_max_id <- max(all_log_vals_id, na.rm = TRUE)
buffer_id   <- (axis_max_id - axis_min_id) * 0.05
limits_id   <- c(axis_min_id - buffer_id, axis_max_id + buffer_id)

cat(sprintf("Shared axis limits: [%.3f, %.3f]\n\n", limits_id[1], limits_id[2]))

# --- Plot A: Conservative log (n=22) ---
plot_scatter_cons_id <- ggplot(
  result_cons_id$data,
  aes(x = dev1_log, y = dev2_log)
) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(method = "lm", se = TRUE,
              color = "gray30",
              fill  = custom_colors["ID"],
              alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50") +
  ggpubr::stat_cor(method = "pearson",
                   aes(label = paste(after_stat(rr.label),
                                     after_stat(p.label),
                                     sep = "~`,`~")),
                   label.x = limits_id[1] + 0.05,
                   label.y = limits_id[2] - 0.1,
                   size = 3.5) +
  coord_cartesian(xlim = limits_id, ylim = limits_id) +
  labs(
    title    = "ID Consistency — Conservative (n=22)",
    subtitle = "All days including battery-failure period",
    x        = "log10(ID1 + 1)",
    y        = "log10(ID2 + 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    panel.border  = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    aspect.ratio  = 1
  )

print(plot_scatter_cons_id)

ggsave("outputs/H1_ID_Consistency_Conservative.png",
       plot   = plot_scatter_cons_id,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ID_Consistency_Conservative.pdf",
       plot   = plot_scatter_cons_id,
       width  = 7, height = 6, device = cairo_pdf)

cat("Conservative scatter plot saved:\n")
cat("   - outputs/H1_ID_Consistency_Conservative.png\n")
cat("   - outputs/H1_ID_Consistency_Conservative.pdf\n\n")

# --- Plot B: Operational log (n=17) ---
plot_scatter_oper_id <- ggplot(
  result_oper_id$data,
  aes(x = dev1_log, y = dev2_log)
) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(method = "lm", se = TRUE,
              color = "gray30",
              fill  = custom_colors["ID"],
              alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50") +
  ggpubr::stat_cor(method = "pearson",
                   aes(label = paste(after_stat(rr.label),
                                     after_stat(p.label),
                                     sep = "~`,`~")),
                   label.x = limits_id[1] + 0.05,
                   label.y = limits_id[2] - 0.1,
                   size = 3.5) +
  coord_cartesian(xlim = limits_id, ylim = limits_id) +
  labs(
    title    = "ID Consistency — Operational (n=17)",
    subtitle = "Battery-failure days excluded (Aug 28 - Sep 1)",
    x        = "log10(ID1 + 1)",
    y        = "log10(ID2 + 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    panel.border  = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    aspect.ratio  = 1
  )

print(plot_scatter_oper_id)

ggsave("outputs/H1_ID_Consistency_Operational.png",
       plot   = plot_scatter_oper_id,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ID_Consistency_Operational.pdf",
       plot   = plot_scatter_oper_id,
       width  = 7, height = 6, device = cairo_pdf)

cat("Operational scatter plot saved:\n")
cat("   - outputs/H1_ID_Consistency_Operational.png\n")
cat("   - outputs/H1_ID_Consistency_Operational.pdf\n\n")

# ==============================================================================
# SECTION 6: TEMPORAL CONCORDANCE PLOT AND ANALYSIS
# Based on OPERATIONAL dataset (n=17) — consistent with main analysis decision
# ==============================================================================

cat("\n=== SECTION 6: TEMPORAL CONCORDANCE (OPERATIONAL, n=17) ===\n\n")

# --- Directional synchronization ---
temporal_concordance_id <- counts_oper_id %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    ID1_change_id = ID1 - dplyr::lag(ID1),
    ID2_change_id = ID2 - dplyr::lag(ID2),
    concordant_id = sign(ID1_change_id) == sign(ID2_change_id)
  ) %>%
  dplyr::filter(!is.na(concordant_id))

n_changes_id     <- nrow(temporal_concordance_id)
n_concordant_id  <- sum(temporal_concordance_id$concordant_id, na.rm = TRUE)
concordance_pct_id <- (n_concordant_id / n_changes_id) * 100

cat(sprintf("Temporal Concordance (operational): %.1f%%\n", concordance_pct_id))
cat(sprintf("  Days analyzed: %d (n-1 lag)\n", n_changes_id))
cat(sprintf("  Concordant days: %d | Discordant: %d\n\n",
            n_concordant_id, n_changes_id - n_concordant_id))

if ((n_changes_id - n_concordant_id) > 0) {
  cat("Days with opposite directional changes:\n")
  discordant_id <- temporal_concordance_id %>%
    dplyr::filter(!concordant_id) %>%
    dplyr::select(Date, ID1, ID2, ID1_change_id, ID2_change_id) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        ID1_change_id > 0 & ID2_change_id < 0 ~ "ID1 up, ID2 down",
        ID1_change_id < 0 & ID2_change_id > 0 ~ "ID1 down, ID2 up",
        TRUE ~ "Other"
      )
    )
  print(discordant_id)
  cat("\n")
}

kendall_test_id <- cor.test(temporal_concordance_id$ID1_change_id,
                            temporal_concordance_id$ID2_change_id,
                            method = "kendall",
                            exact  = FALSE)

cat("Kendall's tau (correlation of directional changes):\n")
cat(sprintf("  tau = %.3f | p = %.4f\n\n",
            kendall_test_id$estimate, kendall_test_id$p.value))

# Day-by-day log-scale divergence
temporal_analysis_id <- counts_oper_id %>%
  dplyr::mutate(
    Difference_log_id = abs(ID1_log_oper_id - ID2_log_oper_id),
    Direction_id = dplyr::case_when(
      ID1_log_oper_id > ID2_log_oper_id ~ "ID1 > ID2",
      ID1_log_oper_id < ID2_log_oper_id ~ "ID2 > ID1",
      TRUE ~ "Equal"
    )
  )

mean_diff_id <- mean(temporal_analysis_id$Difference_log_id, na.rm = TRUE)
sd_diff_id   <- sd(temporal_analysis_id$Difference_log_id,   na.rm = TRUE)
max_diff_id  <- max(temporal_analysis_id$Difference_log_id,  na.rm = TRUE)

cat("Day-to-day divergence (log-scale, operational):\n")
cat(sprintf("  Mean absolute difference: %.3f\n", mean_diff_id))
cat(sprintf("  SD of differences:        %.3f\n", sd_diff_id))
cat(sprintf("  Maximum difference:       %.3f\n\n", max_diff_id))

divergent_days_id <- temporal_analysis_id %>%
  dplyr::filter(Difference_log_id > 0.5) %>%
  dplyr::arrange(desc(Difference_log_id))

if (nrow(divergent_days_id) > 0) {
  cat(sprintf("Days with large divergence (>0.5 log units): %d/%d\n",
              nrow(divergent_days_id), nrow(temporal_analysis_id)))
  print(divergent_days_id %>%
          dplyr::select(Date, ID1, ID2, Difference_log_id, Direction_id) %>%
          head(5))
} else {
  cat("No days with large divergence (> 0.5 log units)\n")
}

cat("\n=== INTERPRETATION ===\n")
if (mean_diff_id < 0.3) {
  cat("Excellent temporal concordance: Mean difference < 0.3 log units\n")
} else if (mean_diff_id < 0.5) {
  cat("Good temporal concordance: Mean difference < 0.5 log units\n")
} else {
  cat("Moderate temporal concordance: Mean difference > 0.5 log units\n")
}
cat("=== TEMPORAL CONCORDANCE COMPLETE ===\n")


# --- Long format for temporal plot ---
id_temporal_id <- counts_oper_id %>%
  dplyr::select(Date, ID1, ID2) %>%
  tidyr::pivot_longer(
    cols      = c(ID1, ID2),
    names_to  = "Device",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(Abundance_log = log10(Abundance + 1))

# --- Temporal plot (log-transformed, operational) ---
plot_temporal_id <- ggplot(
  id_temporal_id,
  aes(x = Date, y = Abundance_log, color = Device, shape = Device)
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
    labels = c("ID1" = "ID1 (Site2 - Maize)",
               "ID2" = "ID2 (Site2 - Maize)")
  ) +
  scale_shape_manual(
    values = c("ID1" = 16, "ID2" = 17),
    labels = c("ID1" = "ID1 (Site2 - Maize)",
               "ID2" = "ID2 (Site2 - Maize)")
  ) +
  labs(
    title    = "Temporal Concordance: ID1 vs ID2 (Operational, n=17)",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       summary(model_log_oper_id)$r.squared,
                       concordance_pct_id,
                       kendall_test_id$estimate,
                       kendall_test_id$p.value),
    x        = "Date",
    y        = "log10(Daily Abundance + 1)",
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

print(plot_temporal_id)

ggsave("outputs/H1_ID_Temporal_Concordance.png",
       plot   = plot_temporal_id,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ID_Temporal_Concordance.pdf",
       plot   = plot_temporal_id,
       width  = 7, height = 6, device = cairo_pdf)

cat("Temporal concordance plot saved:\n")
cat("   - outputs/H1_ID_Temporal_Concordance.png\n")
cat("   - outputs/H1_ID_Temporal_Concordance.pdf\n\n")



# ==============================================================================
# SECTION 6B: TEMPORAL CONCORDANCE — CONSERVATIVE DATASET (n=22)
# Calculates the missing concordance values for PCV documentation
# ==============================================================================

cat("\n=== SECTION 6B: TEMPORAL CONCORDANCE (CONSERVATIVE, n=22) ===\n\n")

# --- Directional concordance (conservative) ---
temporal_concordance_cons_id <- counts_cons_id %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    ID1_change_cons_id = ID1 - dplyr::lag(ID1),
    ID2_change_cons_id = ID2 - dplyr::lag(ID2),
    concordant_cons_id = sign(ID1_change_cons_id) == sign(ID2_change_cons_id)
  ) %>%
  dplyr::filter(!is.na(concordant_cons_id))

n_changes_cons_id      <- nrow(temporal_concordance_cons_id)
n_concordant_cons_id   <- sum(temporal_concordance_cons_id$concordant_cons_id, na.rm = TRUE)
concordance_pct_cons_id <- (n_concordant_cons_id / n_changes_cons_id) * 100

cat(sprintf("Temporal Concordance (conservative): %.1f%%\n", concordance_pct_cons_id))
cat(sprintf("  Days analyzed: %d (n-1 lag)\n", n_changes_cons_id))
cat(sprintf("  Concordant days: %d | Discordant: %d\n\n",
            n_concordant_cons_id, n_changes_cons_id - n_concordant_cons_id))

if ((n_changes_cons_id - n_concordant_cons_id) > 0) {
  cat("Days with opposite directional changes (conservative):\n")
  discordant_cons_id <- temporal_concordance_cons_id %>%
    dplyr::filter(!concordant_cons_id) %>%
    dplyr::select(Date, ID1, ID2, ID1_change_cons_id, ID2_change_cons_id) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        ID1_change_cons_id > 0 & ID2_change_cons_id < 0 ~ "ID1 up, ID2 down",
        ID1_change_cons_id < 0 & ID2_change_cons_id > 0 ~ "ID1 down, ID2 up",
        TRUE ~ "Other"
      )
    )
  print(discordant_cons_id)
  cat("\n")
}

kendall_cons_id <- cor.test(temporal_concordance_cons_id$ID1_change_cons_id,
                            temporal_concordance_cons_id$ID2_change_cons_id,
                            method = "kendall",
                            exact  = FALSE)

cat("Kendall's tau (conservative - correlation of directional changes):\n")
cat(sprintf("  tau = %.3f | p = %.4f\n\n",
            kendall_cons_id$estimate, kendall_cons_id$p.value))

# --- Temporal plot (conservative, log-transformed) ---
id_temporal_cons_id <- counts_cons_id %>%
  dplyr::select(Date, ID1, ID2) %>%
  tidyr::pivot_longer(
    cols      = c(ID1, ID2),
    names_to  = "Device",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(Abundance_log = log10(Abundance + 1))

plot_temporal_cons_id <- ggplot(
  id_temporal_cons_id,
  aes(x = Date, y = Abundance_log, color = Device, shape = Device)
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = custom_colors,
    labels = c("ID1" = "ID1 (Site2 - Maize)",
               "ID2" = "ID2 (Site2 - Maize)")
  ) +
  scale_shape_manual(
    values = c("ID1" = 16, "ID2" = 17),
    labels = c("ID1" = "ID1 (Site2 - Maize)",
               "ID2" = "ID2 (Site2 - Maize)")
  ) +
  labs(
    title    = "Temporal Concordance: ID1 vs ID2 (Conservative, n=22)",
    subtitle = sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                       summary(model_log_cons_id)$r.squared,
                       concordance_pct_cons_id,
                       kendall_cons_id$estimate,
                       kendall_cons_id$p.value),
    x        = "Date",
    y        = "log10(Daily Abundance + 1)",
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

print(plot_temporal_cons_id)

ggsave("outputs/H1_ID_Temporal_Concordance_Conservative.png",
       plot   = plot_temporal_cons_id,
       width  = 7, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ID_Temporal_Concordance_Conservative.pdf",
       plot   = plot_temporal_cons_id,
       width  = 7, height = 6, device = cairo_pdf)

cat("Conservative temporal concordance plot saved:\n")
cat("   - outputs/H1_ID_Temporal_Concordance_Conservative.png\n")
cat("   - outputs/H1_ID_Temporal_Concordance_Conservative.pdf\n\n")

cat("=== SECTION 6B COMPLETE ===\n")

# ------------------------------------------------------------------------------
# COMBINED FIGURES (A: scatter | B: temporal concordance)
# Two versions: operational (n=17) and conservative (n=22)
# ------------------------------------------------------------------------------

cat("--- Generating combined figures (A: scatter | B: temporal) ---\n\n")

# --- Combined: Operational (n=17) ---
plot_scatter_oper_id_A <- plot_scatter_oper_id +
  labs(title = "A) ID Consistency — Operational (n=17)")

plot_temporal_id_B <- plot_temporal_id +
  labs(title = "B) Temporal Concordance: ID1 vs ID2 (Operational, n=17)")

combined_id_oper <- plot_scatter_oper_id_A + plot_temporal_id_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_id_oper)

ggsave("outputs/H1_ID_Combined_Operational.png",
       plot   = combined_id_oper,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ID_Combined_Operational.pdf",
       plot   = combined_id_oper,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved (Operational):\n")
cat("   - outputs/H1_ID_Combined_Operational.png\n")
cat("   - outputs/H1_ID_Combined_Operational.pdf\n\n")

# --- Combined: Conservative (n=22) ---
plot_scatter_cons_id_A <- plot_scatter_cons_id +
  labs(title = "A) ID Consistency — Conservative (n=22)")

plot_temporal_cons_id_B <- plot_temporal_cons_id +
  labs(title = "B) Temporal Concordance: ID1 vs ID2 (Conservative, n=22)")

combined_id_cons <- plot_scatter_cons_id_A + plot_temporal_cons_id_B +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(combined_id_cons)

ggsave("outputs/H1_ID_Combined_Conservative.png",
       plot   = combined_id_cons,
       width  = 14, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H1_ID_Combined_Conservative.pdf",
       plot   = combined_id_cons,
       width  = 14, height = 6, device = cairo_pdf)

cat("Combined figure saved (Conservative):\n")
cat("   - outputs/H1_ID_Combined_Conservative.png\n")
cat("   - outputs/H1_ID_Combined_Conservative.pdf\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("H1 ID ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n")
cat("  [OK] outputs/Table_H1_ID_AllModels_Results.csv (4-model comparison)\n")
cat("  [OK] outputs/H1_ID_Diagnostic_Comparison.png / .pdf (Supplementary)\n")
cat("  [OK] outputs/H1_ID_Consistency_Conservative.png / .pdf\n")
cat("  [OK] outputs/H1_ID_Consistency_Operational.png / .pdf\n")
cat("  [OK] outputs/H1_ID_Temporal_Concordance.png / .pdf\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("  Conservative (n=22): R2 = %.3f%s | Pseudo-R2 robust = %.3f\n",
            summary(model_log_cons_id)$r.squared,
            get_sig_stars_id(summary(model_log_cons_id)$coefficients[2, 4]),
            pseudo_r2_log_cons_id))
cat(sprintf("  Operational  (n=17): R2 = %.3f%s | Pseudo-R2 robust = %.3f\n",
            summary(model_log_oper_id)$r.squared,
            get_sig_stars_id(summary(model_log_oper_id)$coefficients[2, 4]),
            pseudo_r2_log_oper_id))
cat(sprintf("  Temporal concordance (operational): %.1f%% (%d/%d days)\n\n",
            concordance_pct_id, n_concordant_id, n_changes_id))

cat("NOTES:\n")
cat("  - ID3 excluded from all analyses (device failure)\n")
cat("  - Conservative: n=22, all days including battery-failure period\n")
cat("  - Operational: n=17, battery-failure days excluded (Aug 28 - Sep 1)\n")
cat("  - All variables carry _id suffix for cross-script traceability\n\n")

cat("========================================\n")
cat("H1 ID COMPLETE\n")
cat("========================================\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file_id, "\n")