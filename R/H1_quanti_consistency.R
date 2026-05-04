# ==============================================================================
# H1_quanti_consistency.R  — §3.2.1 Quantitative Intra-System Consistency
# FAIRD: §3.2.1.1 (Figs 4–5) + §3.2.1.2 (Table 9 robust suppl.)
# AMMOD: §3.2.1.1 (Fig 6) | Inter-Habitat: §3.2.1.3 (Figs 9–11, Table 10)
# Diagnostic: §2.6 (Fig 3)
# Sections 4–7 (robust summary, inter-habitat, tables) in Part 2.
# ==============================================================================

library(tidyverse)
library(here)
library(car)
library(lmtest)
library(patchwork)
library(scales)
library(MASS)    # load last — masks dplyr::select if not namespaced

here::i_am("R/H1_quanti_consistency.R")

dir.create("outputs/figures",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",   showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console",  showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/H1_quanti_consistency_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("H1 QUANTITATIVE CONSISTENCY — §3.2.1\n")
cat("log10(x+1) | OLS + diagnostics | Spearman & Kendall on daily changes\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source("R/00_load_data.R")
source("R/00_custom_colors.R")


# ==============================================================================
# SECTION 0: SHARED SETUP — DATA, HELPERS, CORE ANALYSIS FUNCTION
# ==============================================================================

etraps_daily  <- data_list$etraps_daily
ammod_biomass <- data_list$ammod_biomass_for_stats  # AMMOD3 failure days pre-excluded

# --- Significance stars ---
sig_stars <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

# --- Core analysis function ---
# Accepts a pre-pivoted wide-format dataframe (Date + device1 + device2 columns).
# Handles both count data (values_fill=0) and biomass data (values_fill=NA +
# pre-filtered). Returns a named list with all regression, diagnostic, and
# non-parametric metrics.
analyze_pair <- function(df_wide, device1, device2) {
  df <- df_wide %>%
    dplyr::mutate(
      dev1_log = log10(.data[[device1]] + 1),
      dev2_log = log10(.data[[device2]] + 1)
    )

  if (nrow(df) < 3) {
    warning(sprintf("Insufficient data for %s vs %s (%d rows)", device1, device2, nrow(df)))
    return(NULL)
  }

  # OLS on log scale
  model_log <- stats::lm(dev2_log ~ dev1_log, data = df)
  model_sum <- summary(model_log)

  # Diagnostics (SW, BP, DW)
  sw_test <- stats::shapiro.test(stats::residuals(model_log))
  bp_test <- car::ncvTest(model_log)
  dw_test <- lmtest::dwtest(model_log)

  # Slope, SE, CI, intercept
  coefs       <- stats::coef(model_sum)
  ci          <- stats::confint(model_log)
  slope       <- coefs["dev1_log", "Estimate"]
  slope_se    <- coefs["dev1_log", "Std. Error"]
  slope_p     <- coefs["dev1_log", "Pr(>|t|)"]
  slope_ci    <- ci["dev1_log", ]
  intercept   <- coefs["(Intercept)", "Estimate"]
  intercept_p <- coefs["(Intercept)", "Pr(>|t|)"]

  # Spearman — two variants:
  # spearman_abs: on absolute log-transformed values → Tables 8, 9, 10 (matches manuscript)
  # spearman_chg: on daily raw changes               → printed alongside Kendall for context
  sp_abs <- stats::cor.test(df$dev1_log, df$dev2_log, method = "spearman", exact = FALSE)

  # Daily changes for non-parametric stats (raw un-logged values)
  changes <- df %>%
    dplyr::arrange(Date) %>%
    dplyr::mutate(
      change1    = .data[[device1]] - dplyr::lag(.data[[device1]]),
      change2    = .data[[device2]] - dplyr::lag(.data[[device2]]),
      concordant = sign(change1) == sign(change2)
    ) %>%
    dplyr::filter(!is.na(concordant))

  n_changes    <- nrow(changes)
  n_concordant <- sum(changes$concordant)

  spearman <- stats::cor.test(changes$change1, changes$change2,
                               method = "spearman", exact = FALSE)
  kendall  <- stats::cor.test(changes$change1, changes$change2,
                               method = "kendall",  exact = FALSE)

  # Robust regression (Huber M-estimation, log scale)
  model_robust <- MASS::rlm(dev2_log ~ dev1_log, data = df, psi = MASS::psi.huber)
  ss_res       <- sum(model_robust$residuals^2)
  ss_tot       <- sum((df$dev2_log - mean(df$dev2_log))^2)
  pseudo_r2    <- 1 - ss_res / ss_tot
  n_low_wt     <- sum(model_robust$w < 0.5)

  list(
    devices         = c(device1, device2),
    n_days          = nrow(df),
    data            = df,
    changes         = changes,
    model           = model_log,
    model_robust    = model_robust,
    r_squared       = model_sum$r.squared,
    adj_r_squared   = model_sum$adj.r.squared,
    slope_p         = slope_p,
    slope           = slope,
    slope_se        = slope_se,
    slope_ci_low    = slope_ci[1],
    slope_ci_high   = slope_ci[2],
    intercept       = intercept,
    intercept_p     = intercept_p,
    shapiro_p       = sw_test$p.value,
    bp_p            = bp_test$p,
    dw_stat         = as.numeric(dw_test$statistic),
    dw_p            = dw_test$p.value,
    spearman_abs    = as.numeric(sp_abs$estimate),   # absolute log values → tables
    spearman_abs_p  = sp_abs$p.value,
    spearman_rho    = as.numeric(spearman$estimate), # daily changes → console
    spearman_p      = spearman$p.value,
    kendall_tau     = as.numeric(kendall$estimate),
    kendall_p       = kendall$p.value,
    concordance_pct = n_concordant / n_changes * 100,
    n_concordant    = n_concordant,
    n_changes       = n_changes,
    pseudo_r2       = pseudo_r2,
    n_low_wt        = n_low_wt
  )
}

# --- Console printer ---
print_pair <- function(res, label) {
  cat(sprintf("\n--- %s (n = %d) ---\n", label, res$n_days))
  cat(sprintf("  R² = %.3f | Adj.R² = %.3f | p = %.4f%s\n",
              res$r_squared, res$adj_r_squared, res$slope_p,
              sig_stars(res$slope_p)))
  cat(sprintf("  Slope = %.3f ± %.3f SE [95%% CI: %.3f, %.3f]\n",
              res$slope, res$slope_se, res$slope_ci_low, res$slope_ci_high))
  cat(sprintf("  Intercept = %.3f (p = %.3f)\n",
              res$intercept, res$intercept_p))
  cat(sprintf("  Shapiro-Wilk p = %.4f | Breusch-Pagan p = %.4f | DW stat = %.3f (p = %.3f)\n",
              res$shapiro_p, res$bp_p, res$dw_stat, res$dw_p))
  cat(sprintf("  Spearman ρ [abs log] = %.3f%s | [daily chg] = %.3f%s\n",
              res$spearman_abs, sig_stars(res$spearman_abs_p),
              res$spearman_rho, sig_stars(res$spearman_p)))
  cat(sprintf("  Kendall τ [raw chg] = %.3f%s\n",
              res$kendall_tau, sig_stars(res$kendall_p)))
  cat(sprintf("  Concordance: %.1f%% (%d/%d days)\n",
              res$concordance_pct, res$n_concordant, res$n_changes))
  cat(sprintf("  Robust pseudo-R² = %.3f | Low-weight days (<0.5): %d/%d\n",
              res$pseudo_r2, res$n_low_wt, res$n_days))
}

# --- Figure helpers ---

# Scatter panel: res$data must contain dev1_log, dev2_log
make_scatter <- function(res, panel_title, x_lab, y_lab, fill_color,
                          shared_limits = NULL) {
  sig <- sig_stars(res$slope_p)
  lab <- sprintf(
    "R² = %.3f%s | Slope = %.3f [%.3f, %.3f]\nSW p = %.3f | BP p = %.3f | DW = %.2f",
    res$r_squared, sig, res$slope,
    res$slope_ci_low, res$slope_ci_high,
    res$shapiro_p, res$bp_p, res$dw_stat
  )
  p <- ggplot2::ggplot(res$data, ggplot2::aes(x = dev1_log, y = dev2_log)) +
    ggplot2::geom_point(size = 3, alpha = 0.7, color = "gray30") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "gray30",
                          fill = fill_color, alpha = 0.3, linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "gray50", linewidth = 0.8) +
    ggplot2::annotate("text", x = Inf, y = -Inf, label = lab,
                      hjust = 1.05, vjust = -0.3, size = 3,
                      color = "gray20", fontface = "italic") +
    ggplot2::labs(title = panel_title, x = x_lab, y = y_lab) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      panel.border     = ggplot2::element_rect(color = "gray80", fill = NA,
                                                linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      aspect.ratio     = 1
    )
  if (!is.null(shared_limits))
    p <- p + ggplot2::coord_cartesian(xlim = shared_limits, ylim = shared_limits)
  p
}

# Temporal concordance panel: df_long must contain Date, Device, value_log
make_temporal <- function(df_long, panel_title, subtitle_text, y_lab,
                           color_vals, shape_vals, label_vals = NULL) {
  scale_labs <- if (!is.null(label_vals)) label_vals else ggplot2::waiver()
  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = Date, y = value_log,
                               color = Device, shape = Device)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.7) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = color_vals, labels = scale_labs) +
    ggplot2::scale_shape_manual(values = shape_vals, labels = scale_labs) +
    ggplot2::labs(title    = panel_title,
                  subtitle = subtitle_text,
                  x = "Date", y = y_lab,
                  color = "Device", shape = "Device") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "gray40"),
      panel.border     = ggplot2::element_rect(color = "gray80", fill = NA,
                                                linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top",
      legend.title     = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

# Build temporal subtitle from a pair result
# Spearman ρ shown is on absolute log values (matching Tables 8/9/10).
temporal_subtitle <- function(res) {
  sprintf(
    "R² = %.3f | Concordance: %.1f%% | Spearman ρ = %.3f%s | Kendall τ = %.3f%s",
    res$r_squared,
    res$concordance_pct,
    res$spearman_abs, sig_stars(res$spearman_abs_p),
    res$kendall_tau,  sig_stars(res$kendall_p)
  )
}

# Pivot pair result to long format for temporal plotting
pair_to_long <- function(df_wide, device1, device2) {
  df_wide %>%
    dplyr::select(Date, dplyr::all_of(c(device1, device2))) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(c(device1, device2)),
      names_to  = "Device",
      values_to = "value"
    ) %>%
    dplyr::mutate(value_log = log10(value + 1))
}


# ==============================================================================
# SECTION 1: DIAGNOSTIC — HETEROSCEDASTICITY (§2.6, Figure 3)
# Representative pair: FAIRD1 vs FAIRD2 (Maize) — justifies log10(x+1) decision
# ==============================================================================

cat("=== SECTION 1: DIAGNOSTIC HETEROSCEDASTICITY (Fig 3) ===\n\n")

# --- Data prep ---
diag_wide <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

cat(sprintf("Diagnostic pair: FAIRD1 vs FAIRD2 (Maize) | n = %d days\n\n",
            nrow(diag_wide)))

# --- Original scale model ---
model_orig    <- stats::lm(FAIRD2 ~ FAIRD1, data = diag_wide)
sw_orig       <- stats::shapiro.test(stats::residuals(model_orig))
bp_orig       <- car::ncvTest(model_orig)
dw_orig       <- lmtest::dwtest(model_orig)
r2_orig       <- summary(model_orig)$r.squared
slope_orig    <- stats::coef(model_orig)[2]
p_orig        <- summary(model_orig)$coefficients[2, 4]
threshold_cooks <- 4 / nrow(diag_wide)
cooks_d       <- stats::cooks.distance(model_orig)
n_influential <- sum(cooks_d > threshold_cooks)

# --- Log-transformed model ---
diag_log <- diag_wide %>%
  dplyr::mutate(FAIRD1_log = log10(FAIRD1 + 1),
                FAIRD2_log = log10(FAIRD2 + 1))

model_log_diag <- stats::lm(FAIRD2_log ~ FAIRD1_log, data = diag_log)
sw_log         <- stats::shapiro.test(stats::residuals(model_log_diag))
bp_log         <- car::ncvTest(model_log_diag)
dw_log         <- lmtest::dwtest(model_log_diag)
r2_log_diag    <- summary(model_log_diag)$r.squared
slope_log_diag <- stats::coef(model_log_diag)[2]
p_log_diag     <- summary(model_log_diag)$coefficients[2, 4]

cat("Original scale:\n")
cat(sprintf("  R² = %.3f | SW p = %.4f (%s) | BP p = %.4f (%s) | DW = %.3f\n",
            r2_orig,
            sw_orig$p.value, ifelse(sw_orig$p.value > 0.05, "normal", "NOT normal"),
            bp_orig$p,      ifelse(bp_orig$p > 0.05, "homoscedastic", "HETEROSCEDASTIC"),
            as.numeric(dw_orig$statistic)))
cat(sprintf("  Cook's D (threshold %.4f): %d/%d influential points\n\n",
            threshold_cooks, n_influential, nrow(diag_wide)))

cat("Log-transformed:\n")
cat(sprintf("  R² = %.3f | SW p = %.4f (%s) | BP p = %.4f (%s) | DW = %.3f\n",
            r2_log_diag,
            sw_log$p.value, ifelse(sw_log$p.value > 0.05, "normal", "NOT normal"),
            bp_log$p,       ifelse(bp_log$p > 0.05, "homoscedastic", "HETEROSCEDASTIC"),
            as.numeric(dw_log$statistic)))
cat(sprintf("\n[DECISION] Use log10(x+1): homoscedasticity restored (BP p = %.4f)\n\n",
            bp_log$p))

# --- Plot data with influential flags from original model ---
plot_diag <- diag_log %>%
  dplyr::mutate(
    cooks_d     = stats::cooks.distance(model_orig),
    influential = cooks_d > threshold_cooks
  )

label_orig <- sprintf(
  "R² = %.3f%s  Slope = %.2f\nSW p = %.3f | BP p = %.3f | DW = %.2f",
  r2_orig,       sig_stars(p_orig),       slope_orig,
  sw_orig$p.value, bp_orig$p, as.numeric(dw_orig$statistic)
)
label_log_d <- sprintf(
  "R² = %.3f%s  Slope = %.2f\nSW p = %.3f | BP p = %.3f | DW = %.2f",
  r2_log_diag,   sig_stars(p_log_diag),   slope_log_diag,
  sw_log$p.value,  bp_log$p,  as.numeric(dw_log$statistic)
)

# Panel A: original scale
p_diag_A <- ggplot2::ggplot(plot_diag,
                             ggplot2::aes(x = FAIRD1, y = FAIRD2)) +
  ggplot2::geom_point(ggplot2::aes(color = influential, size = influential),
                      alpha = 0.7) +
  ggplot2::scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray30"),
    labels = c("TRUE" = "Influential (Cook's D)", "FALSE" = "Normal")
  ) +
  ggplot2::scale_size_manual(
    values = c("TRUE" = 4, "FALSE" = 3),
    labels = c("TRUE" = "Influential (Cook's D)", "FALSE" = "Normal")
  ) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, color = "gray30",
                        fill = custom_colors["FAIRD"], alpha = 0.3, linewidth = 1) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                        color = "gray50", linewidth = 0.8) +
  ggplot2::annotate("text", x = Inf, y = -Inf, label = label_orig,
                    hjust = 1.05, vjust = -0.3, size = 3,
                    color = "gray20", fontface = "italic") +
  ggplot2::labs(title = "A) Original Scale",
                x = "FAIRD1 Daily Abundance",
                y = "FAIRD2 Daily Abundance",
                color = "Point type", size = "Point type") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border     = ggplot2::element_rect(color = "gray80", fill = NA,
                                              linewidth = 0.5),
    legend.position  = "bottom",
    aspect.ratio     = 1
  )

# Panel B: log scale — same influential classification as Panel A
p_diag_B <- ggplot2::ggplot(plot_diag,
                             ggplot2::aes(x = FAIRD1_log, y = FAIRD2_log)) +
  ggplot2::geom_point(ggplot2::aes(color = influential, size = influential),
                      alpha = 0.7) +
  ggplot2::scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray30"),
    labels = c("TRUE" = "Influential (Cook's D)", "FALSE" = "Normal")
  ) +
  ggplot2::scale_size_manual(
    values = c("TRUE" = 4, "FALSE" = 3),
    labels = c("TRUE" = "Influential (Cook's D)", "FALSE" = "Normal")
  ) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, color = "gray30",
                        fill = custom_colors["FAIRD"], alpha = 0.3, linewidth = 1) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                        color = "gray50", linewidth = 0.8) +
  ggplot2::annotate("text", x = Inf, y = -Inf, label = label_log_d,
                    hjust = 1.05, vjust = -0.3, size = 3,
                    color = "gray20", fontface = "italic") +
  ggplot2::labs(title = "B) log₁₀-Transformed Scale",
                x = "log₁₀(FAIRD1 + 1)",
                y = "log₁₀(FAIRD2 + 1)",
                color = "Point type", size = "Point type") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border     = ggplot2::element_rect(color = "gray80", fill = NA,
                                              linewidth = 0.5),
    legend.position  = "bottom",
    aspect.ratio     = 1
  )

fig3 <- p_diag_A + p_diag_B +
  patchwork::plot_layout(guides = "collect", ncol = 2) &
  ggplot2::theme(legend.position = "bottom")

print(fig3)
ggplot2::ggsave("outputs/figures/Fig03_H1_diagnostic_heteroscedasticity.png",
                plot = fig3, width = 12, height = 6, dpi = 300)
cat("[OK] Fig03_H1_diagnostic_heteroscedasticity.png\n\n")


# ==============================================================================
# SECTION 2: INTRA-HABITAT LEVEL — FAIRD (§3.2.1.1, Figs 4–5, Table 9 FAIRD rows)
# Pair 1: FAIRD1 vs FAIRD2 (Maize, n = 22)
# Pair 2: FAIRD3 vs FAIRD4 (Meadow, n = 22)
# ==============================================================================

cat("=== SECTION 2: FAIRD INTRA-HABITAT ===\n\n")

# --- 2a. Data prep ---
faird_maize_wide <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

faird_meadow_wide <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD3", "FAIRD4")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

# --- 2b. Analyses ---
res_faird_maize  <- analyze_pair(faird_maize_wide,  "FAIRD1", "FAIRD2")
res_faird_meadow <- analyze_pair(faird_meadow_wide, "FAIRD3", "FAIRD4")

print_pair(res_faird_maize,  "FAIRD1 vs FAIRD2 (Maize)")
print_pair(res_faird_meadow, "FAIRD3 vs FAIRD4 (Meadow)")
cat("\n")

# --- 2c. Shared axis limits (both FAIRD pairs on same scale) ---
all_faird_log <- c(
  res_faird_maize$data$dev1_log,  res_faird_maize$data$dev2_log,
  res_faird_meadow$data$dev1_log, res_faird_meadow$data$dev2_log
)
buf_f        <- diff(range(all_faird_log)) * 0.05
shared_lim_f <- c(min(all_faird_log) - buf_f, max(all_faird_log) + buf_f)
cat(sprintf("Shared FAIRD axis limits: [%.3f, %.3f]\n\n", shared_lim_f[1], shared_lim_f[2]))

# --- 2d. Fig 4: FAIRD Maize (FAIRD1 vs FAIRD2) ---
p4_scatter <- make_scatter(
  res         = res_faird_maize,
  panel_title = "A) FAIRD Maize — Site 1 vs Site 2",
  x_lab       = "log₁₀(FAIRD1 + 1)",
  y_lab       = "log₁₀(FAIRD2 + 1)",
  fill_color  = custom_colors["FAIRD"],
  shared_limits = shared_lim_f
)

faird_maize_long <- pair_to_long(faird_maize_wide, "FAIRD1", "FAIRD2")

p4_temporal <- make_temporal(
  df_long      = faird_maize_long,
  panel_title  = "B) Temporal Concordance: FAIRD1 vs FAIRD2 (Maize)",
  subtitle_text = temporal_subtitle(res_faird_maize),
  y_lab        = "log₁₀(Abundance + 1)",
  color_vals   = custom_colors,
  shape_vals   = c("FAIRD1" = 16L, "FAIRD2" = 17L),
  label_vals   = c("FAIRD1" = "FAIRD1 (Site 1)", "FAIRD2" = "FAIRD2 (Site 2)")
)

fig4 <- p4_scatter + p4_temporal +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(fig4)
ggplot2::ggsave("outputs/figures/Fig04_H1_FAIRD_Maize.png",
                plot = fig4, width = 14, height = 6, dpi = 300)
cat("[OK] Fig04_H1_FAIRD_Maize.png\n\n")

# --- 2e. Fig 5: FAIRD Meadow (FAIRD3 vs FAIRD4) ---
p5_scatter <- make_scatter(
  res         = res_faird_meadow,
  panel_title = "A) FAIRD Meadow — Site 3 vs Site 4",
  x_lab       = "log₁₀(FAIRD3 + 1)",
  y_lab       = "log₁₀(FAIRD4 + 1)",
  fill_color  = custom_colors["FAIRD"],
  shared_limits = shared_lim_f
)

faird_meadow_long <- pair_to_long(faird_meadow_wide, "FAIRD3", "FAIRD4")

p5_temporal <- make_temporal(
  df_long      = faird_meadow_long,
  panel_title  = "B) Temporal Concordance: FAIRD3 vs FAIRD4 (Meadow)",
  subtitle_text = temporal_subtitle(res_faird_meadow),
  y_lab        = "log₁₀(Abundance + 1)",
  color_vals   = custom_colors,
  shape_vals   = c("FAIRD3" = 16L, "FAIRD4" = 17L),
  label_vals   = c("FAIRD3" = "FAIRD3 (Site 3)", "FAIRD4" = "FAIRD4 (Site 4)")
)

fig5 <- p5_scatter + p5_temporal +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(fig5)
ggplot2::ggsave("outputs/figures/Fig05_H1_FAIRD_Meadow.png",
                plot = fig5, width = 14, height = 6, dpi = 300)
cat("[OK] Fig05_H1_FAIRD_Meadow.png\n\n")

# --- 2f. FAIRD intra-habitat summary to console ---
cat("--- FAIRD Intra-Habitat Summary ---\n")
cat(sprintf("  Maize  (F1 vs F2): R² = %.3f%s | Concordance: %.1f%% | Spearman ρ = %.3f%s | Kendall τ = %.3f%s\n",
            res_faird_maize$r_squared,  sig_stars(res_faird_maize$slope_p),
            res_faird_maize$concordance_pct,
            res_faird_maize$spearman_rho,  sig_stars(res_faird_maize$spearman_p),
            res_faird_maize$kendall_tau,   sig_stars(res_faird_maize$kendall_p)))
cat(sprintf("  Meadow (F3 vs F4): R² = %.3f%s | Concordance: %.1f%% | Spearman ρ = %.3f%s | Kendall τ = %.3f%s\n\n",
            res_faird_meadow$r_squared, sig_stars(res_faird_meadow$slope_p),
            res_faird_meadow$concordance_pct,
            res_faird_meadow$spearman_rho, sig_stars(res_faird_meadow$spearman_p),
            res_faird_meadow$kendall_tau,  sig_stars(res_faird_meadow$kendall_p)))


# ==============================================================================
# SECTION 3: INTRA-HABITAT LEVEL — AMMOD (§3.2.1.1, Fig 6, Table 9 AMMOD row)
# Pair: AMMOD1 vs AMMOD2 (Maize only, n = 12)
# No Meadow pair — AMMOD3 mechanical failure, AMMOD4 complete data loss
# Data source: ammod_biomass_for_stats (AMMOD3 failure days pre-excluded)
# Metric: Live_mass (biomass in mg) | values_fill = NA — genuine missing data
# ==============================================================================

cat("=== SECTION 3: AMMOD INTRA-HABITAT ===\n\n")
cat("Note: Only Maize pair available (AMMOD1 vs AMMOD2, Sites 1–2).\n")
cat("      AMMOD3 excluded (mechanical failure Aug 30–Sep 3).\n")
cat("      AMMOD4 excluded (complete data loss).\n\n")

# --- 3a. Data prep (NA fill — biomass NAs are real missing data, not zeros) ---
ammod_maize_wide <- ammod_biomass %>%
  dplyr::filter(Device %in% c("AMMOD1", "AMMOD2")) %>%
  dplyr::select(Date, Device, Live_mass) %>%
  tidyr::pivot_wider(names_from  = Device,
                     values_from = Live_mass,
                     values_fill = NA) %>%
  dplyr::filter(!is.na(AMMOD1) & !is.na(AMMOD2))

cat(sprintf("AMMOD Maize pair: n = %d days with complete biomass data\n",
            nrow(ammod_maize_wide)))
cat(sprintf("  Date range: %s to %s\n\n",
            min(ammod_maize_wide$Date), max(ammod_maize_wide$Date)))

# --- 3b. Analysis ---
res_ammod_maize <- analyze_pair(ammod_maize_wide, "AMMOD1", "AMMOD2")
print_pair(res_ammod_maize, "AMMOD1 vs AMMOD2 (Maize)")
cat("\n")

# --- 3c. Axis limits (AMMOD biomass scale, standalone) ---
ammod_log_vals  <- c(res_ammod_maize$data$dev1_log, res_ammod_maize$data$dev2_log)
buf_a           <- diff(range(ammod_log_vals)) * 0.08
shared_lim_a    <- c(min(ammod_log_vals) - buf_a, max(ammod_log_vals) + buf_a)

# --- 3d. Fig 6: AMMOD Maize (AMMOD1 vs AMMOD2) ---
p6_scatter <- make_scatter(
  res         = res_ammod_maize,
  panel_title = "A) AMMOD Maize — Site 1 vs Site 2",
  x_lab       = "log₁₀(AMMOD1 Biomass + 1)",
  y_lab       = "log₁₀(AMMOD2 Biomass + 1)",
  fill_color  = custom_colors["AMMOD"],
  shared_limits = shared_lim_a
)

ammod_maize_long <- pair_to_long(ammod_maize_wide, "AMMOD1", "AMMOD2")

p6_temporal <- make_temporal(
  df_long      = ammod_maize_long,
  panel_title  = "B) Temporal Concordance: AMMOD1 vs AMMOD2 (Maize)",
  subtitle_text = temporal_subtitle(res_ammod_maize),
  y_lab        = "log₁₀(Live Biomass + 1)",
  color_vals   = custom_colors,
  shape_vals   = c("AMMOD1" = 16L, "AMMOD2" = 17L),
  label_vals   = c("AMMOD1" = "AMMOD1 (Site 1)", "AMMOD2" = "AMMOD2 (Site 2)")
)

fig6 <- p6_scatter + p6_temporal +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(fig6)
ggplot2::ggsave("outputs/figures/Fig06_H1_AMMOD_Maize.png",
                plot = fig6, width = 14, height = 6, dpi = 300)
cat("[OK] Fig06_H1_AMMOD_Maize.png\n\n")

# --- 3e. AMMOD intra-habitat summary to console ---
cat("--- AMMOD Intra-Habitat Summary ---\n")
cat(sprintf("  Maize (A1 vs A2): R² = %.3f%s | Concordance: %.1f%% | Spearman ρ = %.3f%s | Kendall τ = %.3f%s\n",
            res_ammod_maize$r_squared, sig_stars(res_ammod_maize$slope_p),
            res_ammod_maize$concordance_pct,
            res_ammod_maize$spearman_rho, sig_stars(res_ammod_maize$spearman_p),
            res_ammod_maize$kendall_tau,  sig_stars(res_ammod_maize$kendall_p)))
cat("  Meadow: no pair available (AMMOD4 failed)\n\n")

cat("==============================================================================\n")
cat("PART 1 COMPLETE (Sections 1–3) — continue with Sections 4–7\n")
cat("==============================================================================\n")
cat("Script checkpoint:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 4: ROBUST REGRESSION SUPPLEMENTARY (Table 9 supplement)
# OLS vs robust comparison for all intra-habitat pairs already computed above.
# Quantifies outlier influence: if Pseudo-R² ≈ OLS R², the OLS result is robust.
# ==============================================================================

cat("=== SECTION 4: ROBUST REGRESSION SUPPLEMENTARY ===\n\n")

make_robust_row <- function(res, pair_label, habitat, system) {
  tibble::tibble(
    System      = system,
    Pair        = pair_label,
    Habitat     = habitat,
    n           = res$n_days,
    OLS_R2      = res$r_squared,
    OLS_AdjR2   = res$adj_r_squared,
    Pseudo_R2   = res$pseudo_r2,
    n_low_wt    = res$n_low_wt,
    Delta_R2    = res$pseudo_r2 - res$r_squared,
    SW_p        = res$shapiro_p,
    BP_p        = res$bp_p,
    DW_stat     = res$dw_stat,
    DW_p        = res$dw_p
  )
}

table9_robust <- dplyr::bind_rows(
  make_robust_row(res_faird_maize,  "FAIRD1 vs FAIRD2", "Maize",  "FAIRD"),
  make_robust_row(res_faird_meadow, "FAIRD3 vs FAIRD4", "Meadow", "FAIRD"),
  make_robust_row(res_ammod_maize,  "AMMOD1 vs AMMOD2", "Maize",  "AMMOD")
)

cat("OLS R² vs Pseudo-R² (robust, log scale | low-weight threshold = 0.5):\n\n")
print(table9_robust, n = Inf)

cat("\nInterpretation:\n")
cat("  |ΔR²| < 0.05 → OLS result not driven by outliers\n")
cat("  Low-weight days (<0.5 Huber weight) indicate high-leverage observations\n\n")

readr::write_csv(
  table9_robust %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))),
  "outputs/tables/Table09_supp_OLS_vs_robust.csv"
)
cat("[OK] Table09_supp_OLS_vs_robust.csv\n\n")


# ==============================================================================
# SECTION 5: INTER-HABITAT LEVEL — FAIRD (§3.2.1.3, Fig 9, Table 10 FAIRD row)
# Comparison: mean(FAIRD1+FAIRD2) [Maize] vs mean(FAIRD3+FAIRD4) [Meadow]
# n = 22 days | Metric: mean daily abundance (log10(x+1))
# Note: Averaging approach — pooling invalid (site-pair covariate p = 0.013)
# ==============================================================================

cat("=== SECTION 5: INTER-HABITAT FAIRD ===\n\n")
cat("Note: Averaging approach — site-pair covariate p = 0.013, pooling invalid.\n")
cat("      Compare daily mean(FAIRD1+FAIRD2) [Maize] vs mean(FAIRD3+FAIRD4) [Meadow].\n\n")

# --- 5a. Data prep: daily mean abundance per habitat ---
faird_inter_wide <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1","FAIRD2","FAIRD3","FAIRD4")) %>%
  dplyr::mutate(Habitat = dplyr::case_when(
    Device %in% c("FAIRD1","FAIRD2") ~ "Maize",
    Device %in% c("FAIRD3","FAIRD4") ~ "Meadow"
  )) %>%
  dplyr::group_by(Date, Habitat) %>%
  dplyr::summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from  = Habitat,
                     values_from = mean_abundance,
                     values_fill = 0) %>%
  dplyr::arrange(Date)

cat(sprintf("FAIRD inter-habitat: n = %d days\n", nrow(faird_inter_wide)))
cat(sprintf("  Date range: %s to %s\n\n",
            min(faird_inter_wide$Date), max(faird_inter_wide$Date)))

# --- 5b. Analysis ---
res_faird_inter <- analyze_pair(faird_inter_wide, "Maize", "Meadow")
print_pair(res_faird_inter, "FAIRD Maize vs Meadow (Inter-Habitat)")
cat("\n")

# --- 5c. Axis limits ---
fi_log_vals   <- c(res_faird_inter$data$dev1_log, res_faird_inter$data$dev2_log)
buf_fi        <- diff(range(fi_log_vals)) * 0.08
shared_lim_fi <- c(min(fi_log_vals) - buf_fi, max(fi_log_vals) + buf_fi)

# --- 5d. Fig 9 ---
faird_inter_pal <- setNames(
  c(custom_colors[["FAIRD1"]], custom_colors[["FAIRD3"]]),
  c("Maize", "Meadow")
)

p9_scatter <- make_scatter(
  res           = res_faird_inter,
  panel_title   = "A) FAIRD Inter-Habitat: Maize vs Meadow",
  x_lab         = "log₁₀(Mean Maize Abundance + 1)",
  y_lab         = "log₁₀(Mean Meadow Abundance + 1)",
  fill_color    = custom_colors["FAIRD"],
  shared_limits = shared_lim_fi
)

faird_inter_long <- pair_to_long(faird_inter_wide, "Maize", "Meadow")

p9_temporal <- make_temporal(
  df_long       = faird_inter_long,
  panel_title   = "B) Temporal Concordance: FAIRD Maize vs Meadow",
  subtitle_text = temporal_subtitle(res_faird_inter),
  y_lab         = "log₁₀(Mean Abundance + 1)",
  color_vals    = faird_inter_pal,
  shape_vals    = c("Maize" = 16L, "Meadow" = 17L),
  label_vals    = c("Maize"  = "Maize (mean FAIRD1+2)",
                    "Meadow" = "Meadow (mean FAIRD3+4)")
)

fig9 <- p9_scatter + p9_temporal +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(fig9)
ggplot2::ggsave("outputs/figures/Fig09_H1_FAIRD_InterHabitat.png",
                plot = fig9, width = 14, height = 6, dpi = 300)
cat("[OK] Fig09_H1_FAIRD_InterHabitat.png\n\n")

cat("--- FAIRD Inter-Habitat Summary ---\n")
cat(sprintf("  Maize vs Meadow: R² = %.3f%s | Concordance: %.1f%% | Spearman ρ = %.3f%s | Kendall τ = %.3f%s\n\n",
            res_faird_inter$r_squared, sig_stars(res_faird_inter$slope_p),
            res_faird_inter$concordance_pct,
            res_faird_inter$spearman_rho, sig_stars(res_faird_inter$spearman_p),
            res_faird_inter$kendall_tau,  sig_stars(res_faird_inter$kendall_p)))


# ==============================================================================
# SECTION 6: INTER-HABITAT LEVEL — AMMOD (§3.2.1.3, Figs 10–11, Table 10 AMMOD rows)
# Two separate pairs — no averaging because AMMOD4 failed entirely:
#   Pair 1: AMMOD1 (Site 1, Maize)  vs AMMOD3 (Site 3, Meadow)
#   Pair 2: AMMOD2 (Site 2, Maize)  vs AMMOD3 (Site 3, Meadow)
# Data: data_list$ammod_biomass (raw, NOT _for_stats — NA filter is exclusion mechanism)
# ==============================================================================

cat("=== SECTION 6: INTER-HABITAT AMMOD ===\n\n")
cat("Note: Two separate pairs — averaging impossible (AMMOD4 failed).\n")
cat("      data_list$ammod_biomass used (raw); NA filter excludes invalid days.\n\n")

# --- 6a. Data prep: wide format with AMMOD1, AMMOD2, AMMOD3 ---
ammod_inter_wide <- data_list$ammod_biomass %>%
  dplyr::filter(Device %in% c("AMMOD1","AMMOD2","AMMOD3")) %>%
  dplyr::select(Date, Device, Live_mass) %>%
  tidyr::pivot_wider(names_from  = Device,
                     values_from = Live_mass,
                     values_fill = NA) %>%
  dplyr::arrange(Date)

# Pair 1: AMMOD1 vs AMMOD3 (Site 1 Maize vs Site 3 Meadow)
ammod_a1a3_wide <- ammod_inter_wide %>%
  dplyr::filter(!is.na(AMMOD1) & !is.na(AMMOD3))

# Pair 2: AMMOD2 vs AMMOD3 (Site 2 Maize vs Site 3 Meadow)
ammod_a2a3_wide <- ammod_inter_wide %>%
  dplyr::filter(!is.na(AMMOD2) & !is.na(AMMOD3))

cat(sprintf("AMMOD1 vs AMMOD3: n = %d overlapping biomass days\n", nrow(ammod_a1a3_wide)))
cat(sprintf("AMMOD2 vs AMMOD3: n = %d overlapping biomass days\n\n", nrow(ammod_a2a3_wide)))

# --- 6b. Analysis: Pair 1 (AMMOD1 vs AMMOD3) ---
res_ammod_a1a3 <- analyze_pair(ammod_a1a3_wide, "AMMOD1", "AMMOD3")
print_pair(res_ammod_a1a3, "AMMOD1 vs AMMOD3 (Site 1 Maize vs Site 3 Meadow)")
cat("\n")

# --- 6c. Analysis: Pair 2 (AMMOD2 vs AMMOD3) ---
res_ammod_a2a3 <- analyze_pair(ammod_a2a3_wide, "AMMOD2", "AMMOD3")
print_pair(res_ammod_a2a3, "AMMOD2 vs AMMOD3 (Site 2 Maize vs Site 3 Meadow)")
cat("\n")

# --- 6d. Axis limits (independent scales per pair — different biomass ranges) ---
a1a3_log_vals   <- c(res_ammod_a1a3$data$dev1_log, res_ammod_a1a3$data$dev2_log)
buf_a1a3        <- diff(range(a1a3_log_vals)) * 0.08
shared_lim_a1a3 <- c(min(a1a3_log_vals) - buf_a1a3, max(a1a3_log_vals) + buf_a1a3)

a2a3_log_vals   <- c(res_ammod_a2a3$data$dev1_log, res_ammod_a2a3$data$dev2_log)
buf_a2a3        <- diff(range(a2a3_log_vals)) * 0.08
shared_lim_a2a3 <- c(min(a2a3_log_vals) - buf_a2a3, max(a2a3_log_vals) + buf_a2a3)

# --- 6e. Fig 10: AMMOD1 vs AMMOD3 ---
p10_scatter <- make_scatter(
  res           = res_ammod_a1a3,
  panel_title   = "A) AMMOD Inter-Habitat: Site 1 vs Site 3",
  x_lab         = "log₁₀(AMMOD1 Biomass + 1)  [Maize]",
  y_lab         = "log₁₀(AMMOD3 Biomass + 1)  [Meadow]",
  fill_color    = custom_colors["AMMOD"],
  shared_limits = shared_lim_a1a3
)

ammod_a1a3_long <- pair_to_long(ammod_a1a3_wide, "AMMOD1", "AMMOD3")

p10_temporal <- make_temporal(
  df_long       = ammod_a1a3_long,
  panel_title   = "B) Temporal Concordance: AMMOD1 vs AMMOD3",
  subtitle_text = temporal_subtitle(res_ammod_a1a3),
  y_lab         = "log₁₀(Live Biomass + 1)",
  color_vals    = custom_colors,
  shape_vals    = c("AMMOD1" = 16L, "AMMOD3" = 15L),
  label_vals    = c("AMMOD1" = "AMMOD1 (Site 1, Maize)",
                    "AMMOD3" = "AMMOD3 (Site 3, Meadow)")
)

fig10 <- p10_scatter + p10_temporal +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(fig10)
ggplot2::ggsave("outputs/figures/Fig10_H1_AMMOD_InterHabitat_A1vsA3.png",
                plot = fig10, width = 14, height = 6, dpi = 300)
cat("[OK] Fig10_H1_AMMOD_InterHabitat_A1vsA3.png\n\n")

# --- 6f. Fig 11: AMMOD2 vs AMMOD3 ---
p11_scatter <- make_scatter(
  res           = res_ammod_a2a3,
  panel_title   = "A) AMMOD Inter-Habitat: Site 2 vs Site 3",
  x_lab         = "log₁₀(AMMOD2 Biomass + 1)  [Maize]",
  y_lab         = "log₁₀(AMMOD3 Biomass + 1)  [Meadow]",
  fill_color    = custom_colors["AMMOD"],
  shared_limits = shared_lim_a2a3
)

ammod_a2a3_long <- pair_to_long(ammod_a2a3_wide, "AMMOD2", "AMMOD3")

p11_temporal <- make_temporal(
  df_long       = ammod_a2a3_long,
  panel_title   = "B) Temporal Concordance: AMMOD2 vs AMMOD3",
  subtitle_text = temporal_subtitle(res_ammod_a2a3),
  y_lab         = "log₁₀(Live Biomass + 1)",
  color_vals    = custom_colors,
  shape_vals    = c("AMMOD2" = 17L, "AMMOD3" = 15L),
  label_vals    = c("AMMOD2" = "AMMOD2 (Site 2, Maize)",
                    "AMMOD3" = "AMMOD3 (Site 3, Meadow)")
)

fig11 <- p11_scatter + p11_temporal +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))

print(fig11)
ggplot2::ggsave("outputs/figures/Fig11_H1_AMMOD_InterHabitat_A2vsA3.png",
                plot = fig11, width = 14, height = 6, dpi = 300)
cat("[OK] Fig11_H1_AMMOD_InterHabitat_A2vsA3.png\n\n")

cat("--- AMMOD Inter-Habitat Summary ---\n")
cat(sprintf("  AMMOD1 vs AMMOD3: R² = %.3f%s | Concordance: %.1f%% | Spearman ρ [abs]=%.3f%s | τ=%.3f%s\n",
            res_ammod_a1a3$r_squared, sig_stars(res_ammod_a1a3$slope_p),
            res_ammod_a1a3$concordance_pct,
            res_ammod_a1a3$spearman_abs, sig_stars(res_ammod_a1a3$spearman_abs_p),
            res_ammod_a1a3$kendall_tau,  sig_stars(res_ammod_a1a3$kendall_p)))
cat(sprintf("  AMMOD2 vs AMMOD3: R² = %.3f%s | Concordance: %.1f%% | Spearman ρ [abs]=%.3f%s | τ=%.3f%s\n\n",
            res_ammod_a2a3$r_squared, sig_stars(res_ammod_a2a3$slope_p),
            res_ammod_a2a3$concordance_pct,
            res_ammod_a2a3$spearman_abs, sig_stars(res_ammod_a2a3$spearman_abs_p),
            res_ammod_a2a3$kendall_tau,  sig_stars(res_ammod_a2a3$kendall_p)))


# ==============================================================================
# SECTION 7: SUMMARY TABLES (Tables 8, 9, 10 — §3.2.1)
# Consolidates all H1 quantitative results into publication-ready CSV tables.
#   Table 8: H1 overview — all pairs, compact view
#   Table 9: Intra-Habitat Level (FAIRD Maize, FAIRD Meadow, AMMOD Maize)
#   Table 10: Inter-Habitat Level (FAIRD inter, AMMOD1vA3, AMMOD2vA3)
# ==============================================================================

cat("=== SECTION 7: SUMMARY TABLES ===\n\n")

# Helper: one table row from analyze_pair result
res_to_row <- function(res, pair_label, habitat, level, system) {
  tibble::tibble(
    System          = system,
    Level           = level,
    Pair            = pair_label,
    Habitat         = habitat,
    n               = res$n_days,
    R2              = res$r_squared,
    Adj_R2          = res$adj_r_squared,
    Slope           = res$slope,
    Slope_CI_low    = res$slope_ci_low,
    Slope_CI_high   = res$slope_ci_high,
    Slope_p         = res$slope_p,
    Intercept       = res$intercept,
    Intercept_p     = res$intercept_p,
    Shapiro_p       = res$shapiro_p,
    BP_p            = res$bp_p,
    DW_stat         = res$dw_stat,
    DW_p            = res$dw_p,
    Pseudo_R2       = res$pseudo_r2,
    n_low_wt        = res$n_low_wt,
    Concordance_pct = res$concordance_pct,
    Spearman_rho    = res$spearman_abs,    # absolute log values → manuscript tables
    Spearman_p      = res$spearman_abs_p,
    Kendall_tau     = res$kendall_tau,
    Kendall_p       = res$kendall_p
  )
}

# Display formatter (shared across tables)
format_table <- function(df_raw, include_habitat = TRUE) {
  out <- df_raw %>%
    dplyr::mutate(
      `R² (log)`     = paste0(formatC(R2, digits = 3, format = "f"),
                                    " ", sig_stars(Slope_p)),
      `Slope [95% CI]` = sprintf("%.2f [%.2f, %.2f]",
                                  Slope, Slope_CI_low, Slope_CI_high),
      `Concordance`    = sprintf("%.1f%%", Concordance_pct),
      `Spearman ρ` = paste0(formatC(Spearman_rho, digits = 3, format = "f"),
                                   " ", sig_stars(Spearman_p)),
      `Kendall τ`  = paste0(formatC(Kendall_tau, digits = 3, format = "f"),
                                   " ", sig_stars(Kendall_p)),
      `SW p`           = formatC(Shapiro_p, digits = 3, format = "f"),
      `BP p`           = formatC(BP_p,      digits = 3, format = "f"),
      `DW`             = sprintf("%.2f (p=%.3f)", DW_stat, DW_p),
      `Pseudo-R²` = formatC(Pseudo_R2, digits = 3, format = "f")
    )
  if (include_habitat) {
    out %>% dplyr::select(System, Pair, Habitat, n,
                           `R² (log)`, `Slope [95% CI]`,
                           `Concordance`, `Spearman ρ`, `Kendall τ`,
                           `SW p`, `BP p`, `DW`, `Pseudo-R²`)
  } else {
    out %>% dplyr::select(System, Pair, n,
                           `R² (log)`, `Slope [95% CI]`,
                           `Concordance`, `Spearman ρ`, `Kendall τ`,
                           `SW p`, `BP p`, `DW`, `Pseudo-R²`)
  }
}

# --- Table 9: Intra-Habitat Level ---
table9_raw <- dplyr::bind_rows(
  res_to_row(res_faird_maize,  "FAIRD1 vs FAIRD2", "Maize",  "Intra-Habitat", "FAIRD"),
  res_to_row(res_faird_meadow, "FAIRD3 vs FAIRD4", "Meadow", "Intra-Habitat", "FAIRD"),
  res_to_row(res_ammod_maize,  "AMMOD1 vs AMMOD2", "Maize",  "Intra-Habitat", "AMMOD")
)

table9_display <- format_table(table9_raw, include_habitat = TRUE)

cat("=== TABLE 9: H1 Intra-Habitat Level (§3.2.1.1) ===\n\n")
print(table9_display, n = Inf)
cat("\nSignificance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p>=0.05\n")
cat("R² significance from slope p-value.\n")
cat("Spearman ρ on absolute log values (manuscript tables). Kendall τ on raw daily changes.\n\n")

readr::write_csv(table9_display,
                 "outputs/tables/Table09_H1_intra_habitat.csv")
readr::write_csv(
  table9_raw %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))),
  "outputs/tables/Table09_H1_intra_habitat_numeric.csv"
)
cat("[OK] Table09_H1_intra_habitat.csv\n")
cat("[OK] Table09_H1_intra_habitat_numeric.csv\n\n")

# --- Table 10: Inter-Habitat Level ---
table10_raw <- dplyr::bind_rows(
  res_to_row(res_faird_inter, "Maize mean vs Meadow mean", "Maize/Meadow", "Inter-Habitat", "FAIRD"),
  res_to_row(res_ammod_a1a3,  "AMMOD1 vs AMMOD3",          "Maize/Meadow", "Inter-Habitat", "AMMOD"),
  res_to_row(res_ammod_a2a3,  "AMMOD2 vs AMMOD3",          "Maize/Meadow", "Inter-Habitat", "AMMOD")
)

table10_display <- format_table(table10_raw, include_habitat = FALSE)

cat("=== TABLE 10: H1 Inter-Habitat Level (§3.2.1.3) ===\n\n")
print(table10_display, n = Inf)
cat("\nSignificance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p>=0.05\n")
cat("FAIRD: daily mean(FAIRD1+FAIRD2) [Maize] vs mean(FAIRD3+FAIRD4) [Meadow].\n")
cat("       Averaging used (not pooling) — site-pair covariate p = 0.013.\n")
cat("AMMOD: two separate pairs — no averaging because AMMOD4 failed.\n\n")

readr::write_csv(table10_display,
                 "outputs/tables/Table10_H1_inter_habitat.csv")
readr::write_csv(
  table10_raw %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))),
  "outputs/tables/Table10_H1_inter_habitat_numeric.csv"
)
cat("[OK] Table10_H1_inter_habitat.csv\n")
cat("[OK] Table10_H1_inter_habitat_numeric.csv\n\n")

# --- Table 8: H1 Overview (all levels, compact) ---
table8_raw <- dplyr::bind_rows(
  table9_raw  %>% dplyr::mutate(Level = "Intra-Habitat"),
  table10_raw %>% dplyr::mutate(Level = "Inter-Habitat")
)

table8_display <- table8_raw %>%
  dplyr::mutate(
    `R² (log)`     = paste0(formatC(R2, digits = 3, format = "f"),
                                  " ", sig_stars(Slope_p)),
    `Slope [95% CI]` = sprintf("%.2f [%.2f, %.2f]",
                                Slope, Slope_CI_low, Slope_CI_high),
    `Concordance`    = sprintf("%.1f%%", Concordance_pct),
    `Kendall τ` = paste0(formatC(Kendall_tau, digits = 3, format = "f"),
                                " ", sig_stars(Kendall_p)),
    `Pseudo-R²` = formatC(Pseudo_R2, digits = 3, format = "f")
  ) %>%
  dplyr::select(System, Level, Pair, Habitat, n,
                `R² (log)`, `Slope [95% CI]`,
                `Concordance`, `Kendall τ`, `Pseudo-R²`)

cat("=== TABLE 8: H1 Overview (§3.2.1 all levels) ===\n\n")
print(table8_display, n = Inf)
cat("\nSignificance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p>=0.05\n\n")

readr::write_csv(table8_display,
                 "outputs/tables/Table08_H1_overview.csv")
readr::write_csv(
  table8_raw %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))),
  "outputs/tables/Table08_H1_overview_numeric.csv"
)
cat("[OK] Table08_H1_overview.csv\n")
cat("[OK] Table08_H1_overview_numeric.csv\n\n")


# ==============================================================================
# SCRIPT COMPLETE
# ==============================================================================

cat("==============================================================================\n")
cat("H1_quanti_consistency.R — COMPLETE\n")
cat("Sections 1-7 (Figs 3-11, Tables 8-10) generated successfully.\n")
cat("==============================================================================\n")
cat("\nOutputs summary:\n")
cat("  Figures:  outputs/figures/Fig03 through Fig11\n")
cat("  Tables:   outputs/tables/Table08, Table09 (+ supp), Table10 (+ numeric copies)\n")
cat("  Console:  outputs/console/H1_quanti_consistency_output.txt\n\n")
cat("Script completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

while (sink.number() > 0) sink()
