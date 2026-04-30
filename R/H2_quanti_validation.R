# ==============================================================================
# H2_quanti_validation.R  — §3.3.1 Quantitative Inter-System Agreement
# Site Level: §3.3.1.1 (Figs 16–18) | Intra-Habitat: §3.3.1.2 (Fig 19)
# Overall: §3.3.1.3 (Fig 20) | CCF: §3.3.1.4 (Figs 21–24)
# Detection sensitivity: §3.3.1.5 (Table 14) | Body size: §3.3.1.6 (Table 15)
# All results assembled into Table 13.
# ==============================================================================

library(tidyverse)
library(here)
library(car)
library(lmtest)
library(patchwork)
library(scales)
library(MASS)    # load last — masks dplyr::select if not namespaced

here::i_am("R/H2_quanti_validation.R")

dir.create("outputs/figures",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",   showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console",  showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/H2_quanti_validation_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("H2 QUANTITATIVE VALIDATION — §3.3.1\n")
cat("log10(x+1) | AMMOD=X (predictor), FAIRD=Y (response)\n")
cat("Spearman & Kendall on daily changes | Bootstrap R² (n=1000, seed=42)\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source("R/00_load_data.R")
source("R/00_custom_colors.R")


# ==============================================================================
# SECTION 0: HELPERS AND CORE ANALYSIS FUNCTIONS
# ==============================================================================

etraps_daily      <- data_list$etraps_daily
ammod_biomass     <- data_list$ammod_biomass_for_stats

sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Build per-site analysis data: inner_join FAIRD biomass + AMMOD biomass
build_site_df <- function(site_id) {
  faird_site <- etraps_daily %>%
    dplyr::filter(Device_type == "FAIRD", Site == site_id) %>%
    dplyr::group_by(Date, Site) %>%
    dplyr::summarise(faird_biomass = mean(biomass_mg, na.rm = TRUE), .groups = "drop")

  ammod_site <- ammod_biomass %>%
    dplyr::filter(Site == site_id) %>%
    dplyr::select(Date, Site, Live_mass)

  dplyr::inner_join(faird_site, ammod_site, by = c("Date", "Site")) %>%
    dplyr::arrange(Date) %>%
    dplyr::mutate(
      log_faird = log10(faird_biomass + 1),
      log_ammod = log10(Live_mass + 1)
    )
}

# OLS + full diagnostics for one site df
analyze_site_ols <- function(df, site_id) {
  n <- nrow(df)
  model <- stats::lm(log_faird ~ log_ammod, data = df)
  sm    <- summary(model)

  slope     <- stats::coef(model)[["log_ammod"]]
  intercept <- stats::coef(model)[["(Intercept)"]]
  slope_ci  <- stats::confint(model)["log_ammod", ]
  slope_p   <- sm$coefficients["log_ammod", "Pr(>|t|)"]
  intercept_p <- sm$coefficients["(Intercept)", "Pr(>|t|)"]
  r_sq      <- sm$r.squared
  adj_r_sq  <- sm$adj.r.squared

  # Diagnostics
  sw_p  <- stats::shapiro.test(stats::residuals(model))$p.value
  bp_p  <- car::ncvTest(model)$p
  dw    <- lmtest::dwtest(model)
  dw_stat <- dw$statistic[[1]]
  dw_p    <- dw$p.value

  # Spearman on daily log changes (per CLAUDE.md: on changes, not absolute values)
  delta_log_faird <- diff(df$log_faird)
  delta_log_ammod <- diff(df$log_ammod)
  sp <- stats::cor.test(delta_log_faird, delta_log_ammod, method = "spearman", exact = FALSE)
  spearman_rho <- sp$estimate[[1]]
  spearman_p   <- sp$p.value

  # Kendall and concordance on raw biomass changes (matching run_temporal_sync in source)
  # Log transformation reorders change magnitudes, shifting concordant/discordant counts.
  delta_raw_faird <- diff(df$faird_biomass)
  delta_raw_ammod <- diff(df$Live_mass)
  kt <- stats::cor.test(delta_raw_faird, delta_raw_ammod, method = "kendall", exact = FALSE)
  kendall_tau <- kt$estimate[[1]]
  kendall_p   <- kt$p.value

  # Concordance: sign agreement on raw daily changes
  sign_match  <- sign(delta_raw_faird) == sign(delta_raw_ammod)
  concordance <- mean(sign_match, na.rm = TRUE) * 100

  # Robust regression (MASS::rlm)
  robust_model <- MASS::rlm(log_faird ~ log_ammod, data = df, psi = psi.huber)
  pseudo_r2    <- 1 - sum(robust_model$residuals^2) / sum((df$log_faird - mean(df$log_faird))^2)
  weights_rob  <- robust_model$w
  n_low_wt     <- sum(weights_rob < 0.5)

  list(
    site = site_id, n = n,
    r_squared = r_sq, adj_r_squared = adj_r_sq,
    slope = slope, slope_ci_low = slope_ci[1], slope_ci_high = slope_ci[2],
    slope_p = slope_p, intercept = intercept, intercept_p = intercept_p,
    sw_p = sw_p, bp_p = bp_p, dw_stat = dw_stat, dw_p = dw_p,
    spearman_rho = spearman_rho, spearman_p = spearman_p,
    kendall_tau = kendall_tau, kendall_p = kendall_p,
    concordance_pct = concordance,
    pseudo_r2 = pseudo_r2, n_low_wt = n_low_wt,
    model = model, df = df
  )
}

# Scatter plot: AMMOD (x) vs FAIRD (y), log10(x+1) scale
make_scatter_h2 <- function(df, res, site_label, point_color) {
  r2_label <- sprintf("R² = %.3f%s", res$r_squared, sig_stars(res$slope_p))
  sp_label <- sprintf("ρ = %.3f%s", res$spearman_rho, sig_stars(res$spearman_p))
  n_label  <- sprintf("n = %d days", res$n)

  ggplot2::ggplot(df, ggplot2::aes(x = log_ammod, y = log_faird)) +
    ggplot2::geom_point(color = point_color, size = 2.5, alpha = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black",
                         fill = "grey80", linewidth = 0.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                         linetype = "dashed", color = "grey50", linewidth = 0.6) +
    ggplot2::annotate("text", x = -Inf, y = Inf,
                      label = paste(r2_label, sp_label, n_label, sep = "\n"),
                      hjust = -0.08, vjust = 1.3, size = 3.2) +
    ggplot2::labs(
      title   = site_label,
      x       = "AMMOD log10(Live mass + 1) [mg]",
      y       = "FAIRD log10(Biomass + 1) [mg]"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(face = "bold", size = 11),
      panel.grid.minor  = ggplot2::element_blank()
    )
}

# Temporal trend: FAIRD and AMMOD log biomass over time
make_temporal_h2 <- function(df, site_label, point_color) {
  df_long <- df %>%
    dplyr::select(Date, log_faird, log_ammod) %>%
    tidyr::pivot_longer(c(log_faird, log_ammod),
                        names_to = "System", values_to = "log_biomass") %>%
    dplyr::mutate(
      System = dplyr::recode(System,
                             "log_faird" = "FAIRD",
                             "log_ammod" = "AMMOD")
    )

  sys_colors <- c(
    "FAIRD" = point_color,
    "AMMOD" = custom_colors[["AMMOD1"]]
  )

  ggplot2::ggplot(df_long, ggplot2::aes(x = Date, y = log_biomass,
                                         color = System, group = System)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = sys_colors) +
    ggplot2::labs(
      title = site_label,
      x     = "Date",
      y     = "log10(Biomass + 1) [mg]",
      color = "System"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(face = "bold", size = 11),
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "top"
    )
}


# ==============================================================================
# SECTION 1: SITE LEVEL — §3.3.1.1 (Figs 16–18, Table 13 site rows)
# Source: H2_validation_FAIRD_AMMOD.R
# Site4 excluded: AMMOD4 had mechanical failure (no biomass data).
# inner_join naturally excludes Site4 (no ammod rows).
# ==============================================================================

cat("\n==============================================================================\n")
cat("SECTION 1: SITE LEVEL OLS — §3.3.1.1\n")
cat("==============================================================================\n\n")

# Build per-site data frames
df_site1 <- build_site_df("Site1")
df_site2 <- build_site_df("Site2")
df_site3 <- build_site_df("Site3")

cat(sprintf("Site1: n=%d days | Site2: n=%d days | Site3 (Meadow): n=%d days\n\n",
            nrow(df_site1), nrow(df_site2), nrow(df_site3)))

# Run OLS + diagnostics
res_site1 <- analyze_site_ols(df_site1, "Site1")
res_site2 <- analyze_site_ols(df_site2, "Site2")
res_site3 <- analyze_site_ols(df_site3, "Site3")

# Print site-level results
for (res in list(res_site1, res_site2, res_site3)) {
  cat(sprintf("--- %s (n=%d) ---\n", res$site, res$n))
  cat(sprintf("  R²=%.3f (adj=%.3f)  slope=%.3f [%.3f, %.3f]  slope_p=%.4f %s\n",
              res$r_squared, res$adj_r_squared,
              res$slope, res$slope_ci_low, res$slope_ci_high,
              res$slope_p, sig_stars(res$slope_p)))
  cat(sprintf("  Spearman ρ=%.3f (p=%.4f %s)  Kendall τ=%.3f (p=%.4f %s)\n",
              res$spearman_rho, res$spearman_p, sig_stars(res$spearman_p),
              res$kendall_tau, res$kendall_p, sig_stars(res$kendall_p)))
  cat(sprintf("  Concordance=%.1f%%  Pseudo-R²=%.3f  n_low_wt=%d\n",
              res$concordance_pct, res$pseudo_r2, res$n_low_wt))
  cat(sprintf("  Diagnostics: SW p=%.3f | BP p=%.3f | DW stat=%.3f (p=%.3f)\n\n",
              res$sw_p, res$bp_p, res$dw_stat, res$dw_p))
}

# --- Figure 16: Site1 ---
p16_scatter  <- make_scatter_h2(df_site1, res_site1, "Site 1 — Maize",
                                 custom_colors[["Site1"]])
p16_temporal <- make_temporal_h2(df_site1, "Site 1 — Daily Biomass Trends",
                                  custom_colors[["Site1"]])
fig16 <- p16_scatter / p16_temporal +
  patchwork::plot_annotation(
    title    = "Figure 16 — H2 Site Level: Site 1 (Maize)",
    subtitle = "AMMOD (predictor) vs FAIRD (response) | log10(x+1)",
    theme    = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40")
    )
  )
ggplot2::ggsave("outputs/figures/Fig16_H2_Site1.png", fig16, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: Fig16_H2_Site1.png\n")

# --- Figure 17: Site2 ---
p17_scatter  <- make_scatter_h2(df_site2, res_site2, "Site 2 — Maize",
                                 custom_colors[["Site2"]])
p17_temporal <- make_temporal_h2(df_site2, "Site 2 — Daily Biomass Trends",
                                  custom_colors[["Site2"]])
fig17 <- p17_scatter / p17_temporal +
  patchwork::plot_annotation(
    title    = "Figure 17 — H2 Site Level: Site 2 (Maize)",
    subtitle = "AMMOD (predictor) vs FAIRD (response) | log10(x+1)",
    theme    = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40")
    )
  )
ggplot2::ggsave("outputs/figures/Fig17_H2_Site2.png", fig17, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: Fig17_H2_Site2.png\n")

# --- Figure 18: Site3 (Meadow) ---
p18_scatter  <- make_scatter_h2(df_site3, res_site3, "Site 3 — Meadow",
                                 custom_colors[["Site3"]])
p18_temporal <- make_temporal_h2(df_site3, "Site 3 — Daily Biomass Trends",
                                  custom_colors[["Site3"]])
fig18 <- p18_scatter / p18_temporal +
  patchwork::plot_annotation(
    title    = "Figure 18 — H2 Site Level: Site 3 (Meadow)",
    subtitle = "AMMOD (predictor) vs FAIRD (response) | log10(x+1)",
    theme    = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40")
    )
  )
ggplot2::ggsave("outputs/figures/Fig18_H2_Site3.png", fig18, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: Fig18_H2_Site3.png\n\n")


# ==============================================================================
# SECTION 2: INTRA-HABITAT LEVEL — MAIZE — §3.3.1.2 (Fig 19, Table 13 Maize row)
# Source: H2_validation_aggregated_levels.R
# Pool Site1+Site2 (n=24). Covariate validation: Site covariate p=0.861 → pooling valid.
# Spearman/Kendall: daily changes computed within-site before pooling.
# Bootstrap R²: n_boot=1000, seed=42, quantile CI method.
# ==============================================================================

cat("==============================================================================\n")
cat("SECTION 2: INTRA-HABITAT LEVEL (MAIZE) — §3.3.1.2\n")
cat("==============================================================================\n\n")

df_maize <- dplyr::bind_rows(df_site1, df_site2) %>%
  dplyr::arrange(Site, Date)

cat(sprintf("Maize pooled: n=%d days (Site1=%d, Site2=%d)\n\n",
            nrow(df_maize), nrow(df_site1), nrow(df_site2)))

# Covariate model — validate pooling
covariate_maize <- stats::lm(log_faird ~ log_ammod + Site, data = df_maize)
sm_cov <- summary(covariate_maize)
cat("Covariate model (log_faird ~ log_ammod + Site):\n")
print(sm_cov$coefficients)
site2_p <- tryCatch(sm_cov$coefficients["SiteS2", "Pr(>|t|)"],
                    error = function(e) sm_cov$coefficients[grep("Site2", rownames(sm_cov$coefficients)), "Pr(>|t|)"][1])
# Handle possible row name formats
coeff_names <- rownames(sm_cov$coefficients)
site_coef_rows <- coeff_names[grep("Site", coeff_names)]
if (length(site_coef_rows) > 0) {
  site2_p <- sm_cov$coefficients[site_coef_rows[1], "Pr(>|t|)"]
  cat(sprintf("\nSite covariate (%s) p=%.3f → pooling %s\n\n",
              site_coef_rows[1], site2_p,
              ifelse(site2_p > 0.05, "VALID", "QUESTIONABLE")))
}

# Pooled OLS
model_maize <- stats::lm(log_faird ~ log_ammod, data = df_maize)
sm_maize    <- summary(model_maize)

r2_maize      <- sm_maize$r.squared
slope_maize   <- stats::coef(model_maize)[["log_ammod"]]
slope_ci_m    <- stats::confint(model_maize)["log_ammod", ]
slope_p_m     <- sm_maize$coefficients["log_ammod", "Pr(>|t|)"]
intercept_m   <- stats::coef(model_maize)[["(Intercept)"]]
intercept_p_m <- sm_maize$coefficients["(Intercept)", "Pr(>|t|)"]

sw_maize   <- stats::shapiro.test(stats::residuals(model_maize))$p.value
bp_maize   <- car::ncvTest(model_maize)$p
dw_maize   <- lmtest::dwtest(model_maize)
dw_stat_m  <- dw_maize$statistic[[1]]
dw_p_m     <- dw_maize$p.value

# Pseudo-R² (robust)
robust_maize  <- MASS::rlm(log_faird ~ log_ammod, data = df_maize, psi = psi.huber)
pseudo_r2_m   <- 1 - sum(robust_maize$residuals^2) / sum((df_maize$log_faird - mean(df_maize$log_faird))^2)
n_low_wt_m    <- sum(robust_maize$w < 0.5)

# Spearman/Kendall on within-site daily changes
# Spearman: log changes | Kendall + concordance: raw changes (matching run_temporal_sync)
changes_maize <- df_maize %>%
  dplyr::arrange(Site, Date) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    d_faird_log = log_faird     - dplyr::lag(log_faird),
    d_ammod_log = log_ammod     - dplyr::lag(log_ammod),
    d_faird_raw = faird_biomass - dplyr::lag(faird_biomass),
    d_ammod_raw = Live_mass     - dplyr::lag(Live_mass)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na(d_faird_raw, d_ammod_raw)

sp_m   <- stats::cor.test(changes_maize$d_faird_log, changes_maize$d_ammod_log,
                           method = "spearman", exact = FALSE)
kt_m   <- stats::cor.test(changes_maize$d_faird_raw, changes_maize$d_ammod_raw,
                           method = "kendall",  exact = FALSE)
conc_m <- mean(sign(changes_maize$d_faird_raw) == sign(changes_maize$d_ammod_raw)) * 100

# Bootstrap R²
set.seed(42)
boot_r2_maize <- replicate(1000, {
  idx <- sample(nrow(df_maize), replace = TRUE)
  b   <- stats::lm(log_faird ~ log_ammod, data = df_maize[idx, ])
  summary(b)$r.squared
})
boot_ci_m <- stats::quantile(boot_r2_maize, c(0.025, 0.975))

cat(sprintf("Maize pooled OLS: R²=%.3f%s [Boot 95%% CI: %.3f–%.3f]\n",
            r2_maize, sig_stars(slope_p_m), boot_ci_m[1], boot_ci_m[2]))
cat(sprintf("  slope=%.3f [%.3f, %.3f]  Spearman ρ=%.3f%s  Kendall τ=%.3f%s\n",
            slope_maize, slope_ci_m[1], slope_ci_m[2],
            sp_m$estimate[[1]], sig_stars(sp_m$p.value),
            kt_m$estimate[[1]], sig_stars(kt_m$p.value)))
cat(sprintf("  Concordance=%.1f%%  Pseudo-R²=%.3f  n_low_wt=%d\n",
            conc_m, pseudo_r2_m, n_low_wt_m))
cat(sprintf("  Diagnostics: SW p=%.3f | BP p=%.3f | DW stat=%.3f (p=%.3f)\n\n",
            sw_maize, bp_maize, dw_stat_m, dw_p_m))

# --- Figure 19: Maize pooled scatter ---
site_colors_maize <- c(
  "Site1" = custom_colors[["Site1"]],
  "Site2" = custom_colors[["Site2"]]
)

r2_label_m <- sprintf("R² = %.3f%s [Boot CI: %.3f–%.3f]",
                       r2_maize, sig_stars(slope_p_m), boot_ci_m[1], boot_ci_m[2])
sp_label_m <- sprintf("ρ = %.3f%s", sp_m$estimate[[1]], sig_stars(sp_m$p.value))

fig19 <- ggplot2::ggplot(df_maize, ggplot2::aes(x = log_ammod, y = log_faird, color = Site)) +
  ggplot2::geom_point(size = 2.5, alpha = 0.8) +
  ggplot2::geom_smooth(ggplot2::aes(group = 1),
                       method = "lm", se = TRUE, color = "black",
                       fill = "grey80", linewidth = 0.8) +
  ggplot2::geom_abline(intercept = 0, slope = 1,
                       linetype = "dashed", color = "grey50", linewidth = 0.6) +
  ggplot2::scale_color_manual(values = site_colors_maize) +
  ggplot2::annotate("text", x = -Inf, y = Inf,
                    label = paste(r2_label_m, sp_label_m,
                                  sprintf("n = %d days", nrow(df_maize)), sep = "\n"),
                    hjust = -0.08, vjust = 1.3, size = 3.2) +
  ggplot2::labs(
    title    = "Figure 19 — H2 Intra-Habitat Level: Maize (Sites 1+2 pooled)",
    subtitle = "AMMOD (predictor) vs FAIRD (response) | log10(x+1) | Site covariate p=0.861",
    x        = "AMMOD log10(Live mass + 1) [mg]",
    y        = "FAIRD log10(Biomass + 1) [mg]",
    color    = "Site"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title        = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle     = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor  = ggplot2::element_blank(),
    legend.position   = "top"
  )

ggplot2::ggsave("outputs/figures/Fig19_H2_Maize_Pooled.png", fig19,
                width = 8, height = 6, dpi = 300, bg = "white")
cat("Saved: Fig19_H2_Maize_Pooled.png\n\n")


# ==============================================================================
# SECTION 3: OVERALL LEVEL — §3.3.1.3 (Fig 20, Table 13 Overall row)
# Source: H2_validation_aggregated_levels.R
# All sites pooled (n=32). Site covariate all p>0.64 → pooling valid.
# Primary result per Scherber request: R²=0.321***.
# ==============================================================================

cat("==============================================================================\n")
cat("SECTION 3: OVERALL LEVEL — §3.3.1.3\n")
cat("==============================================================================\n\n")

df_overall <- dplyr::bind_rows(df_site1, df_site2, df_site3) %>%
  dplyr::arrange(Site, Date)

cat(sprintf("Overall pooled: n=%d days (Site1=%d, Site2=%d, Site3=%d)\n\n",
            nrow(df_overall), nrow(df_site1), nrow(df_site2), nrow(df_site3)))

# Covariate model — validate pooling
covariate_overall <- stats::lm(log_faird ~ log_ammod + Site, data = df_overall)
sm_cov_ov <- summary(covariate_overall)
cat("Covariate model (log_faird ~ log_ammod + Site) — Overall:\n")
print(sm_cov_ov$coefficients)
cat("\n")

# Pooled OLS
model_overall  <- stats::lm(log_faird ~ log_ammod, data = df_overall)
sm_overall     <- summary(model_overall)

r2_overall     <- sm_overall$r.squared
slope_ov       <- stats::coef(model_overall)[["log_ammod"]]
slope_ci_ov    <- stats::confint(model_overall)["log_ammod", ]
slope_p_ov     <- sm_overall$coefficients["log_ammod", "Pr(>|t|)"]
intercept_ov   <- stats::coef(model_overall)[["(Intercept)"]]
intercept_p_ov <- sm_overall$coefficients["(Intercept)", "Pr(>|t|)"]

sw_ov   <- stats::shapiro.test(stats::residuals(model_overall))$p.value
bp_ov   <- car::ncvTest(model_overall)$p
dw_ov   <- lmtest::dwtest(model_overall)
dw_stat_ov <- dw_ov$statistic[[1]]
dw_p_ov    <- dw_ov$p.value

# Pseudo-R² (robust)
robust_ov  <- MASS::rlm(log_faird ~ log_ammod, data = df_overall, psi = psi.huber)
pseudo_r2_ov  <- 1 - sum(robust_ov$residuals^2) / sum((df_overall$log_faird - mean(df_overall$log_faird))^2)
n_low_wt_ov   <- sum(robust_ov$w < 0.5)

# Spearman/Kendall on within-site daily changes
# Spearman: log changes | Kendall + concordance: raw changes (matching run_temporal_sync)
changes_overall <- df_overall %>%
  dplyr::arrange(Site, Date) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    d_faird_log = log_faird     - dplyr::lag(log_faird),
    d_ammod_log = log_ammod     - dplyr::lag(log_ammod),
    d_faird_raw = faird_biomass - dplyr::lag(faird_biomass),
    d_ammod_raw = Live_mass     - dplyr::lag(Live_mass)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na(d_faird_raw, d_ammod_raw)

sp_ov   <- stats::cor.test(changes_overall$d_faird_log, changes_overall$d_ammod_log,
                            method = "spearman", exact = FALSE)
kt_ov   <- stats::cor.test(changes_overall$d_faird_raw, changes_overall$d_ammod_raw,
                            method = "kendall",  exact = FALSE)
conc_ov <- mean(sign(changes_overall$d_faird_raw) == sign(changes_overall$d_ammod_raw)) * 100

# Bootstrap R²
set.seed(42)
boot_r2_overall <- replicate(1000, {
  idx <- sample(nrow(df_overall), replace = TRUE)
  b   <- stats::lm(log_faird ~ log_ammod, data = df_overall[idx, ])
  summary(b)$r.squared
})
boot_ci_ov <- stats::quantile(boot_r2_overall, c(0.025, 0.975))

cat(sprintf("Overall OLS: R²=%.3f%s [Boot 95%% CI: %.3f–%.3f]\n",
            r2_overall, sig_stars(slope_p_ov), boot_ci_ov[1], boot_ci_ov[2]))
cat(sprintf("  slope=%.3f [%.3f, %.3f]  Spearman ρ=%.3f%s  Kendall τ=%.3f%s\n",
            slope_ov, slope_ci_ov[1], slope_ci_ov[2],
            sp_ov$estimate[[1]], sig_stars(sp_ov$p.value),
            kt_ov$estimate[[1]], sig_stars(kt_ov$p.value)))
cat(sprintf("  Concordance=%.1f%%  Pseudo-R²=%.3f  n_low_wt=%d\n",
            conc_ov, pseudo_r2_ov, n_low_wt_ov))
cat(sprintf("  Diagnostics: SW p=%.3f | BP p=%.3f | DW stat=%.3f (p=%.3f)\n\n",
            sw_ov, bp_ov, dw_stat_ov, dw_p_ov))

# --- Figure 20: Overall pooled scatter ---
site_colors_ov <- c(
  "Site1" = custom_colors[["Site1"]],
  "Site2" = custom_colors[["Site2"]],
  "Site3" = custom_colors[["Site3"]]
)

r2_label_ov <- sprintf("R² = %.3f%s [Boot CI: %.3f–%.3f]",
                        r2_overall, sig_stars(slope_p_ov), boot_ci_ov[1], boot_ci_ov[2])
sp_label_ov <- sprintf("ρ = %.3f%s", sp_ov$estimate[[1]], sig_stars(sp_ov$p.value))

fig20 <- ggplot2::ggplot(df_overall, ggplot2::aes(x = log_ammod, y = log_faird, color = Site)) +
  ggplot2::geom_point(size = 2.5, alpha = 0.8) +
  ggplot2::geom_smooth(ggplot2::aes(group = 1),
                       method = "lm", se = TRUE, color = "black",
                       fill = "grey80", linewidth = 0.8) +
  ggplot2::geom_abline(intercept = 0, slope = 1,
                       linetype = "dashed", color = "grey50", linewidth = 0.6) +
  ggplot2::scale_color_manual(values = site_colors_ov) +
  ggplot2::annotate("text", x = -Inf, y = Inf,
                    label = paste(r2_label_ov, sp_label_ov,
                                  sprintf("n = %d days", nrow(df_overall)), sep = "\n"),
                    hjust = -0.08, vjust = 1.3, size = 3.2) +
  ggplot2::labs(
    title    = "Figure 20 — H2 Overall Level: All Sites Pooled",
    subtitle = "AMMOD (predictor) vs FAIRD (response) | log10(x+1) | Site covariate all p > 0.64",
    x        = "AMMOD log10(Live mass + 1) [mg]",
    y        = "FAIRD log10(Biomass + 1) [mg]",
    color    = "Site"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title        = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle     = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor  = ggplot2::element_blank(),
    legend.position   = "top"
  )

ggplot2::ggsave("outputs/figures/Fig20_H2_Overall_Pooled.png", fig20,
                width = 8, height = 6, dpi = 300, bg = "white")
cat("Saved: Fig20_H2_Overall_Pooled.png\n\n")


# ==============================================================================
# SECTION 4: CROSS-CORRELATION ANALYSIS — §3.3.1.4 (Figs 21–24, Table 13 CCF)
# Source: H2_validation_cross_corr.R
# Method: Z-scored series, acf(cbind()), CI = 1.96/sqrt(n) (significance threshold).
# Per-site: Figs 21–23. Intra-Habitat Maize summed series: Fig 24 (PRIMARY).
# Fisher's z averaging: descriptive only. Overall Level: footnote only.
# ==============================================================================

cat("==============================================================================\n")
cat("SECTION 4: CROSS-CORRELATION ANALYSIS — §3.3.1.4\n")
cat("==============================================================================\n\n")

# Build time series data: log10(x+1) biomass per site per day (FAIRD only)
ts_data <- dplyr::inner_join(
  etraps_daily %>%
    dplyr::filter(Device_type == "FAIRD") %>%
    dplyr::group_by(Date, Site) %>%
    dplyr::summarise(faird_biomass = mean(biomass_mg, na.rm = TRUE), .groups = "drop"),
  ammod_biomass %>% dplyr::select(Date, Site, Live_mass),
  by = c("Date", "Site")
) %>%
  dplyr::mutate(
    faird_value = log10(faird_biomass + 1),
    ammod_value = log10(Live_mass + 1)
  ) %>%
  dplyr::arrange(Site, Date)

# CCF helper — Z-scored series, lag.max = min(10, floor(n/3))
analyze_ccf <- function(df_site, label) {
  n        <- nrow(df_site)
  z_faird  <- scale(df_site$faird_value)[, 1]
  z_ammod  <- scale(df_site$ammod_value)[, 1]
  lag_max  <- min(10L, floor(n / 3L))
  ci_bound <- stats::qnorm(0.975) / sqrt(n)

  ccf_obj   <- stats::acf(cbind(z_faird, z_ammod), lag.max = lag_max,
                           type = "correlation", plot = FALSE)
  # acf[,1,2] = z_faird → z_ammod cross-correlation; lag 0 at index 1
  ccf_vals  <- as.numeric(ccf_obj$acf[, 1, 2])
  lags      <- as.integer(ccf_obj$lag[, 1, 2])
  lag0_r    <- ccf_vals[lags == 0]

  list(
    label     = label,
    n         = n,
    lag0_r    = lag0_r,
    ci_bound  = ci_bound,
    significant = abs(lag0_r) > ci_bound,
    ccf_vals  = ccf_vals,
    lags      = lags,
    lag_max   = lag_max
  )
}

# CCF bar chart helper
plot_ccf_bars <- function(ccf_res, site_color) {
  df_ccf <- data.frame(lag = ccf_res$lags, r = ccf_res$ccf_vals)
  ci     <- ccf_res$ci_bound
  lag0_r <- ccf_res$lag0_r

  sig_label <- sprintf("lag-0 r = %.3f (%s)\nCI bound = ±%.3f",
                        lag0_r,
                        ifelse(ccf_res$significant, "p < 0.05", "ns"),
                        ci)

  ggplot2::ggplot(df_ccf, ggplot2::aes(x = lag, y = r)) +
    ggplot2::geom_col(fill = site_color, alpha = 0.8, width = 0.7) +
    ggplot2::geom_hline(yintercept =  ci, linetype = "dashed", color = "firebrick", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = -ci, linetype = "dashed", color = "firebrick", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::annotate("text", x = Inf, y = Inf,
                      label = sig_label, hjust = 1.1, vjust = 1.3, size = 3) +
    ggplot2::scale_x_continuous(breaks = df_ccf$lag) +
    ggplot2::labs(x = "Lag (days)", y = "Cross-correlation (r)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# Stacked time series helper (Z-scored)
plot_ts_stack <- function(df_site, site_color) {
  df_long <- df_site %>%
    dplyr::mutate(
      z_faird = scale(faird_value)[, 1],
      z_ammod = scale(ammod_value)[, 1]
    ) %>%
    dplyr::select(Date, z_faird, z_ammod) %>%
    tidyr::pivot_longer(c(z_faird, z_ammod),
                        names_to = "System", values_to = "z") %>%
    dplyr::mutate(System = dplyr::recode(System,
                                          "z_faird" = "FAIRD",
                                          "z_ammod" = "AMMOD"))

  sys_colors <- c("FAIRD" = site_color, "AMMOD" = custom_colors[["AMMOD1"]])

  ggplot2::ggplot(df_long, ggplot2::aes(x = Date, y = z, color = System, group = System)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_color_manual(values = sys_colors) +
    ggplot2::labs(x = "Date", y = "Z-score", color = "System") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )
}

# --- Per-site CCF ---
df_ts_s1 <- ts_data %>% dplyr::filter(Site == "Site1") %>% dplyr::arrange(Date)
df_ts_s2 <- ts_data %>% dplyr::filter(Site == "Site2") %>% dplyr::arrange(Date)
df_ts_s3 <- ts_data %>% dplyr::filter(Site == "Site3") %>% dplyr::arrange(Date)

ccf_s1 <- analyze_ccf(df_ts_s1, "Site1")
ccf_s2 <- analyze_ccf(df_ts_s2, "Site2")
ccf_s3 <- analyze_ccf(df_ts_s3, "Site3")

cat("Per-site CCF (lag-0):\n")
for (cr in list(ccf_s1, ccf_s2, ccf_s3)) {
  cat(sprintf("  %s: r=%.3f  CI±%.3f  %s  (n=%d)\n",
              cr$label, cr$lag0_r, cr$ci_bound,
              ifelse(cr$significant, "SIGNIFICANT", "ns"), cr$n))
}
cat("\n")

# Figure 21 — Site1 CCF
p21_ts  <- plot_ts_stack(df_ts_s1, custom_colors[["Site1"]]) +
  ggplot2::ggtitle("Site 1 — Z-scored Daily Biomass")
p21_ccf <- plot_ccf_bars(ccf_s1, custom_colors[["Site1"]]) +
  ggplot2::ggtitle(sprintf("Site 1 CCF — lag-0 r = %.3f %s",
                             ccf_s1$lag0_r,
                             ifelse(ccf_s1$significant, "(*)", "(ns)")))
fig21 <- p21_ts / p21_ccf +
  patchwork::plot_annotation(title = "Figure 21 — H2 CCF: Site 1 (Maize)",
                              theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
ggplot2::ggsave("outputs/figures/Fig21_H2_CCF_Site1.png", fig21,
                width = 9, height = 7, dpi = 300, bg = "white")
cat("Saved: Fig21_H2_CCF_Site1.png\n")

# Figure 22 — Site2 CCF
p22_ts  <- plot_ts_stack(df_ts_s2, custom_colors[["Site2"]]) +
  ggplot2::ggtitle("Site 2 — Z-scored Daily Biomass")
p22_ccf <- plot_ccf_bars(ccf_s2, custom_colors[["Site2"]]) +
  ggplot2::ggtitle(sprintf("Site 2 CCF — lag-0 r = %.3f %s",
                             ccf_s2$lag0_r,
                             ifelse(ccf_s2$significant, "(*)", "(ns)")))
fig22 <- p22_ts / p22_ccf +
  patchwork::plot_annotation(title = "Figure 22 — H2 CCF: Site 2 (Maize)",
                              theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
ggplot2::ggsave("outputs/figures/Fig22_H2_CCF_Site2.png", fig22,
                width = 9, height = 7, dpi = 300, bg = "white")
cat("Saved: Fig22_H2_CCF_Site2.png\n")

# Figure 23 — Site3 CCF
p23_ts  <- plot_ts_stack(df_ts_s3, custom_colors[["Site3"]]) +
  ggplot2::ggtitle("Site 3 — Z-scored Daily Biomass")
p23_ccf <- plot_ccf_bars(ccf_s3, custom_colors[["Site3"]]) +
  ggplot2::ggtitle(sprintf("Site 3 CCF — lag-0 r = %.3f %s",
                             ccf_s3$lag0_r,
                             ifelse(ccf_s3$significant, "(*)", "(ns)")))
fig23 <- p23_ts / p23_ccf +
  patchwork::plot_annotation(title = "Figure 23 — H2 CCF: Site 3 (Meadow)",
                              theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
ggplot2::ggsave("outputs/figures/Fig23_H2_CCF_Site3.png", fig23,
                width = 9, height = 7, dpi = 300, bg = "white")
cat("Saved: Fig23_H2_CCF_Site3.png\n\n")

# --- Intra-Habitat Maize: summed series (PRIMARY) ---
cat("Intra-Habitat Maize CCF — summed series (PRIMARY):\n")
ts_maize_summed <- ts_data %>%
  dplyr::filter(Site %in% c("Site1", "Site2")) %>%
  dplyr::group_by(Date) %>%
  dplyr::summarise(
    faird_value = sum(faird_value),
    ammod_value = sum(ammod_value),
    n_sites     = dplyr::n(),
    .groups     = "drop"
  ) %>%
  dplyr::filter(n_sites == 2) %>%   # both sites must have data
  dplyr::arrange(Date)

cat(sprintf("  Summed series: n=%d days (days where both Site1 & Site2 have data)\n",
            nrow(ts_maize_summed)))

ccf_maize_summed <- analyze_ccf(ts_maize_summed, "Maize_Summed")
cat(sprintf("  Maize summed CCF: lag-0 r=%.3f  CI±%.3f  %s\n\n",
            ccf_maize_summed$lag0_r, ccf_maize_summed$ci_bound,
            ifelse(ccf_maize_summed$significant, "SIGNIFICANT", "ns")))

# Fisher's z average (descriptive only — no formal CI)
z_s1  <- atanh(ccf_s1$lag0_r)
z_s2  <- atanh(ccf_s2$lag0_r)
r_avg <- tanh(mean(c(z_s1, z_s2)))
cat(sprintf("Fisher's z averaging (descriptive): r_avg=%.3f [atanh(%.3f)+atanh(%.3f))/2 → tanh]\n\n",
            r_avg, ccf_s1$lag0_r, ccf_s2$lag0_r))

# Overall CCF (descriptive, footnote in Table 13)
ccf_overall_desc <- analyze_ccf(
  ts_data %>% dplyr::arrange(Site, Date) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(faird_value = mean(faird_value),
                     ammod_value = mean(ammod_value),
                     .groups = "drop") %>%
    dplyr::arrange(Date),
  "Overall_descriptive"
)
cat(sprintf("Overall CCF (descriptive, excluded as primary — Site3 n=8 inconclusive):\n"))
cat(sprintf("  lag-0 r=%.3f  CI±%.3f  %s\n\n",
            ccf_overall_desc$lag0_r, ccf_overall_desc$ci_bound,
            ifelse(ccf_overall_desc$significant, "SIGNIFICANT", "ns")))

# Figure 24 — Maize summed CCF
p24_ts  <- plot_ts_stack(ts_maize_summed, custom_colors[["Site1"]]) +
  ggplot2::ggtitle("Maize — Z-scored Summed Daily Biomass (Sites 1+2)")
p24_ccf <- plot_ccf_bars(ccf_maize_summed, custom_colors[["Site1"]]) +
  ggplot2::ggtitle(sprintf("Maize Summed CCF — lag-0 r = %.3f %s",
                             ccf_maize_summed$lag0_r,
                             ifelse(ccf_maize_summed$significant, "(*)", "(ns)")))
fig24 <- p24_ts / p24_ccf +
  patchwork::plot_annotation(
    title    = "Figure 24 — H2 CCF: Intra-Habitat Maize (Summed Series, PRIMARY)",
    subtitle = sprintf("n=%d days | both sites present | formal CI valid",
                        nrow(ts_maize_summed)),
    theme    = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey40")
    )
  )
ggplot2::ggsave("outputs/figures/Fig24_H2_CCF_Maize_Summed.png", fig24,
                width = 9, height = 7, dpi = 300, bg = "white")
cat("Saved: Fig24_H2_CCF_Maize_Summed.png\n\n")


# ==============================================================================
# SECTION 5: DETECTION METRIC SENSITIVITY — §3.3.1.5 (Table 14)
# Source: H2_validation_counts_vs_drymass.R
# FIX: original reads intermediate CSVs → replaced with data_list$ objects.
# Comparison: raw abundance counts vs allometric biomass as detection metric.
# ==============================================================================

cat("==============================================================================\n")
cat("SECTION 5: DETECTION METRIC SENSITIVITY — §3.3.1.5\n")
cat("==============================================================================\n\n")

# Use data_list objects — DO NOT read intermediate CSVs
etraps_raw   <- data_list$etraps_daily          # has abundance + biomass_mg columns
ammod_raw    <- data_list$ammod_biomass_for_stats  # has Live_mass + Dry_mass_g columns

# Helper: run OLS + diagnostics, return tidy row
run_metric_ols <- function(df, x_col, y_col, site_id, metric_label) {
  df_fit <- df %>%
    dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]])) %>%
    dplyr::mutate(
      log_x = log10(.data[[x_col]] + 1),
      log_y = log10(.data[[y_col]] + 1)
    )
  if (nrow(df_fit) < 4) return(NULL)

  m    <- stats::lm(log_y ~ log_x, data = df_fit)
  sm   <- summary(m)
  r2   <- sm$r.squared
  sp   <- stats::cor.test(diff(df_fit$log_y), diff(df_fit$log_x),
                           method = "spearman", exact = FALSE)
  sw_p <- stats::shapiro.test(stats::residuals(m))$p.value
  bp_p <- car::ncvTest(m)$p
  dw   <- lmtest::dwtest(m)

  tibble::tibble(
    Site   = site_id,
    Metric = metric_label,
    n      = nrow(df_fit),
    R2     = round(r2, 3),
    Slope  = round(stats::coef(m)[["log_x"]], 3),
    Slope_p = round(sm$coefficients["log_x", "Pr(>|t|)"], 4),
    Spearman_rho = round(sp$estimate[[1]], 3),
    Spearman_p   = round(sp$p.value, 4),
    SW_p   = round(sw_p, 3),
    BP_p   = round(bp_p, 3),
    DW_stat = round(dw$statistic[[1]], 3),
    DW_p   = round(dw$p.value, 3)
  )
}

# Build joined df per site with both count and biomass FAIRD metrics
build_metric_df <- function(site_id) {
  faird_site <- etraps_raw %>%
    dplyr::filter(Device_type == "FAIRD", Site == site_id) %>%
    dplyr::group_by(Date, Site) %>%
    dplyr::summarise(
      faird_count    = sum(abundance, na.rm = TRUE),
      faird_biomass  = mean(biomass_mg, na.rm = TRUE),
      .groups = "drop"
    )

  ammod_site <- ammod_raw %>%
    dplyr::filter(Site == site_id) %>%
    dplyr::select(Date, Site, Live_mass, Dry_mass_g)

  dplyr::inner_join(faird_site, ammod_site, by = c("Date", "Site")) %>%
    dplyr::arrange(Date)
}

table14_rows <- list()

for (site_id in c("Site1", "Site2", "Site3")) {
  df_m <- build_metric_df(site_id)
  if (nrow(df_m) == 0) next

  # Analysis 1: FAIRD counts vs AMMOD dry mass
  row_counts <- run_metric_ols(df_m, "Dry_mass_g", "faird_count", site_id, "Counts vs Dry mass")

  # Analysis 2: FAIRD biomass vs AMMOD live mass (standard analysis)
  row_biomass <- run_metric_ols(df_m, "Live_mass", "faird_biomass", site_id, "Biomass vs Live mass")

  table14_rows[[length(table14_rows) + 1]] <- row_counts
  table14_rows[[length(table14_rows) + 1]] <- row_biomass

  cat(sprintf("%s: counts R²=%.3f | biomass R²=%.3f  ΔR²=%.3f (biomass advantage)\n",
              site_id,
              if (!is.null(row_counts))  row_counts$R2  else NA,
              if (!is.null(row_biomass)) row_biomass$R2 else NA,
              if (!is.null(row_counts) && !is.null(row_biomass))
                row_biomass$R2 - row_counts$R2 else NA))
}

table14_raw <- dplyr::bind_rows(table14_rows)

# Add ΔR² column (biomass minus counts, positive = biomass better)
table14_wide <- table14_raw %>%
  dplyr::select(Site, Metric, R2) %>%
  tidyr::pivot_wider(names_from = Metric, values_from = R2) %>%
  dplyr::mutate(
    Delta_R2 = `Biomass vs Live mass` - `Counts vs Dry mass`
  )

cat("\nTable 14 — Detection Metric Sensitivity:\n")
print(table14_wide)

readr::write_csv(table14_raw,  "outputs/tables/Table14_detection_metric_full.csv")
readr::write_csv(table14_wide, "outputs/tables/Table14_detection_metric_summary.csv")
cat("\nSaved: Table14_detection_metric_full.csv + Table14_detection_metric_summary.csv\n\n")


# ==============================================================================
# SECTION 6: BODY SIZE DISTRIBUTIONS — §3.3.1.6 (Table 15, supplementary figs)
# Source: H2_H3_extra_figures.R
# FAIRD-only: AMMOD has no body length (DNA metabarcoding).
# Size classes: <3mm (absolute limit), 3–5mm (constrained), >5mm (optimal).
# ==============================================================================

cat("==============================================================================\n")
cat("SECTION 6: BODY SIZE DISTRIBUTIONS — §3.3.1.6\n")
cat("==============================================================================\n\n")

etraps_indiv <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD", !is.na(Body_length))

# Per-site size class counts
size_table_rows <- list()

for (site_id in c("Site1", "Site2", "Site3", "Site4")) {
  df_site_sz <- etraps_indiv %>% dplyr::filter(Site == site_id)
  if (nrow(df_site_sz) == 0) next

  n_lt3 <- sum(df_site_sz$Body_length < 3)
  n_3to5 <- sum(df_site_sz$Body_length >= 3 & df_site_sz$Body_length < 5)
  n_gt5  <- sum(df_site_sz$Body_length >= 5)
  n_total <- n_lt3 + n_3to5 + n_gt5

  size_table_rows[[site_id]] <- tibble::tibble(
    Site           = site_id,
    n_total        = n_total,
    n_lt3mm        = n_lt3,
    pct_lt3mm      = round(100 * n_lt3  / n_total, 1),
    n_3to5mm       = n_3to5,
    pct_3to5mm     = round(100 * n_3to5 / n_total, 1),
    n_gt5mm        = n_gt5,
    pct_gt5mm      = round(100 * n_gt5  / n_total, 1)
  )
}

# Overall row
n_lt3_all  <- sum(etraps_indiv$Body_length < 3)
n_3to5_all <- sum(etraps_indiv$Body_length >= 3 & etraps_indiv$Body_length < 5)
n_gt5_all  <- sum(etraps_indiv$Body_length >= 5)
n_tot_all  <- n_lt3_all + n_3to5_all + n_gt5_all

size_table_rows[["Overall"]] <- tibble::tibble(
  Site       = "Overall",
  n_total    = n_tot_all,
  n_lt3mm    = n_lt3_all,
  pct_lt3mm  = round(100 * n_lt3_all  / n_tot_all, 1),
  n_3to5mm   = n_3to5_all,
  pct_3to5mm = round(100 * n_3to5_all / n_tot_all, 1),
  n_gt5mm    = n_gt5_all,
  pct_gt5mm  = round(100 * n_gt5_all  / n_tot_all, 1)
)

table15 <- dplyr::bind_rows(size_table_rows)
cat("Table 15 — FAIRD Body Size Distribution by Site:\n")
print(table15)

readr::write_csv(table15, "outputs/tables/Table15_body_size_distributions.csv")
cat("\nSaved: Table15_body_size_distributions.csv\n\n")

# Supplementary figure — overall body length histogram with thresholds
p_size_hist <- ggplot2::ggplot(etraps_indiv, ggplot2::aes(x = Body_length)) +
  ggplot2::geom_histogram(bins = 40, fill = custom_colors[["FAIRD1"]],
                          color = "white", alpha = 0.85) +
  ggplot2::annotate("rect", xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf,
                    fill = "#E74C3C", alpha = 0.08) +
  ggplot2::annotate("rect", xmin = 3, xmax = 5, ymin = -Inf, ymax = Inf,
                    fill = "#F39C12", alpha = 0.08) +
  ggplot2::geom_vline(xintercept = 3, linetype = "dashed",
                      color = "#E74C3C", linewidth = 1) +
  ggplot2::geom_vline(xintercept = 5, linetype = "dotted",
                      color = "#F39C12", linewidth = 1) +
  ggplot2::annotate("text", x = 3, y = Inf,
                    label = sprintf("Absolute limit 3mm\n(<3mm: %d, %.1f%%)",
                                     n_lt3_all, 100 * n_lt3_all / n_tot_all),
                    hjust = -0.05, vjust = 1.2, size = 3, color = "#E74C3C") +
  ggplot2::annotate("text", x = 5, y = Inf,
                    label = sprintf("Constrained <5mm\n(3–5mm: %d, %.1f%%)",
                                     n_3to5_all, 100 * n_3to5_all / n_tot_all),
                    hjust = -0.05, vjust = 1.2, size = 3, color = "#F39C12") +
  ggplot2::labs(
    title    = "Supplementary — FAIRD Body Length Distribution",
    subtitle = sprintf("n = %d individuals | AMMOD has no body length (DNA metabarcoding)",
                        n_tot_all),
    x        = "Body Length (mm)",
    y        = "Number of Individuals"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title        = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle     = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor  = ggplot2::element_blank()
  )

ggplot2::ggsave("outputs/figures/FigS_H2_BodySize_Distribution.png", p_size_hist,
                width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: FigS_H2_BodySize_Distribution.png\n\n")

# Per-site faceted version
p_size_facet <- ggplot2::ggplot(
  etraps_indiv %>%
    dplyr::mutate(Size_class = dplyr::case_when(
      Body_length < 3              ~ "<3mm",
      Body_length >= 3 & Body_length < 5 ~ "3–5mm",
      Body_length >= 5             ~ "≥5mm"
    )),
  ggplot2::aes(x = Body_length, fill = Size_class)
) +
  ggplot2::geom_histogram(bins = 30, color = "white", alpha = 0.85) +
  ggplot2::scale_fill_manual(values = c("<3mm" = "#E74C3C",
                                         "3–5mm" = "#F39C12",
                                         "≥5mm"  = custom_colors[["FAIRD1"]])) +
  ggplot2::facet_wrap(~ Site, scales = "free_y") +
  ggplot2::geom_vline(xintercept = c(3, 5), linetype = "dashed",
                      color = "grey30", linewidth = 0.5) +
  ggplot2::labs(
    title = "Supplementary — FAIRD Body Length by Site",
    x     = "Body Length (mm)",
    y     = "Number of Individuals",
    fill  = "Size class"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top"
  )

ggplot2::ggsave("outputs/figures/FigS_H2_BodySize_BySite.png", p_size_facet,
                width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: FigS_H2_BodySize_BySite.png\n\n")


# ==============================================================================
# SECTION 7: TABLE 13 ASSEMBLY — all H2 quantitative results
# ==============================================================================

cat("==============================================================================\n")
cat("SECTION 7: TABLE 13 ASSEMBLY\n")
cat("==============================================================================\n\n")

# Helper to build a Table 13 row from a result list
make_t13_row <- function(level, site_or_habitat, n, r2, adj_r2,
                          slope, slope_ci_low, slope_ci_high, slope_p,
                          intercept, intercept_p,
                          sw_p, bp_p, dw_stat, dw_p,
                          pseudo_r2, n_low_wt,
                          concordance_pct,
                          spearman_rho, spearman_p,
                          kendall_tau, kendall_p,
                          boot_ci_low = NA, boot_ci_high = NA,
                          ccf_lag0 = NA, ccf_ci = NA, ccf_sig = NA) {
  tibble::tibble(
    Level          = level,
    Site_Habitat   = site_or_habitat,
    n              = n,
    R2             = round(r2, 3),
    Adj_R2         = round(adj_r2, 3),
    Boot_CI_low    = if (!is.na(boot_ci_low))  round(boot_ci_low,  3) else NA,
    Boot_CI_high   = if (!is.na(boot_ci_high)) round(boot_ci_high, 3) else NA,
    Slope          = round(slope, 3),
    Slope_CI_low   = round(slope_ci_low, 3),
    Slope_CI_high  = round(slope_ci_high, 3),
    Slope_p        = round(slope_p, 4),
    Sig            = sig_stars(slope_p),
    Intercept      = round(intercept, 3),
    Intercept_p    = round(intercept_p, 4),
    Pseudo_R2      = round(pseudo_r2, 3),
    n_low_wt       = n_low_wt,
    Concordance_pct = round(concordance_pct, 1),
    Spearman_rho   = round(spearman_rho, 3),
    Spearman_p     = round(spearman_p, 4),
    Kendall_tau    = round(kendall_tau, 3),
    Kendall_p      = round(kendall_p, 4),
    SW_p           = round(sw_p, 3),
    BP_p           = round(bp_p, 3),
    DW_stat        = round(dw_stat, 3),
    DW_p           = round(dw_p, 3),
    CCF_lag0       = if (!is.na(ccf_lag0)) round(ccf_lag0, 3) else NA,
    CCF_CI_bound   = if (!is.na(ccf_ci))   round(ccf_ci,   3) else NA,
    CCF_sig        = ccf_sig
  )
}

table13_raw <- dplyr::bind_rows(
  # Site Level
  make_t13_row("Site Level", "Site1 (Maize)", res_site1$n,
               res_site1$r_squared, res_site1$adj_r_squared,
               res_site1$slope, res_site1$slope_ci_low, res_site1$slope_ci_high,
               res_site1$slope_p, res_site1$intercept, res_site1$intercept_p,
               res_site1$sw_p, res_site1$bp_p, res_site1$dw_stat, res_site1$dw_p,
               res_site1$pseudo_r2, res_site1$n_low_wt,
               res_site1$concordance_pct,
               res_site1$spearman_rho, res_site1$spearman_p,
               res_site1$kendall_tau, res_site1$kendall_p,
               ccf_lag0 = ccf_s1$lag0_r, ccf_ci = ccf_s1$ci_bound,
               ccf_sig  = ifelse(ccf_s1$significant, "sig", "ns")),

  make_t13_row("Site Level", "Site2 (Maize)", res_site2$n,
               res_site2$r_squared, res_site2$adj_r_squared,
               res_site2$slope, res_site2$slope_ci_low, res_site2$slope_ci_high,
               res_site2$slope_p, res_site2$intercept, res_site2$intercept_p,
               res_site2$sw_p, res_site2$bp_p, res_site2$dw_stat, res_site2$dw_p,
               res_site2$pseudo_r2, res_site2$n_low_wt,
               res_site2$concordance_pct,
               res_site2$spearman_rho, res_site2$spearman_p,
               res_site2$kendall_tau, res_site2$kendall_p,
               ccf_lag0 = ccf_s2$lag0_r, ccf_ci = ccf_s2$ci_bound,
               ccf_sig  = ifelse(ccf_s2$significant, "sig", "ns")),

  make_t13_row("Site Level", "Site3 (Meadow)", res_site3$n,
               res_site3$r_squared, res_site3$adj_r_squared,
               res_site3$slope, res_site3$slope_ci_low, res_site3$slope_ci_high,
               res_site3$slope_p, res_site3$intercept, res_site3$intercept_p,
               res_site3$sw_p, res_site3$bp_p, res_site3$dw_stat, res_site3$dw_p,
               res_site3$pseudo_r2, res_site3$n_low_wt,
               res_site3$concordance_pct,
               res_site3$spearman_rho, res_site3$spearman_p,
               res_site3$kendall_tau, res_site3$kendall_p,
               ccf_lag0 = ccf_s3$lag0_r, ccf_ci = ccf_s3$ci_bound,
               ccf_sig  = ifelse(ccf_s3$significant, "sig", "ns")),

  # Intra-Habitat Level — Maize
  make_t13_row("Intra-Habitat Level", "Maize (Sites 1+2)", nrow(df_maize),
               r2_maize, summary(model_maize)$adj.r.squared,
               slope_maize, slope_ci_m[1], slope_ci_m[2],
               slope_p_m, intercept_m, intercept_p_m,
               sw_maize, bp_maize, dw_stat_m, dw_p_m,
               pseudo_r2_m, n_low_wt_m,
               conc_m,
               sp_m$estimate[[1]], sp_m$p.value,
               kt_m$estimate[[1]], kt_m$p.value,
               boot_ci_low = boot_ci_m[1], boot_ci_high = boot_ci_m[2],
               ccf_lag0 = ccf_maize_summed$lag0_r,
               ccf_ci   = ccf_maize_summed$ci_bound,
               ccf_sig  = ifelse(ccf_maize_summed$significant, "sig", "ns")),

  # Overall Level
  make_t13_row("Overall Level", "All Sites (1+2+3)", nrow(df_overall),
               r2_overall, summary(model_overall)$adj.r.squared,
               slope_ov, slope_ci_ov[1], slope_ci_ov[2],
               slope_p_ov, intercept_ov, intercept_p_ov,
               sw_ov, bp_ov, dw_stat_ov, dw_p_ov,
               pseudo_r2_ov, n_low_wt_ov,
               conc_ov,
               sp_ov$estimate[[1]], sp_ov$p.value,
               kt_ov$estimate[[1]], kt_ov$p.value,
               boot_ci_low = boot_ci_ov[1], boot_ci_high = boot_ci_ov[2],
               ccf_lag0 = ccf_overall_desc$lag0_r,
               ccf_ci   = ccf_overall_desc$ci_bound,
               ccf_sig  = paste0(ifelse(ccf_overall_desc$significant, "sig", "ns"),
                                  " [descriptive only — Site3 n=8]"))
)

cat("Table 13 — H2 Quantitative Validation Summary:\n")
print(table13_raw, width = Inf)

readr::write_csv(table13_raw, "outputs/tables/Table13_H2_quanti_validation.csv")
cat("\nSaved: Table13_H2_quanti_validation.csv\n\n")

# Display-formatted Table 13
table13_display <- table13_raw %>%
  dplyr::mutate(
    R2_display     = sprintf("%.3f%s", R2, Sig),
    Boot_CI        = dplyr::if_else(!is.na(Boot_CI_low),
                                    sprintf("[%.3f–%.3f]", Boot_CI_low, Boot_CI_high),
                                    "—"),
    Slope_display  = sprintf("%.3f [%.3f, %.3f]", Slope, Slope_CI_low, Slope_CI_high),
    Spearman_disp  = sprintf("%.3f (%s)", Spearman_rho, sig_stars(Spearman_p)),
    Kendall_disp   = sprintf("%.3f (%s)", Kendall_tau,  sig_stars(Kendall_p)),
    CCF_display    = dplyr::if_else(!is.na(CCF_lag0),
                                    sprintf("%.3f [±%.3f] %s", CCF_lag0, CCF_CI_bound,
                                            ifelse(CCF_sig == "sig", "*", "ns")),
                                    "—")
  ) %>%
  dplyr::select(Level, Site_Habitat, n,
                R2_display, Boot_CI,
                Slope_display,
                Spearman_disp, Kendall_disp,
                Concordance_pct, Pseudo_R2,
                CCF_display,
                SW_p, BP_p, DW_stat)

cat("\nTable 13 — Display Format:\n")
print(table13_display, width = Inf)

readr::write_csv(table13_display, "outputs/tables/Table13_H2_quanti_validation_display.csv")
cat("Saved: Table13_H2_quanti_validation_display.csv\n\n")

# ==============================================================================
# END
# ==============================================================================

cat("==============================================================================\n")
cat("H2_quanti_validation.R COMPLETE\n")
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("Outputs:\n")
cat("  Figures 16–24: outputs/figures/Fig16–Fig24_*.png\n")
cat("  Supplementary: FigS_H2_BodySize_*.png\n")
cat("  Tables 13–15: outputs/tables/Table13–15_*.csv\n")
cat("==============================================================================\n")

while (sink.number() > 0) sink()
