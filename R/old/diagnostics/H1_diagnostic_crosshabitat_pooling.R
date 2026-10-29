# ==============================================================================
# H1_diagnostic_crosshabitat_pooling.R
# Diagnostic: H1 Cross-Habitat Level using pooling (vs current averaging approach)
#
# Current approach (H1_consistency_FAIRD_Maize_Meadow.R):
#   Average(FAIRD1, FAIRD2) vs Average(FAIRD3, FAIRD4) → n = 22 days
#
# This approach (pooling):
#   FAIRD1↔FAIRD3 and FAIRD2↔FAIRD4 stacked → n = 42 pairs
#   Analogous to H2 pooling (FAIRD1↔AMMOD1, FAIRD2↔AMMOD2 → n=24)
#   Site-pair covariate validates that pooling is appropriate.
#
# Purpose: Exploratory comparison of both approaches before deciding
#          whether to update H1_consistency_FAIRD_Maize_Meadow.R
# ==============================================================================

# ============================================================
# CONSOLE OUTPUT CAPTURE
# ============================================================
script_name <- sub("\\.R$", "", basename(sys.frame(1)$ofile))

output_file <- file.path("outputs",
                         paste0(script_name, "_output_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

while (sink.number() > 0) sink()
sink(output_file, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Script:", paste0(script_name, ".R"), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ============================================================
# LIBRARIES AND DATA
# ============================================================
library(tidyverse)
library(patchwork)
library(car)
library(lmtest)

source("load_data.R")
load("data_cleaned.RData")
source("custom_colors.R")

cat("[OK] Colors loaded: FAIRD =", custom_colors["FAIRD"], "\n\n")

etraps_daily <- data_list$etraps_daily

# ============================================================
# SECTION 1: CURRENT APPROACH — AVERAGING
# ============================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SECTION 1: CURRENT APPROACH — AVERAGING (n = 22 days)\n")
cat("  Average(FAIRD1+FAIRD2)/2 vs Average(FAIRD3+FAIRD4)/2\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

faird_avg <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2", "FAIRD3", "FAIRD4")) %>%
  dplyr::mutate(
    Habitat = dplyr::case_when(
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
    names_from  = Habitat,
    values_from = mean_abundance,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    log_Maize  = log10(Maize + 1),
    log_Meadow = log10(Meadow + 1)
  )

fit_avg     <- lm(log_Meadow ~ log_Maize, data = faird_avg)
sum_avg     <- summary(fit_avg)
spear_avg   <- cor.test(faird_avg$log_Maize, faird_avg$log_Meadow,
                        method = "spearman", exact = FALSE)
pval_avg    <- sum_avg$coefficients["log_Maize", "Pr(>|t|)"]
sig_avg     <- dplyr::case_when(
  pval_avg < 0.001 ~ "***", pval_avg < 0.01 ~ "**",
  pval_avg < 0.05  ~ "*",   TRUE             ~ "ns"
)

# Kendall tau on daily changes
changes_avg <- faird_avg %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    Maize_change  = Maize  - dplyr::lag(Maize),
    Meadow_change = Meadow - dplyr::lag(Meadow),
    concordant    = sign(Maize_change) == sign(Meadow_change)
  ) %>%
  dplyr::filter(!is.na(Maize_change))

concordance_avg <- 100 * sum(changes_avg$concordant) / nrow(changes_avg)
kt_avg          <- cor.test(changes_avg$Maize_change, changes_avg$Meadow_change,
                            method = "kendall", exact = FALSE)

cat(sprintf("n = %d days\n", nrow(faird_avg)))
cat(sprintf("R² = %.3f %s | Adj.R² = %.3f\n",
            sum_avg$r.squared, sig_avg, sum_avg$adj.r.squared))
cat(sprintf("Slope = %.3f [%.3f, %.3f]\n",
            coef(fit_avg)[2], confint(fit_avg)[2,1], confint(fit_avg)[2,2]))
cat(sprintf("Spearman ρ = %.3f (p = %.4f)\n",
            spear_avg$estimate, spear_avg$p.value))
cat(sprintf("Kendall τ (daily changes) = %.3f (p = %.4f)\n",
            kt_avg$estimate, kt_avg$p.value))
cat(sprintf("Concordance = %.1f%%\n", concordance_avg))

# ============================================================
# SECTION 2: NEW APPROACH — POOLING
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SECTION 2: NEW APPROACH — POOLING (n = 42 pairs)\n")
cat("  FAIRD1↔FAIRD3 (n=21) + FAIRD2↔FAIRD4 (n=21)\n")
cat("  Analogous to H2 pooling: FAIRD1↔AMMOD1 + FAIRD2↔AMMOD2\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Build wide format per device
faird_wide <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2", "FAIRD3", "FAIRD4")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

# Create natural cross-habitat pairs (site-matched)
# Pair A: FAIRD1 (Maize Site1) ↔ FAIRD3 (Meadow Site3)
# Pair B: FAIRD2 (Maize Site2) ↔ FAIRD4 (Meadow Site4)
faird_pooled <- dplyr::bind_rows(
  faird_wide %>%
    dplyr::select(Date, Maize = FAIRD1, Meadow = FAIRD3) %>%
    dplyr::mutate(Site_pair = "FAIRD1-FAIRD3"),
  faird_wide %>%
    dplyr::select(Date, Maize = FAIRD2, Meadow = FAIRD4) %>%
    dplyr::mutate(Site_pair = "FAIRD2-FAIRD4")
) %>%
  dplyr::mutate(
    log_Maize  = log10(Maize + 1),
    log_Meadow = log10(Meadow + 1)
  )

cat(sprintf("Pooled dataset: n = %d pairs\n", nrow(faird_pooled)))
cat("n per site-pair:\n")
print(table(faird_pooled$Site_pair))
cat("\n")

# --- Simple OLS (without covariate) — comparable to averaging approach ---
fit_pool    <- lm(log_Meadow ~ log_Maize, data = faird_pooled)
sum_pool    <- summary(fit_pool)
spear_pool  <- cor.test(faird_pooled$log_Maize, faird_pooled$log_Meadow,
                        method = "spearman", exact = FALSE)
pval_pool   <- sum_pool$coefficients["log_Maize", "Pr(>|t|)"]
sig_pool    <- dplyr::case_when(
  pval_pool < 0.001 ~ "***", pval_pool < 0.01 ~ "**",
  pval_pool < 0.05  ~ "*",   TRUE              ~ "ns"
)

# Diagnostics
shapiro_pool <- shapiro.test(residuals(fit_pool))$p.value
bp_pool      <- car::ncvTest(fit_pool)$p
dw_pool      <- lmtest::dwtest(fit_pool)

cat("--- Simple OLS (no covariate) ---\n")
cat(sprintf("R² = %.3f %s | Adj.R² = %.3f\n",
            sum_pool$r.squared, sig_pool, sum_pool$adj.r.squared))
cat(sprintf("Slope = %.3f [%.3f, %.3f]\n",
            coef(fit_pool)[2], confint(fit_pool)[2,1], confint(fit_pool)[2,2]))
cat(sprintf("Spearman ρ = %.3f (p = %.4f)\n",
            spear_pool$estimate, spear_pool$p.value))
cat(sprintf("Diagnostics: Shapiro p=%.3f | BP p=%.3f | DW=%.3f (p=%.4f)\n",
            shapiro_pool, bp_pool,
            as.numeric(dw_pool$statistic), dw_pool$p.value))

# --- Site-pair covariate model — validates pooling ---
cat("\n--- Site-pair covariate model (pooling validation) ---\n")
fit_cov <- lm(log_Meadow ~ log_Maize + Site_pair, data = faird_pooled)
cat("Site_pair covariate:\n")
print(round(summary(fit_cov)$coefficients, 4))
cat("  → If Site_pair p > 0.05: pooling is valid\n\n")

# --- Kendall tau on daily changes (within site-pair, then pooled) ---
changes_pool <- faird_pooled %>%
  dplyr::arrange(Site_pair, Date) %>%
  dplyr::group_by(Site_pair) %>%
  dplyr::mutate(
    Maize_change  = Maize  - dplyr::lag(Maize),
    Meadow_change = Meadow - dplyr::lag(Meadow),
    concordant    = sign(Maize_change) == sign(Meadow_change)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(Maize_change))

concordance_pool <- 100 * sum(changes_pool$concordant) / nrow(changes_pool)
kt_pool          <- cor.test(changes_pool$Maize_change, changes_pool$Meadow_change,
                             method = "kendall", exact = FALSE)

cat(sprintf("Kendall τ (daily changes, n changes = %d) = %.3f (p = %.4f)\n",
            nrow(changes_pool), kt_pool$estimate, kt_pool$p.value))
cat(sprintf("Concordance = %.1f%%\n", concordance_pool))

# ============================================================
# SECTION 3: DIRECT COMPARISON
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SECTION 3: DIRECT COMPARISON — AVERAGING vs POOLING\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("%-30s %10s %10s\n", "Metric", "Averaging", "Pooling"))
cat(strrep("-", 52), "\n")
cat(sprintf("%-30s %10d %10d\n", "n", nrow(faird_avg), nrow(faird_pooled)))
cat(sprintf("%-30s %10s %10s\n", "R²",
            sprintf("%.3f %s", sum_avg$r.squared, sig_avg),
            sprintf("%.3f %s", sum_pool$r.squared, sig_pool)))
cat(sprintf("%-30s %10.3f %10.3f\n", "Adj.R²",
            sum_avg$adj.r.squared, sum_pool$adj.r.squared))
cat(sprintf("%-30s %10.3f %10.3f\n", "Slope",
            coef(fit_avg)[2], coef(fit_pool)[2]))
cat(sprintf("%-30s %10.3f %10.3f\n", "Spearman ρ",
            spear_avg$estimate, spear_pool$estimate))
cat(sprintf("%-30s %10.3f %10.3f\n", "Kendall τ (changes)",
            kt_avg$estimate, kt_pool$estimate))
cat(sprintf("%-30s %9.1f%% %9.1f%%\n", "Concordance",
            concordance_avg, concordance_pool))
cat(strrep("-", 52), "\n\n")

cat("INTERPRETATION:\n")
cat(sprintf("  ΔR²    = %.3f (pooling - averaging)\n",
            sum_pool$r.squared - sum_avg$r.squared))
cat(sprintf("  ΔAdj.R² = %.3f\n",
            sum_pool$adj.r.squared - sum_avg$adj.r.squared))
if (abs(sum_pool$r.squared - sum_avg$r.squared) < 0.05) {
  cat("  → Results are consistent between approaches (ΔR² < 0.05)\n")
  cat("  → Pooling preferred: larger n, no habitat-as-treatment assumption\n")
} else {
  cat("  → Notable difference between approaches — review carefully\n")
}

# ============================================================
# SECTION 4: FIGURE — side-by-side scatter comparison
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SECTION 4: FIGURES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

faird_color <- custom_colors["FAIRD"]

# Panel A: averaging
p_avg <- ggplot(faird_avg, aes(x = log_Maize, y = log_Meadow)) +
  geom_point(size = 3, alpha = 0.7, color = "gray30") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "gray30", fill = faird_color, alpha = 0.25, linewidth = 1) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("R\u00b2 = %.3f %s\nn = %d", sum_avg$r.squared, sig_avg, nrow(faird_avg)),
           hjust = -0.1, vjust = 1.3, size = 3.5, fontface = "italic") +
  labs(title = "A) Averaging approach",
       subtitle = "(FAIRD1+FAIRD2)/2 vs (FAIRD3+FAIRD4)/2",
       x = "log\u2081\u2080(FAIRD Maize + 1)",
       y = "log\u2081\u2080(FAIRD Meadow + 1)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40", size = 9))

# Panel B: pooling
pair_colors <- c("FAIRD1-FAIRD3" = custom_colors["Site1"],
                 "FAIRD2-FAIRD4" = custom_colors["Site2"])

p_pool <- ggplot(faird_pooled, aes(x = log_Maize, y = log_Meadow, color = Site_pair)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE,
              color = "gray30", fill = "gray70", alpha = 0.25, linewidth = 1) +
  scale_color_manual(values = pair_colors, name = "Site pair") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("R\u00b2 = %.3f %s\nn = %d", sum_pool$r.squared, sig_pool, nrow(faird_pooled)),
           hjust = -0.1, vjust = 1.3, size = 3.5, fontface = "italic") +
  labs(title = "B) Pooling approach",
       subtitle = "FAIRD1\u2194FAIRD3 + FAIRD2\u2194FAIRD4",
       x = "log\u2081\u2080(FAIRD Maize + 1)",
       y = "log\u2081\u2080(FAIRD Meadow + 1)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40", size = 9),
        legend.position = "right")

fig_compare <- p_avg + p_pool +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(
    title    = "H1 Cross-Habitat Level: Averaging vs Pooling",
    subtitle = "FAIRD Maize vs Meadow temporal consistency",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

print(fig_compare)

ggsave("outputs/H1_diagnostic_crosshabitat_pooling.png",
       fig_compare, width = 14, height = 6, dpi = 300, bg = "white")
ggsave("outputs/H1_diagnostic_crosshabitat_pooling.pdf",
       fig_compare, width = 14, height = 6, device = cairo_pdf)

cat("[OK] Comparison figure saved to outputs/\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
while (sink.number() > 0) sink()
cat(">> Output saved to:", output_file, "\n")