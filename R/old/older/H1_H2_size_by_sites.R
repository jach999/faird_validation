# ==============================================================================
# EXPLORATORY ANALYSIS: BODY SIZE DISTRIBUTION BY SITE (AMMOD OTU DATA)
# ==============================================================================
# Purpose: Compare daily mean body length of detected taxa across sites
#          to explore whether size-structure differs between maize and meadow
#
# Data source: AMMOD_24h_otu_counts_For_Statistics_trimmed.csv
# Biometric source: Söhlström et al. (2018) — assigned at Family level
#
# Design:
#   - Deduplicate to Family level per site per day (biometrics = family-level)
#   - Daily mean Body_length per site = mean of unique families detected that day
#   - Pairwise Mann-Whitney tests across all 6 site pairs (common days only)
#   - Bonferroni correction for 6 comparisons
#   - Time series plot (4 sites) + violin plot (Maize vs Meadow)
#
# Interpretation scope:
#   "The daily mean body size of taxa detected in meadow sites differed
#    systematically from that in maize sites" — NOT a dominance claim
#
# Author: Juan Chiavassa
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

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

# ==============================================================================
# SECTION 1: LOAD AND PREPARE DATA
# ==============================================================================

cat("=== SECTION 1: DATA PREPARATION ===\n\n")

ammod_otu <- readr::read_csv(
  "ammod_taxonomy_presence_clean.csv",
  show_col_types = FALSE
) %>%
  dplyr::mutate(Date = as.Date(Date))

cat(sprintf("Raw OTU rows loaded: %d\n", nrow(ammod_otu)))
cat(sprintf("Sites present: %s\n", paste(sort(unique(ammod_otu$Site)), collapse = ", ")))
cat(sprintf("Date range: %s to %s\n\n",
            min(ammod_otu$Date), max(ammod_otu$Date)))

# ------------------------------------------------------------------------------
# Deduplicate to Family level per site per day
# Rationale: biometrics from Söhlström et al. (2018) assigned at Family level,
# so multiple species within the same family on the same day carry identical
# Body_length values and must not be treated as independent observations.
# ------------------------------------------------------------------------------

daily_size <- ammod_otu %>%
  dplyr::filter(!is.na(Body_length), !is.na(Family)) %>%
  dplyr::distinct(Date, Site, Family, .keep_all = TRUE) %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    mean_body_length = mean(Body_length, na.rm = TRUE),
    n_families       = dplyr::n(),
    .groups          = "drop"
  ) %>%
  dplyr::mutate(
    Habitat = dplyr::case_when(
      Site %in% c("Site1", "Site2") ~ "Maize",
      Site %in% c("Site3", "Site4") ~ "Meadow"
    )
  )

cat("Daily mean body length (Family-deduplicated) per site:\n")
daily_size %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n_days         = dplyr::n(),
    mean_BL        = round(mean(mean_body_length), 2),
    sd_BL          = round(sd(mean_body_length), 2),
    median_BL      = round(median(mean_body_length), 2),
    mean_n_families = round(mean(n_families), 1),
    .groups = "drop"
  ) %>%
  print()

cat("\n")

# ==============================================================================
# SECTION 2: PAIRWISE MANN-WHITNEY TESTS (common days per pair)
# ==============================================================================

cat("=== SECTION 2: PAIRWISE MANN-WHITNEY TESTS ===\n\n")
cat("Note: Each comparison restricted to days available in BOTH sites\n")
cat("      Bonferroni correction applied for 6 comparisons (alpha = 0.05/6)\n\n")

site_pairs <- list(
  c("Site1", "Site2"),  # intra-maize
  c("Site3", "Site4"),  # intra-meadow
  c("Site1", "Site3"),  # cross-habitat
  c("Site1", "Site4"),  # cross-habitat
  c("Site2", "Site3"),  # cross-habitat
  c("Site2", "Site4")   # cross-habitat
)

bonferroni_alpha <- 0.05 / length(site_pairs)

mw_results <- purrr::map_dfr(site_pairs, function(pair) {
  s1 <- pair[1]
  s2 <- pair[2]
  
  # Restrict to common days
  days_s1 <- daily_size %>% dplyr::filter(Site == s1) %>% dplyr::pull(Date)
  days_s2 <- daily_size %>% dplyr::filter(Site == s2) %>% dplyr::pull(Date)
  common_days <- intersect(days_s1, days_s2)
  
  vals_s1 <- daily_size %>%
    dplyr::filter(Site == s1, Date %in% common_days) %>%
    dplyr::pull(mean_body_length)
  
  vals_s2 <- daily_size %>%
    dplyr::filter(Site == s2, Date %in% common_days) %>%
    dplyr::pull(mean_body_length)
  
  if (length(vals_s1) < 3 || length(vals_s2) < 3) {
    cat(sprintf("  %s vs %s: Insufficient common days (n=%d) — skipped\n",
                s1, s2, length(common_days)))
    return(NULL)
  }
  
  wt <- stats::wilcox.test(vals_s1, vals_s2, exact = FALSE)
  
  hab_s1 <- ifelse(s1 %in% c("Site1", "Site2"), "Maize", "Meadow")
  hab_s2 <- ifelse(s2 %in% c("Site1", "Site2"), "Maize", "Meadow")
  comparison_type <- ifelse(hab_s1 == hab_s2, "Intra-habitat", "Cross-habitat")
  
  tibble::tibble(
    Pair             = paste(s1, "vs", s2),
    Type             = comparison_type,
    n_common_days    = length(common_days),
    median_s1        = round(median(vals_s1), 2),
    median_s2        = round(median(vals_s2), 2),
    W_statistic      = wt$statistic,
    p_value          = wt$p.value,
    p_bonferroni     = min(wt$p.value * length(site_pairs), 1),
    significant_bonf = wt$p.value < bonferroni_alpha
  )
})

cat("Results:\n")
print(mw_results, n = Inf)
cat(sprintf("\nBonferroni-adjusted alpha: %.4f\n\n", bonferroni_alpha))

readr::write_csv(mw_results,
                 "outputs/BodySize_PairwiseMannWhitney.csv")
cat("[OK] outputs/BodySize_PairwiseMannWhitney.csv\n\n")

# ==============================================================================
# SECTION 3: TIME SERIES PLOT (4 sites)
# ==============================================================================

cat("=== SECTION 3: TIME SERIES PLOT ===\n\n")

plot_ts <- ggplot(
  daily_size,
  aes(x = Date, y = mean_body_length,
      color = Site, shape = Site, linetype = Habitat)
) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(
    "Site1" = 16, "Site2" = 17, "Site3" = 15, "Site4" = 18
  )) +
  scale_linetype_manual(values = c("Maize" = "solid", "Meadow" = "dashed")) +
  labs(
    title    = "Daily Mean Body Length of Detected Taxa by Site",
    subtitle = "Family-level deduplication | Biometrics: Söhlström et al. (2018)",
    x        = "Date",
    y        = "Mean Body Length (mm)",
    color    = "Site",
    shape    = "Site",
    linetype = "Habitat"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "gray40"),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(plot_ts)

ggsave("outputs/BodySize_TimeSeries.png",
       plot = plot_ts, width = 12, height = 5, dpi = 300, bg = "white")
ggsave("outputs/BodySize_TimeSeries.pdf",
       plot = plot_ts, width = 12, height = 5, device = cairo_pdf)

cat("[OK] outputs/BodySize_TimeSeries.png/.pdf\n\n")

# ==============================================================================
# SECTION 4: VIOLIN PLOT (Maize vs Meadow, all daily values)
# ==============================================================================

cat("=== SECTION 4: VIOLIN PLOT ===\n\n")

plot_violin <- ggplot(
  daily_size,
  aes(x = Site, y = mean_body_length, fill = Habitat)
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = 16,
               outlier.size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Maize" = unname(custom_colors["Site1"]),
                               "Meadow" = unname(custom_colors["Site3"]))) +
  labs(
    title = "Distribution of Daily Mean Body Length by Site",
    x     = "Site",
    y     = "Mean Body Length (mm)",
    fill  = "Habitat"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

print(plot_violin)

ggsave("outputs/BodySize_Violin.png",
       plot = plot_violin, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("outputs/BodySize_Violin.pdf",
       plot = plot_violin, width = 8, height = 5, device = cairo_pdf)

cat("[OK] outputs/BodySize_Violin.png/.pdf\n\n")

# ==============================================================================
# COMBINED FIGURE
# ==============================================================================

combined_size <- plot_ts + plot_violin +
  patchwork::plot_layout(ncol = 2, widths = c(1.6, 1)) +
  patchwork::plot_annotation(
    title = "Body Size Distribution Across Sites (AMMOD OTU Data)",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

print(combined_size)

ggsave("outputs/BodySize_Combined.png",
       plot = combined_size, width = 16, height = 5, dpi = 300, bg = "white")
ggsave("outputs/BodySize_Combined.pdf",
       plot = combined_size, width = 16, height = 5, device = cairo_pdf)

cat("[OK] outputs/BodySize_Combined.png/.pdf\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("=== ANALYSIS COMPLETE ===\n\n")
cat("IMPORTANT INTERPRETATION NOTE:\n")
cat("  Biometrics assigned at Family level (Söhlström et al. 2018)\n")
cat("  => Daily values reflect mean size of DETECTED FAMILY POOL, not\n")
cat("     individual insect abundance or dominance.\n")
cat("  => Results should be interpreted as exploratory and hypothesis-generating.\n\n")

cat("Files generated:\n")
cat("  [OK] outputs/BodySize_PairwiseMannWhitney.csv\n")
cat("  [OK] outputs/BodySize_TimeSeries.png/.pdf\n")
cat("  [OK] outputs/BodySize_Violin.png/.pdf\n")
cat("  [OK] outputs/BodySize_Combined.png/.pdf\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")