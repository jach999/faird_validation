# ==============================================================================
# H2_quali_validation.R  — §3.3.2 Inter-System Qualitative Agreement
# FAIRD vs AMMOD — Izsak-Price ΔS — Superfamily — 3-day non-overlapping windows
# ==============================================================================

library(tidyverse)
library(here)
library(vegan)
library(effsize)
library(car)
library(patchwork)
library(scales)
library(MASS)

here::i_am("R/H2_quali_validation.R")

dir.create("outputs/figures",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",   showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console",  showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/H2_quali_validation_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("H2 QUALITATIVE VALIDATION — §3.3.2\n")
cat("FAIRD vs AMMOD | Izsak-Price ΔS | Superfamily | 3-day windows\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source("R/00_load_data.R")
source("R/00_custom_colors.R")


# ==============================================================================
# SECTION 0: SHARED TEMPORAL SETUP
# ==============================================================================

date_start           <- as.Date("2023-08-23")
date_end             <- as.Date("2023-09-13")
date_exclude         <- as.Date("2023-09-04")   # bottle change day
ammod4_failure_start <- as.Date("2023-08-28")
ammod4_failure_end   <- as.Date("2023-09-03")
nc                   <- "#N/C"
MIN_INDIV            <- 3
SIZE_3MM             <- 3
SIZE_5MM             <- 5

sites     <- c("Site1", "Site2", "Site3", "Site4")  # Site4 included: OTUs available (no biomass)
all_dates <- seq.Date(date_start, date_end, by = "day")

window_breaks <- seq(date_start, date_end + 2, by = 3)
window_lookup <- tibble::tibble(Date = all_dates) %>%
  dplyr::mutate(
    Window_ID    = as.integer(cut(Date, breaks = window_breaks,
                                  right = FALSE, labels = FALSE)),
    Window_Start = window_breaks[Window_ID],
    Window_End   = pmin(Window_Start + 2, date_end),
    Window_Label = paste0("W", Window_ID, " (",
                          format(Window_Start, "%b %d"), "-",
                          format(Window_End,   "%b %d"), ")")
  )


# ==============================================================================
# SECTION 1: SITE LEVEL (§3.3.2.1)
# ==============================================================================

cat("=== SECTION 1: SITE LEVEL ===\n\n")

# --- 1a. FAIRD data (data_list$faird is already FAIRD-only) ---
faird_data <- data_list$faird %>%
  dplyr::filter(Date >= date_start & Date <= date_end,
                Site %in% sites) %>%
  dplyr::mutate(
    Taxon_SF = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc & Superfamily != "",
      Superfamily, NA_character_
    )
  )

cat("FAIRD individuals:", nrow(faird_data), "\n")
cat("  With Superfamily:", sum(!is.na(faird_data$Taxon_SF)),
    sprintf("(%.1f%%)\n\n", 100 * sum(!is.na(faird_data$Taxon_SF)) / nrow(faird_data)))

# --- 1b. AMMOD taxonomy (presence/absence OTUs, Sites 1-3 only) ---
ammod_data <- data_list$ammod_taxonomy %>%
  dplyr::filter(
    Date >= date_start & Date <= date_end,
    Date != date_exclude,
    Site %in% sites,
    !(Device == "AMMOD4" & Date >= ammod4_failure_start & Date <= ammod4_failure_end)
  ) %>%
  dplyr::mutate(
    Taxon_SF = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc & Superfamily != "",
      Superfamily, NA_character_
    )
  )

cat("AMMOD OTUs:", nrow(ammod_data), "\n")
cat("  With Superfamily:", sum(!is.na(ammod_data$Taxon_SF)),
    sprintf("(%.1f%%)\n\n", 100 * sum(!is.na(ammod_data$Taxon_SF)) / nrow(ammod_data)))

# --- 1c. Build taxonomic distance matrix (Superfamily level) ---
cat("Building taxonomic distance matrix...\n")

faird_sf_taxa <- faird_data %>%
  dplyr::filter(!is.na(Taxon_SF)) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily)
faird_sf_taxa <- faird_sf_taxa[!duplicated(faird_sf_taxa$Taxon_SF), , drop = FALSE]

ammod_sf_taxa <- ammod_data %>%
  dplyr::filter(!is.na(Taxon_SF)) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily)
ammod_sf_taxa <- ammod_sf_taxa[!duplicated(ammod_sf_taxa$Taxon_SF), , drop = FALSE]

all_sf_taxa <- dplyr::bind_rows(faird_sf_taxa, ammod_sf_taxa)
all_sf_taxa <- all_sf_taxa[!duplicated(all_sf_taxa$Taxon_SF), , drop = FALSE]

tax_table_sf <- as.data.frame(all_sf_taxa[, c("Taxon_SF", "Class", "Order",
                                               "Suborder", "Infraorder", "Superfamily")])
rownames(tax_table_sf) <- tax_table_sf$Taxon_SF
tax_table_sf            <- tax_table_sf[, setdiff(colnames(tax_table_sf), "Taxon_SF"), drop = FALSE]

taxon_distances_sf <- vegan::taxa2dist(tax_table_sf, varstep = TRUE, check = TRUE)
dist_matrix_sf     <- as.matrix(taxon_distances_sf)

cat(sprintf("Distance matrix: %d x %d | range: %.4f - %.4f\n\n",
            nrow(dist_matrix_sf), ncol(dist_matrix_sf),
            min(taxon_distances_sf), max(taxon_distances_sf)))

# --- 1d. Izsak-Price ΔS function ---
calculate_ip <- function(taxa_A, taxa_B, dist_mat) {
  taxa_A <- taxa_A[taxa_A %in% rownames(dist_mat)]
  taxa_B <- taxa_B[taxa_B %in% rownames(dist_mat)]
  if (length(taxa_A) == 0 | length(taxa_B) == 0) return(NA_real_)
  min_A_to_B <- sapply(taxa_A, function(a)
    if (a %in% taxa_B) 0 else min(dist_mat[a, taxa_B], na.rm = TRUE))
  min_B_to_A <- sapply(taxa_B, function(b)
    if (b %in% taxa_A) 0 else min(dist_mat[b, taxa_A], na.rm = TRUE))
  1 - (mean(c(min_A_to_B, min_B_to_A)) / max(dist_mat, na.rm = TRUE))
}

# --- 1e. 3-day windowed loop (Sites 1-3) ---
cat(sprintf("3-day windowed Izsak-Price (MIN_INDIV = %d FAIRD individuals per window)\n\n",
            MIN_INDIV))

windowed_results <- list()

for (s in sites) {
  for (w in unique(stats::na.omit(window_lookup$Window_ID))) {
    win_dates <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
    win_label <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>%
      dplyr::pull(Window_Label) %>% unique()

    faird_win    <- faird_data %>%
      dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF))
    faird_sf     <- unique(faird_win$Taxon_SF)
    n_faird_indiv <- nrow(faird_win)

    ammod_win <- ammod_data %>%
      dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF))
    ammod_sf  <- unique(ammod_win$Taxon_SF)

    passes_filter <- (n_faird_indiv >= MIN_INDIV & length(ammod_sf) > 0)
    delta_s       <- NA_real_
    if (passes_filter & length(faird_sf) > 0)
      delta_s <- calculate_ip(faird_sf, ammod_sf, dist_matrix_sf)

    windowed_results[[length(windowed_results) + 1]] <- tibble::tibble(
      Period        = win_label,
      Window_ID     = w,
      Site          = s,
      Habitat       = dplyr::if_else(s %in% c("Site1", "Site2"), "Maize", "Meadow"),
      Delta_S       = delta_s,
      N_taxa_FAIRD  = length(faird_sf),
      N_taxa_AMMOD  = length(ammod_sf),
      N_indiv_FAIRD = n_faird_indiv,
      N_shared      = length(intersect(faird_sf, ammod_sf)),
      Filtered      = !passes_filter
    )
  }
}

windowed_df <- dplyr::bind_rows(windowed_results)
cat("Windowed computations:", nrow(windowed_df), "\n")
cat("Valid:", sum(!is.na(windowed_df$Delta_S)),
    "| Filtered:", sum(windowed_df$Filtered), "\n\n")

# --- 1f. Per-site windowed summary ---
cat("--- Site-level windowed summary ---\n\n")

win_by_site <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n_windows    = dplyr::n(),
    mean_DS      = round(mean(Delta_S) * 100, 1),
    sd_DS        = round(stats::sd(Delta_S) * 100, 1),
    se_DS        = round(stats::sd(Delta_S) / sqrt(dplyr::n()) * 100, 2),
    faird_sf_win = round(mean(N_taxa_FAIRD), 1),
    ammod_sf_win = round(mean(N_taxa_AMMOD), 1),
    shared_win   = round(mean(N_shared), 1),
    .groups = "drop"
  )
print(win_by_site, n = Inf)
cat("\n")

# --- 1g. Site-level temporal figure (Fig 25) ---
plot_valid <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::mutate(DS_pct = Delta_S * 100)

win_breaks_x <- sort(unique(plot_valid$Window_ID))

p_site <- ggplot2::ggplot(plot_valid,
                          ggplot2::aes(x = Window_ID, y = DS_pct, color = Site)) +
  ggplot2::geom_line(ggplot2::aes(group = Site), linewidth = 0.8) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  ggplot2::scale_x_continuous(breaks = win_breaks_x,
                               labels = function(x) paste0("W", x)) +
  ggplot2::scale_color_manual(values = custom_colors,
                               labels = c("Site1" = "Site 1 (Maize)",
                                          "Site2" = "Site 2 (Maize)",
                                          "Site3" = "Site 3 (Meadow)",
                                          "Site4" = "Site 4 (Meadow)")) +
  ggplot2::labs(
    title    = "H2 Qualitative — Taxonomic Similarity: FAIRD vs AMMOD — Site Level",
    subtitle = paste0("Izsak-Price index | 3-day windows | Superfamily level | filter: >= ",
                      MIN_INDIV, " FAIRD individuals"),
    x = "Window", y = "Izsak-Price Similarity (%)",
    color = "Site",
    caption = "Dashed line = 50% similarity reference."
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "gray30"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "bottom",
    legend.title     = ggplot2::element_text(face = "bold")
  )

print(p_site)
ggplot2::ggsave("outputs/figures/Fig25_H2_quali_site_level.png",
                plot = p_site, width = 10, height = 6, dpi = 300)
cat("[OK] Fig25_H2_quali_site_level.png\n\n")


# ==============================================================================
# SECTION 2: INTRA-HABITAT LEVEL (§3.3.2.2)
# ==============================================================================

cat("=== SECTION 2: INTRA-HABITAT LEVEL ===\n\n")

maize_w  <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S), Habitat == "Maize") %>%
  dplyr::pull(Delta_S)
meadow_w <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S), Habitat == "Meadow") %>%
  dplyr::pull(Delta_S)

cat(sprintf("Maize:  %.1f%% +/- %.1f%% (n = %d windows)\n",
            mean(maize_w) * 100, stats::sd(maize_w) * 100, length(maize_w)))
cat(sprintf("Meadow: %.1f%% +/- %.1f%% (n = %d windows)\n\n",
            mean(meadow_w) * 100, stats::sd(meadow_w) * 100, length(meadow_w)))

# Levene's test first (then t-test with correct var.equal)
gd <- data.frame(
  val = c(maize_w, meadow_w),
  grp = factor(c(rep("Maize", length(maize_w)), rep("Meadow", length(meadow_w))))
)
lev_result <- car::leveneTest(val ~ grp, data = gd)
lev_p      <- lev_result$`Pr(>F)`[1]
cat(sprintf("Levene's test: F = %.3f, p = %.4f (%s)\n",
            lev_result$`F value`[1], lev_p,
            ifelse(lev_p >= 0.05, "equal variances assumed", "unequal variances")))

tt <- stats::t.test(maize_w, meadow_w, var.equal = (lev_p >= 0.05))
cat(sprintf("t-test: t(%s) = %.3f, p = %.4f\n",
            ifelse(tt$parameter == round(tt$parameter),
                   as.character(tt$parameter),
                   sprintf("%.1f", tt$parameter)),
            tt$statistic, tt$p.value))

cd <- effsize::cohen.d(maize_w, meadow_w)
cat(sprintf("Cohen's d = %.3f (%s)\n\n", cd$estimate, cd$magnitude))

# Intra-Habitat figure (Fig 26 — violin + boxplot)
violin_df <- data.frame(
  Delta_S_pct = c(maize_w * 100, meadow_w * 100),
  Habitat     = factor(c(rep("Maize", length(maize_w)),
                         rep("Meadow", length(meadow_w))))
)

p_violin <- ggplot2::ggplot(violin_df,
                            ggplot2::aes(x = Habitat, y = Delta_S_pct,
                                         fill = Habitat)) +
  ggplot2::geom_violin(alpha = 0.4, trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.15, alpha = 0.8,
                        outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.08, alpha = 0.6, size = 2,
                       ggplot2::aes(color = Habitat)) +
  ggplot2::scale_fill_manual(values  = custom_colors) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::labs(
    title    = "H2 Qualitative — Intra-Habitat Comparison: FAIRD vs AMMOD",
    subtitle = sprintf(
      "t = %.3f, p = %.4f | Cohen's d = %.3f (%s) | Levene's p = %.4f (%s)",
      tt$statistic, tt$p.value, cd$estimate, cd$magnitude,
      lev_p, ifelse(lev_p >= 0.05, "equal variances", "unequal variances")
    ),
    x = "Habitat", y = "Izsak-Price Similarity (%)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "gray30"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "none"
  )

print(p_violin)
ggplot2::ggsave("outputs/figures/Fig26_H2_quali_intra_habitat.png",
                plot = p_violin, width = 8, height = 6, dpi = 300)
cat("[OK] Fig26_H2_quali_intra_habitat.png\n\n")

# Intra-Habitat temporal figure (Fig 27)
plot_intra <- plot_valid %>%
  dplyr::group_by(Window_ID, Habitat) %>%
  dplyr::summarise(DS_pct = mean(DS_pct), .groups = "drop")

p_intra <- ggplot2::ggplot(plot_intra,
                           ggplot2::aes(x = Window_ID, y = DS_pct,
                                        color = Habitat)) +
  ggplot2::geom_line(ggplot2::aes(group = Habitat), linewidth = 0.8) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  ggplot2::scale_x_continuous(breaks = win_breaks_x,
                               labels = function(x) paste0("W", x)) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::labs(
    title    = "H2 Qualitative — Taxonomic Similarity: FAIRD vs AMMOD — Intra-Habitat Level",
    subtitle = sprintf(
      "Maize: %.1f%% ± %.1f%% (n=%d) | Meadow: %.1f%% ± %.1f%% (n=%d)",
      mean(maize_w) * 100, stats::sd(maize_w) * 100, length(maize_w),
      mean(meadow_w) * 100, stats::sd(meadow_w) * 100, length(meadow_w)
    ),
    x = "Window", y = "Izsak-Price Similarity (%)",
    color = "Habitat",
    caption = "Dashed line = 50% similarity reference."
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "gray30"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "bottom",
    legend.title     = ggplot2::element_text(face = "bold")
  )

print(p_intra)
ggplot2::ggsave("outputs/figures/Fig27_H2_quali_intra_habitat_temporal.png",
                plot = p_intra, width = 10, height = 6, dpi = 300)
cat("[OK] Fig27_H2_quali_intra_habitat_temporal.png\n\n")


# ==============================================================================
# SECTION 3: OVERALL LEVEL (§3.3.2.3)
# ==============================================================================

cat("=== SECTION 3: OVERALL LEVEL ===\n\n")

all_vals <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::pull(Delta_S)

n_all   <- length(all_vals)
mean_all <- mean(all_vals)
sd_all   <- stats::sd(all_vals)
ci_lo    <- (mean_all - 1.96 * sd_all / sqrt(n_all)) * 100
ci_hi    <- (mean_all + 1.96 * sd_all / sqrt(n_all)) * 100

cat(sprintf("Overall ΔS: %.1f%% ± %.1f%% (n = %d, 95%% CI: %.1f–%.1f%%)\n\n",
            mean_all * 100, sd_all * 100, n_all, ci_lo, ci_hi))

plot_overall <- plot_valid %>%
  dplyr::group_by(Window_ID) %>%
  dplyr::summarise(DS_pct = mean(DS_pct), .groups = "drop")

p_overall <- ggplot2::ggplot(plot_overall,
                             ggplot2::aes(x = Window_ID, y = DS_pct)) +
  ggplot2::geom_line(linewidth = 0.8, color = "steelblue") +
  ggplot2::geom_point(size = 3, color = "steelblue") +
  ggplot2::geom_hline(yintercept = 50, linetype = "dashed",
                      color = "gray50") +
  ggplot2::geom_hline(yintercept = mean_all * 100, linetype = "dotted",
                      color = "black", linewidth = 0.6) +
  ggplot2::scale_x_continuous(breaks = win_breaks_x,
                               labels = function(x) paste0("W", x)) +
  ggplot2::labs(
    title    = "H2 Qualitative — Taxonomic Similarity: FAIRD vs AMMOD — Overall Level",
    subtitle = sprintf(
      "Izsak-Price | 3-day windows | Superfamily | Overall mean: %.1f%% ± %.1f%% (n=%d, 95%% CI: %.1f–%.1f%%)",
      mean_all * 100, sd_all * 100, n_all, ci_lo, ci_hi
    ),
    x = "Window", y = "Izsak-Price Similarity (%)",
    caption = "Dashed line = 50% reference. Dotted line = overall mean."
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "gray30"),
    panel.grid.minor = ggplot2::element_blank()
  )

print(p_overall)
ggplot2::ggsave("outputs/figures/Fig28_H2_quali_overall_level.png",
                plot = p_overall, width = 10, height = 6, dpi = 300)
cat("[OK] Fig28_H2_quali_overall_level.png\n\n")

# --- Build and export Table 16 ---
cat("--- TABLE 16: Izsak-Price ΔS Summary ---\n\n")

table16_sites <- win_by_site %>%
  dplyr::rename(
    Level     = Site,
    n_windows = n_windows,
    Mean_DS   = mean_DS,
    SD_DS     = sd_DS
  ) %>%
  dplyr::select(Level, Habitat, n_windows, Mean_DS, SD_DS)

table16_intra <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Habitat) %>%
  dplyr::summarise(
    Level     = paste0(dplyr::first(Habitat), " (Intra-Habitat)"),
    n_windows = dplyr::n(),
    Mean_DS   = round(mean(Delta_S) * 100, 1),
    SD_DS     = round(stats::sd(Delta_S) * 100, 1),
    .groups   = "drop"
  ) %>%
  dplyr::select(Level, Habitat, n_windows, Mean_DS, SD_DS)

table16_overall <- tibble::tibble(
  Level     = "Overall",
  Habitat   = "All",
  n_windows = n_all,
  Mean_DS   = round(mean_all * 100, 1),
  SD_DS     = round(sd_all * 100, 1)
)

table16 <- dplyr::bind_rows(table16_sites, table16_intra, table16_overall)
print(table16, n = Inf)

readr::write_csv(table16, "outputs/tables/Table16_H2_quali_IP_summary.csv")
cat("[OK] Table16_H2_quali_IP_summary.csv\n\n")


# ==============================================================================
# SECTION 4: SIZE-THRESHOLD SENSITIVITY (§3.3.2.4)
# ==============================================================================

cat("=== SECTION 4: SIZE-THRESHOLD SENSITIVITY ===\n\n")
cat("Three AMMOD pools: raw (all OTUs), >=3mm, >=5mm. Same FAIRD pool in all three.\n\n")

# --- AMMOD size-filtered subsets (Body_length column) ---
ammod_base <- data_list$ammod_taxonomy %>%
  dplyr::filter(
    Date >= date_start & Date <= date_end,
    Date != date_exclude,
    Site %in% sites,
    !(Device == "AMMOD4" & Date >= ammod4_failure_start & Date <= ammod4_failure_end)
  ) %>%
  dplyr::mutate(
    Taxon_SF = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc & Superfamily != "",
      Superfamily, NA_character_
    )
  )

ammod_raw_sz <- ammod_base
ammod_3mm    <- ammod_base %>% dplyr::filter(!is.na(Body_length) & Body_length >= SIZE_3MM)
ammod_5mm    <- ammod_base %>% dplyr::filter(!is.na(Body_length) & Body_length >= SIZE_5MM)

cat(sprintf("AMMOD raw OTUs:    %d\n", nrow(ammod_raw_sz)))
cat(sprintf("AMMOD >= 3mm OTUs: %d (%.1f%% of raw)\n",
            nrow(ammod_3mm), nrow(ammod_3mm) / nrow(ammod_raw_sz) * 100))
cat(sprintf("AMMOD >= 5mm OTUs: %d (%.1f%% of raw)\n\n",
            nrow(ammod_5mm), nrow(ammod_5mm) / nrow(ammod_raw_sz) * 100))

# --- Windowed loop: three parallel ΔS computations (reuse dist_matrix_sf and calculate_ip) ---
sz_results <- list()

for (s in sites) {
  for (w in unique(stats::na.omit(window_lookup$Window_ID))) {
    win_dates <- window_lookup %>% dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)

    faird_win     <- faird_data %>%
      dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF))
    faird_sf      <- unique(faird_win$Taxon_SF)
    n_faird_indiv <- nrow(faird_win)

    sf_raw <- unique((ammod_raw_sz %>%
                        dplyr::filter(Date %in% win_dates, Site == s,
                                      !is.na(Taxon_SF)))$Taxon_SF)
    sf_3mm <- unique((ammod_3mm %>%
                        dplyr::filter(Date %in% win_dates, Site == s,
                                      !is.na(Taxon_SF)))$Taxon_SF)
    sf_5mm <- unique((ammod_5mm %>%
                        dplyr::filter(Date %in% win_dates, Site == s,
                                      !is.na(Taxon_SF)))$Taxon_SF)

    passes <- (n_faird_indiv >= MIN_INDIV)

    ds_raw <- NA_real_; ds_3mm <- NA_real_; ds_5mm <- NA_real_
    if (passes & length(faird_sf) > 0) {
      if (length(sf_raw) > 0) ds_raw <- calculate_ip(faird_sf, sf_raw, dist_matrix_sf)
      if (length(sf_3mm) > 0) ds_3mm <- calculate_ip(faird_sf, sf_3mm, dist_matrix_sf)
      if (length(sf_5mm) > 0) ds_5mm <- calculate_ip(faird_sf, sf_5mm, dist_matrix_sf)
    }

    sz_results[[length(sz_results) + 1]] <- tibble::tibble(
      Site          = s,
      Window_ID     = w,
      Habitat       = dplyr::if_else(s %in% c("Site1", "Site2"), "Maize", "Meadow"),
      DS_raw        = ds_raw,
      DS_3mm        = ds_3mm,
      DS_5mm        = ds_5mm,
      Taxa_AMMOD_raw = length(sf_raw),
      Taxa_AMMOD_3mm = length(sf_3mm),
      Taxa_AMMOD_5mm = length(sf_5mm)
    )
  }
}

res_df <- dplyr::bind_rows(sz_results) %>%
  dplyr::filter(!is.na(DS_raw) & !is.na(DS_3mm) & !is.na(DS_5mm))

cat(sprintf("Windows with all three ΔS computed: %d\n\n", nrow(res_df)))

# --- Table 17: size-threshold summary ---
cat("--- TABLE 17: Size-Threshold Sensitivity ---\n\n")

table17_sites <- res_df %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n_windows        = dplyr::n(),
    Mean_DS_raw      = round(mean(DS_raw) * 100, 1),
    Mean_DS_3mm      = round(mean(DS_3mm) * 100, 1),
    Mean_DS_5mm      = round(mean(DS_5mm) * 100, 1),
    Delta_raw_to_3mm = round((mean(DS_3mm) - mean(DS_raw)) * 100, 1),
    Delta_3mm_to_5mm = round((mean(DS_5mm) - mean(DS_3mm)) * 100, 1),
    Taxa_raw         = round(mean(Taxa_AMMOD_raw), 1),
    Taxa_3mm         = round(mean(Taxa_AMMOD_3mm), 1),
    Taxa_5mm         = round(mean(Taxa_AMMOD_5mm), 1),
    .groups = "drop"
  )

table17_all <- res_df %>%
  dplyr::summarise(
    Site = "All sites", Habitat = "-",
    n_windows        = dplyr::n(),
    Mean_DS_raw      = round(mean(DS_raw) * 100, 1),
    Mean_DS_3mm      = round(mean(DS_3mm) * 100, 1),
    Mean_DS_5mm      = round(mean(DS_5mm) * 100, 1),
    Delta_raw_to_3mm = round((mean(DS_3mm) - mean(DS_raw)) * 100, 1),
    Delta_3mm_to_5mm = round((mean(DS_5mm) - mean(DS_3mm)) * 100, 1),
    Taxa_raw         = round(mean(Taxa_AMMOD_raw), 1),
    Taxa_3mm         = round(mean(Taxa_AMMOD_3mm), 1),
    Taxa_5mm         = round(mean(Taxa_AMMOD_5mm), 1)
  )

table17 <- dplyr::bind_rows(table17_sites, table17_all)
print(table17, n = Inf)
cat("\nNOTE: Delta_raw_to_3mm = DS_3mm - DS_raw (pp). Positive = filtering >=3mm improves concordance.\n")
cat("      Delta_3mm_to_5mm = DS_5mm - DS_3mm (pp). Same interpretation for stricter threshold.\n\n")

readr::write_csv(table17, "outputs/tables/Table17_H2_quali_size_threshold.csv")
cat("[OK] Table17_H2_quali_size_threshold.csv\n\n")

# --- Sub-3mm Superfamily diagnostic ---
cat("--- Sub-3mm Superfamily diagnostic ---\n\n")
cat("For each site × window: do sub-3mm OTUs add NEW Superfamilies to the raw pool?\n")
cat("(SF present in raw but absent in >=3mm fraction)\n\n")

diag_results <- list()

for (s in sites) {
  for (w in unique(stats::na.omit(window_lookup$Window_ID))) {
    win_dates <- window_lookup %>% dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)

    sf_r <- unique((ammod_raw_sz %>%
                      dplyr::filter(Date %in% win_dates, Site == s,
                                    !is.na(Taxon_SF)))$Taxon_SF)
    sf_m <- unique((ammod_3mm %>%
                      dplyr::filter(Date %in% win_dates, Site == s,
                                    !is.na(Taxon_SF)))$Taxon_SF)

    sf_sub3mm_only <- setdiff(sf_r, sf_m)

    diag_results[[length(diag_results) + 1]] <- tibble::tibble(
      Site                  = s,
      Window_ID             = w,
      Habitat               = dplyr::if_else(s %in% c("Site1", "Site2"), "Maize", "Meadow"),
      N_SF_raw              = length(sf_r),
      N_SF_3mm              = length(sf_m),
      N_SF_sub3mm_exclusive = length(sf_sub3mm_only),
      Pct_new_SF            = dplyr::if_else(
        length(sf_r) > 0,
        round(length(sf_sub3mm_only) / length(sf_r) * 100, 1),
        NA_real_
      ),
      SF_new_names = paste(sf_sub3mm_only, collapse = "; ")
    )
  }
}

diag_df <- dplyr::bind_rows(diag_results) %>%
  dplyr::filter(N_SF_raw > 0)

diag_summary <- diag_df %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n_windows           = dplyr::n(),
    mean_SF_raw         = round(mean(N_SF_raw), 1),
    mean_SF_3mm         = round(mean(N_SF_3mm), 1),
    mean_new_SF         = round(mean(N_SF_sub3mm_exclusive), 2),
    max_new_SF          = max(N_SF_sub3mm_exclusive),
    windows_with_new_SF = sum(N_SF_sub3mm_exclusive > 0),
    mean_pct_new        = round(mean(Pct_new_SF, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(diag_summary, n = Inf)

diag_overall <- diag_df %>%
  dplyr::summarise(
    n_windows           = dplyr::n(),
    mean_SF_raw         = round(mean(N_SF_raw), 1),
    mean_SF_3mm         = round(mean(N_SF_3mm), 1),
    mean_new_SF         = round(mean(N_SF_sub3mm_exclusive), 2),
    max_new_SF          = max(N_SF_sub3mm_exclusive),
    windows_with_new_SF = sum(N_SF_sub3mm_exclusive > 0),
    mean_pct_new        = round(mean(Pct_new_SF, na.rm = TRUE), 1)
  )
cat("\nOverall:\n")
print(diag_overall)

all_new_sf <- diag_df %>%
  dplyr::filter(SF_new_names != "") %>%
  dplyr::pull(SF_new_names) %>%
  strsplit("; ") %>%
  unlist() %>%
  unique() %>%
  sort()

cat("\nUnique Superfamilies appearing ONLY in sub-3mm fraction (across all windows):\n")
if (length(all_new_sf) == 0) {
  cat("  None — sub-3mm fraction adds no new Superfamilies in any window.\n")
} else {
  cat(paste0("  ", all_new_sf, "\n"), sep = "")
}
cat("\n")


# ==============================================================================
# SECTION 5: TAXONOMIC BREADTH — FAIRD ONLY (§3.3.2.5)
# ==============================================================================

cat("=== SECTION 5: TAXONOMIC BREADTH (FAIRD) ===\n\n")

# --- 5a. Table 18: Diptera dominance by Habitat (Ambient) and Site ---
cat("--- TABLE 18: Diptera dominance by Habitat and Site ---\n\n")

# By habitat (Ambient = Maize / Meadow) — #N/C excluded: not a valid taxon
faird_by_ambient <- data_list$faird %>%
  dplyr::filter(is.na(Order) | Order != nc) %>%
  dplyr::group_by(Ambient, Order) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Ambient) %>%
  dplyr::mutate(
    total = sum(n),
    pct   = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Ambient, dplyr::desc(pct))

cat("--- Order composition by Habitat ---\n")
print(faird_by_ambient, n = 30)

cat("\n--- Diptera % by Habitat ---\n")
faird_by_ambient %>%
  dplyr::filter(Order == "Diptera") %>%
  dplyr::select(Ambient, n, total, pct) %>%
  dplyr::rename(Diptera_n = n, Total_n = total, Diptera_pct = pct) %>%
  print()

# By site — #N/C excluded
faird_by_site <- data_list$faird %>%
  dplyr::filter(is.na(Order) | Order != nc) %>%
  dplyr::group_by(Site, Ambient, Order) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    total = sum(n),
    pct   = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Ambient, Site, dplyr::desc(pct))

cat("\n--- Order composition by Site ---\n")
print(faird_by_site, n = 40)

cat("\n--- Diptera % by Site ---\n")
faird_by_site %>%
  dplyr::filter(Order == "Diptera") %>%
  dplyr::select(Site, Ambient, n, total, pct) %>%
  dplyr::rename(Diptera_n = n, Total_n = total, Diptera_pct = pct) %>%
  print()

# Overall reference — #N/C excluded
cat("\n--- Diptera % overall (all sites, #N/C excluded) ---\n")
diptera_overall <- data_list$faird %>%
  dplyr::filter(is.na(Order) | Order != nc) %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(total = sum(n), pct = round(n / total * 100, 1)) %>%
  dplyr::filter(Order == "Diptera")
print(diptera_overall)

# Export Table 18
table18 <- faird_by_ambient %>%
  dplyr::filter(Order == "Diptera") %>%
  dplyr::select(Ambient, n, total, pct) %>%
  dplyr::rename(Habitat = Ambient, Diptera_n = n, Total_n = total, Diptera_pct = pct)

readr::write_csv(table18, "outputs/tables/Table18_H2_quali_diptera_dominance.csv")
cat("\n[OK] Table18_H2_quali_diptera_dominance.csv\n\n")

# --- 5b. Table 19: Taxonomic resolution by Order ---
cat("--- TABLE 19: Taxonomic resolution by Order ---\n\n")

is_valid <- function(x) !is.na(x) & x != "#N/C" & x != ""

# Filter #N/C at the source: these individuals have no valid Order-level classification
# and must not appear as a row in Table 19 or inflate the Overall totals.
resolution_data <- data_list$faird %>%
  dplyr::filter(is.na(Order) | Order != nc) %>%
  dplyr::mutate(
    finest_level = dplyr::case_when(
      is_valid(Genus)       ~ "Genus",
      is_valid(Family)      ~ "Family",
      is_valid(Superfamily) ~ "Superfamily",
      is_valid(Infraorder)  ~ "Infraorder",
      is_valid(Suborder)    ~ "Suborder",
      is_valid(Order)       ~ "Order",
      TRUE                  ~ "Class_or_unresolved"
    ),
    finest_level = factor(finest_level,
                          levels = c("Genus", "Family", "Superfamily",
                                     "Infraorder", "Suborder", "Order",
                                     "Class_or_unresolved"))
  )

# By Order (full detail)
resolution_by_order <- resolution_data %>%
  dplyr::group_by(Order, finest_level) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Order) %>%
  dplyr::mutate(total = sum(n), pct = round(n / total * 100, 1)) %>%
  dplyr::ungroup()

resolution_table_order <- resolution_by_order %>%
  dplyr::select(Order, finest_level, n, pct, total) %>%
  tidyr::pivot_wider(
    names_from  = finest_level,
    values_from = c(n, pct),
    values_fill = list(n = 0L, pct = 0.0)
  ) %>%
  dplyr::select(
    Order, total,
    dplyr::any_of(c(
      "n_Genus",               "pct_Genus",
      "n_Family",              "pct_Family",
      "n_Superfamily",         "pct_Superfamily",
      "n_Infraorder",          "pct_Infraorder",
      "n_Suborder",            "pct_Suborder",
      "n_Order",               "pct_Order",
      "n_Class_or_unresolved", "pct_Class_or_unresolved"
    ))
  ) %>%
  dplyr::arrange(dplyr::desc(total))

cat("--- Resolution by Order (detailed) ---\n")
print(resolution_table_order, n = 50)

# Condensed (PCV14-style: Genus_ID / Family_only / Superfamily_or_lower)
resolution_condensed <- resolution_data %>%
  dplyr::mutate(
    resolution_class = dplyr::case_when(
      finest_level == "Genus"  ~ "Genus_ID",
      finest_level == "Family" ~ "Family_only",
      TRUE                     ~ "Superfamily_or_lower"
    )
  ) %>%
  dplyr::group_by(Order, resolution_class) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Order) %>%
  dplyr::mutate(total = sum(n), pct = round(n / total * 100, 1)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    names_from  = resolution_class,
    values_from = c(n, pct),
    values_fill = list(n = 0L, pct = 0.0)
  ) %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::select(
    Order, total,
    dplyr::any_of(c(
      "pct_Genus_ID",             "n_Genus_ID",
      "pct_Family_only",          "n_Family_only",
      "pct_Superfamily_or_lower", "n_Superfamily_or_lower"
    ))
  )

cat("\n--- Resolution condensed by Order ---\n")
print(resolution_condensed, n = 50)

# Overall condensed
resolution_overall <- resolution_data %>%
  dplyr::mutate(
    resolution_class = dplyr::case_when(
      finest_level == "Genus"  ~ "Genus_ID",
      finest_level == "Family" ~ "Family_only",
      TRUE                     ~ "Superfamily_or_lower"
    )
  ) %>%
  dplyr::group_by(resolution_class) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(total = sum(n), pct = round(n / total * 100, 1)) %>%
  tidyr::pivot_wider(
    names_from  = resolution_class,
    values_from = c(n, pct),
    values_fill = list(n = 0L, pct = 0.0)
  )

cat("\n--- Overall resolution (all Orders combined) ---\n")
print(resolution_overall)

readr::write_csv(resolution_table_order, "outputs/tables/Table19_H2_quali_resolution_by_order.csv")
readr::write_csv(resolution_condensed,   "outputs/tables/Table19_H2_quali_resolution_condensed.csv")
cat("\n[OK] Table19_H2_quali_resolution_by_order.csv\n")
cat("[OK] Table19_H2_quali_resolution_condensed.csv\n\n")


# ==============================================================================
# WRAP UP
# ==============================================================================

cat("==============================================================================\n")
cat("H2_quali_validation.R COMPLETE\n")
cat("==============================================================================\n\n")

cat("FIGURES:\n")
cat("  outputs/figures/Fig25_H2_quali_site_level.png\n")
cat("  outputs/figures/Fig26_H2_quali_intra_habitat.png\n")
cat("  outputs/figures/Fig27_H2_quali_intra_habitat_temporal.png\n")
cat("  outputs/figures/Fig28_H2_quali_overall_level.png\n\n")

cat("TABLES:\n")
cat("  outputs/tables/Table16_H2_quali_IP_summary.csv\n")
cat("  outputs/tables/Table17_H2_quali_size_threshold.csv\n")
cat("  outputs/tables/Table18_H2_quali_diptera_dominance.csv\n")
cat("  outputs/tables/Table19_H2_quali_resolution_by_order.csv\n")
cat("  outputs/tables/Table19_H2_quali_resolution_condensed.csv\n\n")

cat("Script finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

while (sink.number() > 0) sink()
