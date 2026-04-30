# ==============================================================================
# H1: QUALITATIVE INTRA-SYSTEM CONSISTENCY — §3.2.2
# ==============================================================================
# Metrics:
#   FAIRD — Morisita-Horn (1 - vegdist(method="horn")):
#            abundance-based, appropriate for FAIRD count data
#   AMMOD — Izsak-Price delta-S (presence/absence + taxonomic distance):
#            presence/absence-based, required for AMMOD OTU data
#
# Temporal resolution: daily (unfiltered) + 3-day windows (filtered >= 3)
#
# Pairs:
#   Intra-habitat FAIRD: FAIRD1<->FAIRD2 (Maize), FAIRD3<->FAIRD4 (Meadow)
#   Intra-habitat AMMOD: AMMOD1<->AMMOD2 (Maize), AMMOD3<->AMMOD4 (Meadow)
#   Inter-habitat FAIRD: Pool(FAIRD1+2) vs Pool(FAIRD3+4)
#   Inter-habitat AMMOD: Pool(AMMOD1+2) vs Pool(AMMOD3+4)
#
# Outputs:
#   Figures 12-15  -> outputs/figures/
#   Tables 11-12   -> outputs/tables/
#   Console log    -> outputs/console/
#
# Author: Juan Chiavassa & Claude
# ==============================================================================

library(tidyverse)
library(scales)
library(vegan)     # vegdist (Morisita-Horn), taxa2dist (Izsak-Price)
library(effsize)   # cohen.d()
library(car)       # leveneTest()
library(here)

# Anchor here to the project root regardless of where this script is executed from.
# here::i_am() walks upward from the script location until it finds a directory
# that contains R/H1_quali_consistency.R — that directory becomes the project root.
here::i_am("R/H1_quali_consistency.R")

# ==============================================================================
# CONSOLE OUTPUT
# ==============================================================================

dir.create("outputs/figures",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",   showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console",  showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/H1_quali_consistency_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("H1: QUALITATIVE INTRA-SYSTEM CONSISTENCY — §3.2.2\n")
cat("Morisita-Horn (FAIRD) | Izsak-Price (AMMOD)\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source("R/00_load_data.R")
source("R/00_custom_colors.R")


# ==============================================================================
# SECTION 0: SHARED TEMPORAL SETUP
# ==============================================================================

cat("=== SECTION 0: SHARED TEMPORAL SETUP ===\n\n")

nc <- "#N/C"

date_start        <- as.Date("2023-08-23")
date_end          <- as.Date("2023-09-13")
date_exclude      <- as.Date("2023-09-04")   # bottle change day (AMMOD only)
ammod4_fail_start <- as.Date("2023-08-28")
ammod4_fail_end   <- as.Date("2023-09-03")

all_dates_faird <- seq.Date(date_start, date_end, by = "day")
all_dates_ammod <- all_dates_faird[all_dates_faird != date_exclude]

MIN_INDIV <- 3   # FAIRD: minimum individuals per device per window
MIN_OTUS  <- 3   # AMMOD: minimum OTUs per device per window

window_breaks <- seq(date_start, date_end + 2, by = 3)
window_lookup <- tibble::tibble(Date = all_dates_faird) %>%
  dplyr::mutate(
    Window_ID    = as.integer(cut(Date, breaks = window_breaks,
                                  right = FALSE, labels = FALSE)),
    Window_Start = window_breaks[Window_ID],
    Window_End   = pmin(Window_Start + 2, date_end),
    Window_Label = paste0("W", Window_ID, " (",
                          format(Window_Start, "%b %d"), "-",
                          format(Window_End,   "%b %d"), ")")
  )

cat("3-day windows defined:\n")
window_lookup %>%
  dplyr::group_by(Window_ID, Window_Label) %>%
  dplyr::summarise(n_days = dplyr::n(), .groups = "drop") %>%
  print(n = Inf)
cat("\n")


# ==============================================================================
# SECTION 1: FAIRD INTRA-HABITAT — MORISITA-HORN (§3.2.2.1)
# ==============================================================================
# Source: H1_consistency_FAIRD_taxonomic.R, Sections 2-8
# Pairs:  FAIRD1<->FAIRD2 (Maize), FAIRD3<->FAIRD4 (Meadow)
# ==============================================================================

cat("=== SECTION 1: FAIRD INTRA-HABITAT — MORISITA-HORN (§3.2.2.1) ===\n\n")

# --- 1.1 Build taxon columns and attach windows ---

faird_raw <- data_list$faird

cat("FAIRD individuals loaded:", nrow(faird_raw), "\n")
cat("Devices:", paste(sort(unique(faird_raw$Device)), collapse = ", "), "\n")
cat("Date range:", as.character(min(faird_raw$Date)),
    "to", as.character(max(faird_raw$Date)), "\n\n")

faird_taxo <- faird_raw %>%
  dplyr::mutate(
    Taxon_Genus       = dplyr::if_else(!is.na(Genus)       & Genus       != nc,
                                       Genus,       NA_character_),
    Taxon_Superfamily = dplyr::if_else(!is.na(Superfamily) & Superfamily != nc,
                                       Superfamily, NA_character_),
    Taxon_Order       = dplyr::if_else(!is.na(Order)       & Order       != nc,
                                       Order,       NA_character_)
  ) %>%
  dplyr::left_join(window_lookup, by = "Date")

n_total <- nrow(faird_taxo)
cat("Taxonomic coverage (n =", n_total, "individuals):\n")
cat(sprintf("  Order:       %d (%.1f%%)\n",
            sum(!is.na(faird_taxo$Taxon_Order)),
            100 * sum(!is.na(faird_taxo$Taxon_Order))       / n_total))
cat(sprintf("  Superfamily: %d (%.1f%%)\n",
            sum(!is.na(faird_taxo$Taxon_Superfamily)),
            100 * sum(!is.na(faird_taxo$Taxon_Superfamily)) / n_total))
cat(sprintf("  Genus:       %d (%.1f%%)\n",
            sum(!is.na(faird_taxo$Taxon_Genus)),
            100 * sum(!is.na(faird_taxo$Taxon_Genus))       / n_total))
cat("\n")

# --- 1.2 Morisita-Horn helper functions ---

compute_mh <- function(data, dev_a, dev_b, taxon_col) {
  day_data <- data %>%
    dplyr::filter(Device %in% c(dev_a, dev_b), !is.na(.data[[taxon_col]])) %>%
    dplyr::group_by(Device, .data[[taxon_col]]) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  if (length(unique(day_data$Device)) < 2) return(NA_real_)

  comm_matrix <- day_data %>%
    tidyr::pivot_wider(names_from  = dplyr::all_of(taxon_col),
                       values_from = count, values_fill = 0) %>%
    dplyr::arrange(Device) %>%
    dplyr::select(-Device)

  if (any(rowSums(comm_matrix) == 0)) return(NA_real_)
  1 - as.numeric(vegan::vegdist(as.matrix(comm_matrix), method = "horn"))
}

get_counts_mh <- function(data, device_name, taxon_col) {
  sub <- data %>%
    dplyr::filter(Device == device_name, !is.na(.data[[taxon_col]]))
  list(n_indiv = nrow(sub), n_taxa = length(unique(sub[[taxon_col]])))
}

# --- 1.3 Daily Morisita-Horn (unfiltered) ---

cat("--- 1.3 Daily Morisita-Horn ---\n\n")

faird_pairs    <- list(
  list(dev_a = "FAIRD1", dev_b = "FAIRD2", habitat = "Maize"),
  list(dev_a = "FAIRD3", dev_b = "FAIRD4", habitat = "Meadow")
)
tax_levels     <- c("Taxon_Order", "Taxon_Superfamily", "Taxon_Genus")
tax_labels     <- c("Order",       "Superfamily",        "Genus")

faird_daily <- list()

for (p in seq_along(faird_pairs)) {
  pair <- faird_pairs[[p]]
  cat(sprintf("  %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))

  for (lev in seq_along(tax_levels)) {
    tax_col   <- tax_levels[lev]
    tax_label <- tax_labels[lev]

    for (d in all_dates_faird) {
      date_val <- as.Date(d, origin = "1970-01-01")
      day_data <- faird_taxo %>% dplyr::filter(Date == date_val)
      mh_sim   <- compute_mh(day_data, pair$dev_a, pair$dev_b, tax_col)
      counts_a <- get_counts_mh(day_data, pair$dev_a, tax_col)
      counts_b <- get_counts_mh(day_data, pair$dev_b, tax_col)

      faird_daily[[length(faird_daily) + 1]] <- tibble::tibble(
        Resolution = "Daily",
        Pair       = paste0(pair$dev_a, " vs ", pair$dev_b),
        Habitat    = pair$habitat,
        Period     = as.character(date_val),
        Tax_Level  = tax_label,
        MH_Sim     = mh_sim,
        N_indiv_A  = counts_a$n_indiv,
        N_indiv_B  = counts_b$n_indiv,
        N_taxa_A   = counts_a$n_taxa,
        N_taxa_B   = counts_b$n_taxa,
        Filtered   = FALSE
      )
    }
  }
}

faird_daily_df <- dplyr::bind_rows(faird_daily)
cat(sprintf("\nDaily MH: %d rows | Valid: %d\n\n",
            nrow(faird_daily_df), sum(!is.na(faird_daily_df$MH_Sim))))

# --- 1.4 3-day windowed Morisita-Horn (>= MIN_INDIV individuals) ---

cat(sprintf("--- 1.4 Windowed MH (>= %d individuals per device) ---\n\n", MIN_INDIV))

faird_windowed <- list()

for (p in seq_along(faird_pairs)) {
  pair <- faird_pairs[[p]]
  cat(sprintf("  %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))

  for (lev in seq_along(tax_levels)) {
    tax_col   <- tax_levels[lev]
    tax_label <- tax_labels[lev]

    for (w in unique(stats::na.omit(faird_taxo$Window_ID))) {
      win_data  <- faird_taxo %>% dplyr::filter(Window_ID == w)
      win_label <- unique(win_data$Window_Label)[1]
      counts_a  <- get_counts_mh(win_data, pair$dev_a, tax_col)
      counts_b  <- get_counts_mh(win_data, pair$dev_b, tax_col)

      passes <- counts_a$n_indiv >= MIN_INDIV & counts_b$n_indiv >= MIN_INDIV
      mh_sim <- if (passes) compute_mh(win_data, pair$dev_a, pair$dev_b, tax_col) else NA_real_

      faird_windowed[[length(faird_windowed) + 1]] <- tibble::tibble(
        Resolution = "3-day window",
        Pair       = paste0(pair$dev_a, " vs ", pair$dev_b),
        Habitat    = pair$habitat,
        Period     = win_label,
        Tax_Level  = tax_label,
        MH_Sim     = mh_sim,
        N_indiv_A  = counts_a$n_indiv,
        N_indiv_B  = counts_b$n_indiv,
        N_taxa_A   = counts_a$n_taxa,
        N_taxa_B   = counts_b$n_taxa,
        Filtered   = !passes
      )
    }
  }
}

faird_windowed_df <- dplyr::bind_rows(faird_windowed)
cat(sprintf("\nWindowed MH: %d rows | Valid: %d | Filtered: %d\n\n",
            nrow(faird_windowed_df),
            sum(!is.na(faird_windowed_df$MH_Sim)),
            sum(faird_windowed_df$Filtered)))

cat("Filter impact by pair and level:\n")
faird_windowed_df %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_windows    = dplyr::n(),
    n_valid      = sum(!is.na(MH_Sim)),
    n_filtered   = sum(Filtered),
    pct_retained = round(100 * n_valid / n_windows, 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")), Habitat) %>%
  print(n = Inf)
cat("\n")

# --- 1.5 Summary statistics ---

cat("--- 1.5 Daily summary (unfiltered) ---\n\n")

faird_daily_summary <- faird_daily_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_days    = dplyr::n(),
    mean_MH   = round(mean(MH_Sim) * 100, 1),
    sd_MH     = round(sd(MH_Sim)   * 100, 1),
    median_MH = round(median(MH_Sim) * 100, 1),
    min_MH    = round(min(MH_Sim)  * 100, 1),
    max_MH    = round(max(MH_Sim)  * 100, 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")), Habitat)
print(faird_daily_summary, n = Inf, width = 200)

cat("\n--- Windowed summary (filtered) ---\n\n")

faird_windowed_summary <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_windows    = dplyr::n(),
    mean_MH      = round(mean(MH_Sim)   * 100, 1),
    sd_MH        = round(sd(MH_Sim)     * 100, 1),
    median_MH    = round(median(MH_Sim) * 100, 1),
    min_MH       = round(min(MH_Sim)    * 100, 1),
    max_MH       = round(max(MH_Sim)    * 100, 1),
    mean_indiv_A = round(mean(N_indiv_A), 1),
    mean_indiv_B = round(mean(N_indiv_B), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")), Habitat)
print(faird_windowed_summary, n = Inf, width = 200)
cat("\n")

# --- 1.6 Temporal trend (windowed, Superfamily) ---

cat("--- 1.6 Temporal trend (windowed, Superfamily) ---\n\n")

faird_sfam_trend <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily") %>%
  dplyr::mutate(Window_Num = as.numeric(gsub("W(\\d+).*", "\\1", Period)))

if (nrow(faird_sfam_trend) >= 5) {
  trend_all <- lm(MH_Sim ~ Window_Num, data = faird_sfam_trend)
  ts        <- summary(trend_all)
  cat(sprintf("Overall: slope = %.4f/window (%.2f pp/win) | R2 = %.3f | p = %.4f | n = %d\n\n",
              coef(trend_all)[2], coef(trend_all)[2] * 100,
              ts$r.squared, coef(ts)[2, 4], nrow(faird_sfam_trend)))

  for (pair_name in unique(faird_sfam_trend$Pair)) {
    pd <- faird_sfam_trend %>% dplyr::filter(Pair == pair_name)
    if (nrow(pd) >= 4) {
      m <- lm(MH_Sim ~ Window_Num, data = pd)
      s <- summary(m)
      cat(sprintf("  %s: slope = %.4f (%.2f pp/win) | R2 = %.3f | p = %.4f | n = %d\n",
                  pair_name, coef(m)[2], coef(m)[2] * 100,
                  s$r.squared, coef(s)[2, 4], nrow(pd)))
    }
  }
  cat("\n")
}

# --- 1.7 Habitat comparison (windowed, Superfamily) ---

cat("--- 1.7 Habitat comparison (windowed, Superfamily) ---\n\n")

maize_mh  <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily", Habitat == "Maize") %>%
  dplyr::pull(MH_Sim)
meadow_mh <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily", Habitat == "Meadow") %>%
  dplyr::pull(MH_Sim)

if (length(maize_mh) >= 3 && length(meadow_mh) >= 3) {
  grp_df   <- data.frame(
    val = c(maize_mh, meadow_mh),
    grp = factor(c(rep("Maize", length(maize_mh)), rep("Meadow", length(meadow_mh))))
  )
  levene_p <- car::leveneTest(val ~ grp, data = grp_df)$`Pr(>F)`[1]
  tt       <- t.test(maize_mh, meadow_mh, var.equal = levene_p >= 0.05)
  cd       <- effsize::cohen.d(maize_mh, meadow_mh)

  cat(sprintf("Levene p = %.4f -> %s\n",
              levene_p, ifelse(levene_p < 0.05, "Welch's t-test", "Student's t-test")))
  cat(sprintf("Maize:  %.1f%% +/- %.1f%% (n = %d windows)\n",
              mean(maize_mh) * 100, sd(maize_mh) * 100, length(maize_mh)))
  cat(sprintf("Meadow: %.1f%% +/- %.1f%% (n = %d windows)\n",
              mean(meadow_mh) * 100, sd(meadow_mh) * 100, length(meadow_mh)))
  cat(sprintf("t = %.3f | df = %.1f | p = %.4f\n",
              tt$statistic, tt$parameter, tt$p.value))
  cat(sprintf("Cohen's d = %.3f (%s)\n\n", cd$estimate, cd$magnitude))
}

# --- 1.8 Figure 12: Windowed Superfamily temporal ---

cat("--- 1.8 Figure 12 ---\n")

f12_data <- faird_windowed_df %>%
  dplyr::filter(Tax_Level == "Superfamily", !is.na(MH_Sim)) %>%
  dplyr::mutate(
    MH_pct     = MH_Sim * 100,
    Window_Num = as.numeric(gsub("W(\\d+).*", "\\1", Period))
  )

f12_stats <- f12_data %>%
  dplyr::group_by(Habitat) %>%
  dplyr::summarise(mean_MH = round(mean(MH_pct), 1),
                   sd_MH   = round(sd(MH_pct),   1),
                   n       = dplyr::n(), .groups = "drop")

subtitle_f12 <- f12_stats %>%
  dplyr::mutate(lab = sprintf("%s: %.1f%% ± %.1f%% (n=%d)",
                              Habitat, mean_MH, sd_MH, n)) %>%
  dplyr::pull(lab) %>%
  paste(collapse = " | ")

win_labs_f12 <- f12_data %>%
  dplyr::distinct(Window_Num, Period) %>%
  dplyr::arrange(Window_Num) %>%
  dplyr::mutate(label = gsub("W\\d+ \\((.*)\\)", "\\1", Period))

fig12 <- ggplot(f12_data, aes(x = Window_Num, y = MH_pct, color = Habitat)) +
  geom_line(aes(group = Pair), linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_x_continuous(breaks = win_labs_f12$Window_Num,
                     labels = win_labs_f12$label) +
  scale_color_manual(
    values = custom_colors,
    labels = c("Maize"  = "FAIRD1 vs FAIRD2 (Maize)",
               "Meadow" = "FAIRD3 vs FAIRD4 (Meadow)")
  ) +
  labs(
    title    = "H1 Qualitative: FAIRD Intra-Habitat Consistency — Superfamily Level",
    subtitle = subtitle_f12,
    x        = "3-Day Window",
    y        = "Morisita-Horn Similarity (%)",
    color    = "Pair",
    caption  = "Dashed line = 50% similarity reference."
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

ggsave("outputs/figures/Fig12_H1_quali_FAIRD_intra_MH_windowed_Superfamily.png",
       plot = fig12, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/figures/Fig12_H1_quali_FAIRD_intra_MH_windowed_Superfamily.png\n\n")

# --- 1.9 Figure 13: Multi-level boxplot (windowed) ---

cat("--- 1.9 Figure 13 ---\n")

f13_data <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::mutate(
    MH_pct     = MH_Sim * 100,
    Tax_Level  = factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")),
    Pair_label = dplyr::case_when(
      Pair == "FAIRD1 vs FAIRD2" ~ "FAIRD1 vs FAIRD2\n(Maize)",
      Pair == "FAIRD3 vs FAIRD4" ~ "FAIRD3 vs FAIRD4\n(Meadow)",
      TRUE ~ Pair
    )
  )

fig13 <- ggplot(f13_data, aes(x = Tax_Level, y = MH_pct, fill = Habitat)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  facet_wrap(~Pair_label) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = custom_colors) +
  labs(
    title    = "H1 Qualitative: FAIRD Intra-Habitat Consistency — By Taxonomic Level",
    subtitle = paste0("Morisita-Horn — 3-Day Windows — filtered: >= ",
                      MIN_INDIV, " individuals per device per window"),
    x        = "Taxonomic Level",
    y        = "Morisita-Horn Similarity (%)",
    fill     = "Habitat"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("outputs/figures/Fig13_H1_quali_FAIRD_intra_MH_windowed_by_level.png",
       plot = fig13, width = 9, height = 6, dpi = 300)
cat("[OK] outputs/figures/Fig13_H1_quali_FAIRD_intra_MH_windowed_by_level.png\n\n")


# ==============================================================================
# SECTION 2: AMMOD INTRA-HABITAT — IZSAK-PRICE (§3.2.2.1)
# ==============================================================================
# Source: H1_consistency_AMMOD_taxonomic.R, Sections 1-7b
# Pairs:  AMMOD1<->AMMOD2 (Maize), AMMOD3<->AMMOD4 (Meadow)
# ==============================================================================

cat("=== SECTION 2: AMMOD INTRA-HABITAT — IZSAK-PRICE (§3.2.2.1) ===\n\n")

# --- 2.1 Load AMMOD taxonomy data ---

ammod_raw <- data_list$ammod_taxonomy %>%
  dplyr::filter(
    Date >= date_start & Date <= date_end,
    Date != date_exclude,
    !(Device == "AMMOD4" & Date >= ammod4_fail_start & Date <= ammod4_fail_end)
  ) %>%
  dplyr::select(Date, Device, Site, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus, Species) %>%
  dplyr::mutate(
    Taxon_SF      = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc & Superfamily != "",
      Superfamily, NA_character_),
    Taxon_Genus   = dplyr::if_else(
      !is.na(Genus)       & Genus       != nc & Genus       != "",
      Genus,       NA_character_),
    Taxon_Species = dplyr::if_else(
      !is.na(Species)     & Species     != nc & Species     != "",
      Species,     NA_character_)
  )

cat("AMMOD OTUs loaded:", nrow(ammod_raw), "\n")
cat(sprintf("  With Superfamily: %d (%.1f%%)\n",
            sum(!is.na(ammod_raw$Taxon_SF)),
            100 * sum(!is.na(ammod_raw$Taxon_SF)) / nrow(ammod_raw)))
cat(sprintf("  With Genus:       %d (%.1f%%)\n",
            sum(!is.na(ammod_raw$Taxon_Genus)),
            100 * sum(!is.na(ammod_raw$Taxon_Genus)) / nrow(ammod_raw)))
cat(sprintf("  With Species:     %d (%.1f%%)\n",
            sum(!is.na(ammod_raw$Taxon_Species)),
            100 * sum(!is.na(ammod_raw$Taxon_Species)) / nrow(ammod_raw)))

cat("\nPer-device coverage:\n")
ammod_raw %>%
  dplyr::group_by(Device) %>%
  dplyr::summarise(
    n_otus = dplyr::n(),
    n_sf   = sum(!is.na(Taxon_SF)),
    pct_sf = round(100 * n_sf / dplyr::n(), 1),
    n_days = length(unique(Date)),
    .groups = "drop"
  ) %>%
  print()
cat("\nNote: AMMOD4 reduced to 13 valid taxo days (mechanical failure Aug 28-Sep 3).\n\n")

# --- 2.2 Taxonomic distance matrices ---

cat("--- 2.2 Building taxonomic distance matrices ---\n\n")

build_dist_matrix <- function(data, taxon_col, hier_cols) {
  taxa_df <- data %>%
    dplyr::filter(!is.na(.data[[taxon_col]])) %>%
    dplyr::select(dplyr::all_of(c(taxon_col, hier_cols)))
  taxa_df   <- taxa_df[!duplicated(taxa_df[[taxon_col]]), , drop = FALSE]
  tax_table <- as.data.frame(taxa_df[, hier_cols, drop = FALSE])
  rownames(tax_table) <- taxa_df[[taxon_col]]
  as.matrix(vegan::taxa2dist(tax_table, varstep = TRUE, check = TRUE))
}

hier_sf  <- c("Class", "Order", "Suborder", "Infraorder", "Superfamily")
hier_gen <- c("Class", "Order", "Suborder", "Infraorder",
              "Superfamily", "Family", "Subfamily", "Genus")
hier_sp  <- c("Class", "Order", "Suborder", "Infraorder",
              "Superfamily", "Family", "Subfamily", "Genus", "Species")

dist_sf  <- build_dist_matrix(ammod_raw, "Taxon_SF",      hier_sf)
dist_gen <- build_dist_matrix(ammod_raw, "Taxon_Genus",   hier_gen)
dist_sp  <- build_dist_matrix(ammod_raw, "Taxon_Species", hier_sp)

cat(sprintf("Superfamily matrix: %d x %d | max dist = %.4f\n",
            nrow(dist_sf),  ncol(dist_sf),  max(dist_sf)))
cat(sprintf("Genus matrix:       %d x %d | max dist = %.4f\n",
            nrow(dist_gen), ncol(dist_gen), max(dist_gen)))
cat(sprintf("Species matrix:     %d x %d | max dist = %.4f\n\n",
            nrow(dist_sp),  ncol(dist_sp),  max(dist_sp)))

# --- 2.3 Izsak-Price function ---

calculate_ip <- function(taxa_a, taxa_b, dist_matrix) {
  taxa_a <- taxa_a[taxa_a %in% rownames(dist_matrix)]
  taxa_b <- taxa_b[taxa_b %in% rownames(dist_matrix)]
  if (length(taxa_a) == 0 | length(taxa_b) == 0) return(NA_real_)

  min_a <- sapply(taxa_a, function(a) {
    if (a %in% taxa_b) 0 else min(dist_matrix[a, taxa_b], na.rm = TRUE)
  })
  min_b <- sapply(taxa_b, function(b) {
    if (b %in% taxa_a) 0 else min(dist_matrix[b, taxa_a], na.rm = TRUE)
  })
  1 - (mean(c(min_a, min_b)) / max(dist_matrix, na.rm = TRUE))
}

# --- 2.4 Daily Izsak-Price (all levels) ---

cat("--- 2.4 Daily Izsak-Price (Species / Genus / Superfamily) ---\n\n")

ammod_pairs <- list(
  list(dev_a = "AMMOD1", dev_b = "AMMOD2", habitat = "Maize"),
  list(dev_a = "AMMOD3", dev_b = "AMMOD4", habitat = "Meadow")
)

ammod_tax_levels <- list(
  list(col = "Taxon_Species",  label = "Species",     dmat = dist_sp),
  list(col = "Taxon_Genus",    label = "Genus",        dmat = dist_gen),
  list(col = "Taxon_SF",       label = "Superfamily",  dmat = dist_sf)
)

ammod_daily <- list()

for (p in seq_along(ammod_pairs)) {
  pair <- ammod_pairs[[p]]
  cat(sprintf("  %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))

  for (lev in ammod_tax_levels) {
    for (d in all_dates_ammod) {
      date_val <- as.Date(d, origin = "1970-01-01")

      taxa_a <- ammod_raw %>%
        dplyr::filter(Date == date_val, Device == pair$dev_a,
                      !is.na(.data[[lev$col]])) %>%
        dplyr::pull(.data[[lev$col]]) %>% unique()
      taxa_b <- ammod_raw %>%
        dplyr::filter(Date == date_val, Device == pair$dev_b,
                      !is.na(.data[[lev$col]])) %>%
        dplyr::pull(.data[[lev$col]]) %>% unique()

      delta_s <- if (length(taxa_a) > 0 & length(taxa_b) > 0)
        calculate_ip(taxa_a, taxa_b, lev$dmat) else NA_real_

      ammod_daily[[length(ammod_daily) + 1]] <- tibble::tibble(
        Resolution = "Daily",
        Pair       = paste0(pair$dev_a, " vs ", pair$dev_b),
        Habitat    = pair$habitat,
        Date       = date_val,
        Period     = as.character(date_val),
        Tax_Level  = lev$label,
        Delta_S    = delta_s,
        N_taxa_A   = length(taxa_a),
        N_taxa_B   = length(taxa_b),
        N_shared   = length(intersect(taxa_a, taxa_b))
      )
    }
  }
}

ammod_daily_df <- dplyr::bind_rows(ammod_daily)
cat(sprintf("\nDaily IP: %d rows | Valid: %d\n\n",
            nrow(ammod_daily_df), sum(!is.na(ammod_daily_df$Delta_S))))

# --- 2.5 3-day windowed Izsak-Price (Superfamily only, >= MIN_OTUS) ---

cat(sprintf("--- 2.5 Windowed IP (Superfamily, >= %d OTUs per device) ---\n\n", MIN_OTUS))

ammod_windowed <- list()

for (p in seq_along(ammod_pairs)) {
  pair <- ammod_pairs[[p]]
  cat(sprintf("  %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))

  for (w in unique(window_lookup$Window_ID)) {
    win_dates <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
    win_dates <- win_dates[win_dates != date_exclude]
    win_label <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>% dplyr::pull(Window_Label) %>% unique()

    sf_a <- ammod_raw %>%
      dplyr::filter(Date %in% win_dates, Device == pair$dev_a,
                    !is.na(Taxon_SF)) %>%
      dplyr::pull(Taxon_SF) %>% unique()
    sf_b <- ammod_raw %>%
      dplyr::filter(Date %in% win_dates, Device == pair$dev_b,
                    !is.na(Taxon_SF)) %>%
      dplyr::pull(Taxon_SF) %>% unique()

    passes  <- length(sf_a) >= MIN_OTUS & length(sf_b) >= MIN_OTUS
    delta_s <- if (passes) calculate_ip(sf_a, sf_b, dist_sf) else NA_real_

    ammod_windowed[[length(ammod_windowed) + 1]] <- tibble::tibble(
      Resolution = "3-day window",
      Pair       = paste0(pair$dev_a, " vs ", pair$dev_b),
      Habitat    = pair$habitat,
      Period     = win_label,
      Window_ID  = w,
      Tax_Level  = "Superfamily",
      Delta_S    = delta_s,
      N_taxa_A   = length(sf_a),
      N_taxa_B   = length(sf_b),
      N_shared   = length(intersect(sf_a, sf_b)),
      Filtered   = !passes
    )
  }
}

ammod_windowed_df <- dplyr::bind_rows(ammod_windowed)
cat(sprintf("\nWindowed IP: %d rows | Valid: %d | Filtered: %d\n\n",
            nrow(ammod_windowed_df),
            sum(!is.na(ammod_windowed_df$Delta_S)),
            sum(ammod_windowed_df$Filtered)))

cat("Filter impact:\n")
ammod_windowed_df %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    n_windows    = dplyr::n(),
    n_valid      = sum(!is.na(Delta_S)),
    n_filtered   = sum(Filtered),
    pct_retained = round(100 * n_valid / n_windows, 1),
    .groups = "drop"
  ) %>%
  print()
cat("\n")

# --- 2.6 Summary statistics ---

cat("--- 2.6 Daily summary (all levels) ---\n\n")

ammod_daily_summary <- ammod_daily_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_days      = dplyr::n(),
    mean_DS     = round(mean(Delta_S)   * 100, 1),
    sd_DS       = round(sd(Delta_S)     * 100, 1),
    median_DS   = round(median(Delta_S) * 100, 1),
    min_DS      = round(min(Delta_S)    * 100, 1),
    max_DS      = round(max(Delta_S)    * 100, 1),
    mean_taxa_A = round(mean(N_taxa_A), 1),
    mean_taxa_B = round(mean(N_taxa_B), 1),
    mean_shared = round(mean(N_shared), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Pair,
                 factor(Tax_Level, levels = c("Species", "Genus", "Superfamily")))
print(ammod_daily_summary, n = Inf, width = 200)

cat("\n--- Windowed summary (Superfamily) ---\n\n")

ammod_windowed_summary <- ammod_windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    n_windows   = dplyr::n(),
    mean_DS     = round(mean(Delta_S)   * 100, 1),
    sd_DS       = round(sd(Delta_S)     * 100, 1),
    median_DS   = round(median(Delta_S) * 100, 1),
    min_DS      = round(min(Delta_S)    * 100, 1),
    max_DS      = round(max(Delta_S)    * 100, 1),
    mean_taxa_A = round(mean(N_taxa_A), 1),
    mean_taxa_B = round(mean(N_taxa_B), 1),
    mean_shared = round(mean(N_shared), 1),
    .groups = "drop"
  )
print(ammod_windowed_summary, n = Inf, width = 200)
cat("\n")

# --- 2.7 Figure 14: Daily Species temporal ---

cat("--- 2.7 Figure 14 ---\n")

f14_data <- ammod_daily_df %>%
  dplyr::filter(Tax_Level == "Species", !is.na(Delta_S)) %>%
  dplyr::mutate(DS_pct = Delta_S * 100)

f14_stats <- f14_data %>%
  dplyr::group_by(Habitat) %>%
  dplyr::summarise(mean_ds = round(mean(DS_pct), 1),
                   sd_ds   = round(sd(DS_pct),   1),
                   n       = dplyr::n(), .groups = "drop")

subtitle_f14 <- f14_stats %>%
  dplyr::mutate(lab = sprintf("%s: %.1f%% ± %.1f%% (n=%d)",
                              Habitat, mean_ds, sd_ds, n)) %>%
  dplyr::pull(lab) %>% paste(collapse = " | ")

fig14 <- ggplot(f14_data, aes(x = Date, y = DS_pct, color = Habitat)) +
  geom_line(aes(group = Pair), linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = custom_colors,
    labels = c("Maize"  = "AMMOD1 vs AMMOD2 (Maize)",
               "Meadow" = "AMMOD3 vs AMMOD4 (Meadow)")
  ) +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  labs(
    title    = "H1 Qualitative: AMMOD Intra-Habitat Consistency — Species Level",
    subtitle = subtitle_f14,
    x        = "Date",
    y        = "Izsak-Price Similarity (%)",
    color    = "Pair",
    caption  = "Dashed line = 50% similarity reference."
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

ggsave("outputs/figures/Fig14_H1_quali_AMMOD_intra_IP_daily_Species.png",
       plot = fig14, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/figures/Fig14_H1_quali_AMMOD_intra_IP_daily_Species.png\n\n")

# --- 2.8 Figure 15: Windowed Superfamily temporal ---

cat("--- 2.8 Figure 15 ---\n")

f15_data <- ammod_windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::mutate(
    DS_pct     = Delta_S * 100,
    Window_Num = as.numeric(gsub("W(\\d+).*", "\\1", Period))
  )

f15_stats <- f15_data %>%
  dplyr::group_by(Habitat) %>%
  dplyr::summarise(mean_ds = round(mean(DS_pct), 1),
                   sd_ds   = round(sd(DS_pct),   1),
                   n       = dplyr::n(), .groups = "drop")

subtitle_f15 <- f15_stats %>%
  dplyr::mutate(lab = sprintf("%s: %.1f%% ± %.1f%% (n=%d)",
                              Habitat, mean_ds, sd_ds, n)) %>%
  dplyr::pull(lab) %>% paste(collapse = " | ")

win_labs_f15 <- f15_data %>%
  dplyr::distinct(Window_Num, Period) %>%
  dplyr::arrange(Window_Num) %>%
  dplyr::mutate(label = gsub("W\\d+ \\((.*)\\)", "\\1", Period))

fig15 <- ggplot(f15_data, aes(x = Window_Num, y = DS_pct, color = Habitat)) +
  geom_line(aes(group = Pair), linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_x_continuous(breaks = win_labs_f15$Window_Num,
                     labels = win_labs_f15$label) +
  scale_color_manual(
    values = custom_colors,
    labels = c("Maize"  = "AMMOD1 vs AMMOD2 (Maize)",
               "Meadow" = "AMMOD3 vs AMMOD4 (Meadow)")
  ) +
  labs(
    title    = "H1 Qualitative: AMMOD Intra-Habitat Consistency — Superfamily Level",
    subtitle = subtitle_f15,
    x        = "3-Day Window",
    y        = "Izsak-Price Similarity (%)",
    color    = "Pair",
    caption  = paste0("Dashed line = 50% reference. Filtered: >= ",
                      MIN_OTUS, " OTUs per device per window.")
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

ggsave("outputs/figures/Fig15_H1_quali_AMMOD_intra_IP_windowed_Superfamily.png",
       plot = fig15, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/figures/Fig15_H1_quali_AMMOD_intra_IP_windowed_Superfamily.png\n\n")


# ==============================================================================
# SECTION 3: FAIRD INTER-HABITAT — MORISITA-HORN (§3.2.2.2)
# ==============================================================================
# Source: H1_consistency_FAIRD_taxonomic.R, Section 8.5
# Comparison: Pool(FAIRD1+FAIRD2) vs Pool(FAIRD3+FAIRD4)
# Reuses compute_mh() via Device column rename trick from source script
# ==============================================================================

cat("=== SECTION 3: FAIRD INTER-HABITAT — MORISITA-HORN (§3.2.2.2) ===\n\n")

faird_pooled <- faird_taxo %>%
  dplyr::mutate(
    Pool = dplyr::case_when(
      Device %in% c("FAIRD1", "FAIRD2") ~ "Maize_Pool",
      Device %in% c("FAIRD3", "FAIRD4") ~ "Meadow_Pool",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Pool)) %>%
  dplyr::select(-Device) %>%
  dplyr::rename(Device = Pool)

cat("Processing: Maize_Pool vs Meadow_Pool | Superfamily | 3-day windows\n\n")

faird_inter <- list()

for (w in unique(stats::na.omit(faird_pooled$Window_ID))) {
  win_data  <- faird_pooled %>% dplyr::filter(Window_ID == w)
  win_label <- unique(win_data$Window_Label)[1]
  counts_a  <- get_counts_mh(win_data, "Maize_Pool",  "Taxon_Superfamily")
  counts_b  <- get_counts_mh(win_data, "Meadow_Pool", "Taxon_Superfamily")

  passes <- counts_a$n_indiv >= MIN_INDIV & counts_b$n_indiv >= MIN_INDIV
  mh_sim <- if (passes)
    compute_mh(win_data, "Maize_Pool", "Meadow_Pool", "Taxon_Superfamily") else NA_real_

  faird_inter[[length(faird_inter) + 1]] <- tibble::tibble(
    Resolution = "3-day window",
    Pair       = "Maize_Pool vs Meadow_Pool",
    Habitat    = "Inter-Habitat",
    Period     = win_label,
    Tax_Level  = "Superfamily",
    MH_Sim     = mh_sim,
    N_indiv_A  = counts_a$n_indiv,
    N_indiv_B  = counts_b$n_indiv,
    N_taxa_A   = counts_a$n_taxa,
    N_taxa_B   = counts_b$n_taxa,
    Filtered   = !passes
  )
}

faird_inter_df <- dplyr::bind_rows(faird_inter)
cat(sprintf("Windows: %d | Valid: %d | Filtered: %d\n\n",
            nrow(faird_inter_df),
            sum(!is.na(faird_inter_df$MH_Sim)),
            sum(faird_inter_df$Filtered)))

faird_inter_summary <- faird_inter_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::summarise(
    n_windows         = dplyr::n(),
    mean_MH           = round(mean(MH_Sim)   * 100, 1),
    sd_MH             = round(sd(MH_Sim)     * 100, 1),
    median_MH         = round(median(MH_Sim) * 100, 1),
    min_MH            = round(min(MH_Sim)    * 100, 1),
    max_MH            = round(max(MH_Sim)    * 100, 1),
    mean_indiv_Maize  = round(mean(N_indiv_A), 1),
    mean_indiv_Meadow = round(mean(N_indiv_B), 1)
  )

cat("--- FAIRD Inter-Habitat Summary (MH, %) ---\n\n")
print(faird_inter_summary, width = 200)
cat("\n")


# ==============================================================================
# SECTION 4: AMMOD INTER-HABITAT — IZSAK-PRICE (§3.2.2.2)
# ==============================================================================
# Source: H1_consistency_AMMOD_taxonomic.R, Sections 5.5 + 6.5
# Comparison: Pool(AMMOD1+AMMOD2) vs Pool(AMMOD3+AMMOD4)
# ==============================================================================

cat("=== SECTION 4: AMMOD INTER-HABITAT — IZSAK-PRICE (§3.2.2.2) ===\n\n")

# --- 4.1 Daily IP (pooled habitats, all levels) ---

cat("--- 4.1 Daily IP (pooled, all levels) ---\n\n")
cat("Processing: Maize_Pool (AMMOD1+2) vs Meadow_Pool (AMMOD3+4)\n\n")

ammod_inter_daily <- list()

for (lev in ammod_tax_levels) {
  for (d in all_dates_ammod) {
    date_val <- as.Date(d, origin = "1970-01-01")

    taxa_maize <- ammod_raw %>%
      dplyr::filter(Date == date_val, Device %in% c("AMMOD1", "AMMOD2"),
                    !is.na(.data[[lev$col]])) %>%
      dplyr::pull(.data[[lev$col]]) %>% unique()

    taxa_meadow <- ammod_raw %>%
      dplyr::filter(Date == date_val, Device %in% c("AMMOD3", "AMMOD4"),
                    !is.na(.data[[lev$col]])) %>%
      dplyr::pull(.data[[lev$col]]) %>% unique()

    delta_s <- if (length(taxa_maize) > 0 & length(taxa_meadow) > 0)
      calculate_ip(taxa_maize, taxa_meadow, lev$dmat) else NA_real_

    ammod_inter_daily[[length(ammod_inter_daily) + 1]] <- tibble::tibble(
      Resolution = "Daily",
      Pair       = "Maize_Pool vs Meadow_Pool",
      Habitat    = "Inter-Habitat",
      Date       = date_val,
      Period     = as.character(date_val),
      Tax_Level  = lev$label,
      Delta_S    = delta_s,
      N_taxa_A   = length(taxa_maize),
      N_taxa_B   = length(taxa_meadow),
      N_shared   = length(intersect(taxa_maize, taxa_meadow))
    )
  }
}

ammod_inter_daily_df <- dplyr::bind_rows(ammod_inter_daily)
cat(sprintf("Daily inter-habitat IP: %d rows | Valid: %d\n\n",
            nrow(ammod_inter_daily_df),
            sum(!is.na(ammod_inter_daily_df$Delta_S))))

cat("--- Inter-Habitat Daily Summary ---\n\n")
ammod_inter_daily_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Tax_Level) %>%
  dplyr::summarise(
    n_days      = dplyr::n(),
    mean_DS     = round(mean(Delta_S)   * 100, 1),
    sd_DS       = round(sd(Delta_S)     * 100, 1),
    median_DS   = round(median(Delta_S) * 100, 1),
    mean_taxa_A = round(mean(N_taxa_A), 1),
    mean_taxa_B = round(mean(N_taxa_B), 1),
    mean_shared = round(mean(N_shared), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Species", "Genus", "Superfamily"))) %>%
  print(n = Inf, width = 200)
cat("\n")

# --- 4.2 Windowed IP (Superfamily, pooled habitats) ---

cat(sprintf("--- 4.2 Windowed IP (Superfamily, pooled, >= %d OTUs) ---\n\n", MIN_OTUS))

ammod_inter_windowed <- list()

for (w in unique(window_lookup$Window_ID)) {
  win_dates <- window_lookup %>%
    dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
  win_dates <- win_dates[win_dates != date_exclude]
  win_label <- window_lookup %>%
    dplyr::filter(Window_ID == w) %>% dplyr::pull(Window_Label) %>% unique()

  sf_maize <- ammod_raw %>%
    dplyr::filter(Date %in% win_dates, Device %in% c("AMMOD1", "AMMOD2"),
                  !is.na(Taxon_SF)) %>%
    dplyr::pull(Taxon_SF) %>% unique()
  sf_meadow <- ammod_raw %>%
    dplyr::filter(Date %in% win_dates, Device %in% c("AMMOD3", "AMMOD4"),
                  !is.na(Taxon_SF)) %>%
    dplyr::pull(Taxon_SF) %>% unique()

  passes  <- length(sf_maize) >= MIN_OTUS & length(sf_meadow) >= MIN_OTUS
  delta_s <- if (passes) calculate_ip(sf_maize, sf_meadow, dist_sf) else NA_real_

  ammod_inter_windowed[[length(ammod_inter_windowed) + 1]] <- tibble::tibble(
    Resolution = "3-day window",
    Pair       = "Maize_Pool vs Meadow_Pool",
    Habitat    = "Inter-Habitat",
    Period     = win_label,
    Window_ID  = w,
    Tax_Level  = "Superfamily",
    Delta_S    = delta_s,
    N_taxa_A   = length(sf_maize),
    N_taxa_B   = length(sf_meadow),
    N_shared   = length(intersect(sf_maize, sf_meadow)),
    Filtered   = !passes
  )
}

ammod_inter_windowed_df <- dplyr::bind_rows(ammod_inter_windowed)
cat(sprintf("Windowed inter-habitat IP: %d rows | Valid: %d | Filtered: %d\n\n",
            nrow(ammod_inter_windowed_df),
            sum(!is.na(ammod_inter_windowed_df$Delta_S)),
            sum(ammod_inter_windowed_df$Filtered)))

ammod_inter_summary <- ammod_inter_windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::summarise(
    n_windows        = dplyr::n(),
    mean_DS          = round(mean(Delta_S)   * 100, 1),
    sd_DS            = round(sd(Delta_S)     * 100, 1),
    median_DS        = round(median(Delta_S) * 100, 1),
    min_DS           = round(min(Delta_S)    * 100, 1),
    max_DS           = round(max(Delta_S)    * 100, 1),
    mean_taxa_Maize  = round(mean(N_taxa_A), 1),
    mean_taxa_Meadow = round(mean(N_taxa_B), 1),
    mean_shared      = round(mean(N_shared), 1)
  )

cat("--- AMMOD Inter-Habitat Summary (IP, %) ---\n\n")
print(ammod_inter_summary, width = 200)
cat("\n")


# ==============================================================================
# SECTION 5: INTRA vs INTER-HABITAT STATISTICAL TESTS
# ==============================================================================
# Source: H1_H2_H3_taxonomic_size_analysis.R, lines 366-450
# Wilcoxon rank-sum + Cohen's d on windowed values (same metric within system)
# ==============================================================================

cat("=== SECTION 5: INTRA vs INTER-HABITAT STATISTICAL TESTS ===\n\n")

# --- 5.1 FAIRD: Intra vs Inter (MH, windowed Superfamily) ---

cat("--- 5.1 FAIRD: Intra vs Inter-Habitat (MH, windowed) ---\n\n")

faird_intra_vals <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily") %>%
  dplyr::pull(MH_Sim)
faird_inter_vals <- faird_inter_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::pull(MH_Sim)

if (length(faird_intra_vals) >= 3 && length(faird_inter_vals) >= 3) {
  wt_f <- stats::wilcox.test(faird_intra_vals, faird_inter_vals, exact = FALSE)
  cd_f <- effsize::cohen.d(faird_intra_vals, faird_inter_vals)
  cat(sprintf("FAIRD Intra: n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(faird_intra_vals),
              mean(faird_intra_vals) * 100, median(faird_intra_vals) * 100))
  cat(sprintf("FAIRD Inter: n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(faird_inter_vals),
              mean(faird_inter_vals) * 100, median(faird_inter_vals) * 100))
  cat(sprintf("Wilcoxon W = %.1f | p = %.4f%s\n",
              wt_f$statistic, wt_f$p.value,
              ifelse(wt_f$p.value < 0.05, " *", " ns")))
  cat(sprintf("Cohen's d = %.3f (%s)\n\n", cd_f$estimate, cd_f$magnitude))
}

# --- 5.2 AMMOD: Intra vs Inter (IP, windowed Superfamily) ---

cat("--- 5.2 AMMOD: Intra vs Inter-Habitat (IP, windowed) ---\n\n")

ammod_intra_vals <- ammod_windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>% dplyr::pull(Delta_S)
ammod_inter_vals <- ammod_inter_windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>% dplyr::pull(Delta_S)

if (length(ammod_intra_vals) >= 3 && length(ammod_inter_vals) >= 3) {
  wt_a <- stats::wilcox.test(ammod_intra_vals, ammod_inter_vals, exact = FALSE)
  cd_a <- effsize::cohen.d(ammod_intra_vals, ammod_inter_vals)
  cat(sprintf("AMMOD Intra: n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(ammod_intra_vals),
              mean(ammod_intra_vals) * 100, median(ammod_intra_vals) * 100))
  cat(sprintf("AMMOD Inter: n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(ammod_inter_vals),
              mean(ammod_inter_vals) * 100, median(ammod_inter_vals) * 100))
  cat(sprintf("Wilcoxon W = %.1f | p = %.4f%s\n",
              wt_a$statistic, wt_a$p.value,
              ifelse(wt_a$p.value < 0.05, " *", " ns")))
  cat(sprintf("Cohen's d = %.3f (%s)\n\n", cd_a$estimate, cd_a$magnitude))
}

# --- 5.3 FAIRD: Intra-Maize vs Intra-Meadow (MH, windowed Superfamily) ---

cat("--- 5.3 FAIRD: Intra-Maize vs Intra-Meadow (MH, windowed) ---\n\n")

faird_maize_w  <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily", Habitat == "Maize") %>%
  dplyr::pull(MH_Sim)
faird_meadow_w <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily", Habitat == "Meadow") %>%
  dplyr::pull(MH_Sim)

if (length(faird_maize_w) >= 3 && length(faird_meadow_w) >= 3) {
  wt_fm <- stats::wilcox.test(faird_maize_w, faird_meadow_w, exact = FALSE)
  cd_fm <- effsize::cohen.d(faird_maize_w, faird_meadow_w)
  cat(sprintf("FAIRD1 vs FAIRD2 (Maize):  n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(faird_maize_w),
              mean(faird_maize_w) * 100, median(faird_maize_w) * 100))
  cat(sprintf("FAIRD3 vs FAIRD4 (Meadow): n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(faird_meadow_w),
              mean(faird_meadow_w) * 100, median(faird_meadow_w) * 100))
  cat(sprintf("Wilcoxon W = %.1f | p = %.4f%s\n",
              wt_fm$statistic, wt_fm$p.value,
              ifelse(wt_fm$p.value < 0.05, " *", " ns")))
  cat(sprintf("Cohen's d = %.3f (%s)\n\n", cd_fm$estimate, cd_fm$magnitude))
}

# --- 5.4 AMMOD: Intra-Maize vs Intra-Meadow (IP, windowed Superfamily) ---

cat("--- 5.4 AMMOD: Intra-Maize vs Intra-Meadow (IP, windowed) ---\n\n")

ammod_maize_w  <- ammod_windowed_df %>%
  dplyr::filter(!is.na(Delta_S), Habitat == "Maize") %>% dplyr::pull(Delta_S)
ammod_meadow_w <- ammod_windowed_df %>%
  dplyr::filter(!is.na(Delta_S), Habitat == "Meadow") %>% dplyr::pull(Delta_S)

if (length(ammod_maize_w) >= 3 && length(ammod_meadow_w) >= 3) {
  wt_am <- stats::wilcox.test(ammod_maize_w, ammod_meadow_w, exact = FALSE)
  cd_am <- effsize::cohen.d(ammod_maize_w, ammod_meadow_w)
  cat(sprintf("AMMOD1 vs AMMOD2 (Maize):  n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(ammod_maize_w),
              mean(ammod_maize_w) * 100, median(ammod_maize_w) * 100))
  cat(sprintf("AMMOD3 vs AMMOD4 (Meadow): n=%d | mean=%.1f%% | median=%.1f%%\n",
              length(ammod_meadow_w),
              mean(ammod_meadow_w) * 100, median(ammod_meadow_w) * 100))
  cat(sprintf("Wilcoxon W = %.1f | p = %.4f%s\n",
              wt_am$statistic, wt_am$p.value,
              ifelse(wt_am$p.value < 0.05, " *", " ns")))
  cat(sprintf("Cohen's d = %.3f (%s)\n\n", cd_am$estimate, cd_am$magnitude))
}


# ==============================================================================
# SECTION 6: CONSOLIDATED TABLES (Tables 11-12)
# ==============================================================================

cat("=== SECTION 6: CONSOLIDATED TABLES ===\n\n")

# --- Table 11: All pairs — windowed Superfamily summary ---

cat("--- Table 11: Windowed Superfamily similarity (all pairs, all levels) ---\n\n")

t11_faird_intra <- faird_windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily") %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    System     = "FAIRD",
    Level      = "Intra-Habitat",
    Metric     = "Morisita-Horn",
    n_windows  = dplyr::n(),
    mean_sim   = round(mean(MH_Sim)   * 100, 1),
    sd_sim     = round(sd(MH_Sim)     * 100, 1),
    median_sim = round(median(MH_Sim) * 100, 1),
    min_sim    = round(min(MH_Sim)    * 100, 1),
    max_sim    = round(max(MH_Sim)    * 100, 1),
    .groups = "drop"
  )

t11_faird_inter <- faird_inter_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::summarise(
    Pair       = "Maize_Pool vs Meadow_Pool",
    Habitat    = "Inter-Habitat",
    System     = "FAIRD",
    Level      = "Inter-Habitat",
    Metric     = "Morisita-Horn",
    n_windows  = dplyr::n(),
    mean_sim   = round(mean(MH_Sim)   * 100, 1),
    sd_sim     = round(sd(MH_Sim)     * 100, 1),
    median_sim = round(median(MH_Sim) * 100, 1),
    min_sim    = round(min(MH_Sim)    * 100, 1),
    max_sim    = round(max(MH_Sim)    * 100, 1)
  )

t11_ammod_intra <- ammod_windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    System     = "AMMOD",
    Level      = "Intra-Habitat",
    Metric     = "Izsak-Price",
    n_windows  = dplyr::n(),
    mean_sim   = round(mean(Delta_S)   * 100, 1),
    sd_sim     = round(sd(Delta_S)     * 100, 1),
    median_sim = round(median(Delta_S) * 100, 1),
    min_sim    = round(min(Delta_S)    * 100, 1),
    max_sim    = round(max(Delta_S)    * 100, 1),
    .groups = "drop"
  )

t11_ammod_inter <- ammod_inter_windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::summarise(
    Pair       = "Maize_Pool vs Meadow_Pool",
    Habitat    = "Inter-Habitat",
    System     = "AMMOD",
    Level      = "Inter-Habitat",
    Metric     = "Izsak-Price",
    n_windows  = dplyr::n(),
    mean_sim   = round(mean(Delta_S)   * 100, 1),
    sd_sim     = round(sd(Delta_S)     * 100, 1),
    median_sim = round(median(Delta_S) * 100, 1),
    min_sim    = round(min(Delta_S)    * 100, 1),
    max_sim    = round(max(Delta_S)    * 100, 1)
  )

table_11 <- dplyr::bind_rows(t11_faird_intra, t11_faird_inter,
                              t11_ammod_intra, t11_ammod_inter) %>%
  dplyr::select(System, Level, Pair, Habitat, Metric,
                n_windows, mean_sim, sd_sim, median_sim, min_sim, max_sim) %>%
  dplyr::arrange(System, Level, Habitat)

print(table_11, n = Inf, width = 200)
readr::write_csv(table_11, "outputs/tables/Table11_H1_quali_taxonomic_similarity.csv")
cat("\n[OK] outputs/tables/Table11_H1_quali_taxonomic_similarity.csv\n\n")

# --- Table 12: AMMOD detailed — daily IP all levels, intra + inter ---

cat("--- Table 12: AMMOD daily IP all levels (intra + inter) ---\n\n")

table_12 <- dplyr::bind_rows(
  ammod_daily_df       %>% dplyr::mutate(Scope = "Intra-Habitat"),
  ammod_inter_daily_df %>% dplyr::mutate(Scope = "Inter-Habitat")
) %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Scope, Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_days      = dplyr::n(),
    mean_DS     = round(mean(Delta_S)   * 100, 1),
    sd_DS       = round(sd(Delta_S)     * 100, 1),
    median_DS   = round(median(Delta_S) * 100, 1),
    min_DS      = round(min(Delta_S)    * 100, 1),
    max_DS      = round(max(Delta_S)    * 100, 1),
    mean_taxa_A = round(mean(N_taxa_A), 1),
    mean_taxa_B = round(mean(N_taxa_B), 1),
    mean_shared = round(mean(N_shared), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Scope, Pair,
                 factor(Tax_Level, levels = c("Species", "Genus", "Superfamily")))

print(table_12, n = Inf, width = 200)
readr::write_csv(table_12, "outputs/tables/Table12_H1_quali_AMMOD_IP_all_levels.csv")
cat("\n[OK] outputs/tables/Table12_H1_quali_AMMOD_IP_all_levels.csv\n\n")


# ==============================================================================
# WRAP UP
# ==============================================================================

cat("==============================================================================\n")
cat("Script completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Outputs:\n")
cat("  Fig 12 — FAIRD intra-habitat MH windowed Superfamily\n")
cat("  Fig 13 — FAIRD intra-habitat MH windowed by level\n")
cat("  Fig 14 — AMMOD intra-habitat IP daily Species\n")
cat("  Fig 15 — AMMOD intra-habitat IP windowed Superfamily\n")
cat("  Table 11 — outputs/tables/Table11_H1_quali_taxonomic_similarity.csv\n")
cat("  Table 12 — outputs/tables/Table12_H1_quali_AMMOD_IP_all_levels.csv\n")
cat("  Console  — outputs/console/H1_quali_consistency_output.txt\n")
cat("==============================================================================\n")

while (sink.number() > 0) sink()
