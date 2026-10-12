# ==============================================================================
# H2: INTER-SYSTEM QUALITATIVE VALIDATION — SUPERFAMILY LEVEL
# ==============================================================================
# Companion to H3_qualitative_FAIRD_AMMOD_size.R (Genus level)
# Computes Izsak-Price qualitative similarity at SUPERFAMILY level
# with daily + 3-day windowed (filtered >= 3 FAIRD individuals) analyses.
#
# Rationale: Superfamily (~68% FAIRD coverage) reduces artefacts from
# variable classification depth across devices/observers while preserving
# ecologically meaningful groupings (Muscoidea, Syrphoidea, Oestroidea).
#
# Author: Juan Chiavassa & Claude
# Date: March 2026
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(vegan)       # taxa2dist()
library(lubridate)
library(effsize)

# ==============================================================================
# CONSOLE OUTPUT CAPTURE
# ==============================================================================
script_name <- "H2_qualitative_FAIRD_AMMOD_superfamily"

output_file <- file.path("outputs",
                         paste0(format(Sys.time(), "%Y%m%d_%H%M%S"),
                                "_", script_name, ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)
while (sink.number() > 0) sink()
sink(output_file, split = TRUE)

cat("==============================================================================\n")
cat("H3: IZSAK-PRICE SIMILARITY AT SUPERFAMILY LEVEL\n")
cat("Daily + 3-day windows (filtered >= 3 FAIRD individuals)\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: LOAD AND PREPARE DATA
# ==============================================================================

cat("=== SECTION 1: LOAD DATA ===\n\n")

source("load_data.R")
source("custom_colors.R")

# --- Temporal scope (same as H3_qualitative_FAIRD_AMMOD_size.R) ---
date_start <- as.Date("2023-08-23")
date_end   <- as.Date("2023-09-13")
date_exclude <- as.Date("2023-09-04")
ammod4_failure_start <- as.Date("2023-08-28")
ammod4_failure_end   <- as.Date("2023-09-03")

nc <- "#N/C"

# --- FAIRD data: assign Superfamily ---
faird_data <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD",
                Date >= date_start & Date <= date_end) %>%
  dplyr::select(Date, Device, Site, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Genus) %>%
  dplyr::mutate(
    Taxon_SF = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc & Superfamily != "",
      Superfamily, NA_character_
    )
  )

cat("FAIRD loaded:", nrow(faird_data), "individuals\n")
cat("  With Superfamily:", sum(!is.na(faird_data$Taxon_SF)),
    sprintf("(%.1f%%)\n", 100 * sum(!is.na(faird_data$Taxon_SF)) / nrow(faird_data)))

# --- AMMOD data: collapse to Superfamily presence/absence ---
ammod_data <- data_list$ammod_taxonomy %>%
  dplyr::filter(
    Date >= date_start & Date <= date_end,
    Date != date_exclude,
    !(Device == "AMMOD4" & Date >= ammod4_failure_start & Date <= ammod4_failure_end)
  ) %>%
  dplyr::select(Date, Device, Site, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Genus) %>%
  dplyr::mutate(
    Taxon_SF = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc & Superfamily != "",
      Superfamily, NA_character_
    )
  )

cat("AMMOD loaded:", nrow(ammod_data), "OTUs\n")
cat("  With Superfamily:", sum(!is.na(ammod_data$Taxon_SF)),
    sprintf("(%.1f%%)\n\n", 100 * sum(!is.na(ammod_data$Taxon_SF)) / nrow(ammod_data)))


# ==============================================================================
# SECTION 2: BUILD TAXONOMIC DISTANCE MATRIX (SUPERFAMILY LEVEL)
# ==============================================================================

cat("=== SECTION 2: DISTANCE MATRIX (SUPERFAMILY) ===\n\n")

# Collect unique Superfamilies with their hierarchy
faird_sf_taxa <- faird_data %>%
  dplyr::filter(!is.na(Taxon_SF)) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  dplyr::distinct(Taxon_SF, .keep_all = TRUE)

ammod_sf_taxa <- ammod_data %>%
  dplyr::filter(!is.na(Taxon_SF)) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  dplyr::distinct(Taxon_SF, .keep_all = TRUE)

all_sf_taxa <- dplyr::bind_rows(faird_sf_taxa, ammod_sf_taxa) %>%
  dplyr::distinct(Taxon_SF, .keep_all = TRUE)

cat("Unique Superfamilies: FAIRD =", nrow(faird_sf_taxa),
    "| AMMOD =", nrow(ammod_sf_taxa),
    "| Combined =", nrow(all_sf_taxa), "\n")
cat("Shared:", length(intersect(faird_sf_taxa$Taxon_SF, ammod_sf_taxa$Taxon_SF)), "\n\n")

# Build hierarchy table for taxa2dist
tax_table_sf <- all_sf_taxa %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  as.data.frame()

rownames(tax_table_sf) <- tax_table_sf$Taxon_SF
tax_table_sf <- tax_table_sf %>% dplyr::select(-Taxon_SF)

taxon_distances_sf <- vegan::taxa2dist(tax_table_sf, varstep = TRUE, check = TRUE)
dist_matrix_sf <- as.matrix(taxon_distances_sf)

cat("Distance matrix:", nrow(dist_matrix_sf), "x", ncol(dist_matrix_sf), "\n")
cat("Distance range:", round(min(taxon_distances_sf), 4), "-",
    round(max(taxon_distances_sf), 4), "\n\n")


# ==============================================================================
# SECTION 3: IZSAK-PRICE FUNCTION (same as main H3 script)
# ==============================================================================

calculate_qualitative_similarity <- function(taxa_A, taxa_B, dist_matrix) {
  taxa_A <- taxa_A[taxa_A %in% rownames(dist_matrix)]
  taxa_B <- taxa_B[taxa_B %in% rownames(dist_matrix)]
  
  if (length(taxa_A) == 0 | length(taxa_B) == 0) return(NA)
  
  min_distances_A_to_B <- sapply(taxa_A, function(a) {
    if (a %in% taxa_B) 0
    else min(dist_matrix[a, taxa_B], na.rm = TRUE)
  })
  
  min_distances_B_to_A <- sapply(taxa_B, function(b) {
    if (b %in% taxa_A) 0
    else min(dist_matrix[b, taxa_A], na.rm = TRUE)
  })
  
  avg_min_dist <- mean(c(min_distances_A_to_B, min_distances_B_to_A))
  max_dist <- max(dist_matrix, na.rm = TRUE)
  
  return(1 - (avg_min_dist / max_dist))
}


# ==============================================================================
# SECTION 4: DAILY IZSAK-PRICE (SUPERFAMILY)
# ==============================================================================

cat("=== SECTION 4: DAILY SIMILARITY ===\n\n")

sites <- c("Site1", "Site2", "Site3", "Site4")
all_dates <- seq.Date(date_start, date_end, by = "day")

daily_results <- list()

for (s in sites) {
  for (d in all_dates) {
    date_val <- as.Date(d, origin = "1970-01-01")
    
    faird_sf <- faird_data %>%
      dplyr::filter(Date == date_val, Site == s, !is.na(Taxon_SF)) %>%
      dplyr::pull(Taxon_SF) %>% unique()
    
    ammod_sf <- ammod_data %>%
      dplyr::filter(Date == date_val, Site == s, !is.na(Taxon_SF)) %>%
      dplyr::pull(Taxon_SF) %>% unique()
    
    n_faird_indiv <- faird_data %>%
      dplyr::filter(Date == date_val, Site == s, !is.na(Taxon_SF)) %>% nrow()
    
    delta_s <- NA_real_
    if (length(faird_sf) > 0 & length(ammod_sf) > 0) {
      delta_s <- calculate_qualitative_similarity(faird_sf, ammod_sf, dist_matrix_sf)
    }
    
    daily_results[[length(daily_results) + 1]] <- tibble::tibble(
      Resolution    = "Daily",
      Date          = date_val,
      Site          = s,
      Habitat       = dplyr::if_else(s %in% c("Site1", "Site2"), "Maize", "Meadow"),
      Delta_S       = delta_s,
      N_taxa_FAIRD  = length(faird_sf),
      N_taxa_AMMOD  = length(ammod_sf),
      N_indiv_FAIRD = n_faird_indiv,
      N_shared      = length(intersect(faird_sf, ammod_sf))
    )
  }
}

daily_df <- dplyr::bind_rows(daily_results)
cat("Daily computations:", nrow(daily_df), "\n")
cat("Valid Delta_S:", sum(!is.na(daily_df$Delta_S)), "\n\n")


# ==============================================================================
# SECTION 5: 3-DAY WINDOWED IZSAK-PRICE (FILTERED >= 3 FAIRD INDIVIDUALS)
# ==============================================================================

cat("=== SECTION 5: 3-DAY WINDOWED SIMILARITY ===\n\n")

MIN_INDIV <- 3
cat(sprintf("Minimum FAIRD individuals per window: %d\n\n", MIN_INDIV))

# Build windows
window_breaks <- seq(date_start, date_end + 2, by = 3)
window_lookup <- tibble::tibble(Date = all_dates) %>%
  dplyr::mutate(
    Window_ID = as.integer(cut(Date, breaks = window_breaks, right = FALSE,
                               labels = FALSE)),
    Window_Start = window_breaks[Window_ID],
    Window_End   = pmin(Window_Start + 2, date_end),
    Window_Label = paste0("W", Window_ID, " (",
                          format(Window_Start, "%b %d"), "-",
                          format(Window_End, "%b %d"), ")")
  )

windowed_results <- list()

for (s in sites) {
  for (w in unique(na.omit(window_lookup$Window_ID))) {
    win_dates <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
    win_label <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>% dplyr::pull(Window_Label) %>%
      unique()
    
    # Pool FAIRD across window (presence/absence of Superfamilies)
    faird_win <- faird_data %>%
      dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF))
    faird_sf <- unique(faird_win$Taxon_SF)
    n_faird_indiv <- nrow(faird_win)
    
    # Pool AMMOD across window (already presence/absence)
    ammod_win <- ammod_data %>%
      dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF))
    ammod_sf <- unique(ammod_win$Taxon_SF)
    
    passes_filter <- (n_faird_indiv >= MIN_INDIV & length(ammod_sf) > 0)
    
    delta_s <- NA_real_
    if (passes_filter & length(faird_sf) > 0) {
      delta_s <- calculate_qualitative_similarity(faird_sf, ammod_sf, dist_matrix_sf)
    }
    
    windowed_results[[length(windowed_results) + 1]] <- tibble::tibble(
      Resolution    = "3-day window",
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


# ==============================================================================
# SECTION 6: SUMMARY STATISTICS
# ==============================================================================

cat("=== SECTION 6: SUMMARY STATISTICS ===\n\n")

# --- 6a: Daily summary by site ---
cat("--- 6a: Daily summary (Superfamily, unfiltered) ---\n\n")

daily_by_site <- daily_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n_days       = dplyr::n(),
    mean_DS      = round(mean(Delta_S) * 100, 1),
    sd_DS        = round(sd(Delta_S) * 100, 1),
    faird_sf_day = round(mean(N_taxa_FAIRD), 1),
    ammod_sf_day = round(mean(N_taxa_AMMOD), 1),
    exact_match  = round(mean(N_shared), 1),
    .groups = "drop"
  )
print(daily_by_site, n = Inf)

# Grand mean
daily_grand <- daily_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::summarise(
    Site = "All", Habitat = "--",
    n_days = dplyr::n(),
    mean_DS = round(mean(Delta_S) * 100, 1),
    sd_DS = round(sd(Delta_S) * 100, 1),
    faird_sf_day = round(mean(N_taxa_FAIRD), 1),
    ammod_sf_day = round(mean(N_taxa_AMMOD), 1),
    exact_match = round(mean(N_shared), 1)
  )
print(daily_grand)

# --- 6b: Windowed summary by site ---
cat("\n--- 6b: 3-day windowed summary (filtered >= 3) ---\n\n")

win_by_site <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n_windows    = dplyr::n(),
    mean_DS      = round(mean(Delta_S) * 100, 1),
    sd_DS        = round(sd(Delta_S) * 100, 1),
    faird_sf_win = round(mean(N_taxa_FAIRD), 1),
    ammod_sf_win = round(mean(N_taxa_AMMOD), 1),
    exact_match  = round(mean(N_shared), 1),
    .groups = "drop"
  )
print(win_by_site, n = Inf)

win_grand <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::summarise(
    Site = "All", Habitat = "--",
    n_windows = dplyr::n(),
    mean_DS = round(mean(Delta_S) * 100, 1),
    sd_DS = round(sd(Delta_S) * 100, 1),
    faird_sf_win = round(mean(N_taxa_FAIRD), 1),
    ammod_sf_win = round(mean(N_taxa_AMMOD), 1),
    exact_match = round(mean(N_shared), 1)
  )
print(win_grand)

# --- Intra-Habitat and Overall mean individuals per window ---
cat("\n--- FAIRD individuals per site-window (for Table 15) ---\n\n")

indiv_by_habitat <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Habitat) %>%
  dplyr::summarise(
    n_site_windows    = dplyr::n(),
    mean_indiv_FAIRD  = round(mean(N_indiv_FAIRD), 1),
    mean_SF_FAIRD     = round(mean(N_taxa_FAIRD), 1),
    mean_SF_AMMOD     = round(mean(N_taxa_AMMOD), 1),
    mean_shared       = round(mean(N_shared), 1),
    .groups = "drop"
  )
print(indiv_by_habitat)

indiv_overall <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::summarise(
    n_site_windows    = dplyr::n(),
    mean_indiv_FAIRD  = round(mean(N_indiv_FAIRD), 1),
    mean_SF_FAIRD     = round(mean(N_taxa_FAIRD), 1),
    mean_SF_AMMOD     = round(mean(N_taxa_AMMOD), 1),
    mean_shared       = round(mean(N_shared), 1)
  )
cat("\nOverall:\n")
print(indiv_overall)

# --- 6c: Habitat comparison (windowed) ---
cat("\n--- 6c: Habitat comparison (windowed) ---\n\n")

maize_w <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S), Habitat == "Maize") %>%
  dplyr::pull(Delta_S)
meadow_w <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S), Habitat == "Meadow") %>%
  dplyr::pull(Delta_S)

if (length(maize_w) >= 3 & length(meadow_w) >= 3) {
  gd <- data.frame(val = c(maize_w, meadow_w),
                   grp = factor(c(rep("Maize", length(maize_w)),
                                  rep("Meadow", length(meadow_w)))))
  lev_p <- car::leveneTest(val ~ grp, data = gd)$`Pr(>F)`[1]
  tt <- t.test(maize_w, meadow_w, var.equal = (lev_p >= 0.05))
  cd <- effsize::cohen.d(maize_w, meadow_w)
  
  cat(sprintf("Maize:  %.1f%% +/- %.1f%% (n=%d)\n",
              mean(maize_w)*100, sd(maize_w)*100, length(maize_w)))
  cat(sprintf("Meadow: %.1f%% +/- %.1f%% (n=%d)\n",
              mean(meadow_w)*100, sd(meadow_w)*100, length(meadow_w)))
  cat(sprintf("t = %.3f, p = %.4f, Cohen's d = %.3f (%s)\n\n",
              tt$statistic, tt$p.value, cd$estimate, cd$magnitude))
}

# --- 6d: Intra-Habitat and Overall (windowed) ---
cat("--- 6d: Intra-Habitat and Overall (windowed) ---\n\n")

for (hab in c("Maize", "Meadow")) {
  vals <- windowed_df %>%
    dplyr::filter(!is.na(Delta_S), Habitat == hab) %>%
    dplyr::pull(Delta_S)
  cat(sprintf("Intra-Habitat %s: %.1f%% +/- %.1f%% (n=%d)\n",
              hab, mean(vals)*100, sd(vals)*100, length(vals)))
}

all_vals <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::pull(Delta_S)
cat(sprintf("Overall: %.1f%% +/- %.1f%% (n=%d, 95%% CI: %.1f-%.1f%%)\n\n",
            mean(all_vals)*100, sd(all_vals)*100, length(all_vals),
            (mean(all_vals) - 1.96*sd(all_vals)/sqrt(length(all_vals)))*100,
            (mean(all_vals) + 1.96*sd(all_vals)/sqrt(length(all_vals)))*100))


# ==============================================================================
# SECTION 7: FIGURES
# ==============================================================================

cat("=== SECTION 7: FIGURES ===\n\n")

# --- 7a: Windowed temporal by site ---
cat("Creating windowed temporal plot...\n")

plot_win <- windowed_df %>%
  dplyr::filter(!Filtered) %>%
  dplyr::mutate(DS_pct = Delta_S * 100)

p_win_temporal <- ggplot(plot_win,
                         aes(x = Window_ID, y = DS_pct, color = Site)) +
  geom_line(aes(group = Site), linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = sort(unique(plot_win$Window_ID)),
                     labels = function(x) paste0("W", x)) +
  scale_color_manual(values = custom_colors,
                     labels = c("Site1" = "Site 1 (Maize)",
                                "Site2" = "Site 2 (Maize)",
                                "Site3" = "Site 3 (Meadow)",
                                "Site4" = "Site 4 (Meadow)")) +
  labs(
    title    = "H2 Qualitative - Taxonomic Similarity: FAIRD vs AMMOD  — Site Level",
    subtitle = paste0("Izsak-Price index — 3-Day Windows — Superfamily level filtered: >= ", MIN_INDIV, " FAIRD individuals"),
    x = "Window", y = "Izsak-Price Similarity (%)",
    color = "Site",
    caption = "Dashed line = 50% similarity reference."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(p_win_temporal)
ggsave("outputs/H3_superfamily_windowed_site.png",
       plot = p_win_temporal, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/H3_superfamily_windowed_site.png\n\n")

# --- 7b: Habitat boxplot (windowed) ---
cat("Creating habitat boxplot...\n")

plot_box <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::mutate(DS_pct = Delta_S * 100)

p_habitat_box <- ggplot(plot_box, aes(x = Habitat, y = DS_pct, fill = Habitat)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "H3: FAIRD-AMMOD Similarity by Habitat (Superfamily, 3-Day Windows)",
    x = "Habitat", y = "Izsak-Price Similarity (%)", fill = "Habitat"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p_habitat_box)
ggsave("outputs/H3_superfamily_windowed_habitat_boxplot.png",
       plot = p_habitat_box, width = 7, height = 6, dpi = 300)
cat("[OK] outputs/H3_superfamily_windowed_habitat_boxplot.png\n\n")


# ==============================================================================
# SECTION 8: EXPORT
# ==============================================================================

cat("=== SECTION 8: EXPORT ===\n\n")

readr::write_csv(daily_df, "outputs/H3_superfamily_daily_results.csv")
readr::write_csv(windowed_df, "outputs/H3_superfamily_windowed_results.csv")
readr::write_csv(dplyr::bind_rows(daily_by_site, daily_grand),
                 "outputs/H3_superfamily_daily_summary.csv")
readr::write_csv(dplyr::bind_rows(win_by_site, win_grand),
                 "outputs/H3_superfamily_windowed_summary.csv")
cat("[OK] All outputs exported\n\n")


# ==============================================================================
# WRAP UP
# ==============================================================================

cat("==============================================================================\n")
cat("Script completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

while (sink.number() > 0) sink()