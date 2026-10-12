# ==============================================================================
# H1: INTRA-SYSTEM TAXONOMIC CONSISTENCY — AMMOD REPLICATES
# ==============================================================================
# Purpose: Assess whether replicate AMMOD traps detect the same taxa
#          on the same days, using Izsak-Price similarity at Superfamily level.
#
# Metric:  Izsak-Price taxonomic similarity (presence/absence)
#          AMMOD data is qualitative (metabarcoding OTUs) — no abundance info.
#
# Level:   SUPERFAMILY (consistent with FAIRD H1 taxonomic and H3 analyses)
#
# Temporal resolution:
#   - DAILY:     original resolution (no filtering)
#   - 3-DAY:     non-overlapping windows, filtered >= 3 OTUs with Superfamily
#                per device per window
#
# Pairs:   AMMOD1 <-> AMMOD2 (Maize, Sites 1-2) — 21 valid taxo days each
#          AMMOD3 <-> AMMOD4 (Meadow, Sites 3-4) — 21 and 13 valid taxo days
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
script_name <- "H1_consistency_AMMOD_taxonomic"

output_file <- file.path("outputs",
                         paste0(format(Sys.time(), "%Y%m%d_%H%M%S"),
                                "_", script_name, ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)
while (sink.number() > 0) sink()
sink(output_file, split = TRUE)

cat("==============================================================================\n")
cat("H1: INTRA-SYSTEM TAXONOMIC CONSISTENCY — AMMOD REPLICATES\n")
cat("Izsak-Price similarity at Superfamily level\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: LOAD DATA
# ==============================================================================

cat("=== SECTION 1: LOAD DATA ===\n\n")

source("load_data.R")
source("custom_colors.R")

nc <- "#N/C"

# --- Temporal scope ---
date_start <- as.Date("2023-08-23")
date_end   <- as.Date("2023-09-13")
date_exclude <- as.Date("2023-09-04")  # Bottle change day
ammod4_failure_start <- as.Date("2023-08-28")
ammod4_failure_end   <- as.Date("2023-09-03")

# --- AMMOD data: extract and assign Superfamily ---
ammod_data <- data_list$ammod_taxonomy %>%
  dplyr::filter(
    Date >= date_start & Date <= date_end,
    Date != date_exclude,
    !(Device == "AMMOD4" & Date >= ammod4_failure_start & Date <= ammod4_failure_end)
  ) %>%
  dplyr::select(Date, Device, Site, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus, Species) %>%
  dplyr::mutate(
    Taxon_SF = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc & Superfamily != "",
      Superfamily, NA_character_
    ),
    Taxon_Genus = dplyr::if_else(
      !is.na(Genus) & Genus != nc & Genus != "",
      Genus, NA_character_
    ),
    Taxon_Species = dplyr::if_else(
      !is.na(Species) & Species != nc & Species != "",
      Species, NA_character_
    )
  )

cat("AMMOD loaded:", nrow(ammod_data), "OTUs\n")
cat("  With Superfamily:", sum(!is.na(ammod_data$Taxon_SF)),
    sprintf("(%.1f%%)\n", 100 * sum(!is.na(ammod_data$Taxon_SF)) / nrow(ammod_data)))
cat("  With Genus:", sum(!is.na(ammod_data$Taxon_Genus)),
    sprintf("(%.1f%%)\n", 100 * sum(!is.na(ammod_data$Taxon_Genus)) / nrow(ammod_data)))
cat("  With Species:", sum(!is.na(ammod_data$Taxon_Species)),
    sprintf("(%.1f%%)\n", 100 * sum(!is.na(ammod_data$Taxon_Species)) / nrow(ammod_data)))

# Per-device summary
cat("\nPer-device coverage:\n")
dev_summary <- ammod_data %>%
  dplyr::group_by(Device) %>%
  dplyr::summarise(
    n_otus = dplyr::n(),
    n_sf   = sum(!is.na(Taxon_SF)),
    pct_sf = round(100 * sum(!is.na(Taxon_SF)) / dplyr::n(), 1),
    n_days = length(unique(Date)),
    .groups = "drop"
  )
print(dev_summary)

# AMMOD4 exclusion details
cat("\nNote: AMMOD4 has reduced valid days (13 taxo) due to mechanical failure.\n")
cat("      AMMOD4 Sep 10 excluded (contaminated sample R2_4_06).\n\n")


# ==============================================================================
# SECTION 2: TAXONOMIC DISTANCE MATRICES (SPECIES, GENUS, SUPERFAMILY)
# ==============================================================================

cat("=== SECTION 2: DISTANCE MATRICES ===\n\n")

cat("--- Superfamily-level distance matrix ---\n")

# Collect unique Superfamilies with hierarchy
ammod_sf_taxa <- ammod_data %>%
  dplyr::filter(!is.na(Taxon_SF)) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  dplyr::distinct(Taxon_SF, .keep_all = TRUE)

cat("Unique AMMOD Superfamilies:", nrow(ammod_sf_taxa), "\n")

# Build hierarchy table for taxa2dist
tax_table_sf <- ammod_sf_taxa %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  as.data.frame()

rownames(tax_table_sf) <- tax_table_sf$Taxon_SF
tax_table_sf <- tax_table_sf %>% dplyr::select(-Taxon_SF)

taxon_distances_sf <- vegan::taxa2dist(tax_table_sf, varstep = TRUE, check = TRUE)
dist_matrix_sf <- as.matrix(taxon_distances_sf)

cat("Distance matrix:", nrow(dist_matrix_sf), "x", ncol(dist_matrix_sf), "\n")
cat("Distance range:", round(min(taxon_distances_sf), 4), "-",
    round(max(taxon_distances_sf), 4), "\n\n")

# --- Genus-level distance matrix ---
cat("--- Genus-level distance matrix ---\n")

ammod_genus_taxa <- ammod_data %>%
  dplyr::filter(!is.na(Taxon_Genus)) %>%
  dplyr::select(Taxon_Genus, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus) %>%
  dplyr::distinct(Taxon_Genus, .keep_all = TRUE)

cat("Unique AMMOD Genera:", nrow(ammod_genus_taxa), "\n")

tax_table_gen <- ammod_genus_taxa %>%
  dplyr::select(Taxon_Genus, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus) %>%
  as.data.frame()
rownames(tax_table_gen) <- tax_table_gen$Taxon_Genus
tax_table_gen <- tax_table_gen %>% dplyr::select(-Taxon_Genus)

taxon_distances_gen <- vegan::taxa2dist(tax_table_gen, varstep = TRUE, check = TRUE)
dist_matrix_gen <- as.matrix(taxon_distances_gen)
cat("Distance matrix:", nrow(dist_matrix_gen), "x", ncol(dist_matrix_gen), "\n\n")

# --- Species-level distance matrix ---
cat("--- Species-level distance matrix ---\n")

ammod_species_taxa <- ammod_data %>%
  dplyr::filter(!is.na(Taxon_Species)) %>%
  dplyr::select(Taxon_Species, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus, Species) %>%
  dplyr::distinct(Taxon_Species, .keep_all = TRUE)

cat("Unique AMMOD Species:", nrow(ammod_species_taxa), "\n")

tax_table_sp <- ammod_species_taxa %>%
  dplyr::select(Taxon_Species, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus, Species) %>%
  as.data.frame()
rownames(tax_table_sp) <- tax_table_sp$Taxon_Species
tax_table_sp <- tax_table_sp %>% dplyr::select(-Taxon_Species)

taxon_distances_sp <- vegan::taxa2dist(tax_table_sp, varstep = TRUE, check = TRUE)
dist_matrix_sp <- as.matrix(taxon_distances_sp)
cat("Distance matrix:", nrow(dist_matrix_sp), "x", ncol(dist_matrix_sp), "\n\n")


# ==============================================================================
# SECTION 3: IZSAK-PRICE FUNCTION
# ==============================================================================

calculate_taxonomic_similarity <- function(taxa_A, taxa_B, dist_matrix) {
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
# SECTION 4: DEFINE PAIRS AND WINDOWS
# ==============================================================================

cat("=== SECTION 4: DEFINE PAIRS AND WINDOWS ===\n\n")

pairs <- list(
  list(dev_a = "AMMOD1", dev_b = "AMMOD2", habitat = "Maize"),
  list(dev_a = "AMMOD3", dev_b = "AMMOD4", habitat = "Meadow")
)

all_dates <- seq.Date(date_start, date_end, by = "day")

# 3-day windows (same as FAIRD H1)
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
  ) %>%
  dplyr::filter(!is.na(Window_ID))

cat("Windows defined:\n")
print(window_lookup %>%
        dplyr::group_by(Window_ID, Window_Label) %>%
        dplyr::summarise(n_days = dplyr::n(), .groups = "drop"),
      n = Inf)
cat("\n")


# ==============================================================================
# SECTION 5: COMPUTE DAILY IZSAK-PRICE
# ==============================================================================

cat("=== SECTION 5: DAILY IZSAK-PRICE (ALL LEVELS) ===\n\n")

MIN_OTUS <- 3
cat(sprintf("Minimum OTUs per device (for windowed analysis): %d\n\n", MIN_OTUS))

# Define taxonomic levels for daily analysis
tax_levels <- list(
  list(col = "Taxon_Species",  label = "Species",     dmat = dist_matrix_sp),
  list(col = "Taxon_Genus",    label = "Genus",        dmat = dist_matrix_gen),
  list(col = "Taxon_SF",       label = "Superfamily",  dmat = dist_matrix_sf)
)

daily_results <- list()

for (p in seq_along(pairs)) {
  pair <- pairs[[p]]
  cat(sprintf("Processing: %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))
  
  for (lev in tax_levels) {
    for (d in all_dates) {
      date_val <- as.Date(d, origin = "1970-01-01")
      if (date_val == date_exclude) next
      
      taxa_a <- ammod_data %>%
        dplyr::filter(Date == date_val, Device == pair$dev_a,
                      !is.na(.data[[lev$col]])) %>%
        dplyr::pull(.data[[lev$col]]) %>% unique()
      
      taxa_b <- ammod_data %>%
        dplyr::filter(Date == date_val, Device == pair$dev_b,
                      !is.na(.data[[lev$col]])) %>%
        dplyr::pull(.data[[lev$col]]) %>% unique()
      
      delta_s <- NA_real_
      if (length(taxa_a) > 0 & length(taxa_b) > 0) {
        delta_s <- calculate_taxonomic_similarity(taxa_a, taxa_b, lev$dmat)
      }
      
      daily_results[[length(daily_results) + 1]] <- tibble::tibble(
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

daily_df <- dplyr::bind_rows(daily_results)
cat("\nDaily computations:", nrow(daily_df), "\n")
cat("Valid Delta_S:", sum(!is.na(daily_df$Delta_S)), "\n\n")

# ==============================================================================
# SECTION 5.5: INTER-HABITAT DAILY IZSAK-PRICE (ALL LEVELS)
# ==============================================================================

cat("=== SECTION 5.5: INTER-HABITAT DAILY IZSAK-PRICE (ALL LEVELS) ===\n\n")

inter_daily_results <- list()

cat("Processing Inter-Habitat: Maize_Pool (AMMOD1+2) <-> Meadow_Pool (AMMOD3+4)\n")

for (lev in tax_levels) {
  for (d in all_dates) {
    date_val <- as.Date(d, origin = "1970-01-01")
    # Skip excluded dates (e.g., bottle change day)
    if (date_val == date_exclude) next
    
    # Pool presences for Maize (AMMOD 1 and 2)
    taxa_maize <- ammod_data %>%
      dplyr::filter(Date == date_val, Device %in% c("AMMOD1", "AMMOD2"),
                    !is.na(.data[[lev$col]])) %>%
      dplyr::pull(.data[[lev$col]]) %>% unique()
    
    # Pool presences for Meadow (AMMOD 3 and 4)
    taxa_meadow <- ammod_data %>%
      dplyr::filter(Date == date_val, Device %in% c("AMMOD3", "AMMOD4"),
                    !is.na(.data[[lev$col]])) %>%
      dplyr::pull(.data[[lev$col]]) %>% unique()
    
    delta_s <- NA_real_
    if (length(taxa_maize) > 0 & length(taxa_meadow) > 0) {
      delta_s <- calculate_taxonomic_similarity(taxa_maize, taxa_meadow, lev$dmat)
    }
    
    inter_daily_results[[length(inter_daily_results) + 1]] <- tibble::tibble(
      Resolution = "Daily",
      Pair       = "Maize vs Meadow",
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

inter_daily_df <- dplyr::bind_rows(inter_daily_results)

cat("\nInter-Habitat daily computations:", nrow(inter_daily_df), "\n")
cat("Valid Delta_S:", sum(!is.na(inter_daily_df$Delta_S)), "\n\n")

# Print the summary
cat("--- Inter-Habitat Daily Summary ---\n\n")
inter_daily_summary <- inter_daily_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Pair, Tax_Level) %>%
  dplyr::summarise(
    n_days   = dplyr::n(),
    mean_DS  = round(mean(Delta_S) * 100, 1),
    sd_DS    = round(sd(Delta_S) * 100, 1),
    median_DS = round(median(Delta_S) * 100, 1),
    min_DS   = round(min(Delta_S) * 100, 1),
    max_DS   = round(max(Delta_S) * 100, 1),
    mean_taxa_Maize  = round(mean(N_taxa_A), 1),
    mean_taxa_Meadow = round(mean(N_taxa_B), 1),
    mean_shared      = round(mean(N_shared), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Species", "Genus", "Superfamily")))

print(inter_daily_summary, n = Inf, width = 200)

# Export specific results
readr::write_csv(inter_daily_df, "outputs/H1_AMMOD_taxonomic_inter_habitat_daily_results.csv")
cat("\n[OK] outputs/H1_AMMOD_taxonomic_inter_habitat_daily_results.csv\n\n")

# ==============================================================================
# SECTION 6: COMPUTE 3-DAY WINDOWED IZSAK-PRICE (FILTERED)
# ==============================================================================

cat("=== SECTION 6: 3-DAY WINDOWED IZSAK-PRICE ===\n\n")

windowed_results <- list()

for (p in seq_along(pairs)) {
  pair <- pairs[[p]]
  cat(sprintf("Processing: %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))
  
  for (w in unique(window_lookup$Window_ID)) {
    win_dates <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
    win_label <- window_lookup %>%
      dplyr::filter(Window_ID == w) %>% dplyr::pull(Window_Label) %>%
      unique()
    
    # Remove excluded dates from window
    win_dates <- win_dates[win_dates != date_exclude]
    
    sf_a <- ammod_data %>%
      dplyr::filter(Date %in% win_dates, Device == pair$dev_a, !is.na(Taxon_SF)) %>%
      dplyr::pull(Taxon_SF) %>% unique()
    
    sf_b <- ammod_data %>%
      dplyr::filter(Date %in% win_dates, Device == pair$dev_b, !is.na(Taxon_SF)) %>%
      dplyr::pull(Taxon_SF) %>% unique()
    
    passes_filter <- (length(sf_a) >= MIN_OTUS & length(sf_b) >= MIN_OTUS)
    
    delta_s <- NA_real_
    if (passes_filter) {
      delta_s <- calculate_taxonomic_similarity(sf_a, sf_b, dist_matrix_sf)
    }
    
    windowed_results[[length(windowed_results) + 1]] <- tibble::tibble(
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
      Filtered   = !passes_filter
    )
  }
}

windowed_df <- dplyr::bind_rows(windowed_results)
cat("\nWindowed computations:", nrow(windowed_df), "\n")
cat("Valid:", sum(!is.na(windowed_df$Delta_S)),
    "| Filtered:", sum(windowed_df$Filtered), "\n\n")

# Filter impact
cat("--- Filter impact ---\n")
filter_impact <- windowed_df %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    n_windows    = dplyr::n(),
    n_valid      = sum(!is.na(Delta_S)),
    n_filtered   = sum(Filtered),
    pct_retained = round(100 * sum(!is.na(Delta_S)) / dplyr::n(), 1),
    .groups = "drop"
  )
print(filter_impact)
cat("\n")

# ==============================================================================
# SECTION 6.5: INTER-HABITAT IZSAK-PRICE (MAIZE VS MEADOW)
# ==============================================================================

cat("=== SECTION 6.5: INTER-HABITAT IZSAK-PRICE (MAIZE VS MEADOW) ===\n\n")

inter_habitat_results <- list()

cat("Processing Inter-Habitat: Maize_Pool (AMMOD1+2) <-> Meadow_Pool (AMMOD3+4) (Superfamily Level)\n")

for (w in unique(window_lookup$Window_ID)) {
  win_dates <- window_lookup %>%
    dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
  win_label <- window_lookup %>%
    dplyr::filter(Window_ID == w) %>% dplyr::pull(Window_Label) %>%
    unique()
  
  # Remove excluded dates (e.g., bottle change day)
  win_dates <- win_dates[win_dates != date_exclude]
  
  # Pool presences for Maize (AMMOD 1 and 2)
  sf_maize <- ammod_data %>%
    dplyr::filter(Date %in% win_dates, Device %in% c("AMMOD1", "AMMOD2"), !is.na(Taxon_SF)) %>%
    dplyr::pull(Taxon_SF) %>% unique()
  
  # Pool presences for Meadow (AMMOD 3 and 4)
  sf_meadow <- ammod_data %>%
    dplyr::filter(Date %in% win_dates, Device %in% c("AMMOD3", "AMMOD4"), !is.na(Taxon_SF)) %>%
    dplyr::pull(Taxon_SF) %>% unique()
  
  # Apply filter: minimum 3 OTUs per HABITAT per window
  passes_filter <- (length(sf_maize) >= MIN_OTUS & length(sf_meadow) >= MIN_OTUS)
  
  delta_s <- NA_real_
  if (passes_filter) {
    delta_s <- calculate_taxonomic_similarity(sf_maize, sf_meadow, dist_matrix_sf)
  }
  
  inter_habitat_results[[length(inter_habitat_results) + 1]] <- tibble::tibble(
    Resolution = "3-day window",
    Pair       = "Maize vs Meadow",
    Habitat    = "Inter-Habitat",
    Period     = win_label,
    Window_ID  = w,
    Tax_Level  = "Superfamily",
    Delta_S    = delta_s,
    N_taxa_A   = length(sf_maize),
    N_taxa_B   = length(sf_meadow),
    N_shared   = length(intersect(sf_maize, sf_meadow)),
    Filtered   = !passes_filter
  )
}

inter_habitat_df <- dplyr::bind_rows(inter_habitat_results)

cat("\nInter-Habitat windowed computations:", nrow(inter_habitat_df), "\n")
cat("Valid Delta_S (passed filter):", sum(!is.na(inter_habitat_df$Delta_S)), "\n")
cat("Filtered out (< ", MIN_OTUS, " OTUs):", sum(inter_habitat_df$Filtered), "\n\n")

# Print the summary
cat("--- Inter-Habitat Summary (Filtered >= 3 OTUs) ---\n\n")
summary_inter <- inter_habitat_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::summarise(
    Pair = "Maize vs Meadow",
    n_windows = dplyr::n(),
    mean_DS   = round(mean(Delta_S) * 100, 1),
    sd_DS     = round(sd(Delta_S) * 100, 1),
    median_DS = round(median(Delta_S) * 100, 1),
    min_DS    = round(min(Delta_S) * 100, 1),
    max_DS    = round(max(Delta_S) * 100, 1),
    mean_taxa_Maize  = round(mean(N_taxa_A), 1),
    mean_taxa_Meadow = round(mean(N_taxa_B), 1),
    mean_shared      = round(mean(N_shared), 1),
    .groups = "drop"
  )
print(summary_inter, width = 200)

# Export specific results
readr::write_csv(inter_habitat_df, "outputs/H1_AMMOD_taxonomic_inter_habitat_results.csv")
cat("\n[OK] outputs/H1_AMMOD_taxonomic_inter_habitat_results.csv\n\n")

# ==============================================================================
# SECTION 7: SUMMARY STATISTICS
# ==============================================================================

cat("=== SECTION 7: SUMMARY STATISTICS ===\n\n")

# --- 7a: Daily ---
cat("--- 7a: Daily summary ---\n\n")

daily_summary <- daily_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_days   = dplyr::n(),
    mean_DS  = round(mean(Delta_S) * 100, 1),
    sd_DS    = round(sd(Delta_S) * 100, 1),
    median_DS = round(median(Delta_S) * 100, 1),
    min_DS   = round(min(Delta_S) * 100, 1),
    max_DS   = round(max(Delta_S) * 100, 1),
    mean_taxa_A = round(mean(N_taxa_A), 1),
    mean_taxa_B = round(mean(N_taxa_B), 1),
    mean_shared = round(mean(N_shared), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Pair, factor(Tax_Level, levels = c("Species", "Genus", "Superfamily")))
print(daily_summary, n = Inf, width = 200)

# --- 7b: Windowed ---
cat("\n--- 7b: Windowed summary (filtered >= 3 OTUs) ---\n\n")

win_summary <- windowed_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    n_windows = dplyr::n(),
    mean_DS   = round(mean(Delta_S) * 100, 1),
    sd_DS     = round(sd(Delta_S) * 100, 1),
    median_DS = round(median(Delta_S) * 100, 1),
    min_DS    = round(min(Delta_S) * 100, 1),
    max_DS    = round(max(Delta_S) * 100, 1),
    mean_taxa_A = round(mean(N_taxa_A), 1),
    mean_taxa_B = round(mean(N_taxa_B), 1),
    mean_shared = round(mean(N_shared), 1),
    .groups = "drop"
  )
print(win_summary, width = 200)

readr::write_csv(win_summary, "outputs/H1_AMMOD_taxonomic_windowed_summary.csv")
cat("\n[OK] outputs/H1_AMMOD_taxonomic_windowed_summary.csv\n\n")

# --- 7c: Comparison with FAIRD H1 (context) ---
cat("--- 7c: AMMOD vs FAIRD reference context ---\n\n")
cat("AMMOD consistency serves as the REFERENCE LEVEL for H1 taxonomic.\n")
cat("FAIRD H1 taxonomic results (Superfamily, windowed, M-H):\n")
cat("  FAIRD1 vs FAIRD2 (Maize):  88.6% +/- 17.3% (n=5)\n")
cat("  FAIRD3 vs FAIRD4 (Meadow): 82.9% +/- 13.6% (n=7)\n")
cat("NOTE: Different metrics (M-H vs Izsak-Price) — values not directly comparable.\n")
cat("      Both are [0,1] similarity indices, but M-H weighs by abundance\n")
cat("      while Izsak-Price weighs by taxonomic distance.\n\n")


# ==============================================================================
# SECTION 8: FIGURES
# ==============================================================================

cat("=== SECTION 8: FIGURES ===\n\n")

# --- 8a: Intra-Habitat Level — Daily Species temporal ---
cat("Creating Intra-Habitat Level temporal plot (Daily, Species)...\n")

plot_intra_daily <- daily_df %>%
  dplyr::filter(Tax_Level == "Species", !is.na(Delta_S)) %>%
  dplyr::mutate(DS_pct = Delta_S * 100)

# Compute summary stats per pair for subtitle
intra_stats <- plot_intra_daily %>%
  dplyr::group_by(Habitat) %>%
  dplyr::summarise(
    mean_ds = round(mean(DS_pct), 1),
    sd_ds   = round(sd(DS_pct), 1),
    n       = dplyr::n(),
    .groups = "drop"
  )

p_intra <- ggplot(plot_intra_daily,
                  aes(x = Date, y = DS_pct, color = Habitat)) +
  geom_line(aes(group = Pair), linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = custom_colors,
                     labels = c("Maize"  = "AMMOD1 vs AMMOD2 (Maize)",
                                "Meadow" = "AMMOD3 vs AMMOD4 (Meadow)")) +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  labs(
    title    = "H1 Qualitative - Taxonomic Consistency: AMMOD Replicates — Intra-Habitat Level",
    subtitle = sprintf("Izsak-Price index — Daily — Species level | Maize: %.1f%% +/- %.1f%% (n=%d) | Meadow: %.1f%% +/- %.1f%% (n=%d)",
                       intra_stats$mean_ds[intra_stats$Habitat == "Maize"],
                       intra_stats$sd_ds[intra_stats$Habitat == "Maize"],
                       intra_stats$n[intra_stats$Habitat == "Maize"],
                       intra_stats$mean_ds[intra_stats$Habitat == "Meadow"],
                       intra_stats$sd_ds[intra_stats$Habitat == "Meadow"],
                       intra_stats$n[intra_stats$Habitat == "Meadow"]),
    x = "Date", y = "Izsak-Price Similarity (%)",
    color = "Pair",
    caption = "Dashed line = 50% similarity reference."
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

print(p_intra)
ggsave("outputs/H1_AMMOD_taxonomic_daily_species_intra_habitat.png",
       plot = p_intra, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/H1_AMMOD_taxonomic_daily_species_intra_habitat.png\n\n")


# --- 8b: Inter-Habitat Level — Daily Species temporal ---
cat("Creating Inter-Habitat Level temporal plot (Daily, Species)...\n")

plot_inter_daily <- inter_daily_df %>%
  dplyr::filter(Tax_Level == "Species", !is.na(Delta_S)) %>%
  dplyr::mutate(DS_pct = Delta_S * 100)

inter_mean <- mean(plot_inter_daily$DS_pct)

p_inter <- ggplot(plot_inter_daily,
                  aes(x = Date, y = DS_pct)) +
  geom_line(linewidth = 0.8, color = custom_colors["AMMOD"]) +
  geom_point(size = 2, color = custom_colors["AMMOD"]) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = inter_mean,
             linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  labs(
    title    = "H1 Qualitative - Taxonomic Consistency: AMMOD Replicates — Inter-Habitat Level",
    subtitle = sprintf("Izsak-Price index — Daily — Species level | Maize vs Meadow: %.1f%% +/- %.1f%% (n=%d)",
                       inter_mean, sd(plot_inter_daily$DS_pct),
                       nrow(plot_inter_daily)),
    x = "Date", y = "Izsak-Price Similarity (%)",
    caption = "Dashed line = 50% reference. Dotted line = mean."
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

print(p_inter)
ggsave("outputs/H1_AMMOD_taxonomic_daily_species_inter_habitat.png",
       plot = p_inter, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/H1_AMMOD_taxonomic_daily_species_inter_habitat.png\n\n")


# --- 8c: Inter-Habitat Level — Daily all taxonomic levels ---
cat("Creating Inter-Habitat Level multi-level temporal plot (Daily)...\n")

plot_inter_all <- inter_daily_df %>%
  dplyr::filter(!is.na(Delta_S)) %>%
  dplyr::mutate(
    DS_pct = Delta_S * 100,
    Tax_Level = factor(Tax_Level, levels = c("Species", "Genus", "Superfamily"))
  )

# Summary stats per level for subtitle
inter_all_stats <- plot_inter_all %>%
  dplyr::group_by(Tax_Level) %>%
  dplyr::summarise(
    mean_ds = round(mean(DS_pct), 1),
    sd_ds   = round(sd(DS_pct), 1),
    n       = dplyr::n(),
    .groups = "drop"
  )

level_colors <- c(
  "Species"     = unname(custom_colors["AMMOD"]),
  "Genus"       =  "darkorchid1",
  "Superfamily" = "darkred"
)

p_inter_all <- ggplot(plot_inter_all,
                      aes(x = Date, y = DS_pct, color = Tax_Level)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  scale_color_manual(values = level_colors) +
  labs(
    title    = "H1 Qualitative - Taxonomic Consistency: AMMOD Replicates — Inter-Habitat Level",
    subtitle = sprintf("Izsak-Price index — Daily — Maize vs Meadow | Species: %.1f%% | Genus: %.1f%% | Superfamily: %.1f%%",
                       inter_all_stats$mean_ds[inter_all_stats$Tax_Level == "Species"],
                       inter_all_stats$mean_ds[inter_all_stats$Tax_Level == "Genus"],
                       inter_all_stats$mean_ds[inter_all_stats$Tax_Level == "Superfamily"]),
    x = "Date", y = "Izsak-Price Similarity (%)",
    color = "Taxonomic Level",
    caption = "Dashed line = 50% similarity reference."
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

print(p_inter_all)
ggsave("outputs/H1_AMMOD_taxonomic_daily_all_levels_inter_habitat.png",
       plot = p_inter_all, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/H1_AMMOD_taxonomic_daily_all_levels_inter_habitat.png\n\n")


# --- 8d: Daily detail per window (table print) ---
cat("--- Detailed windowed results ---\n\n")
windowed_detail <- windowed_df %>%
  dplyr::mutate(DS_pct = round(Delta_S * 100, 1)) %>%
  dplyr::select(Pair, Period, N_taxa_A, N_taxa_B, N_shared, DS_pct, Filtered)
print(windowed_detail, n = Inf)


# ==============================================================================
# SECTION 9: EXPORT
# ==============================================================================

cat("\n=== SECTION 9: EXPORT ===\n\n")

readr::write_csv(daily_df, "outputs/H1_AMMOD_taxonomic_daily_results.csv")
readr::write_csv(windowed_df, "outputs/H1_AMMOD_taxonomic_windowed_results.csv")
cat("[OK] All outputs exported\n\n")


# ==============================================================================
# WRAP UP
# ==============================================================================

cat("==============================================================================\n")
cat("Script completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

while (sink.number() > 0) sink()