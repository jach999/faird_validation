# ==============================================================================
# H2 BRIDGE EXPLORATORY: Biomass Ratio vs Taxonomic Similarity (Unified)
# ==============================================================================
# Description: Evaluates if taxonomic mismatch drives biomass mismatch by
# correlating the AMMOD/FAIRD biomass ratio with the Izsak-Price taxonomic
# similarity across 3-day windows.
# Computed hierarchically at Site, Habitat, and Overall levels.
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(vegan)
library(patchwork) # For combining plots

script_name <- "H2_bridge_exploratory_unified"
output_file <- file.path("outputs/console", paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", script_name, "_output.txt"))
if (!dir.exists("outputs/console")) dir.create("outputs/console", recursive = TRUE)
if (!dir.exists("outputs/figures")) dir.create("outputs/figures", recursive = TRUE)

while (sink.number() > 0) sink()
sink(output_file, split = TRUE)

cat("==============================================================================\n")
cat("H2 BRIDGE: BIOMASS RATIO VS TAXONOMIC SIMILARITY (EXPLORATORY)\n")
cat("Hierarchical evaluation: Site, Habitat, and Overall Levels\n")
cat("==============================================================================\n\n")

source("R/00_load_data.R")

# --- 1. SETUP PARAMETERS ---
date_start   <- as.Date("2023-08-23")
date_end     <- as.Date("2023-09-03") # AMMOD biomass limit
date_exclude <- as.Date("2023-09-04")
nc           <- "#N/C"
MIN_INDIV    <- 3

window_breaks <- seq(date_start, as.Date("2023-09-13") + 2, by = 3)
window_lookup <- tibble::tibble(Date = seq.Date(date_start, as.Date("2023-09-13"), by = "day")) %>%
  dplyr::mutate(Window_ID = as.integer(cut(Date, breaks = window_breaks, right = FALSE, labels = FALSE)))

# --- 2. PREPARE RAW DATA (BIOMASS & TAXONOMY) ---
cat("Preparing raw datasets...\n")

# Biomass
bio_raw <- data_list$ammod_biomass_for_stats %>%
  dplyr::mutate(Device_type = "AMMOD", Biomass_mg = Live_mass) %>%
  dplyr::select(Date, Site, Device_type, Biomass_mg) %>%
  dplyr::bind_rows(
    data_list$etraps_daily %>%
      dplyr::filter(Device_type == "FAIRD") %>%
      dplyr::select(Date, Site, Device_type, biomass_mg) %>%
      dplyr::rename(Biomass_mg = biomass_mg)
  ) %>%
  dplyr::filter(Date >= date_start & Date <= date_end) %>%
  dplyr::left_join(window_lookup, by = "Date") %>%
  dplyr::filter(Date != date_exclude) %>%
  dplyr::mutate(Habitat = dplyr::if_else(Site %in% c("Site1", "Site2"), "Maize", "Meadow"))

# Taxonomy
faird_tax <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD", Date >= date_start & Date <= date_end) %>%
  dplyr::mutate(Taxon_SF = dplyr::if_else(!is.na(Superfamily) & Superfamily != nc & Superfamily != "", Superfamily, NA_character_),
                Habitat = dplyr::if_else(Site %in% c("Site1", "Site2"), "Maize", "Meadow"))

ammod_tax <- data_list$ammod_taxonomy %>%
  dplyr::filter(Date >= date_start & Date <= date_end, Date != date_exclude,
                !(Device == "AMMOD4" & Date >= as.Date("2023-08-28") & Date <= as.Date("2023-09-03"))) %>%
  dplyr::mutate(Taxon_SF = dplyr::if_else(!is.na(Superfamily) & Superfamily != nc & Superfamily != "", Superfamily, NA_character_),
                Habitat = dplyr::if_else(Site %in% c("Site1", "Site2"), "Maize", "Meadow"))

# Distance Matrix
all_sf_taxa <- dplyr::bind_rows(faird_tax, ammod_tax) %>%
  dplyr::filter(!is.na(Taxon_SF)) %>%
  dplyr::distinct(Taxon_SF, .keep_all = TRUE) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  as.data.frame()
rownames(all_sf_taxa) <- all_sf_taxa$Taxon_SF
all_sf_taxa <- all_sf_taxa %>% dplyr::select(-Taxon_SF)
dist_matrix_sf <- as.matrix(vegan::taxa2dist(all_sf_taxa, varstep = TRUE, check = TRUE))

calc_ip <- function(taxa_A, taxa_B, dist_mat) {
  taxa_A <- taxa_A[taxa_A %in% rownames(dist_mat)]
  taxa_B <- taxa_B[taxa_B %in% rownames(dist_mat)]
  if (length(taxa_A) == 0 | length(taxa_B) == 0) return(NA)
  min_A_to_B <- sapply(taxa_A, function(a) if (a %in% taxa_B) 0 else min(dist_mat[a, taxa_B], na.rm = TRUE))
  min_B_to_A <- sapply(taxa_B, function(b) if (b %in% taxa_A) 0 else min(dist_mat[b, taxa_A], na.rm = TRUE))
  return(1 - (mean(c(min_A_to_B, min_B_to_A)) / max(dist_mat, na.rm = TRUE)))
}

# ==============================================================================
# 3. HELPER FUNCTION TO PROCESS A GIVEN LEVEL
# ==============================================================================
process_level <- function(group_vars) {
  # Biomass aggregation
  bio_agg <- bio_raw %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "Window_ID", "Device_type")))) %>%
    dplyr::summarise(Sum_Biomass = sum(Biomass_mg, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Device_type, values_from = Sum_Biomass) %>%
    dplyr::filter(!is.na(AMMOD) & !is.na(FAIRD) & FAIRD > 0) %>%
    dplyr::mutate(Biomass_Ratio = AMMOD / FAIRD)
  
  # Taxonomic aggregation
  taxo_results <- list()
  group_combinations <- bio_agg %>% dplyr::select(dplyr::all_of(c(group_vars, "Window_ID"))) %>% dplyr::distinct()
  
  for (i in 1:nrow(group_combinations)) {
    curr_row <- group_combinations[i, ]
    w <- curr_row$Window_ID
    win_dates <- window_lookup %>% dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
    
    # Filter by group dynamically
    f_win <- faird_tax %>% dplyr::filter(Date %in% win_dates, !is.na(Taxon_SF))
    a_win <- ammod_tax %>% dplyr::filter(Date %in% win_dates, !is.na(Taxon_SF))
    
    for (v in group_vars) {
      f_win <- f_win[f_win[[v]] == curr_row[[v]], ]
      a_win <- a_win[a_win[[v]] == curr_row[[v]], ]
    }
    
    f_sf <- unique(f_win$Taxon_SF)
    a_sf <- unique(a_win$Taxon_SF)
    n_f <- nrow(f_win)
    
    ds <- NA_real_
    if (n_f >= MIN_INDIV & length(a_sf) > 0 & length(f_sf) > 0) {
      ds <- calc_ip(f_sf, a_sf, dist_matrix_sf)
    }
    
    res_row <- curr_row
    res_row$Delta_S_pct <- ds * 100
    taxo_results[[length(taxo_results) + 1]] <- res_row
  }
  
  taxo_agg <- dplyr::bind_rows(taxo_results)
  
  # Merge
  final_df <- bio_agg %>%
    dplyr::inner_join(taxo_agg, by = c(group_vars, "Window_ID")) %>%
    dplyr::filter(!is.na(Delta_S_pct))
  
  return(final_df)
}

# ==============================================================================
# 4. RUN ANALYSES FOR ALL LEVELS
# ==============================================================================
cat("Running computations for Site, Habitat, and Overall levels...\n\n")

df_site    <- process_level(c("Site"))
df_habitat <- process_level(c("Habitat"))
df_overall <- process_level(character(0)) # Empty vector for overall (no grouping)

# Correlation tests
cor_site    <- cor.test(df_site$Biomass_Ratio,    df_site$Delta_S_pct,    method = "spearman")
cor_habitat <- cor.test(df_habitat$Biomass_Ratio, df_habitat$Delta_S_pct, method = "spearman")
cor_overall <- cor.test(df_overall$Biomass_Ratio, df_overall$Delta_S_pct, method = "spearman")

# Print results
cat("--- SITE LEVEL ---\n")
cat(sprintf("n = %d | rho = %.3f | p = %.4f\n\n", nrow(df_site), cor_site$estimate, cor_site$p.value))

cat("--- HABITAT LEVEL ---\n")
cat(sprintf("n = %d | rho = %.3f | p = %.4f\n\n", nrow(df_habitat), cor_habitat$estimate, cor_habitat$p.value))

cat("--- OVERALL LEVEL ---\n")
cat(sprintf("n = %d | rho = %.3f | p = %.4f\n\n", nrow(df_overall), cor_overall$estimate, cor_overall$p.value))

# ==============================================================================
# 5. GENERATE UNIFIED PLOT
# ==============================================================================
cat("Generating unified plot...\n")

site_cols <- if(exists("custom_colors")) custom_colors else c("Site1"="#E69F00", "Site2"="#56B4E9", "Site3"="#009E73", "Site4"="#F0E442")

p1 <- ggplot(df_site, aes(x = Biomass_Ratio, y = Delta_S_pct)) +
  geom_point(aes(color = Site), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  scale_x_log10() + scale_color_manual(values = site_cols) +
  labs(title = "Site Level", subtitle = sprintf("\u03c1 = %.3f (p=%.3f)", cor_site$estimate, cor_site$p.value),
       x = "Biomass Ratio (log)", y = "Taxonomic Similarity (%)") + theme_minimal()

p2 <- ggplot(df_habitat, aes(x = Biomass_Ratio, y = Delta_S_pct)) +
  geom_point(aes(color = Habitat), size = 4, alpha = 0.9) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  scale_x_log10() + scale_color_manual(values = c("Maize"="#E69F00", "Meadow"="#009E73")) +
  labs(title = "Habitat Level", subtitle = sprintf("\u03c1 = %.3f (p=%.3f)", cor_habitat$estimate, cor_habitat$p.value),
       x = "Biomass Ratio (log)", y = "") + theme_minimal()

p3 <- ggplot(df_overall, aes(x = Biomass_Ratio, y = Delta_S_pct)) +
  geom_point(color = "steelblue", size = 5, alpha = 0.9) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  scale_x_log10() +
  labs(title = "Overall Level", subtitle = sprintf("\u03c1 = %.3f (p=%.3f)", cor_overall$estimate, cor_overall$p.value),
       x = "Biomass Ratio (log)", y = "") + theme_minimal()

# Combine with patchwork
unified_plot <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(
    title = "Exploratory Analysis: Does Taxonomic Mismatch drive Biomass Mismatch?",
    subtitle = "Spearman correlations across different spatial aggregation levels (3-day windows)"
  )

print(unified_plot)
ggsave("outputs/figures/H2_bridge_unified_exploratory.png", plot = unified_plot, width = 14, height = 5, dpi = 300, bg = "white")
cat("[OK] Saved plot to: outputs/figures/H2_bridge_unified_exploratory.png\n")

cat("\n==============================================================================\n")
cat("Script completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")
while (sink.number() > 0) sink()