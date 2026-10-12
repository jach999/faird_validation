# =============================================================================
# BrayCurtis_FAIRD_Similarity.R
# Abundance-based taxonomic similarity (Bray-Curtis) between FAIRD units
# and between habitats at genus level.
#
# PURPOSE:
#   Complement the qualitative Izsak-Price presence/absence similarity (H3)
#   with a quantitative, abundance-weighted similarity metric. Bray-Curtis
#   dissimilarity accounts for relative genus-level abundance differences
#   between FAIRD replicates — something Izsak-Price cannot capture.
#
# ANALYSES:
#   1. Daily pairwise Bray-Curtis similarity between FAIRD replicates
#      - Intra-habitat: FAIRD1 vs FAIRD2 (Maize), FAIRD3 vs FAIRD4 (Meadow)
#      - Cross-habitat: Pool(FAIRD1+FAIRD2) vs Pool(FAIRD3+FAIRD4)
#   2. Overall (full-period) pooled Bray-Curtis between habitats
#   3. Summary statistics and visualisation
#
# INPUT:  data_cleaned.RData (data_list$etraps)
#         analysis_functions.R (get_best_classification)
#         custom_colors.R
#
# OUTPUT: outputs/BrayCurtis_daily_pairs.csv
#         outputs/BrayCurtis_pooled_habitats.csv
#         outputs/BrayCurtis_violin_daily.png/.pdf
#         outputs/BrayCurtis_heatmap_pooled.png/.pdf
#
# AUTHOR: [Author initials] — generated via Claude (March 2026)
# =============================================================================

library(tidyverse)
library(vegan)       # vegdist() for Bray-Curtis
library(ggplot2)
library(patchwork)

# --- 0. Load project files ---------------------------------------------------

load("data_cleaned.RData")            # produces data_list
source("analysis_functions.R")        # get_best_classification()
source("custom_colors.R")             # project colour palette

etraps_raw <- data_list$etraps        # individual-level e-trap records

# --- 1. Prepare FAIRD genus-level abundance data -----------------------------

# Apply get_best_classification() to obtain the most resolved taxonomic level
# available per individual. For Bray-Curtis we keep genus-level resolutions
# only (consistent with Izsak-Price analysis in H3).

faird_raw <- etraps_raw %>%
  filter(Device_type == "FAIRD") %>%
  rowwise() %>%
  mutate(best_taxon = get_best_classification(cur_data())) %>%
  ungroup() %>%
  filter(!is.na(best_taxon))

# Retain only records classified to Genus level
# (get_best_classification returns the finest resolved level label;
#  adapt the filter below if the function returns a level indicator separately)
faird_genus <- faird_raw %>%
  filter(Taxon_level == "Genus")      # ← adjust column name if needed

# NOTE: If get_best_classification() does not return a Taxon_level column,
# replace the filter above with:
#   filter(str_detect(best_taxon, "^[A-Z][a-z]+$"))  # simple genus regex
# or filter by the genus column directly:
#   filter(!is.na(Genus))             # ← adjust to actual column name

# Daily abundance matrix: Date × Genus, per Device
faird_daily <- faird_genus %>%
  group_by(Device, Date, best_taxon) %>%
  summarise(n = n(), .groups = "drop") %>%
  rename(Genus = best_taxon)

# Add Ambient (habitat) label for each device
device_habitat <- tibble(
  Device  = c("FAIRD1", "FAIRD2", "FAIRD3", "FAIRD4"),
  Ambient = c("Maize",  "Maize",  "Meadow", "Meadow")
)

faird_daily <- faird_daily %>%
  left_join(device_habitat, by = "Device")

# =============================================================================
# ANALYSIS 1 — Daily pairwise Bray-Curtis similarity
# =============================================================================
# For each day, compute Bray-Curtis similarity (= 1 − dissimilarity) between:
#   Pair A: FAIRD1 vs FAIRD2  (intra-Maize)
#   Pair B: FAIRD3 vs FAIRD4  (intra-Meadow)
#   Pair C: Pool(FAIRD1+2) vs Pool(FAIRD3+4)  (cross-habitat)

compute_daily_bc <- function(data, device_a, device_b, pair_label) {
  
  # Find days when BOTH devices have at least one detection
  days_a <- data %>% filter(Device == device_a) %>% pull(Date) %>% unique()
  days_b <- data %>% filter(Device == device_b) %>% pull(Date) %>% unique()
  shared_days <- intersect(days_a, days_b)
  
  if (length(shared_days) == 0) {
    message("No shared detection days for pair: ", pair_label)
    return(tibble())
  }
  
  results <- map_dfr(shared_days, function(d) {
    
    vec_a <- data %>%
      filter(Device == device_a, Date == d) %>%
      select(Genus, n) %>%
      pivot_wider(names_from = Genus, values_from = n, values_fill = 0)
    
    vec_b <- data %>%
      filter(Device == device_b, Date == d) %>%
      select(Genus, n) %>%
      pivot_wider(names_from = Genus, values_from = n, values_fill = 0)
    
    # Align to same genus columns
    all_genera <- union(names(vec_a), names(vec_b))
    vec_a <- vec_a %>% bind_rows(tibble()) %>%
      select(any_of(all_genera)) %>%
      mutate(across(everything(), ~ replace_na(., 0)))
    vec_b <- vec_b %>% bind_rows(tibble()) %>%
      select(any_of(all_genera)) %>%
      mutate(across(everything(), ~ replace_na(., 0)))
    
    # Fill missing genera with 0
    for (g in setdiff(all_genera, names(vec_a))) vec_a[[g]] <- 0
    for (g in setdiff(all_genera, names(vec_b))) vec_b[[g]] <- 0
    vec_a <- vec_a[, all_genera]
    vec_b <- vec_b[, all_genera]
    
    mat <- rbind(as.numeric(vec_a), as.numeric(vec_b))
    
    # Bray-Curtis dissimilarity → convert to similarity
    bc_dissim <- as.numeric(vegdist(mat, method = "bray"))
    bc_sim    <- 1 - bc_dissim
    
    tibble(
      Date       = d,
      Pair       = pair_label,
      Device_A   = device_a,
      Device_B   = device_b,
      BC_similarity   = bc_sim,
      BC_dissimilarity = bc_dissim,
      n_genera_A = sum(as.numeric(vec_a) > 0),
      n_genera_B = sum(as.numeric(vec_b) > 0)
    )
  })
  
  results
}

# --- Pair A: FAIRD1 vs FAIRD2 (Maize) ----------------------------------------
bc_pair_maize <- compute_daily_bc(faird_daily, "FAIRD1", "FAIRD2",
                                  "Intra-Maize (FAIRD1 vs FAIRD2)")

# --- Pair B: FAIRD3 vs FAIRD4 (Meadow) ---------------------------------------
bc_pair_meadow <- compute_daily_bc(faird_daily, "FAIRD3", "FAIRD4",
                                   "Intra-Meadow (FAIRD3 vs FAIRD4)")

# --- Pair C: Cross-habitat pooled -------------------------------------------
# Pool daily counts across FAIRD1+FAIRD2 (Maize) and FAIRD3+FAIRD4 (Meadow)
faird_pooled <- faird_daily %>%
  mutate(Pool = if_else(Ambient == "Maize", "Pool_Maize", "Pool_Meadow")) %>%
  group_by(Pool, Date, Genus) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  rename(Device = Pool)

bc_pair_cross <- compute_daily_bc(faird_pooled,
                                  "Pool_Maize", "Pool_Meadow",
                                  "Cross-habitat (Maize vs Meadow)")

# --- Combine daily results ---------------------------------------------------
bc_daily_all <- bind_rows(bc_pair_maize, bc_pair_meadow, bc_pair_cross)

# Add habitat classification for summary
bc_daily_all <- bc_daily_all %>%
  mutate(Comparison = case_when(
    str_detect(Pair, "Intra-Maize")  ~ "Intra-Maize",
    str_detect(Pair, "Intra-Meadow") ~ "Intra-Meadow",
    str_detect(Pair, "Cross")        ~ "Cross-habitat"
  ))

# --- Summary statistics ------------------------------------------------------
bc_summary <- bc_daily_all %>%
  group_by(Comparison) %>%
  summarise(
    n_days        = n(),
    mean_BC_sim   = mean(BC_similarity, na.rm = TRUE),
    sd_BC_sim     = sd(BC_similarity,   na.rm = TRUE),
    median_BC_sim = median(BC_similarity, na.rm = TRUE),
    min_BC_sim    = min(BC_similarity,  na.rm = TRUE),
    max_BC_sim    = max(BC_similarity,  na.rm = TRUE),
    .groups = "drop"
  )

print("=== Daily Bray-Curtis Similarity Summary ===")
print(bc_summary)

# --- Intra vs Cross-habitat comparison (Wilcoxon) ----------------------------
bc_intra <- bc_daily_all %>%
  filter(Comparison %in% c("Intra-Maize", "Intra-Meadow")) %>%
  pull(BC_similarity)

bc_cross <- bc_daily_all %>%
  filter(Comparison == "Cross-habitat") %>%
  pull(BC_similarity)

wilcox_result <- wilcox.test(bc_intra, bc_cross, exact = FALSE)
cohen_d_val   <- (mean(bc_intra, na.rm=TRUE) - mean(bc_cross, na.rm=TRUE)) /
  sqrt((var(bc_intra, na.rm=TRUE) + var(bc_cross, na.rm=TRUE)) / 2)

cat("\nIntra vs Cross-habitat Bray-Curtis:\n")
cat("  Intra mean:", round(mean(bc_intra, na.rm=TRUE), 3), "\n")
cat("  Cross mean:", round(mean(bc_cross, na.rm=TRUE), 3), "\n")
cat("  Wilcoxon p:", round(wilcox_result$p.value, 4), "\n")
cat("  Cohen's d: ", round(cohen_d_val, 3), "\n")

# --- Maize vs Meadow intra-habitat comparison (Wilcoxon) --------------------
bc_maize_vals  <- bc_daily_all %>% filter(Comparison == "Intra-Maize")  %>% pull(BC_similarity)
bc_meadow_vals <- bc_daily_all %>% filter(Comparison == "Intra-Meadow") %>% pull(BC_similarity)

wilcox_habitat <- wilcox.test(bc_maize_vals, bc_meadow_vals, exact = FALSE)
cohen_d_hab    <- (mean(bc_meadow_vals, na.rm=TRUE) - mean(bc_maize_vals, na.rm=TRUE)) /
  sqrt((var(bc_meadow_vals, na.rm=TRUE) + var(bc_maize_vals, na.rm=TRUE)) / 2)

cat("\nIntra-Maize vs Intra-Meadow Bray-Curtis:\n")
cat("  Maize  mean:", round(mean(bc_maize_vals,  na.rm=TRUE), 3), "\n")
cat("  Meadow mean:", round(mean(bc_meadow_vals, na.rm=TRUE), 3), "\n")
cat("  Wilcoxon p: ", round(wilcox_habitat$p.value, 4), "\n")
cat("  Cohen's d:  ", round(cohen_d_hab, 3), "\n")

# =============================================================================
# ANALYSIS 2 — Pooled full-period Bray-Curtis between habitats
# =============================================================================
# Aggregate all FAIRD detections over entire study period per habitat
# → single abundance vector per habitat → one Bray-Curtis value

faird_pooled_total <- faird_genus %>%
  left_join(device_habitat, by = "Device") %>%
  group_by(Ambient, Genus) %>%
  summarise(n = n(), .groups = "drop")

mat_habitat <- faird_pooled_total %>%
  pivot_wider(names_from = Genus, values_from = n, values_fill = 0) %>%
  column_to_rownames("Ambient")

bc_habitat_pooled <- 1 - as.numeric(vegdist(mat_habitat, method = "bray"))

cat("\n=== Pooled Full-Period Bray-Curtis (Maize vs Meadow) ===\n")
cat("  BC similarity:", round(bc_habitat_pooled, 3), "\n")
cat("  BC dissimilarity:", round(1 - bc_habitat_pooled, 3), "\n")

# Per-device pooled (for heatmap)
mat_per_device <- faird_genus %>%
  left_join(device_habitat, by = "Device") %>%
  group_by(Device, Genus) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = n, values_fill = 0) %>%
  column_to_rownames("Device")

bc_device_dist <- vegdist(mat_per_device, method = "bray")
bc_device_sim  <- 1 - as.matrix(bc_device_dist)

cat("\n=== Full-Period Pairwise Bray-Curtis Similarity Matrix (all FAIRD units) ===\n")
print(round(bc_device_sim, 3))

# =============================================================================
# FIGURES
# =============================================================================

dir.create("outputs", showWarnings = FALSE)

# Colour mapping consistent with project palette
comp_colors <- c(
  "Intra-Maize"   = unname(project_colors["FAIRD"]),   # springgreen4
  "Intra-Meadow"  = unname(project_colors["AMMOD"]),   # dodgerblue
  "Cross-habitat" = "grey50"
)

# --- Figure A: Violin + boxplot of daily BC similarity by comparison type ---
p_violin <- ggplot(bc_daily_all,
                   aes(x = Comparison, y = BC_similarity, fill = Comparison)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = 16, outlier.size = 1.5,
               fill = "white", alpha = 0.8) +
  scale_fill_manual(values = comp_colors) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Daily Bray-Curtis Similarity between FAIRD replicates",
    subtitle = "Abundance-weighted genus-level similarity",
    x        = NULL,
    y        = "Bray-Curtis Similarity"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold"),
        axis.text.x   = element_text(size = 10))

# --- Figure B: Full-period pairwise heatmap ----------------------------------
bc_long <- as.data.frame(bc_device_sim) %>%
  rownames_to_column("Device_A") %>%
  pivot_longer(-Device_A, names_to = "Device_B", values_to = "BC_sim")

p_heatmap <- ggplot(bc_long,
                    aes(x = Device_A, y = Device_B, fill = BC_sim)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(BC_sim, 2)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#a8ddb5", high = "darkgreen",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "BC\nSimilarity") +
  labs(
    title    = "Full-period Bray-Curtis Similarity (FAIRD units)",
    subtitle = "Pooled across entire study period (Aug 23 – Sep 13)",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title  = element_text(face = "bold"),
        axis.text   = element_text(size = 11))

# --- Combined figure ---------------------------------------------------------
p_combined <- p_violin + p_heatmap +
  plot_annotation(tag_levels = "A")

ggsave("outputs/BrayCurtis_combined.png", p_combined,
       width = 12, height = 6, dpi = 300)
ggsave("outputs/BrayCurtis_combined.pdf", p_combined,
       width = 12, height = 6)

# =============================================================================
# SAVE RESULTS
# =============================================================================

write_csv(bc_daily_all,
          "outputs/BrayCurtis_daily_pairs.csv")

write_csv(
  as.data.frame(bc_device_sim) %>% rownames_to_column("Device"),
  "outputs/BrayCurtis_pooled_device_matrix.csv"
)

write_csv(
  tibble(
    Comparison         = "Maize_vs_Meadow_pooled",
    BC_similarity      = bc_habitat_pooled,
    BC_dissimilarity   = 1 - bc_habitat_pooled
  ),
  "outputs/BrayCurtis_pooled_habitats.csv"
)

cat("\n✅ BrayCurtis_FAIRD_Similarity.R completed.\n")
cat("   Outputs saved to outputs/\n")