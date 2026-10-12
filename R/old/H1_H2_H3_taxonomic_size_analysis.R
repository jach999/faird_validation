# ==============================================================================
# TAXONOMIC CONSISTENCY ANALYSIS: INTRA- AND CROSS-HABITAT IZSAK-PRICE SIMILARITY
# ==============================================================================
# Purpose: Evaluate whether devices of the same type detect consistent
#          taxonomic composition within habitats (intra-habitat replicates)
#          and whether composition differs between habitats (cross-habitat).
#          Uses the same Izsak-Price (2001) metric as H3 for methodological
#          consistency.
#
# Pairs analyzed:
#   Intra-habitat:    FAIRD1 vs FAIRD2 (Maize) | FAIRD3 vs FAIRD4 (Meadow)
#                     AMMOD1 vs AMMOD2 (Maize)  | AMMOD3 vs AMMOD4 (Meadow)
#   Cross-habitat:    Pool(FAIRD1+FAIRD2) vs Pool(FAIRD3+FAIRD4)
#                     Pool(AMMOD1+AMMOD2) vs Pool(AMMOD3+AMMOD4)
#   Site-level cross: FAIRD1 vs FAIRD3 | FAIRD1 vs FAIRD4
#                     FAIRD2 vs FAIRD3 | FAIRD2 vs FAIRD4
#                     AMMOD1 vs AMMOD3 | AMMOD1 vs AMMOD4
#                     AMMOD2 vs AMMOD3 | AMMOD2 vs AMMOD4
#
# Statistics: Wilcoxon rank-sum + Cohen's d (intra vs cross per method);
#             site-level cross comparisons to identify which maize site
#             drives intra-maize heterogeneity (edge/core investigation)
#
# Author: Juan Chiavassa
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(vegan)     # taxa2dist()
library(effsize)   # cohen.d()
#install.packages("colorspace")
library(colorspace)

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

source("custom_colors.R")
source("load_data.R")
load("data_cleaned.RData")

# ==============================================================================
# SECTION 1: DATA PREPARATION (same scope as H3)
# ==============================================================================

cat("=== SECTION 1: DATA PREPARATION ===\n\n")

SIZE_THRESHOLD  <- 3  # mm — same as H3
date_start      <- as.Date("2023-08-23")
date_end        <- as.Date("2023-09-13")
date_exclude    <- as.Date("2023-09-04")   # bottle change
ammod4_fail_start <- as.Date("2023-08-28")
ammod4_fail_end   <- as.Date("2023-09-03")

all_dates <- seq.Date(date_start, date_end, by = "day")
all_dates <- all_dates[all_dates != date_exclude]

cat(sprintf("Date range: %s to %s (%d days, Sep 4 excluded)\n",
            date_start, date_end, length(all_dates)))
cat(sprintf("Size threshold: >= %d mm (consistent with H3)\n", SIZE_THRESHOLD))
cat(sprintf("AMMOD4 failure excluded: %s to %s\n\n",
            ammod4_fail_start, ammod4_fail_end))

# Helper: build Taxon_ID from lowest available taxonomic level
build_taxon_id <- function(df) {
  df %>%
    dplyr::mutate(
      Taxon_ID = dplyr::case_when(
        !is.na(Genus)       & Genus       != "" ~ Genus,
        !is.na(Family)      & Family      != "" ~ Family,
        !is.na(Superfamily) & Superfamily != "" ~ Superfamily,
        !is.na(Infraorder)  & Infraorder  != "" ~ Infraorder,
        !is.na(Suborder)    & Suborder    != "" ~ Suborder,
        !is.na(Order)       & Order       != "" ~ Order,
        !is.na(Class)       & Class       != "" ~ Class,
        TRUE ~ "Unknown"
      )
    )
}

# --- FAIRD data ---
faird_data <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD", Date %in% all_dates) %>%
  build_taxon_id() %>%
  dplyr::filter(Taxon_ID != "Unknown")

cat(sprintf("FAIRD records loaded: %d\n", nrow(faird_data)))
cat("FAIRD devices:", paste(sort(unique(faird_data$Device)), collapse = ", "), "\n\n")

# --- AMMOD data ---
ammod_data <- data_list$ammod_taxonomy %>%
  dplyr::filter(
    Date %in% all_dates,
    !is.na(Body_length),
    Body_length >= SIZE_THRESHOLD
  ) %>%
  dplyr::filter(
    !(Device == "AMMOD4" &
        Date >= ammod4_fail_start &
        Date <= ammod4_fail_end)
  ) %>%
  build_taxon_id() %>%
  dplyr::filter(Taxon_ID != "Unknown")

cat(sprintf("AMMOD records loaded (>= %dmm): %d\n", SIZE_THRESHOLD, nrow(ammod_data)))
cat("AMMOD devices:", paste(sort(unique(ammod_data$Device)), collapse = ", "), "\n\n")

# ==============================================================================
# SECTION 2: BUILD TAXONOMIC DISTANCE MATRIX (same method as H3)
# ==============================================================================

cat("=== SECTION 2: TAXONOMIC DISTANCE MATRIX ===\n\n")

tax_cols <- c("Taxon_ID", "Class", "Order", "Suborder", "Infraorder",
              "Superfamily", "Family", "Subfamily", "Genus", "Species")

faird_taxa <- faird_data %>%
  dplyr::select(dplyr::all_of(tax_cols)) %>%
  dplyr::distinct()

ammod_taxa <- ammod_data %>%
  dplyr::select(dplyr::all_of(tax_cols)) %>%
  dplyr::distinct()

all_taxa <- dplyr::bind_rows(
  faird_taxa %>% dplyr::mutate(Source = "FAIRD"),
  ammod_taxa %>% dplyr::mutate(Source = "AMMOD")
) %>%
  dplyr::distinct(Taxon_ID, .keep_all = TRUE)

cat(sprintf("Total unique taxa in distance matrix: %d\n", nrow(all_taxa)))
cat(sprintf("  FAIRD: %d | AMMOD: %d | Shared: %d\n\n",
            nrow(faird_taxa), nrow(ammod_taxa),
            length(intersect(faird_taxa$Taxon_ID, ammod_taxa$Taxon_ID))))

taxonomy_table <- all_taxa %>%
  dplyr::select(-Source, -Taxon_ID) %>%
  as.data.frame()

rownames(taxonomy_table) <- all_taxa$Taxon_ID

taxon_distances <- vegan::taxa2dist(taxonomy_table, varstep = TRUE, check = TRUE)
dist_matrix     <- as.matrix(taxon_distances)

cat(sprintf("Distance matrix: %d x %d | range [%.4f, %.4f]\n\n",
            nrow(dist_matrix), ncol(dist_matrix),
            min(taxon_distances), max(taxon_distances)))

# ==============================================================================
# SECTION 3: IZSAK-PRICE FUNCTION (identical to H3)
# ==============================================================================

calculate_taxonomic_similarity <- function(taxa_A, taxa_B, dist_matrix) {
  taxa_A <- taxa_A[taxa_A %in% rownames(dist_matrix)]
  taxa_B <- taxa_B[taxa_B %in% rownames(dist_matrix)]
  if (length(taxa_A) == 0 | length(taxa_B) == 0) return(NA)
  
  min_A_to_B <- sapply(taxa_A, function(a) {
    if (a %in% taxa_B) 0L
    else min(dist_matrix[a, taxa_B], na.rm = TRUE)
  })
  
  min_B_to_A <- sapply(taxa_B, function(b) {
    if (b %in% taxa_A) 0L
    else min(dist_matrix[b, taxa_A], na.rm = TRUE)
  })
  
  avg_min  <- mean(c(min_A_to_B, min_B_to_A))
  max_dist <- max(dist_matrix, na.rm = TRUE)
  1 - (avg_min / max_dist)
}

# ==============================================================================
# SECTION 4: DEFINE PAIRS AND CALCULATE DAILY SIMILARITY
# ==============================================================================

cat("=== SECTION 4: DAILY IZSAK-PRICE SIMILARITY BY PAIR ===\n\n")

# Helper: get taxa vector for a device on a date
get_taxa <- function(data, device, date) {
  data %>%
    dplyr::filter(Device == device, Date == date) %>%
    dplyr::pull(Taxon_ID) %>%
    unique()
}

# Helper: get taxa vector for a pool of devices on a date
get_taxa_pool <- function(data, devices, date) {
  data %>%
    dplyr::filter(Device %in% devices, Date == date) %>%
    dplyr::pull(Taxon_ID) %>%
    unique()
}

# Pair definitions
pairs <- list(
  # --- Intra-habitat ---
  list(label = "FAIRD1 vs FAIRD2", type = "Intra-habitat",
       habitat = "Maize", method = "FAIRD",
       data = faird_data, mode = "single",
       dev_A = "FAIRD1", dev_B = "FAIRD2"),
  
  list(label = "FAIRD3 vs FAIRD4", type = "Intra-habitat",
       habitat = "Meadow", method = "FAIRD",
       data = faird_data, mode = "single",
       dev_A = "FAIRD3", dev_B = "FAIRD4"),
  
  list(label = "AMMOD1 vs AMMOD2", type = "Intra-habitat",
       habitat = "Maize", method = "AMMOD",
       data = ammod_data, mode = "single",
       dev_A = "AMMOD1", dev_B = "AMMOD2"),
  
  list(label = "AMMOD3 vs AMMOD4", type = "Intra-habitat",
       habitat = "Meadow", method = "AMMOD",
       data = ammod_data, mode = "single",
       dev_A = "AMMOD3", dev_B = "AMMOD4"),
  
  # --- Cross-habitat (pooled) ---
  list(label = "FAIRD Maize vs Meadow", type = "Cross-habitat",
       habitat = "Cross", method = "FAIRD",
       data = faird_data, mode = "pool",
       pool_A = c("FAIRD1", "FAIRD2"),
       pool_B = c("FAIRD3", "FAIRD4")),
  
  list(label = "AMMOD Maize vs Meadow", type = "Cross-habitat",
       habitat = "Cross", method = "AMMOD",
       data = ammod_data, mode = "pool",
       pool_A = c("AMMOD1", "AMMOD2"),
       pool_B = c("AMMOD3", "AMMOD4")),
  
  # --- Site-level cross (individual maize site vs individual meadow site) ---
  # FAIRD
  list(label = "FAIRD1 vs FAIRD3", type = "Site-level cross",
       habitat = "S1-S3", method = "FAIRD",
       data = faird_data, mode = "single",
       dev_A = "FAIRD1", dev_B = "FAIRD3"),
  
  list(label = "FAIRD1 vs FAIRD4", type = "Site-level cross",
       habitat = "S1-S4", method = "FAIRD",
       data = faird_data, mode = "single",
       dev_A = "FAIRD1", dev_B = "FAIRD4"),
  
  list(label = "FAIRD2 vs FAIRD3", type = "Site-level cross",
       habitat = "S2-S3", method = "FAIRD",
       data = faird_data, mode = "single",
       dev_A = "FAIRD2", dev_B = "FAIRD3"),
  
  list(label = "FAIRD2 vs FAIRD4", type = "Site-level cross",
       habitat = "S2-S4", method = "FAIRD",
       data = faird_data, mode = "single",
       dev_A = "FAIRD2", dev_B = "FAIRD4"),
  
  # AMMOD
  list(label = "AMMOD1 vs AMMOD3", type = "Site-level cross",
       habitat = "S1-S3", method = "AMMOD",
       data = ammod_data, mode = "single",
       dev_A = "AMMOD1", dev_B = "AMMOD3"),
  
  list(label = "AMMOD1 vs AMMOD4", type = "Site-level cross",
       habitat = "S1-S4", method = "AMMOD",
       data = ammod_data, mode = "single",
       dev_A = "AMMOD1", dev_B = "AMMOD4"),
  
  list(label = "AMMOD2 vs AMMOD3", type = "Site-level cross",
       habitat = "S2-S3", method = "AMMOD",
       data = ammod_data, mode = "single",
       dev_A = "AMMOD2", dev_B = "AMMOD3"),
  
  list(label = "AMMOD2 vs AMMOD4", type = "Site-level cross",
       habitat = "S2-S4", method = "AMMOD",
       data = ammod_data, mode = "single",
       dev_A = "AMMOD2", dev_B = "AMMOD4")
)

# Calculate daily similarity for each pair
all_results <- purrr::map_dfr(pairs, function(p) {
  cat(sprintf("Processing: %s\n", p$label))
  
  purrr::map_dfr(all_dates, function(d) {
    if (p$mode == "single") {
      # Skip days where dev_B is AMMOD4 and within failure period
      if (p$dev_B == "AMMOD4" &&
          d >= ammod4_fail_start && d <= ammod4_fail_end) {
        return(NULL)
      }
      taxa_A <- get_taxa(p$data, p$dev_A, d)
      taxa_B <- get_taxa(p$data, p$dev_B, d)
    } else {
      # Pool mode: exclude AMMOD4 from pool_B during failure period
      pool_B_active <- if ("AMMOD4" %in% p$pool_B &&
                           d >= ammod4_fail_start && d <= ammod4_fail_end) {
        p$pool_B[p$pool_B != "AMMOD4"]
      } else {
        p$pool_B
      }
      # If pool_B is now empty after exclusion, skip this day
      if (length(pool_B_active) == 0) return(NULL)
      
      taxa_A <- get_taxa_pool(p$data, p$pool_A, d)
      taxa_B <- get_taxa_pool(p$data, pool_B_active, d)
    }
    
    delta_s <- calculate_taxonomic_similarity(taxa_A, taxa_B, dist_matrix)
    
    if (is.na(delta_s)) return(NULL)
    
    tibble::tibble(
      Date    = d,
      Pair    = p$label,
      Type    = p$type,
      Habitat = p$habitat,
      Method  = p$method,
      Delta_S = delta_s,
      n_A     = length(taxa_A),
      n_B     = length(taxa_B)
    )
  })
})

cat(sprintf("\nTotal daily observations computed: %d\n\n", nrow(all_results)))

# Summary per pair
pair_summary <- all_results %>%
  dplyr::group_by(Pair, Type, Habitat, Method) %>%
  dplyr::summarise(
    n_days    = dplyr::n(),
    mean_DS   = round(mean(Delta_S, na.rm = TRUE), 4),
    sd_DS     = round(sd(Delta_S, na.rm = TRUE), 4),
    median_DS = round(median(Delta_S, na.rm = TRUE), 4),
    min_DS    = round(min(Delta_S, na.rm = TRUE), 4),
    max_DS    = round(max(Delta_S, na.rm = TRUE), 4),
    .groups   = "drop"
  )

cat("=== SUMMARY PER PAIR ===\n\n")
print(pair_summary, n = Inf)
cat("\n")

readr::write_csv(pair_summary, "outputs/TaxConsistency_Summary.csv")
readr::write_csv(all_results,  "outputs/TaxConsistency_DailyResults.csv")
cat("[OK] outputs/TaxConsistency_Summary.csv\n")
cat("[OK] outputs/TaxConsistency_DailyResults.csv\n\n")

# ==============================================================================
# SECTION 5: STATISTICAL TESTS — INTRA vs CROSS HABITAT PER METHOD
# ==============================================================================

cat("=== SECTION 5: WILCOXON + COHEN'S D (Intra vs Cross per Method) ===\n\n")
cat("Comparison: Are intra-habitat replicates more similar than cross-habitat?\n\n")

for (method in c("FAIRD", "AMMOD")) {
  cat(sprintf("--- %s ---\n", method))
  
  intra_vals <- all_results %>%
    dplyr::filter(Method == method, Type == "Intra-habitat") %>%
    dplyr::pull(Delta_S)
  
  cross_vals <- all_results %>%
    dplyr::filter(Method == method, Type == "Cross-habitat") %>%
    dplyr::pull(Delta_S)
  
  if (length(intra_vals) < 3 || length(cross_vals) < 3) {
    cat("  Insufficient data for test\n\n")
    next
  }
  
  wt <- stats::wilcox.test(intra_vals, cross_vals, exact = FALSE)
  cd <- effsize::cohen.d(intra_vals, cross_vals)
  
  cat(sprintf("  Intra-habitat: n=%d | mean=%.4f | median=%.4f\n",
              length(intra_vals), mean(intra_vals), median(intra_vals)))
  cat(sprintf("  Cross-habitat: n=%d | mean=%.4f | median=%.4f\n",
              length(cross_vals), mean(cross_vals), median(cross_vals)))
  cat(sprintf("  Wilcoxon W = %.1f | p = %.4f%s\n",
              wt$statistic, wt$p.value,
              ifelse(wt$p.value < 0.05, " *", " ns")))
  cat(sprintf("  Cohen's d = %.3f (%s)\n\n",
              cd$estimate, cd$magnitude))
}

# ==============================================================================
# SECTION 5B: STATISTICAL TESTS — INTRA-MAIZE vs INTRA-MEADOW PER METHOD
# ==============================================================================

cat("=== SECTION 5B: WILCOXON + COHEN'S D (Intra-Maize vs Intra-Meadow) ===\n\n")
cat("Comparison: Does within-habitat taxonomic similarity differ between habitat types?\n")
cat("Pairs: FAIRD1 vs FAIRD2 (Maize) vs FAIRD3 vs FAIRD4 (Meadow)\n")
cat("       AMMOD1 vs AMMOD2 (Maize) vs AMMOD3 vs AMMOD4 (Meadow)\n\n")

intra_pairs <- list(
  list(method = "FAIRD",
       maize_label  = "FAIRD1 vs FAIRD2",
       meadow_label = "FAIRD3 vs FAIRD4"),
  list(method = "AMMOD",
       maize_label  = "AMMOD1 vs AMMOD2",
       meadow_label = "AMMOD3 vs AMMOD4")
)

for (p in intra_pairs) {
  cat(sprintf("--- %s ---\n", p$method))
  
  maize_vals <- all_results %>%
    dplyr::filter(Pair == p$maize_label) %>%
    dplyr::pull(Delta_S)
  
  meadow_vals <- all_results %>%
    dplyr::filter(Pair == p$meadow_label) %>%
    dplyr::pull(Delta_S)
  
  if (length(maize_vals) < 3 || length(meadow_vals) < 3) {
    cat("  Insufficient data for test\n\n")
    next
  }
  
  wt <- stats::wilcox.test(maize_vals, meadow_vals, exact = FALSE)
  cd <- effsize::cohen.d(maize_vals, meadow_vals)
  
  cat(sprintf("  Intra-Maize  (%s): n=%d | mean=%.4f | median=%.4f\n",
              p$maize_label,  length(maize_vals),
              mean(maize_vals),  median(maize_vals)))
  cat(sprintf("  Intra-Meadow (%s): n=%d | mean=%.4f | median=%.4f\n",
              p$meadow_label, length(meadow_vals),
              mean(meadow_vals), median(meadow_vals)))
  cat(sprintf("  Wilcoxon W = %.1f | p = %.4f%s\n",
              wt$statistic, wt$p.value,
              ifelse(wt$p.value < 0.05, " *", " ns")))
  cat(sprintf("  Cohen's d = %.3f (%s)\n\n",
              cd$estimate, cd$magnitude))
}



# ==============================================================================
# SECTION 5C: SITE-LEVEL CROSS COMPARISONS — EDGE/CORE INVESTIGATION
# ==============================================================================

cat("=== SECTION 5C: SITE-LEVEL CROSS COMPARISONS (Edge/Core Investigation) ===\n\n")
cat("Question: Which maize site (Site1 or Site2) is more similar to meadow sites?\n")
cat("If Site1 drives intra-maize heterogeneity, Site1 should show lower ΔS\n")
cat("with meadow sites than Site2 (or vice versa).\n\n")

site_cross_summary <- all_results %>%
  dplyr::filter(Type == "Site-level cross") %>%
  dplyr::group_by(Pair, Method, Habitat) %>%
  dplyr::summarise(
    n_days    = dplyr::n(),
    mean_DS   = round(mean(Delta_S, na.rm = TRUE), 4),
    sd_DS     = round(sd(Delta_S, na.rm = TRUE), 4),
    median_DS = round(median(Delta_S, na.rm = TRUE), 4),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(Method, Pair)

cat("--- Summary: mean ΔS per site-level cross pair ---\n\n")
print(site_cross_summary, n = Inf)
cat("\n")

# For each method: compare Site1-meadow vs Site2-meadow similarity
for (method in c("FAIRD", "AMMOD")) {
  cat(sprintf("--- %s: Site1 vs meadow compared to Site2 vs meadow ---\n", method))
  
  dev1  <- if (method == "FAIRD") "FAIRD1" else "AMMOD1"
  dev2  <- if (method == "FAIRD") "FAIRD2" else "AMMOD2"
  dev3  <- if (method == "FAIRD") "FAIRD3" else "AMMOD3"
  dev4  <- if (method == "FAIRD") "FAIRD4" else "AMMOD4"
  
  lab_1_3 <- sprintf("%s vs %s", dev1, dev3)
  lab_1_4 <- sprintf("%s vs %s", dev1, dev4)
  lab_2_3 <- sprintf("%s vs %s", dev2, dev3)
  lab_2_4 <- sprintf("%s vs %s", dev2, dev4)
  
  s1_meadow <- all_results %>%
    dplyr::filter(Pair %in% c(lab_1_3, lab_1_4)) %>%
    dplyr::pull(Delta_S)
  
  s2_meadow <- all_results %>%
    dplyr::filter(Pair %in% c(lab_2_3, lab_2_4)) %>%
    dplyr::pull(Delta_S)
  
  if (length(s1_meadow) < 3 || length(s2_meadow) < 3) {
    cat("  Insufficient data\n\n"); next
  }
  
  wt <- stats::wilcox.test(s1_meadow, s2_meadow, exact = FALSE)
  cd <- effsize::cohen.d(s1_meadow, s2_meadow)
  
  cat(sprintf("  Site1 vs meadow (pairs: %s, %s): n=%d | mean=%.4f | median=%.4f\n",
              lab_1_3, lab_1_4, length(s1_meadow),
              mean(s1_meadow), median(s1_meadow)))
  cat(sprintf("  Site2 vs meadow (pairs: %s, %s): n=%d | mean=%.4f | median=%.4f\n",
              lab_2_3, lab_2_4, length(s2_meadow),
              mean(s2_meadow), median(s2_meadow)))
  cat(sprintf("  Wilcoxon W = %.1f | p = %.4f%s\n",
              wt$statistic, wt$p.value,
              ifelse(wt$p.value < 0.05, " *", " ns")))
  cat(sprintf("  Cohen's d = %.3f (%s)\n", cd$estimate, cd$magnitude))
  cat(sprintf("  Interpretation: %s is more similar to meadow\n\n",
              ifelse(mean(s1_meadow) > mean(s2_meadow), "Site1", "Site2")))
}

# ==============================================================================
# SECTION 6: PLOTS
# ==============================================================================

cat("=== SECTION 6: PLOTS ===\n\n")

# --- 6A: Time series of daily Delta_S per pair ---

# Color logic: habitat color for intra-habitat pairs, method color for cross-habitat
pair_colors <- c(
  "FAIRD1 vs FAIRD2"      = unname(custom_colors["Maize"]),
  "FAIRD3 vs FAIRD4"      = unname(custom_colors["Meadow"]),
  "AMMOD1 vs AMMOD2"      = unname(custom_colors["Maize"]),
  "AMMOD3 vs AMMOD4"      = unname(custom_colors["Meadow"]),
  "FAIRD Maize vs Meadow" = unname(custom_colors["FAIRD"]),
  "AMMOD Maize vs Meadow" = unname(custom_colors["AMMOD"])
)

plot_ts <- ggplot(
  all_results,
  aes(x = Date, y = Delta_S, color = Pair, linetype = Type)
) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 1.8, alpha = 0.8) +
  facet_wrap(~ Method, ncol = 1) +
  scale_color_manual(values = pair_colors) +
  scale_linetype_manual(values = c("Intra-habitat" = "solid",
                                   "Cross-habitat"  = "dashed")) +
  labs(
    title    = "Daily Taxonomic Similarity (Izsak-Price ΔS) by Pair",
    subtitle = sprintf("Size threshold: >= %d mm | Same metric as H3", SIZE_THRESHOLD),
    x        = "Date",
    y        = "Izsak-Price ΔS",
    color    = "Pair",
    linetype = "Comparison type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "gray40"),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold")
  )

# --- 6B: Violin plot — Intra vs Cross per Method ---
plot_violin <- ggplot(
  all_results,
  aes(x = Type, y = Delta_S, fill = Method)
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = 16,
               outlier.size = 2, alpha = 0.8) +
  facet_wrap(~ Method, ncol = 2) +
  scale_fill_manual(values = c(
    "FAIRD" = unname(custom_colors["FAIRD"]),
    "AMMOD" = unname(custom_colors["AMMOD"])
  )) +
  labs(
    title = "Taxonomic Similarity: Intra- vs Cross-Habitat",
    x     = "Comparison type",
    y     = "Izsak-Price ΔS",
    fill  = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    strip.text       = element_text(face = "bold")
  )

# --- 6C: Violin plot — Intra-Maize vs Intra-Meadow per Method ---
intra_habitat_data <- all_results %>%
  dplyr::filter(Type == "Intra-habitat") %>%
  dplyr::mutate(
    Habitat_label = dplyr::case_when(
      Pair %in% c("FAIRD1 vs FAIRD2", "AMMOD1 vs AMMOD2") ~ "Maize",
      Pair %in% c("FAIRD3 vs FAIRD4", "AMMOD3 vs AMMOD4") ~ "Meadow",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Habitat_label))

plot_violin_habitat <- ggplot(
  intra_habitat_data,
  aes(x = Habitat_label, y = Delta_S, fill = Method)
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = 16,
               outlier.size = 2, alpha = 0.8) +
  facet_wrap(~ Method, ncol = 2) +
  scale_fill_manual(values = c(
    "FAIRD" = unname(custom_colors["FAIRD"]),
    "AMMOD" = unname(custom_colors["AMMOD"])
  )) +
  labs(
    title    = "Intra-Habitat Similarity: Maize vs Meadow",
    subtitle = "Intra-habitat pairs only",
    x        = "Habitat",
    y        = "Izsak-Price ΔS",
    fill     = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "gray40"),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    strip.text       = element_text(face = "bold")
  )

# --- 6D: Heatmap — mean ΔS between all site pairs (AMMOD only, unbiased) ---
heatmap_data <- all_results %>%
  dplyr::filter(Method == "AMMOD") %>%
  dplyr::group_by(Pair, Type) %>%
  dplyr::summarise(mean_DS = mean(Delta_S, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    Site_A = stringr::str_extract(Pair, "AMMOD\\d"),
    Site_B = stringr::str_extract(Pair, "AMMOD\\d$")
  ) %>%
  dplyr::filter(!is.na(Site_A), !is.na(Site_B)) %>%
  # Mirror the matrix (symmetric)
  dplyr::bind_rows(
    dplyr::rename(., Site_A = Site_B, Site_B = Site_A)
  ) %>%
  dplyr::distinct(Site_A, Site_B, .keep_all = TRUE) %>%
  # Add self-similarity = 1
  dplyr::bind_rows(
    tibble::tibble(
      Site_A = c("AMMOD1","AMMOD2","AMMOD3","AMMOD4"),
      Site_B = c("AMMOD1","AMMOD2","AMMOD3","AMMOD4"),
      mean_DS = 1, Type = "Self"
    )
  ) %>%
  dplyr::mutate(
    Site_A = factor(Site_A, levels = c("AMMOD1","AMMOD2","AMMOD3","AMMOD4")),
    Site_B = factor(Site_B, levels = c("AMMOD4","AMMOD3","AMMOD2","AMMOD1"))
  )

plot_heatmap <- ggplot(heatmap_data,
                       aes(x = Site_A, y = Site_B, fill = mean_DS)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(Site_A == Site_B, "—",
                               sprintf("%.3f", mean_DS))),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#1a3a5c", "#2E86AB", "#A8DADC", "#F1FAEE"),
    limits = c(0.60, 1.0),
    name   = "Mean ΔS"
  ) +
  scale_x_discrete(labels = c("AMMOD1\n(Maize S1)", "AMMOD2\n(Maize S2)",
                              "AMMOD3\n(Meadow S3)", "AMMOD4\n(Meadow S4)")) +
  scale_y_discrete(labels = c("AMMOD4\n(Meadow S4)", "AMMOD3\n(Meadow S3)",
                              "AMMOD2\n(Maize S2)", "AMMOD1\n(Maize S1)")) +
  labs(
    title    = "Pairwise Taxonomic Similarity Between Sites (AMMOD)",
    subtitle = "Mean Izsak-Price ΔS | Higher = more similar community",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "gray40"),
    axis.text        = element_text(size = 9),
    panel.grid       = element_blank(),
    legend.position  = "right"
  )

# --- Combined ---
combined_plot <- (plot_ts + plot_violin + plot_violin_habitat) / plot_heatmap +
  patchwork::plot_layout(heights = c(2, 1))

print(combined_plot)

ggsave("outputs/TaxConsistency_Combined.png",
       plot = combined_plot, width = 18, height = 14,
       dpi = 300, bg = "white")
ggsave("outputs/TaxConsistency_Combined.pdf",
       plot = combined_plot, width = 18, height = 14,
       device = cairo_pdf)

cat("[OK] outputs/TaxConsistency_Combined.png/.pdf\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("=== ANALYSIS COMPLETE ===\n\n")
cat("INTERPRETATION SCOPE:\n")
cat("  ΔS reflects similarity in the DETECTED TAXON POOL per device per day.\n")
cat("  AMMOD counts are metabarcoding-based (presence only, not abundance).\n")
cat("  Results are exploratory and hypothesis-generating.\n\n")
cat("Files generated:\n")
cat("  [OK] outputs/TaxConsistency_Summary.csv\n")
cat("  [OK] outputs/TaxConsistency_DailyResults.csv\n")
cat("  [OK] outputs/TaxConsistency_Combined.png/.pdf\n\n")

# ==============================================================================
# SECTION 7: SUMMARY TABLES — BODY SIZE AND TAXONOMIC SIMILARITY BY SITE,
#            DEVICE TYPE, AND HABITAT
# ==============================================================================

cat("=== SECTION 7: SUMMARY TABLES ===\n\n")
cat("Metrics reported:\n")
cat("  - Mean body length (mm) +/- SD: AMMOD only, family-level deduplication\n")
cat("    (Söhlström et al. 2018; reflects detected family pool, not abundance)\n")
cat("  - Delta BL vs meadow: maize site mean BL minus mean BL of both meadow sites\n")
cat("  - Intra-Habitat Taxonomy DeltaS +/- SD: within-habitat pair, daily values\n")
cat("    (FAIRD: intra-habitat pair; AMMOD: intra-habitat pair)\n")
cat("  - Taxonomy DeltaS vs meadow +/- SD: mean daily DeltaS of each maize site\n")
cat("    against both meadow sites combined (AMMOD only; site-level cross pairs)\n\n")

# ------------------------------------------------------------------------------
# 7.1  Body size — daily mean BL per site
#   AMMOD: family-level deduplication (Söhlström pool; each family = 1 obs/day)
#   FAIRD: abundance-weighted mean BL (each row = 1 individual detection, no dedup)
# ------------------------------------------------------------------------------

# AMMOD — family-deduplicated
daily_bl_ammod <- ammod_data %>%
  dplyr::filter(!is.na(Body_length), !is.na(Family), Family != "") %>%
  dplyr::distinct(Date, Device, Family, .keep_all = TRUE) %>%
  dplyr::mutate(
    Site = dplyr::case_when(
      Device == "AMMOD1" ~ "Site1",
      Device == "AMMOD2" ~ "Site2",
      Device == "AMMOD3" ~ "Site3",
      Device == "AMMOD4" ~ "Site4"
    )
  ) %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    mean_bl = mean(Body_length, na.rm = TRUE),
    n_obs   = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Habitat = dplyr::if_else(Site %in% c("Site1", "Site2"), "Maize", "Meadow"),
    Method  = "AMMOD"
  )

# FAIRD — abundance-weighted (each row = 1 individual detection)
daily_bl_faird <- faird_data %>%
  dplyr::filter(!is.na(Body_length)) %>%
  dplyr::mutate(
    Site = dplyr::case_when(
      Device == "FAIRD1" ~ "Site1",
      Device == "FAIRD2" ~ "Site2",
      Device == "FAIRD3" ~ "Site3",
      Device == "FAIRD4" ~ "Site4"
    )
  ) %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    mean_bl = mean(Body_length, na.rm = TRUE),
    n_obs   = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Habitat = dplyr::if_else(Site %in% c("Site1", "Site2"), "Maize", "Meadow"),
    Method  = "FAIRD"
  )

cat(sprintf("AMMOD daily BL observations: %d | FAIRD daily BL observations: %d\n",
            nrow(daily_bl_ammod), nrow(daily_bl_faird)))

# Site-level summary per method
site_bl_ammod <- daily_bl_ammod %>%
  dplyr::group_by(Site, Habitat, Method) %>%
  dplyr::summarise(
    n_days_bl = dplyr::n(),
    mean_BL   = mean(mean_bl, na.rm = TRUE),
    sd_BL     = sd(mean_bl,   na.rm = TRUE),
    .groups   = "drop"
  )

site_bl_faird <- daily_bl_faird %>%
  dplyr::group_by(Site, Habitat, Method) %>%
  dplyr::summarise(
    n_days_bl = dplyr::n(),
    mean_BL   = mean(mean_bl, na.rm = TRUE),
    sd_BL     = sd(mean_bl,   na.rm = TRUE),
    .groups   = "drop"
  )

site_bl <- dplyr::bind_rows(site_bl_ammod, site_bl_faird)

cat("\nSite-level mean body length by method:\n")
print(site_bl, n = Inf)
cat("\n")

# Delta BL vs meadow per method: maize site mean BL minus grand mean of meadow sites
meadow_mean_ammod <- site_bl_ammod %>%
  dplyr::filter(Habitat == "Meadow") %>%
  dplyr::summarise(grand = mean(mean_BL)) %>% dplyr::pull(grand)

meadow_mean_faird <- site_bl_faird %>%
  dplyr::filter(Habitat == "Meadow") %>%
  dplyr::summarise(grand = mean(mean_BL)) %>% dplyr::pull(grand)

site_bl <- site_bl %>%
  dplyr::mutate(
    meadow_ref = dplyr::if_else(Method == "AMMOD", meadow_mean_ammod, meadow_mean_faird),
    delta_BL_vs_meadow = dplyr::if_else(
      Habitat == "Maize",
      round(mean_BL - meadow_ref, 3),
      NA_real_
    )
  ) %>%
  dplyr::select(-meadow_ref)

# ------------------------------------------------------------------------------
# 7.2  Intra-habitat taxonomy DeltaS — assign pair values to component sites
# ------------------------------------------------------------------------------

# Map each site to its intra-habitat pair label (per method)
site_to_pair <- tibble::tibble(
  Site   = c("Site1","Site2","Site3","Site4","Site1","Site2","Site3","Site4"),
  Method = c("AMMOD","AMMOD","AMMOD","AMMOD","FAIRD","FAIRD","FAIRD","FAIRD"),
  Pair   = c("AMMOD1 vs AMMOD2","AMMOD1 vs AMMOD2",
             "AMMOD3 vs AMMOD4","AMMOD3 vs AMMOD4",
             "FAIRD1 vs FAIRD2","FAIRD1 vs FAIRD2",
             "FAIRD3 vs FAIRD4","FAIRD3 vs FAIRD4")
)

intra_ds_vals <- all_results %>%
  dplyr::filter(Type == "Intra-habitat") %>%
  dplyr::group_by(Pair, Method) %>%
  dplyr::summarise(
    mean_intra_DS = mean(Delta_S, na.rm = TRUE),
    sd_intra_DS   = sd(Delta_S,   na.rm = TRUE),
    n_intra       = dplyr::n(),
    .groups       = "drop"
  )

site_intra_ds <- site_to_pair %>%
  dplyr::left_join(intra_ds_vals, by = c("Pair", "Method"))

# ------------------------------------------------------------------------------
# 7.3  Taxonomy DeltaS vs meadow — per maize site, AMMOD only
# ------------------------------------------------------------------------------

# Site1 vs meadow = daily values from pairs AMMOD1 vs AMMOD3 and AMMOD1 vs AMMOD4
# Site2 vs meadow = daily values from pairs AMMOD2 vs AMMOD3 and AMMOD2 vs AMMOD4

cross_meadow_ammod <- all_results %>%
  dplyr::filter(Type == "Site-level cross", Method == "AMMOD") %>%
  dplyr::mutate(
    Maize_site = dplyr::case_when(
      Pair %in% c("AMMOD1 vs AMMOD3", "AMMOD1 vs AMMOD4") ~ "Site1",
      Pair %in% c("AMMOD2 vs AMMOD3", "AMMOD2 vs AMMOD4") ~ "Site2"
    )
  ) %>%
  dplyr::filter(!is.na(Maize_site)) %>%
  dplyr::group_by(Maize_site) %>%
  dplyr::summarise(
    mean_cross_DS = mean(Delta_S, na.rm = TRUE),
    sd_cross_DS   = sd(Delta_S,   na.rm = TRUE),
    n_cross       = dplyr::n(),
    .groups       = "drop"
  ) %>%
  dplyr::rename(Site = Maize_site)

cross_meadow_faird <- all_results %>%
  dplyr::filter(Type == "Site-level cross", Method == "FAIRD") %>%
  dplyr::mutate(
    Maize_site = dplyr::case_when(
      Pair %in% c("FAIRD1 vs FAIRD3", "FAIRD1 vs FAIRD4") ~ "Site1",
      Pair %in% c("FAIRD2 vs FAIRD3", "FAIRD2 vs FAIRD4") ~ "Site2"
    )
  ) %>%
  dplyr::filter(!is.na(Maize_site)) %>%
  dplyr::group_by(Maize_site) %>%
  dplyr::summarise(
    mean_cross_DS = mean(Delta_S, na.rm = TRUE),
    sd_cross_DS   = sd(Delta_S,   na.rm = TRUE),
    n_cross       = dplyr::n(),
    .groups       = "drop"
  ) %>%
  dplyr::rename(Site = Maize_site)

# ==============================================================================
# TABLE 1: BY SITE
# ==============================================================================

cat("--- TABLE 1: BY SITE ---\n\n")

table1_ammod <- site_intra_ds %>%
  dplyr::filter(Method == "AMMOD") %>%
  dplyr::left_join(
    site_bl %>% dplyr::filter(Method == "AMMOD") %>%
      dplyr::select(Site, mean_BL, sd_BL, delta_BL_vs_meadow),
    by = "Site"
  ) %>%
  dplyr::left_join(
    cross_meadow_ammod %>% dplyr::select(Site, mean_cross_DS, sd_cross_DS, n_cross),
    by = "Site"
  ) %>%
  dplyr::mutate(
    Habitat      = dplyr::if_else(Site %in% c("Site1","Site2"), "Maize", "Meadow"),
    Method_label = "AMMOD",
    BL_report    = sprintf("%.3f +/- %.3f", mean_BL, sd_BL),
    dBL_report   = dplyr::if_else(!is.na(delta_BL_vs_meadow),
                                  sprintf("%+.3f", delta_BL_vs_meadow), "—"),
    intra_DS_report = sprintf("%.3f +/- %.3f (n=%d)", mean_intra_DS, sd_intra_DS, n_intra),
    cross_DS_report = dplyr::if_else(!is.na(mean_cross_DS),
                                     sprintf("%.3f +/- %.3f (n=%d)", mean_cross_DS,
                                             sd_cross_DS, n_cross), "—")
  ) %>%
  dplyr::select(Site, Habitat, Method_label,
                BL_report, dBL_report,
                intra_DS_report, cross_DS_report)

table1_faird <- site_intra_ds %>%
  dplyr::filter(Method == "FAIRD") %>%
  dplyr::left_join(
    site_bl %>% dplyr::filter(Method == "FAIRD") %>%
      dplyr::select(Site, mean_BL, sd_BL, delta_BL_vs_meadow),
    by = "Site"
  ) %>%
  dplyr::left_join(
    cross_meadow_faird %>% dplyr::select(Site, mean_cross_DS, sd_cross_DS, n_cross),
    by = "Site"
  ) %>%
  dplyr::mutate(
    Habitat      = dplyr::if_else(Site %in% c("Site1","Site2"), "Maize", "Meadow"),
    Method_label = "FAIRD",
    BL_report    = sprintf("%.3f +/- %.3f", mean_BL, sd_BL),
    dBL_report   = dplyr::if_else(!is.na(delta_BL_vs_meadow),
                                  sprintf("%+.3f", delta_BL_vs_meadow), "—"),
    intra_DS_report = sprintf("%.3f +/- %.3f (n=%d)", mean_intra_DS, sd_intra_DS, n_intra),
    cross_DS_report = dplyr::if_else(!is.na(mean_cross_DS),
                                     sprintf("%.3f +/- %.3f (n=%d)", mean_cross_DS,
                                             sd_cross_DS, n_cross), "—")
  ) %>%
  dplyr::select(Site, Habitat, Method_label,
                BL_report, dBL_report,
                intra_DS_report, cross_DS_report)

table1 <- dplyr::bind_rows(table1_ammod, table1_faird) %>%
  dplyr::arrange(Site, Method_label) %>%
  dplyr::rename(
    Method                     = Method_label,
    `Mean BL +/- SD (mm)`      = BL_report,
    `Delta BL vs meadow`       = dBL_report,
    `Intra-Habitat Tax. DeltaS`= intra_DS_report,
    `Tax. DeltaS vs meadow`    = cross_DS_report
  )

print(table1, n = Inf)
cat("\nNote: Intra-Habitat Tax. DeltaS is a pair-level metric — Site1 and Site2 share\n")
cat("the same value (they form the intra-maize pair); same applies to Site3 and Site4.\n")
cat("AMMOD BL: family-deduplicated daily mean (Soehlstroem et al. 2018).\n")
cat("FAIRD BL: abundance-weighted daily mean (each detection = 1 individual).\n")
cat("Delta BL vs meadow = site mean BL minus grand mean BL of Sites 3+4 (per method).\n\n")

readr::write_csv(table1, "outputs/CommunityStructure_Table1_BySite.csv")
cat("[OK] outputs/CommunityStructure_Table1_BySite.csv\n\n")

# ==============================================================================
# TABLE 2: BY DEVICE TYPE
# ==============================================================================

cat("--- TABLE 2: BY DEVICE TYPE ---\n\n")

table2_ammod <- tibble::tibble(
  Method = "AMMOD",
  `Mean BL +/- SD (mm)` = sprintf("%.3f +/- %.3f",
                                  mean(site_bl_ammod$mean_BL), mean(site_bl_ammod$sd_BL)),
  `Intra-Habitat Tax. DeltaS (Maize)` = sprintf("%.3f +/- %.3f",
                                                intra_ds_vals %>% dplyr::filter(Pair == "AMMOD1 vs AMMOD2") %>% dplyr::pull(mean_intra_DS),
                                                intra_ds_vals %>% dplyr::filter(Pair == "AMMOD1 vs AMMOD2") %>% dplyr::pull(sd_intra_DS)),
  `Intra-Habitat Tax. DeltaS (Meadow)` = sprintf("%.3f +/- %.3f",
                                                 intra_ds_vals %>% dplyr::filter(Pair == "AMMOD3 vs AMMOD4") %>% dplyr::pull(mean_intra_DS),
                                                 intra_ds_vals %>% dplyr::filter(Pair == "AMMOD3 vs AMMOD4") %>% dplyr::pull(sd_intra_DS)),
  `Tax. DeltaS vs meadow (Site1)` = sprintf("%.3f +/- %.3f",
                                            cross_meadow_ammod %>% dplyr::filter(Site == "Site1") %>% dplyr::pull(mean_cross_DS),
                                            cross_meadow_ammod %>% dplyr::filter(Site == "Site1") %>% dplyr::pull(sd_cross_DS)),
  `Tax. DeltaS vs meadow (Site2)` = sprintf("%.3f +/- %.3f",
                                            cross_meadow_ammod %>% dplyr::filter(Site == "Site2") %>% dplyr::pull(mean_cross_DS),
                                            cross_meadow_ammod %>% dplyr::filter(Site == "Site2") %>% dplyr::pull(sd_cross_DS))
)

table2_faird <- tibble::tibble(
  Method = "FAIRD",
  `Mean BL +/- SD (mm)` = sprintf("%.3f +/- %.3f",
                                  mean(site_bl_faird$mean_BL), mean(site_bl_faird$sd_BL)),
  `Intra-Habitat Tax. DeltaS (Maize)` = sprintf("%.3f +/- %.3f",
                                                intra_ds_vals %>% dplyr::filter(Pair == "FAIRD1 vs FAIRD2") %>% dplyr::pull(mean_intra_DS),
                                                intra_ds_vals %>% dplyr::filter(Pair == "FAIRD1 vs FAIRD2") %>% dplyr::pull(sd_intra_DS)),
  `Intra-Habitat Tax. DeltaS (Meadow)` = sprintf("%.3f +/- %.3f",
                                                 intra_ds_vals %>% dplyr::filter(Pair == "FAIRD3 vs FAIRD4") %>% dplyr::pull(mean_intra_DS),
                                                 intra_ds_vals %>% dplyr::filter(Pair == "FAIRD3 vs FAIRD4") %>% dplyr::pull(sd_intra_DS)),
  `Tax. DeltaS vs meadow (Site1)` = sprintf("%.3f +/- %.3f",
                                            cross_meadow_faird %>% dplyr::filter(Site == "Site1") %>% dplyr::pull(mean_cross_DS),
                                            cross_meadow_faird %>% dplyr::filter(Site == "Site1") %>% dplyr::pull(sd_cross_DS)),
  `Tax. DeltaS vs meadow (Site2)` = sprintf("%.3f +/- %.3f",
                                            cross_meadow_faird %>% dplyr::filter(Site == "Site2") %>% dplyr::pull(mean_cross_DS),
                                            cross_meadow_faird %>% dplyr::filter(Site == "Site2") %>% dplyr::pull(sd_cross_DS))
)

table2 <- dplyr::bind_rows(table2_ammod, table2_faird)
print(table2, n = Inf)
cat("\n")

readr::write_csv(table2, "outputs/CommunityStructure_Table2_ByDevice.csv")
cat("[OK] outputs/CommunityStructure_Table2_ByDevice.csv\n\n")

# ==============================================================================
# TABLE 3: BY HABITAT
# ==============================================================================

cat("--- TABLE 3: BY HABITAT ---\n\n")

# Body size by habitat — per method
bl_by_habitat <- dplyr::bind_rows(
  daily_bl_ammod %>% dplyr::group_by(Habitat, Method) %>%
    dplyr::summarise(mean_BL = mean(mean_bl, na.rm = TRUE),
                     sd_BL   = sd(mean_bl,   na.rm = TRUE), .groups = "drop"),
  daily_bl_faird %>% dplyr::group_by(Habitat, Method) %>%
    dplyr::summarise(mean_BL = mean(mean_bl, na.rm = TRUE),
                     sd_BL   = sd(mean_bl,   na.rm = TRUE), .groups = "drop")
) %>%
  dplyr::mutate(BL_report = sprintf("%.3f +/- %.3f", mean_BL, sd_BL))

# Intra-habitat DeltaS by habitat — aggregate across both pairs per habitat
intra_ds_by_habitat <- all_results %>%
  dplyr::filter(Type == "Intra-habitat") %>%
  dplyr::group_by(Method, Habitat) %>%
  dplyr::summarise(
    mean_DS = mean(Delta_S, na.rm = TRUE),
    sd_DS   = sd(Delta_S,   na.rm = TRUE),
    n       = dplyr::n(),
    .groups = "drop"
  )

table3 <- tidyr::expand_grid(
  Habitat = c("Maize", "Meadow"),
  Method  = c("AMMOD", "FAIRD")
) %>%
  dplyr::left_join(
    bl_by_habitat %>% dplyr::select(Habitat, Method, BL_report),
    by = c("Habitat", "Method")
  ) %>%
  dplyr::left_join(
    intra_ds_by_habitat %>% dplyr::mutate(
      DS_report = sprintf("%.3f +/- %.3f (n=%d)", mean_DS, sd_DS, n)
    ) %>% dplyr::select(Habitat, Method, DS_report),
    by = c("Habitat", "Method")
  ) %>%
  dplyr::rename(
    `Mean BL +/- SD (mm)`       = BL_report,
    `Intra-Habitat Tax. DeltaS` = DS_report
  )

print(table3, n = Inf)
cat("\n")

readr::write_csv(table3, "outputs/CommunityStructure_Table3_ByHabitat.csv")
cat("[OK] outputs/CommunityStructure_Table3_ByHabitat.csv\n\n")

cat("=== SECTION 7 COMPLETE (internal tables — both methods) ===\n\n")

# ==============================================================================
# SECTION 7B: PUBLICATION TABLES — AMMOD ONLY
# ==============================================================================
# Rationale: AMMOD and FAIRD body length metrics are NOT comparable.
#   AMMOD BL = mean of unique families detected per day (presence/absence pool,
#              Söhlström et al. 2018) → reflects taxonomic size diversity
#   FAIRD BL = abundance-weighted mean of individual detections → dominated
#              by small abundant flies; systematically lower due to small-insect
#              abundance bias (not a community size descriptor)
# Publication tables therefore report AMMOD BL only.
# FAIRD BL retained in internal tables (Section 7) and PCV.
# ==============================================================================

cat("=== SECTION 7B: PUBLICATION TABLES (AMMOD only) ===\n\n")

# --- Publication Table 1: By Site (AMMOD only) ---
cat("--- PUB TABLE 1: BY SITE (AMMOD only) ---\n\n")

pub_table1 <- table1 %>%
  dplyr::filter(Method == "AMMOD") %>%
  dplyr::select(-Method)

print(pub_table1, n = Inf)
cat("\nNote: Intra-Habitat Tax. DeltaS is a pair-level metric shared by both sites\n")
cat("within the same habitat (Site1+Site2 form the intra-maize pair; Site3+Site4\n")
cat("form the intra-meadow pair).\n")
cat("BL: family-deduplicated daily mean (Sohlstrom et al. 2018).\n")
cat("Delta BL vs meadow = site mean BL minus grand mean BL of Sites 3+4.\n\n")

readr::write_csv(pub_table1, "outputs/PUB_CommunityStructure_Table1_BySite_AMMOD.csv")
cat("[OK] outputs/PUB_CommunityStructure_Table1_BySite_AMMOD.csv\n\n")

# --- Publication Table 2: By Device (AMMOD only) ---
cat("--- PUB TABLE 2: BY DEVICE TYPE (AMMOD only) ---\n\n")

pub_table2 <- table2 %>%
  dplyr::filter(Method == "AMMOD") %>%
  dplyr::select(-Method)

print(pub_table2, n = Inf)
cat("\n")

readr::write_csv(pub_table2, "outputs/PUB_CommunityStructure_Table2_ByDevice_AMMOD.csv")
cat("[OK] outputs/PUB_CommunityStructure_Table2_ByDevice_AMMOD.csv\n\n")

# --- Publication Table 3: By Habitat (AMMOD only) ---
cat("--- PUB TABLE 3: BY HABITAT (AMMOD only) ---\n\n")

pub_table3 <- table3 %>%
  dplyr::filter(Method == "AMMOD") %>%
  dplyr::select(-Method)

print(pub_table3, n = Inf)
cat("\n")

readr::write_csv(pub_table3, "outputs/PUB_CommunityStructure_Table3_ByHabitat_AMMOD.csv")
cat("[OK] outputs/PUB_CommunityStructure_Table3_ByHabitat_AMMOD.csv\n\n")

cat("=== SECTION 7B COMPLETE ===\n\n")

# ==============================================================================
# SECTION 7B2: PUBLICATION TABLE 14 — BOTH METHODS, ΔS ONLY
# ==============================================================================
# Table 14 in manuscript (§3.3.6.2): combines AMMOD and FAIRD for taxonomic
# similarity columns. BL excluded for FAIRD (not comparable to AMMOD BL;
# see Section 7B rationale). Both methods included for ΔS metrics.
# ==============================================================================

cat("=== SECTION 7B2: PUB TABLE 14 (Both methods, DeltaS only) ===\n\n")

# AMMOD rows: BL + intra-habitat DS + cross-meadow DS
pub_table14_ammod <- table1 %>%
  dplyr::filter(Method == "AMMOD") %>%
  dplyr::select(Site, Habitat, Method,
                `Mean BL +/- SD (mm)`,
                `Intra-Habitat Tax. DeltaS`,
                `Tax. DeltaS vs meadow`)

# FAIRD rows: BL excluded (not comparable), intra-habitat DS + cross-meadow DS
pub_table14_faird <- table1 %>%
  dplyr::filter(Method == "FAIRD") %>%
  dplyr::select(Site, Habitat, Method,
                `Intra-Habitat Tax. DeltaS`,
                `Tax. DeltaS vs meadow`) %>%
  dplyr::mutate(`Mean BL +/- SD (mm)` = "—") %>%
  dplyr::select(Site, Habitat, Method,
                `Mean BL +/- SD (mm)`,
                `Intra-Habitat Tax. DeltaS`,
                `Tax. DeltaS vs meadow`)

pub_table14 <- dplyr::bind_rows(pub_table14_ammod, pub_table14_faird) %>%
  dplyr::arrange(Site, Method)

cat("--- PUB TABLE 14: COMBINED (both methods) ---\n\n")
print(pub_table14, n = Inf)
cat("\nNotes:\n")
cat("  BL reported for AMMOD only (family-deduplicated; Soehlstroem et al. 2018).\n")
cat("  FAIRD BL (abundance-weighted individual detections) is not comparable\n")
cat("  to AMMOD BL and is therefore excluded from this table.\n")
cat("  Intra-Habitat Tax. DeltaS: pair-level metric shared by both sites\n")
cat("  within the same habitat (Site1+Site2 = intra-maize pair; Site3+Site4 = intra-meadow pair).\n")
cat("  Tax. DeltaS vs meadow: mean daily DeltaS of each maize site against\n")
cat("  both meadow sites combined (site-level cross pairs).\n\n")

# Habitat-level summary (both methods) for manuscript Table 14 footer
cat("--- Habitat-level summary (for Table 14 footer) ---\n\n")

pub_table14_habitat <- dplyr::bind_rows(
  # AMMOD habitat level
  table3 %>%
    dplyr::filter(Method == "AMMOD") %>%
    dplyr::rename(`Intra-Habitat Tax. DeltaS` = `Intra-Habitat Tax. DeltaS`) %>%
    dplyr::mutate(`Tax. DeltaS vs meadow (Site1)` =
                    sprintf("%.3f +/- %.3f",
                            cross_meadow_ammod %>%
                              dplyr::filter(Site == "Site1") %>%
                              dplyr::pull(mean_cross_DS),
                            cross_meadow_ammod %>%
                              dplyr::filter(Site == "Site1") %>%
                              dplyr::pull(sd_cross_DS)),
                  `Tax. DeltaS vs meadow (Site2)` =
                    sprintf("%.3f +/- %.3f",
                            cross_meadow_ammod %>%
                              dplyr::filter(Site == "Site2") %>%
                              dplyr::pull(mean_cross_DS),
                            cross_meadow_ammod %>%
                              dplyr::filter(Site == "Site2") %>%
                              dplyr::pull(sd_cross_DS))),
  # FAIRD habitat level
  table3 %>%
    dplyr::filter(Method == "FAIRD") %>%
    dplyr::mutate(`Tax. DeltaS vs meadow (Site1)` =
                    sprintf("%.3f +/- %.3f",
                            cross_meadow_faird %>%
                              dplyr::filter(Site == "Site1") %>%
                              dplyr::pull(mean_cross_DS),
                            cross_meadow_faird %>%
                              dplyr::filter(Site == "Site1") %>%
                              dplyr::pull(sd_cross_DS)),
                  `Tax. DeltaS vs meadow (Site2)` =
                    sprintf("%.3f +/- %.3f",
                            cross_meadow_faird %>%
                              dplyr::filter(Site == "Site2") %>%
                              dplyr::pull(mean_cross_DS),
                            cross_meadow_faird %>%
                              dplyr::filter(Site == "Site2") %>%
                              dplyr::pull(sd_cross_DS)))
) %>%
  dplyr::arrange(Method, Habitat)

print(pub_table14_habitat, n = Inf)
cat("\n")

readr::write_csv(pub_table14,         "outputs/PUB_Table14_CommunityStructure_BothMethods.csv")
readr::write_csv(pub_table14_habitat, "outputs/PUB_Table14_CommunityStructure_HabitatLevel.csv")
cat("[OK] outputs/PUB_Table14_CommunityStructure_BothMethods.csv\n")
cat("[OK] outputs/PUB_Table14_CommunityStructure_HabitatLevel.csv\n\n")

# Statistical tests summary for manuscript footnote
cat("--- Statistical tests (for Table 14 footnote) ---\n\n")

for (method in c("FAIRD", "AMMOD")) {
  dev1 <- if (method == "FAIRD") "FAIRD1" else "AMMOD1"
  dev2 <- if (method == "FAIRD") "FAIRD2" else "AMMOD2"
  dev3 <- if (method == "FAIRD") "FAIRD3" else "AMMOD3"
  dev4 <- if (method == "FAIRD") "FAIRD4" else "AMMOD4"
  
  s1_meadow <- all_results %>%
    dplyr::filter(Pair %in% c(sprintf("%s vs %s", dev1, dev3),
                              sprintf("%s vs %s", dev1, dev4))) %>%
    dplyr::pull(Delta_S)
  
  s2_meadow <- all_results %>%
    dplyr::filter(Pair %in% c(sprintf("%s vs %s", dev2, dev3),
                              sprintf("%s vs %s", dev2, dev4))) %>%
    dplyr::pull(Delta_S)
  
  intra_vals <- all_results %>%
    dplyr::filter(Method == method, Type == "Intra-habitat") %>%
    dplyr::pull(Delta_S)
  
  cross_vals <- all_results %>%
    dplyr::filter(Method == method, Type == "Site-level cross") %>%
    dplyr::pull(Delta_S)
  
  wt_s1s2  <- stats::wilcox.test(s1_meadow, s2_meadow, exact = FALSE)
  cd_s1s2  <- effsize::cohen.d(s1_meadow, s2_meadow)
  wt_ic    <- stats::wilcox.test(intra_vals, cross_vals, exact = FALSE)
  cd_ic    <- effsize::cohen.d(intra_vals, cross_vals)
  
  cat(sprintf("%s — Site1 vs Site2 similarity to meadow:\n", method))
  cat(sprintf("  Site1 vs meadow: n=%d | mean=%.3f +/- %.3f\n",
              length(s1_meadow), mean(s1_meadow), sd(s1_meadow)))
  cat(sprintf("  Site2 vs meadow: n=%d | mean=%.3f +/- %.3f\n",
              length(s2_meadow), mean(s2_meadow), sd(s2_meadow)))
  cat(sprintf("  Wilcoxon W = %.1f | p = %.4f%s | Cohen's d = %.3f (%s)\n\n",
              wt_s1s2$statistic, wt_s1s2$p.value,
              ifelse(wt_s1s2$p.value < 0.05, " *", " ns"),
              cd_s1s2$estimate, cd_s1s2$magnitude))
  
  cat(sprintf("%s — Intra-habitat vs cross-habitat similarity:\n", method))
  cat(sprintf("  Intra-habitat: n=%d | mean=%.3f +/- %.3f\n",
              length(intra_vals), mean(intra_vals), sd(intra_vals)))
  cat(sprintf("  Cross-habitat: n=%d | mean=%.3f +/- %.3f\n",
              length(cross_vals), mean(cross_vals), sd(cross_vals)))
  cat(sprintf("  Wilcoxon W = %.1f | p = %.4f%s | Cohen's d = %.3f (%s)\n\n",
              wt_ic$statistic, wt_ic$p.value,
              ifelse(wt_ic$p.value < 0.05, " *", " ns"),
              cd_ic$estimate, cd_ic$magnitude))
}

cat("=== SECTION 7B2 COMPLETE ===\n\n")

# ==============================================================================
# SECTION 7C: BODY SIZE — PAIRWISE MANN-WHITNEY AND PLOTS (AMMOD only)
# ==============================================================================
# Formal test of whether body size differs between any pair of sites.
# Equivalent to H1_H2_size_by_sites.R Section 2-4 (now consolidated here).
# AMMOD only: family-deduplicated daily mean BL (Söhlström et al. 2018).
# ==============================================================================

cat("=== SECTION 7C: BODY SIZE PAIRWISE TESTS AND PLOTS (AMMOD) ===\n\n")
cat("Pairwise Mann-Whitney tests across all 6 site pairs (common days only)\n")
cat("Bonferroni correction for 6 comparisons (alpha = 0.05/6 = 0.0083)\n\n")

site_pairs_bl <- list(
  c("Site1", "Site2"),  # intra-maize
  c("Site3", "Site4"),  # intra-meadow
  c("Site1", "Site3"),  # cross-habitat
  c("Site1", "Site4"),  # cross-habitat
  c("Site2", "Site3"),  # cross-habitat
  c("Site2", "Site4")   # cross-habitat
)

bonferroni_alpha <- 0.05 / length(site_pairs_bl)

mw_bl_results <- purrr::map_dfr(site_pairs_bl, function(pair) {
  s1 <- pair[1]; s2 <- pair[2]
  
  days_s1 <- daily_bl_ammod %>% dplyr::filter(Site == s1) %>% dplyr::pull(Date)
  days_s2 <- daily_bl_ammod %>% dplyr::filter(Site == s2) %>% dplyr::pull(Date)
  common_days <- intersect(days_s1, days_s2)
  
  vals_s1 <- daily_bl_ammod %>%
    dplyr::filter(Site == s1, Date %in% common_days) %>%
    dplyr::pull(mean_bl)
  vals_s2 <- daily_bl_ammod %>%
    dplyr::filter(Site == s2, Date %in% common_days) %>%
    dplyr::pull(mean_bl)
  
  if (length(vals_s1) < 3 || length(vals_s2) < 3) {
    cat(sprintf("  %s vs %s: Insufficient common days (n=%d) — skipped\n",
                s1, s2, length(common_days)))
    return(NULL)
  }
  
  wt <- stats::wilcox.test(vals_s1, vals_s2, exact = FALSE)
  hab_s1 <- dplyr::if_else(s1 %in% c("Site1","Site2"), "Maize", "Meadow")
  hab_s2 <- dplyr::if_else(s2 %in% c("Site1","Site2"), "Maize", "Meadow")
  
  tibble::tibble(
    Pair             = paste(s1, "vs", s2),
    Type             = dplyr::if_else(hab_s1 == hab_s2, "Intra-habitat", "Cross-habitat"),
    n_common_days    = length(common_days),
    median_s1        = round(median(vals_s1), 3),
    median_s2        = round(median(vals_s2), 3),
    W_statistic      = wt$statistic,
    p_value          = round(wt$p.value, 4),
    p_bonferroni     = round(min(wt$p.value * length(site_pairs_bl), 1), 4),
    significant_bonf = wt$p.value < bonferroni_alpha
  )
})

cat("Results:\n")
print(mw_bl_results, n = Inf)
cat(sprintf("\nBonferroni-adjusted alpha: %.4f\n\n", bonferroni_alpha))

readr::write_csv(mw_bl_results, "outputs/BodySize_PairwiseMannWhitney_AMMOD.csv")
cat("[OK] outputs/BodySize_PairwiseMannWhitney_AMMOD.csv\n\n")

# --- Time series plot: daily mean BL by site (AMMOD) ---
site_colors_bl <- c(
  "Site1" = unname(custom_colors["Site1"]),
  "Site2" = unname(custom_colors["Site2"]),
  "Site3" = unname(custom_colors["Site3"]),
  "Site4" = unname(custom_colors["Site4"])
)

plot_bl_ts <- ggplot(
  daily_bl_ammod,
  aes(x = Date, y = mean_bl, color = Site, shape = Site, linetype = Habitat)
) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = site_colors_bl) +
  scale_shape_manual(values = c(
    "Site1" = 16, "Site2" = 17, "Site3" = 15, "Site4" = 18
  )) +
  scale_linetype_manual(values = c("Maize" = "solid", "Meadow" = "dashed")) +
  labs(
    title    = "Daily Mean Body Length of Detected Taxa by Site (AMMOD)",
    subtitle = "Family-level deduplication | Biometrics: Söhlström et al. (2018)",
    x        = "Date",
    y        = "Mean Body Length (mm)",
    color    = "Site", shape = "Site", linetype = "Habitat"
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

# --- Violin plot: BL distribution by site (AMMOD) ---
plot_bl_violin <- ggplot(
  daily_bl_ammod %>%
    dplyr::mutate(Habitat = dplyr::if_else(
      Site %in% c("Site1","Site2"), "Maize", "Meadow")),
  aes(x = Site, y = mean_bl, fill = Habitat)
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = 16,
               outlier.size = 2, alpha = 0.8) +
  scale_fill_manual(values = c(
    "Maize"  = unname(custom_colors["Maize"]),
    "Meadow" = unname(custom_colors["Meadow"])
  )) +
  labs(
    title = "Distribution of Daily Mean Body Length by Site (AMMOD)",
    x = "Site", y = "Mean Body Length (mm)", fill = "Habitat"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

combined_bl <- plot_bl_ts + plot_bl_violin +
  patchwork::plot_layout(ncol = 2, widths = c(1.6, 1)) +
  patchwork::plot_annotation(
    title = "Body Size Distribution Across Sites (AMMOD)",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

print(combined_bl)

ggsave("outputs/BodySize_Combined_AMMOD.png",
       plot = combined_bl, width = 16, height = 5, dpi = 300, bg = "white")
ggsave("outputs/BodySize_Combined_AMMOD.pdf",
       plot = combined_bl, width = 16, height = 5, device = cairo_pdf)

cat("[OK] outputs/BodySize_Combined_AMMOD.png/.pdf\n\n")
cat("=== SECTION 7C COMPLETE ===\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")