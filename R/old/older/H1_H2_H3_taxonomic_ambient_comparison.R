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

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")