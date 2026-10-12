# ==============================================================================
# H1: INTRA-SYSTEM TAXONOMIC CONSISTENCY — FAIRD REPLICATES
# ==============================================================================
# Purpose: Assess whether replicate FAIR-Device units detect the same taxa
#          in the same proportions on the same days.
#
# Metric:  Morisita-Horn similarity (1 - vegdist(method="horn"))
#          Computed at three fixed taxonomic levels in parallel:
#            - SUPERFAMILY (primary, ~68% of individuals)
#            - ORDER       (control,  ~97% of individuals)
#            - GENUS        (fine,    ~43% of individuals; consistent with H3)
#
# Temporal resolution:
#   - DAILY:     original daily resolution (full dataset, no filtering)
#   - 3-DAY:     non-overlapping 3-day windows, filtered to >= 3 individuals
#                per device per window at the corresponding taxonomic level
#
# Pairs:   FAIRD1 <-> FAIRD2 (Maize, Sites 1-2)
#          FAIRD3 <-> FAIRD4 (Meadow, Sites 3-4)
#
# Author: Juan Chiavassa & Claude
# Date: March 2026
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(vegan)       # vegdist(method = "horn") for Morisita-Horn
library(lubridate)
library(effsize)     # cohen.d()

# ==============================================================================
# CONSOLE OUTPUT CAPTURE
# ==============================================================================
script_name <- "H1_consistency_FAIRD_taxonomic"

output_file <- file.path("outputs",
                         paste0(format(Sys.time(), "%Y%m%d_%H%M%S"),
                                "_", script_name, ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

while (sink.number() > 0) sink()

sink(output_file, split = TRUE)
cat("==============================================================================\n")
cat("H1: INTRA-SYSTEM TAXONOMIC CONSISTENCY — FAIRD REPLICATES\n")
cat("Morisita-Horn similarity: Daily + 3-day windows\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: LOAD DATA
# ==============================================================================

cat("=== SECTION 1: LOAD DATA ===\n\n")

source("load_data.R")
source("custom_colors.R")

faird_data <- data_list$faird

cat("\nFAIRD data loaded:", nrow(faird_data), "individuals\n")
cat("Devices:", paste(sort(unique(faird_data$Device)), collapse = ", "), "\n")
cat("Date range:", as.character(min(faird_data$Date)), "to",
    as.character(max(faird_data$Date)), "\n\n")


# ==============================================================================
# SECTION 2: BUILD TAXON IDS AND DEFINE WINDOWS
# ==============================================================================

cat("=== SECTION 2: BUILD TAXON IDS AND DEFINE WINDOWS ===\n\n")

nc <- "#N/C"

# Add taxonomic level columns
faird_taxo <- faird_data %>%
  dplyr::mutate(
    Taxon_Genus = dplyr::if_else(
      !is.na(Genus) & Genus != nc, Genus, NA_character_
    ),
    Taxon_Superfamily = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != nc, Superfamily, NA_character_
    ),
    Taxon_Order = dplyr::if_else(
      !is.na(Order) & Order != nc, Order, NA_character_
    )
  )

# Coverage diagnostics
n_total <- nrow(faird_taxo)
n_genus <- sum(!is.na(faird_taxo$Taxon_Genus))
n_sfam  <- sum(!is.na(faird_taxo$Taxon_Superfamily))
n_order <- sum(!is.na(faird_taxo$Taxon_Order))

cat("Taxonomic level coverage (n =", n_total, "individuals):\n")
cat(sprintf("  Genus:       %d (%.1f%%)\n", n_genus, 100 * n_genus / n_total))
cat(sprintf("  Superfamily: %d (%.1f%%)\n", n_sfam,  100 * n_sfam  / n_total))
cat(sprintf("  Order:       %d (%.1f%%)\n", n_order, 100 * n_order / n_total))

# Per-device breakdown
cat("\nCoverage by device:\n")
device_coverage <- faird_taxo %>%
  dplyr::group_by(Device) %>%
  dplyr::summarise(
    n_total       = dplyr::n(),
    pct_Genus     = round(100 * sum(!is.na(Taxon_Genus))       / dplyr::n(), 1),
    pct_Superfam  = round(100 * sum(!is.na(Taxon_Superfamily)) / dplyr::n(), 1),
    pct_Order     = round(100 * sum(!is.na(Taxon_Order))       / dplyr::n(), 1),
    .groups = "drop"
  )
print(device_coverage)

# --- Define 3-day windows ---
experiment_start <- as.Date("2023-08-23")
experiment_end   <- as.Date("2023-09-13")
all_dates <- seq(experiment_start, experiment_end, by = "day")

# Build window lookup table
window_breaks <- seq(experiment_start, experiment_end + 2, by = 3)
window_lookup <- tibble::tibble(Date = all_dates) %>%
  dplyr::mutate(
    Window_ID = as.integer(cut(Date, breaks = window_breaks, right = FALSE,
                               labels = FALSE)),
    Window_Start = window_breaks[Window_ID],
    Window_End   = pmin(Window_Start + 2, experiment_end),
    Window_Label = paste0("W", Window_ID, " (",
                          format(Window_Start, "%b %d"), "-",
                          format(Window_End, "%b %d"), ")")
  )

# Attach window info to data
faird_taxo <- faird_taxo %>%
  dplyr::left_join(window_lookup, by = "Date")

cat("\n3-day window definitions:\n")
window_summary <- window_lookup %>%
  dplyr::group_by(Window_ID, Window_Label) %>%
  dplyr::summarise(n_days = dplyr::n(), .groups = "drop")
print(window_summary, n = Inf)
cat("\n")


# ==============================================================================
# SECTION 3: MORISITA-HORN COMPUTATION FUNCTION
# ==============================================================================

# Shared function for both daily and windowed computation
compute_mh <- function(data, dev_a, dev_b, taxon_col) {
  # Aggregate counts per taxon for each device
  day_data <- data %>%
    dplyr::filter(Device %in% c(dev_a, dev_b),
                  !is.na(.data[[taxon_col]])) %>%
    dplyr::group_by(Device, .data[[taxon_col]]) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  devs_present <- unique(day_data$Device)
  if (length(devs_present) < 2) return(NA_real_)
  
  comm_matrix <- day_data %>%
    tidyr::pivot_wider(
      names_from  = dplyr::all_of(taxon_col),
      values_from = count,
      values_fill = 0
    ) %>%
    dplyr::arrange(Device) %>%
    dplyr::select(-Device)
  
  row_sums <- rowSums(comm_matrix)
  if (any(row_sums == 0)) return(NA_real_)
  
  mh_dissim <- vegan::vegdist(as.matrix(comm_matrix), method = "horn")
  return(1 - as.numeric(mh_dissim))
}

# Get individual/taxa counts per device from a data subset
get_counts <- function(data, device_name, taxon_col) {
  sub <- data %>%
    dplyr::filter(Device == device_name, !is.na(.data[[taxon_col]]))
  list(
    n_indiv = nrow(sub),
    n_taxa  = length(unique(sub[[taxon_col]]))
  )
}


# ==============================================================================
# SECTION 4: COMPUTE DAILY M-H (FULL, NO FILTER)
# ==============================================================================

cat("=== SECTION 4: COMPUTE DAILY M-H (FULL, NO FILTER) ===\n\n")

pairs <- list(
  list(dev_a = "FAIRD1", dev_b = "FAIRD2", habitat = "Maize"),
  list(dev_a = "FAIRD3", dev_b = "FAIRD4", habitat = "Meadow")
)

tax_levels <- c("Taxon_Order", "Taxon_Superfamily", "Taxon_Genus")
tax_labels <- c("Order", "Superfamily", "Genus")

daily_results <- list()

for (p in seq_along(pairs)) {
  pair <- pairs[[p]]
  cat(sprintf("Processing pair: %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))
  
  for (lev in seq_along(tax_levels)) {
    tax_col   <- tax_levels[lev]
    tax_label <- tax_labels[lev]
    
    for (d in all_dates) {
      date_val <- as.Date(d, origin = "1970-01-01")
      day_data <- faird_taxo %>% dplyr::filter(Date == date_val)
      
      mh_sim   <- compute_mh(day_data, pair$dev_a, pair$dev_b, tax_col)
      counts_a <- get_counts(day_data, pair$dev_a, tax_col)
      counts_b <- get_counts(day_data, pair$dev_b, tax_col)
      
      daily_results[[length(daily_results) + 1]] <- tibble::tibble(
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

daily_df <- dplyr::bind_rows(daily_results)
cat("\nDaily computations:", nrow(daily_df), "\n")
cat("Valid M-H values:", sum(!is.na(daily_df$MH_Sim)), "\n\n")


# ==============================================================================
# SECTION 5: COMPUTE 3-DAY WINDOWED M-H (FILTERED >= 3 INDIVIDUALS)
# ==============================================================================

cat("=== SECTION 5: COMPUTE 3-DAY WINDOWED M-H (FILTERED >= 3 INDIV.) ===\n\n")

MIN_INDIV <- 3
cat(sprintf("Minimum individuals per device per window: %d\n\n", MIN_INDIV))

windowed_results <- list()

for (p in seq_along(pairs)) {
  pair <- pairs[[p]]
  cat(sprintf("Processing pair: %s <-> %s (%s)\n", pair$dev_a, pair$dev_b, pair$habitat))
  
  for (lev in seq_along(tax_levels)) {
    tax_col   <- tax_levels[lev]
    tax_label <- tax_labels[lev]
    
    for (w in unique(na.omit(faird_taxo$Window_ID))) {
      win_data <- faird_taxo %>% dplyr::filter(Window_ID == w)
      win_label <- unique(win_data$Window_Label)[1]
      
      counts_a <- get_counts(win_data, pair$dev_a, tax_col)
      counts_b <- get_counts(win_data, pair$dev_b, tax_col)
      
      # Apply minimum individual filter
      passes_filter <- (counts_a$n_indiv >= MIN_INDIV &
                          counts_b$n_indiv >= MIN_INDIV)
      
      if (passes_filter) {
        mh_sim <- compute_mh(win_data, pair$dev_a, pair$dev_b, tax_col)
      } else {
        mh_sim <- NA_real_
      }
      
      windowed_results[[length(windowed_results) + 1]] <- tibble::tibble(
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
        Filtered   = !passes_filter
      )
    }
  }
}

windowed_df <- dplyr::bind_rows(windowed_results)

cat("\nWindowed computations:", nrow(windowed_df), "\n")
cat("Valid M-H values (passed filter):", sum(!is.na(windowed_df$MH_Sim)), "\n")
cat("Filtered out (< ", MIN_INDIV, " indiv.):", sum(windowed_df$Filtered), "\n\n")

# Detailed filter impact per pair and level
cat("--- Filter impact by pair and level ---\n")
filter_impact <- windowed_df %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_windows     = dplyr::n(),
    n_valid       = sum(!is.na(MH_Sim)),
    n_filtered    = sum(Filtered),
    pct_retained  = round(100 * sum(!is.na(MH_Sim)) / dplyr::n(), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")),
                 Habitat)
print(filter_impact, n = Inf)
cat("\n")


# ==============================================================================
# SECTION 6: COMBINE ALL RESULTS
# ==============================================================================

results_all <- dplyr::bind_rows(daily_df, windowed_df)


# ==============================================================================
# SECTION 7: SUMMARY STATISTICS
# ==============================================================================

cat("=== SECTION 7: SUMMARY STATISTICS ===\n\n")

# --- 7a: Daily (unfiltered) summary ---
cat("--- 7a: Daily (unfiltered) summary by pair and level ---\n\n")

summary_daily <- daily_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_days    = dplyr::n(),
    mean_MH   = round(mean(MH_Sim) * 100, 1),
    sd_MH     = round(sd(MH_Sim) * 100, 1),
    median_MH = round(median(MH_Sim) * 100, 1),
    min_MH    = round(min(MH_Sim) * 100, 1),
    max_MH    = round(max(MH_Sim) * 100, 1),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")),
                 Habitat)
print(summary_daily, n = Inf, width = 200)

# --- 7b: 3-day windowed (filtered) summary ---
cat("\n--- 7b: 3-day windowed (filtered >= 3 indiv.) summary ---\n\n")

summary_windowed <- windowed_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::group_by(Pair, Habitat, Tax_Level) %>%
  dplyr::summarise(
    n_windows = dplyr::n(),
    mean_MH   = round(mean(MH_Sim) * 100, 1),
    sd_MH     = round(sd(MH_Sim) * 100, 1),
    median_MH = round(median(MH_Sim) * 100, 1),
    min_MH    = round(min(MH_Sim) * 100, 1),
    max_MH    = round(max(MH_Sim) * 100, 1),
    mean_indiv_A = round(mean(N_indiv_A), 1),
    mean_indiv_B = round(mean(N_indiv_B), 1),
    mean_taxa_A  = round(mean(N_taxa_A), 1),
    mean_taxa_B  = round(mean(N_taxa_B), 1),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")),
                 Habitat)
print(summary_windowed, n = Inf, width = 200)

readr::write_csv(summary_windowed,
                 "outputs/H1_taxonomic_MH_summary_windowed.csv")
cat("\n[OK] outputs/H1_taxonomic_MH_summary_windowed.csv\n\n")


# --- 7c: Comparison table: daily vs windowed at Superfamily ---
cat("--- 7c: Daily vs Windowed comparison (Superfamily) ---\n\n")

comp_daily <- daily_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily") %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    Resolution = "Daily (unfiltered)",
    n = dplyr::n(),
    mean_MH = round(mean(MH_Sim) * 100, 1),
    sd_MH   = round(sd(MH_Sim) * 100, 1),
    median_MH = round(median(MH_Sim) * 100, 1),
    .groups = "drop"
  )

comp_windowed <- windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily") %>%
  dplyr::group_by(Pair, Habitat) %>%
  dplyr::summarise(
    Resolution = "3-day (filtered >= 3)",
    n = dplyr::n(),
    mean_MH = round(mean(MH_Sim) * 100, 1),
    sd_MH   = round(sd(MH_Sim) * 100, 1),
    median_MH = round(median(MH_Sim) * 100, 1),
    .groups = "drop"
  )

comparison_table <- dplyr::bind_rows(comp_daily, comp_windowed) %>%
  dplyr::arrange(Pair, Resolution)

print(comparison_table, n = Inf, width = 200)

readr::write_csv(comparison_table,
                 "outputs/H1_taxonomic_MH_daily_vs_windowed.csv")
cat("\n[OK] outputs/H1_taxonomic_MH_daily_vs_windowed.csv\n\n")


# --- 7d: Habitat comparison tests (windowed, Superfamily) ---
cat("--- 7d: Habitat comparison (windowed, Superfamily) ---\n\n")

maize_w <- windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily",
                Habitat == "Maize") %>%
  dplyr::pull(MH_Sim)

meadow_w <- windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily",
                Habitat == "Meadow") %>%
  dplyr::pull(MH_Sim)

if (length(maize_w) >= 3 & length(meadow_w) >= 3) {
  group_data <- data.frame(
    val = c(maize_w, meadow_w),
    grp = factor(c(rep("Maize", length(maize_w)),
                   rep("Meadow", length(meadow_w))))
  )
  levene_p <- car::leveneTest(val ~ grp, data = group_data)$`Pr(>F)`[1]
  cat(sprintf("Levene's test p = %.4f\n", levene_p))
  
  if (levene_p < 0.05) {
    tt <- t.test(maize_w, meadow_w, var.equal = FALSE)
    cat("Using Welch's t-test (unequal variances)\n")
  } else {
    tt <- t.test(maize_w, meadow_w, var.equal = TRUE)
    cat("Using Student's t-test (equal variances)\n")
  }
  cat(sprintf("Maize:  %.1f%% +/- %.1f%% (n=%d windows)\n",
              mean(maize_w) * 100, sd(maize_w) * 100, length(maize_w)))
  cat(sprintf("Meadow: %.1f%% +/- %.1f%% (n=%d windows)\n",
              mean(meadow_w) * 100, sd(meadow_w) * 100, length(meadow_w)))
  cat(sprintf("t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
  
  cd <- effsize::cohen.d(maize_w, meadow_w)
  cat(sprintf("Cohen's d = %.3f (%s effect)\n\n", cd$estimate, cd$magnitude))
} else {
  cat("Insufficient data for habitat comparison\n\n")
}


# ==============================================================================
# SECTION 8: TEMPORAL TREND ANALYSIS (WINDOWED, SUPERFAMILY)
# ==============================================================================

cat("=== SECTION 8: TEMPORAL TREND ANALYSIS ===\n\n")

cat("--- Temporal trend: 3-day windows, Superfamily level ---\n\n")

sfam_windowed <- windowed_df %>%
  dplyr::filter(!is.na(MH_Sim), Tax_Level == "Superfamily") %>%
  dplyr::mutate(
    Window_Num = as.numeric(gsub("W(\\d+).*", "\\1", Period))
  )

if (nrow(sfam_windowed) >= 5) {
  trend_model <- lm(MH_Sim ~ Window_Num, data = sfam_windowed)
  trend_summ  <- summary(trend_model)
  
  cat(sprintf("Overall: slope = %.4f per window (%.2f pp/window)\n",
              coef(trend_model)[2], coef(trend_model)[2] * 100))
  cat(sprintf("  R2 = %.3f, p = %.4f, n = %d windows\n\n",
              trend_summ$r.squared, coef(trend_summ)[2, 4], nrow(sfam_windowed)))
  
  for (pair_name in unique(sfam_windowed$Pair)) {
    pair_data <- sfam_windowed %>% dplyr::filter(Pair == pair_name)
    if (nrow(pair_data) >= 4) {
      m <- lm(MH_Sim ~ Window_Num, data = pair_data)
      s <- summary(m)
      cat(sprintf("  %s: slope = %.4f (%.2f pp/win), R2 = %.3f, p = %.4f, n = %d\n",
                  pair_name, coef(m)[2], coef(m)[2] * 100,
                  s$r.squared, coef(s)[2, 4], nrow(pair_data)))
    }
  }
}

cat("\n")

# ==============================================================================
# SECTION 8.5: INTER-HABITAT TAXONOMIC CONSISTENCY (MAIZE VS MEADOW)
# ==============================================================================

cat("=== SECTION 8.5: INTER-HABITAT TAXONOMIC CONSISTENCY (MAIZE VS MEADOW) ===\n\n")

# 1. Create the pooled dataset
pooled_taxo <- faird_taxo %>%
  dplyr::mutate(
    Pool = dplyr::case_when(
      Device %in% c("FAIRD1", "FAIRD2") ~ "Maize_Pool",
      Device %in% c("FAIRD3", "FAIRD4") ~ "Meadow_Pool",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Pool)) %>%
  # Rename Pool to Device to re-use existing functions
  dplyr::select(-Device) %>%
  dplyr::rename(Device = Pool) 

inter_habitat_results <- list()
tax_col <- "Taxon_Superfamily"
tax_label <- "Superfamily"

cat("Processing Inter-Habitat: Maize_Pool <-> Meadow_Pool (Superfamily Level)\n")

# 2. Calculate Morisita-Horn in 3-day windows
for (w in unique(na.omit(pooled_taxo$Window_ID))) {
  win_data <- pooled_taxo %>% dplyr::filter(Window_ID == w)
  win_label <- unique(win_data$Window_Label)[1]
  
  counts_a <- get_counts(win_data, "Maize_Pool", tax_col)
  counts_b <- get_counts(win_data, "Meadow_Pool", tax_col)
  
  # Apply filter: minimum 3 individuals per HABITAT per window
  passes_filter <- (counts_a$n_indiv >= MIN_INDIV & counts_b$n_indiv >= MIN_INDIV)
  
  if (passes_filter) {
    mh_sim <- compute_mh(win_data, "Maize_Pool", "Meadow_Pool", tax_col)
  } else {
    mh_sim <- NA_real_
  }
  
  inter_habitat_results[[length(inter_habitat_results) + 1]] <- tibble::tibble(
    Resolution = "3-day window",
    Pair       = "Maize vs Meadow",
    Habitat    = "Inter-Habitat",
    Period     = win_label,
    Tax_Level  = tax_label,
    MH_Sim     = mh_sim,
    N_indiv_A  = counts_a$n_indiv, # Maize
    N_indiv_B  = counts_b$n_indiv, # Meadow
    N_taxa_A   = counts_a$n_taxa,
    N_taxa_B   = counts_b$n_taxa,
    Filtered   = !passes_filter
  )
}

inter_habitat_df <- dplyr::bind_rows(inter_habitat_results)

cat("\nInter-Habitat windowed computations:", nrow(inter_habitat_df), "\n")
cat("Valid M-H values (passed filter):", sum(!is.na(inter_habitat_df$MH_Sim)), "\n")
cat("Filtered out (< ", MIN_INDIV, " indiv.):", sum(inter_habitat_df$Filtered), "\n\n")

# 3. Print the summary
cat("--- Inter-Habitat Summary (Filtered >= 3 indiv.) ---\n\n")
summary_inter <- inter_habitat_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::summarise(
    Pair = "Maize vs Meadow",
    n_windows = dplyr::n(),
    mean_MH   = round(mean(MH_Sim) * 100, 1),
    sd_MH     = round(sd(MH_Sim) * 100, 1),
    median_MH = round(median(MH_Sim) * 100, 1),
    min_MH    = round(min(MH_Sim) * 100, 1),
    max_MH    = round(max(MH_Sim) * 100, 1),
    mean_indiv_Maize = round(mean(N_indiv_A), 1),
    mean_indiv_Meadow = round(mean(N_indiv_B), 1),
    .groups = "drop"
  )

print(summary_inter, width = 200)

# 4. Export specific results
readr::write_csv(inter_habitat_df, "outputs/H1_taxonomic_MH_inter_habitat_results.csv")
cat("\n[OK] outputs/H1_taxonomic_MH_inter_habitat_results.csv\n\n")

# ==============================================================================
# SECTION 9: FIGURES
# ==============================================================================

cat("=== SECTION 9: FIGURES ===\n\n")

# --- 9a: Windowed Superfamily temporal plot (PRIMARY FIGURE) ---
cat("Creating windowed Superfamily temporal plot...\n")

plot_w_sfam <- windowed_df %>%
  dplyr::filter(Tax_Level == "Superfamily", !is.na(MH_Sim)) %>%
  dplyr::mutate(
    MH_pct = MH_Sim * 100,
    Window_Num = as.numeric(gsub("W(\\d+).*", "\\1", Period))
  )

p_windowed_sfam <- ggplot(plot_w_sfam,
                          aes(x = Window_Num, y = MH_pct, color = Habitat)) +
  geom_line(aes(group = Pair), linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50",
             linewidth = 0.5) +
  scale_x_continuous(
    breaks = unique(plot_w_sfam$Window_Num),
    labels = function(x) {
      labs <- unique(plot_w_sfam[, c("Window_Num", "Period")])
      labs <- labs[order(labs$Window_Num), ]
      gsub("W\\d+ \\((.*)\\)", "\\1", labs$Period[match(x, labs$Window_Num)])
    }
  ) +
  scale_color_manual(values = custom_colors,
                     labels = c("Maize"  = "FAIRD1 vs FAIRD2 (Maize)",
                                "Meadow" = "FAIRD3 vs FAIRD4 (Meadow)")) +
  labs(
    title    = "H1 Qualitative - Taxonomic Consistency: FAIRD Replicates — Intra-Habitat Level",
    subtitle = paste0("Morisita-Horn index — 3-Day Windows — Superfamily level filtered: >= ",
                      MIN_INDIV, " individuals per device per window"),
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

print(p_windowed_sfam)

ggsave("outputs/H1_taxonomic_MH_windowed_Superfamily.png",
       plot = p_windowed_sfam, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/H1_taxonomic_MH_windowed_Superfamily.png\n\n")


# --- 9b: Inter-Habitat Level — Windowed Superfamily temporal ---
cat("Creating Inter-Habitat Level temporal plot...\n")

plot_inter <- inter_habitat_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::mutate(
    MH_pct = MH_Sim * 100,
    Window_Num = as.numeric(gsub("W(\\d+).*", "\\1", Period))
  )

inter_mean <- mean(plot_inter$MH_pct)

p_inter_habitat <- ggplot(plot_inter,
                          aes(x = Window_Num, y = MH_pct)) +
  geom_line(linewidth = 0.8, color = custom_colors["FAIRD"]) +
  geom_point(size = 3, color = custom_colors["FAIRD"]) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50",
             linewidth = 0.5) +
  geom_hline(yintercept = inter_mean,
             linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_x_continuous(
    breaks = unique(plot_inter$Window_Num),
    labels = function(x) {
      labs <- unique(plot_inter[, c("Window_Num", "Period")])
      labs <- labs[order(labs$Window_Num), ]
      gsub("W\\d+ \\((.*)\\)", "\\1", labs$Period[match(x, labs$Window_Num)])
    }
  ) +
  labs(
    title    = "H1 Qualitative - Taxonomic Consistency: FAIRD Replicates — Inter-Habitat Level",
    subtitle = sprintf("Morisita-Horn index — 3-Day Windows — Superfamily level | Maize vs Meadow: %.1f%% +/- %.1f%% (n=%d)",
                       mean(plot_inter$MH_pct), sd(plot_inter$MH_pct),
                       nrow(plot_inter)),
    x        = "3-Day Window",
    y        = "Morisita-Horn Similarity (%)",
    caption  = "Dashed line = 50% reference. Dotted line = mean."
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

print(p_inter_habitat)

ggsave("outputs/H1_taxonomic_MH_windowed_inter_habitat.png",
       plot = p_inter_habitat, width = 10, height = 6, dpi = 300)
cat("[OK] outputs/H1_taxonomic_MH_windowed_inter_habitat.png\n\n")


# --- 9c: Multi-level boxplot comparison (windowed) ---
cat("Creating multi-level boxplot (windowed)...\n")

plot_w_levels <- windowed_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::mutate(
    MH_pct = MH_Sim * 100,
    Tax_Level = factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")),
    Pair_label = dplyr::case_when(
      Pair == "FAIRD1 vs FAIRD2" ~ "FAIRD1 vs FAIRD2\n(Maize)",
      Pair == "FAIRD3 vs FAIRD4" ~ "FAIRD3 vs FAIRD4\n(Meadow)",
      TRUE ~ Pair
    )
  )

p_w_levels <- ggplot(plot_w_levels,
                     aes(x = Tax_Level, y = MH_pct, fill = Habitat)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  facet_wrap(~Pair_label) +
  scale_fill_manual(values = custom_colors) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  labs(
    title    = "H1 Qualitative - Taxonomic Consistency: FAIRD Replicates — by Taxonomic Level",
    subtitle = paste0("Morisita-Horn index — 3-Day Windows — filtered: >= ",
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

print(p_w_levels)

ggsave("outputs/H1_taxonomic_MH_windowed_by_level_boxplot.png",
       plot = p_w_levels, width = 9, height = 6, dpi = 300)
cat("[OK] outputs/H1_taxonomic_MH_windowed_by_level_boxplot.png\n\n")


# --- 9d: Windowed all levels temporal (faceted by pair) ---
cat("Creating multi-level temporal plot (windowed)...\n")

plot_w_all <- windowed_df %>%
  dplyr::filter(!is.na(MH_Sim)) %>%
  dplyr::mutate(
    MH_pct = MH_Sim * 100,
    Window_Num = as.numeric(gsub("W(\\d+).*", "\\1", Period)),
    Tax_Level = factor(Tax_Level, levels = c("Order", "Superfamily", "Genus")),
    Pair_label = dplyr::case_when(
      Pair == "FAIRD1 vs FAIRD2" ~ "FAIRD1 vs FAIRD2 (Maize)",
      Pair == "FAIRD3 vs FAIRD4" ~ "FAIRD3 vs FAIRD4 (Meadow)",
      TRUE ~ Pair
    )
  )

level_colors <- c(
  "Order"       = "steelblue",
  "Superfamily" = "darkorange",
  "Genus"       = "forestgreen"
)

p_w_temporal_all <- ggplot(plot_w_all,
                           aes(x = Window_Num, y = MH_pct, color = Tax_Level)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50",
             linewidth = 0.4) +
  facet_wrap(~Pair_label, ncol = 1) +
  scale_x_continuous(
    breaks = sort(unique(plot_w_all$Window_Num)),
    labels = function(x) paste0("W", x)
  ) +
  scale_color_manual(values = level_colors) +
  labs(
    title    = "H1 Taxonomic Consistency: Temporal M-H by Level (3-Day Windows)",
    subtitle = paste0("Filtered: >= ", MIN_INDIV, " individuals per device per window"),
    x        = "Window",
    y        = "Morisita-Horn Similarity (%)",
    color    = "Taxonomic Level",
    caption  = "Dashed line = 50%. Missing points = filtered windows."
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

print(p_w_temporal_all)

ggsave("outputs/H1_taxonomic_MH_windowed_temporal_all_levels.png",
       plot = p_w_temporal_all, width = 10, height = 8, dpi = 300)
cat("[OK] outputs/H1_taxonomic_MH_windowed_temporal_all_levels.png\n\n")


# ==============================================================================
# SECTION 10: EXPORT FULL RESULTS
# ==============================================================================

cat("=== SECTION 10: EXPORT FULL RESULTS ===\n\n")

readr::write_csv(daily_df, "outputs/H1_taxonomic_MH_daily_results.csv")
cat("[OK] outputs/H1_taxonomic_MH_daily_results.csv\n")

readr::write_csv(windowed_df, "outputs/H1_taxonomic_MH_windowed_results.csv")
cat("[OK] outputs/H1_taxonomic_MH_windowed_results.csv\n")

cat("\nTotal rows: daily =", nrow(daily_df),
    "| windowed =", nrow(windowed_df), "\n\n")


# ==============================================================================
# WRAP UP
# ==============================================================================

cat("==============================================================================\n")
cat("Script completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

while (sink.number() > 0) sink()