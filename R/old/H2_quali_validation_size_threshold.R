# ==============================================================================
# H2 QUALITATIVE: SIZE THRESHOLD SENSITIVITY (raw vs 3mm vs 5mm) - SUPERFAMILY & 3-DAY WINDOWS
# ==============================================================================
# Three AMMOD pools are compared against the same FAIRD Superfamily pool:
#   raw  : all AMMOD OTUs (no body-length filter)
#   3mm  : AMMOD OTUs with Body_length >= 3 mm
#   5mm  : AMMOD OTUs with Body_length >= 5 mm
# The raw vs 3mm delta quantifies how much concordance changes when the AMMOD
# pool is restricted to the FAIRD-detectable size fraction.
# ==============================================================================
library(dplyr)
library(tidyr)
library(vegan)

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


source("load_data.R")

# --- 1. CONFIGURATION & TEMPORAL SCOPE ---
SIZE_3MM <- 3
SIZE_5MM <- 5
date_start <- as.Date("2023-08-23")
date_end   <- as.Date("2023-09-13")
date_exclude <- as.Date("2023-09-04")
ammod4_failure_start <- as.Date("2023-08-28")
ammod4_failure_end   <- as.Date("2023-09-03")
nc <- "#N/C"
MIN_INDIV <- 3

cat("=== SIZE THRESHOLD SENSITIVITY (SUPERFAMILY, 3-DAY WINDOWS) ===\n\n")

# --- 2. PREPARE FAIRD DATA ---
faird_data <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD", Date >= date_start & Date <= date_end) %>%
  dplyr::mutate(Taxon_SF = dplyr::if_else(!is.na(Superfamily) & Superfamily != nc & Superfamily != "", Superfamily, NA_character_))

# --- 3. PREPARE AMMOD DATA (raw, 3mm, and 5mm) ---
ammod_base <- data_list$ammod_taxonomy %>%
  dplyr::filter(Date >= date_start & Date <= date_end, Date != date_exclude,
                !(Device == "AMMOD4" & Date >= ammod4_failure_start & Date <= ammod4_failure_end)) %>%
  dplyr::mutate(Taxon_SF = dplyr::if_else(!is.na(Superfamily) & Superfamily != nc & Superfamily != "", Superfamily, NA_character_))

# raw: all OTUs regardless of body length (includes rows with NA Body_length)
ammod_raw <- ammod_base
# 3mm and 5mm: require non-missing Body_length and minimum size
ammod_3mm <- ammod_base %>% dplyr::filter(!is.na(Body_length) & Body_length >= SIZE_3MM)
ammod_5mm <- ammod_base %>% dplyr::filter(!is.na(Body_length) & Body_length >= SIZE_5MM)

cat("AMMOD raw OTUs:     ", nrow(ammod_raw), "\n")
cat("AMMOD >= 3mm OTUs:  ", nrow(ammod_3mm), "\n")
cat("AMMOD >= 5mm OTUs:  ", nrow(ammod_5mm), "\n")
cat(sprintf("  3mm retains %.1f%% of raw OTUs\n", nrow(ammod_3mm) / nrow(ammod_raw) * 100))
cat(sprintf("  5mm retains %.1f%% of raw OTUs\n\n", nrow(ammod_5mm) / nrow(ammod_raw) * 100))

# --- 4. BUILD DISTANCE MATRIX (SUPERFAMILY) ---
all_sf_taxa <- dplyr::bind_rows(faird_data, ammod_base) %>%
  dplyr::filter(!is.na(Taxon_SF)) %>%
  dplyr::distinct(Taxon_SF, .keep_all = TRUE) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  as.data.frame()

rownames(all_sf_taxa) <- all_sf_taxa$Taxon_SF
all_sf_taxa <- all_sf_taxa %>% dplyr::select(-Taxon_SF)

dist_matrix_sf <- as.matrix(vegan::taxa2dist(all_sf_taxa, varstep = TRUE, check = TRUE))

calculate_taxonomic_similarity <- function(taxa_A, taxa_B, dist_matrix) {
  taxa_A <- taxa_A[taxa_A %in% rownames(dist_matrix)]
  taxa_B <- taxa_B[taxa_B %in% rownames(dist_matrix)]
  if (length(taxa_A) == 0 | length(taxa_B) == 0) return(NA)
  min_A_to_B <- sapply(taxa_A, function(a) if(a %in% taxa_B) 0 else min(dist_matrix[a, taxa_B], na.rm = TRUE))
  min_B_to_A <- sapply(taxa_B, function(b) if(b %in% taxa_A) 0 else min(dist_matrix[b, taxa_A], na.rm = TRUE))
  return(1 - (mean(c(min_A_to_B, min_B_to_A)) / max(dist_matrix, na.rm = TRUE)))
}

# --- 5. WINDOWING AND COMPUTATION ---
all_dates <- seq.Date(date_start, date_end, by = "day")
window_breaks <- seq(date_start, date_end + 2, by = 3)
window_lookup <- tibble::tibble(Date = all_dates) %>%
  dplyr::mutate(Window_ID = as.integer(cut(Date, breaks = window_breaks, right = FALSE, labels = FALSE)))

sites <- c("Site1", "Site2", "Site3", "Site4")
results <- list()

for (s in sites) {
  for (w in unique(na.omit(window_lookup$Window_ID))) {
    win_dates <- window_lookup %>% dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
    
    faird_win <- faird_data %>% dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF))
    faird_sf <- unique(faird_win$Taxon_SF)
    n_faird_indiv <- nrow(faird_win)
    
    ammod_sf_raw <- unique((ammod_raw %>% dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF)))$Taxon_SF)
    ammod_sf_3mm <- unique((ammod_3mm %>% dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF)))$Taxon_SF)
    ammod_sf_5mm <- unique((ammod_5mm %>% dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF)))$Taxon_SF)
    
    passes_filter <- (n_faird_indiv >= MIN_INDIV)
    
    ds_raw <- NA_real_
    ds_3mm <- NA_real_
    ds_5mm <- NA_real_
    
    if (passes_filter & length(faird_sf) > 0) {
      if (length(ammod_sf_raw) > 0) ds_raw <- calculate_taxonomic_similarity(faird_sf, ammod_sf_raw, dist_matrix_sf)
      if (length(ammod_sf_3mm) > 0) ds_3mm <- calculate_taxonomic_similarity(faird_sf, ammod_sf_3mm, dist_matrix_sf)
      if (length(ammod_sf_5mm) > 0) ds_5mm <- calculate_taxonomic_similarity(faird_sf, ammod_sf_5mm, dist_matrix_sf)
    }
    
    results[[length(results) + 1]] <- tibble::tibble(
      Site = s, Window_ID = w, 
      Habitat = dplyr::if_else(s %in% c("Site1", "Site2"), "Maize", "Meadow"),
      DS_raw = ds_raw, DS_3mm = ds_3mm, DS_5mm = ds_5mm,
      Taxa_AMMOD_raw = length(ammod_sf_raw),
      Taxa_AMMOD_3mm = length(ammod_sf_3mm),
      Taxa_AMMOD_5mm = length(ammod_sf_5mm)
    )
  }
}

res_df <- dplyr::bind_rows(results) %>% dplyr::filter(!is.na(DS_raw) & !is.na(DS_3mm) & !is.na(DS_5mm))

# --- 6. SUMMARY TABLE (Table 16 format) ---
# Delta_raw_to_3mm: DS_3mm - DS_raw  (positive = filtering to 3mm improves concordance)
# Delta_3mm_to_5mm: DS_5mm - DS_3mm  (positive = filtering to 5mm improves concordance)
cat("--- TABLE 16: THRESHOLD SENSITIVITY (raw vs 3mm vs 5mm, SUPERFAMILY, 3-DAY) ---\n\n")

table16_by_site <- res_df %>%
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

table16_all <- res_df %>%
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

final_table <- dplyr::bind_rows(table16_by_site, table16_all)
print(final_table, n = Inf)

cat("\nNOTE: Delta_raw_to_3mm = DS_3mm - DS_raw (pp).",
    "Positive values indicate that restricting AMMOD to >=3mm improves concordance.\n")
cat("      Delta_3mm_to_5mm = DS_5mm - DS_3mm (pp).",
    "Same interpretation for the stricter threshold.\n")

# ==============================================================================
# SECTION 7: DIAGNOSTIC — DO SUB-3MM OTUS ADD NEW SUPERFAMILIES?
# ==============================================================================
# For each site × window, compute:
#   SF_raw    : Superfamilies in unfiltered AMMOD pool
#   SF_3mm    : Superfamilies in >=3mm AMMOD pool
#   SF_sub3mm : Superfamilies present ONLY in the sub-3mm fraction
#               (i.e., in raw but not in >=3mm)
#   SF_new    : SF_sub3mm taxa NOT found in SF_3mm (genuinely new Superfamilies)
#
# If SF_new is consistently 0 or near-0, the sub-3mm fraction adds no new
# Superfamilies — confirming that the >=3mm filter cannot improve concordance
# at the Superfamily level regardless of its effect on OTU count.
# ==============================================================================

cat("\n--- SECTION 7: SUB-3MM SUPERFAMILY DIAGNOSTIC ---\n\n")

diag_results <- list()

for (s in sites) {
  for (w in unique(na.omit(window_lookup$Window_ID))) {
    win_dates <- window_lookup %>% dplyr::filter(Window_ID == w) %>% dplyr::pull(Date)
    
    sf_raw <- unique((ammod_raw %>%
                        dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF)))$Taxon_SF)
    sf_3mm <- unique((ammod_3mm %>%
                        dplyr::filter(Date %in% win_dates, Site == s, !is.na(Taxon_SF)))$Taxon_SF)
    
    sf_sub3mm_only <- setdiff(sf_raw, sf_3mm)  # SF present in raw but absent in >=3mm
    sf_new         <- length(sf_sub3mm_only)    # = genuinely new SF from sub-3mm
    
    diag_results[[length(diag_results) + 1]] <- tibble::tibble(
      Site          = s,
      Window_ID     = w,
      Habitat       = dplyr::if_else(s %in% c("Site1", "Site2"), "Maize", "Meadow"),
      N_SF_raw      = length(sf_raw),
      N_SF_3mm      = length(sf_3mm),
      N_SF_sub3mm_exclusive = sf_new,
      Pct_new_SF    = dplyr::if_else(
        length(sf_raw) > 0,
        round(sf_new / length(sf_raw) * 100, 1),
        NA_real_
      ),
      SF_new_names  = paste(sf_sub3mm_only, collapse = "; ")
    )
  }
}

diag_df <- dplyr::bind_rows(diag_results) %>%
  dplyr::filter(N_SF_raw > 0)

cat("Windows with at least one AMMOD Superfamily detected:", nrow(diag_df), "\n\n")

# Summary by site
cat("--- Exclusive sub-3mm Superfamilies by site ---\n\n")
diag_summary <- diag_df %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n_windows          = dplyr::n(),
    mean_SF_raw        = round(mean(N_SF_raw), 1),
    mean_SF_3mm        = round(mean(N_SF_3mm), 1),
    mean_new_SF        = round(mean(N_SF_sub3mm_exclusive), 2),
    max_new_SF         = max(N_SF_sub3mm_exclusive),
    windows_with_new_SF = sum(N_SF_sub3mm_exclusive > 0),
    mean_pct_new       = round(mean(Pct_new_SF, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(diag_summary, n = Inf)

# Overall
cat("\n--- Overall ---\n")
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
print(diag_overall)

# List all unique exclusive SF names found
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
  cat(paste0("  ", all_new_sf), sep = "\n")
}


# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")