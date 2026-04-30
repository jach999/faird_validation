# ==============================================================================
# H2 DIVERGENCE DIAGNOSTIC: AMMOD-EXCLUSIVE TAXA ON HIGH-RESIDUAL DAYS
# ==============================================================================
# Purpose: Identify the taxonomic composition of what AMMOD detected but FAIRD
#          missed on days where the H2 log-log biomass model shows the largest
#          negative residuals (FAIRD underperformance relative to AMMOD).
#
# Rationale: Days with large negative residuals are the primary drivers of
#            weak R² in the H2 validation. Understanding what AMMOD captured
#            on those days (that FAIRD didn't) may reveal whether the gap is
#            driven by: (a) large/heavy insects entering AMMOD by chance
#            (lottery effect), (b) small swarming species below FAIRD's
#            detection threshold, or (c) other compositional differences.
#
# Output: A table of AMMOD-exclusive Superfamilies on high-residual days,
#         with full taxonomy down to Species for manual expert review.
#
# Author: JCA / Claude
# Date: 2026-03-28
# ==============================================================================

# --- Setup ---
library(tidyverse)
library(lubridate)

while (sink.number() > 0) sink()

output_file <- paste0("outputs/H2_divergence_diagnostic_output_",
                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
sink(output_file, split = TRUE)

cat("==============================================================================\n")
cat("H2 DIVERGENCE DIAGNOSTIC: AMMOD-EXCLUSIVE TAXA ON HIGH-RESIDUAL DAYS\n")
cat("==============================================================================\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Source data loading ---
source("load_data.R")
source("custom_colors.R")

# ==============================================================================
# SECTION 1: RECONSTRUCT H2 LOG-LOG RESIDUALS
# ==============================================================================

cat("\n=== SECTION 1: H2 LOG-LOG RESIDUALS ===\n\n")

# --- Prepare paired daily biomass (FAIRD vs AMMOD) per site ---
# Site-device mapping: Site1=FAIRD1/AMMOD1, Site2=FAIRD2/AMMOD2, Site3=FAIRD3/AMMOD3

faird_daily <- data_list$etraps_daily %>%
  dplyr::filter(Device_type == "FAIRD",
                Device %in% c("FAIRD1", "FAIRD2", "FAIRD3")) %>%
  dplyr::select(Date, Device, Site, biomass_mg) %>%
  dplyr::rename(FAIRD_biomass = biomass_mg)

ammod_daily <- data_list$ammod_biomass_for_stats %>%
  dplyr::filter(Device %in% c("AMMOD1", "AMMOD2", "AMMOD3"),
                !is.na(Live_mass)) %>%
  dplyr::select(Date, Device, Site, Live_mass) %>%
  dplyr::rename(AMMOD_biomass = Live_mass)

# Merge by Date + Site
paired_biomass <- dplyr::inner_join(
  faird_daily %>% dplyr::select(Date, Site, FAIRD_biomass),
  ammod_daily %>% dplyr::select(Date, Site, AMMOD_biomass),
  by = c("Date", "Site")
) %>%
  dplyr::mutate(
    log_FAIRD = log10(FAIRD_biomass + 1),
    log_AMMOD = log10(AMMOD_biomass + 1),
    Habitat   = dplyr::if_else(Site %in% c("Site1", "Site2"), "Maize", "Meadow")
  )

cat("Paired biomass observations:\n")
cat("  Total:", nrow(paired_biomass), "\n")
paired_biomass %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  print()
cat("\n")

# --- Fit log-log models per site and extract residuals ---
sites_h2 <- c("Site1", "Site2", "Site3")

residual_results <- list()

for (s in sites_h2) {
  site_data <- paired_biomass %>% dplyr::filter(Site == s)
  
  if (nrow(site_data) < 4) {
    cat(sprintf("  %s: Skipped (n=%d, too few observations)\n", s, nrow(site_data)))
    next
  }
  
  model <- lm(log_FAIRD ~ log_AMMOD, data = site_data)
  site_data$residual    <- residuals(model)
  site_data$fitted      <- fitted(model)
  site_data$std_residual <- rstandard(model)
  
  cat(sprintf("  %s (n=%d): R²=%.3f, Residual SD=%.3f\n",
              s, nrow(site_data),
              summary(model)$r.squared,
              sigma(model)))
  
  residual_results[[s]] <- site_data
}

all_residuals <- dplyr::bind_rows(residual_results)

cat(sprintf("\nTotal observations with residuals: %d\n", nrow(all_residuals)))

# ==============================================================================
# SECTION 2: IDENTIFY HIGH-DIVERGENCE DAYS
# ==============================================================================

cat("\n=== SECTION 2: HIGH-DIVERGENCE DAYS ===\n\n")

# Criterion: standardized residuals < -1 (approximately bottom ~16%)
# These are days where FAIRD biomass was substantially lower than predicted
# by AMMOD biomass through the site-specific log-log model.

RESIDUAL_THRESHOLD <- -1.0

divergent_days <- all_residuals %>%
  dplyr::filter(std_residual < RESIDUAL_THRESHOLD) %>%
  dplyr::arrange(std_residual)

cat(sprintf("Threshold: standardized residual < %.1f\n", RESIDUAL_THRESHOLD))
cat(sprintf("Days meeting criterion: %d / %d (%.1f%%)\n\n",
            nrow(divergent_days), nrow(all_residuals),
            100 * nrow(divergent_days) / nrow(all_residuals)))

cat("--- Divergent days detail ---\n\n")
divergent_days %>%
  dplyr::select(Date, Site, Habitat, FAIRD_biomass, AMMOD_biomass,
                log_FAIRD, log_AMMOD, residual, std_residual) %>%
  dplyr::mutate(
    FAIRD_biomass = round(FAIRD_biomass, 1),
    AMMOD_biomass = round(AMMOD_biomass, 1),
    log_FAIRD     = round(log_FAIRD, 3),
    log_AMMOD     = round(log_AMMOD, 3),
    residual      = round(residual, 3),
    std_residual  = round(std_residual, 3)
  ) %>%
  print(n = Inf)

# Also show full residual distribution for context
cat("\n--- Full residual distribution by site ---\n\n")
all_residuals %>%
  dplyr::group_by(Site, Habitat) %>%
  dplyr::summarise(
    n             = dplyr::n(),
    mean_std_res  = round(mean(std_residual), 3),
    sd_std_res    = round(sd(std_residual), 3),
    min_std_res   = round(min(std_residual), 3),
    max_std_res   = round(max(std_residual), 3),
    n_below_neg1  = sum(std_residual < RESIDUAL_THRESHOLD),
    n_above_pos1  = sum(std_residual > abs(RESIDUAL_THRESHOLD)),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# ==============================================================================
# SECTION 3: EXTRACT AMMOD-EXCLUSIVE SUPERFAMILIES ON DIVERGENT DAYS
# ==============================================================================

cat("\n=== SECTION 3: AMMOD-EXCLUSIVE SUPERFAMILIES ON DIVERGENT DAYS ===\n\n")

# Define non-classified marker
nc <- "not_classified"

# --- Helper: extract unique Superfamilies for a system on a given Date × Site ---
get_superfamilies <- function(data, date_val, site_val) {
  data %>%
    dplyr::filter(Date == date_val, Site == site_val,
                  !is.na(Superfamily), Superfamily != "", Superfamily != nc) %>%
    dplyr::pull(Superfamily) %>%
    unique()
}

# --- Prepare FAIRD individual-level data with Site mapping ---
faird_indiv <- data_list$faird %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2", "FAIRD3"))

# --- Prepare AMMOD taxonomy data ---
ammod_taxo <- data_list$ammod_taxonomy %>%
  dplyr::filter(Device %in% c("AMMOD1", "AMMOD2", "AMMOD3"))

# --- Process each divergent day ---
exclusive_results <- list()

for (i in seq_len(nrow(divergent_days))) {
  d   <- divergent_days$Date[i]
  s   <- divergent_days$Site[i]
  hab <- divergent_days$Habitat[i]
  res <- divergent_days$std_residual[i]
  
  # Get Superfamilies detected by each system
  faird_sfs <- get_superfamilies(faird_indiv, d, s)
  ammod_sfs <- get_superfamilies(ammod_taxo, d, s)
  
  # AMMOD-exclusive Superfamilies
  exclusive_sfs <- setdiff(ammod_sfs, faird_sfs)
  shared_sfs    <- intersect(ammod_sfs, faird_sfs)
  
  cat(sprintf("Day %s | %s (%s) | std.res=%.3f\n", d, s, hab, res))
  cat(sprintf("  FAIRD SFs: %d | AMMOD SFs: %d | Shared: %d | AMMOD-exclusive: %d\n",
              length(faird_sfs), length(ammod_sfs),
              length(shared_sfs), length(exclusive_sfs)))
  
  if (length(faird_sfs) > 0) {
    cat(sprintf("  FAIRD:     %s\n", paste(sort(faird_sfs), collapse = ", ")))
  } else {
    cat("  FAIRD:     [none detected]\n")
  }
  
  if (length(exclusive_sfs) > 0) {
    cat(sprintf("  Exclusive: %s\n", paste(sort(exclusive_sfs), collapse = ", ")))
  } else {
    cat("  Exclusive: [none — all AMMOD SFs also in FAIRD]\n")
  }
  cat("\n")
  
  # Extract full taxonomy for AMMOD-exclusive SFs
  if (length(exclusive_sfs) > 0) {
    excl_taxa <- ammod_taxo %>%
      dplyr::filter(Date == d, Site == s,
                    Superfamily %in% exclusive_sfs) %>%
      dplyr::select(Date, Site, Class, Order, Suborder, Infraorder,
                    Superfamily, Family, Subfamily, Genus, Species,
                    Body_length, Body_width, Live_mass) %>%
      dplyr::mutate(
        Std_Residual = round(res, 3),
        Habitat      = hab
      ) %>%
      dplyr::arrange(Order, Superfamily, Family, Genus, Species)
    
    exclusive_results[[length(exclusive_results) + 1]] <- excl_taxa
  }
}

# ==============================================================================
# SECTION 4: COMPREHENSIVE OUTPUT TABLE
# ==============================================================================

cat("\n=== SECTION 4: COMPREHENSIVE AMMOD-EXCLUSIVE TAXA TABLE ===\n\n")

if (length(exclusive_results) > 0) {
  excl_all <- dplyr::bind_rows(exclusive_results) %>%
    dplyr::arrange(Date, Site, Order, Superfamily, Family, Genus, Species)
  
  cat(sprintf("Total AMMOD-exclusive OTUs across all divergent days: %d\n", nrow(excl_all)))
  cat(sprintf("Unique Superfamilies: %d\n",
              dplyr::n_distinct(excl_all$Superfamily)))
  cat(sprintf("Unique Families: %d\n",
              dplyr::n_distinct(excl_all$Family[!is.na(excl_all$Family) &
                                                  excl_all$Family != "" &
                                                  excl_all$Family != nc])))
  cat(sprintf("Unique Genera: %d\n",
              dplyr::n_distinct(excl_all$Genus[!is.na(excl_all$Genus) &
                                                 excl_all$Genus != "" &
                                                 excl_all$Genus != nc])))
  cat(sprintf("Unique Species: %d\n\n",
              dplyr::n_distinct(excl_all$Species[!is.na(excl_all$Species) &
                                                   excl_all$Species != "" &
                                                   excl_all$Species != nc])))
  
  # --- Print full table ---
  cat("--- FULL TABLE: AMMOD-EXCLUSIVE OTUs ON HIGH-RESIDUAL DAYS ---\n\n")
  print(excl_all, n = Inf, width = Inf)
  
  # --- Summary by Superfamily ---
  cat("\n\n--- SUMMARY BY SUPERFAMILY (across all divergent days) ---\n\n")
  sf_summary <- excl_all %>%
    dplyr::group_by(Order, Superfamily) %>%
    dplyr::summarise(
      n_otus            = dplyr::n(),
      n_days_appearing  = dplyr::n_distinct(Date),
      n_sites           = dplyr::n_distinct(Site),
      n_families        = dplyr::n_distinct(Family[!is.na(Family) &
                                                     Family != "" & Family != nc]),
      n_genera          = dplyr::n_distinct(Genus[!is.na(Genus) &
                                                    Genus != "" & Genus != nc]),
      n_species         = dplyr::n_distinct(Species[!is.na(Species) &
                                                      Species != "" & Species != nc]),
      mean_BL_mm        = round(mean(Body_length, na.rm = TRUE), 2),
      min_BL_mm         = round(min(Body_length, na.rm = TRUE), 2),
      max_BL_mm         = round(max(Body_length, na.rm = TRUE), 2),
      example_genera    = paste(unique(Genus[!is.na(Genus) &
                                               Genus != "" & Genus != nc]),
                                collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Order, Superfamily)
  
  print(sf_summary, n = Inf, width = Inf)
  
  # --- Save CSV for manual review ---
  csv_path <- paste0("outputs/H2_divergence_AMMOD_exclusive_taxa_",
                     format(Sys.time(), "%Y%m%d"), ".csv")
  readr::write_csv(excl_all, csv_path)
  cat(sprintf("\n[OK] Full table saved to: %s\n", csv_path))
  
  csv_summary_path <- paste0("outputs/H2_divergence_AMMOD_exclusive_summary_",
                             format(Sys.time(), "%Y%m%d"), ".csv")
  readr::write_csv(sf_summary, csv_summary_path)
  cat(sprintf("[OK] Summary saved to: %s\n", csv_summary_path))
  
} else {
  cat("No AMMOD-exclusive Superfamilies found on divergent days.\n")
  cat("This would mean FAIRD detected all Superfamilies that AMMOD detected\n")
  cat("on these days — the divergence is driven by shared taxa, not missing ones.\n")
}

# ==============================================================================
# SECTION 5: ALSO EXAMINE HIGH POSITIVE RESIDUALS (COMPLEMENTARY)
# ==============================================================================

cat("\n=== SECTION 5: HIGH POSITIVE RESIDUALS (FAIRD 'OVERPERFORMANCE') ===\n\n")
cat("For completeness: days where FAIRD captured MORE biomass than expected.\n")
cat("These may reveal days where a large insect entered FAIRD but not AMMOD.\n\n")

positive_divergent <- all_residuals %>%
  dplyr::filter(std_residual > abs(RESIDUAL_THRESHOLD)) %>%
  dplyr::arrange(dplyr::desc(std_residual))

cat(sprintf("Days with std.residual > +%.1f: %d\n\n", abs(RESIDUAL_THRESHOLD),
            nrow(positive_divergent)))

if (nrow(positive_divergent) > 0) {
  positive_divergent %>%
    dplyr::select(Date, Site, Habitat, FAIRD_biomass, AMMOD_biomass,
                  std_residual) %>%
    dplyr::mutate(
      FAIRD_biomass = round(FAIRD_biomass, 1),
      AMMOD_biomass = round(AMMOD_biomass, 1),
      std_residual  = round(std_residual, 3)
    ) %>%
    print(n = Inf)
  
  # For positive residuals: what did FAIRD detect that AMMOD didn't?
  cat("\n--- FAIRD-exclusive Superfamilies on positive-residual days ---\n\n")
  
  for (i in seq_len(nrow(positive_divergent))) {
    d   <- positive_divergent$Date[i]
    s   <- positive_divergent$Site[i]
    res <- positive_divergent$std_residual[i]
    
    faird_sfs <- get_superfamilies(faird_indiv, d, s)
    ammod_sfs <- get_superfamilies(ammod_taxo, d, s)
    
    faird_excl <- setdiff(faird_sfs, ammod_sfs)
    
    cat(sprintf("Day %s | %s | std.res=+%.3f\n", d, s, res))
    cat(sprintf("  FAIRD SFs: %d | AMMOD SFs: %d | FAIRD-exclusive: %d\n",
                length(faird_sfs), length(ammod_sfs), length(faird_excl)))
    
    if (length(faird_excl) > 0) {
      # List the FAIRD-exclusive taxa with body length
      faird_excl_taxa <- faird_indiv %>%
        dplyr::filter(Date == d, Site == s,
                      Superfamily %in% faird_excl) %>%
        dplyr::select(Superfamily, Family, Genus, Species, Body_length) %>%
        dplyr::arrange(Superfamily, Family, Genus)
      
      cat(sprintf("  FAIRD-exclusive SFs: %s\n",
                  paste(sort(faird_excl), collapse = ", ")))
      cat("  Taxa detail:\n")
      print(faird_excl_taxa, n = Inf)
    } else {
      cat("  FAIRD-exclusive: [none]\n")
    }
    cat("\n")
  }
}

# ==============================================================================
# SECTION 6: CONTEXT — ALL RESIDUALS TABLE
# ==============================================================================

cat("\n=== SECTION 6: FULL RESIDUALS TABLE (for reference) ===\n\n")

all_residuals %>%
  dplyr::select(Date, Site, Habitat, FAIRD_biomass, AMMOD_biomass,
                log_FAIRD, log_AMMOD, residual, std_residual) %>%
  dplyr::mutate(
    FAIRD_biomass = round(FAIRD_biomass, 1),
    AMMOD_biomass = round(AMMOD_biomass, 1),
    log_FAIRD     = round(log_FAIRD, 3),
    log_AMMOD     = round(log_AMMOD, 3),
    residual      = round(residual, 3),
    std_residual  = round(std_residual, 3),
    flag          = dplyr::case_when(
      std_residual < RESIDUAL_THRESHOLD  ~ "<<< NEGATIVE",
      std_residual > abs(RESIDUAL_THRESHOLD) ~ ">>> POSITIVE",
      TRUE ~ ""
    )
  ) %>%
  dplyr::arrange(Site, Date) %>%
  print(n = Inf, width = Inf)

# ==============================================================================
# CLOSE LOG
# ==============================================================================
sink()
cat(">> Output saved to:", output_file, "\n")