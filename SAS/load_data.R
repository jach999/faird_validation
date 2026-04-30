# ==============================================================================
# DATA LOADING SCRIPT - FAIR-Device Validation Project
# ==============================================================================
# Description: Loads and prepares all data files for analysis
# Date: 2025-10-10
# Updated: 2025-10-28 (JCF/Gemini) - Corrected completeness check
# ==============================================================================

# --- Load required libraries ---
library(tidyverse)
library(lubridate)
library(here)

cat("\n=== LOADING E-TRAPS DATA ===\n")

# ==============================================================================
# 1. LOAD E-TRAPS DATA (FAIRD AND ID)
# ==============================================================================

# --- Load FAIRD data ---
faird_data <- read_csv("source_tables/FAIRD_For_Statistics.csv",
                       show_col_types = FALSE) %>%
  dplyr::mutate(
    # Dates are already parsed correctly by read_csv
    Date = as.Date(Checkin),
    
    # Ensure Device_type is consistent
    Device_type = "FAIRD"
  ) %>%
  # Remove rows with no taxonomic information
  dplyr::filter(!is.na(Order))

cat("FAIRD data loaded:", nrow(faird_data), "individuals\n")
cat("Date range:", min(faird_data$Date), "to", max(faird_data$Date), "\n")
cat("Devices:", paste(unique(faird_data$Device), collapse = ", "), "\n")

# --- Load Insect Detect data ---
id_data <- read_csv("source_tables/ID_For_Statistics.csv",
                    show_col_types = FALSE) %>%
  dplyr::mutate(
    # Dates are already parsed correctly by read_csv
    Date = as.Date(Checkin),
    
    # Ensure Device_type is consistent
    Device_type = "ID"
  ) %>%
  # Remove rows with no taxonomic information
  dplyr::filter(!is.na(Order))

cat("ID data loaded:", nrow(id_data), "individuals\n")
cat("Date range:", min(id_data$Date), "to", max(id_data$Date), "\n")
cat("Devices:", paste(unique(id_data$Device), collapse = ", "), "\n")

# Check for ID3 (should be missing/minimal)
id3_count <- id_data %>% dplyr::filter(Device == "ID3") %>% nrow()
if(id3_count == 0) {
  cat("âš ï¸  ID3 has no data (device failure as expected)\n")
} else {
  cat("âš ï¸  ID3 has", id3_count, "records (unexpected - verify data)\n")
}

# --- Combine e-traps data ---
etraps_data <- bind_rows(faird_data, id_data) %>%
  dplyr::arrange(Date, Checkin)

cat("\nTotal e-traps individuals:", nrow(etraps_data), "\n")

# ==============================================================================
# 2. LOAD AMMOD TAXONOMIC DATA (METABARCODING)
# ==============================================================================

cat("\n=== LOADING AMMOD TAXONOMIC DATA ===\n")

ammod_taxonomy <- read_csv("source_tables/AMMOD_24h_OTUs_For_Statistics.csv",
                           show_col_types = FALSE,
                           na = c("", "NA", "#N/D", "#N/A")) %>%
  dplyr::mutate(
    # Date already comes as Date object from read_csv
    Device_type = "AMMOD"
  ) %>%
  # Remove rows with no taxonomic information
  dplyr::filter(!is.na(Order))

cat("AMMOD taxonomy loaded:", nrow(ammod_taxonomy), "OTUs/taxa\n")
cat("Date range:", min(ammod_taxonomy$Date), "to", max(ammod_taxonomy$Date), "\n")
cat("Devices:", paste(unique(ammod_taxonomy$Device), collapse = ", "), "\n")

# CRITICAL NOTE: Count column interpretation
cat("\nâš ï¸  IMPORTANT: AMMOD taxonomy is QUALITATIVE (presence/absence)\n")
cat("   'Count' column = DNA sequence reads (NOT individual insect counts)\n")
cat("   Use ONLY for presence/absence analyses\n")

# Convert to presence/absence for analysis
ammod_taxonomy_presence <- ammod_taxonomy %>%
  dplyr::mutate(Present = 1) %>%
  distinct(Date, Device, Device_type, Site, Ambient, Order, Family, Genus, .keep_all = TRUE) #%>%
#dplyr::select(-Count)  # Remove Count to avoid confusion

cat("   Converted to presence/absence format\n")

# ==============================================================================
# 2.5. TRUTH CALENDAR: VALIDITY FILTERS FOR ALL DEVICES
# ==============================================================================
# Centralised validity definition - all downstream scripts use this via data_list.
# This avoids each script redefining its own ad-hoc filters.

cat("\n=== BUILDING TRUTH CALENDAR ===\n")

experiment_start <- as.Date("2023-08-23")
experiment_end   <- as.Date("2023-09-13")

all_devices <- c("FAIRD1","FAIRD2","FAIRD3","FAIRD4",
                 "ID1","ID2","ID3","ID4",
                 "AMMOD1","AMMOD2","AMMOD3","AMMOD4")

effort_base <- tidyr::expand_grid(
  Device = all_devices,
  Date   = seq(experiment_start, experiment_end, by = "day")
) %>%
  dplyr::mutate(
    # --- Field validity: was the device physically operational? ---
    Field_Valid = dplyr::case_when(
      Device == "AMMOD4" & Date >= "2023-08-28" & Date <= "2023-09-03" ~ FALSE,  # Jammed bottle
      stringr::str_detect(Device, "AMMOD") & Date == "2023-09-04"      ~ FALSE,  # Bottle change day
      Device == "ID3"                                                   ~ FALSE,  # Corrupted SD card
      TRUE ~ TRUE
    ),
    # --- Taxonomy validity: field-valid AND lab taxonomy processing OK? ---
    Taxo_Valid = dplyr::case_when(
      !Field_Valid                                            ~ FALSE,
      Device == "AMMOD4" & Date == "2023-09-10"             ~ FALSE,  # Contaminated sample R2_4_06
      TRUE ~ TRUE
    ),
    # --- Biomass validity: field-valid AND gravimetric processing available? ---
    Bio_Valid = dplyr::case_when(
      !Field_Valid                                                                  ~ FALSE,
      Device == "AMMOD4"                                                           ~ FALSE,  # No biomass (Batch 1 & 2 ATL)
      stringr::str_detect(Device, "AMMOD") & Date >= "2023-09-05"                 ~ FALSE,  # All AMMOD Batch 2 ATL
      Device == "AMMOD3" & Date >= "2023-08-31" & Date <= "2023-09-03"            ~ FALSE,  # AMMOD3 partial Batch 1 ATL
      TRUE ~ TRUE
    )
  )

effort_summary <- effort_base %>%
  dplyr::group_by(Device) %>%
  dplyr::summarise(
    Field_Days_Effort = sum(Field_Valid),
    Valid_Taxo_Days   = sum(Taxo_Valid),
    Valid_Bio_Days    = sum(Bio_Valid),
    .groups = "drop"
  )

cat("Truth Calendar built for", length(all_devices), "devices over",
    as.numeric(experiment_end - experiment_start + 1), "days.\n")
print(effort_summary)

# Apply Taxo_Valid filter to ammod_taxonomy_presence
# This is the single authoritative filtered version used by all analysis scripts.
n_before <- nrow(ammod_taxonomy_presence)

ammod_taxonomy_presence <- ammod_taxonomy_presence %>%
  dplyr::left_join(
    effort_base %>% dplyr::select(Device, Date, Taxo_Valid),
    by = c("Device", "Date")
  ) %>%
  dplyr::filter(Taxo_Valid == TRUE) %>%
  dplyr::select(-Taxo_Valid)

n_after <- nrow(ammod_taxonomy_presence)
cat(sprintf("\nammod_taxonomy_presence filtered by Truth Calendar:\n"))
cat(sprintf("   Before: %d rows  |  After: %d rows  |  Removed: %d rows\n",
            n_before, n_after, n_before - n_after))
cat("   Removed rows correspond to AMMOD4 jammed-bottle period (Aug 28-Sep 3)\n")
cat("   and contaminated sample R2_4_06 (AMMOD4, Sep 10).\n")

# ==============================================================================
# 3. LOAD AMMOD BIOMASS DATA (GRAVIMETRIC)
# ==============================================================================

cat("\n=== LOADING AMMOD BIOMASS DATA ===\n")

ammod_biomass_raw <- read_csv("source_tables/AMMOD_24h_live_mass_For_Statistics.csv",
                              show_col_types = FALSE)

cat("Raw columns:", paste(colnames(ammod_biomass_raw), collapse = ", "), "\n")

# Process AMMOD biomass data
ammod_biomass <- ammod_biomass_raw %>%
  dplyr::mutate(
    # Date already comes as Date object from read_csv
    Device_type = "AMMOD",
    
    # Column is already named Live_mass (not Live_mass_mg)
    # Verify conversion: Live_mass should = Dry_mass_g * 1000 / 0.25 = Dry_mass_g * 4000
    Calculated_live_mass = Dry_mass_g * 4000,
    Conversion_check = abs(Live_mass - Calculated_live_mass) < 1
  ) %>%
  # Rename column with typo
  rename(DevicexAmbient = DevicexAmbien)

# Check conversion consistency
conversion_issues <- ammod_biomass %>%
  dplyr::filter(!is.na(Live_mass), !Conversion_check)

if(nrow(conversion_issues) > 0) {
  cat("âš ï¸  Warning: Found", nrow(conversion_issues), 
      "rows where conversion formula doesn't match\n")
  print(conversion_issues %>% dplyr::select(Date, Device, Dry_mass_g, Live_mass, Calculated_live_mass))
}

# Remove helper columns
ammod_biomass <- ammod_biomass %>%
  dplyr::select(-Calculated_live_mass, -Conversion_check)

cat("AMMOD biomass loaded:", nrow(ammod_biomass), "daily measurements\n")
cat("Date range:", min(ammod_biomass$Date, na.rm = TRUE), "to", 
    max(ammod_biomass$Date, na.rm = TRUE), "\n")
cat("Devices:", paste(unique(ammod_biomass$Device), collapse = ", "), "\n")
# ==============================================================================
# 4. DATA QUALITY CHECKS & DOCUMENTATION
# ==============================================================================

cat("\n=== AMMOD DATA QUALITY ASSESSMENT ===\n")

# Check completeness by device
ammod_completeness <- ammod_biomass %>%
  dplyr::group_by(Device, Site, Ambient) %>%
  dplyr::summarise(
    total_days = dplyr::n(),
    days_with_biomass = sum(!is.na(Live_mass)),
    days_missing = sum(is.na(Live_mass)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Device)

cat("\nData completeness by device:\n")
print(ammod_completeness)

# Detailed checks for incomplete devices
incomplete_devices <- ammod_completeness %>%
  dplyr::filter(days_missing > 0)

if(nrow(incomplete_devices) > 0) {
  cat("\nâš ï¸  DEVICES WITH INCOMPLETE DATA:\n")
  
  for(i in 1:nrow(incomplete_devices)) {
    device <- incomplete_devices$Device[i]
    site <- incomplete_devices$Site[i]
    ambient <- incomplete_devices$Ambient[i]
    missing <- incomplete_devices$days_missing[i]
    available <- incomplete_devices$days_with_biomass[i]
    
    cat(sprintf("\n%s (%s, %s):\n", device, site, ambient))
    cat(sprintf("  â€¢ Days available: %d\n", available))
    cat(sprintf("  â€¢ Days missing: %d\n", missing))
    
    # Device-specific explanations
    if(device == "AMMOD3") {
      cat("  â€¢ Status: Probes processing method changed in cursu (ATL) -> Starting from day 9, biomass not available\n")
      cat("  â€¢ Impact: Site3 for H2 validation (FAIRD-AMMOD) has reduced matching days (n=8)\n")
      cat("  â€¢ Note: OTU data available for H3 taxonomic analysis\n")
      
    } else if(device == "AMMOD4") {
      # Check which days are missing
      ammod4_dates <- ammod_biomass %>%
        dplyr::filter(Device == "AMMOD4") %>%
        dplyr::mutate(has_biomass = !is.na(Live_mass)) %>%
        dplyr::select(Date, has_biomass)
      
      missing_dates <- ammod4_dates %>% 
        dplyr::filter(!has_biomass) %>% 
        dplyr::pull(Date)
      
      cat("  â€¢ Status: Mechanical failure - rotation mechanism jammed (Aug 28-Sep 3)\n")
      cat(sprintf("  â€¢ Missing dates: %s to %s\n", 
                  min(missing_dates), max(missing_dates)))
      cat("  â€¢ Explanation: Bottle rotation system jammed on 28-Aug (bottle 6)\n")
      cat("  â€¢              Bottle 6 accumulated 7 days (28-Aug to 3-Sep)\n")
      cat("  â€¢              All biomass unavailable (ATL processing for entire Batch 1)\n")
      cat("  â€¢ Taxonomy: Batch 1 n=5 (Aug 23-27 valid) + Batch 2 n=9 (Sep 5-13) = 14 days\n")
      cat("  â€¢           Days Aug 28-Sep 3 excluded from H3 (accumulated sample)\n")
      cat("  â€¢ Impact: AMMOD4 excluded from H1 and H2 (no biomass data)\n")
      cat("  â€¢         AMMOD4 partially available for H3 (14 days taxonomy only)\n")
    }
  }
  cat("\n")
}

# Summary for analyses

cat("SUMMARY FOR STATISTICAL ANALYSES:\n")

# Create version with only valid measurements (NAs removed)
ammod_biomass_original <- ammod_biomass  # Keep all data including NAs
ammod_biomass_for_stats <- ammod_biomass %>%
  dplyr::filter(!is.na(Live_mass))

# Create version with only valid measurements (NAs removed)
ammod_biomass_for_stats <- ammod_biomass %>%
  dplyr::filter(!is.na(Live_mass))

by_site_summary <- ammod_biomass_for_stats %>%
  dplyr::group_by(Site, Device, Ambient) %>%
  dplyr::summarise(
    n_days = dplyr::n(),
    mean_biomass = mean(Live_mass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Site)

print(by_site_summary)

cat("\nTotal valid measurements:", nrow(ammod_biomass_for_stats), "\n")
cat("  â€¢ Site1 (Maize): n =", sum(ammod_biomass_for_stats$Site == "Site1"), "\n")
cat("  â€¢ Site2 (Maize): n =", sum(ammod_biomass_for_stats$Site == "Site2"), "\n")
cat("  â€¢ Site3 (Meadow): n =", sum(ammod_biomass_for_stats$Site == "Site3"), "\n")
cat("  â€¢ Site4 (Meadow): n =", sum(ammod_biomass_for_stats$Site == "Site4"), "\n")

cat("\nâœ… H2 FAIRD-AMMOD VALIDATION:\n")
cat("   Available sites: Site1 (n=12), Site2 (n=12), Site3 (n=8)\n")
cat("   Total sample size: n=32 (12+12+8)\n")

cat("\nâš ï¸  H1 AMMOD-AMMOD CONSISTENCY:\n")
ammod_pairs_summary <- ammod_biomass_for_stats %>%
  dplyr::group_by(Ambient) %>%
  dplyr::summarise(
    devices = paste(unique(Device), collapse = ", "),
    n_devices = dplyr::n_distinct(Device),
    total_obs = dplyr::n(),
    .groups = "drop"
  )

for(i in 1:nrow(ammod_pairs_summary)) {
  habitat <- ammod_pairs_summary$Ambient[i]
  n_dev <- ammod_pairs_summary$n_devices[i]
  devs <- ammod_pairs_summary$devices[i]
  n_obs <- ammod_pairs_summary$total_obs[i]
  
  if(n_dev == 2) {
    cat(sprintf("   %s pair: %s (n=%d total) âœ…\n", habitat, devs, n_obs))
  } else {
    cat(sprintf("   %s: Only %d device with biomass (%s) âŒ\n", 
                habitat, n_dev, devs))
  }
}

cat("\n")


# ==============================================================================
# 5. AGGREGATE E-TRAPS DATA TO DAILY BIOMASS
# ==============================================================================

cat("\n=== AGGREGATING E-TRAPS TO DAILY BIOMASS ===\n")

etraps_daily_biomass <- etraps_data %>%
  dplyr::group_by(Date, Device, Device_type, Site, Ambient) %>%
  dplyr::summarise(
    # Abundance (count of individuals)
    abundance = n(),
    
    # Biomass (sum of individual live masses)
    biomass_mg = sum(Live_mass, na.rm = TRUE),
    
    # Taxonomic richness
    n_orders = n_distinct(Order, na.rm = TRUE),
    n_families = n_distinct(Family, na.rm = TRUE),
    n_genera = n_distinct(Genus, na.rm = TRUE),
    
    # Mean body size
    mean_body_length = mean(Body_length, na.rm = TRUE),
    mean_body_width = mean(Body_width, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("E-traps daily aggregation complete\n")


# ==============================================================================
# 5.5. ++ NEW SECTION (CORRECTED): DATA COMPLETENESS ++
# ==============================================================================
# Purpose: Ensure all e-trap devices have one row for every day of the 
#          experiment, filling missing days with 0.

cat("\n=== ENFORCING DATA COMPLETENESS (E-TRAPS) ===\n")

# 1. Define the full date range *of the e-trap experiment*
#    We derive this from the loaded e-trap data
full_date_range <- seq.Date(min(etraps_data$Date, na.rm = TRUE), 
                            max(etraps_data$Date, na.rm = TRUE), 
                            by = "day")

# 2. Get all unique device metadata (excluding known failures like ID3)
device_metadata <- etraps_data %>%
  dplyr::filter(Device != "ID3") %>%
  dplyr::distinct(Device, Device_type, Site, Ambient)

cat(sprintf("Creating full data grid for %d devices over %d days...\n", 
            nrow(device_metadata), length(full_date_range)))

# 3. Create a "master grid" of all possible device-day combinations
master_grid <- tidyr::expand_grid(
  Date = full_date_range,
  device_metadata
)

# 4. Join the aggregated data to the master grid
#    This will create NA rows for any missing day/device combo
etraps_daily_biomass_complete <- master_grid %>%
  dplyr::left_join(
    etraps_daily_biomass,
    # Specify all columns in the 'device_metadata' to join by
    by = c("Date", "Device", "Device_type", "Site", "Ambient") 
  ) %>%
  
  # 5. Replace NAs in key numeric columns with 0 (true zeros)
  dplyr::mutate(
    abundance = replace_na(abundance, 0),
    biomass_mg = replace_na(biomass_mg, 0),
    n_orders = replace_na(n_orders, 0),
    n_families = replace_na(n_families, 0),
    n_genera = replace_na(n_genera, 0)
    # Note: mean_body_length/width are left as NA, as 0 is not meaningful
  )

# 6. Overwrite the old dataframe
etraps_daily_biomass <- etraps_daily_biomass_complete

cat("âœ… E-traps data completeness enforced.\n")
cat("   Missing days filled with abundance = 0 and biomass_mg = 0.\n")


# Check data days (NOW should be complete)
devices_to_check <- device_metadata$Device

for(dev in devices_to_check) {
  n_days <- etraps_daily_biomass %>% 
    dplyr::filter(Device == dev) %>% 
    distinct(Date) %>% 
    nrow()
  
  # This should now match length(full_date_range)
  cat(sprintf("%s: %d days of data (out of %d)\n", 
              dev, n_days, length(full_date_range)))
}


# ==============================================================================
# 6. CREATE TAXONOMIC PRESENCE/ABSENCE MATRICES
# ==============================================================================

cat("\n=== CREATING PRESENCE/ABSENCE MATRICES ===\n")

# E-traps presence/absence (daily)
etraps_presence_daily <- etraps_data %>%
  dplyr::group_by(Date, Device, Device_type, Site, Ambient, Order) %>%
  dplyr::summarise(present = 1, .groups = "drop")

# AMMOD presence/absence (daily) - using processed presence data
ammod_presence_daily <- ammod_taxonomy_presence %>%
  dplyr::select(Date, Device, Device_type, Site, Ambient, Order, Present) %>%
  rename(present = Present)

cat("Presence/absence matrices created\n")
cat("E-traps presence records:", nrow(etraps_presence_daily), "\n")
cat("AMMOD presence records:", nrow(ammod_presence_daily), "\n")

# ==============================================================================
# 7. CREATE ANALYSIS-READY DATASETS
# ==============================================================================

cat("\n=== CREATING ANALYSIS-READY DATASETS ===\n")

# Get overlapping date range
faird_date_range <- range(etraps_daily_biomass %>% 
                            dplyr::filter(Device_type == "FAIRD") %>% 
                            pull(Date))
ammod_date_range <- range(ammod_biomass_for_stats$Date, na.rm = TRUE)

overlap_start <- max(faird_date_range[1], ammod_date_range[1])
overlap_end <- min(faird_date_range[2], ammod_date_range[2])

cat("Overlapping date range for FAIRD-AMMOD comparison:\n")
cat("  Start:", as.character(overlap_start), "\n")
cat("  End:", as.character(overlap_end), "\n")
cat("  Duration:", as.numeric(overlap_end - overlap_start + 1), "days\n")

# Filter to overlapping dates
# NOTE: etraps_daily_biomass is already complete, this just subsets the dates
etraps_daily_overlap <- etraps_daily_biomass %>%
  dplyr::filter(Date >= overlap_start, Date <= overlap_end)

ammod_biomass_overlap <- ammod_biomass_for_stats %>%
  dplyr::filter(Date >= overlap_start, Date <= overlap_end)

# ==============================================================================
# 8. EXPORT CLEANED DATASETS FOR ANALYSIS
# ==============================================================================

cat("\n=== EXPORTING CLEANED DATASETS ===\n")

# Create a list with all datasets
data_list <- list(
  # Individual-level data
  faird = faird_data,
  id = id_data,
  etraps = etraps_data,
  
  # Daily aggregated data (NOW COMPLETE)
  etraps_daily = etraps_daily_biomass,
  etraps_daily_overlap = etraps_daily_overlap,
  
  # AMMOD data (not entangled)
  ammod_taxonomy = ammod_taxonomy_presence,   # Already filtered by Taxo_Valid (Truth Calendar)
  ammod_biomass_for_stats = ammod_biomass_for_stats,
  ammod_biomass_original = ammod_biomass_original,
  ammod_biomass = ammod_biomass_for_stats,    # Default to stats version
  ammod_biomass_overlap = ammod_biomass_overlap,
  
  # Presence/absence
  etraps_presence = etraps_presence_daily,
  ammod_presence = ammod_presence_daily,
  
  # Truth Calendar (centralised validity filters - use these in all analysis scripts)
  effort_base    = effort_base,     # Row-level validity flags (Field/Taxo/Bio) per device-day
  effort_summary = effort_summary,  # Aggregated valid-day counts per device
  
  # Metadata
  overlap_dates = list(start = overlap_start, end = overlap_end),
  experiment_dates = list(start = experiment_start, end = experiment_end)
)

# Save as RData for easy loading
save(data_list, file = "data_cleaned.RData")
cat("Cleaned data saved to: data_cleaned.RData\n")

# Also save individual CSVs
write_csv(etraps_daily_biomass, "etraps_daily_biomass_clean.csv")
write_csv(ammod_biomass_for_stats, "ammod_biomass_for_stats.csv")
write_csv(ammod_taxonomy_presence, "ammod_taxonomy_presence_clean.csv")
cat("Daily biomass and taxonomy files saved\n")

# ==============================================================================
# 9. SUMMARY STATISTICS
# ==============================================================================

cat("\n=== DATA SUMMARY ===\n")

summary_stats <- tibble(
  Dataset = c("FAIRD", "ID", "E-traps Total", "AMMOD Taxonomy", "AMMOD Biomass"),
  N_rows = c(
    nrow(faird_data),
    nrow(id_data),
    nrow(etraps_data),
    nrow(ammod_taxonomy_presence),
    nrow(ammod_biomass_for_stats)
  ),
  Date_range = c(
    paste(min(faird_data$Date), "to", max(faird_data$Date)),
    paste(min(id_data$Date), "to", max(id_data$Date)),
    paste(min(etraps_data$Date), "to", max(etraps_data$Date)),
    paste(min(ammod_taxonomy_presence$Date), "to", max(ammod_taxonomy_presence$Date)),
    paste(min(ammod_biomass_for_stats$Date, na.rm=TRUE), "to", 
          max(ammod_biomass_for_stats$Date, na.rm=TRUE))
  ),
  N_devices = c(
    n_distinct(faird_data$Device),
    n_distinct(id_data$Device),
    n_distinct(etraps_data$Device),
    n_distinct(ammod_taxonomy_presence$Device),
    n_distinct(ammod_biomass_for_stats$Device)
  )
)

print(summary_stats)

# Taxonomic summary
cat("\n--- Taxonomic Summary ---\n")
cat("E-traps:\n")
cat("  Orders:", n_distinct(etraps_data$Order, na.rm = TRUE), "\n")
cat("  Families:", n_distinct(etraps_data$Family, na.rm = TRUE), "\n")
cat("  Genera:", n_distinct(etraps_data$Genus, na.rm = TRUE), "\n")

cat("AMMOD:\n")
cat("  Orders:", n_distinct(ammod_taxonomy_presence$Order, na.rm = TRUE), "\n")
cat("  Families:", n_distinct(ammod_taxonomy_presence$Family, na.rm = TRUE), "\n")
cat("  Genera:", n_distinct(ammod_taxonomy_presence$Genus, na.rm = TRUE), "\n")

# Sites availability summary
cat("\n--- Sites Availability ---\n")
site_summary <- bind_rows(
  etraps_daily_biomass %>% 
    dplyr::group_by(Site, Device_type) %>% 
    dplyr::summarise(n_days = n_distinct(Date), .groups = "drop") %>%
    dplyr::mutate(Dataset = "E-traps"),
  ammod_biomass_for_stats %>%
    dplyr::filter(!is.na(Live_mass)) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(n_days = n_distinct(Date), .groups = "drop") %>%
    dplyr::mutate(Device_type = "AMMOD", Dataset = "AMMOD Biomass")
) %>%
  dplyr::arrange(Site, Device_type)

print(site_summary)

cat("\n=== DATA LOADING COMPLETE ===\n")
cat("Ready for analysis!\n")
cat("Access datasets via: data_list$faird, data_list$etraps_daily, etc.\n")