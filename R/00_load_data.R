# ==============================================================================
# DATA LOADING — FAIR-Device Validation Project
# ==============================================================================
# All file paths use here::here() so the script resolves correctly whether
# sourced from the project root via .Rproj or run directly from R/.
# ==============================================================================

library(tidyverse)
library(here)

# Anchor here to the project root regardless of where this script is executed from.
# here::i_am() walks upward from the script location until it finds a directory
# that contains R/00_load_data.R — that directory becomes the project root.
here::i_am("R/00_load_data.R")

cat("\n=== LOADING E-TRAPS DATA ===\n")

# ==============================================================================
# 1. LOAD E-TRAPS DATA (FAIRD AND ID)
# ==============================================================================

faird_data <- readr::read_csv(
  here::here("source_tables", "FAIRD_For_Statistics.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Date        = as.Date(Checkin),
    Device_type = "FAIRD"
  ) %>%
  dplyr::filter(!is.na(Order))

cat("FAIRD loaded:", nrow(faird_data), "individuals\n")
cat("Date range:", as.character(min(faird_data$Date)),
    "to", as.character(max(faird_data$Date)), "\n")
cat("Devices:", paste(sort(unique(faird_data$Device)), collapse = ", "), "\n")

id_data <- readr::read_csv(
  here::here("source_tables", "ID_For_Statistics.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Date        = as.Date(Checkin),
    Device_type = "ID"
  ) %>%
  dplyr::filter(!is.na(Order))

cat("ID loaded:", nrow(id_data), "individuals\n")
cat("Date range:", as.character(min(id_data$Date)),
    "to", as.character(max(id_data$Date)), "\n")
cat("Devices:", paste(sort(unique(id_data$Device)), collapse = ", "), "\n")

id3_count <- id_data %>% dplyr::filter(Device == "ID3") %>% nrow()
if (id3_count == 0) {
  cat("ID3 has no data (device failure as expected)\n")
} else {
  cat("WARNING: ID3 has", id3_count, "records (unexpected - verify data)\n")
}

etraps_data <- dplyr::bind_rows(faird_data, id_data) %>%
  dplyr::arrange(Date, Checkin)

cat("\nTotal e-traps individuals:", nrow(etraps_data), "\n")


# ==============================================================================
# 2. LOAD AMMOD TAXONOMIC DATA (METABARCODING)
# ==============================================================================

cat("\n=== LOADING AMMOD TAXONOMIC DATA ===\n")

ammod_taxonomy <- readr::read_csv(
  here::here("source_tables", "AMMOD_24h_OTUs_For_Statistics.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "#N/D", "#N/A")
) %>%
  dplyr::mutate(Device_type = "AMMOD") %>%
  dplyr::filter(!is.na(Order))

cat("AMMOD taxonomy loaded:", nrow(ammod_taxonomy), "OTUs\n")
cat("Date range:", as.character(min(ammod_taxonomy$Date)),
    "to", as.character(max(ammod_taxonomy$Date)), "\n")
cat("Devices:", paste(sort(unique(ammod_taxonomy$Device)), collapse = ", "), "\n")

cat("\nIMPORTANT: AMMOD taxonomy is QUALITATIVE (presence/absence)\n")
cat("  'Count' = DNA sequence reads, NOT individual insect counts.\n")
cat("  Use ONLY for presence/absence analyses.\n")

# Presence/absence deduplication — Count column intentionally discarded
ammod_taxonomy_presence <- ammod_taxonomy %>%
  dplyr::mutate(Present = 1) %>%
  dplyr::distinct(Date, Device, Device_type, Site, Ambient,
                  Order, Family, Genus, .keep_all = TRUE)

cat("Converted to presence/absence format\n")


# ==============================================================================
# 2.5. TRUTH CALENDAR: CENTRALISED VALIDITY FLAGS
# ==============================================================================

cat("\n=== BUILDING TRUTH CALENDAR ===\n")

experiment_start <- as.Date("2023-08-23")
experiment_end   <- as.Date("2023-09-13")

all_devices <- c("FAIRD1", "FAIRD2", "FAIRD3", "FAIRD4",
                 "ID1",    "ID2",    "ID3",    "ID4",
                 "AMMOD1", "AMMOD2", "AMMOD3", "AMMOD4")

effort_base <- tidyr::expand_grid(
  Device = all_devices,
  Date   = seq(experiment_start, experiment_end, by = "day")
) %>%
  dplyr::mutate(
    Field_Valid = dplyr::case_when(
      Device == "AMMOD4" & Date >= "2023-08-28" & Date <= "2023-09-03" ~ FALSE,
      stringr::str_detect(Device, "AMMOD") & Date == "2023-09-04"      ~ FALSE,
      Device == "ID3"                                                   ~ FALSE,
      TRUE ~ TRUE
    ),
    Taxo_Valid = dplyr::case_when(
      !Field_Valid                                        ~ FALSE,
      Device == "AMMOD4" & Date == "2023-09-10"         ~ FALSE,
      TRUE ~ TRUE
    ),
    Bio_Valid = dplyr::case_when(
      !Field_Valid                                                              ~ FALSE,
      Device == "AMMOD4"                                                       ~ FALSE,
      stringr::str_detect(Device, "AMMOD") & Date >= "2023-09-05"             ~ FALSE,
      Device == "AMMOD3" & Date >= "2023-08-31" & Date <= "2023-09-03"        ~ FALSE,
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

# Apply Taxo_Valid to ammod_taxonomy_presence — single authoritative filter
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
cat(sprintf("  Before: %d | After: %d | Removed: %d\n",
            n_before, n_after, n_before - n_after))
cat("  Removed: AMMOD4 jammed-bottle period (Aug 28-Sep 3) +",
    "contaminated sample R2_4_06 (AMMOD4, Sep 10).\n")


# ==============================================================================
# 3. LOAD AMMOD BIOMASS DATA (GRAVIMETRIC)
# ==============================================================================

cat("\n=== LOADING AMMOD BIOMASS DATA ===\n")

ammod_biomass_raw <- readr::read_csv(
  here::here("source_tables", "AMMOD_24h_live_mass_For_Statistics.csv"),
  show_col_types = FALSE
)

cat("Raw columns:", paste(colnames(ammod_biomass_raw), collapse = ", "), "\n")

ammod_biomass <- ammod_biomass_raw %>%
  dplyr::mutate(
    Device_type          = "AMMOD",
    Calculated_live_mass = Dry_mass_g * 4000,
    Conversion_check     = abs(Live_mass - Calculated_live_mass) < 1
  ) %>%
  dplyr::rename(DevicexAmbient = DevicexAmbien)

conversion_issues <- ammod_biomass %>%
  dplyr::filter(!is.na(Live_mass), !Conversion_check)

if (nrow(conversion_issues) > 0) {
  cat("WARNING:", nrow(conversion_issues),
      "rows where conversion formula does not match\n")
  print(conversion_issues %>%
          dplyr::select(Date, Device, Dry_mass_g, Live_mass, Calculated_live_mass))
}

ammod_biomass <- ammod_biomass %>%
  dplyr::select(-Calculated_live_mass, -Conversion_check)

cat("AMMOD biomass loaded:", nrow(ammod_biomass), "daily measurements\n")
cat("Date range:", as.character(min(ammod_biomass$Date, na.rm = TRUE)),
    "to", as.character(max(ammod_biomass$Date, na.rm = TRUE)), "\n")
cat("Devices:", paste(sort(unique(ammod_biomass$Device)), collapse = ", "), "\n")


# ==============================================================================
# 4. AMMOD DATA QUALITY CHECKS
# ==============================================================================

cat("\n=== AMMOD DATA QUALITY ASSESSMENT ===\n")

ammod_completeness <- ammod_biomass %>%
  dplyr::group_by(Device, Site, Ambient) %>%
  dplyr::summarise(
    total_days        = dplyr::n(),
    days_with_biomass = sum(!is.na(Live_mass)),
    days_missing      = sum(is.na(Live_mass)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Device)

cat("\nCompleteness by device:\n")
print(ammod_completeness)

incomplete_devices <- ammod_completeness %>% dplyr::filter(days_missing > 0)

if (nrow(incomplete_devices) > 0) {
  cat("\nDEVICES WITH INCOMPLETE DATA:\n")

  for (i in seq_len(nrow(incomplete_devices))) {
    device    <- incomplete_devices$Device[i]
    site      <- incomplete_devices$Site[i]
    ambient   <- incomplete_devices$Ambient[i]
    missing   <- incomplete_devices$days_missing[i]
    available <- incomplete_devices$days_with_biomass[i]

    cat(sprintf("\n%s (%s, %s):\n", device, site, ambient))
    cat(sprintf("  Days available: %d | Days missing: %d\n", available, missing))

    if (device == "AMMOD3") {
      cat("  Status: Processing method changed in cursu (ATL) from day 9.\n")
      cat("  Impact: Site3 for H2 has reduced matching days (n=8).\n")
      cat("  Note:   OTU data available for H1 qualitative analysis.\n")
    } else if (device == "AMMOD4") {
      missing_dates <- ammod_biomass %>%
        dplyr::filter(Device == "AMMOD4", is.na(Live_mass)) %>%
        dplyr::pull(Date)
      cat(sprintf("  Status: Mechanical failure — rotation jammed (Aug 28-Sep 3).\n"))
      cat(sprintf("  Missing dates: %s to %s\n",
                  as.character(min(missing_dates)),
                  as.character(max(missing_dates))))
      cat("  All biomass unavailable (ATL Batch 1). Excluded from H2.\n")
    }
  }
  cat("\n")
}

ammod_biomass_original  <- ammod_biomass
ammod_biomass_for_stats <- ammod_biomass %>% dplyr::filter(!is.na(Live_mass))

by_site_summary <- ammod_biomass_for_stats %>%
  dplyr::group_by(Site, Device, Ambient) %>%
  dplyr::summarise(
    n_days       = dplyr::n(),
    mean_biomass = mean(Live_mass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Site)

cat("\nSummary for statistical analyses:\n")
print(by_site_summary)
cat("\nTotal valid measurements:", nrow(ammod_biomass_for_stats), "\n")
cat("  Site1 (Maize):  n =", sum(ammod_biomass_for_stats$Site == "Site1"), "\n")
cat("  Site2 (Maize):  n =", sum(ammod_biomass_for_stats$Site == "Site2"), "\n")
cat("  Site3 (Meadow): n =", sum(ammod_biomass_for_stats$Site == "Site3"), "\n")
cat("  Site4 (Meadow): n =", sum(ammod_biomass_for_stats$Site == "Site4"), "\n")


# ==============================================================================
# 5. AGGREGATE E-TRAPS TO DAILY BIOMASS
# ==============================================================================

cat("\n=== AGGREGATING E-TRAPS TO DAILY BIOMASS ===\n")

etraps_daily_biomass <- etraps_data %>%
  dplyr::group_by(Date, Device, Device_type, Site, Ambient) %>%
  dplyr::summarise(
    abundance         = dplyr::n(),
    biomass_mg        = sum(Live_mass,    na.rm = TRUE),
    n_orders          = dplyr::n_distinct(Order,  na.rm = TRUE),
    n_families        = dplyr::n_distinct(Family, na.rm = TRUE),
    n_genera          = dplyr::n_distinct(Genus,  na.rm = TRUE),
    mean_body_length  = mean(Body_length, na.rm = TRUE),
    mean_body_width   = mean(Body_width,  na.rm = TRUE),
    .groups = "drop"
  )

cat("E-traps daily aggregation complete\n")


# ==============================================================================
# 5.5. ENFORCE DATA COMPLETENESS (zero-capture days as real observations)
# ==============================================================================

cat("\n=== ENFORCING DATA COMPLETENESS (E-TRAPS) ===\n")

full_date_range <- seq.Date(min(etraps_data$Date, na.rm = TRUE),
                            max(etraps_data$Date, na.rm = TRUE),
                            by = "day")

device_metadata <- etraps_data %>%
  dplyr::filter(Device != "ID3") %>%
  dplyr::distinct(Device, Device_type, Site, Ambient)

cat(sprintf("Full grid: %d devices x %d days\n",
            nrow(device_metadata), length(full_date_range)))

master_grid <- tidyr::expand_grid(Date = full_date_range, device_metadata)

etraps_daily_biomass <- master_grid %>%
  dplyr::left_join(etraps_daily_biomass,
                   by = c("Date", "Device", "Device_type", "Site", "Ambient")) %>%
  dplyr::mutate(
    abundance  = tidyr::replace_na(abundance,  0),
    biomass_mg = tidyr::replace_na(biomass_mg, 0),
    n_orders   = tidyr::replace_na(n_orders,   0),
    n_families = tidyr::replace_na(n_families, 0),
    n_genera   = tidyr::replace_na(n_genera,   0)
  )

cat("Zero-capture days filled (abundance = 0, biomass_mg = 0).\n")
cat("mean_body_length/width left as NA on zero-capture days.\n\n")

for (dev in device_metadata$Device) {
  n_days <- etraps_daily_biomass %>%
    dplyr::filter(Device == dev) %>%
    dplyr::distinct(Date) %>%
    nrow()
  cat(sprintf("  %s: %d days (expected %d)\n",
              dev, n_days, length(full_date_range)))
}


# ==============================================================================
# 6. PRESENCE/ABSENCE MATRICES
# ==============================================================================

cat("\n=== CREATING PRESENCE/ABSENCE MATRICES ===\n")

etraps_presence_daily <- etraps_data %>%
  dplyr::group_by(Date, Device, Device_type, Site, Ambient, Order) %>%
  dplyr::summarise(present = 1, .groups = "drop")

ammod_presence_daily <- ammod_taxonomy_presence %>%
  dplyr::select(Date, Device, Device_type, Site, Ambient, Order, Present) %>%
  dplyr::rename(present = Present)

cat("E-traps presence records:", nrow(etraps_presence_daily), "\n")
cat("AMMOD presence records:",   nrow(ammod_presence_daily), "\n")


# ==============================================================================
# 7. ANALYSIS-READY DATASETS
# ==============================================================================

cat("\n=== CREATING ANALYSIS-READY DATASETS ===\n")

faird_date_range <- range(etraps_daily_biomass %>%
                            dplyr::filter(Device_type == "FAIRD") %>%
                            dplyr::pull(Date))
ammod_date_range <- range(ammod_biomass_for_stats$Date, na.rm = TRUE)

overlap_start <- max(faird_date_range[1], ammod_date_range[1])
overlap_end   <- min(faird_date_range[2], ammod_date_range[2])

cat("FAIRD-AMMOD overlapping date range:\n")
cat("  Start:", as.character(overlap_start), "\n")
cat("  End:",   as.character(overlap_end),   "\n")
cat("  Days:",  as.numeric(overlap_end - overlap_start + 1), "\n")

etraps_daily_overlap  <- etraps_daily_biomass %>%
  dplyr::filter(Date >= overlap_start, Date <= overlap_end)
ammod_biomass_overlap <- ammod_biomass_for_stats %>%
  dplyr::filter(Date >= overlap_start, Date <= overlap_end)


# ==============================================================================
# 8. PACKAGE INTO data_list AND SAVE
# ==============================================================================

cat("\n=== EXPORTING CLEANED DATASETS ===\n")

data_list <- list(
  faird                  = faird_data,
  id                     = id_data,
  etraps                 = etraps_data,
  etraps_daily           = etraps_daily_biomass,
  etraps_daily_overlap   = etraps_daily_overlap,
  ammod_taxonomy         = ammod_taxonomy_presence,
  ammod_biomass_for_stats= ammod_biomass_for_stats,
  ammod_biomass_original = ammod_biomass_original,
  ammod_biomass          = ammod_biomass_for_stats,
  ammod_biomass_overlap  = ammod_biomass_overlap,
  etraps_presence        = etraps_presence_daily,
  ammod_presence         = ammod_presence_daily,
  effort_base            = effort_base,
  effort_summary         = effort_summary,
  overlap_dates          = list(start = overlap_start, end = overlap_end),
  experiment_dates       = list(start = experiment_start, end = experiment_end)
)

save(data_list, file = here::here("data_cleaned.RData"))
cat("data_list saved to: data_cleaned.RData\n")


# ==============================================================================
# 9. SUMMARY
# ==============================================================================

cat("\n=== DATA SUMMARY ===\n")

summary_stats <- tibble::tibble(
  Dataset  = c("FAIRD", "ID", "E-traps Total",
               "AMMOD Taxonomy", "AMMOD Biomass"),
  N_rows   = c(nrow(faird_data), nrow(id_data), nrow(etraps_data),
               nrow(ammod_taxonomy_presence), nrow(ammod_biomass_for_stats)),
  N_devices = c(dplyr::n_distinct(faird_data$Device),
                dplyr::n_distinct(id_data$Device),
                dplyr::n_distinct(etraps_data$Device),
                dplyr::n_distinct(ammod_taxonomy_presence$Device),
                dplyr::n_distinct(ammod_biomass_for_stats$Device))
)
print(summary_stats)

cat("\n--- Taxonomic richness ---\n")
cat("E-traps: Orders =",   dplyr::n_distinct(etraps_data$Order,  na.rm = TRUE),
    "| Families =", dplyr::n_distinct(etraps_data$Family, na.rm = TRUE),
    "| Genera =",   dplyr::n_distinct(etraps_data$Genus,  na.rm = TRUE), "\n")
cat("AMMOD:   Orders =",   dplyr::n_distinct(ammod_taxonomy_presence$Order,  na.rm = TRUE),
    "| Families =", dplyr::n_distinct(ammod_taxonomy_presence$Family, na.rm = TRUE),
    "| Genera =",   dplyr::n_distinct(ammod_taxonomy_presence$Genus,  na.rm = TRUE), "\n")

cat("\n--- Site availability ---\n")
dplyr::bind_rows(
  etraps_daily_biomass %>%
    dplyr::group_by(Site, Device_type) %>%
    dplyr::summarise(n_days = dplyr::n_distinct(Date), .groups = "drop") %>%
    dplyr::mutate(Dataset = "E-traps"),
  ammod_biomass_for_stats %>%
    dplyr::filter(!is.na(Live_mass)) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(n_days = dplyr::n_distinct(Date), .groups = "drop") %>%
    dplyr::mutate(Device_type = "AMMOD", Dataset = "AMMOD Biomass")
) %>%
  dplyr::arrange(Site, Device_type) %>%
  print()

cat("\n=== DATA LOADING COMPLETE ===\n")
cat("Access datasets via: data_list$faird, data_list$etraps_daily, etc.\n")
