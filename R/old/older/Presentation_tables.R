# ==============================================================================
# SCRIPT: PRESENTATION_TABLES_V13
# Objective: Pre-filter raw data to discard contamination/invalid days, 
# apply correct AMMOD Batch 2 biomass rules, and sum total 'n' for summary.
# NOTE: Truth Calendar (effort_base, effort_summary) is now centralised in
#       load_data.R and accessed via data_list. Do NOT redefine it here.
# ==============================================================================

library(tidyverse)

# --- 0. GLOBAL SETTINGS ---
start_date <- as.Date("2023-08-23")
end_date   <- as.Date("2023-09-13")

# --- 1. LOAD DATA (includes Truth Calendar from load_data.R) ---
source("load_data.R")
load("data_cleaned.RData")

# Retrieve Truth Calendar from data_list (defined once in load_data.R)
effort_base    <- data_list$effort_base
effort_summary <- data_list$effort_summary

# --- 2. AGGREGATE RAW DETECTIONS (PRE-FILTERED BY TRUTH CALENDAR) ---
# effort_base comes from data_list (defined in load_data.R)

# A. E-traps (Biomass & Taxonomy separated to apply specific valid filters)
etraps_bio_clean <- etraps_daily_biomass %>%
  filter(Date >= start_date & Date <= end_date) %>%
  left_join(effort_base %>% select(Device, Date, Bio_Valid), by = c("Device", "Date")) %>%
  filter(Bio_Valid == TRUE) %>% # ONLY keep valid biomass days
  group_by(Device, Device_type) %>%
  summarise(
    Biomass_g         = sum(biomass_mg, na.rm = TRUE) / 1000,
    Bio_Days_Detected = n_distinct(Date[!is.na(biomass_mg) & biomass_mg > 0]),
    .groups = 'drop'
  )

etraps_tax_clean <- etraps_daily_biomass %>%
  filter(Date >= start_date & Date <= end_date) %>%
  left_join(effort_base %>% select(Device, Date, Taxo_Valid), by = c("Device", "Date")) %>%
  filter(Taxo_Valid == TRUE) %>% # ONLY keep valid taxonomy days
  group_by(Device, Device_type) %>%
  summarise(
    Individuals       = sum(abundance, na.rm = TRUE),
    Tax_Days_Detected = n_distinct(Date[!is.na(abundance) & abundance > 0]),
    .groups = 'drop'
  )

etraps_clean <- full_join(etraps_bio_clean, etraps_tax_clean, by = c("Device", "Device_type"))

# B. AMMOD
ammod_bio_clean <- ammod_biomass_for_stats %>%
  filter(Date >= start_date & Date <= end_date) %>%
  left_join(effort_base %>% select(Device, Date, Bio_Valid), by = c("Device", "Date")) %>%
  filter(Bio_Valid == TRUE) %>% # Drops Batch 2 and AMMOD3/4 ATL periods completely
  group_by(Device, Device_type) %>%
  summarise(
    Biomass_g         = sum(Live_mass, na.rm = TRUE) / 1000,
    Bio_Days_Detected = n_distinct(Date[!is.na(Live_mass) & Live_mass > 0]),
    .groups = 'drop'
  )

ammod_tax_clean <- ammod_taxonomy_presence %>%
  filter(Date >= start_date & Date <= end_date) %>%
  left_join(effort_base %>% select(Device, Date, Taxo_Valid), by = c("Device", "Date")) %>%
  filter(Taxo_Valid == TRUE) %>% # Drops the contaminated days (Aug 28-Sep 3) BEFORE counting!
  group_by(Device, Device_type) %>%
  summarise(
    Tax_Days_Detected = n_distinct(Date),
    .groups = 'drop'
  )

ammod_clean <- full_join(ammod_bio_clean, ammod_tax_clean, by = c("Device", "Device_type")) %>%
  mutate(Individuals = NA)

# Combine Raw Data for all devices
all_raw_data <- bind_rows(etraps_clean, ammod_clean)

# --- 3. MERGE DETAIL TABLE (Table 1) ---
final_detail_table <- effort_summary %>%
  left_join(all_raw_data, by = "Device") %>%
  # Fill NAs for devices that had zero detections or complete failures
  mutate(across(c(Biomass_g, Individuals, Bio_Days_Detected, Tax_Days_Detected), ~replace_na(., 0))) %>%
  mutate(
    Device_type = case_when(
      is.na(Device_type) & str_detect(Device, "FAIRD") ~ "FAIRD",
      is.na(Device_type) & str_detect(Device, "ID")    ~ "ID",
      is.na(Device_type) & str_detect(Device, "AMMOD") ~ "AMMOD",
      TRUE ~ Device_type
    ),
    # Calculate Success Rates based ONLY on valid lab days
    Taxo_Success_Pct = if_else(Valid_Taxo_Days > 0, round((Tax_Days_Detected / Valid_Taxo_Days) * 100, 1), 0),
    Bio_Success_Pct  = if_else(Valid_Bio_Days > 0, round((Bio_Days_Detected / Valid_Bio_Days) * 100, 1), 0)
  ) %>%
  select(Device, Device_type, Field_Days_Effort, 
         Valid_Taxo_Days, Tax_Days_Detected, Taxo_Success_Pct,
         Valid_Bio_Days, Bio_Days_Detected, Bio_Success_Pct,
         Biomass_g, Individuals) %>%
  arrange(Device_type, Device)

# --- 4. CREATE SUMMARY TYPE TABLE (Table 2) ---
summary_type_table <- final_detail_table %>%
  group_by(Device_type) %>%
  summarise(
    `Nr. Devices`               = n_distinct(Device),
    `Total Biomass (g)`         = round(sum(Biomass_g), 2),
    `Total Individuals`         = sum(Individuals),
    # Now calculating the SUM of all Valid Days to represent total 'n'
    `Total Valid Taxo Days (n)` = sum(Valid_Taxo_Days),
    `Total Valid Bio Days (n)`  = sum(Valid_Bio_Days),
    .groups = 'drop'
  )

# Calculate Taxonomic Richness STRICTLY on Valid_Taxo_Days
taxa_richness <- bind_rows(
  faird_data %>% mutate(Type = "FAIRD"),
  id_data %>% mutate(Type = "ID"),
  ammod_taxonomy_presence %>% mutate(Type = "AMMOD")
) %>%
  filter(Date >= start_date & Date <= end_date) %>%
  # Join with our Truth Calendar to exclude taxa from invalid days (lab errors/contamination)
  left_join(effort_base %>% select(Device, Date, Taxo_Valid), by = c("Device", "Date")) %>%
  filter(Taxo_Valid == TRUE) %>%
  group_by(Type) %>%
  summarise(
    `Unique Orders`   = n_distinct(Order, na.rm = TRUE),
    `Unique Families` = n_distinct(Family, na.rm = TRUE),
    `Unique Genera`   = n_distinct(Genus, na.rm = TRUE),
    .groups = 'drop'
  )

final_summary_table <- summary_type_table %>%
  left_join(taxa_richness, by = c("Device_type" = "Type"))

# --- OUTPUTS ---
cat("\n=== TABLE 1: DEVICE-LEVEL EFFORT VS SUCCESS ===\n")
print(as.data.frame(final_detail_table))

cat("\n=== TABLE 2: GENERAL CHARACTERIZATION ===\n")
print(as.data.frame(final_summary_table))

write_csv(final_detail_table, "outputs/Table_Device_Detail_V12.csv")
write_csv(final_summary_table, "outputs/Table_Summary_Type_V12.csv")