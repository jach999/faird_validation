# -----------------------------------------------------------------------------
# Script: prepare_sas_data.R
# Description: Generates two clean, pre-filtered CSV files for SAS.
#              1. Biomass Validation (Days 1-12 only, Sites 1-3 only)
#              2. Abundance Comparison (FAIRD vs ID, all days, all sites)
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)

# 1. Load master data
# Ensure 'data_cleaned.RData' is in your working directory
load("data_cleaned.RData")

# 2. Prepare the Master Dataframe
# This combines all necessary columns and calculates log transformations
df_master <- bind_rows(
  # E-traps part
  data_list$etraps_daily %>%
    select(Date, Site, Ambient, Device_type, biomass_mg, abundance) %>%
    filter(!(Device_type == "ID" & Site == "Site3")), # Remove failed device
  
  # AMMOD part
  data_list$ammod_biomass_original %>%
    select(Date, Site, Ambient, Live_mass) %>%
    rename(biomass_mg = Live_mass) %>%
    mutate(Device_type = "AMMOD", abundance = NA) # AMMOD has no daily counts
) %>%
  mutate(
    # Transformations
    log_biomass = log10(biomass_mg + 1),
    log_abundance = log10(abundance + 1),
    # Numeric day (1 to 22)
    day_num = as.numeric(Date - min(Date) + 1)
  ) %>%
  # Clean column names for SAS (lowercase)
  select(date = Date, site = Site, ambient = Ambient, 
         device_type = Device_type, log_biomass, log_abundance, day_num)

# -----------------------------------------------------------------------------
# 3. EXPORT FILE A: BIOMASS VALIDATION (STRICT COMPARISON)
# -----------------------------------------------------------------------------
# Criteria:
# - Only Days 1-12 (Overlap period)
# - Only Sites 1, 2, 3 (Exclude Site 4 because AMMOD4 has no biomass data)
# This ensures a balanced design for the interaction test.

df_biomass_strict <- df_master %>%
  filter(day_num <= 12) %>%
  filter(site != "Site4") %>% 
  filter(!is.na(log_biomass)) 

write_csv(df_biomass_strict, "sas_tables/sas_final_biomass_12day.csv", na = ".")
cat("✅ Generated: sas_tables/sas_final_biomass_12day.csv (Rows:", nrow(df_biomass_strict), ")\n")
cat("   Note: Site 4 excluded to ensure balanced validation design.\n")

# -----------------------------------------------------------------------------
# 4. EXPORT FILE B: ABUNDANCE COMPARISON (FAIRD vs ID - FULL PERIOD)
# -----------------------------------------------------------------------------
# Criteria:
# - All days (1-22)
# - Exclude AMMOD entirely (no daily counts)
# - Keep all sites (Site 4 is valid for FAIRD vs ID)

df_counts_full <- df_master %>%
  filter(device_type != "AMMOD") %>%
  filter(site != "Site3") %>% # ID3 has failed
  filter(!is.na(log_abundance))

write_csv(df_counts_full, "sas_tables/sas_final_counts_22day.csv", na = ".")
cat("✅ Generated: sas_tables/sas_final_counts_22day.csv (Rows:", nrow(df_counts_full), ")\n")