# -----------------------------------------------------------------------------
# Script: FAIRD_ID_AMMOD_prepare_sas_data.R
# Description: Generates two clean, pre-filtered CSV files for SAS.
#              1. Biomass Validation (Days 1-12 only, Sites 1-3 only)
#              2. Abundance Comparison (FAIRD vs ID, all days, balanced sites)
# -----------------------------------------------------------------------------

# Load required packages
library(dplyr)
library(readr)
library(here) # Library for robust path management

# DIAGNOSTIC: Check where 'here' thinks the root is
cat("📂 Project root detected at:", here::here(), "\n")

# 1. Load master data
# 'here::here()' ensures R looks for the file in the project root directory,
# making the script portable across different computers.
load(here::here("data_cleaned.RData"))

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

# Export using 'here' to save in the project root directory
write_csv(df_biomass_strict, here::here("FAIRD_ID_AMMOD_sas_biomass_12day.csv"), na = ".")

cat("✅ Generated:", here::here("FAIRD_ID_AMMOD_sas_biomass_12day.csv"), "\n")
cat("   Rows:", nrow(df_biomass_strict), "\n")
cat("   Note: Site 4 excluded to ensure balanced validation design.\n")

# -----------------------------------------------------------------------------
# 4. EXPORT FILE B: ABUNDANCE COMPARISON (FAIRD vs ID - FULL PERIOD)
# -----------------------------------------------------------------------------
# Criteria:
# - All days (1-22)
# - Exclude AMMOD entirely (no daily counts)
# - Exclude Site 3 entirely (because ID3 failed, avoiding unbalanced comparison)

df_counts_balanced <- df_master %>%
  filter(device_type != "AMMOD") %>%
  filter(site != "Site3") %>%  # CRITICAL: Ensures balanced design (Sites 1, 2, 4 only)
  filter(!is.na(log_abundance))

# Export using 'here' to save in the project root directory
write_csv(df_counts_balanced, here::here("FAIRD_ID_AMMOD_sas_counts_22day.csv"), na = ".")

cat("✅ Generated:", here::here("FAIRD_ID_AMMOD_sas_counts_22day.csv"), "\n")
cat("   Rows:", nrow(df_counts_balanced), "\n")
cat("   Note: Site 3 excluded to ensure balanced paired comparison (ID3 missing).\n")