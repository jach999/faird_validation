#!/usr/bin/env Rscript
# ==============================================================================
# DEVICE COMPARISON SUMMARY TABLE
# Generate comprehensive comparison table across all device types
# ==============================================================================
# NOTE: All data is sourced exclusively through load_data.R and data_list.
#       Do NOT load raw CSVs directly — the Truth Calendar filter (Taxo_Valid)
#       is applied in load_data.R and must be inherited by this script.
# ==============================================================================

# Load required libraries with explicit namespaces
library(readr)
library(dplyr)
library(tidyr)

# ============================================================
# CONSOLE OUTPUT CAPTURE
# ============================================================
script_name <- sub("\\.R$", "", basename(sys.frame(1)$ofile))

output_file <- file.path("outputs",
                         paste0(script_name, "_output_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

sink(output_file, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Script:", paste0(script_name, ".R"), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Source custom colors if needed for output
source("custom_colors.R")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data via load_data.R...\n")

# Single authoritative source — Truth Calendar applied inside load_data.R
source("load_data.R")
load("data_cleaned.RData")

# Extract datasets from data_list
faird_data   <- data_list$faird            # Individual-level FAIRD records
id_data      <- data_list$id               # Individual-level ID records
ammod_biomass   <- data_list$ammod_biomass_for_stats  # Valid biomass measurements only
ammod_taxonomy  <- data_list$ammod_taxonomy           # Presence/absence, Truth Calendar applied (2023 rows)

cat("✓ Data loaded successfully\n\n")

# ==============================================================================
# 2. CALCULATE BIOMASS TOTALS
# ==============================================================================

cat("Calculating total biomass per device type...\n")

# FAIRD total biomass
faird_biomass <- faird_data %>%
  dplyr::summarise(
    # Convert sum from mg to g
    total_biomass_g = sum(Live_mass, na.rm = TRUE) / 1000,
    n_individuals = dplyr::n(),
    .groups = "drop"
  )

# ID total biomass
id_biomass <- id_data %>%
  dplyr::summarise(
    # Convert sum from mg to g
    total_biomass_g = sum(Live_mass, na.rm = TRUE) / 1000,
    n_individuals = dplyr::n(),
    .groups = "drop"
  )

# AMMOD total biomass — Live_mass column is in mg, convert to g
# ammod_biomass_for_stats already excludes NA rows (invalid days)
ammod_biomass_total <- ammod_biomass %>%
  dplyr::summarise(
    total_biomass_g = sum(Live_mass, na.rm = TRUE) / 1000,
    n_days = dplyr::n_distinct(Date),
    .groups = "drop"
  )

cat("✓ Biomass calculated\n")
cat(sprintf("  FAIRD: %.2f g (from %d individuals)\n", 
            faird_biomass$total_biomass_g, faird_biomass$n_individuals))
cat(sprintf("  ID: %.2f g (from %d individuals)\n", 
            id_biomass$total_biomass_g, id_biomass$n_individuals))
cat(sprintf("  AMMOD: %.2f g (across %d days)\n\n", 
            ammod_biomass_total$total_biomass_g, ammod_biomass_total$n_days))

# ==============================================================================
# 3. CALCULATE TAXONOMIC RICHNESS
# ==============================================================================

cat("Calculating taxonomic richness per device type...\n")

# Function to count unique taxa at different levels
count_unique_taxa <- function(data) {
  data %>%
    dplyr::summarise(
      unique_orders = dplyr::n_distinct(Order[!is.na(Order) & Order != ""]),
      unique_families = dplyr::n_distinct(Family[!is.na(Family) & Family != ""]),
      unique_genera = dplyr::n_distinct(Genus[!is.na(Genus) & Genus != ""]),
      .groups = "drop"
    )
}

# FAIRD taxonomic richness
faird_richness <- count_unique_taxa(faird_data)

# ID taxonomic richness
id_richness <- count_unique_taxa(id_data)

# AMMOD taxonomic richness
# Note: AMMOD data is presence/absence at OTU level
ammod_richness <- count_unique_taxa(ammod_taxonomy)

cat("✓ Taxonomic richness calculated\n")
cat(sprintf("  FAIRD: %d Orders | %d Families | %d Genera\n",
            faird_richness$unique_orders,
            faird_richness$unique_families,
            faird_richness$unique_genera))
cat(sprintf("  ID: %d Orders | %d Families | %d Genera\n",
            id_richness$unique_orders,
            id_richness$unique_families,
            id_richness$unique_genera))
cat(sprintf("  AMMOD: %d Orders | %d Families | %d Genera\n\n",
            ammod_richness$unique_orders,
            ammod_richness$unique_families,
            ammod_richness$unique_genera))

# ==============================================================================
# 4. CALCULATE VALID MONITORING DAYS (from Truth Calendar)
# ==============================================================================

cat("Calculating valid monitoring days from Truth Calendar...\n")

# effort_summary comes from data_list (built in load_data.R)
effort_summary <- data_list$effort_summary

# Sum valid days per device type
valid_days_summary <- effort_summary %>%
  dplyr::mutate(Device_type = dplyr::case_when(
    stringr::str_detect(Device, "FAIRD") ~ "FAIRD",
    stringr::str_detect(Device, "^ID")   ~ "ID",
    stringr::str_detect(Device, "AMMOD") ~ "AMMOD"
  )) %>%
  dplyr::group_by(Device_type) %>%
  dplyr::summarise(
    n_devices        = dplyr::n(),
    total_taxo_days  = sum(Valid_Taxo_Days),
    total_bio_days   = sum(Valid_Bio_Days),
    .groups = "drop"
  )

# Convenience vectors for assembly below
faird_valid  <- valid_days_summary %>% dplyr::filter(Device_type == "FAIRD")
id_valid     <- valid_days_summary %>% dplyr::filter(Device_type == "ID")
ammod_valid  <- valid_days_summary %>% dplyr::filter(Device_type == "AMMOD")

cat("✓ Valid days calculated\n")
cat(sprintf("  FAIRD: %d taxo days | %d bio days\n",
            faird_valid$total_taxo_days, faird_valid$total_bio_days))
cat(sprintf("  ID:    %d taxo days | %d bio days\n",
            id_valid$total_taxo_days, id_valid$total_bio_days))
cat(sprintf("  AMMOD: %d taxo days | %d bio days\n\n",
            ammod_valid$total_taxo_days, ammod_valid$total_bio_days))

# ==============================================================================
# 5. CREATE SUMMARY TABLE
# ==============================================================================

cat("Creating summary table...\n")

summary_table <- tibble::tibble(
  Device_Type = c("FAIRD", "Insect Detect", "AMMOD"),
  Nr_Devices = c(
    faird_valid$n_devices,
    id_valid$n_devices,
    ammod_valid$n_devices
  ),
  Est_Live_Biomass_g = c(
    faird_biomass$total_biomass_g,
    id_biomass$total_biomass_g,
    ammod_biomass_total$total_biomass_g
  ),
  N_Individuals = c(
    faird_biomass$n_individuals,
    id_biomass$n_individuals,
    NA  # AMMOD: metabarcoding, no individual counts
  ),
  Total_Valid_Taxo_Days = c(
    faird_valid$total_taxo_days,
    id_valid$total_taxo_days,
    ammod_valid$total_taxo_days
  ),
  Total_Valid_Bio_Days = c(
    faird_valid$total_bio_days,
    id_valid$total_bio_days,
    ammod_valid$total_bio_days
  ),
  Unique_Orders = c(
    faird_richness$unique_orders,
    id_richness$unique_orders,
    ammod_richness$unique_orders
  ),
  Unique_Families = c(
    faird_richness$unique_families,
    id_richness$unique_families,
    ammod_richness$unique_families
  ),
  Unique_Genera = c(
    faird_richness$unique_genera,
    id_richness$unique_genera,
    ammod_richness$unique_genera
  )
)

# ==============================================================================
# 6. DISPLAY RESULTS
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n")
cat("DEVICE COMPARISON SUMMARY TABLE\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n\n")

# Create formatted display table
display_table <- summary_table %>%
  dplyr::mutate(
    Device                    = Device_Type,
    `Nr. Devices`             = Nr_Devices,
    `Total Biomass (g)`       = sprintf("%.2f", Est_Live_Biomass_g),
    `Total Individuals`       = ifelse(is.na(N_Individuals), "N/A (DNA)",
                                       sprintf("%d", as.integer(N_Individuals))),
    `Total Valid Taxo Days (n)` = Total_Valid_Taxo_Days,
    `Total Valid Bio Days (n)`  = Total_Valid_Bio_Days,
    `Unique Orders`           = Unique_Orders,
    `Unique Families`         = Unique_Families,
    `Unique Genera`           = Unique_Genera
  ) %>%
  dplyr::select(Device, `Nr. Devices`, `Total Biomass (g)`, `Total Individuals`,
                `Total Valid Taxo Days (n)`, `Total Valid Bio Days (n)`,
                `Unique Orders`, `Unique Families`, `Unique Genera`)

print(display_table, width = Inf)

cat("\n")

# ==============================================================================
# 7. SAVE OUTPUT
# ==============================================================================

# Save raw data
output_file <- "outputs/Device_Summary_Comparison.csv"
readr::write_csv(summary_table, output_file)
cat(sprintf("✓ Raw data saved to: %s\n", output_file))

# Save formatted version for presentation
formatted_file <- "outputs/Device_Summary_Formatted.csv"
readr::write_csv(display_table, formatted_file)
cat(sprintf("✓ Formatted table saved to: %s\n\n", formatted_file))

# ==============================================================================
# 7. ADDITIONAL INSIGHTS
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("KEY INSIGHTS\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Biomass per valid bio day
cat("Average Biomass per Valid Bio Day:\n")
cat(sprintf("  FAIRD: %.3f g/day\n",
            summary_table$Est_Live_Biomass_g[1] / summary_table$Total_Valid_Bio_Days[1]))
cat(sprintf("  ID: %.3f g/day\n",
            summary_table$Est_Live_Biomass_g[2] / summary_table$Total_Valid_Bio_Days[2]))
cat(sprintf("  AMMOD: %.3f g/day\n\n",
            summary_table$Est_Live_Biomass_g[3] / summary_table$Total_Valid_Bio_Days[3]))

# Taxonomic efficiency (genera per gram)
cat("Taxonomic Efficiency (Unique Genera per gram of biomass):\n")
cat(sprintf("  FAIRD: %.2f genera/g\n", 
            summary_table$Unique_Genera[1] / summary_table$Est_Live_Biomass_g[1]))
cat(sprintf("  ID: %.2f genera/g\n", 
            summary_table$Unique_Genera[2] / summary_table$Est_Live_Biomass_g[2]))
cat(sprintf("  AMMOD: %.2f genera/g\n\n", 
            summary_table$Unique_Genera[3] / summary_table$Est_Live_Biomass_g[3]))

# Individuals per genus (for e-traps only)
cat("Sampling Intensity (Individuals per unique genus):\n")
cat(sprintf("  FAIRD: %.1f individuals/genus\n", 
            summary_table$N_Individuals[1] / summary_table$Unique_Genera[1]))
cat(sprintf("  ID: %.1f individuals/genus\n", 
            summary_table$N_Individuals[2] / summary_table$Unique_Genera[2]))
cat("  AMMOD: N/A (DNA metabarcoding, no individual counts)\n\n")

# Diversity comparison
cat("Taxonomic Diversity Comparison:\n")
cat(sprintf("  FAIRD captures %.1f%% of AMMOD's genus-level diversity\n",
            100 * summary_table$Unique_Genera[1] / summary_table$Unique_Genera[3]))
cat(sprintf("  ID captures %.1f%% of AMMOD's genus-level diversity\n",
            100 * summary_table$Unique_Genera[2] / summary_table$Unique_Genera[3]))
cat(sprintf("  FAIRD is %.1f× more diverse than ID at genus level\n\n",
            summary_table$Unique_Genera[1] / summary_table$Unique_Genera[2]))

cat("═══════════════════════════════════════════════════════════════════\n")
cat("SCRIPT COMPLETED SUCCESSFULLY\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")