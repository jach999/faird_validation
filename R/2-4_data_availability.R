# ==============================================================================
# 2-4_data_availability.R
# Manuscript Section 2.4 (Data availability and exclusions): Table 1
# ------------------------------------------------------------------------------
# Counts only, no models.
#
# 1. Incident log: every data loss of the campaign, with unit, data stream,
#    period, cause and treatment (field and laboratory records). Unit-days
#    per incident are counted from the calendar. Incidents can overlap (for
#    example AMMOD4 biomass), so their unit-days do not add up to the loss.
# 2. Valid unit-days per system and data stream, from the validity calendar
#    `effort_base` built in 00_load_data.R (Field_Valid, Taxo_Valid,
#    Bio_Valid, Field_Valid_Operational_ID). These are the denominators of
#    Table 2 (3-2a_composition_orders.R) and the sample sizes of Section 3.
# 3. Checks: the incident log must reproduce every FALSE cell of the validity
#    calendar, and the calendar must agree with the data actually present
#    (AMMOD biomass rows, AMMOD taxonomy records, camera records).
#
# The ID power outages are days of partial, unknown recording effort, not
# losses: they stay valid in Field_Valid and are removed only in the runs
# with the window excluded (Field_Valid_Operational_ID).
# ==============================================================================

required_packages <- c("tidyverse", "here")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)

source(here::here("R", "00_load_data.R"))

SCRIPT  <- "2-4_data_availability"
TAB_DIR <- here::here("outputs", "tables")
CON_DIR <- here::here("outputs", "console")
for (d in c(TAB_DIR, CON_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, paste0(SCRIPT, "_output.txt")), split = TRUE)

cat("==============================================================================\n")
cat("SECTION 2.4, DATA AVAILABILITY AND EXCLUSIONS\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cal <- data_list$effort_base %>%
  dplyr::mutate(Date = as.Date(Date),
                System = stringr::str_remove(Device, "[0-9]+$"))
n_days_campaign <- dplyr::n_distinct(cal$Date)
units_per_system <- cal %>% dplyr::distinct(System, Device) %>% dplyr::count(System, name = "n_units")


# ==============================================================================
# SECTION 1: INCIDENT LOG (Table 1)
# ==============================================================================

d <- function(x) as.Date(x)
incidents <- tibble::tribble(
  ~System, ~Units,                       ~Stream,    ~Start,          ~End,            ~Cause_type,  ~Cause,                                                            ~Treatment,
  "AMMOD", "AMMOD1,AMMOD2,AMMOD3,AMMOD4", "all",      d("2023-09-04"), d("2023-09-04"), "Protocol",   "Bottle change, no sample collected",                              "No sample",
  "AMMOD", "AMMOD1,AMMOD2,AMMOD3,AMMOD4", "biomass",  d("2023-09-05"), d("2023-09-13"), "Laboratory", "ATL buffer extraction destroys the material for weighing",        "Lost",
  "AMMOD", "AMMOD3",                      "biomass",  d("2023-08-31"), d("2023-09-03"), "Laboratory", "ATL buffer extraction",                                           "Lost",
  "AMMOD", "AMMOD4",                      "biomass",  d("2023-08-23"), d("2023-09-03"), "Laboratory", "ATL buffer extraction",                                           "Lost",
  "AMMOD", "AMMOD4",                      "all",      d("2023-08-28"), d("2023-09-03"), "Device",     "Bottle rotation jammed, seven days in one bottle",                "Excluded",
  "AMMOD", "AMMOD4",                      "taxonomy", d("2023-09-10"), d("2023-09-10"), "Laboratory", "Sample discarded",                                                "Excluded",
  "ID",    "ID3",                         "all",      d("2023-08-23"), d("2023-09-13"), "Operation",  "Corrupted SD card, detected after the campaign",                  "Lost",
  "ID",    "ID1,ID2,ID4",                 "all",      d("2023-08-28"), d("2023-09-01"), "Power",      "Insufficient solar charging during a cold front",                 "Partial operation; analysed with and without"
)

# Expand to one row per unit-day and incident.
incident_days <- incidents %>%
  dplyr::mutate(incident = dplyr::row_number()) %>%
  tidyr::separate_rows(Units, sep = ",") %>%
  dplyr::rename(Device = Units) %>%
  dplyr::mutate(Date = purrr::map2(Start, End, ~ seq(.x, .y, by = "day"))) %>%
  tidyr::unnest(Date)

table1 <- incident_days %>%
  dplyr::group_by(incident, System, Stream, Start, End, Cause_type, Cause, Treatment) %>%
  dplyr::summarise(Units = paste(sort(unique(Device)), collapse = ", "),
                   Unit_days = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(incident) %>%
  dplyr::select(System, Units, Stream, Start, End, Unit_days, Cause_type, Cause, Treatment)

cat("=== TABLE 1: INCIDENTS ===\n")
cat(sprintf("Planned unit-days per system: %d units x %d days\n",
            max(units_per_system$n_units), n_days_campaign))
print(table1, n = Inf, width = Inf)


# ==============================================================================
# SECTION 2: VALID UNIT-DAYS PER SYSTEM AND DATA STREAM
# ==============================================================================

valid_by_system <- cal %>%
  dplyr::group_by(System) %>%
  dplyr::summarise(
    Planned            = dplyr::n(),
    Field_valid        = sum(Field_Valid),
    Taxonomy_valid     = sum(Taxo_Valid),
    Biomass_valid      = sum(Bio_Valid),
    Full_effort_valid  = sum(Field_Valid_Operational_ID),
    .groups = "drop")

valid_by_unit <- cal %>%
  dplyr::group_by(System, Device) %>%
  dplyr::summarise(
    Field_valid       = sum(Field_Valid),
    Taxonomy_valid    = sum(Taxo_Valid),
    Biomass_valid     = sum(Bio_Valid),
    Full_effort_valid = sum(Field_Valid_Operational_ID),
    .groups = "drop")

cat("\n=== VALID UNIT-DAYS PER SYSTEM ===\n")
cat("Biomass_valid for AMMOD: days with a weighed sample; for the camera systems:\n")
cat("  days with valid recordings (biomass estimated from detections).\n")
cat("Full_effort_valid: Field_valid minus the ID power-outage window.\n\n")
print(valid_by_system, n = Inf, width = Inf)
cat("\n=== VALID UNIT-DAYS PER UNIT ===\n")
print(valid_by_unit, n = Inf, width = Inf)


# ==============================================================================
# SECTION 3: CHECKS
# ==============================================================================

cat("\n=== CHECKS ===\n")

# 3a. Incident log vs validity calendar.
lost <- function(streams, exclude_types = "Power") {
  incident_days %>%
    dplyr::filter(Stream %in% streams, !(Cause_type %in% exclude_types)) %>%
    dplyr::distinct(Device, Date) %>%
    dplyr::mutate(lost = TRUE)
}
cmp_flag <- function(flag, streams, exclude_types = "Power") {
  cal %>%
    dplyr::select(Device, Date, valid = dplyr::all_of(flag)) %>%
    dplyr::left_join(lost(streams, exclude_types), by = c("Device", "Date")) %>%
    dplyr::mutate(lost = dplyr::coalesce(lost, FALSE)) %>%
    dplyr::filter(valid == lost)
}
checks <- list(
  Field_Valid                = cmp_flag("Field_Valid", "all"),
  Taxo_Valid                 = cmp_flag("Taxo_Valid", c("all", "taxonomy")),
  Bio_Valid                  = cmp_flag("Bio_Valid", c("all", "biomass")),
  Field_Valid_Operational_ID = cmp_flag("Field_Valid_Operational_ID", "all", character(0))
)
for (nm in names(checks)) {
  n_bad <- nrow(checks[[nm]])
  cat(sprintf("%s %s: %d unit-days where the incident log and the calendar disagree\n",
              ifelse(n_bad == 0, "[OK]  ", "[WARN]"), nm, n_bad))
  if (n_bad > 0) print(checks[[nm]], n = Inf)
}

# 3b. Calendar vs data present.
ammod_bio_days <- data_list$ammod_biomass_for_stats %>%
  dplyr::distinct(Device, Date = as.Date(Date))
bio_cal_ammod <- cal %>% dplyr::filter(System == "AMMOD", Bio_Valid) %>%
  dplyr::distinct(Device, Date)
cat(sprintf("%s AMMOD biomass: %d weighed unit-days in the data, %d valid in the calendar, %d in both\n",
            ifelse(nrow(dplyr::setdiff(ammod_bio_days, bio_cal_ammod)) == 0 &&
                     nrow(dplyr::setdiff(bio_cal_ammod, ammod_bio_days)) == 0, "[OK]  ", "[WARN]"),
            nrow(ammod_bio_days), nrow(bio_cal_ammod),
            nrow(dplyr::intersect(ammod_bio_days, bio_cal_ammod))))

ammod_taxo_days <- data_list$ammod_taxonomy %>%
  dplyr::distinct(Device, Date = as.Date(Date))
taxo_cal_ammod <- cal %>% dplyr::filter(System == "AMMOD", Taxo_Valid) %>%
  dplyr::distinct(Device, Date)
cat(sprintf("%s AMMOD taxonomy: %d unit-days with at least one record, %d valid in the calendar; records outside valid days: %d\n",
            ifelse(nrow(dplyr::setdiff(ammod_taxo_days, taxo_cal_ammod)) == 0, "[OK]  ", "[WARN]"),
            nrow(ammod_taxo_days), nrow(taxo_cal_ammod),
            nrow(dplyr::setdiff(ammod_taxo_days, taxo_cal_ammod))))

cam_days <- data_list$etraps_daily %>% dplyr::distinct(Device, Date = as.Date(Date))
cam_cal  <- cal %>% dplyr::filter(System %in% c("FAIRD", "ID"), Field_Valid) %>%
  dplyr::distinct(Device, Date)
cat(sprintf("%s Camera systems: %d unit-days in the daily series, %d valid in the calendar\n",
            ifelse(nrow(dplyr::setdiff(cam_days, cam_cal)) == 0 &&
                     nrow(dplyr::setdiff(cam_cal, cam_days)) == 0, "[OK]  ", "[WARN]"),
            nrow(cam_days), nrow(cam_cal)))


# ==============================================================================
# SECTION 4: EXPORT
# ==============================================================================

readr::write_csv(table1, file.path(TAB_DIR, paste0(SCRIPT, "_table1_incidents.csv")))
readr::write_csv(valid_by_system, file.path(TAB_DIR, paste0(SCRIPT, "_valid_by_system.csv")))
readr::write_csv(valid_by_unit, file.path(TAB_DIR, paste0(SCRIPT, "_valid_by_unit.csv")))
cat("\n[OK] Tables written to", TAB_DIR, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

while (sink.number() > 0) sink()
