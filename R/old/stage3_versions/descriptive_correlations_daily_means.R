# ==============================================================================
# descriptive_correlations_daily_means.R
# Simple correlations of daily system means, biomass and abundance (Stage 3)
# ------------------------------------------------------------------------------
# Question: what does a simple correlation give for each comparison of the
# V6 manuscript, and does camera biomass carry more noise than camera
# abundance? Descriptive only: no models, no tests, no p-values.
#
# Place in V6: Results 3.1 (values in parentheses next to the model-based
# results), Methods 2.5 (definition). At the repository reorganisation this
# script becomes a section of the quantity / concordance file.
#
# Definition (identical to TOTAL in Scherber_temporal_concordance_biomass_4.R):
#   daily mean per unit, averaged within block, then across blocks; one
#   series per unit; correlation between the two series over the dates on
#   which both units have a value in the same block (block-day pairing).
#   "Unit" is a system, a half of the FAIRD network, or a single device in
#   the Maize-only references.
#
# Measures:
#   biomass   (mg live mass): all three systems.
#   abundance (individuals or ID tracks per day): FAIRD and ID only; AMMOD
#             taxonomy is presence/absence and has no abundance.
# Coefficients: Pearson on the raw scale (as TOTAL) and Spearman (rank
# based, not dominated by peak days; comparable with earlier field studies).
#
# Comparisons, windows and exclusions: copied from the concordance script.
#   A_12d      = AMMOD biomass window (Aug 23 to Sep 3), full ID series.
#   B_22d_full = Aug 23 to Sep 13, all days.
#   B_22d_excl = Aug 23 to Sep 13, ID power-outage window (Aug 28 to Sep 1)
#                removed for BOTH units of every comparison.
#
# Check: biomass Pearson must reproduce r_total of
# outputs/tables/Scherber_concordance_biomass.csv (spline df 3 rows). The
# script compares them and reports the largest difference.
# ==============================================================================

library(tidyverse)
library(here)

source(here::here("R", "00_load_data.R"))

dir.create(here::here("outputs", "console"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("outputs", "tables"),  showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(here::here("outputs", "console", "descriptive_correlations_daily_means_output.txt"),
     split = TRUE)

cat("==============================================================================\n")
cat("SIMPLE CORRELATIONS OF DAILY SYSTEM MEANS, BIOMASS AND ABUNDANCE\n")
cat("Descriptive only: no models, no tests\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: DATA
# ==============================================================================

etraps_daily <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::transmute(Date, Device, Device_type, Site,
                   biomass = biomass_mg, abundance = as.numeric(abundance))

ammod_daily <- data_list$ammod_biomass_for_stats %>%
  dplyr::transmute(Date, Device, Device_type = "AMMOD", Site,
                   biomass = Live_mass, abundance = NA_real_)

daily_all <- dplyr::bind_rows(etraps_daily, ammod_daily) %>%
  dplyr::mutate(
    Date   = as.Date(Date),
    Device = as.character(Device),
    Site   = as.character(Site),
    Block  = ifelse(Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")
  )

dates_12 <- seq(min(ammod_daily$Date), max(ammod_daily$Date), by = "day")
dates_22 <- seq(min(etraps_daily$Date), max(etraps_daily$Date), by = "day")

adverse_dates <- data_list$effort_base %>%
  dplyr::filter(Field_Valid, !Field_Valid_Operational_ID) %>%
  dplyr::distinct(Date) %>%
  dplyr::pull(Date) %>%
  as.Date()

windows <- list(
  A_12d      = list(dates = dates_12, drop = NULL),
  B_22d_full = list(dates = dates_22, drop = NULL),
  B_22d_excl = list(dates = dates_22, drop = adverse_dates)
)

cat(sprintf("A_12d: %s to %s (%d dates)\n", min(dates_12), max(dates_12), length(dates_12)))
cat(sprintf("B_22d: %s to %s (%d dates)\n", min(dates_22), max(dates_22), length(dates_22)))
cat(sprintf("Excluded window (B_22d_excl, both units): %s\n\n",
            paste(format(sort(adverse_dates)), collapse = ", ")))


# ==============================================================================
# SECTION 2: COMPARISONS (same set and ids as the concordance script)
# ==============================================================================

sys_dev <- list(AMMOD = paste0("AMMOD", 1:4), FAIRD = paste0("FAIRD", 1:4),
                ID = paste0("ID", 1:4))

cmp_between <- function(a, b, window, version) {
  list(id = sprintf("%s-%s|%s|%s", a, b, window, version), type = "between",
       comparison = sprintf("%s vs %s", a, b), systems = c(a, b),
       lab_A = a, dev_A = sys_dev[[a]], lab_B = b, dev_B = sys_dev[[b]],
       window = window, version = version)
}

cmp_internal <- function(system, split, lab_A, dev_A, lab_B, dev_B, window, version) {
  list(id = sprintf("%s|%s|%s|%s", system, split, window, version), type = "internal",
       comparison = sprintf("%s: %s vs %s", system,
                            paste(dev_A, collapse = "+"), paste(dev_B, collapse = "+")),
       systems = system, lab_A = lab_A, dev_A = dev_A, lab_B = lab_B, dev_B = dev_B,
       window = window, version = version)
}

comparisons <- list()
for (p in list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))) {
  for (v in c("two_block", "maize_only")) {
    comparisons <- c(comparisons, list(cmp_between(p[1], p[2], "A_12d", v)))
  }
}
for (w in c("B_22d_full", "B_22d_excl")) {
  for (v in c("two_block", "maize_only")) {
    comparisons <- c(comparisons, list(cmp_between("FAIRD", "ID", w, v)))
  }
}
for (w in c("A_12d", "B_22d_full", "B_22d_excl")) {
  comparisons <- c(comparisons, list(
    cmp_internal("FAIRD", "S1", "FAIRD_h13", c("FAIRD1", "FAIRD3"),
                 "FAIRD_h24", c("FAIRD2", "FAIRD4"), w, "two_block"),
    cmp_internal("FAIRD", "S2", "FAIRD_h14", c("FAIRD1", "FAIRD4"),
                 "FAIRD_h23", c("FAIRD2", "FAIRD3"), w, "two_block")))
  maize_systems <- if (w == "A_12d") c("FAIRD", "AMMOD", "ID") else c("FAIRD", "ID")
  for (s in maize_systems) {
    comparisons <- c(comparisons, list(
      cmp_internal(s, "maize", paste0(s, "1"), paste0(s, "1"),
                   paste0(s, "2"), paste0(s, "2"), w, "maize_only")))
  }
}

cat(sprintf("%d comparisons defined.\n\n", length(comparisons)))


# ==============================================================================
# SECTION 3: CORRELATION OF DAILY UNIT MEANS
# ==============================================================================

# daily_mean_r(): one comparison, one measure. Returns n dates, Pearson and
# Spearman r between the two unit series, or NULL when the measure does not
# exist for one of the systems (abundance with AMMOD).
daily_mean_r <- function(cmp, measure) {
  if (measure == "abundance" && "AMMOD" %in% cmp$systems) return(NULL)
  w   <- windows[[cmp$window]]
  map <- dplyr::bind_rows(
    tibble::tibble(Device = cmp$dev_A, Unit = cmp$lab_A),
    tibble::tibble(Device = cmp$dev_B, Unit = cmp$lab_B)
  )
  d <- daily_all %>%
    dplyr::inner_join(map, by = "Device") %>%
    dplyr::rename(value = dplyr::all_of(measure)) %>%
    dplyr::filter(!is.na(value), Date %in% w$dates)
  if (!is.null(w$drop)) d <- d %>% dplyr::filter(!(Date %in% w$drop))
  if (cmp$version == "maize_only") d <- d %>% dplyr::filter(Block == "BLOCK1")
  means <- d %>%
    dplyr::group_by(Block, Date) %>%
    dplyr::filter(dplyr::n_distinct(Unit) == 2) %>%
    dplyr::group_by(Unit, Block, Date) %>%
    dplyr::summarise(value = mean(value), .groups = "drop") %>%
    dplyr::group_by(Unit, Date) %>%
    dplyr::summarise(value = mean(value), .groups = "drop")
  wide <- means %>%
    tidyr::pivot_wider(names_from = Unit, values_from = value) %>%
    tidyr::drop_na()
  a <- wide[[cmp$lab_A]]
  b <- wide[[cmp$lab_B]]
  tibble::tibble(
    id = cmp$id, type = cmp$type, comparison = cmp$comparison,
    window = cmp$window, version = cmp$version, measure = measure,
    n_dates  = nrow(wide),
    pearson  = stats::cor(a, b, method = "pearson"),
    spearman = stats::cor(a, b, method = "spearman")
  )
}

res <- dplyr::bind_rows(lapply(comparisons, function(cmp) {
  dplyr::bind_rows(daily_mean_r(cmp, "biomass"), daily_mean_r(cmp, "abundance"))
}))

# FAIRD halves: mean of the two splits, as in the concordance script summary.
halves <- res %>%
  dplyr::filter(type == "internal", version == "two_block") %>%
  dplyr::group_by(window, measure) %>%
  dplyr::summarise(pearson = mean(pearson), spearman = mean(spearman),
                   n_dates = min(n_dates), n_splits = dplyr::n(), .groups = "drop")


# ==============================================================================
# SECTION 4: REPRODUCTION CHECK AGAINST THE CONCORDANCE SCRIPT
# ==============================================================================

cat("=== CHECK: biomass Pearson vs r_total of the concordance script ===\n")
ref_file <- here::here("outputs", "tables", "Scherber_concordance_biomass.csv")
if (file.exists(ref_file)) {
  ref <- readr::read_csv(ref_file, show_col_types = FALSE) %>%
    dplyr::filter(spline_df == 3) %>%
    dplyr::select(id, r_total)
  chk <- res %>%
    dplyr::filter(measure == "biomass") %>%
    dplyr::inner_join(ref, by = "id") %>%
    dplyr::mutate(diff = abs(pearson - r_total))
  n_missing <- sum(res$measure == "biomass") - nrow(chk)
  max_diff  <- max(chk$diff)
  cat(sprintf("Matched %d of %d biomass comparisons; largest difference %.2e.\n",
              nrow(chk), sum(res$measure == "biomass"), max_diff))
  if (n_missing == 0 && max_diff < 1e-6) {
    cat("[OK] Biomass Pearson reproduces r_total for every comparison.\n\n")
  } else {
    cat("[WARN] Biomass Pearson does not reproduce r_total; rows with a difference:\n")
    print(chk %>% dplyr::filter(diff >= 1e-6) %>%
            dplyr::select(id, pearson, r_total, diff), n = Inf)
    cat("\n")
  }
} else {
  cat("[WARN] Reference table not found:", ref_file, "\n\n")
}


# ==============================================================================
# SECTION 5: RESULTS
# ==============================================================================

fmt <- function(x) x %>% dplyr::mutate(dplyr::across(c(pearson, spearman), ~ round(.x, 3)))

side_by_side <- res %>%
  tidyr::pivot_wider(id_cols = c(id, type, comparison, window, version),
                     names_from = measure, values_from = c(n_dates, pearson, spearman))

cat("=== BETWEEN SYSTEMS (two blocks = main text; Maize only = secondary) ===\n")
print(side_by_side %>% dplyr::filter(type == "between") %>%
        dplyr::select(-id, -type) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with(c("pearson", "spearman")),
                                    ~ round(.x, 3))),
      n = Inf, width = Inf)

cat("\n=== WITHIN SYSTEMS ===\n")
print(side_by_side %>% dplyr::filter(type == "internal") %>%
        dplyr::select(-id, -type) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with(c("pearson", "spearman")),
                                    ~ round(.x, 3))),
      n = Inf, width = Inf)

cat("\n=== FAIRD HALVES OVER BOTH BLOCKS, MEAN OF THE TWO SPLITS ===\n")
print(fmt(halves), n = Inf, width = Inf)

cat("\nReading notes (definitions, not results):\n")
cat("  - Pearson on the raw scale is dominated by days of high values (the seasonal\n")
cat("    peak); Spearman weights all days equally.\n")
cat("  - Biomass vs abundance within the camera systems: a lower biomass r than\n")
cat("    abundance r indicates noise added by the mass estimates.\n")


# ==============================================================================
# SECTION 6: EXPORT
# ==============================================================================

out_file <- here::here("outputs", "tables", "descriptive_correlations_daily_means.csv")
readr::write_csv(res, out_file)
cat("\n[OK] Table written:", out_file, "\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

while (sink.number() > 0) sink()
