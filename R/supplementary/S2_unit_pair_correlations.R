# ==============================================================================
# supplementary/S2_unit_pair_correlations.R
# Supplementary Material S2: conventional correlations between paired units
# ------------------------------------------------------------------------------
# Correlations between two individual units on log10(x + 1)-transformed daily
# values, as used in earlier field comparisons of monitoring devices. Given
# for comparison with that literature only: the manuscript bases its
# inference on the system-level models of Section 3.1, and single units
# cannot separate a system from the peculiarities of a device.
#
# Pairs:
#   Within systems (two units of the same habitat): FAIRD1 vs FAIRD2,
#     FAIRD3 vs FAIRD4, ID1 vs ID2 (22 days; and without the ID power-outage
#     window, 17 days) in abundance and biomass; AMMOD1 vs AMMOD2 in biomass
#     (12 days).
#   Between FAIRD and AMMOD (co-located units, biomass): Site1 and Site2
#     (12 days), Site3 (8 days); site pairs stacked for Maize and for all
#     three sites.
# Statistics per pair: n, Pearson r and R2 on log10(x + 1), OLS slope of the
# second unit on the first (with 95% CI), Spearman rho (rank based, so the
# same on the raw and the log scale).
# Other comparisons of this kind (Kendall tau on daily changes, robust
# regression, cross-correlation) are in the frozen V5 state of the
# repository (tag v5-supplementary).
# ==============================================================================

required_packages <- c("tidyverse", "here")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)

source(here::here("R", "00_load_data.R"))

SCRIPT  <- "S2_unit_pair_correlations"
TAB_DIR <- here::here("outputs", "tables")
CON_DIR <- here::here("outputs", "console")
for (dd in c(TAB_DIR, CON_DIR)) dir.create(dd, showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, paste0(SCRIPT, "_output.txt")), split = TRUE)

cat("==============================================================================\n")
cat("SUPPLEMENTARY S2: CORRELATIONS BETWEEN PAIRED UNITS, log10(x + 1)\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: DAILY VALUES PER UNIT
# ==============================================================================

daily <- dplyr::bind_rows(
  data_list$etraps_daily %>%
    dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
    dplyr::transmute(Date = as.Date(Date), Device, Site,
                     abundance = as.numeric(abundance), biomass = biomass_mg),
  data_list$ammod_biomass_for_stats %>%
    dplyr::transmute(Date = as.Date(Date), Device, Site,
                     abundance = NA_real_, biomass = Live_mass)
)

outage_dates <- data_list$effort_base %>%
  dplyr::filter(Field_Valid, !Field_Valid_Operational_ID) %>%
  dplyr::distinct(Date) %>%
  dplyr::pull(Date) %>%
  as.Date()
cat(sprintf("ID power-outage window (%d dates): %s\n\n", length(outage_dates),
            paste(format(sort(outage_dates)), collapse = ", ")))


# ==============================================================================
# SECTION 2: PAIR STATISTICS
# ==============================================================================

# pair_values(): daily values of two units on their common dates.
pair_values <- function(dev_x, dev_y, measure, drop_dates = NULL) {
  x <- daily %>% dplyr::filter(Device == dev_x) %>%
    dplyr::select(Date, x = dplyr::all_of(measure))
  y <- daily %>% dplyr::filter(Device == dev_y) %>%
    dplyr::select(Date, y = dplyr::all_of(measure))
  out <- dplyr::inner_join(x, y, by = "Date") %>% tidyr::drop_na()
  if (!is.null(drop_dates)) out <- out %>% dplyr::filter(!(Date %in% drop_dates))
  out
}

# pair_stats(): statistics on log10(x + 1); x and y may stack several pairs.
pair_stats <- function(v, label, group, measure, window) {
  lx  <- log10(v$x + 1)
  ly  <- log10(v$y + 1)
  fit <- stats::lm(ly ~ lx)
  ci  <- stats::confint(fit)["lx", ]
  tibble::tibble(
    group = group, pair = label, measure = measure, window = window,
    n = nrow(v),
    pearson_log = stats::cor(lx, ly),
    R2_log      = summary(fit)$r.squared,
    p_log       = summary(fit)$coefficients["lx", "Pr(>|t|)"],
    slope       = unname(stats::coef(fit)["lx"]),
    slope_lo    = unname(ci[1]),
    slope_hi    = unname(ci[2]),
    spearman    = stats::cor(v$x, v$y, method = "spearman")
  )
}

rows <- list()

# Within systems.
within_pairs <- list(c("FAIRD1", "FAIRD2"), c("FAIRD3", "FAIRD4"), c("ID1", "ID2"))
for (pp in within_pairs) {
  for (m in c("abundance", "biomass")) {
    rows[[length(rows) + 1]] <- pair_stats(pair_values(pp[1], pp[2], m),
                                           paste(pp, collapse = " vs "), "within", m, "22 days")
    if (startsWith(pp[1], "ID")) {
      rows[[length(rows) + 1]] <- pair_stats(pair_values(pp[1], pp[2], m, outage_dates),
                                             paste(pp, collapse = " vs "), "within", m,
                                             "22 days, outage window excluded")
    }
  }
}
rows[[length(rows) + 1]] <- pair_stats(pair_values("AMMOD1", "AMMOD2", "biomass"),
                                       "AMMOD1 vs AMMOD2", "within", "biomass", "12 days")

# Between FAIRD and the co-located AMMOD (AMMOD on x, FAIRD on y).
site_pairs <- list(Site1 = c("AMMOD1", "FAIRD1"), Site2 = c("AMMOD2", "FAIRD2"),
                   Site3 = c("AMMOD3", "FAIRD3"))
site_values <- lapply(site_pairs, function(pp) pair_values(pp[1], pp[2], "biomass"))
for (s in names(site_pairs)) {
  rows[[length(rows) + 1]] <- pair_stats(site_values[[s]],
                                         paste(rev(site_pairs[[s]]), collapse = " vs "),
                                         "between", "biomass", sprintf("%s, AMMOD days", s))
}
rows[[length(rows) + 1]] <- pair_stats(dplyr::bind_rows(site_values[c("Site1", "Site2")]),
                                       "FAIRD vs AMMOD, Site1 and Site2 stacked",
                                       "between", "biomass", "Maize, AMMOD days")
rows[[length(rows) + 1]] <- pair_stats(dplyr::bind_rows(site_values),
                                       "FAIRD vs AMMOD, Sites 1 to 3 stacked",
                                       "between", "biomass", "All sites, AMMOD days")

res <- dplyr::bind_rows(rows)

cat("=== CORRELATIONS BETWEEN PAIRED UNITS ===\n")
cat("x = first unit of the pair (AMMOD for the between-system rows), y = second.\n\n")
print(res %>% dplyr::mutate(dplyr::across(tidyselect::where(is.numeric) & !dplyr::all_of("n"),
                                          ~ round(.x, 3))),
      n = Inf, width = Inf)
cat("\n  NOTES:\n")
cat("  - Stacked rows combine site pairs; their observations are not independent.\n")
cat("  - Spearman rho is rank based and does not depend on the transformation.\n")


# ==============================================================================
# SECTION 3: EXPORT
# ==============================================================================

readr::write_csv(res, file.path(TAB_DIR, paste0(SCRIPT, ".csv")))
cat("\n[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, ".csv")), "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

while (sink.number() > 0) sink()
