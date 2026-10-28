# ==============================================================================
# supplementary/S4_abundance_cameras.R
# Supplementary Material S4: daily abundance of the camera systems
# ------------------------------------------------------------------------------
# Descriptive only, no models. Abundance (individuals per day; ID: tracks per
# day) exists for FAIRD and ID only: AMMOD taxonomy is presence/absence.
#
# Table: mean abundance per unit-day by system, over all 22 days and without
# the ID power-outage window (removed for both systems), and by unit.
# Figure: daily abundance of each unit and daily system means (within block,
# then across blocks), linear axis; outage-window days drawn with open
# symbols (days of partial, unknown ID recording effort).
# Correlations of daily abundance between and within the camera systems are
# in 3-1b_quantity_concordance.R (simple correlations of daily means).
# ==============================================================================

required_packages <- c("tidyverse", "here", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)

source(here::here("R", "00_load_data.R"))
source(here::here("R", "00_custom_colors.R"))

SCRIPT  <- "S4_abundance_cameras"
FIG_DIR <- here::here("outputs", "figures")
TAB_DIR <- here::here("outputs", "tables")
CON_DIR <- here::here("outputs", "console")
for (d in c(FIG_DIR, TAB_DIR, CON_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, paste0(SCRIPT, "_output.txt")), split = TRUE)

cat("==============================================================================\n")
cat("SUPPLEMENTARY S4: DAILY ABUNDANCE OF THE CAMERA SYSTEMS\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: DATA
# ==============================================================================

outage_dates <- data_list$effort_base %>%
  dplyr::filter(Field_Valid, !Field_Valid_Operational_ID) %>%
  dplyr::distinct(Date) %>%
  dplyr::pull(Date) %>%
  as.Date()

abund <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::transmute(Date = as.Date(Date), Device, Device_type, Site,
                   abundance = as.numeric(abundance),
                   Block  = ifelse(Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2"),
                   Effort = ifelse(Date %in% outage_dates, "outage window", "full effort"))

cat(sprintf("%d unit-days (%s); outage window: %s\n\n", nrow(abund),
            paste(sort(unique(abund$Device)), collapse = ", "),
            paste(format(sort(outage_dates)), collapse = ", ")))


# ==============================================================================
# SECTION 2: MEAN ABUNDANCE PER UNIT-DAY
# ==============================================================================

by_system <- dplyr::bind_rows(
  abund %>% dplyr::mutate(Days = "22 days"),
  abund %>% dplyr::filter(Effort == "full effort") %>%
    dplyr::mutate(Days = "outage window excluded")
) %>%
  dplyr::group_by(Days, Device_type) %>%
  dplyr::summarise(n_units = dplyr::n_distinct(Device), n_unit_days = dplyr::n(),
                   total = sum(abundance), mean_per_unit_day = mean(abundance),
                   median_per_unit_day = stats::median(abundance), .groups = "drop")

by_unit <- abund %>%
  dplyr::group_by(Device_type, Device, Site) %>%
  dplyr::summarise(n_days = dplyr::n(), total = sum(abundance),
                   mean_per_day = mean(abundance), max_per_day = max(abundance),
                   zero_days = sum(abundance == 0), .groups = "drop")

cat("=== MEAN ABUNDANCE PER UNIT-DAY, BY SYSTEM ===\n")
print(by_system %>% dplyr::mutate(dplyr::across(c(mean_per_unit_day, median_per_unit_day),
                                                ~ round(.x, 1))), n = Inf, width = Inf)
cat("\n=== BY UNIT (22 days) ===\n")
print(by_unit %>% dplyr::mutate(mean_per_day = round(mean_per_day, 1)), n = Inf, width = Inf)


# ==============================================================================
# SECTION 3: FIGURE
# ==============================================================================

system_daily <- abund %>%
  dplyr::group_by(Device_type, Block, Date) %>%
  dplyr::summarise(m = mean(abundance), .groups = "drop") %>%
  dplyr::group_by(Device_type, Date) %>%
  dplyr::summarise(mean_abundance = mean(m), .groups = "drop")

q98   <- stats::quantile(abund$abundance, 0.98, na.rm = TRUE)
step  <- 10^floor(log10(q98)) / 2
y_cap <- ceiling(q98 / step) * step

p <- ggplot2::ggplot() +
  ggplot2::annotate("rect", xmin = min(outage_dates) - 0.5, xmax = max(outage_dates) + 0.5,
                    ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  ggplot2::geom_point(data = abund,
                      ggplot2::aes(x = Date, y = abundance, color = Device_type, shape = Effort),
                      alpha = 0.6, size = 1.8) +
  ggplot2::geom_line(data = system_daily,
                     ggplot2::aes(x = Date, y = mean_abundance, color = Device_type),
                     linewidth = 1) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = c("full effort" = 16, "outage window" = 1)) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
  ggplot2::labs(x = "Date", y = "Abundance (individuals per unit and day)",
                color = "System", shape = NULL,
                caption = sprintf(paste0("Lines: daily system means (within block, then across blocks). ",
                                         "Shaded: ID power-outage window. View capped at %s ",
                                         "(rounded 98th percentile)."),
                                  scales::label_comma()(y_cap))) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position  = "top",
                 plot.caption     = ggplot2::element_text(size = 8, color = "grey50"))

print(p)
f_fig <- file.path(FIG_DIR, paste0(SCRIPT, "_daily_abundance.png"))
ggplot2::ggsave(f_fig, p, width = 10, height = 6, dpi = 300, bg = "white")
cat("\n[OK] Saved:", f_fig, "\n")


# ==============================================================================
# SECTION 4: EXPORT
# ==============================================================================

readr::write_csv(by_system, file.path(TAB_DIR, paste0(SCRIPT, "_by_system.csv")))
readr::write_csv(by_unit,   file.path(TAB_DIR, paste0(SCRIPT, "_by_unit.csv")))
cat("[OK] Tables written to", TAB_DIR, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

while (sink.number() > 0) sink()
