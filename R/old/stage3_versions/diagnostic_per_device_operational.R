# ==============================================================================
# diagnostic_per_device_operational.R
# Throwaway diagnostic: is the central ID peak (late August) that FAIRD does not
# follow present at ALL co-located sites, or driven by a single device?
# Splits the operational-period abundance by Site, plotting FAIRD and ID side by
# side within each site. FAIRD and ID coexist only at Site1, Site2, Site4
# (Site3 has no ID), so those are the comparable panels.
# Not part of the pipeline; run once, look, discard.
# ==============================================================================

library(tidyverse)
library(splines)
library(ggplot2)

source("R/00_load_data.R")
source("R/00_custom_colors.R")

# --- Operational-period FAIRD + ID daily abundance ----------------------------
# Same construction as the main script operational run: apply the ID battery
# flag, keep only rows valid for the operational window.
abund_op <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, abundance) %>%
  dplyr::left_join(
    data_list$effort_base %>%
      dplyr::select(Device, Date, Field_Valid_Operational_ID),
    by = c("Device", "Date")
  ) %>%
  dplyr::filter(Field_Valid_Operational_ID %in% TRUE)

# Keep only sites where FAIRD and ID coexist (both device types present).
sites_both <- abund_op %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(n_types = dplyr::n_distinct(Device_type), .groups = "drop") %>%
  dplyr::filter(n_types == 2) %>%
  dplyr::pull(Site)

cat("=== SITES WITH BOTH FAIRD AND ID (operational) ===\n")
cat(sprintf("  %s\n", paste(sort(sites_both), collapse = ", ")))

abund_op_both <- abund_op %>% dplyr::filter(Site %in% sites_both)

# --- Numeric summary: where does each device peak? ----------------------------
# For each device, the date of its maximum daily abundance and that value.
# If the ID peak is at one device only, this table shows it before the plot.
cat("\n=== PER-DEVICE PEAK (operational period) ===\n\n")
peak_tbl <- abund_op_both %>%
  dplyr::group_by(Site, Device_type, Device) %>%
  dplyr::summarise(
    peak_date  = Date[which.max(abundance)],
    peak_value = max(abundance),
    mean_abund = round(mean(abundance), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Site, Device_type)
print(peak_tbl)

# Mean abundance in the late-August central window (the peak region), per
# device, to quantify who is high there. Window chosen around the visible
# ID hump (Aug 27 to Sep 1).
central_start <- as.Date("2023-08-27")
central_end   <- as.Date("2023-09-01")
cat(sprintf("\n=== MEAN ABUNDANCE IN CENTRAL WINDOW (%s to %s) ===\n\n",
            central_start, central_end))
central_tbl <- abund_op_both %>%
  dplyr::filter(Date >= central_start, Date <= central_end) %>%
  dplyr::group_by(Site, Device_type, Device) %>%
  dplyr::summarise(mean_central = round(mean(abundance), 1),
                   max_central  = max(abundance),
                   .groups = "drop") %>%
  dplyr::arrange(Site, Device_type)
print(central_tbl)

# --- Faceted plot: FAIRD vs ID per site ---------------------------------------
# One panel per site. If the ID-above-FAIRD central divergence appears in all
# panels, it is a system-level pattern; if in one, it is device-specific.
abund_op_both <- abund_op_both %>%
  dplyr::mutate(Date.num = as.numeric(Date))

p <- ggplot2::ggplot(abund_op_both,
                     ggplot2::aes(x = Date, y = abundance, color = Device_type)) +
  ggplot2::geom_point(alpha = 0.6, size = 1.8) +
  ggplot2::geom_smooth(
    ggplot2::aes(fill = Device_type),
    method = "lm", formula = y ~ splines::bs(x, 3),
    se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  ggplot2::facet_wrap(~ Site, ncol = 1, scales = "fixed") +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_fill_manual(values = custom_colors) +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::labs(
    title    = "Per-device FAIRD vs ID abundance (operational period)",
    subtitle = "One panel per co-located site. Is the central ID peak in all sites or one?",
    x = "Date", y = "Abundance (individuals / day)",
    color = "Device_type", fill = "Device_type",
    caption = "Diagnostic only. y-axis capped at 100 for readability."
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
    strip.text       = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top"
  )

out_file <- "outputs/figures/diagnostic_per_device_operational.png"
print(p)
ggplot2::ggsave(out_file, p, width = 9, height = 10, dpi = 200, bg = "white")
cat(sprintf("\n[OK] Saved: %s\n", out_file))
