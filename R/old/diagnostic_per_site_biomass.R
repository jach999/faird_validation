# ==============================================================================
# diagnostic_per_site_biomass.R
# Per-site decomposition of the daily biomass trajectories in the AMMOD overlap
# window (Batch 1). Biomass analogue of diagnostic_per_device_operational.R.
#
# Question: is the aggregate pattern of the temporal biomass model
# (Scherber_temporal_biomass.R) present at every site, or is it the
# superposition of unrelated site-level dynamics? With n = 1 device per type
# per site the aggregate model cannot answer this by itself.
#
# Data treatment matches the main script (conservative run: ID battery-failure
# days enter as zeros). Those days are FLAGGED via Field_Valid_Operational_ID
# and shaded in the figures: a dip of the ID series inside that window is a
# power failure, not a biological signal.
#
# Outputs (diagnostic only, not part of the manuscript pipeline):
#   console: daily values per site, window profile per device, AMMOD/FAIRD
#            daily ratio per site, direction concordance per site
#   figures: (A) raw biomass, dual axis, one panel per site
#            (B) shape only: each device relative to its own mean over the
#                site-days on which all its co-located devices are valid
# ==============================================================================

library(tidyverse)
library(splines)
library(ggplot2)

source("R/00_load_data.R")
source("R/00_custom_colors.R")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console", showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/diagnostic_per_site_biomass_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("PER-SITE BIOMASS DIAGNOSTIC (AMMOD overlap window, conservative ID run)\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Helpers (identical to Scherber_temporal_biomass.R) -----------------------
nice_ceiling <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1 || !is.finite(x) || x <= 0) return(NA_real_)
  step <- 10^floor(log10(x)) / 2
  ceiling(x / step) * step
}

# --- Data: same construction as the main script -------------------------------
etraps_bio <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, biomass_mg) %>%
  dplyr::rename(Biomass_mg = biomass_mg)

ammod_bio <- data_list$ammod_biomass_for_stats %>%
  dplyr::mutate(Device_type = "AMMOD", Biomass_mg = Live_mass) %>%
  dplyr::select(Date, Device, Device_type, Site, Biomass_mg)

date_range <- range(ammod_bio$Date)

bio_site <- dplyr::bind_rows(etraps_bio, ammod_bio) %>%
  dplyr::filter(Date >= date_range[1], Date <= date_range[2],
                !is.na(Biomass_mg)) %>%
  dplyr::left_join(
    data_list$effort_base %>%
      dplyr::select(Device, Date, Field_Valid_Operational_ID),
    by = c("Device", "Date")
  ) %>%
  dplyr::mutate(
    ID_power_failure = Device_type == "ID" &
      !(Field_Valid_Operational_ID %in% TRUE),
    Device_type = factor(Device_type, levels = c("AMMOD", "FAIRD", "ID")),
    Site        = factor(Site)
  ) %>%
  dplyr::arrange(Site, Device_type, Date)

cat(sprintf("Overlap window: %s to %s | %d observations\n",
            date_range[1], date_range[2], nrow(bio_site)))
cat("\nDevices present per site:\n")
print(bio_site %>%
        dplyr::count(Site, Device_type, Device, name = "n_days"),
      n = Inf)

# Verify the flag before using it: dates flagged as ID power failure.
fail_dates <- bio_site %>%
  dplyr::filter(ID_power_failure) %>%
  dplyr::pull(Date)
cat("\nID power-failure days flagged (Field_Valid_Operational_ID not TRUE):\n")
print(bio_site %>%
        dplyr::filter(ID_power_failure) %>%
        dplyr::group_by(Device) %>%
        dplyr::summarise(n = dplyr::n(), first = min(Date), last = max(Date),
                         biomass_sum = sum(Biomass_mg), .groups = "drop"))
has_failure <- length(fail_dates) > 0

# --- Scale factor for the dual axis (same definition as main script) ----------
paired_af <- bio_site %>%
  dplyr::filter(Device_type %in% c("AMMOD", "FAIRD")) %>%
  dplyr::select(Site, Date, Device_type, Biomass_mg) %>%
  dplyr::mutate(Device_type = as.character(Device_type)) %>%
  tidyr::pivot_wider(names_from = Device_type, values_from = Biomass_mg) %>%
  dplyr::filter(!is.na(AMMOD), !is.na(FAIRD))
scale_factor_ammod <- signif(sum(paired_af$AMMOD) / sum(paired_af$FAIRD), 2)
cat(sprintf("\nDual-axis scale factor (sum AMMOD / sum FAIRD, co-located site-days): %s\n",
            scale_factor_ammod))


# ==============================================================================
# TABLE 1: daily values per site (wide)
# ==============================================================================
cat("\n=== TABLE 1: DAILY BIOMASS (mg) PER SITE ===\n")
cat("ID_fail = TRUE marks ID battery-failure days (ID value is a forced zero).\n\n")

daily_wide <- bio_site %>%
  dplyr::mutate(Biomass_mg = round(Biomass_mg, 1)) %>%
  dplyr::select(Site, Date, Device_type, Biomass_mg) %>%
  tidyr::pivot_wider(names_from = Device_type, values_from = Biomass_mg,
                     names_sort = TRUE) %>%
  dplyr::left_join(
    bio_site %>%
      dplyr::group_by(Site, Date) %>%
      dplyr::summarise(ID_fail = any(ID_power_failure), .groups = "drop"),
    by = c("Site", "Date")
  ) %>%
  dplyr::arrange(Site, Date)
print(daily_wide, n = Inf)


# ==============================================================================
# TABLE 2: window profile (early / central / late thirds of the window)
# ==============================================================================
# The window is split into three equal blocks of consecutive dates. For each
# device: mean over all days and mean over valid days only (ID failure days
# removed), plus the ratio of each later window to the early window. A
# fall-and-recovery shape appears as central/early < 1 and late/early > central.
all_dates <- sort(unique(bio_site$Date))
window_lookup <- tibble::tibble(
  Date   = all_dates,
  Window = cut(seq_along(all_dates), breaks = 3,
               labels = c("early", "central", "late"))
)
cat("\n=== TABLE 2: WINDOW PROFILE PER DEVICE ===\n\nWindow definition:\n")
print(window_lookup %>%
        dplyr::group_by(Window) %>%
        dplyr::summarise(from = min(Date), to = max(Date), n_days = dplyr::n(),
                         .groups = "drop"))

window_long <- bio_site %>%
  dplyr::left_join(window_lookup, by = "Date") %>%
  dplyr::group_by(Site, Device_type, Device, Window) %>%
  dplyr::summarise(
    n_days     = dplyr::n(),
    n_failure  = sum(ID_power_failure),
    mean_all   = mean(Biomass_mg),
    mean_valid = if (any(!ID_power_failure)) mean(Biomass_mg[!ID_power_failure])
                 else NA_real_,
    .groups = "drop"
  )
cat("\nPer device and window (n_failure = ID battery days inside the window):\n")
print(window_long %>%
        dplyr::mutate(dplyr::across(c(mean_all, mean_valid), ~ round(.x, 1))),
      n = Inf)

window_ratio <- window_long %>%
  dplyr::select(Site, Device_type, Device, Window, mean_valid) %>%
  tidyr::pivot_wider(names_from = Window, values_from = mean_valid) %>%
  dplyr::mutate(
    central_vs_early = dplyr::if_else(early > 0, central / early, NA_real_),
    late_vs_early    = dplyr::if_else(early > 0, late / early, NA_real_)
  )
cat("\nRatios to the early window (valid days only):\n")
print(window_ratio %>%
        dplyr::mutate(dplyr::across(c(early, central, late), ~ round(.x, 1)),
                      dplyr::across(c(central_vs_early, late_vs_early), ~ round(.x, 2))),
      n = Inf)


# ==============================================================================
# TABLE 3: AMMOD/FAIRD daily ratio per site
# ==============================================================================
# Under the log link of the Tweedie model, parallel trajectories mean a
# constant ratio. The Spearman rho of the ratio against date is a descriptive
# drift indicator per site (positive: AMMOD grows relative to FAIRD). No
# p-values: n is 8 to 12 days per site.
cat("\n=== TABLE 3: AMMOD/FAIRD DAILY RATIO PER SITE ===\n\n")
ratio_site <- paired_af %>%
  dplyr::mutate(ratio = dplyr::if_else(FAIRD > 0, AMMOD / FAIRD, NA_real_)) %>%
  dplyr::arrange(Site, Date)

print(ratio_site %>%
        dplyr::mutate(dplyr::across(c(AMMOD, FAIRD), ~ round(.x, 1)),
                      ratio = round(ratio, 1)),
      n = Inf)

ratio_summary <- ratio_site %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    n_days        = dplyr::n(),
    n_FAIRD_zero  = sum(FAIRD == 0),
    ratio_overall = sum(AMMOD) / sum(FAIRD),
    ratio_median  = stats::median(ratio, na.rm = TRUE),
    ratio_min     = min(ratio, na.rm = TRUE),
    ratio_max     = max(ratio, na.rm = TRUE),
    rho_ratio_vs_date = suppressWarnings(
      stats::cor(as.numeric(Date), ratio, method = "spearman",
                 use = "complete.obs")),
    .groups = "drop"
  )
cat("\nSummary per site:\n")
print(ratio_summary %>%
        dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ round(.x, 2))))


# ==============================================================================
# TABLE 4: direction concordance of consecutive daily changes, per site
# ==============================================================================
# Percentage of consecutive day-pairs in which two co-located systems change in
# the same direction. Intervals touching an ID failure day and intervals with a
# zero change in either series are excluded.
cat("\n=== TABLE 4: DIRECTION CONCORDANCE PER SITE ===\n\n")
diffs <- bio_site %>%
  dplyr::arrange(Site, Device_type, Date) %>%
  dplyr::group_by(Site, Device_type) %>%
  dplyr::mutate(
    consecutive     = as.numeric(Date - dplyr::lag(Date)) == 1,
    delta           = Biomass_mg - dplyr::lag(Biomass_mg),
    touches_failure = ID_power_failure | dplyr::lag(ID_power_failure, default = FALSE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(consecutive %in% TRUE) %>%
  dplyr::select(Site, Date, Device_type, delta, touches_failure)

type_pairs <- list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))
concordance <- dplyr::bind_rows(lapply(type_pairs, function(pr) {
  a <- diffs %>% dplyr::filter(Device_type == pr[1])
  b <- diffs %>% dplyr::filter(Device_type == pr[2])
  dplyr::inner_join(a, b, by = c("Site", "Date"), suffix = c("_a", "_b")) %>%
    dplyr::filter(!touches_failure_a, !touches_failure_b,
                  delta_a != 0, delta_b != 0) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      Pair = paste(pr, collapse = " vs "),
      n_intervals = dplyr::n(),
      pct_same_direction = round(100 * mean(sign(delta_a) == sign(delta_b)), 1),
      .groups = "drop"
    )
}))
print(concordance %>% dplyr::arrange(Pair, Site), n = Inf)


# ==============================================================================
# FIGURE A: raw biomass, dual axis, one panel per site
# ==============================================================================
# Smooths are per-device quasi-Poisson GLMs on bs(Date, 3) with a log link,
# the closest descriptive analogue of the Tweedie spline model. No CI ribbons:
# with free y-scales, unstable edge intervals would set the panel range.
k <- scale_factor_ammod
plot_a <- bio_site %>%
  dplyr::mutate(y_plot = dplyr::if_else(Device_type == "AMMOD",
                                        Biomass_mg / k, Biomass_mg))

p_a <- ggplot2::ggplot(plot_a, ggplot2::aes(x = Date, y = y_plot,
                                            color = Device_type))
if (has_failure) {
  p_a <- p_a + ggplot2::annotate(
    "rect", xmin = min(fail_dates) - 0.5, xmax = max(fail_dates) + 0.5,
    ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.12
  )
}
p_a <- p_a +
  ggplot2::geom_line(alpha = 0.35, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(shape = ID_power_failure), size = 2, alpha = 0.8) +
  ggplot2::geom_smooth(
    method = "glm", formula = y ~ splines::bs(x, 3),
    method.args = list(family = stats::quasipoisson(link = "log")),
    se = FALSE, linewidth = 1
  ) +
  ggplot2::facet_wrap(~ Site, ncol = 1, scales = "free_y") +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                              name = "ID power failure") +
  ggplot2::scale_y_continuous(
    name     = "FAIRD / ID Biomass (mg)",
    labels   = scales::label_comma(),
    sec.axis = ggplot2::sec_axis(~ . * k,
                                 name   = sprintf("AMMOD Biomass (mg), axis = left x %s", k),
                                 labels = scales::label_comma())
  ) +
  ggplot2::labs(
    title    = "Per-site daily biomass, AMMOD overlap window",
    subtitle = "One panel per site. Is the aggregate trajectory pattern present at every site?",
    x = "Date", color = "Device_type",
    caption  = paste0("Diagnostic only. Linear axes, free y per panel. ",
                      "Shaded band: ID battery-failure days (affects ID only; ",
                      "open circles = forced zeros). Lines: log-link spline smooths.")
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title         = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle      = ggplot2::element_text(size = 10, color = "grey40"),
    plot.caption       = ggplot2::element_text(size = 8, color = "grey50"),
    strip.text         = ggplot2::element_text(face = "bold"),
    panel.grid.minor   = ggplot2::element_blank(),
    legend.position    = "top",
    axis.title.y.right = ggplot2::element_text(color = "red"),
    axis.text.y.right  = ggplot2::element_text(color = "red")
  )

out_a <- "outputs/figures/diagnostic_per_site_biomass_dualaxis.png"
print(p_a)
ggplot2::ggsave(out_a, p_a, width = 9, height = 12, dpi = 200, bg = "white")
cat(sprintf("\n[OK] Saved: %s\n", out_a))


# ==============================================================================
# FIGURE B: shape only, each device relative to its own mean
# ==============================================================================
# Removes the magnitude difference between systems without a log axis and
# without an arbitrary scale factor. Each device is divided by its mean over
# the COMMON VALID DAYS of its site: days on which every device present at
# that site has a value and no ID failure occurs. Using the same days for all
# devices of a site keeps the normalisations comparable (AMMOD3 has 8 days, ID
# has forced zeros).
n_types_site <- bio_site %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(n_types = dplyr::n_distinct(Device_type), .groups = "drop")

common_days <- bio_site %>%
  dplyr::group_by(Site, Date) %>%
  dplyr::summarise(n_present = dplyr::n_distinct(Device_type),
                   any_fail  = any(ID_power_failure), .groups = "drop") %>%
  dplyr::left_join(n_types_site, by = "Site") %>%
  dplyr::filter(n_present == n_types, !any_fail) %>%
  dplyr::select(Site, Date)

cat("\n=== FIGURE B: COMMON VALID DAYS USED FOR NORMALISATION ===\n\n")
print(common_days %>%
        dplyr::group_by(Site) %>%
        dplyr::summarise(n_days = dplyr::n(), from = min(Date), to = max(Date),
                         .groups = "drop"))

device_ref <- bio_site %>%
  dplyr::semi_join(common_days, by = c("Site", "Date")) %>%
  dplyr::group_by(Device) %>%
  dplyr::summarise(ref_mean = mean(Biomass_mg), .groups = "drop")

plot_b <- bio_site %>%
  dplyr::left_join(device_ref, by = "Device") %>%
  dplyr::mutate(rel = dplyr::if_else(ref_mean > 0, Biomass_mg / ref_mean, NA_real_))

y_cap_b <- nice_ceiling(stats::quantile(plot_b$rel, 0.98, na.rm = TRUE))

p_b <- ggplot2::ggplot(plot_b, ggplot2::aes(x = Date, y = rel, color = Device_type))
if (has_failure) {
  p_b <- p_b + ggplot2::annotate(
    "rect", xmin = min(fail_dates) - 0.5, xmax = max(fail_dates) + 0.5,
    ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.12
  )
}
p_b <- p_b +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  ggplot2::geom_line(alpha = 0.35, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(shape = ID_power_failure), size = 2, alpha = 0.8) +
  ggplot2::geom_smooth(
    method = "glm", formula = y ~ splines::bs(x, 3),
    method.args = list(family = stats::quasipoisson(link = "log")),
    se = FALSE, linewidth = 1
  ) +
  ggplot2::facet_wrap(~ Site, ncol = 1) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                              name = "ID power failure") +
  ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  ggplot2::coord_cartesian(ylim = c(0, y_cap_b)) +
  ggplot2::labs(
    title    = "Per-site biomass trajectory SHAPE (each device / its own mean)",
    subtitle = "1 = device mean over the common valid days of its site. Linear axis.",
    x = "Date", y = "Relative biomass", color = "Device_type",
    caption  = sprintf(paste0("Diagnostic only. View capped at %s (rounded 98th percentile). ",
                              "Shaded band: ID battery-failure days, open circles = forced zeros."),
                       y_cap_b)
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
    plot.caption     = ggplot2::element_text(size = 8, color = "grey50"),
    strip.text       = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top"
  )

out_b <- "outputs/figures/diagnostic_per_site_biomass_relative.png"
print(p_b)
ggplot2::ggsave(out_b, p_b, width = 9, height = 12, dpi = 200, bg = "white")
cat(sprintf("[OK] Saved: %s\n", out_b))

cat("\n==============================================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================================\n")

while (sink.number() > 0) sink()
