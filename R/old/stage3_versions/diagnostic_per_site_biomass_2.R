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
# Data treatment matches the main script (conservative run). Days inside the
# ID battery window (Field_Valid_Operational_ID not TRUE) are FLAGGED and
# shaded in the figures. The data show that ID kept recording on most of those
# days (non-zero biomass), so they are days of PARTIAL and unknown sampling
# effort, not forced zeros. ID values inside the window are therefore not
# comparable with full-effort days, in either direction.
#
# Outputs (diagnostic only, not part of the manuscript pipeline):
#   console: daily values per site, window profile per device, daily ratio
#            per pair of systems (site, block, intra-FAIRD check), direction
#            concordance per site
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
cat("ID battery-window days are flagged as PARTIAL-EFFORT days, not as zeros.\n")
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
    ID_battery_window = Device_type == "ID" &
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

# Verify the flag before using it: dates flagged as ID battery window.
window_dates <- bio_site %>%
  dplyr::filter(ID_battery_window) %>%
  dplyr::pull(Date)
cat("\nID battery-window days flagged (Field_Valid_Operational_ID not TRUE):\n")
print(bio_site %>%
        dplyr::filter(ID_battery_window) %>%
        dplyr::group_by(Device) %>%
        dplyr::summarise(n = dplyr::n(), first = min(Date), last = max(Date),
                         n_zero_days = sum(Biomass_mg == 0),
                         biomass_sum = sum(Biomass_mg), .groups = "drop"))
cat("  n_zero_days < n means the device recorded on some flagged days.\n")
has_window <- length(window_dates) > 0

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
cat("ID_window = TRUE marks ID battery-window days (partial, unknown effort).\n\n")

daily_wide <- bio_site %>%
  dplyr::mutate(Biomass_mg = round(Biomass_mg, 1)) %>%
  dplyr::select(Site, Date, Device_type, Biomass_mg) %>%
  tidyr::pivot_wider(names_from = Device_type, values_from = Biomass_mg,
                     names_sort = TRUE) %>%
  dplyr::left_join(
    bio_site %>%
      dplyr::group_by(Site, Date) %>%
      dplyr::summarise(ID_window = any(ID_battery_window), .groups = "drop"),
    by = c("Site", "Date")
  ) %>%
  dplyr::arrange(Site, Date)
print(daily_wide, n = Inf)


# ==============================================================================
# TABLE 2: window profile (early / central / late thirds of the window)
# ==============================================================================
# The window is split into three equal blocks of consecutive dates. For each
# device: mean over all days and mean over valid days only (ID battery-window
# days removed), plus the ratio of each later window to the early window. A
# fall-and-recovery shape appears as central/early < 1 and late/early > central.
# n_valid shows how many full-effort days support mean_valid: for ID the
# central and late windows rest on very few days.
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
    n_window   = sum(ID_battery_window),
    n_valid    = sum(!ID_battery_window),
    mean_all   = mean(Biomass_mg),
    mean_valid = if (any(!ID_battery_window)) mean(Biomass_mg[!ID_battery_window])
                 else NA_real_,
    .groups = "drop"
  )
cat("\nPer device and window (n_window = ID battery-window days, n_valid = full-effort days):\n")
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
# TABLE 3: daily ratio per pair of systems
# ==============================================================================
# Under the log link of the Tweedie model, parallel trajectories mean a
# constant ratio. Three readings, matching the pair models of the main script:
#   3a, per site (co-located site-days): AMMOD/FAIRD, AMMOD/ID, FAIRD/ID.
#   3b, per block (PAIRED_BLOCK, the primary pairing of the main script):
#       ratio of DEVICE MEANS per block-day, so that unequal device numbers
#       per type (Meadow: AMMOD3 vs FAIRD3 + FAIRD4) do not bias it.
#   3c, intra-system FAIRD ratio between the two sites of each block
#       (FAIRD1/FAIRD2, FAIRD3/FAIRD4): internal check of the assumption behind
#       PAIRED_BLOCK that sites within a block share their temporal dynamics.
#       Diagnostic only, not for the manuscript.
# Summaries over all days and over days outside the ID battery window. The
# Spearman rho of the ratio against date is a descriptive drift indicator
# (positive: numerator grows relative to denominator). No p-values: n is 5 to
# 12 days per unit.
cat("\n=== TABLE 3: DAILY RATIO PER PAIR OF SYSTEMS ===\n")

ratio_pairs <- list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))
id_window_dates <- unique(window_dates)

bio_block <- bio_site %>%
  dplyr::mutate(
    Block       = dplyr::if_else(Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2"),
    Device_type = as.character(Device_type)
  )

# 3a: co-located site-days.
ratio_3a <- dplyr::bind_rows(lapply(ratio_pairs, function(pr) {
  bio_block %>%
    dplyr::filter(Device_type %in% pr) %>%
    dplyr::select(Site, Date, Device_type, Biomass_mg) %>%
    tidyr::pivot_wider(names_from = Device_type, values_from = Biomass_mg) %>%
    dplyr::filter(!is.na(.data[[pr[1]]]), !is.na(.data[[pr[2]]])) %>%
    dplyr::transmute(Level = "3a_site", Pair = paste(pr, collapse = "/"),
                     Unit = as.character(Site), Date,
                     num = .data[[pr[1]]], den = .data[[pr[2]]])
}))

# 3b: block-days, ratio of device means.
ratio_3b <- dplyr::bind_rows(lapply(ratio_pairs, function(pr) {
  bio_block %>%
    dplyr::filter(Device_type %in% pr) %>%
    dplyr::group_by(Block, Date, Device_type) %>%
    dplyr::summarise(mean_mg = mean(Biomass_mg), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Device_type, values_from = mean_mg) %>%
    dplyr::filter(!is.na(.data[[pr[1]]]), !is.na(.data[[pr[2]]])) %>%
    dplyr::transmute(Level = "3b_block", Pair = paste(pr, collapse = "/"),
                     Unit = Block, Date,
                     num = .data[[pr[1]]], den = .data[[pr[2]]])
}))

# 3c: FAIRD between the two sites of each block.
ratio_3c <- bio_block %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::mutate(role = dplyr::if_else(as.character(Site) %in% c("Site1", "Site3"),
                                      "num", "den")) %>%
  dplyr::select(Block, Date, role, Device, Biomass_mg) %>%
  dplyr::mutate(Device = as.character(Device)) %>%
  tidyr::pivot_wider(names_from = role, values_from = c(Device, Biomass_mg)) %>%
  dplyr::filter(!is.na(Biomass_mg_num), !is.na(Biomass_mg_den)) %>%
  dplyr::transmute(Level = "3c_intra_FAIRD", Pair = paste(Device_num, Device_den, sep = "/"),
                   Unit = Block, Date, num = Biomass_mg_num, den = Biomass_mg_den)

# Days on which the block has AMMOD biomass: the days that enter the AMMOD
# pairs under PAIRED_BLOCK (Meadow: AMMOD3 days only).
ammod_block_days <- bio_block %>%
  dplyr::filter(Device_type == "AMMOD") %>%
  dplyr::distinct(Block, Date) %>%
  dplyr::mutate(ammod_day = TRUE)

ratio_all <- dplyr::bind_rows(ratio_3a, ratio_3b, ratio_3c) %>%
  dplyr::mutate(
    ratio     = dplyr::if_else(den > 0, num / den, NA_real_),
    in_window = Date %in% id_window_dates,
    Block     = dplyr::if_else(Unit %in% c("Site1", "Site2", "BLOCK1"), "BLOCK1", "BLOCK2")
  ) %>%
  dplyr::left_join(ammod_block_days, by = c("Block", "Date")) %>%
  dplyr::mutate(ammod_day = dplyr::coalesce(ammod_day, FALSE)) %>%
  dplyr::arrange(Level, Pair, Unit, Date)

summarise_ratio <- function(d, subset_label) {
  d %>%
    dplyr::group_by(Level, Pair, Unit) %>%
    dplyr::summarise(
      Subset            = subset_label,
      n_days            = dplyr::n(),
      n_den_zero        = sum(den == 0),
      ratio_overall     = sum(num) / sum(den),
      ratio_median      = stats::median(ratio, na.rm = TRUE),
      ratio_min         = suppressWarnings(min(ratio, na.rm = TRUE)),
      ratio_max         = suppressWarnings(max(ratio, na.rm = TRUE)),
      rho_ratio_vs_date = suppressWarnings(
        stats::cor(as.numeric(Date), ratio, method = "spearman",
                   use = "complete.obs")),
      .groups = "drop"
    )
}

ratio_summary <- dplyr::bind_rows(
  summarise_ratio(ratio_all, "all_days"),
  summarise_ratio(ratio_all %>% dplyr::filter(!in_window), "outside_ID_window"),
  summarise_ratio(ratio_all %>% dplyr::filter(Level == "3c_intra_FAIRD", ammod_day),
                  "AMMOD_biomass_days")
) %>%
  dplyr::arrange(Level, Pair, Unit, Subset)

cat("\nSummary per pair and unit (3a site, 3b block, 3c intra-FAIRD check):\n")
print(ratio_summary %>%
        dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ round(.x, 2))),
      n = Inf)

cat("\nDaily ratios, 3b (block level) and 3c (intra-FAIRD check):\n")
print(ratio_all %>%
        dplyr::filter(Level != "3a_site") %>%
        dplyr::select(Level, Pair, Unit, Date, num, den, ratio, in_window, ammod_day) %>%
        dplyr::mutate(dplyr::across(c(num, den), ~ round(.x, 1)),
                      ratio = round(ratio, 2)),
      n = Inf)
cat("  in_window = ID battery-window day; ammod_day = block has AMMOD biomass.\n")

# ==============================================================================
# TABLE 4: direction concordance of consecutive daily changes, per site
# ==============================================================================
# Percentage of consecutive day-pairs in which two co-located systems change in
# the same direction. Intervals touching an ID battery-window day and intervals
# with a zero change in either series are excluded.
cat("\n=== TABLE 4: DIRECTION CONCORDANCE PER SITE ===\n\n")
diffs <- bio_site %>%
  dplyr::arrange(Site, Device_type, Date) %>%
  dplyr::group_by(Site, Device_type) %>%
  dplyr::mutate(
    consecutive     = as.numeric(Date - dplyr::lag(Date)) == 1,
    delta           = Biomass_mg - dplyr::lag(Biomass_mg),
    touches_window = ID_battery_window | dplyr::lag(ID_battery_window, default = FALSE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(consecutive %in% TRUE) %>%
  dplyr::select(Site, Date, Device_type, delta, touches_window)

type_pairs <- list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))
concordance <- dplyr::bind_rows(lapply(type_pairs, function(pr) {
  a <- diffs %>% dplyr::filter(Device_type == pr[1])
  b <- diffs %>% dplyr::filter(Device_type == pr[2])
  dplyr::inner_join(a, b, by = c("Site", "Date"), suffix = c("_a", "_b")) %>%
    dplyr::filter(!touches_window_a, !touches_window_b,
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
if (has_window) {
  p_a <- p_a + ggplot2::annotate(
    "rect", xmin = min(window_dates) - 0.5, xmax = max(window_dates) + 0.5,
    ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.12
  )
}
p_a <- p_a +
  ggplot2::geom_line(alpha = 0.35, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(shape = ID_battery_window), size = 2, alpha = 0.8) +
  ggplot2::geom_smooth(
    method = "glm", formula = y ~ splines::bs(x, 3),
    method.args = list(family = stats::quasipoisson(link = "log")),
    se = FALSE, linewidth = 1
  ) +
  ggplot2::facet_wrap(~ Site, ncol = 1, scales = "free_y") +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                              name = "ID battery window") +
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
                      "Shaded band: ID battery window (affects ID only; ",
                      "open circles = ID days of partial, unknown effort). Lines: log-link spline smooths.")
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
# that site has a value and no ID battery-window day occurs. Using the same
# days for all devices of a site keeps the normalisations comparable (AMMOD3
# has 8 days, ID has partial-effort days).
n_types_site <- bio_site %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(n_types = dplyr::n_distinct(Device_type), .groups = "drop")

common_days <- bio_site %>%
  dplyr::group_by(Site, Date) %>%
  dplyr::summarise(n_present = dplyr::n_distinct(Device_type),
                   any_window  = any(ID_battery_window), .groups = "drop") %>%
  dplyr::left_join(n_types_site, by = "Site") %>%
  dplyr::filter(n_present == n_types, !any_window) %>%
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
if (has_window) {
  p_b <- p_b + ggplot2::annotate(
    "rect", xmin = min(window_dates) - 0.5, xmax = max(window_dates) + 0.5,
    ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.12
  )
}
p_b <- p_b +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  ggplot2::geom_line(alpha = 0.35, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(shape = ID_battery_window), size = 2, alpha = 0.8) +
  ggplot2::geom_smooth(
    method = "glm", formula = y ~ splines::bs(x, 3),
    method.args = list(family = stats::quasipoisson(link = "log")),
    se = FALSE, linewidth = 1
  ) +
  ggplot2::facet_wrap(~ Site, ncol = 1) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                              name = "ID battery window") +
  ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  ggplot2::coord_cartesian(ylim = c(0, y_cap_b)) +
  ggplot2::labs(
    title    = "Per-site biomass trajectory SHAPE (each device / its own mean)",
    subtitle = "1 = device mean over the common valid days of its site. Linear axis.",
    x = "Date", y = "Relative biomass", color = "Device_type",
    caption  = sprintf(paste0("Diagnostic only. View capped at %s (rounded 98th percentile). ",
                              "Shaded band: ID battery window, open circles = ID days of partial, unknown effort."),
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
