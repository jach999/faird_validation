# ==============================================================================
# supplementary/S3_site_diagnostics.R
# Supplementary Material S3: per-site diagnostics
# ------------------------------------------------------------------------------
# With one unit per system and site, a system-level pattern can be the
# superposition of unrelated site dynamics. These diagnostics show each site
# separately; they are descriptive and enter no test.
#
# Part A, biomass, 12 days with AMMOD biomass (all ID days kept): daily values
#   per site, profile of the ID power-outage window per unit, daily ratio per
#   pair of systems (by site, by block, and between the two FAIRD units of
#   each block), direction concordance per site. Days inside the outage
#   window are flagged as days of partial, unknown recording effort (ID kept
#   recording on most of them), not as zeros.
#   Figures: (A) biomass on a dual axis, one panel per site; (B) each unit
#   relative to its own mean over the site-days on which all co-located
#   units are valid.
# Part B, abundance of FAIRD and ID, all days with the ID power-outage window
#   flagged for both systems: per-unit peak (all days and full-effort days),
#   mean in the central window of the season, one panel per site with both
#   camera systems (Site1, Site2, Site4; Site3 has no ID data).
# Part C, ID detection classes (Pre-classification) per unit: daily counts per
#   class on the outage-window days against the same class on full-effort days,
#   and class shares per day. Shows whether the high counts of some
#   outage-window days sit in one class or are spread over all classes. ID body
#   length and live mass are class means, so classes cannot separate individuals.
# ==============================================================================

required_packages <- c("tidyverse", "here", "splines", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)
library(splines)
library(ggplot2)

source(here::here("R", "00_load_data.R"))
source(here::here("R", "00_custom_colors.R"))

SCRIPT  <- "S3_site_diagnostics"
FIG_DIR <- here::here("outputs", "figures")
CON_DIR <- here::here("outputs", "console")
TAB_DIR <- here::here("outputs", "tables")
for (d in c(FIG_DIR, CON_DIR, TAB_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Central window of the season for Part B (late-August peak of ID abundance).
CENTRAL_START <- as.Date("2023-08-27")
CENTRAL_END   <- as.Date("2023-09-01")

# Part C: number of most frequent ID classes shown separately (display
# parameter); the remaining classes are pooled as "other classes".
N_TOP_CLASSES <- 3

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, paste0(SCRIPT, "_output.txt")), split = TRUE)

cat("==============================================================================\n")
cat("SUPPLEMENTARY S3: PER-SITE DIAGNOSTICS\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("=== PART A: BIOMASS BY SITE (12 days with AMMOD biomass) ===\n\n")

# --- Helpers (identical to R/00_trajectory_functions.R) ------------------------
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

out_a <- file.path(FIG_DIR, paste0(SCRIPT, "_biomass_dualaxis.png"))
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

out_b <- file.path(FIG_DIR, paste0(SCRIPT, "_biomass_relative.png"))
print(p_b)
ggplot2::ggsave(out_b, p_b, width = 9, height = 12, dpi = 200, bg = "white")
cat(sprintf("[OK] Saved: %s\n", out_b))


# ==============================================================================
# PART B: ABUNDANCE OF THE CAMERA SYSTEMS BY SITE
# ==============================================================================

# FAIRD and ID daily abundance, all days. Days of the ID power-outage window
# (Field_Valid_Operational_ID FALSE for ID) are flagged as partial effort; the
# same dates are flagged for FAIRD, so that any exclusion is symmetric.
outage_dates <- data_list$effort_base %>%
  dplyr::filter(Field_Valid, !Field_Valid_Operational_ID) %>%
  dplyr::distinct(Date) %>%
  dplyr::pull(Date) %>%
  as.Date()

abund <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, abundance) %>%
  dplyr::mutate(Date = as.Date(Date),
                Effort = ifelse(Date %in% outage_dates, "outage window", "full effort"))

# Sites where FAIRD and ID coexist.
sites_both <- abund %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(n_types = dplyr::n_distinct(Device_type), .groups = "drop") %>%
  dplyr::filter(n_types == 2) %>%
  dplyr::pull(Site)

cat("\n\n=== PART B: ABUNDANCE OF THE CAMERA SYSTEMS BY SITE ===\n\n")
cat(sprintf("Sites with both FAIRD and ID: %s\n", paste(sort(sites_both), collapse = ", ")))
cat(sprintf("Outage window (flagged for both systems): %s\n",
            paste(format(sort(outage_dates)), collapse = ", ")))

abund_both <- abund %>% dplyr::filter(Site %in% sites_both)

# --- Peak per unit, all days and full-effort days only ------------------------
cat("\n=== PEAK AND MEAN PER UNIT ===\n\n")
peak_tbl <- dplyr::bind_rows(
  abund_both %>% dplyr::mutate(Days = "all days"),
  abund_both %>% dplyr::filter(Effort == "full effort") %>% dplyr::mutate(Days = "full effort")
) %>%
  dplyr::group_by(Days, Site, Device_type, Device) %>%
  dplyr::summarise(
    n_days     = dplyr::n(),
    peak_date  = Date[which.max(abundance)],
    peak_value = max(abundance),
    mean_abund = round(mean(abundance), 1),
    .groups = "drop") %>%
  dplyr::arrange(Days, Site, Device_type)
print(peak_tbl, n = Inf)

# --- Mean in the central window of the season ---------------------------------
cat(sprintf("\n=== MEAN ABUNDANCE IN THE CENTRAL WINDOW (%s to %s), ALL DAYS ===\n\n",
            CENTRAL_START, CENTRAL_END))
central_tbl <- abund_both %>%
  dplyr::filter(Date >= CENTRAL_START, Date <= CENTRAL_END) %>%
  dplyr::group_by(Site, Device_type, Device) %>%
  dplyr::summarise(n_days         = dplyr::n(),
                   n_outage_days  = sum(Effort == "outage window"),
                   mean_central   = round(mean(abundance), 1),
                   max_central    = max(abundance),
                   .groups = "drop") %>%
  dplyr::arrange(Site, Device_type)
print(central_tbl, n = Inf)
cat("  Outage-window days are days of partial, unknown ID recording effort.\n")

# --- One panel per site ---------------------------------------------------------
# Points: all days (open symbols inside the outage window). Smooth: cubic
# spline fitted to the full-effort days of both systems.
abund_cap <- nice_ceiling(stats::quantile(abund_both$abundance, 0.98, na.rm = TRUE))

p <- ggplot2::ggplot(abund_both,
                     ggplot2::aes(x = Date, y = abundance, color = Device_type)) +
  ggplot2::geom_point(ggplot2::aes(shape = Effort), alpha = 0.7, size = 1.8) +
  ggplot2::geom_smooth(
    data = abund_both %>% dplyr::filter(Effort == "full effort"),
    ggplot2::aes(fill = Device_type),
    method = "lm", formula = y ~ splines::bs(x, 3),
    se = TRUE, alpha = 0.15, linewidth = 1) +
  ggplot2::facet_wrap(~ Site, ncol = 1, scales = "fixed") +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_fill_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = c("full effort" = 16, "outage window" = 1)) +
  ggplot2::coord_cartesian(ylim = c(0, abund_cap)) +
  ggplot2::labs(
    title    = "FAIRD and ID daily abundance by site",
    subtitle = "Sites with both camera systems; spline over full-effort days of both systems",
    x = "Date", y = "Abundance (individuals per day)",
    color = "System", fill = "System", shape = NULL,
    caption = sprintf("View capped at %s (rounded 98th percentile).", abund_cap)) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
    strip.text       = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top")

out_file <- file.path(FIG_DIR, paste0(SCRIPT, "_abundance_by_site.png"))
print(p)
ggplot2::ggsave(out_file, p, width = 9, height = 10, dpi = 200, bg = "white")
cat(sprintf("\n[OK] Saved: %s\n", out_file))


# ==============================================================================
# PART C: ID DETECTION CLASSES ON OUTAGE-WINDOW DAYS
# ==============================================================================

cat("\n\n=== PART C: ID DETECTION CLASSES ON OUTAGE-WINDOW DAYS ===\n\n")

id_class <- data_list$id %>%
  dplyr::mutate(Date = as.Date(Date), Class_ID = `Pre-classification`)

top_classes <- id_class %>%
  dplyr::count(Class_ID, sort = TRUE) %>%
  dplyr::slice_head(n = N_TOP_CLASSES) %>%
  dplyr::pull(Class_ID)
class_levels <- c(top_classes, "other classes")
cat(sprintf("Most frequent ID classes (all units, all days): %s; the rest pooled as 'other classes'.\n",
            paste(top_classes, collapse = ", ")))

class_counts <- id_class %>%
  dplyr::mutate(Class_ID = ifelse(Class_ID %in% top_classes, Class_ID, "other classes")) %>%
  dplyr::count(Device, Date, Class_ID, name = "n")

# Unit-day grid from the daily abundance table, so that zero days are kept.
id_grid <- abund %>%
  dplyr::filter(Device_type == "ID") %>%
  dplyr::distinct(Device, Date, Effort)

id_class_daily <- id_grid %>%
  tidyr::crossing(Class_ID = class_levels) %>%
  dplyr::left_join(class_counts, by = c("Device", "Date", "Class_ID")) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0L),
                Class_ID = factor(Class_ID, levels = class_levels)) %>%
  dplyr::arrange(Device, Date, Class_ID)

# Check: class counts add up to the daily abundance of each unit.
check_tot <- id_class_daily %>%
  dplyr::group_by(Device, Date) %>%
  dplyr::summarise(n_classes = sum(n), .groups = "drop") %>%
  dplyr::left_join(abund %>% dplyr::select(Device, Date, abundance), by = c("Device", "Date"))
n_mismatch <- sum(check_tot$n_classes != check_tot$abundance, na.rm = TRUE)
cat(sprintf("Check: class counts equal daily abundance on %d of %d unit-days %s\n\n",
            nrow(check_tot) - n_mismatch, nrow(check_tot),
            ifelse(n_mismatch == 0, "[OK]", "[MISMATCH]")))

# --- Table C1: counts per class, outage-window days vs full-effort days ---------
cat("--- Table C1: detections per class and day ---\n")
cat("full_* = over full-effort days; day_MMDD = count on that outage-window day;\n")
cat("n_above_full_max = outage-window days on which the count exceeds full_max.\n\n")

full_stats <- id_class_daily %>%
  dplyr::filter(Effort == "full effort") %>%
  dplyr::group_by(Device, Class_ID) %>%
  dplyr::summarise(full_n_days = dplyr::n(),
                   full_mean   = round(mean(n), 1),
                   full_median = stats::median(n),
                   full_max    = max(n),
                   .groups = "drop")

outage_long <- id_class_daily %>%
  dplyr::filter(Effort == "outage window") %>%
  dplyr::left_join(full_stats %>% dplyr::select(Device, Class_ID, full_max),
                   by = c("Device", "Class_ID"))

outage_above <- outage_long %>%
  dplyr::group_by(Device, Class_ID) %>%
  dplyr::summarise(n_above_full_max = sum(n > full_max), .groups = "drop")

outage_wide <- outage_long %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Day = paste0("day_", format(Date, "%m%d"))) %>%
  dplyr::select(Device, Class_ID, Day, n) %>%
  tidyr::pivot_wider(names_from = Day, values_from = n)

table_c1 <- full_stats %>%
  dplyr::left_join(outage_wide,  by = c("Device", "Class_ID")) %>%
  dplyr::left_join(outage_above, by = c("Device", "Class_ID")) %>%
  dplyr::arrange(Device, Class_ID)
print(as.data.frame(table_c1), row.names = FALSE)

# --- Table C2: class shares per day ---------------------------------------------
cat("\n--- Table C2: class shares (%) per outage-window day and over full-effort days ---\n\n")
table_c2 <- dplyr::bind_rows(
  id_class_daily %>%
    dplyr::filter(Effort == "outage window") %>%
    dplyr::mutate(Period = format(Date, "%Y-%m-%d")),
  id_class_daily %>%
    dplyr::filter(Effort == "full effort") %>%
    dplyr::mutate(Period = "full-effort days, pooled")
) %>%
  dplyr::group_by(Device, Period, Class_ID) %>%
  dplyr::summarise(n = sum(n), .groups = "drop") %>%
  dplyr::group_by(Device, Period) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(share_pct = ifelse(total > 0, round(100 * n / total, 1), NA_real_)) %>%
  dplyr::select(Device, Period, total, Class_ID, share_pct) %>%
  tidyr::pivot_wider(names_from = Class_ID, values_from = share_pct, names_prefix = "pct_") %>%
  dplyr::arrange(Device, Period)
print(as.data.frame(table_c2), row.names = FALSE)

readr::write_csv(table_c1, file.path(TAB_DIR, paste0(SCRIPT, "_id_class_counts.csv")))
readr::write_csv(table_c2, file.path(TAB_DIR, paste0(SCRIPT, "_id_class_shares.csv")))
cat(sprintf("\n[OK] Tables written to %s\n", TAB_DIR))


cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
while (sink.number() > 0) sink()
