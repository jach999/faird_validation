# ==============================================================================
# diagnostic_yaxis_percentiles.R
# Quick, throwaway diagnostic: inspect the abundance distribution for FAIRD + ID
# to choose a common y-axis upper limit for the temporal plots.
# Not part of the analysis pipeline; run once, read the numbers, then discard.
# ==============================================================================

source("R/00_load_data.R")

# Same base data the main script uses: FAIRD + ID daily abundance.
abund <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::pull(abundance)

cat("=== ABUNDANCE DISTRIBUTION (FAIRD + ID, all 22 days) ===\n\n")
cat(sprintf("N observations : %d\n", length(abund)))
cat(sprintf("Min / Max      : %g / %g\n", min(abund), max(abund)))
cat(sprintf("Mean / Median  : %.2f / %g\n", mean(abund), stats::median(abund)))

# Candidate percentiles for the y-axis cap. The upper tail is what matters:
# a cap near p95-p99 keeps the curves readable while excluding only a handful
# of extreme points from the visible panel (the fit still uses every point).
probs <- c(0.90, 0.95, 0.975, 0.99, 1.00)
qs <- stats::quantile(abund, probs = probs, names = FALSE)

cat("\n=== UPPER-TAIL PERCENTILES (candidates for y-axis cap) ===\n\n")
for (i in seq_along(probs)) {
  cat(sprintf("  p%-5s : %6.1f\n", format(probs[i] * 100), qs[i]))
}

# How many points would fall outside a given cap (i.e. be clipped from view).
cat("\n=== POINTS ABOVE CANDIDATE CAPS (clipped from view, fit unaffected) ===\n\n")
for (cap in c(80, 100, 120, 150)) {
  n_above <- sum(abund > cap)
  cat(sprintf("  cap = %3d : %2d of %d points above (%.1f%%)\n",
              cap, n_above, length(abund), 100 * n_above / length(abund)))
}

# Same, split by Device_type, so we see whether the clipped points are mostly
# ID (expected, given the battery-window spike) or spread across both.
cat("\n=== POINTS ABOVE 100, BY DEVICE TYPE ===\n\n")
print(
  data_list$etraps_daily %>%
    dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
    dplyr::group_by(Device_type) %>%
    dplyr::summarise(
      n            = dplyr::n(),
      above_100    = sum(abundance > 100),
      max_abund    = max(abundance),
      .groups = "drop"
    )
)
