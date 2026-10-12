###Quick diagnostic scripts###

# Check how many days survive different thresholds
# For H1 (FAIRD vs FAIRD at Superfamily level)
faird_daily_sfam <- faird_data %>%
  dplyr::mutate(
    Taxon_Sfam = dplyr::if_else(
      !is.na(Superfamily) & Superfamily != "#N/C", 
      Superfamily, NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Taxon_Sfam)) %>%
  dplyr::group_by(Date, Device) %>%
  dplyr::summarise(n_indiv = dplyr::n(), .groups = "drop")

# For each threshold, how many valid pair-days remain?
for (thresh in c(1, 2, 3, 5)) {
  cat(sprintf("\n=== Threshold: >= %d individuals ===\n", thresh))
  
  # FAIRD pairs for H1
  for (pair in list(c("FAIRD1","FAIRD2"), c("FAIRD3","FAIRD4"))) {
    valid_days <- faird_daily_sfam %>%
      dplyr::filter(Device %in% pair, n_indiv >= thresh) %>%
      dplyr::group_by(Date) %>%
      dplyr::summarise(n_devs = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n_devs == 2) %>%
      nrow()
    cat(sprintf("  %s vs %s: %d valid days\n", pair[1], pair[2], valid_days))
  }
}

# For H3 (FAIRD genera per day — how many days have >= threshold)
cat("\n=== H3: FAIRD genera/day at Genus level ===\n")
faird_daily_genus <- faird_data %>%
  dplyr::filter(!is.na(Genus) & Genus != "#N/C") %>%
  dplyr::group_by(Date, Device, Site) %>%
  dplyr::summarise(
    n_indiv = dplyr::n(),
    n_genera = length(unique(Genus)),
    .groups = "drop"
  )

for (thresh in c(1, 2, 3, 5)) {
  n_valid <- faird_daily_genus %>%
    dplyr::filter(n_indiv >= thresh) %>%
    nrow()
  cat(sprintf("  >= %d individuals: %d site-days (of %d total)\n",
              thresh, n_valid, nrow(faird_daily_genus)))
}


# Compare body length of AMMOD taxa shared vs exclusive to AMMOD
# at Superfamily level

source("load_data.R")
nc <- "#N/C"

# FAIRD Superfamilies (unique across all devices, all days)
faird_sf <- data_list$faird %>%
  dplyr::filter(!is.na(Superfamily) & Superfamily != nc) %>%
  dplyr::pull(Superfamily) %>% unique()

cat("FAIRD unique Superfamilies:", length(faird_sf), "\n")

# AMMOD taxa with body length — use data_list$ammod_taxonomy directly
ammod_with_bl <- data_list$ammod_taxonomy %>%
  dplyr::filter(!is.na(Superfamily) & Superfamily != nc,
                !is.na(Body_length)) %>%
  dplyr::mutate(
    Detection = dplyr::if_else(
      Superfamily %in% faird_sf, "Shared", "AMMOD-exclusive"
    )
  )

cat("AMMOD OTUs with Body_length:", nrow(ammod_with_bl), "\n\n")

# Summary
size_comparison <- ammod_with_bl %>%
  dplyr::group_by(Detection) %>%
  dplyr::summarise(
    n_otus = dplyr::n(),
    n_superfamilies = length(unique(Superfamily)),
    mean_BL = round(mean(Body_length, na.rm = TRUE), 2),
    sd_BL = round(sd(Body_length, na.rm = TRUE), 2),
    median_BL = round(median(Body_length, na.rm = TRUE), 2),
    pct_below_3mm = round(100 * mean(Body_length < 3, na.rm = TRUE), 1),
    .groups = "drop"
  )

print(size_comparison)

# Statistical test
shared_bl <- ammod_with_bl %>%
  dplyr::filter(Detection == "Shared") %>%
  dplyr::pull(Body_length)
exclusive_bl <- ammod_with_bl %>%
  dplyr::filter(Detection == "AMMOD-exclusive") %>%
  dplyr::pull(Body_length)

cat(sprintf("\nWilcoxon test: p = %.6f\n",
            wilcox.test(shared_bl, exclusive_bl)$p.value))
cat(sprintf("Cohen's d = %.3f\n",
            effsize::cohen.d(shared_bl, exclusive_bl)$estimate))