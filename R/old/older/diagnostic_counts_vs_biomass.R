# ==============================================================================
# DIAGNOSTIC: Why are correlations different?
# ==============================================================================

library(tidyverse)
load("data_cleaned.RData")

# Extract data
etraps_daily <- data_list$etraps_daily

# ==============================================================================
# COMPARISON 1: Abundance (counts) vs Biomass
# ==============================================================================

cat("\n=== COMPARISON 1: Counts vs Biomass ===\n\n")

# ----- COUNTS (Abundance) -----
# Prepare data for FAIRD Maize (primary analysis)
counts_data <- etraps_daily %>%
  dplyr::filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  dplyr::select(Date, Device, abundance) %>%
  tidyr::pivot_wider(names_from = Device, values_from = abundance, values_fill = 0)

cor_counts <- cor.test(counts_data$FAIRD1, counts_data$FAIRD2)
lm_counts <- lm(FAIRD2 ~ FAIRD1, data = counts_data)

cat("--- FAIRD1 vs FAIRD2 using COUNTS (abundance) ---\n")
cat("Sample size:", nrow(counts_data), "days\n")
cat("Pearson r =", round(cor_counts$estimate, 3), "\n")
cat("R² =", round(summary(lm_counts)$r.squared, 3), "\n")
cat("p-value =", format.pval(cor_counts$p.value, digits = 3), "\n\n")

# ----- BIOMASS -----
biomass_data <- etraps_daily %>%
  filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  select(Date, Device, biomass_mg) %>%
  pivot_wider(names_from = Device, values_from = biomass_mg, values_fill = 0)

cor_biomass <- cor.test(biomass_data$FAIRD1, biomass_data$FAIRD2)
lm_biomass <- lm(FAIRD2 ~ FAIRD1, data = biomass_data)

cat("--- FAIRD1 vs FAIRD2 using BIOMASS (mg) ---\n")
cat("Sample size:", nrow(biomass_data), "days\n")
cat("Pearson r =", round(cor_biomass$estimate, 3), "\n")
cat("R² =", round(summary(lm_biomass)$r.squared, 3), "\n")
cat("p-value =", format.pval(cor_biomass$p.value, digits = 3), "\n\n")

# ==============================================================================
# COMPARISON 2: Original vs Log-transformed
# ==============================================================================

cat("\n=== COMPARISON 2: Original vs Log-transformed ===\n\n")

# ----- ORIGINAL (no transform) -----
cat("--- BIOMASS - ORIGINAL scale ---\n")
cat("R² =", round(summary(lm_biomass)$r.squared, 3), "\n\n")

# ----- LOG-TRANSFORMED -----
biomass_log <- biomass_data %>%
  mutate(
    FAIRD1_log = log10(FAIRD1 + 1),
    FAIRD2_log = log10(FAIRD2 + 1)
  )

lm_biomass_log <- lm(FAIRD2_log ~ FAIRD1_log, data = biomass_log)

cat("--- BIOMASS - LOG10(x+1) transformed ---\n")
cat("R² =", round(summary(lm_biomass_log)$r.squared, 3), "\n\n")

# ==============================================================================
# COMPARISON 3: Different handling of zeros/NAs
# ==============================================================================

cat("\n=== COMPARISON 3: Handling of zeros/NAs ===\n\n")

# ----- With zeros (values_fill = 0) -----
biomass_with_zeros <- etraps_daily %>%
  filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  select(Date, Device, biomass_mg) %>%
  pivot_wider(names_from = Device, values_from = biomass_mg, values_fill = 0)

lm_with_zeros <- lm(FAIRD2 ~ FAIRD1, data = biomass_with_zeros)

cat("--- With explicit zeros (values_fill = 0) ---\n")
cat("Sample size:", nrow(biomass_with_zeros), "days\n")
cat("R² =", round(summary(lm_with_zeros)$r.squared, 3), "\n\n")

# ----- Without zeros (drop_na) -----
biomass_no_na <- etraps_daily %>%
  filter(Device %in% c("FAIRD1", "FAIRD2")) %>%
  select(Date, Device, biomass_mg) %>%
  pivot_wider(names_from = Device, values_from = biomass_mg) %>%
  drop_na()

lm_no_na <- lm(FAIRD2 ~ FAIRD1, data = biomass_no_na)

cat("--- Without NAs (drop_na) ---\n")
cat("Sample size:", nrow(biomass_no_na), "days\n")
cat("R² =", round(summary(lm_no_na)$r.squared, 3), "\n\n")

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

cat("\n=== SUMMARY: All Combinations ===\n\n")

summary_table <- tibble(
  Method = c(
    "Counts (abundance), original",
    "Counts (abundance), log10",
    "Biomass, original",
    "Biomass, log10",
    "Biomass, with zeros",
    "Biomass, drop NA"
  ),
  R_squared = c(
    summary(lm(FAIRD2 ~ FAIRD1, data = counts_data))$r.squared,
    summary(lm(log10(FAIRD2+1) ~ log10(FAIRD1+1), data = counts_data))$r.squared,
    summary(lm_biomass)$r.squared,
    summary(lm_biomass_log)$r.squared,
    summary(lm_with_zeros)$r.squared,
    summary(lm_no_na)$r.squared
  ),
  n_days = c(
    nrow(counts_data),
    nrow(counts_data),
    nrow(biomass_data),
    nrow(biomass_log),
    nrow(biomass_with_zeros),
    nrow(biomass_no_na)
  )
) %>%
  mutate(R_squared = round(R_squared, 3))

print(summary_table)

cat("\n✅ Diagnosis complete!\n")
cat("Compare these results with your alternative script.\n")