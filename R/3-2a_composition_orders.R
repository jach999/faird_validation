# ==============================================================================
# 3-2a_composition_orders.R
# Manuscript Section 3.2 (Composition): orders detected by each system
# ------------------------------------------------------------------------------
# Figure 4 and Table 2.
#
# Record unit: one record per order present per unit-day, for every system.
# AMMOD taxonomy (metabarcoding) is deduplicated in 00_load_data.R to one row
# per Date x Device x Order x Family x Genus; it is reduced here to one row per
# Date x Device x Order, the unit of the camera systems. Abundance enters for
# no system, so the model describes detection frequency, not the share of
# individuals.
#
# Model: nnet::multinom(Order ~ bs(Date.num, 3) * Device_type) on the six most
# frequent orders (by name, records not assigned to an order excluded); terms
# tested by likelihood-ratio tests (car::Anova, Type II). The model has no
# random effect for unit: the daily records of all units enter as independent
# observations, so the tests are descriptive.
#
# Table 2: number of unit-days on which each order was detected, per system.
# Denominators (valid taxonomic unit-days) come from 2-4_data_availability.R.
# ==============================================================================

required_packages <- c("tidyverse", "here", "nnet", "car", "effects", "splines",
                       "lattice", "plyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)
library(nnet)
library(splines)
library(effects)

source(here::here("R", "00_load_data.R"))

SCRIPT  <- "3-2a_composition_orders"
FIG_DIR <- here::here("outputs", "figures")
TAB_DIR <- here::here("outputs", "tables")
CON_DIR <- here::here("outputs", "console")
for (d in c(FIG_DIR, TAB_DIR, CON_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, paste0(SCRIPT, "_output.txt")), split = TRUE)

cat("==============================================================================\n")
cat("SECTION 3.2, COMPOSITION: ORDERS DETECTED BY EACH SYSTEM\n")
cat("Multinomial model: Order ~ bs(Date.num, 3) * Device_type\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: PRESENCE RECORDS, ONE PER ORDER PER UNIT-DAY
# ==============================================================================

ammod_presence_order <- ammod_presence_daily[
  !duplicated(ammod_presence_daily[, c("Date", "Device", "Order")]), ]

cat("AMMOD presence records, genus level:", nrow(ammod_presence_daily),
    "| order level:", nrow(ammod_presence_order), "\n")
cat("Camera presence records (order level):", nrow(etraps_presence_daily), "\n\n")

presence <- plyr::rbind.fill(etraps_presence_daily, ammod_presence_order)
presence$Order <- factor(presence$Order)
levels(presence$Order)[levels(presence$Order) == "#N/C"] <- "Indet."
presence$Date     <- as.Date(presence$Date)
presence$Date.num <- as.numeric(presence$Date)

order_table_all <- sort(table(presence$Order), decreasing = TRUE)
cat("Order frequency, all systems (unit-days with the order present):\n")
print(order_table_all)

cat("\nTable 2 source: unit-days with the order present, by system:\n")
order_by_system <- table(presence$Order, presence$Device_type)
print(order_by_system)


# ==============================================================================
# SECTION 2: MULTINOMIAL MODEL, SIX MOST FREQUENT ORDERS
# ==============================================================================

toporders <- head(setdiff(names(order_table_all), "Indet."), 6)
cat("\nSelected orders:\n")
print(toporders)

subdata <- subset(presence, Order %in% toporders)
subdata$Order <- droplevels(subdata$Order)

mmod1 <- nnet::multinom(Order ~ bs(Date.num, 3) * Device_type, data = subdata,
                        MaxNWts = 50000, maxit = 100000)

cat("\n--- Likelihood-ratio tests (car::Anova, Type II) ---\n")
anova_mmod <- car::Anova(mmod1)
print(anova_mmod)


# ==============================================================================
# SECTION 3: FIGURE 4
# ==============================================================================

# levels2dates() (date tick labels) is called unqualified, as in the original
# figure code; it must be available in the session (effects package).
eff_date <- effects::Effect("Date.num", mmod1)

f_fig4 <- file.path(FIG_DIR, paste0(SCRIPT, "_fig4_order_composition.pdf"))
grDevices::pdf(f_fig4, width = 10, height = 7)
print(
  plot(effects::allEffects(mmod1, xlevels = 100), style = "stacked",
       lattice = list(key.args = list(columns = 3)),
       axes = list(x = list(Date.num = list(
         lab   = "Date",
         ticks = list(at = levels2dates(eff_date, "Date.num", "1970-01-01")),
         rotate = 45))))
)
grDevices::dev.off()
cat("\n[OK] Saved:", f_fig4, "\n")


# ==============================================================================
# SECTION 4: EXPORT
# ==============================================================================

table2 <- as.data.frame.matrix(order_by_system) %>%
  tibble::rownames_to_column("Order") %>%
  dplyr::arrange(dplyr::desc(AMMOD))
readr::write_csv(table2, file.path(TAB_DIR, paste0(SCRIPT, "_table2_order_unit_days.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_table2_order_unit_days.csv")), "\n")

anova_df <- as.data.frame(anova_mmod) %>% tibble::rownames_to_column("Term")
readr::write_csv(anova_df, file.path(TAB_DIR, paste0(SCRIPT, "_multinomial_LR_tests.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_multinomial_LR_tests.csv")), "\n")

cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

while (sink.number() > 0) sink()
