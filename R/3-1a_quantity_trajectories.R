# ==============================================================================
# 3-1a_quantity_trajectories.R
# Manuscript Section 3.1 (Quantity): biomass trajectories of the three systems
# ------------------------------------------------------------------------------
# Daily biomass of AMMOD (gravimetric live mass), FAIRD and ID (allometric
# live mass) over the 12 days with AMMOD biomass (Aug 23 to Sep 3), and of
# FAIRD and ID over 22 days.
#
# Parts:
#   1. Figure model and Figure 3:
#        biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site / Device),
#        Tweedie, log link. Figure 3 shows its predictions on linear axes,
#        FAIRD and ID on the left axis, AMMOD on the right axis divided by a
#        data-derived factor k (sum AMMOD / sum FAIRD over co-located
#        site-days, 2 significant digits).
#   2. Test model: the figure model plus Block (BLOCK1 = Maize, Sites 1-2;
#        BLOCK2 = Meadow, Sites 3-4) as a fixed term and AR1 within unit.
#        Block is fixed because a two-level random effect collapses to a
#        singular fit (robustness check in S1). Cascade AR1 -> no AR1 -> no
#        Block, each level gated by tmb_is_valid(). Type III Wald test of the
#        Date x Device_type interaction for the three systems (6 df).
#   3. Pair models: the test model on each pair of systems (3 df).
#        Pairing PAIRED_BLOCK: block-days on which both systems have a value
#        in the same block. A_12d for the three pairs (primary family, Holm
#        over the three LRT p-values); FAIRD vs ID also over 22 days, all days
#        and with the ID power-outage window removed for both systems
#        (complementary). Daily ratio of unit means per pair (drift).
#   4. Internal parallelism: the same pair model on two units of one system.
#        FAIRD halves over both blocks ({1,3} vs {2,4}; {1,4} vs {2,3}); Maize
#        unit pairs for FAIRD, AMMOD and ID (Supplementary Material S2).
#   5. Figures: Figure 3, one figure per pair model, one per internal model.
#
# Under the log link a non-significant interaction means a constant ratio
# between trajectories; a non-significant test is not evidence of it.
# Sensitivity analyses (other pairings, identity link, Block random, Gaussian
# cross-check) are in supplementary/S1_pair_model_sensitivity.R.
# ==============================================================================

required_packages <- c("tidyverse", "here", "glmmTMB", "car", "effects",
                       "splines", "RColorBrewer", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)
library(glmmTMB)
library(splines)
library(effects)

source(here::here("R", "00_load_data.R"))
source(here::here("R", "00_custom_colors.R"))
source(here::here("R", "00_trajectory_functions.R"))

SCRIPT  <- "3-1a_quantity_trajectories"
FIG_DIR <- here::here("outputs", "figures")
TAB_DIR <- here::here("outputs", "tables")
CON_DIR <- here::here("outputs", "console")
for (d in c(FIG_DIR, TAB_DIR, CON_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, paste0(SCRIPT, "_output.txt")), split = TRUE)

cat("==============================================================================\n")
cat("SECTION 3.1, QUANTITY: BIOMASS TRAJECTORIES OF THE THREE SYSTEMS\n")
cat("Figure model: bs(Date,3) * Device_type + (1|Site/Device), Tweedie\n")
cat("Test model:   figure model + Block (fixed) + AR1 within unit\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: DATA
# ==============================================================================

cat("=== DATA ===\n\n")

# FAIRD and ID daily biomass (zero-capture days completed as real zeros in
# 00_load_data.R); AMMOD daily live mass.
etraps_bio <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, biomass_mg) %>%
  dplyr::rename(Biomass_mg = biomass_mg)

ammod_bio <- data_list$ammod_biomass_for_stats %>%
  dplyr::mutate(Device_type = "AMMOD", Biomass_mg = Live_mass) %>%
  dplyr::select(Date, Device, Device_type, Site, Biomass_mg)

bio_combined <- dplyr::bind_rows(etraps_bio, ammod_bio) %>%
  dplyr::arrange(Device_type, Site, Date) %>%
  dplyr::mutate(
    Date.num    = as.numeric(Date),
    Device_type = factor(Device_type, levels = c("AMMOD", "FAIRD", "ID")),
    Site        = factor(Site),
    Device      = factor(Device)
  )

# 12-day window with AMMOD biomass.
date_range  <- range(ammod_bio$Date)
bio_overlap <- bio_combined %>%
  dplyr::filter(Date >= date_range[1] & Date <= date_range[2])
bio_overlap$Block <- factor(
  ifelse(bio_overlap$Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2"))
bio_overlap$times <- factor(bio_overlap$Date)

# 22-day series of the camera systems.
bio_22 <- bio_combined %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::mutate(
    Block = factor(ifelse(Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")),
    times = factor(Date)
  )

# ID power-outage window, from the validity calendar (removed for both
# systems of a comparison).
adverse_dates <- data_list$effort_base %>%
  dplyr::filter(Field_Valid, !Field_Valid_Operational_ID) %>%
  dplyr::distinct(Date) %>%
  dplyr::arrange(Date) %>%
  dplyr::pull(Date)

time_levels_12 <- levels(bio_overlap$times)
time_levels_22 <- levels(bio_22$times)

cat(sprintf("12-day window: %s to %s, %d observations, %d dates\n",
            date_range[1], date_range[2], nrow(bio_overlap),
            dplyr::n_distinct(bio_overlap$Date)))
cat(sprintf("22-day camera series: %d observations\n", nrow(bio_22)))
cat(sprintf("ID power-outage window (%d dates): %s\n", length(adverse_dates),
            paste(format(adverse_dates), collapse = ", ")))

cat("\nObservations per system (12-day window):\n")
print(bio_overlap %>%
        dplyr::group_by(Device_type) %>%
        dplyr::summarise(n = dplyr::n(), n_units = dplyr::n_distinct(Device),
                         mean_biomass_mg = round(mean(Biomass_mg), 1), .groups = "drop"))
cat("\nObservations per unit (12-day window):\n")
print(bio_overlap %>% dplyr::count(Block, Device_type, Device, Site, name = "n"), n = Inf)


# ==============================================================================
# SECTION 2: FIGURE MODEL AND SCALE FACTOR k
# ==============================================================================

cat("\n\n=== FIGURE MODEL (Figure 3) ===\n")

m_fig <- glmmTMB::glmmTMB(Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site / Device),
                          family = glmmTMB::tweedie(), data = bio_overlap)
cat(sprintf("Valid fit (pdHess and finite AIC): %s\n", tmb_is_valid(m_fig)))
print(summary(m_fig))

# Scale factor k: sum of AMMOD over sum of FAIRD biomass, co-located
# site-days only, 2 significant digits.
ratio_info <- ammod_faird_ratio(bio_overlap)
k_ammod    <- signif(ratio_info$overall, 2)
cat(sprintf("\nDual-axis scale factor (sum AMMOD / sum FAIRD, %d co-located site-days): %s\n",
            ratio_info$n_site_days, k_ammod))


# ==============================================================================
# SECTION 3: TEST MODEL, THREE SYSTEMS (6-df interaction)
# ==============================================================================

cat("\n\n=== TEST MODEL, THREE SYSTEMS ===\n")

B1 <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + Block +
    (1 | Site / Device) + ar1(times + 0 | Device),
  bio_overlap, "Test model B1: Block fixed + (1|Site/Device) + AR1")
B0 <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + Block + (1 | Site / Device),
  bio_overlap, "Model B0: Block fixed + (1|Site/Device), no AR1")
B_simple <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site / Device),
  bio_overlap, "Model B_simple: (1|Site/Device), no Block, no AR1")

cat("\n--- AIC comparison (valid fits only) ---\n")
tmb_models <- list(B1 = B1, B0 = B0, B_simple = B_simple)
tmb_models <- tmb_models[vapply(tmb_models, tmb_is_valid, logical(1))]
if (length(tmb_models) >= 1) {
  aic_tbl <- data.frame(
    Model = names(tmb_models),
    df    = vapply(tmb_models, function(m) attr(logLik(m), "df"), numeric(1)),
    AIC   = vapply(tmb_models, AIC, numeric(1)),
    row.names = NULL)
  aic_tbl <- aic_tbl[order(aic_tbl$AIC), ]
  aic_tbl$dAIC <- aic_tbl$AIC - min(aic_tbl$AIC)
  print(aic_tbl, row.names = FALSE)
}

if (tmb_is_valid(B1)) {
  m_test <- B1; m_test_name <- "B1 (Block fixed, Site/Device, AR1)"
} else if (tmb_is_valid(B0)) {
  m_test <- B0; m_test_name <- "B0 (Block fixed, Site/Device, no AR1)"
  cat("[WARNING] B1 not valid -> test model B0.\n")
} else if (tmb_is_valid(B_simple)) {
  m_test <- B_simple; m_test_name <- "B_simple (Site/Device, no Block, no AR1)"
  cat("[WARNING] B1 and B0 not valid -> test model B_simple.\n")
} else {
  m_test <- NULL; m_test_name <- NA_character_
  cat("[ERROR] no valid test model.\n")
}
cat(sprintf("\nTest model used: %s\n", ifelse(is.na(m_test_name), "NONE", m_test_name)))

p_interact <- NA_real_
anova_test <- NULL
if (!is.null(m_test)) {
  cat("\n--- Type III ANOVA, test model ---\n")
  anova_test <- tryCatch(car::Anova(m_test, type = "III"), error = function(e) {
    cat("  [SKIP] Type III ANOVA failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(anova_test)) {
    print(anova_test)
    if (interaction_term %in% rownames(anova_test)) {
      p_interact <- anova_test[interaction_term, "Pr(>Chisq)"]
      cat(sprintf("\n  Date x Device_type interaction (6 df): p = %.4f\n", p_interact))
    }
  }
  cat("\n--- Fixed-effect coefficients, test model ---\n")
  print(glmmTMB::fixef(m_test)$cond)
  cat("\n--- Variance components, test model ---\n")
  print(glmmTMB::VarCorr(m_test))
}


# ==============================================================================
# SECTION 4: PAIR MODELS (3-df interaction), PAIRED_BLOCK
# ==============================================================================

cat("\n\n=== PAIR MODELS: AMMOD-FAIRD, AMMOD-ID, FAIRD-ID ===\n")

pair_types <- list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))
pair_runs  <- list()

for (tp in pair_types) {
  r <- run_pair(tp, bio_overlap, time_levels_12, "Block", "A_12d", "full", "primary")
  pair_runs[[r$label]] <- r
}
r <- run_pair(c("FAIRD", "ID"), bio_22, time_levels_22, "Block",
              "B_22d", "full", "complementary")
pair_runs[[r$label]] <- r
r <- run_pair(c("FAIRD", "ID"), bio_22, time_levels_22, "Block",
              "B_22d", "window_excluded", "complementary", drop_dates = adverse_dates)
pair_runs[[r$label]] <- r

pair_tbl <- dplyr::bind_rows(lapply(pair_runs, `[[`, "row"))

is_family <- pair_tbl$Analysis == "A_12d" & pair_tbl$Role == "primary"
n_family  <- sum(is_family)
if (anyNA(pair_tbl$p_interaction_LRT[is_family])) {
  cat("\n[WARNING] at least one primary pair has no LRT p-value; Holm keeps\n")
  cat("          n = 3 comparisons so the correction is not weakened.\n")
}
pair_tbl$p_LRT_Holm  <- NA_real_
pair_tbl$p_Wald_Holm <- NA_real_
pair_tbl$p_LRT_Holm[is_family] <- stats::p.adjust(
  pair_tbl$p_interaction_LRT[is_family], method = "holm", n = n_family)
pair_tbl$p_Wald_Holm[is_family] <- stats::p.adjust(
  pair_tbl$p_interaction_Wald[is_family], method = "holm", n = n_family)

omnibus_row <- data.frame(
  Pair = "ALL (3 systems)", Analysis = "A_12d", Window = "full", Pairing = "FULL",
  Role = "omnibus", Variant = "OMNIBUS (3 systems)",
  n_obs = nrow(bio_overlap), n_days = dplyr::n_distinct(bio_overlap$Date),
  n_devices = dplyr::n_distinct(bio_overlap$Device), Model = m_test_name,
  Chisq_interaction = if (!is.null(anova_test) && interaction_term %in% rownames(anova_test))
    anova_test[interaction_term, "Chisq"] else NA_real_,
  df_interaction = if (!is.null(anova_test) && interaction_term %in% rownames(anova_test))
    anova_test[interaction_term, "Df"] else NA_real_,
  p_interaction_Wald = p_interact,
  stringsAsFactors = FALSE)
pair_tbl <- dplyr::bind_rows(omnibus_row, pair_tbl)

cat("\n------------------------------------------------------------------------\n")
cat("DATE x DEVICE_TYPE INTERACTION: THREE SYSTEMS AND PAIR MODELS\n")
cat("------------------------------------------------------------------------\n")
print(pair_tbl %>%
        dplyr::select(Pair, Analysis, Window, Pairing, Role, n_obs, n_days,
                      n_devices, Model, df_interaction, p_interaction_Wald,
                      p_interaction_LRT, LRT_level, p_Wald_at_LRT_level,
                      p_LRT_Holm, p_Wald_Holm) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with("p_"), ~ round(.x, 4))),
      row.names = FALSE)
cat("\n  NOTES:\n")
cat("  - Primary family = the three A_12d PAIRED_BLOCK rows; Holm on the LRT is\n")
cat("    the declared correction, Holm on Wald is reported alongside.\n")
cat("  - Compare p_interaction_LRT with p_Wald_at_LRT_level (same cascade level).\n")
cat("  - B_22d rows are complementary and outside the Holm family.\n")

readr::write_csv(pair_tbl, file.path(TAB_DIR, paste0(SCRIPT, "_pair_models.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_pair_models.csv")), "\n")

cat("\n--- Daily ratio of unit means per pair (PAIRED_BLOCK data) ---\n")
for (r in pair_runs) {
  rr <- system_ratio(r$data, r$types[1], r$types[2])
  cat(sprintf("\n  %s: overall ratio %s/%s = %.2f, Spearman rho(ratio, date) = %.2f\n",
              r$label, r$types[1], r$types[2], rr$overall,
              suppressWarnings(stats::cor(as.numeric(rr$daily$Date), rr$daily$ratio,
                                          method = "spearman", use = "complete.obs"))))
  print(rr$daily %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 2))),
        n = Inf)
}


# ==============================================================================
# SECTION 5: INTERNAL PARALLELISM (two units of one system)
# ==============================================================================

cat("\n\n=== INTERNAL CONSISTENCY: PARALLELISM WITHIN EACH SYSTEM ===\n")

faird_splits <- list(
  S1 = list("FAIRD1+FAIRD3" = c("FAIRD1", "FAIRD3"), "FAIRD2+FAIRD4" = c("FAIRD2", "FAIRD4")),
  S2 = list("FAIRD1+FAIRD4" = c("FAIRD1", "FAIRD4"), "FAIRD2+FAIRD3" = c("FAIRD2", "FAIRD3"))
)
maize_units <- function(s) stats::setNames(list(paste0(s, "1"), paste0(s, "2")),
                                           c(paste0(s, "1"), paste0(s, "2")))

internal_windows <- list(
  list(dat = bio_overlap, levels = time_levels_12, analysis = "A_12d", window = "full",
       drop = NULL, systems = c("FAIRD", "AMMOD", "ID")),
  list(dat = bio_22, levels = time_levels_22, analysis = "B_22d", window = "full",
       drop = NULL, systems = c("FAIRD", "ID")),
  list(dat = bio_22, levels = time_levels_22, analysis = "B_22d", window = "window_excluded",
       drop = adverse_dates, systems = c("FAIRD", "ID"))
)

internal_runs <- list()
for (iw in internal_windows) {
  for (sp in names(faird_splits)) {
    r <- run_internal("FAIRD", faird_splits[[sp]], iw$dat, iw$levels, iw$analysis,
                      iw$window, "two_block", sp, drop_dates = iw$drop)
    internal_runs[[r$label]] <- r
  }
  for (s in iw$systems) {
    r <- run_internal(s, maize_units(s), iw$dat, iw$levels, iw$analysis, iw$window,
                      "maize_only", "maize", maize_only = TRUE, drop_dates = iw$drop)
    internal_runs[[r$label]] <- r
  }
}

internal_tbl <- dplyr::bind_rows(lapply(internal_runs, `[[`, "row"))

cat("\n------------------------------------------------------------------------\n")
cat("INTERNAL CONSISTENCY: DATE x UNIT INTERACTION WITHIN EACH SYSTEM\n")
cat("------------------------------------------------------------------------\n")
print(internal_tbl %>%
        dplyr::select(System, Comparison, Analysis, Window, Version, n_obs, n_days,
                      Model, df_interaction, p_interaction_Wald, p_interaction_LRT,
                      LRT_level, p_Wald_at_LRT_level) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with("p_"), ~ round(.x, 4))),
      row.names = FALSE)
cat("\n  NOTES:\n")
cat("  - Same model and LRT as the pair models; constant ratio (log link) only.\n")
cat("  - two_block rows: FAIRD halves, one unit per block in each half (main text).\n")
cat("  - maize_only rows compare single units (Supplementary Material S2); they\n")
cat("    may use the last-resort levels AR1_noRE or fixed (see LRT_level), and a\n")
cat("    p-value at level fixed is less conservative.\n")
cat("  - Unadjusted p-values: references for the pair models, not a family.\n")

readr::write_csv(internal_tbl, file.path(TAB_DIR, paste0(SCRIPT, "_internal_parallelism.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_internal_parallelism.csv")), "\n")


# ==============================================================================
# SECTION 6: FIGURES (linear axes)
# ==============================================================================

cat("\n\n=== FIGURES ===\n\n")

# --- Figure 3: figure model, dual axis ----------------------------------------
# Predictions from effects::Effect() on a fine date grid, response scale (mg).
date_grid <- seq(min(bio_overlap$Date.num), max(bio_overlap$Date.num), length.out = 100)
eff_fig   <- effects::Effect(c("Device_type", "Date.num"), m_fig,
                             xlevels = list(Date.num = date_grid))
pred_fig  <- as.data.frame(eff_fig)
pred_fig$Date <- as.Date(pred_fig$Date.num, origin = "1970-01-01")

is_ammod_pred <- pred_fig$Device_type == "AMMOD"
pred_fig$fit_plot   <- ifelse(is_ammod_pred, pred_fig$fit   / k_ammod, pred_fig$fit)
pred_fig$lower_plot <- ifelse(is_ammod_pred, pred_fig$lower / k_ammod, pred_fig$lower)
pred_fig$upper_plot <- ifelse(is_ammod_pred, pred_fig$upper / k_ammod, pred_fig$upper)

obs_fig <- bio_overlap %>%
  dplyr::mutate(y_plot = dplyr::if_else(Device_type == "AMMOD", Biomass_mg / k_ammod, Biomass_mg))

# View capped at the rounded 98th percentile of the plotted points
# (coord_cartesian: nothing is dropped from lines or ribbons).
q98   <- stats::quantile(obs_fig$y_plot, 0.98, na.rm = TRUE)
y_cap <- ceiling(q98 / 10^floor(log10(q98))) * 10^floor(log10(q98))

system_colors <- stats::setNames(RColorBrewer::brewer.pal(n = 3, name = "Set1"),
                                 c("AMMOD", "FAIRD", "ID"))

p_fig3 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(data = pred_fig,
                       ggplot2::aes(x = Date, ymin = lower_plot, ymax = upper_plot,
                                    fill = Device_type), alpha = 0.15) +
  ggplot2::geom_line(data = pred_fig,
                     ggplot2::aes(x = Date, y = fit_plot, color = Device_type), linewidth = 1) +
  ggplot2::geom_point(data = obs_fig,
                      ggplot2::aes(x = Date, y = y_plot, color = Device_type),
                      alpha = 0.5, size = 1.8) +
  ggplot2::scale_color_manual(values = system_colors) +
  ggplot2::scale_fill_manual(values = system_colors) +
  ggplot2::scale_y_continuous(
    name     = "FAIRD / ID biomass (mg)",
    labels   = scales::label_comma(),
    sec.axis = ggplot2::sec_axis(~ . * k_ammod,
                                 name   = sprintf("AMMOD biomass (mg), axis = left x %s", k_ammod),
                                 labels = scales::label_comma())) +
  ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
  ggplot2::labs(x = "Date", color = "System", fill = "System",
                caption = sprintf("View capped at %s (rounded 98th percentile of plotted points).",
                                  scales::label_comma()(y_cap))) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position  = "top",
                 plot.caption     = ggplot2::element_text(size = 8, color = "grey50"))

print(p_fig3)
f_fig3 <- file.path(FIG_DIR, paste0(SCRIPT, "_fig3_biomass_dualaxis.png"))
ggplot2::ggsave(f_fig3, p_fig3, width = 10, height = 7, dpi = 300, bg = "white")
cat("[OK] Saved:", f_fig3, "\n")

# --- One figure per pair model --------------------------------------------------
# Pairs with AMMOD on a dual axis whose factor is the ratio of unit means over
# the pair data; FAIRD vs ID on a single linear axis.
for (r in pair_runs) {
  if (is.null(r$res$model)) {
    cat(sprintf("\n[SKIP] Pair figure %s: no valid pair model.\n", r$label))
    next
  }
  has_ammod <- "AMMOD" %in% r$types
  other     <- setdiff(r$types, "AMMOD")[1]
  k_pair    <- if (has_ammod) signif(system_ratio(r$data, "AMMOD", other)$overall, 2) else 1
  row       <- r$row
  ref_block <- levels(r$data$Block)[1]
  pred_pair <- predict_fixed(r$res$model, r$data, ref_block)
  window_note <- if (r$window == "window_excluded") {
    sprintf(" | excluded: %s to %s", format(min(adverse_dates)), format(max(adverse_dates)))
  } else ""
  axis_note <- if (has_ammod) sprintf(" | AMMOD right axis, scale factor %s (data-derived)", k_pair) else ""
  p_pair <- make_dual_plot(
    pred      = pred_pair,
    obs       = r$data,
    k         = k_pair,
    title     = sprintf("Biomass trajectories: %s (pair model, %s, %s)",
                        sub("-", " vs ", row$Pair), row$Analysis, row$Window),
    subtitle  = sprintf(
      "%s | %s | predictions at %s | Date x Device_type p = %.4f (Wald, %s df), %.4f (LRT, level %s)%s%s",
      r$res$spec$name, row$Pairing, ref_block, row$p_interaction_Wald, row$df_interaction,
      row$p_interaction_LRT, row$LRT_level, window_note, axis_note),
    left_name = if (has_ammod) sprintf("%s Biomass (mg)", other) else "Biomass (mg)",
    dual      = has_ammod)
  print(p_pair)
  fig_file <- file.path(FIG_DIR, sprintf("%s_pair_%s_%s_%s.png", SCRIPT,
                                         gsub("-", "_", row$Pair), row$Analysis, row$Window))
  ggplot2::ggsave(fig_file, p_pair, width = 10, height = 7, dpi = 300, bg = "white")
  cat(sprintf("[OK] Saved: %s\n", fig_file))
}

# --- One figure per internal-consistency model ----------------------------------
for (r in internal_runs) {
  if (is.null(r$res$model)) {
    cat(sprintf("\n[SKIP] Internal figure %s: no valid model.\n", r$label))
    next
  }
  row       <- r$row
  ref_block <- levels(r$data$Block)[1]
  pred_int  <- predict_fixed(r$res$model, r$data, ref_block)
  window_note <- if (r$window == "window_excluded") {
    sprintf(" | excluded: %s to %s", format(min(adverse_dates)), format(max(adverse_dates)))
  } else ""
  p_int <- make_internal_plot(
    pred     = pred_int,
    obs      = r$data,
    units    = r$units,
    system   = r$system,
    title    = sprintf("Internal consistency, %s: %s (%s, %s)", r$system, row$Comparison,
                       row$Analysis, row$Window),
    subtitle = sprintf(
      "%s | %s | predictions at %s | Date x Unit p = %.4f (Wald, %s df), %.4f (LRT, level %s)%s",
      r$res$spec$name, r$version, ref_block, row$p_interaction_Wald, row$df_interaction,
      row$p_interaction_LRT, row$LRT_level, window_note))
  print(p_int)
  fig_file <- file.path(FIG_DIR, sprintf("%s_internal_%s_%s_%s_%s.png", SCRIPT, r$system,
                                         gsub("[^A-Za-z0-9]+", "_", row$Comparison),
                                         row$Analysis, row$Window))
  ggplot2::ggsave(fig_file, p_int, width = 10, height = 7, dpi = 300, bg = "white")
  cat(sprintf("[OK] Saved: %s\n", fig_file))
}


# ==============================================================================
# SECTION 7: EXPORT
# ==============================================================================

cat("\n=== EXPORT ===\n\n")

if (!is.null(anova_test)) {
  anova_df <- as.data.frame(anova_test) %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::mutate(Model = m_test_name)
  readr::write_csv(anova_df, file.path(TAB_DIR, paste0(SCRIPT, "_test_model_ANOVA_III.csv")))
  cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_test_model_ANOVA_III.csv")), "\n")
}

fig_coef <- as.data.frame(summary(m_fig)$coefficients$cond) %>%
  tibble::rownames_to_column("Term")
readr::write_csv(fig_coef, file.path(TAB_DIR, paste0(SCRIPT, "_figure_model_coefficients.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_figure_model_coefficients.csv")), "\n")

cat("\n==============================================================================\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

while (sink.number() > 0) sink()
