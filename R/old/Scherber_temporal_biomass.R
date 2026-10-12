# ==============================================================================
# Scherber_temporal_biomass.R — Temporal biomass analysis (mixed models)
# ------------------------------------------------------------------------------
# Reconstructs the in-office exploration of temporal biomass trajectories
# as an executable, defensive pipeline.
#
# Statistical intent (unchanged from Scherber edits):
#   - Block structure: Site1+Site2 = BLOCK1 (Maize), Site3+Site4 = BLOCK2 (Meadow)
#   - Random slope over time, nested Block/Site/Device
#   - Temporal autocorrelation (AR1) on consecutive days within device
#   - Model comparison by AIC (B0 without AR1 vs B1 with AR1)
#   - Variance-component inspection
#
# Defensive design:
#   - Every model is wrapped in tryCatch(); a non-converging model is logged
#     and the script continues.
#   - Plots run on `model_for_plots`, selected by cascade B1 -> B0 -> B_simple.
#   - A convergence banner reports which model feeds the figures.
#
# NOTE: package installation is NOT performed here. If a package is missing,
#       the script logs it and skips the dependent step. Install once, separately
#       (see the availability check in SECTION 0).
# ==============================================================================

library(tidyverse)
library(glmmTMB)
library(car)
library(splines)
library(ggplot2)
library(nlme)

source("R/00_load_data.R")
source("R/00_custom_colors.R")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console", showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/Scherber_temporal_biomass_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("SCHERBER TEMPORAL BIOMASS ANALYSIS (executable / defensive reconstruction)\n")
cat("Mixed models: bs(Date, 3) * Device_type + Block/Site/Device nesting + AR1\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 0: PACKAGE AVAILABILITY CHECK (no installation performed)
# ==============================================================================

cat("=== OPTIONAL PACKAGE AVAILABILITY ===\n\n")

# Optional packages used only for variance-component inspection / RE plots.
optional_pkgs <- c("mixedup", "sjPlot", "ape")
pkg_available <- stats::setNames(logical(length(optional_pkgs)), optional_pkgs)

for (pkg in optional_pkgs) {
  pkg_available[pkg] <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-8s: %s\n", pkg,
              ifelse(pkg_available[pkg], "available",
                     "NOT installed (dependent step will be skipped)")))
}
cat("\nTo enable skipped steps, install once, separately, e.g.:\n")
cat('  install.packages(c("sjPlot", "ape"))\n')
cat('  remotes::install_github("m-clark/mixedup")\n\n')


# ==============================================================================
# SECTION 1: DATA PREPARATION
# ==============================================================================

cat("=== DATA PREPARATION ===\n\n")

# FAIRD and ID daily biomass from etraps_daily (zero-capture days already
# filled as real zeros upstream -> Tweedie family is appropriate).
etraps_bio <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, biomass_mg) %>%
  dplyr::rename(Biomass_mg = biomass_mg)

# AMMOD daily biomass — rename to match structure
ammod_bio <- data_list$ammod_biomass_for_stats %>%
  dplyr::mutate(
    Device      = Device,
    Device_type = "AMMOD",
    Biomass_mg  = Live_mass
  ) %>%
  dplyr::select(Date, Device, Device_type, Site, Biomass_mg)

# Combine into a single dataset
bio_combined <- dplyr::bind_rows(etraps_bio, ammod_bio) %>%
  dplyr::arrange(Device_type, Site, Date) %>%
  dplyr::mutate(
    Date.num    = as.numeric(Date),
    Device_type = factor(Device_type, levels = c("AMMOD", "FAIRD", "ID")),
    Site        = factor(Site),
    Device      = factor(Device)
  )

# Restrict to overlapping date range (AMMOD only has Batch 1)
date_range  <- range(ammod_bio$Date)
bio_overlap <- bio_combined %>%
  dplyr::filter(Date >= date_range[1] & Date <= date_range[2])

# --- Design factors (Block structure) ------------------------------------------
# Block: Maize (Sites 1-2) vs Meadow (Sites 3-4), treated as a blocking factor.
bio_overlap$Block <- factor(
  ifelse(bio_overlap$Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")
)
# times: Date as a factor, required for the AR1 structure in glmmTMB.
bio_overlap$times <- factor(bio_overlap$Date)

cat(sprintf("Combined dataset: %d observations\n", nrow(bio_combined)))
cat(sprintf("Overlap period (%s to %s): %d observations\n",
            date_range[1], date_range[2], nrow(bio_overlap)))
cat(sprintf("Device types: %s\n",
            paste(levels(bio_combined$Device_type), collapse = ", ")))
cat(sprintf("Blocks: BLOCK1 = Site1+Site2 (Maize), BLOCK2 = Site3+Site4 (Meadow)\n"))
cat(sprintf("Total biomass variance: %.3g mg^2 (large, due to unit of measurement)\n",
            var(bio_overlap$Biomass_mg, na.rm = TRUE)))

cat("\nObservations per Device_type (overlap period):\n")
print(bio_overlap %>%
        dplyr::group_by(Device_type) %>%
        dplyr::summarise(n = dplyr::n(),
                         n_devices    = dplyr::n_distinct(Device),
                         mean_biomass = round(mean(Biomass_mg, na.rm = TRUE), 1),
                         .groups = "drop"))

cat("\nObservations per Device (overlap period):\n")
print(bio_overlap %>%
        dplyr::group_by(Block, Device_type, Device, Site) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop"))


# ==============================================================================
# SECTION 2: MODEL FITTING (each model defended with tryCatch)
# ==============================================================================

cat("\n=== MODEL FITTING ===\n\n")

# Track convergence status for the plot-model cascade and the final banner.
converged <- c(BX = FALSE, B0 = FALSE, B1 = FALSE, B_simple = FALSE)

# ------------------------------------------------------------------------------
# Model BX: nlme::lme — random slope over time, Block/Site/Device, AR1
# Purpose: nlme allows easy inspection of random effects (BLUPs, augPred,
#          variance components). Convergence is uncertain given the thin design
#          (12 days, 1 device per type per site), hence the tryCatch.
# ------------------------------------------------------------------------------
cat("--- Model BX: nlme::lme (random slope, Block/Site/Device, corAR1) ---\n")
cat("Formula: Biomass_mg ~ bs(Date.num, 3) * Device_type,",
    "random = ~Date.num | Block/Site/Device, correlation = corAR1()\n\n")

BX <- tryCatch(
  nlme::lme(
    Biomass_mg ~ bs(Date.num, 3) * Device_type,
    random      = ~ Date.num | Block / Site / Device,
    correlation = nlme::corAR1(),
    data        = bio_overlap,
    control     = list(opt = "optim")
  ),
  error = function(e) {
    cat("  [FAILED] nlme::lme did not converge:\n    ", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(BX)) {
  converged["BX"] <- TRUE
  cat("  [OK] BX converged.\n\n")
  cat("--- BX summary ---\n"); print(summary(BX))
  cat("\n--- BX Type II ANOVA (car::Anova) ---\n"); print(car::Anova(BX))

  # Random-effects inspection (BLUPs = fixed effects + random deviations)
  cat("\n--- BX random effects (BLUPs via coef) ---\n")
  print(coef(BX))

  # Variance components plot (requires ape::varcomp)
  if (pkg_available["ape"]) {
    vc_plot_file <- "outputs/figures/Scherber_BX_varcomp.png"
    tryCatch({
      grDevices::png(vc_plot_file, width = 900, height = 700, res = 120)
      plot(ape::varcomp(BX))
      grDevices::dev.off()
      cat(sprintf("  [OK] Saved variance-components plot: %s\n", vc_plot_file))
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      cat("  [SKIP] ape::varcomp plot failed:", conditionMessage(e), "\n")
    })
  } else {
    cat("  [SKIP] ape not installed -> variance-components plot skipped.\n")
  }
} 


# ------------------------------------------------------------------------------
# Model B0: glmmTMB Tweedie, Block/Site/Device intercepts, NO AR1
# ------------------------------------------------------------------------------
cat("\n--- Model B0: glmmTMB Tweedie, (1 | Block/Site/Device), no AR1 ---\n")

B0 <- tryCatch(
  glmmTMB::glmmTMB(
    Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Block / Site / Device),
    family = glmmTMB::tweedie(),
    data   = bio_overlap
  ),
  error = function(e) {
    cat("  [FAILED] B0 did not converge:\n    ", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(B0)) { converged["B0"] <- TRUE; cat("  [OK] B0 converged.\n") }


# ------------------------------------------------------------------------------
# Model B1: B0 + AR1 temporal autocorrelation within device
# This is the model Scherber promoted to `B1` and that feeds the figures.
# ------------------------------------------------------------------------------
cat("\n--- Model B1: B0 + ar1(times + 0 | Device) ---\n")

B1 <- tryCatch(
  glmmTMB::glmmTMB(
    Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Block / Site / Device) +
      ar1(times + 0 | Device),
    family = glmmTMB::tweedie(),
    data   = bio_overlap
  ),
  error = function(e) {
    cat("  [FAILED] B1 did not converge:\n    ", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(B1)) { converged["B1"] <- TRUE; cat("  [OK] B1 converged.\n") }


# ------------------------------------------------------------------------------
# Model B_simple: fallback — original simple Tweedie, (1 | Site/Device), no AR1
# Guarantees the plot cascade always has a model to fall back on.
# ------------------------------------------------------------------------------
cat("\n--- Model B_simple (fallback): glmmTMB Tweedie, (1 | Site/Device) ---\n")

B_simple <- tryCatch(
  glmmTMB::glmmTMB(
    Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site / Device),
    family = glmmTMB::tweedie(),
    data   = bio_overlap
  ),
  error = function(e) {
    cat("  [FAILED] B_simple did not converge:\n    ", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(B_simple)) { converged["B_simple"] <- TRUE; cat("  [OK] B_simple converged.\n") }


# --- AIC comparison of the glmmTMB models that converged ----------------------
cat("\n--- AIC comparison (converged glmmTMB models only) ---\n")
tmb_models <- Filter(Negate(is.null), list(B0 = B0, B1 = B1, B_simple = B_simple))
if (length(tmb_models) >= 1) {
  # Build the table by name to avoid passing list names as AIC() arguments.
  aic_tbl <- data.frame(
    Model = names(tmb_models),
    df    = vapply(tmb_models, function(m) attr(logLik(m), "df"), numeric(1)),
    AIC   = vapply(tmb_models, AIC, numeric(1)),
    row.names = NULL
  )
  aic_tbl <- aic_tbl[order(aic_tbl$AIC), ]
  aic_tbl$dAIC <- aic_tbl$AIC - min(aic_tbl$AIC)
  print(aic_tbl, row.names = FALSE)
  if (nrow(aic_tbl) >= 2) {
    cat(sprintf("\n  Best by AIC: %s (dAIC = 0)\n", aic_tbl$Model[1]))
  }
} else {
  cat("  No glmmTMB model converged; AIC comparison skipped.\n")
}


# ==============================================================================
# SECTION 2b: PLOT-MODEL CASCADE (B1 -> B0 -> B_simple) + CONVERGENCE BANNER
# ==============================================================================

cat("\n=== PLOT-MODEL SELECTION ===\n\n")

if (!is.null(B1)) {
  model_for_plots <- B1
  model_for_plots_name <- "B1 (Block/Site/Device + AR1)"
} else if (!is.null(B0)) {
  model_for_plots <- B0
  model_for_plots_name <- "B0 (Block/Site/Device, no AR1)"
  cat("[WARNING] B1 did not converge -> falling back to B0 for figures.\n")
} else if (!is.null(B_simple)) {
  model_for_plots <- B_simple
  model_for_plots_name <- "B_simple (Site/Device, no AR1)"
  cat("[WARNING] Neither B1 nor B0 converged -> falling back to B_simple for figures.\n")
} else {
  model_for_plots <- NULL
  model_for_plots_name <- NA_character_
  cat("[ERROR] No glmmTMB model converged -> figures cannot be generated.\n")
}

cat("\n------------------------------------------------------------------------\n")
cat("CONVERGENCE SUMMARY\n")
cat(sprintf("  BX  (nlme lme, random slope + AR1) : %s\n", ifelse(converged["BX"], "converged", "FAILED")))
cat(sprintf("  B0  (glmmTMB, no AR1)              : %s\n", ifelse(converged["B0"], "converged", "FAILED")))
cat(sprintf("  B1  (glmmTMB + AR1)                : %s\n", ifelse(converged["B1"], "converged", "FAILED")))
cat(sprintf("  B_simple (fallback)                : %s\n", ifelse(converged["B_simple"], "converged", "FAILED")))
cat(sprintf("  --> Model feeding the figures      : %s\n",
            ifelse(is.na(model_for_plots_name), "NONE", model_for_plots_name)))
cat("------------------------------------------------------------------------\n")


# --- Type III ANOVA + interaction p on the plotting model ---------------------
p_interact <- NA_real_
if (!is.null(model_for_plots)) {
  cat("\n--- Type III ANOVA (plotting model) ---\n")
  anova_plot <- tryCatch(car::Anova(model_for_plots, type = "III"),
                         error = function(e) {
                           cat("  [SKIP] Type III ANOVA failed:", conditionMessage(e), "\n")
                           NULL
                         })
  if (!is.null(anova_plot)) {
    print(anova_plot)
    if ("bs(Date.num, 3):Device_type" %in% rownames(anova_plot)) {
      p_interact <- anova_plot["bs(Date.num, 3):Device_type", "Pr(>Chisq)"]
      cat(sprintf("\n  Date x Device_type interaction: p = %.4f %s\n", p_interact,
                  ifelse(p_interact < 0.05,
                         "-> Trajectories NOT parallel",
                         "-> Trajectories parallel")))
    }
  }

  # Variance-component inspection on the plotting model (optional packages)
  if (pkg_available["mixedup"]) {
    cat("\n--- Variance components (mixedup::extract_vc) ---\n")
    tryCatch(print(mixedup::extract_vc(model_for_plots)),
             error = function(e) cat("  [SKIP] mixedup failed:", conditionMessage(e), "\n"))
  } else {
    cat("\n--- Variance components (glmmTMB::VarCorr fallback) ---\n")
    print(glmmTMB::VarCorr(model_for_plots))
  }

  if (pkg_available["sjPlot"]) {
    re_plot_file <- "outputs/figures/Scherber_plotmodel_RE.png"
    tryCatch({
      re_plot <- sjPlot::plot_model(model_for_plots, type = "re", transform = NULL)
      ggplot2::ggsave(re_plot_file, re_plot, width = 9, height = 7, dpi = 150, bg = "white")
      cat(sprintf("  [OK] Saved random-effects plot: %s\n", re_plot_file))
    }, error = function(e) cat("  [SKIP] sjPlot RE plot failed:", conditionMessage(e), "\n"))
  } else {
    cat("  [SKIP] sjPlot not installed -> random-effects plot skipped.\n")
  }
}


# ==============================================================================
# SECTION 3: TEMPORAL TREND PLOTS
# (only if a plotting model is available)
# ==============================================================================

if (is.null(model_for_plots)) {
  cat("\n=== TEMPORAL PLOTS SKIPPED (no converged model) ===\n")
} else {

  cat("\n\n=== GENERATING TEMPORAL PLOTS ===\n\n")
  cat(sprintf("Figures generated from: %s\n\n", model_for_plots_name))

  # Prediction grid (fixed effects only)
  date_seq <- seq(min(bio_overlap$Date.num), max(bio_overlap$Date.num),
                  length.out = 200)
  pred_grid <- expand.grid(
    Date.num    = date_seq,
    Device_type = levels(bio_overlap$Device_type)
  ) %>%
    dplyr::mutate(Site = NA, Device = NA, Block = NA, times = NA)

  pred_grid$fit <- predict(model_for_plots, newdata = pred_grid,
                           type = "response", re.form = NA, se.fit = FALSE)

  pred_link <- predict(model_for_plots, newdata = pred_grid,
                       type = "link", re.form = NA, se.fit = TRUE)

  pred_grid$fit_link <- pred_link$fit
  pred_grid$se_link  <- pred_link$se.fit
  linkinv <- model_for_plots$modelInfo$family$linkinv
  pred_grid$lower <- linkinv(pred_grid$fit_link - 1.96 * pred_grid$se_link)
  pred_grid$upper <- linkinv(pred_grid$fit_link + 1.96 * pred_grid$se_link)

  pred_grid$Date <- as.Date(pred_grid$Date.num, origin = "1970-01-01")
  bio_overlap$Date_plot <- bio_overlap$Date

  # Project palette (custom_colors: FAIRD = blue, ID = magenta, AMMOD = red)
  interaction_label <- ifelse(is.na(p_interact), "not available",
                              sprintf("%.4f", p_interact))

  # --- Plot 1: all device types, log10 y-axis ---
  p_temporal <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_grid,
      ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = Device_type),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = pred_grid,
      ggplot2::aes(x = Date, y = fit, color = Device_type),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = bio_overlap,
      ggplot2::aes(x = Date_plot, y = Biomass_mg, color = Device_type),
      alpha = 0.5, size = 1.8
    ) +
    ggplot2::geom_rug(
      data = bio_overlap,
      ggplot2::aes(x = Date_plot),
      sides = "b", alpha = 0.3, length = ggplot2::unit(0.02, "npc")
    ) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_log10(
      labels = scales::label_comma(),
      breaks = c(10, 50, 100, 250, 500, 1000, 2000, 5000)
    ) +
    ggplot2::labs(
      title    = "Temporal Biomass Trends by Device Type",
      subtitle = sprintf("Model: %s | Date x Device interaction p = %s | log10 y-axis",
                         model_for_plots_name, interaction_label),
      x = "Date", y = "Biomass (mg), log scale",
      color = "Device_type", fill = "Device_type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )

  print(p_temporal)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass.png", p_temporal,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass.png\n")

  # --- Plot 2: FAIRD vs AMMOD only, log10 y-axis ---
  pred_FA     <- pred_grid  %>% dplyr::filter(Device_type %in% c("AMMOD", "FAIRD"))
  bio_FA_plot <- bio_overlap %>% dplyr::filter(Device_type %in% c("AMMOD", "FAIRD"))

  p_temporal_FA <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_FA,
      ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = Device_type),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = pred_FA,
      ggplot2::aes(x = Date, y = fit, color = Device_type),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = bio_FA_plot,
      ggplot2::aes(x = Date_plot, y = Biomass_mg, color = Device_type),
      alpha = 0.5, size = 1.8
    ) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_log10(
      labels = scales::label_comma(),
      breaks = c(10, 50, 100, 250, 500, 1000, 2000, 5000)
    ) +
    ggplot2::labs(
      title    = "Temporal Biomass Trends: FAIRD vs AMMOD",
      subtitle = sprintf("Model: %s | log10 y-axis", model_for_plots_name),
      x = "Date", y = "Biomass (mg), log scale",
      color = "Device_type", fill = "Device_type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )

  print(p_temporal_FA)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_FA_only.png", p_temporal_FA,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_FA_only.png\n")

  # --- Plot 3: linear y-axis version ---
  p_temporal_nolog <- p_temporal +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(
      y = "Biomass (mg)",
      subtitle = sprintf("Model: %s | Date x Device interaction p = %s | linear y-axis",
                         model_for_plots_name, interaction_label)
    )
  print(p_temporal_nolog)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_linear.png", p_temporal_nolog,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_linear.png\n")

  # --- Plot 4: dual-axis (AMMOD scaled onto a secondary axis) ---
  scale_factor <- 10
  bio_plot <- bio_overlap %>%
    dplyr::mutate(
      Biomass_plot = dplyr::if_else(Device_type == "AMMOD",
                                    Biomass_mg / scale_factor, Biomass_mg)
    )
  pred_plot <- pred_grid %>%
    dplyr::mutate(
      fit_plot   = ifelse(Device_type == "AMMOD", fit   / scale_factor, fit),
      lower_plot = ifelse(Device_type == "AMMOD", lower / scale_factor, lower),
      upper_plot = ifelse(Device_type == "AMMOD", upper / scale_factor, upper)
    )
  y_max <- 1200

  p_dual2 <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_plot,
      ggplot2::aes(x = Date, ymin = lower_plot, ymax = upper_plot, fill = Device_type),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = pred_plot,
      ggplot2::aes(x = Date, y = fit_plot, color = Device_type),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = bio_plot,
      ggplot2::aes(x = Date, y = Biomass_plot, color = Device_type),
      alpha = 0.5, size = 1.8
    ) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_continuous(
      name   = "FAIRD / ID Biomass (mg)",
      labels = scales::label_comma(),
      limits = c(0, y_max),
      sec.axis = ggplot2::sec_axis(~ . * scale_factor,
                                   name   = "AMMOD Biomass (mg)",
                                   labels = scales::label_comma())
    ) +
    ggplot2::labs(
      title    = "Temporal Biomass Trends by Device Type",
      subtitle = sprintf("Model: %s | AMMOD on right axis (x%d)",
                         model_for_plots_name, scale_factor),
      x = "Date", color = "Device_type", fill = "Device_type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle      = ggplot2::element_text(size = 10, color = "grey40"),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "top",
      axis.title.y.left  = ggplot2::element_text(color = "grey30"),
      axis.title.y.right = ggplot2::element_text(color = "red"),
      axis.text.y.right  = ggplot2::element_text(color = "red")
    )

  print(p_dual2)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_dualaxis.png", p_dual2,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_dualaxis.png\n")
}


# ==============================================================================
# SECTION 4: EXPORT
# ==============================================================================

cat("\n=== EXPORT ===\n\n")

if (!is.null(model_for_plots) && exists("anova_plot") && !is.null(anova_plot)) {
  anova_df <- as.data.frame(anova_plot) %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::mutate(Model = model_for_plots_name)
  readr::write_csv(anova_df, "outputs/tables/Scherber_glmmTMB_ANOVA_III.csv")
  cat("[OK] Saved: outputs/tables/Scherber_glmmTMB_ANOVA_III.csv\n")
} else {
  cat("[SKIP] No ANOVA table to export (no converged plotting model).\n")
}

cat("\n==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================================\n")

while (sink.number() > 0) sink()
