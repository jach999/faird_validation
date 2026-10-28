# ==============================================================================
# Scherber_temporal_biomass.R — Temporal biomass analysis (production version)
# ------------------------------------------------------------------------------
# Compares the daily biomass trajectories of the three monitoring systems
# (AMMOD gravimetric, FAIRD, Insect Detect) over the AMMOD overlap window.
#
# Statistical specification (settled 2026-09-29):
#   - Fixed part: bs(Date.num, 3) * Device_type + Block
#   - Block (BLOCK1 = Maize/Sites1-2, BLOCK2 = Meadow/Sites3-4) enters as a
#     FIXED additive term, not as a random effect. Rationale: Block has only
#     2 levels; fitted as a random effect its variance collapses
#     (Std.Dev. = 0.018, variance = 3.3e-04, i.e. a singular fit), which
#     mis-calibrates every other standard error. The Block-random alternative
#     is still fitted below as a documented robustness check: it gives the same
#     interaction conclusion (p = 0.0052 vs 0.0063), so the manuscript claim
#     does not depend on this choice.
#   - Random part: (1 | Site/Device)
#   - Temporal autocorrelation: ar1(times + 0 | Device)
#   - Model comparison by AIC among VALID fits only
#
# Defensive design:
#   - tmb_is_valid() (reused verbatim from Scherber_temporal_abundance_FAIRD_ID.R)
#     is the convergence criterion, not !is.null(). A glmmTMB fit can be
#     returned without erroring yet be numerically degenerate; tryCatch does
#     not catch that, so pdHess and finite AIC are checked explicitly.
#   - Plot-model cascade B1 -> B0 -> B_simple, each gated by tmb_is_valid().
#   - A convergence banner reports which model actually feeds the figures.
#
# NOTE: package installation is NOT performed here. Missing optional packages
#       are logged and the dependent step is skipped.
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
cat("SCHERBER TEMPORAL BIOMASS ANALYSIS (production version)\n")
cat("Primary: bs(Date,3) * Device_type + Block + (1|Site/Device) + AR1\n")
cat("Block is a FIXED term (2 levels); Block-random kept as robustness check.\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 0: PACKAGE AVAILABILITY CHECK (no installation performed)
# ==============================================================================

cat("=== OPTIONAL PACKAGE AVAILABILITY ===\n\n")

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
# tmb_is_valid(): strict convergence check for a glmmTMB fit.
# Reused verbatim from Scherber_temporal_abundance_FAIRD_ID.R so the validity
# criterion is identical across the project. Returns TRUE only if the fit
# exists, its sdreport Hessian is positive-definite, and its AIC is finite.
# ==============================================================================

tmb_is_valid <- function(m) {
  if (is.null(m)) return(FALSE)
  pd_hess <- isTRUE(tryCatch(m$sdr$pdHess, error = function(e) FALSE))
  aic_ok  <- tryCatch(is.finite(AIC(m)), error = function(e) FALSE)
  pd_hess && aic_ok
}

# Convenience: fit a glmmTMB model, trap errors, report validity in one line.
fit_tmb <- function(formula, data, label) {
  cat(sprintf("\n--- %s ---\n", label))
  m <- tryCatch(
    glmmTMB::glmmTMB(formula, family = glmmTMB::tweedie(), data = data),
    error = function(e) {
      cat("  [FAILED] model errored:\n    ", conditionMessage(e), "\n")
      NULL
    }
  )
  if (tmb_is_valid(m)) {
    cat("  [VALID] positive-definite Hessian, finite AIC.\n")
  } else if (!is.null(m)) {
    cat("  [INVALID] fit returned but numerically degenerate (pdHess/AIC check failed).\n")
  }
  m
}


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

# AMMOD daily biomass, renamed to match structure.
ammod_bio <- data_list$ammod_biomass_for_stats %>%
  dplyr::mutate(
    Device      = Device,
    Device_type = "AMMOD",
    Biomass_mg  = Live_mass
  ) %>%
  dplyr::select(Date, Device, Device_type, Site, Biomass_mg)

# Combine into a single dataset.
bio_combined <- dplyr::bind_rows(etraps_bio, ammod_bio) %>%
  dplyr::arrange(Device_type, Site, Date) %>%
  dplyr::mutate(
    Date.num    = as.numeric(Date),
    Device_type = factor(Device_type, levels = c("AMMOD", "FAIRD", "ID")),
    Site        = factor(Site),
    Device      = factor(Device)
  )

# Restrict to the overlapping date range (AMMOD only has Batch 1: 12 days).
date_range  <- range(ammod_bio$Date)
bio_overlap <- bio_combined %>%
  dplyr::filter(Date >= date_range[1] & Date <= date_range[2])

# Design factors.
# Block: Maize (Sites 1-2) vs Meadow (Sites 3-4), a blocking factor.
bio_overlap$Block <- factor(
  ifelse(bio_overlap$Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")
)
# times: Date as a factor, required for the AR1 structure in glmmTMB.
bio_overlap$times <- factor(bio_overlap$Date)

cat(sprintf("Combined dataset: %d observations\n", nrow(bio_combined)))
cat(sprintf("Overlap period (%s to %s): %d observations\n",
            date_range[1], date_range[2], nrow(bio_overlap)))
cat(sprintf("Distinct days: %d\n", dplyr::n_distinct(bio_overlap$Date)))
cat("Blocks: BLOCK1 = Site1+Site2 (Maize), BLOCK2 = Site3+Site4 (Meadow)\n")
cat(sprintf("Block levels: %d | Site levels: %d | Device levels: %d\n",
            dplyr::n_distinct(bio_overlap$Block),
            dplyr::n_distinct(bio_overlap$Site),
            dplyr::n_distinct(bio_overlap$Device)))

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
# SECTION 2: MODEL FITTING
# ==============================================================================

cat("\n=== MODEL FITTING ===\n")

# ------------------------------------------------------------------------------
# Model B1 (PRIMARY): Block fixed, (1|Site/Device), AR1 within device.
# ------------------------------------------------------------------------------
B1 <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + Block +
    (1 | Site / Device) + ar1(times + 0 | Device),
  bio_overlap,
  "Model B1 (PRIMARY): Block fixed + (1|Site/Device) + AR1"
)

# ------------------------------------------------------------------------------
# Model B0: same fixed part, no AR1. Isolates the contribution of the temporal
# autocorrelation structure.
# ------------------------------------------------------------------------------
B0 <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + Block + (1 | Site / Device),
  bio_overlap,
  "Model B0: Block fixed + (1|Site/Device), no AR1"
)

# ------------------------------------------------------------------------------
# Model B_ranef (ROBUSTNESS CHECK, not the reported model):
# Block as a random effect, the specification sketched during the Scherber
# meeting. Retained so the manuscript can state that the interaction conclusion
# is invariant to this choice, and so the Block variance collapse is on record.
# ------------------------------------------------------------------------------
B_ranef <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type +
    (1 | Block / Site / Device) + ar1(times + 0 | Device),
  bio_overlap,
  "Model B_ranef (robustness check): Block RANDOM + AR1"
)

# ------------------------------------------------------------------------------
# Model B_simple: fallback, no Block at all. Guarantees the plot cascade
# always has something to fall back on.
# ------------------------------------------------------------------------------
B_simple <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site / Device),
  bio_overlap,
  "Model B_simple (fallback): (1|Site/Device), no Block, no AR1"
)

# ------------------------------------------------------------------------------
# Model BX: nlme::lme cross-check, Gaussian, random slope over time.
# Purpose: nlme exposes BLUPs and variance components directly. Block is fixed
# here too, for consistency with the primary specification. Convergence is
# uncertain given the thin design (12 days, 1 device per type per site).
# ------------------------------------------------------------------------------
cat("\n--- Model BX (cross-check): nlme::lme, Gaussian, random slope + corAR1 ---\n")
cat("Formula: Biomass_mg ~ bs(Date.num,3) * Device_type + Block,",
    "random = ~Date.num | Site/Device, correlation = corAR1()\n")

BX <- tryCatch(
  nlme::lme(
    Biomass_mg ~ bs(Date.num, 3) * Device_type + Block,
    random      = ~ Date.num | Site / Device,
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
  cat("  [OK] BX converged.\n\n")
  cat("--- BX summary ---\n"); print(summary(BX))
  cat("\n--- BX Type II ANOVA (car::Anova) ---\n"); print(car::Anova(BX))
  cat("\n--- BX random effects (BLUPs via coef) ---\n"); print(coef(BX))

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


# ==============================================================================
# SECTION 2b: AIC COMPARISON (valid glmmTMB fits only)
# ==============================================================================

cat("\n\n--- AIC comparison (VALID glmmTMB fits only) ---\n")

tmb_models <- list(B1 = B1, B0 = B0, B_ranef = B_ranef, B_simple = B_simple)
tmb_models <- tmb_models[vapply(tmb_models, tmb_is_valid, logical(1))]

if (length(tmb_models) >= 1) {
  aic_tbl <- data.frame(
    Model = names(tmb_models),
    df    = vapply(tmb_models, function(m) attr(logLik(m), "df"), numeric(1)),
    AIC   = vapply(tmb_models, AIC, numeric(1)),
    row.names = NULL
  )
  aic_tbl <- aic_tbl[order(aic_tbl$AIC), ]
  aic_tbl$dAIC <- aic_tbl$AIC - min(aic_tbl$AIC)
  print(aic_tbl, row.names = FALSE)
  cat("\n  NOTE: B_ranef is NOT AIC-comparable with the Block-fixed models as a\n")
  cat("  selection criterion (a factor treated as fixed vs random is not nested\n")
  cat("  in the usual sense). It is listed for reference only.\n")
} else {
  cat("  No valid glmmTMB fit; AIC comparison skipped.\n")
}


# ==============================================================================
# SECTION 2c: PLOT-MODEL CASCADE (B1 -> B0 -> B_simple) + CONVERGENCE BANNER
# ==============================================================================

cat("\n=== PLOT-MODEL SELECTION ===\n\n")

if (tmb_is_valid(B1)) {
  model_for_plots      <- B1
  model_for_plots_name <- "B1 (Block fixed, Site/Device, AR1)"
} else if (tmb_is_valid(B0)) {
  model_for_plots      <- B0
  model_for_plots_name <- "B0 (Block fixed, Site/Device, no AR1)"
  cat("[WARNING] B1 not valid -> falling back to B0 for figures.\n")
} else if (tmb_is_valid(B_simple)) {
  model_for_plots      <- B_simple
  model_for_plots_name <- "B_simple (Site/Device, no Block, no AR1)"
  cat("[WARNING] Neither B1 nor B0 valid -> falling back to B_simple for figures.\n")
} else {
  model_for_plots      <- NULL
  model_for_plots_name <- NA_character_
  cat("[ERROR] No valid glmmTMB fit -> figures cannot be generated.\n")
}

cat("\n------------------------------------------------------------------------\n")
cat("CONVERGENCE SUMMARY (tmb_is_valid: pdHess + finite AIC)\n")
cat(sprintf("  B1  (PRIMARY: Block fixed + AR1)   : %s\n",
            ifelse(tmb_is_valid(B1), "VALID", "INVALID")))
cat(sprintf("  B0  (Block fixed, no AR1)          : %s\n",
            ifelse(tmb_is_valid(B0), "VALID", "INVALID")))
cat(sprintf("  B_ranef (Block random, check)      : %s\n",
            ifelse(tmb_is_valid(B_ranef), "VALID", "INVALID")))
cat(sprintf("  B_simple (fallback)                : %s\n",
            ifelse(tmb_is_valid(B_simple), "VALID", "INVALID")))
cat(sprintf("  BX  (nlme cross-check)             : %s\n",
            ifelse(is.null(BX), "FAILED", "converged")))
cat(sprintf("  --> Model feeding the figures      : %s\n",
            ifelse(is.na(model_for_plots_name), "NONE", model_for_plots_name)))
cat("------------------------------------------------------------------------\n")


# ==============================================================================
# SECTION 2d: TYPE III ANOVA + ROBUSTNESS COMPARISON
# ==============================================================================

interaction_term <- "bs(Date.num, 3):Device_type"
p_interact <- NA_real_

get_p <- function(m, term) {
  if (!tmb_is_valid(m)) return(NA_real_)
  av <- tryCatch(car::Anova(m, type = "III"), error = function(e) NULL)
  if (is.null(av) || !(term %in% rownames(av))) return(NA_real_)
  av[term, "Pr(>Chisq)"]
}

if (!is.null(model_for_plots)) {
  cat("\n--- Type III ANOVA (plotting model) ---\n")
  anova_plot <- tryCatch(car::Anova(model_for_plots, type = "III"),
                         error = function(e) {
                           cat("  [SKIP] Type III ANOVA failed:", conditionMessage(e), "\n")
                           NULL
                         })
  if (!is.null(anova_plot)) {
    print(anova_plot)
    if (interaction_term %in% rownames(anova_plot)) {
      p_interact <- anova_plot[interaction_term, "Pr(>Chisq)"]
      cat(sprintf("\n  Date x Device_type interaction: p = %.4f %s\n", p_interact,
                  ifelse(p_interact < 0.05,
                         "-> Trajectories NOT parallel",
                         "-> Trajectories parallel")))
    }
  }

  cat("\n--- Fixed-effect coefficients (plotting model) ---\n")
  print(glmmTMB::fixef(model_for_plots)$cond)

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

# --- Robustness: does the interaction conclusion depend on the Block choice? ---
cat("\n--- ROBUSTNESS: Block fixed vs Block random ---\n")

p_fixed  <- get_p(B1,      interaction_term)
p_random <- get_p(B_ranef, interaction_term)

block_var <- NA_real_
if (tmb_is_valid(B_ranef)) {
  vc <- tryCatch(glmmTMB::VarCorr(B_ranef)$cond, error = function(e) NULL)
  if (!is.null(vc) && "Block" %in% names(vc)) {
    block_var <- as.numeric(attr(vc$Block, "stddev"))^2
  }
}

robustness <- data.frame(
  Specification = c("Block FIXED (reported)", "Block RANDOM (check)"),
  p_interaction = c(p_fixed, p_random),
  Block_variance = c(NA_real_, block_var),
  stringsAsFactors = FALSE
)
print(robustness, row.names = FALSE)

cat("\n  Interpretation: a Block variance at or near zero confirms that the\n")
cat("  2-level random effect is a singular fit, and that the fixed\n")
cat("  specification is the better-calibrated one to report.\n")

readr::write_csv(robustness,
                 "outputs/tables/Scherber_biomass_block_robustness.csv")
cat("[OK] Saved: outputs/tables/Scherber_biomass_block_robustness.csv\n")


# ==============================================================================
# SECTION 3: TEMPORAL TREND PLOTS
# ==============================================================================

if (is.null(model_for_plots)) {
  cat("\n=== TEMPORAL PLOTS SKIPPED (no valid model) ===\n")
} else {

  cat("\n\n=== GENERATING TEMPORAL PLOTS ===\n\n")
  cat(sprintf("Figures generated from: %s\n\n", model_for_plots_name))

  # Prediction grid (fixed effects only). Block is a fixed term in the primary
  # model, so it must carry a real level: predictions are made at the reference
  # block (BLOCK1 / Maize) and this is stated in the subtitle.
  ref_block <- levels(bio_overlap$Block)[1]

  date_seq <- seq(min(bio_overlap$Date.num), max(bio_overlap$Date.num),
                  length.out = 200)
  pred_grid <- expand.grid(
    Date.num    = date_seq,
    Device_type = levels(bio_overlap$Device_type)
  ) %>%
    dplyr::mutate(
      Site   = NA,
      Device = NA,
      times  = NA,
      Block  = factor(ref_block, levels = levels(bio_overlap$Block))
    )

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

  interaction_label <- ifelse(is.na(p_interact), "not available",
                              sprintf("%.4f", p_interact))

  subtitle_base <- sprintf(
    "Model: %s | predictions at %s | Date x Device interaction p = %s",
    model_for_plots_name, ref_block, interaction_label
  )

  # --- Plot 1: all device types, log10 y-axis --------------------------------
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
      subtitle = paste0(subtitle_base, " | log10 y-axis"),
      x = "Date", y = "Biomass (mg), log scale",
      color = "Device_type", fill = "Device_type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 9, color = "grey40"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )

  print(p_temporal)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass.png", p_temporal,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass.png\n")

  # --- Plot 2: FAIRD vs AMMOD only, log10 y-axis -----------------------------
  # This is the pair that carries the validation argument; ID is dropped so the
  # two systems being compared are not visually crowded.
  pred_FA     <- pred_grid   %>% dplyr::filter(Device_type %in% c("AMMOD", "FAIRD"))
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
      subtitle = paste0(subtitle_base, " | log10 y-axis"),
      x = "Date", y = "Biomass (mg), log scale",
      color = "Device_type", fill = "Device_type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 9, color = "grey40"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )

  print(p_temporal_FA)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_FA_only.png", p_temporal_FA,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_FA_only.png\n")

  # --- Plot 3: linear y-axis version -----------------------------------------
  p_temporal_nolog <- p_temporal +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(
      y = "Biomass (mg)",
      subtitle = paste0(subtitle_base, " | linear y-axis")
    )
  print(p_temporal_nolog)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_linear.png", p_temporal_nolog,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_linear.png\n")

  # --- Plot 4: dual-axis (AMMOD scaled onto a secondary axis) ----------------
  # AMMOD biomass is roughly an order of magnitude above the camera systems;
  # the secondary axis makes the SHAPE of the trajectories comparable without
  # a log transform. The scale factor is cosmetic, not a model quantity.
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
      subtitle = sprintf("%s | AMMOD on right axis (x%d)",
                         subtitle_base, scale_factor),
      x = "Date", color = "Device_type", fill = "Device_type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle      = ggplot2::element_text(size = 9, color = "grey40"),
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
  readr::write_csv(anova_df, "outputs/tables/Scherber_temporal_biomass_ANOVA_III.csv")
  cat("[OK] Saved: outputs/tables/Scherber_temporal_biomass_ANOVA_III.csv\n")
} else {
  cat("[SKIP] No ANOVA table to export (no valid plotting model).\n")
}

cat("\n==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================================\n")

while (sink.number() > 0) sink()
