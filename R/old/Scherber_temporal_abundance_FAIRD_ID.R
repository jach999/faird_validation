# ==============================================================================
# Scherber_temporal_abundance_FAIRD_ID.R
# Temporal abundance analysis of homologous cameras (FAIRD vs ID), mixed models
# ------------------------------------------------------------------------------
# Companion to Scherber_temporal_biomass.R. Same defensive architecture,
# but restricted to the two camera-based systems (FAIRD, ID) and using daily
# abundance (individual counts per day) as the response instead of biomass.
#
# Why cameras only:
#   Pairing the systems keeps each comparison on a solid temporal base. The
#   strict three-system operational overlap collapses to 7 non-consecutive days
#   with a central gap (Aug 28-Sep 1), which is unusable for a spline. The
#   FAIRD-ID pair, by contrast, retains 22 days (conservative) or 17 days
#   (operational). AMMOD-ID is dropped for insufficient overlap and lies outside
#   the H1/H2 framework anyway.
#
# Two runs on the SAME pipeline:
#   - conservative : all 22 days, no operational exclusion.
#   - operational  : 17 days, applying Field_Valid_Operational_ID from
#                    effort_base (removes the ID1/ID2/ID4 battery window,
#                    Aug 28-Sep 1). ID3 is already absent upstream.
#
# Central output:
#   The Date x Device_type interaction is compared side by side between the
#   conservative and operational runs. If FAIRD-ID trajectory non-parallelism
#   survives removal of the ID power artefact, it is a real divergence; if it
#   disappears, it was driven by the battery episode.
#
# Statistical intent (mirrors Scherber edits):
#   - Block structure: Site1+Site2 = BLOCK1 (Maize), Site3+Site4 = BLOCK2 (Meadow)
#   - Spline over time interacted with system: bs(Date.num, 3) * Device_type
#   - Temporal autocorrelation (AR1) on consecutive days within device
#   - Tweedie family (zero-capture days are real zeros)
#   - Model comparison by AIC; defensive fallback cascade for the figures
#
# Defensive design:
#   - Every model is wrapped in tryCatch(); a non-converging model is logged
#     and the pipeline continues.
#   - Plots run on model_for_plots, selected by cascade B1 -> B0 -> B_simple.
#   - A convergence banner reports which model feeds the figures, per run.
#
# NOTE: no package installation is performed here. A missing optional package
#       is logged and its dependent step is skipped.
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
sink("outputs/console/Scherber_temporal_abundance_FAIRD_ID_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("SCHERBER TEMPORAL ABUNDANCE ANALYSIS (FAIRD vs ID, camera pair)\n")
cat("Mixed models: bs(Date, 3) * Device_type + Block/Site/Device nesting + AR1\n")
cat("Two runs compared: conservative (22 days) vs operational (17 days)\n")
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
# SECTION 1: BUILD THE BASE CAMERA-PAIR DATASET (FAIRD + ID)
# ==============================================================================

cat("=== DATA PREPARATION (FAIRD + ID daily abundance) ===\n\n")

# Daily abundance from etraps_daily. Zero-capture days are already filled as
# real zeros upstream (00_load_data.R), so the Tweedie family is appropriate.
etraps_abund <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, abundance)

# Attach the operational flag by Device + Date. effort_base spans the full
# 22-day calendar, so the join is done on the raw abundance frame before any
# date restriction. A left join keeps every abundance row; the flag column is
# what the operational run filters on later.
etraps_abund <- etraps_abund %>%
  dplyr::left_join(
    data_list$effort_base %>%
      dplyr::select(Device, Date, Field_Valid_Operational_ID),
    by = c("Device", "Date")
  )

# Guard: report any abundance rows that failed to match the effort calendar
# (should be zero; a non-zero count signals a Device/Date mismatch to inspect).
n_unmatched <- sum(is.na(etraps_abund$Field_Valid_Operational_ID))
if (n_unmatched > 0) {
  cat(sprintf("[WARNING] %d abundance rows did not match effort_base on Device+Date.\n",
              n_unmatched))
  cat("  These rows will be dropped in the operational run (NA flag).\n\n")
} else {
  cat("[OK] All abundance rows matched effort_base on Device+Date.\n\n")
}

cat(sprintf("Base camera-pair dataset: %d device-day observations\n", nrow(etraps_abund)))
cat(sprintf("Device types: %s\n",
            paste(sort(unique(etraps_abund$Device_type)), collapse = ", ")))
cat(sprintf("Date span: %s to %s\n",
            min(etraps_abund$Date), max(etraps_abund$Date)))


# ==============================================================================
# SECTION 2: PIPELINE FUNCTION (fit + compare, run once per scenario)
# ==============================================================================
#
# run_temporal_analysis() takes an already-prepared data frame for one scenario
# and returns a list with the chosen plotting model, its name, the Date x
# Device_type interaction p-values (BX Type II and best valid glmmTMB Type III),
# the ANOVA tables, and the prediction grid. Both scenarios go through identical
# code so the interaction p-values are directly comparable.
# ==============================================================================

# ------------------------------------------------------------------------------
# tmb_is_valid(): strict convergence check for a glmmTMB fit.
# A model can be returned by glmmTMB() without erroring yet still be numerically
# degenerate (non-positive-definite Hessian, NA AIC). tryCatch does not catch
# this, so it is checked explicitly here. This is what caught the false
# convergence of B1 (AR1) in the operational run, where the 5-day ID gap
# collapses the ar1(times + 0 | Device) structure.
# Returns TRUE only if the fit exists, its sdreport Hessian is positive-definite,
# and its AIC is finite.
# ------------------------------------------------------------------------------
tmb_is_valid <- function(m) {
  if (is.null(m)) return(FALSE)
  pd_hess <- isTRUE(tryCatch(m$sdr$pdHess, error = function(e) FALSE))
  aic_ok  <- tryCatch(is.finite(AIC(m)), error = function(e) FALSE)
  pd_hess && aic_ok
}

run_temporal_analysis <- function(dat, run_label) {

  cat("\n==============================================================================\n")
  cat(sprintf("RUN: %s\n", toupper(run_label)))
  cat("==============================================================================\n\n")

  # --- Design factors -----------------------------------------------------------
  dat <- dat %>%
    dplyr::arrange(Device_type, Site, Date) %>%
    dplyr::mutate(
      Date.num    = as.numeric(Date),
      Device_type = factor(Device_type, levels = c("FAIRD", "ID")),
      Site        = factor(Site),
      Device      = factor(Device)
    )

  # Block: Maize (Sites 1-2) vs Meadow (Sites 3-4), a blocking factor.
  dat$Block <- factor(
    ifelse(dat$Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")
  )
  # times: Date as a factor, required for the AR1 structure in glmmTMB.
  dat$times <- factor(dat$Date)

  n_days <- dplyr::n_distinct(dat$Date)
  cat(sprintf("Observations: %d | distinct days: %d\n", nrow(dat), n_days))
  cat("Blocks: BLOCK1 = Site1+Site2 (Maize), BLOCK2 = Site3+Site4 (Meadow)\n")

  cat("\nObservations per Device_type:\n")
  print(dat %>%
          dplyr::group_by(Device_type) %>%
          dplyr::summarise(n = dplyr::n(),
                           n_devices     = dplyr::n_distinct(Device),
                           mean_abund    = round(mean(abundance, na.rm = TRUE), 2),
                           .groups = "drop"))

  cat("\nObservations per Device:\n")
  print(dat %>%
          dplyr::group_by(Block, Device_type, Device, Site) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop"))

  # --- Model fitting (each defended with tryCatch) ----------------------------
  cat("\n--- MODEL FITTING ---\n\n")
  converged <- c(BX = FALSE, B0 = FALSE, B1 = FALSE, B_simple = FALSE)

  # Model BX: nlme::lme, random slope over time, Block/Site/Device, AR1.
  # nlme is kept for random-effects inspection; convergence is uncertain on a
  # thin design (few days, 1 device per type per site), hence the tryCatch.
  cat("--- Model BX: nlme::lme (random slope, Block/Site/Device, corAR1) ---\n")
  BX <- tryCatch(
    nlme::lme(
      abundance ~ bs(Date.num, 3) * Device_type,
      random      = ~ Date.num | Block / Site / Device,
      correlation = nlme::corAR1(),
      data        = dat,
      control     = list(opt = "optim")
    ),
    error = function(e) {
      cat("  [FAILED] nlme::lme did not converge:\n    ", conditionMessage(e), "\n")
      NULL
    }
  )
  # BX interaction p-value, extracted from the Type II ANOVA. This is the
  # robust reference for the central comparison: BX (gaussian lme) fits cleanly
  # in both runs, unlike the Tweedie AR1 model B1.
  p_interact_BX <- NA_real_
  bx_anova      <- NULL
  if (!is.null(BX)) {
    converged["BX"] <- TRUE
    cat("  [OK] BX converged.\n")
    cat("\n--- BX Type II ANOVA (car::Anova) ---\n")
    bx_anova <- tryCatch(car::Anova(BX),
                         error = function(e) {
                           cat("  [SKIP] BX ANOVA failed:", conditionMessage(e), "\n"); NULL
                         })
    if (!is.null(bx_anova)) {
      print(bx_anova)
      if ("bs(Date.num, 3):Device_type" %in% rownames(bx_anova)) {
        p_interact_BX <- bx_anova["bs(Date.num, 3):Device_type", "Pr(>Chisq)"]
        cat(sprintf("\n  BX Date x Device_type interaction (Type II): p = %.4f %s\n",
                    p_interact_BX,
                    ifelse(p_interact_BX < 0.05,
                           "-> Trajectories NOT parallel",
                           "-> Trajectories parallel")))
      }
    }
    if (pkg_available["ape"]) {
      vc_plot_file <- sprintf("outputs/figures/Scherber_abund_%s_BX_varcomp.png", run_label)
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

  # Model B0: glmmTMB Tweedie, Block/Site/Device intercepts, no AR1.
  cat("\n--- Model B0: glmmTMB Tweedie, (1 | Block/Site/Device), no AR1 ---\n")
  B0 <- tryCatch(
    glmmTMB::glmmTMB(
      abundance ~ bs(Date.num, 3) * Device_type + (1 | Block / Site / Device),
      family = glmmTMB::tweedie(),
      data   = dat
    ),
    error = function(e) {
      cat("  [FAILED] B0 did not converge:\n    ", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(B0)) {
    if (tmb_is_valid(B0)) {
      converged["B0"] <- TRUE; cat("  [OK] B0 converged (valid).\n")
    } else {
      cat("  [INVALID] B0 returned a fit but failed the strict check ",
          "(non-positive-definite Hessian or NA AIC) -> discarded.\n", sep = "")
      B0 <- NULL
    }
  }

  # Model B1: B0 + AR1 temporal autocorrelation within device.
  # Note: ar1(times + 0 | Device) is written without namespace on purpose,
  # ar1() is a glmmTMB formula constructor resolved inside the model call.
  cat("\n--- Model B1: B0 + ar1(times + 0 | Device) ---\n")
  B1 <- tryCatch(
    glmmTMB::glmmTMB(
      abundance ~ bs(Date.num, 3) * Device_type + (1 | Block / Site / Device) +
        ar1(times + 0 | Device),
      family = glmmTMB::tweedie(),
      data   = dat
    ),
    error = function(e) {
      cat("  [FAILED] B1 did not converge:\n    ", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(B1)) {
    if (tmb_is_valid(B1)) {
      converged["B1"] <- TRUE; cat("  [OK] B1 converged (valid).\n")
    } else {
      cat("  [INVALID] B1 returned a fit but failed the strict check.\n")
      cat("            Numerical instability: the ar1(times + 0 | Device) structure\n")
      cat("            is degenerate here (non-positive-definite Hessian, NA AIC).\n")
      cat("            This is the false-convergence case; B1 is discarded and does\n")
      cat("            NOT enter the AIC comparison or the plot cascade.\n")
      B1 <- NULL
    }
  }

  # Model B_simple: fallback, simple Tweedie, (1 | Site/Device), no AR1.
  cat("\n--- Model B_simple (fallback): glmmTMB Tweedie, (1 | Site/Device) ---\n")
  B_simple <- tryCatch(
    glmmTMB::glmmTMB(
      abundance ~ bs(Date.num, 3) * Device_type + (1 | Site / Device),
      family = glmmTMB::tweedie(),
      data   = dat
    ),
    error = function(e) {
      cat("  [FAILED] B_simple did not converge:\n    ", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(B_simple)) {
    if (tmb_is_valid(B_simple)) {
      converged["B_simple"] <- TRUE; cat("  [OK] B_simple converged (valid).\n")
    } else {
      cat("  [INVALID] B_simple returned a fit but failed the strict check ",
          "(non-positive-definite Hessian or NA AIC) -> discarded.\n", sep = "")
      B_simple <- NULL
    }
  }

  # --- AIC comparison of converged glmmTMB models -----------------------------
  cat("\n--- AIC comparison (converged glmmTMB models only) ---\n")
  tmb_models <- Filter(Negate(is.null), list(B0 = B0, B1 = B1, B_simple = B_simple))
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
    if (nrow(aic_tbl) >= 2) {
      cat(sprintf("\n  Best by AIC: %s (dAIC = 0)\n", aic_tbl$Model[1]))
    }
  } else {
    cat("  No glmmTMB model converged; AIC comparison skipped.\n")
  }

  # --- Plot-model cascade (B1 -> B0 -> B_simple) + convergence banner ---------
  # model_for_plots_short is the bare label (B1 / B0 / B_simple) so the central
  # comparison can show WHICH valid glmmTMB fed each Type III cell.
  if (!is.null(B1)) {
    model_for_plots <- B1
    model_for_plots_name  <- "B1 (Block/Site/Device + AR1)"
    model_for_plots_short <- "B1"
  } else if (!is.null(B0)) {
    model_for_plots <- B0
    model_for_plots_name  <- "B0 (Block/Site/Device, no AR1)"
    model_for_plots_short <- "B0"
    cat("\n[WARNING] B1 not valid -> falling back to B0 for figures.\n")
  } else if (!is.null(B_simple)) {
    model_for_plots <- B_simple
    model_for_plots_name  <- "B_simple (Site/Device, no AR1)"
    model_for_plots_short <- "B_simple"
    cat("\n[WARNING] Neither B1 nor B0 valid -> falling back to B_simple for figures.\n")
  } else {
    model_for_plots <- NULL
    model_for_plots_name  <- NA_character_
    model_for_plots_short <- NA_character_
    cat("\n[ERROR] No valid glmmTMB model -> figures cannot be generated for this run.\n")
  }

  cat("\n------------------------------------------------------------------------\n")
  cat(sprintf("CONVERGENCE SUMMARY (%s run)\n", run_label))
  cat(sprintf("  BX  (nlme lme, random slope + AR1) : %s\n", ifelse(converged["BX"], "converged", "FAILED")))
  cat(sprintf("  B0  (glmmTMB, no AR1)              : %s\n", ifelse(converged["B0"], "converged", "FAILED")))
  cat(sprintf("  B1  (glmmTMB + AR1)                : %s\n", ifelse(converged["B1"], "converged", "FAILED")))
  cat(sprintf("  B_simple (fallback)                : %s\n", ifelse(converged["B_simple"], "converged", "FAILED")))
  cat(sprintf("  --> Model feeding the figures      : %s\n",
              ifelse(is.na(model_for_plots_name), "NONE", model_for_plots_name)))
  cat("------------------------------------------------------------------------\n")

  # --- Type III ANOVA + interaction p-value on the plotting model -------------
  p_interact <- NA_real_
  anova_plot <- NULL
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

    # Variance-component inspection on the plotting model.
    if (pkg_available["mixedup"]) {
      cat("\n--- Variance components (mixedup::extract_vc) ---\n")
      tryCatch(print(mixedup::extract_vc(model_for_plots)),
               error = function(e) cat("  [SKIP] mixedup failed:", conditionMessage(e), "\n"))
    } else {
      cat("\n--- Variance components (glmmTMB::VarCorr fallback) ---\n")
      print(glmmTMB::VarCorr(model_for_plots))
    }

    if (pkg_available["sjPlot"]) {
      re_plot_file <- sprintf("outputs/figures/Scherber_abund_%s_RE.png", run_label)
      tryCatch({
        re_plot <- sjPlot::plot_model(model_for_plots, type = "re", transform = NULL)
        ggplot2::ggsave(re_plot_file, re_plot, width = 9, height = 7, dpi = 150, bg = "white")
        cat(sprintf("  [OK] Saved random-effects plot: %s\n", re_plot_file))
      }, error = function(e) cat("  [SKIP] sjPlot RE plot failed:", conditionMessage(e), "\n"))
    } else {
      cat("  [SKIP] sjPlot not installed -> random-effects plot skipped.\n")
    }
  }

  # --- Prediction grid (fixed effects only) -----------------------------------
  pred_grid <- NULL
  if (!is.null(model_for_plots)) {
    date_seq <- seq(min(dat$Date.num), max(dat$Date.num), length.out = 200)
    pred_grid <- expand.grid(
      Date.num    = date_seq,
      Device_type = levels(dat$Device_type)
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
    pred_grid$Date  <- as.Date(pred_grid$Date.num, origin = "1970-01-01")
  }

  # Return everything the caller needs for the cross-run comparison and plots.
  # p_interact_tmb / model_for_plots_short: Type III p-value from the best VALID
  #   glmmTMB, plus the bare label (B0 / B_simple / ...) that produced it.
  # p_interact_BX: Type II p-value from the robust nlme::lme (BX).
  list(
    run_label             = run_label,
    n_days                = n_days,
    data                  = dat,
    model_for_plots       = model_for_plots,
    model_for_plots_name  = model_for_plots_name,
    model_for_plots_short = model_for_plots_short,
    p_interact_tmb        = p_interact,
    p_interact_BX         = p_interact_BX,
    anova_plot            = anova_plot,
    bx_anova              = bx_anova,
    pred_grid             = pred_grid,
    converged             = converged
  )
}


# ==============================================================================
# SECTION 3: BUILD THE TWO SCENARIO DATASETS AND RUN THE PIPELINE
# ==============================================================================

# Conservative: all 22 days, no operational exclusion.
dat_conservative <- etraps_abund %>%
  dplyr::select(Date, Device, Device_type, Site, abundance)

# Operational: 17 days, apply the ID battery-window flag. NA flags (unmatched
# rows, if any) are excluded by the strict TRUE filter.
dat_operational <- etraps_abund %>%
  dplyr::filter(Field_Valid_Operational_ID %in% TRUE) %>%
  dplyr::select(Date, Device, Device_type, Site, abundance)

cat(sprintf("\nConservative scenario: %d rows, %d distinct days\n",
            nrow(dat_conservative), dplyr::n_distinct(dat_conservative$Date)))
cat(sprintf("Operational scenario : %d rows, %d distinct days\n",
            nrow(dat_operational), dplyr::n_distinct(dat_operational$Date)))

res_conservative <- run_temporal_analysis(dat_conservative, "conservative")
res_operational  <- run_temporal_analysis(dat_operational,  "operational")


# ==============================================================================
# SECTION 4: CENTRAL COMPARISON (interaction p-value, conservative vs operational)
# ==============================================================================

cat("\n==============================================================================\n")
cat("CENTRAL COMPARISON: Date x Device_type interaction across runs\n")
cat("Two analytical frames side by side (robustness to family and ANOVA type):\n")
cat("  - BX      : nlme::lme (gaussian), Type II  ANOVA. Robust reference,\n")
cat("              fits cleanly in both runs.\n")
cat("  - glmmTMB : best VALID Tweedie model, Type III ANOVA. Label in parentheses\n")
cat("              shows which model (B0 / B1 / B_simple) fed each cell; B1 is\n")
cat("              dropped when its AR1 structure is numerically degenerate.\n")
cat("==============================================================================\n\n")

fmt_p <- function(p) ifelse(is.na(p), "not available", sprintf("%.4f", p))
verdict <- function(p) {
  if (is.na(p)) return("not available")
  if (p < 0.05) "NOT parallel" else "parallel"
}
# glmmTMB cell: "p (Model)" so a different valid model per run is visible.
fmt_p_labelled <- function(p, short) {
  if (is.na(p)) return(ifelse(is.na(short), "not available",
                              sprintf("not available (%s)", short)))
  sprintf("%.4f (%s)", p, short)
}

comparison_tbl <- data.frame(
  Run                = c(res_conservative$run_label, res_operational$run_label),
  Days               = c(res_conservative$n_days, res_operational$n_days),
  BX_p_TypeII        = c(res_conservative$p_interact_BX, res_operational$p_interact_BX),
  BX_verdict         = c(verdict(res_conservative$p_interact_BX),
                         verdict(res_operational$p_interact_BX)),
  glmmTMB_p_TypeIII  = c(fmt_p_labelled(res_conservative$p_interact_tmb,
                                        res_conservative$model_for_plots_short),
                         fmt_p_labelled(res_operational$p_interact_tmb,
                                        res_operational$model_for_plots_short)),
  glmmTMB_verdict    = c(verdict(res_conservative$p_interact_tmb),
                         verdict(res_operational$p_interact_tmb)),
  stringsAsFactors = FALSE
)
print(comparison_tbl, row.names = FALSE)

cat("\nInterpretation:\n")
cat(sprintf("  Conservative (%d days): BX (Type II) p = %s -> %s | glmmTMB %s (Type III) p = %s -> %s\n",
            res_conservative$n_days,
            fmt_p(res_conservative$p_interact_BX), verdict(res_conservative$p_interact_BX),
            ifelse(is.na(res_conservative$model_for_plots_short), "NONE",
                   res_conservative$model_for_plots_short),
            fmt_p(res_conservative$p_interact_tmb), verdict(res_conservative$p_interact_tmb)))
cat(sprintf("  Operational  (%d days): BX (Type II) p = %s -> %s | glmmTMB %s (Type III) p = %s -> %s\n",
            res_operational$n_days,
            fmt_p(res_operational$p_interact_BX), verdict(res_operational$p_interact_BX),
            ifelse(is.na(res_operational$model_for_plots_short), "NONE",
                   res_operational$model_for_plots_short),
            fmt_p(res_operational$p_interact_tmb), verdict(res_operational$p_interact_tmb)))

cat("\n  Reading the table:\n")
cat("  - BX is the robust anchor. If its interaction is non-significant in BOTH\n")
cat("    runs, FAIRD and ID follow PARALLEL temporal trajectories (they may still\n")
cat("    differ in level via the Device_type main effect).\n")
cat("  - A significant glmmTMB interaction that BX does not corroborate, and that\n")
cat("    coincides with B1 being dropped for AR1 degeneracy, should be read as\n")
cat("    numerical instability of the AR1 term, not as a real divergence.\n")
cat("  - Agreement between BX and a VALID glmmTMB (B0 / B_simple) across both runs\n")
cat("    is the strongest evidence, robust to family and to ANOVA type.\n")

readr::write_csv(comparison_tbl,
                 "outputs/tables/Scherber_abundance_FAIRD_ID_interaction_comparison.csv")
cat("\n[OK] Saved: outputs/tables/Scherber_abundance_FAIRD_ID_interaction_comparison.csv\n")


# ==============================================================================
# SECTION 5: TEMPORAL TREND PLOTS (one figure per run)
# ==============================================================================

# custom_colors provides FAIRD = blue, ID = magenta.
make_temporal_plot <- function(res) {
  if (is.null(res$model_for_plots) || is.null(res$pred_grid)) {
    cat(sprintf("\n[SKIP] No plotting model for the %s run -> figure skipped.\n",
                res$run_label))
    return(invisible(NULL))
  }

  interaction_label_tmb <- fmt_p(res$p_interact_tmb)
  interaction_label_bx  <- fmt_p(res$p_interact_BX)
  dat_pts <- res$data

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = res$pred_grid,
      ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = Device_type),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = res$pred_grid,
      ggplot2::aes(x = Date, y = fit, color = Device_type),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = dat_pts,
      ggplot2::aes(x = Date, y = abundance, color = Device_type),
      alpha = 0.5, size = 1.8
    ) +
    ggplot2::geom_rug(
      data = dat_pts,
      ggplot2::aes(x = Date),
      sides = "b", alpha = 0.3, length = ggplot2::unit(0.02, "npc")
    ) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::labs(
      title    = sprintf("Temporal Abundance Trends: FAIRD vs ID (%s run, %d days)",
                         res$run_label, res$n_days),
      subtitle = sprintf("Curves: %s | interaction p: %s (Type III) / %s (BX, Type II)",
                         res$model_for_plots_name, interaction_label_tmb,
                         interaction_label_bx),
      x = "Date", y = "Abundance (individuals / day)",
      color = "Device_type", fill = "Device_type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )

  out_file <- sprintf("outputs/figures/Scherber_temporal_abundance_FAIRD_ID_%s.png",
                      res$run_label)
  print(p)
  ggplot2::ggsave(out_file, p, width = 10, height = 7, dpi = 300, bg = "white")
  cat(sprintf("[OK] Saved: %s\n", out_file))

  # log10 companion (abundance + 1 to keep zero-capture days on a log axis).
  p_log <- p +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::labs(
      y = "Abundance (individuals / day), log scale",
      subtitle = sprintf("Curves: %s | interaction p: %s (Type III) / %s (BX, Type II) | log10 y-axis",
                         res$model_for_plots_name, interaction_label_tmb,
                         interaction_label_bx)
    )
  out_file_log <- sprintf("outputs/figures/Scherber_temporal_abundance_FAIRD_ID_%s_log10.png",
                          res$run_label)
  tryCatch({
    print(p_log)
    ggplot2::ggsave(out_file_log, p_log, width = 10, height = 7, dpi = 300, bg = "white")
    cat(sprintf("[OK] Saved: %s\n", out_file_log))
  }, error = function(e) cat("  [SKIP] log10 figure failed:", conditionMessage(e), "\n"))

  invisible(p)
}

cat("\n=== GENERATING TEMPORAL PLOTS ===\n\n")
make_temporal_plot(res_conservative)
make_temporal_plot(res_operational)


# ==============================================================================
# SECTION 6: EXPORT ANOVA TABLES
# ==============================================================================

cat("\n=== EXPORT ===\n\n")

export_anova <- function(res) {
  if (!is.null(res$model_for_plots) && !is.null(res$anova_plot)) {
    anova_df <- as.data.frame(res$anova_plot) %>%
      tibble::rownames_to_column("Term") %>%
      dplyr::mutate(Run = res$run_label, Model = res$model_for_plots_name)
    out_file <- sprintf("outputs/tables/Scherber_abundance_FAIRD_ID_%s_ANOVA_III.csv",
                        res$run_label)
    readr::write_csv(anova_df, out_file)
    cat(sprintf("[OK] Saved: %s\n", out_file))
  } else {
    cat(sprintf("[SKIP] No ANOVA table to export for the %s run.\n", res$run_label))
  }
}

export_anova(res_conservative)
export_anova(res_operational)

cat("\n==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================================\n")

while (sink.number() > 0) sink()
