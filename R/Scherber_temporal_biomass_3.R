# ==============================================================================
# Scherber_temporal_biomass.R, temporal biomass analysis (production version)
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
#   - Pair model AMMOD vs FAIRD (Section 2e, restores former Model B2 of
#     glmmTMB.R): same specification as B1 on the two-system subset, so the
#     Date x Device_type interaction is tested with 3 df instead of the 6 df of
#     the three-system omnibus test. Fitted on two data variants (FULL and
#     PAIRED, see Section 2e), plus a post hoc sensitivity refit of PAIRED
#     without the first days of the series. The interaction is tested by Wald
#     and by LRT, the LRT falling back down the cascade when the additive
#     model is not valid at the reported level.
#
# Presentation (Scherber request): no logarithmic axes. Figures use a linear
# y-axis; the dual-axis figures show AMMOD on a secondary axis divided by a
# scale factor DERIVED FROM THE DATA (ratio of AMMOD to FAIRD biomass over
# co-located site-days), not a hardcoded constant. The y view is capped at the
# rounded 98th percentile with coord_cartesian(), which trims the view without
# dropping data from the fit or the ribbons. The model itself keeps the log
# link of the Tweedie family; only the display scale is linear.
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
cat("Pair model AMMOD vs FAIRD (3-df interaction) in Section 2e.\n")
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

# nice_ceiling(): round a positive number UP to a readable axis limit
# (1, 1.5, 2, 2.5, ... times a power of ten). Used for the 98th-percentile cap.
nice_ceiling <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1 || !is.finite(x) || x <= 0) return(NA_real_)
  step <- 10^floor(log10(x)) / 2
  ceiling(x / step) * step
}

# ammod_faird_ratio(): data-derived scale factor for the dual-axis figures.
# Uses only site-days where BOTH systems have a value at the same site, so the
# ratio is not distorted by sites or days present for one system only
# (Site4 FAIRD, Site3 FAIRD after the AMMOD3 ATL loss). Returns the ratio of
# summed biomass (equivalent to a ratio of means over the same site-days) and
# the daily pooled ratio series, which shows whether the ratio drifts over time.
ammod_faird_ratio <- function(dat) {
  paired <- dat %>%
    dplyr::filter(Device_type %in% c("AMMOD", "FAIRD"), !is.na(Biomass_mg)) %>%
    dplyr::select(Site, Date, Device_type, Biomass_mg) %>%
    dplyr::mutate(Device_type = as.character(Device_type)) %>%
    tidyr::pivot_wider(names_from = Device_type, values_from = Biomass_mg) %>%
    dplyr::filter(!is.na(AMMOD), !is.na(FAIRD))
  daily <- paired %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(
      n_sites   = dplyr::n(),
      AMMOD_sum = sum(AMMOD),
      FAIRD_sum = sum(FAIRD),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(ratio = dplyr::if_else(FAIRD_sum > 0,
                                         AMMOD_sum / FAIRD_sum, NA_real_))
  list(
    overall = sum(paired$AMMOD) / sum(paired$FAIRD),
    daily   = daily,
    n_site_days = nrow(paired)
  )
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
# SECTION 2e: PAIR MODEL, AMMOD vs FAIRD (restores former Model B2)
# ------------------------------------------------------------------------------
# The omnibus test above (6 df) says the three trajectories are not all
# parallel, but not WHICH pair departs. The validation argument rests on AMMOD
# vs FAIRD, so this pair is tested on its own with the same specification as
# B1 (Block fixed, (1|Site/Device), AR1, Tweedie). The interaction has 3 df.
#
# Two data variants, because system composition differs by site and day:
#   FULL   : all AMMOD and FAIRD observations in the overlap window. FAIRD
#            contributes Site4 and the Site3 days after the AMMOD3 ATL loss,
#            AMMOD does not; part of a trajectory difference can therefore come
#            from site composition rather than from the systems.
#   PAIRED : only site-days where BOTH systems have a value (co-located, same
#            day). Composition is identical for the two systems, so the
#            interaction cannot be produced by site imbalance.
#
# Each variant runs its own cascade AR1 -> no AR1 -> no Block, gated by
# tmb_is_valid(). The interaction is tested by Type III Wald chi-square and,
# as a small-sample check, by a likelihood-ratio test against the additive
# model with the same random and AR1 structure.
# ==============================================================================

cat("\n\n=== PAIR MODEL: AMMOD vs FAIRD ===\n")

pair_specs <- list(
  AR1 = list(
    full = "Biomass_mg ~ bs(Date.num, 3) * Device_type + Block + (1 | Site / Device) + ar1(times + 0 | Device)",
    add  = "Biomass_mg ~ bs(Date.num, 3) + Device_type + Block + (1 | Site / Device) + ar1(times + 0 | Device)",
    name = "Block fixed + (1|Site/Device) + AR1"
  ),
  noAR1 = list(
    full = "Biomass_mg ~ bs(Date.num, 3) * Device_type + Block + (1 | Site / Device)",
    add  = "Biomass_mg ~ bs(Date.num, 3) + Device_type + Block + (1 | Site / Device)",
    name = "Block fixed + (1|Site/Device), no AR1"
  ),
  simple = list(
    full = "Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site / Device)",
    add  = "Biomass_mg ~ bs(Date.num, 3) + Device_type + (1 | Site / Device)",
    name = "(1|Site/Device), no Block, no AR1"
  )
)

fit_pair_cascade <- function(dat, variant) {
  cat(sprintf("\n##### Variant %s #####\n", variant))
  for (lvl in names(pair_specs)) {
    spec <- pair_specs[[lvl]]
    m <- fit_tmb(stats::as.formula(spec$full), dat,
                 sprintf("Pair %s | %s", variant, spec$name))
    if (tmb_is_valid(m)) {
      if (lvl != "AR1") {
        cat(sprintf("  [WARNING] higher cascade levels invalid -> using '%s'.\n", lvl))
      }
      return(list(model = m, level = lvl, spec = spec))
    }
  }
  cat("  [ERROR] no valid fit anywhere in the cascade.\n")
  list(model = NULL, level = NA_character_, spec = NULL)
}

# pair_lrt(): likelihood-ratio test of the interaction (full vs additive fixed
# part, identical random and AR1 structure). Non-positive-definite Hessian
# warnings printed by glmmTMB during this step come from fits that
# tmb_is_valid() then rejects; they are expected and handled. Tried first at the cascade level
# of the reported model; if the full or the additive fit is not valid there,
# it moves down the cascade and refits BOTH models at the lower level, so the
# two models compared always share their structure. The Wald p-value at the
# level where the LRT was computed is returned too, so Wald and LRT are
# compared like for like.
pair_lrt <- function(res, dat, variant) {
  empty <- list(p_lrt = NA_real_, p_wald = NA_real_, level = NA_character_)
  if (is.null(res$model)) return(empty)
  lvls  <- names(pair_specs)
  start <- match(res$level, lvls)
  for (lvl in lvls[start:length(lvls)]) {
    spec   <- pair_specs[[lvl]]
    m_full <- if (lvl == res$level) res$model else tryCatch(
      glmmTMB::glmmTMB(stats::as.formula(spec$full),
                       family = glmmTMB::tweedie(), data = dat),
      error = function(e) NULL)
    m_add <- tryCatch(
      glmmTMB::glmmTMB(stats::as.formula(spec$add),
                       family = glmmTMB::tweedie(), data = dat),
      error = function(e) NULL)
    if (tmb_is_valid(m_full) && tmb_is_valid(m_add)) {
      # LRT computed directly from the log-likelihoods. stats::anova() is not
      # used because anova.glmmTMB compares the deparsed name of the data
      # argument in each model call, and the two fits are created by different
      # functions (data vs dat), which triggers a false "not the same data"
      # error. Both models are ML fits on the same data frame by construction.
      ll_full <- stats::logLik(m_full)
      ll_add  <- stats::logLik(m_add)
      lr_stat <- 2 * (as.numeric(ll_full) - as.numeric(ll_add))
      lr_df   <- attr(ll_full, "df") - attr(ll_add, "df")
      p_lrt   <- stats::pchisq(lr_stat, df = lr_df, lower.tail = FALSE)
      lrt_tbl <- data.frame(
        Model  = c("additive", "interaction"),
        df     = c(attr(ll_add, "df"), attr(ll_full, "df")),
        logLik = round(c(as.numeric(ll_add), as.numeric(ll_full)), 3),
        AIC    = round(c(AIC(m_add), AIC(m_full)), 2),
        Chisq  = c(NA, round(lr_stat, 3)),
        Chi_df = c(NA, lr_df),
        p      = c(NA, p_lrt)
      )
      cat(sprintf("\n--- LRT: interaction vs additive, pair %s, level '%s' ---\n",
                  variant, lvl))
      if (lvl != res$level) {
        cat("  [NOTE] LRT computed at a LOWER cascade level than the reported Wald\n")
        cat("         test, because the additive model was not valid at that level.\n")
      }
      print(lrt_tbl, row.names = FALSE)
      p_wald <- tryCatch({
        av <- car::Anova(m_full, type = "III")
        av[interaction_term, "Pr(>Chisq)"]
      }, error = function(e) NA_real_)
      cat(sprintf("  Same level: Wald p = %.4f | LRT p = %.4f\n", p_wald, p_lrt))
      return(list(p_lrt = p_lrt, p_wald = p_wald, level = lvl))
    }
    cat(sprintf("  [INFO] LRT not possible at level '%s' (full: %s, additive: %s).\n",
                lvl,
                ifelse(tmb_is_valid(m_full), "valid", "INVALID"),
                ifelse(tmb_is_valid(m_add),  "valid", "INVALID")))
  }
  cat("  [SKIP] no cascade level with both models valid -> LRT not computed.\n")
  empty
}

summarise_pair <- function(res, dat, variant) {
  out <- data.frame(
    Variant            = variant,
    n_obs              = nrow(dat),
    n_devices          = dplyr::n_distinct(dat$Device),
    Model              = ifelse(is.null(res$spec), NA_character_, res$spec$name),
    Chisq_interaction  = NA_real_,
    df_interaction     = NA_real_,
    p_interaction_Wald = NA_real_,
    p_Device_type      = NA_real_,
    p_spline           = NA_real_,
    p_Block            = NA_real_,
    p_interaction_LRT  = NA_real_,
    LRT_level          = NA_character_,
    p_Wald_at_LRT_level = NA_real_,
    AIC                = NA_real_,
    stringsAsFactors   = FALSE
  )
  if (is.null(res$model)) return(out)
  m <- res$model
  out$AIC <- AIC(m)

  av <- tryCatch(car::Anova(m, type = "III"), error = function(e) {
    cat("  [SKIP] Type III ANOVA failed:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(av)) {
    cat(sprintf("\n--- Type III ANOVA, pair %s ---\n", variant)); print(av)
    if (interaction_term %in% rownames(av)) {
      out$Chisq_interaction  <- av[interaction_term, "Chisq"]
      out$df_interaction     <- av[interaction_term, "Df"]
      out$p_interaction_Wald <- av[interaction_term, "Pr(>Chisq)"]
    }
    if ("Device_type" %in% rownames(av))     out$p_Device_type <- av["Device_type", "Pr(>Chisq)"]
    if ("bs(Date.num, 3)" %in% rownames(av)) out$p_spline      <- av["bs(Date.num, 3)", "Pr(>Chisq)"]
    if ("Block" %in% rownames(av))           out$p_Block       <- av["Block", "Pr(>Chisq)"]
  }

  # Likelihood-ratio test of the interaction, with cascade fallback.
  lrt_res <- pair_lrt(res, dat, variant)
  out$p_interaction_LRT       <- lrt_res$p_lrt
  out$LRT_level               <- lrt_res$level
  out$p_Wald_at_LRT_level     <- lrt_res$p_wald

  cat(sprintf("\n--- Fixed-effect coefficients, pair %s ---\n", variant))
  print(glmmTMB::fixef(m)$cond)
  cat(sprintf("\n--- Variance components, pair %s ---\n", variant))
  print(glmmTMB::VarCorr(m))
  out
}

# --- Data variants ------------------------------------------------------------
bio_FA_full <- bio_overlap %>%
  dplyr::filter(Device_type %in% c("AMMOD", "FAIRD"), !is.na(Biomass_mg)) %>%
  droplevels()

bio_FA_paired <- bio_FA_full %>%
  dplyr::group_by(Site, Date) %>%
  dplyr::filter(dplyr::n_distinct(Device_type) == 2) %>%
  dplyr::ungroup() %>%
  droplevels()

cat("\nObservations per Device, FULL variant:\n")
print(bio_FA_full %>%
        dplyr::count(Block, Device_type, Device, Site, name = "n"))
cat("\nObservations per Device, PAIRED variant:\n")
print(bio_FA_paired %>%
        dplyr::count(Block, Device_type, Device, Site, name = "n"))

# --- Fit both variants --------------------------------------------------------
pair_full   <- fit_pair_cascade(bio_FA_full,   "FULL")
pair_paired <- fit_pair_cascade(bio_FA_paired, "PAIRED")

# --- Sensitivity (POST HOC): PAIRED without the first days of the series ------
# The per-site diagnostic (diagnostic_per_site_biomass.R) showed the
# AMMOD/FAIRD ratio is high on the first days at every site and then drops.
# This refit asks whether the interaction persists once those days are
# removed. It is a post hoc sensitivity check, chosen after seeing the data,
# and must be reported as such; it does not replace the PAIRED model.
n_drop_initial  <- 3
initial_dates   <- sort(unique(bio_FA_paired$Date))[seq_len(n_drop_initial)]
bio_FA_paired_trim <- bio_FA_paired %>%
  dplyr::filter(!(Date %in% initial_dates)) %>%
  droplevels()
sens_label <- sprintf("PAIRED_SENS_drop_first%d", n_drop_initial)

cat(sprintf("\nSensitivity variant %s: dates removed = %s\n", sens_label,
            paste(format(initial_dates), collapse = ", ")))
cat(sprintf("  Remaining: %d observations, %d days (%s to %s)\n",
            nrow(bio_FA_paired_trim), dplyr::n_distinct(bio_FA_paired_trim$Date),
            min(bio_FA_paired_trim$Date), max(bio_FA_paired_trim$Date)))

pair_sens <- fit_pair_cascade(bio_FA_paired_trim, sens_label)

pair_tbl <- dplyr::bind_rows(
  summarise_pair(pair_full,   bio_FA_full,        "FULL"),
  summarise_pair(pair_paired, bio_FA_paired,      "PAIRED"),
  summarise_pair(pair_sens,   bio_FA_paired_trim, sens_label)
)

# Omnibus row for side-by-side reading (6-df test from the plotting model).
omnibus_row <- data.frame(
  Variant            = "OMNIBUS (3 systems)",
  n_obs              = nrow(bio_overlap),
  n_devices          = dplyr::n_distinct(bio_overlap$Device),
  Model              = model_for_plots_name,
  Chisq_interaction  = if (exists("anova_plot") && !is.null(anova_plot) &&
                           interaction_term %in% rownames(anova_plot))
                         anova_plot[interaction_term, "Chisq"] else NA_real_,
  df_interaction     = if (exists("anova_plot") && !is.null(anova_plot) &&
                           interaction_term %in% rownames(anova_plot))
                         anova_plot[interaction_term, "Df"] else NA_real_,
  p_interaction_Wald = p_interact,
  stringsAsFactors   = FALSE
)
pair_tbl <- dplyr::bind_rows(omnibus_row, pair_tbl)

cat("\n------------------------------------------------------------------------\n")
cat("DATE x DEVICE_TYPE INTERACTION: OMNIBUS vs AMMOD-FAIRD PAIR\n")
cat("------------------------------------------------------------------------\n")
print(pair_tbl %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with("p_"), ~ round(.x, 4)),
                      Chisq_interaction = round(Chisq_interaction, 2),
                      AIC = round(AIC, 2)),
      row.names = FALSE)
cat("\n  NOTE: Wald (Type III) and LRT p-values test the same hypothesis. Compare\n")
cat("  p_interaction_LRT with p_Wald_at_LRT_level (same cascade level, given in\n")
cat("  LRT_level). A large gap means the Wald approximation is unreliable at this n.\n")
cat(sprintf("  %s is a POST HOC sensitivity check, not a model choice.\n", sens_label))

readr::write_csv(pair_tbl, "outputs/tables/Scherber_biomass_pair_AMMOD_FAIRD.csv")
cat("[OK] Saved: outputs/tables/Scherber_biomass_pair_AMMOD_FAIRD.csv\n")

# --- AMMOD/FAIRD ratio: dual-axis scale factor and drift over time ------------
# Under the log link, parallel trajectories mean a CONSTANT AMMOD/FAIRD ratio.
# The daily pooled ratio below is the raw-scale reading of the interaction:
# a stable ratio is parallelism, a trend in the ratio is non-parallelism.
ratio_info <- ammod_faird_ratio(bio_overlap)
scale_factor_ammod <- signif(ratio_info$overall, 2)

cat("\n--- AMMOD/FAIRD biomass ratio over co-located site-days ---\n")
cat(sprintf("  Site-days used: %d\n", ratio_info$n_site_days))
cat(sprintf("  Overall ratio (sum AMMOD / sum FAIRD): %.2f\n", ratio_info$overall))
cat(sprintf("  Scale factor used for dual-axis figures: %s\n", scale_factor_ammod))
cat("\n  Daily pooled ratio (sites with both systems that day):\n")
print(ratio_info$daily %>%
        dplyr::mutate(dplyr::across(c(AMMOD_sum, FAIRD_sum), ~ round(.x, 1)),
                      ratio = round(ratio, 1)),
      n = Inf)

# Variant shown in the pair figure (Section 3). PAIRED is the default because
# both systems are observed on exactly the same site-days.
pair_for_plot         <- pair_paired
pair_for_plot_data    <- bio_FA_paired
pair_for_plot_variant <- "PAIRED"
if (is.null(pair_for_plot$model) && !is.null(pair_full$model)) {
  pair_for_plot         <- pair_full
  pair_for_plot_data    <- bio_FA_full
  pair_for_plot_variant <- "FULL"
  cat("[WARNING] PAIRED pair model not valid -> pair figure uses FULL.\n")
}


# ==============================================================================
# SECTION 3: TEMPORAL TREND PLOTS (linear y-axis only, no log axes)
# ------------------------------------------------------------------------------
# Figures:
#   1. Three systems, single linear axis (raw reference view).
#   2. Three systems, dual axis: AMMOD divided by scale_factor_ammod on the
#      right axis, FAIRD and ID on the left axis.
#   3. AMMOD vs FAIRD, dual axis, predictions from the PAIR model (Section 2e).
# The two former log10 figures (Scherber_temporal_biomass.png and
# Scherber_temporal_biomass_FA_only.png) are no longer produced.
# ==============================================================================

# predict_fixed(): fixed-effects predictions with a 95% CI built on the link
# scale and back-transformed. Block is a fixed term, so predictions are made at
# the reference block (BLOCK1 / Maize); this is stated in the subtitles.
predict_fixed <- function(m, dat, ref_block) {
  date_seq <- seq(min(dat$Date.num), max(dat$Date.num), length.out = 200)
  grid <- expand.grid(
    Date.num    = date_seq,
    Device_type = levels(dat$Device_type)
  ) %>%
    dplyr::mutate(
      Site   = NA,
      Device = NA,
      times  = NA,
      Block  = factor(ref_block, levels = levels(dat$Block))
    )
  grid$fit <- predict(m, newdata = grid, type = "response",
                      re.form = NA, se.fit = FALSE)
  pl <- predict(m, newdata = grid, type = "link", re.form = NA, se.fit = TRUE)
  linkinv    <- m$modelInfo$family$linkinv
  grid$lower <- linkinv(pl$fit - 1.96 * pl$se.fit)
  grid$upper <- linkinv(pl$fit + 1.96 * pl$se.fit)
  grid$Date  <- as.Date(grid$Date.num, origin = "1970-01-01")
  grid
}

# make_dual_plot(): linear dual-axis figure. AMMOD (observed, fit, CI) is
# divided by k and read on the right axis; other systems on the left axis.
# The view is capped at the rounded 98th percentile of the plotted points via
# coord_cartesian(), so nothing is dropped from lines or ribbons.
make_dual_plot <- function(pred, obs, k, title, subtitle, left_name) {
  obs_s <- obs %>%
    dplyr::mutate(y_plot = dplyr::if_else(Device_type == "AMMOD",
                                          Biomass_mg / k, Biomass_mg))
  pred_s <- pred %>%
    dplyr::mutate(
      fit_plot   = dplyr::if_else(Device_type == "AMMOD", fit   / k, fit),
      lower_plot = dplyr::if_else(Device_type == "AMMOD", lower / k, lower),
      upper_plot = dplyr::if_else(Device_type == "AMMOD", upper / k, upper)
    )
  y_cap <- nice_ceiling(stats::quantile(obs_s$y_plot, 0.98, na.rm = TRUE))

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_s,
      ggplot2::aes(x = Date, ymin = lower_plot, ymax = upper_plot, fill = Device_type),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = pred_s,
      ggplot2::aes(x = Date, y = fit_plot, color = Device_type),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = obs_s,
      ggplot2::aes(x = Date, y = y_plot, color = Device_type),
      alpha = 0.5, size = 1.8
    ) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_continuous(
      name     = left_name,
      labels   = scales::label_comma(),
      sec.axis = ggplot2::sec_axis(~ . * k,
                                   name   = sprintf("AMMOD Biomass (mg), axis = left x %s", k),
                                   labels = scales::label_comma())
    ) +
    ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
    ggplot2::labs(
      title = title, subtitle = subtitle, x = "Date",
      color = "Device_type", fill = "Device_type",
      caption = sprintf("View capped at %s (rounded 98th percentile of plotted points).",
                        scales::label_comma()(y_cap))
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle      = ggplot2::element_text(size = 9, color = "grey40"),
      plot.caption       = ggplot2::element_text(size = 8, color = "grey50"),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "top",
      axis.title.y.left  = ggplot2::element_text(color = "grey30"),
      axis.title.y.right = ggplot2::element_text(color = "red"),
      axis.text.y.right  = ggplot2::element_text(color = "red")
    )
}

if (is.null(model_for_plots)) {
  cat("\n=== TEMPORAL PLOTS SKIPPED (no valid model) ===\n")
} else {

  cat("\n\n=== GENERATING TEMPORAL PLOTS ===\n\n")
  cat(sprintf("Three-system figures generated from: %s\n", model_for_plots_name))

  ref_block <- levels(bio_overlap$Block)[1]
  pred_grid <- predict_fixed(model_for_plots, bio_overlap, ref_block)

  interaction_label <- ifelse(is.na(p_interact), "not available",
                              sprintf("%.4f", p_interact))
  subtitle_base <- sprintf(
    "Model: %s | predictions at %s | Date x Device_type p = %s (6 df)",
    model_for_plots_name, ref_block, interaction_label
  )

  # --- Plot 1: three systems, single linear axis -----------------------------
  y_cap_lin <- nice_ceiling(stats::quantile(bio_overlap$Biomass_mg, 0.98, na.rm = TRUE))

  p_linear <- ggplot2::ggplot() +
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
      ggplot2::aes(x = Date, y = Biomass_mg, color = Device_type),
      alpha = 0.5, size = 1.8
    ) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::coord_cartesian(ylim = c(0, y_cap_lin)) +
    ggplot2::labs(
      title    = "Temporal Biomass Trends by Device Type",
      subtitle = paste0(subtitle_base, " | linear y-axis"),
      x = "Date", y = "Biomass (mg)",
      color = "Device_type", fill = "Device_type",
      caption = sprintf("View capped at %s (rounded 98th percentile).",
                        scales::label_comma()(y_cap_lin))
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 9, color = "grey40"),
      plot.caption     = ggplot2::element_text(size = 8, color = "grey50"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )

  print(p_linear)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_linear.png", p_linear,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_linear.png\n")

  # --- Plot 2: three systems, dual axis --------------------------------------
  p_dual <- make_dual_plot(
    pred      = pred_grid,
    obs       = bio_overlap,
    k         = scale_factor_ammod,
    title     = "Temporal Biomass Trends by Device Type",
    subtitle  = sprintf("%s | AMMOD on right axis, scale factor %s (data-derived)",
                        subtitle_base, scale_factor_ammod),
    left_name = "FAIRD / ID Biomass (mg)"
  )
  print(p_dual)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_dualaxis.png", p_dual,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_dualaxis.png\n")
}

# --- Plot 3: AMMOD vs FAIRD, dual axis, from the PAIR model -------------------
if (is.null(pair_for_plot$model)) {
  cat("\n[SKIP] Pair figure: no valid pair model.\n")
} else {
  pair_row <- pair_tbl %>% dplyr::filter(Variant == pair_for_plot_variant)
  pred_pair <- predict_fixed(pair_for_plot$model, pair_for_plot_data,
                             levels(pair_for_plot_data$Block)[1])
  p_pair <- make_dual_plot(
    pred      = pred_pair,
    obs       = pair_for_plot_data,
    k         = scale_factor_ammod,
    title     = "Temporal Biomass Trends: AMMOD vs FAIRD (pair model)",
    subtitle  = sprintf(
      "Pair model %s, variant %s | predictions at %s | Date x Device_type p = %.4f (Wald, %d df), %.4f (LRT, level %s)",
      pair_for_plot$spec$name, pair_for_plot_variant,
      levels(pair_for_plot_data$Block)[1],
      pair_row$p_interaction_Wald, as.integer(pair_row$df_interaction),
      pair_row$p_interaction_LRT, pair_row$LRT_level
    ),
    left_name = "FAIRD Biomass (mg)"
  )
  print(p_pair)
  ggplot2::ggsave("outputs/figures/Scherber_temporal_biomass_FA_dualaxis.png", p_pair,
                  width = 10, height = 7, dpi = 300, bg = "white")
  cat("[OK] Saved: outputs/figures/Scherber_temporal_biomass_FA_dualaxis.png\n")
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
