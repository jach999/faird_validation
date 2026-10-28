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
#   - Pair models (Section 2e), one generic function for AMMOD-FAIRD,
#     AMMOD-ID and FAIRD-ID: same specification as B1 on each two-system
#     subset, so the Date x Device_type interaction is tested with 3 df instead
#     of the 6 df of the omnibus test. Primary pairing PAIRED_BLOCK (block-days
#     with both systems), sensitivities PAIRED_SITE and FULL; Holm over the
#     three primary 12-day LRTs. FAIRD-ID is also fitted over 22 days, full vs
#     adverse window excluded (complementary). Wald and LRT for every model,
#     the LRT falling back down the cascade when the additive model is not
#     valid at the reported level.
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
cat("Pair models AMMOD-FAIRD, AMMOD-ID, FAIRD-ID (3-df interaction) in Section 2e.\n")
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
# SECTION 2e: PAIR MODELS, AMMOD-FAIRD, AMMOD-ID AND FAIRD-ID
# ------------------------------------------------------------------------------
# The omnibus test above (6 df) says the three trajectories are not all
# parallel, but not WHICH pair departs. Each pair is tested on its own with the
# specification of B1 (Block fixed, (1|Site/Device), AR1, Tweedie), so the
# Date x Device_type interaction has 3 df. One generic function (run_pair)
# serves the three pairs; there is no pair-specific code.
#
# Pairing criterion (pair_by), decided 2026-10-01:
#   PAIRED_BLOCK (PRIMARY): only block-days on which BOTH systems have a value
#        in the same block, even at different sites. Follows the block framing
#        requested by Scherber: sites within a block are pseudoreplicates, and
#        Device as a random effect absorbs unequal device numbers per type.
#        This gives AMMOD-ID a Meadow component (AMMOD3 vs ID4).
#   PAIRED_SITE (sensitivity): only co-located site-days (former PAIRED).
#   FULL (sensitivity): all observations of the two systems.
# A sensitivity variant whose data are identical to PAIRED_BLOCK is not
# refitted (e.g. FAIRD-ID FULL).
#
# Two analysis blocks:
#   A_12d: 12-day AMMOD window, full ID series, the three pairs. Holm
#          correction over the three PAIRED_BLOCK LRT p-values (primary
#          family); a Wald-based Holm column is given alongside.
#   B_22d: FAIRD-ID over 22 days, full series vs adverse window excluded,
#          PAIRED_BLOCK. Complementary analysis, parallel to the abundance
#          analysis. Excluded dates come from the validity calendar
#          (Field_Valid_Operational_ID), not from hardcoded dates, and are
#          removed for BOTH systems, so the exclusion is symmetric. Inside the
#          12-day window the exclusion would leave 7 days (5 + 2), too few for
#          a spline, so it is applied to the 22-day series only.
#
# Time levels: after subsetting, `times` is rebuilt with ALL dates of the
# analysis window as levels. droplevels() would otherwise delete dates absent
# for every device (e.g. the excluded window), and ar1() would then treat the
# days on either side of the gap as consecutive without any warning. Unobserved
# levels are latent AR1 states with no data, which is the intended treatment.
# ==============================================================================

cat("\n\n=== PAIR MODELS: AMMOD-FAIRD, AMMOD-ID, FAIRD-ID ===\n")

pairing_label <- c(Block = "PAIRED_BLOCK", Site = "PAIRED_SITE", none = "FULL")

# make_pair_specs(): cascade of specifications for one pair. With two Block
# levels it is the cascade used so far (AR1 -> no AR1 -> no Block), formula
# strings unchanged. With a single Block level (e.g. AMMOD-ID PAIRED_SITE,
# Maize only) Block is not estimable and is left out of every level, so AR1 is
# not lost for a reason unrelated to AR1.
make_pair_specs <- function(has_block) {
  re <- "(1 | Site / Device)"
  ar <- "ar1(times + 0 | Device)"
  mk <- function(rhs, name) list(
    full = paste0("Biomass_mg ~ bs(Date.num, 3) * Device_type", rhs),
    add  = paste0("Biomass_mg ~ bs(Date.num, 3) + Device_type", rhs),
    name = name
  )
  if (has_block) {
    list(
      AR1    = mk(paste0(" + Block + ", re, " + ", ar), "Block fixed + (1|Site/Device) + AR1"),
      noAR1  = mk(paste0(" + Block + ", re),            "Block fixed + (1|Site/Device), no AR1"),
      simple = mk(paste0(" + ", re),                    "(1|Site/Device), no Block, no AR1")
    )
  } else {
    list(
      AR1    = mk(paste0(" + ", re, " + ", ar), "(1|Site/Device) + AR1, single Block"),
      simple = mk(paste0(" + ", re),            "(1|Site/Device), single Block, no AR1")
    )
  }
}

# build_pair_data(): subset to two systems, optionally drop dates, apply the
# pairing criterion, and rebuild `times` with the full window levels.
build_pair_data <- function(dat, types, pair_by, time_levels, drop_dates = NULL) {
  out <- dat %>%
    dplyr::filter(Device_type %in% types, !is.na(Biomass_mg))
  if (!is.null(drop_dates)) {
    out <- out %>% dplyr::filter(!(Date %in% drop_dates))
  }
  if (pair_by %in% c("Site", "Block")) {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(pair_by)), Date) %>%
      dplyr::filter(dplyr::n_distinct(Device_type) == 2) %>%
      dplyr::ungroup()
  }
  out <- droplevels(out)
  out$times <- factor(as.character(out$Date), levels = time_levels)
  out
}

pair_keys <- function(d) sort(paste(d$Device, d$Date))

fit_pair_cascade <- function(dat, label, specs) {
  cat(sprintf("\n##### %s #####\n", label))
  for (lvl in names(specs)) {
    spec <- specs[[lvl]]
    m <- fit_tmb(stats::as.formula(spec$full), dat,
                 sprintf("%s | %s", label, spec$name))
    if (tmb_is_valid(m)) {
      if (lvl != names(specs)[1]) {
        cat(sprintf("  [WARNING] higher cascade levels invalid -> using '%s'.\n", lvl))
      }
      return(list(model = m, level = lvl, spec = spec, specs = specs))
    }
  }
  cat("  [ERROR] no valid fit anywhere in the cascade.\n")
  list(model = NULL, level = NA_character_, spec = NULL, specs = specs)
}

# pair_lrt(): likelihood-ratio test of the interaction (full vs additive fixed
# part, identical random and AR1 structure). Non-positive-definite Hessian
# warnings printed by glmmTMB during this step come from fits that
# tmb_is_valid() then rejects; they are expected and handled. Tried first at
# the cascade level of the reported model; if the full or the additive fit is
# not valid there, it moves down the cascade and refits BOTH models at the
# lower level, so the two models compared always share their structure. The
# Wald p-value at the level where the LRT was computed is returned too.
pair_lrt <- function(res, dat, label) {
  empty <- list(p_lrt = NA_real_, p_wald = NA_real_, level = NA_character_)
  if (is.null(res$model)) return(empty)
  lvls  <- names(res$specs)
  start <- match(res$level, lvls)
  for (lvl in lvls[start:length(lvls)]) {
    spec   <- res$specs[[lvl]]
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
      # argument in each model call, and the fits are created by different
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
      cat(sprintf("\n--- LRT: interaction vs additive, %s, level '%s' ---\n", label, lvl))
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

summarise_pair <- function(res, dat, label) {
  out <- data.frame(
    Variant             = label,
    n_obs               = nrow(dat),
    n_days              = dplyr::n_distinct(dat$Date),
    n_devices           = dplyr::n_distinct(dat$Device),
    Model               = ifelse(is.null(res$spec), NA_character_, res$spec$name),
    Chisq_interaction   = NA_real_,
    df_interaction      = NA_real_,
    p_interaction_Wald  = NA_real_,
    p_Device_type       = NA_real_,
    p_spline            = NA_real_,
    p_Block             = NA_real_,
    p_interaction_LRT   = NA_real_,
    LRT_level           = NA_character_,
    p_Wald_at_LRT_level = NA_real_,
    AIC                 = NA_real_,
    stringsAsFactors    = FALSE
  )
  if (is.null(res$model)) return(out)
  m <- res$model
  out$AIC <- AIC(m)

  av <- tryCatch(car::Anova(m, type = "III"), error = function(e) {
    cat("  [SKIP] Type III ANOVA failed:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(av)) {
    cat(sprintf("\n--- Type III ANOVA, %s ---\n", label)); print(av)
    if (interaction_term %in% rownames(av)) {
      out$Chisq_interaction  <- av[interaction_term, "Chisq"]
      out$df_interaction     <- av[interaction_term, "Df"]
      out$p_interaction_Wald <- av[interaction_term, "Pr(>Chisq)"]
    }
    if ("Device_type" %in% rownames(av))     out$p_Device_type <- av["Device_type", "Pr(>Chisq)"]
    if ("bs(Date.num, 3)" %in% rownames(av)) out$p_spline      <- av["bs(Date.num, 3)", "Pr(>Chisq)"]
    if ("Block" %in% rownames(av))           out$p_Block       <- av["Block", "Pr(>Chisq)"]
  }

  lrt_res <- pair_lrt(res, dat, label)
  out$p_interaction_LRT   <- lrt_res$p_lrt
  out$LRT_level           <- lrt_res$level
  out$p_Wald_at_LRT_level <- lrt_res$p_wald

  cat(sprintf("\n--- Fixed-effect coefficients, %s ---\n", label))
  print(glmmTMB::fixef(m)$cond)
  cat(sprintf("\n--- Variance components, %s ---\n", label))
  print(glmmTMB::VarCorr(m))
  out
}

# run_pair(): build data, fit the cascade, summarise. One call per
# pair x analysis x window x pairing.
run_pair <- function(types, dat, time_levels, pair_by, analysis, window_label,
                     role, drop_dates = NULL) {
  pair_name <- paste(types, collapse = "-")
  label <- sprintf("%s | %s | %s | %s", pair_name, analysis, window_label,
                   pairing_label[[pair_by]])
  pdat <- build_pair_data(dat, types, pair_by, time_levels, drop_dates)
  cat(sprintf("\n\n>>> PAIR %s\n", label))
  cat(sprintf("  %d observations, %d days (%s to %s)\n", nrow(pdat),
              dplyr::n_distinct(pdat$Date), min(pdat$Date), max(pdat$Date)))
  cat("  Observations per device:\n")
  print(pdat %>% dplyr::count(Block, Device_type, Device, Site, name = "n"))
  has_block <- dplyr::n_distinct(pdat$Block) == 2
  if (!has_block) {
    cat("  [NOTE] single Block level in these data -> Block left out of the model.\n")
  }
  res <- fit_pair_cascade(pdat, label, make_pair_specs(has_block))
  row <- summarise_pair(res, pdat, label) %>%
    dplyr::mutate(Pair = pair_name, Analysis = analysis, Window = window_label,
                  Pairing = pairing_label[[pair_by]], Role = role, .before = 1)
  list(res = res, data = pdat, row = row, types = types, label = label,
       analysis = analysis, window = window_label, role = role)
}

# system_ratio(): ratio of DEVICE MEANS of two systems over the pair data.
# Means rather than sums, because under PAIRED_BLOCK the device numbers per
# type can differ within a block-day (Meadow: AMMOD3 vs FAIRD3 + FAIRD4).
# The daily series is the raw-scale reading of the interaction: under the log
# link, parallel trajectories mean a constant ratio.
system_ratio <- function(pdat, num, den) {
  daily <- pdat %>%
    dplyr::filter(Device_type %in% c(num, den)) %>%
    dplyr::mutate(Device_type = as.character(Device_type)) %>%
    dplyr::group_by(Date, Device_type) %>%
    dplyr::summarise(mean_mg = mean(Biomass_mg), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Device_type, values_from = mean_mg) %>%
    dplyr::filter(!is.na(.data[[num]]), !is.na(.data[[den]])) %>%
    dplyr::mutate(ratio = dplyr::if_else(.data[[den]] > 0,
                                         .data[[num]] / .data[[den]], NA_real_))
  overall <- mean(pdat$Biomass_mg[pdat$Device_type == num]) /
    mean(pdat$Biomass_mg[pdat$Device_type == den])
  list(overall = overall, daily = daily)
}

# --- Block A: 12-day window, full ID series -----------------------------------
pair_types     <- list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))
time_levels_12 <- levels(bio_overlap$times)
pair_runs      <- list()

for (tp in pair_types) {
  primary <- run_pair(tp, bio_overlap, time_levels_12, "Block",
                      "A_12d", "full", "primary")
  pair_runs[[primary$label]] <- primary
  if (identical(tp, c("AMMOD", "FAIRD"))) af_primary <- primary
  for (pb in c("Site", "none")) {
    cand <- build_pair_data(bio_overlap, tp, pb, time_levels_12)
    if (identical(pair_keys(cand), pair_keys(primary$data))) {
      cat(sprintf("\n[SKIP] %s, %s: data identical to PAIRED_BLOCK, not refitted.\n",
                  paste(tp, collapse = "-"), pairing_label[[pb]]))
      next
    }
    r <- run_pair(tp, bio_overlap, time_levels_12, pb, "A_12d", "full", "sensitivity")
    pair_runs[[r$label]] <- r
  }
}

# --- Sensitivity (POST HOC): AMMOD-FAIRD primary without the first days -------
# The per-site diagnostic showed the AMMOD/FAIRD ratio is high on the first
# days at every site and then drops. Post hoc, chosen after seeing the data,
# reported as such. Applied to AMMOD-FAIRD only (not to the new pairs).
n_drop_initial <- 3
initial_dates  <- sort(unique(af_primary$data$Date))[seq_len(n_drop_initial)]
cat(sprintf("\nPost hoc sensitivity: dates removed = %s\n",
            paste(format(initial_dates), collapse = ", ")))
r <- run_pair(c("AMMOD", "FAIRD"), bio_overlap, time_levels_12, "Block", "A_12d",
              sprintf("drop_first%d", n_drop_initial), "sensitivity_posthoc",
              drop_dates = initial_dates)
pair_runs[[r$label]] <- r

# --- Block B: FAIRD-ID over 22 days, full vs adverse window excluded ----------
adverse_dates <- data_list$effort_base %>%
  dplyr::filter(Field_Valid, !Field_Valid_Operational_ID) %>%
  dplyr::distinct(Date) %>%
  dplyr::arrange(Date) %>%
  dplyr::pull(Date)
cat(sprintf("\nAdverse window from the validity calendar (%d dates): %s\n",
            length(adverse_dates), paste(format(adverse_dates), collapse = ", ")))

bio_22 <- bio_combined %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::mutate(
    Block = factor(ifelse(Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")),
    times = factor(Date)
  )
time_levels_22 <- levels(bio_22$times)

r <- run_pair(c("FAIRD", "ID"), bio_22, time_levels_22, "Block",
              "B_22d", "full", "complementary")
pair_runs[[r$label]] <- r
r <- run_pair(c("FAIRD", "ID"), bio_22, time_levels_22, "Block",
              "B_22d", "window_excluded", "complementary",
              drop_dates = adverse_dates)
pair_runs[[r$label]] <- r

# --- Results table, Holm correction over the primary family -------------------
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

# Omnibus row for side-by-side reading (6-df test from the plotting model).
omnibus_row <- data.frame(
  Pair               = "ALL (3 systems)",
  Analysis           = "A_12d",
  Window             = "full",
  Pairing            = "FULL",
  Role               = "omnibus",
  Variant            = "OMNIBUS (3 systems)",
  n_obs              = nrow(bio_overlap),
  n_days             = dplyr::n_distinct(bio_overlap$Date),
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
cat("DATE x DEVICE_TYPE INTERACTION: OMNIBUS AND PAIR MODELS\n")
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
cat("  - AMMOD-ID PAIRED_SITE is Maize only (Sites 1-2, Block not estimable).\n")
cat("  - drop_first rows are a POST HOC sensitivity check, not a model choice.\n")
cat("  - B_22d rows are complementary and outside the Holm family.\n")

readr::write_csv(pair_tbl, "outputs/tables/Scherber_biomass_pairs.csv")
cat("[OK] Saved: outputs/tables/Scherber_biomass_pairs.csv\n")

# --- Ratios per primary pair: drift over time ---------------------------------
cat("\n--- Daily ratio of device means per pair (PAIRED_BLOCK data) ---\n")
for (r in pair_runs) {
  if (!(r$role %in% c("primary", "complementary"))) next
  rr <- system_ratio(r$data, r$types[1], r$types[2])
  cat(sprintf("\n  %s: overall ratio %s/%s = %.2f, Spearman rho(ratio, date) = %.2f\n",
              r$label, r$types[1], r$types[2], rr$overall,
              suppressWarnings(stats::cor(as.numeric(rr$daily$Date), rr$daily$ratio,
                                          method = "spearman", use = "complete.obs"))))
  print(rr$daily %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 2))),
        n = Inf)
}

# --- AMMOD/FAIRD ratio for the three-system dual-axis figure (unchanged) ------
ratio_info <- ammod_faird_ratio(bio_overlap)
scale_factor_ammod <- signif(ratio_info$overall, 2)

cat("\n--- AMMOD/FAIRD biomass ratio over co-located site-days ---\n")
cat(sprintf("  Site-days used: %d\n", ratio_info$n_site_days))
cat(sprintf("  Overall ratio (sum AMMOD / sum FAIRD): %.2f\n", ratio_info$overall))
cat(sprintf("  Scale factor used for the three-system dual-axis figure: %s\n",
            scale_factor_ammod))

# ==============================================================================
# SECTION 3: TEMPORAL TREND PLOTS (linear y-axis only, no log axes)
# ------------------------------------------------------------------------------
# Figures:
#   1. Three systems, single linear axis (raw reference view).
#   2. Three systems, dual axis: AMMOD divided by scale_factor_ammod on the
#      right axis, FAIRD and ID on the left axis.
#   3. One figure per pair model (Section 2e): the three primary A_12d pairs
#      and the two complementary B_22d FAIRD-ID runs. Dual axis for pairs with
#      AMMOD (data-derived factor), single linear axis for FAIRD-ID.
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
make_dual_plot <- function(pred, obs, k, title, subtitle, left_name, dual = TRUE) {
  # dual = FALSE: single linear axis, no AMMOD rescaling (pairs without AMMOD).
  if (!dual) k <- 1
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
      sec.axis = if (dual) ggplot2::sec_axis(~ . * k,
                                   name   = sprintf("AMMOD Biomass (mg), axis = left x %s", k),
                                   labels = scales::label_comma()) else ggplot2::waiver()
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

# --- Plot 3: one figure per pair model ----------------------------------------
# Primary pair models of block A (PAIRED_BLOCK, 12 days) and the two
# complementary FAIRD-ID models of block B (22 days, full vs window excluded).
# Pairs with AMMOD use a dual axis whose scale factor is the ratio of device
# means over the pair data (system_ratio), derived from the data. FAIRD-ID is
# drawn on a single linear axis: both systems are e-traps in mg.
fig_runs <- Filter(function(r) r$role %in% c("primary", "complementary"), pair_runs)

for (r in fig_runs) {
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
    title     = sprintf("Temporal Biomass Trends: %s (pair model, %s, %s)",
                        sub("-", " vs ", row$Pair), row$Analysis, row$Window),
    subtitle  = sprintf(
      "%s | %s | predictions at %s | Date x Device_type p = %.4f (Wald, %s df), %.4f (LRT, level %s)%s%s",
      r$res$spec$name, row$Pairing, ref_block,
      row$p_interaction_Wald, row$df_interaction,
      row$p_interaction_LRT, row$LRT_level, window_note, axis_note
    ),
    left_name = if (has_ammod) sprintf("%s Biomass (mg)", other) else "Biomass (mg)",
    dual      = has_ammod
  )
  print(p_pair)
  fig_file <- sprintf("outputs/figures/Scherber_temporal_biomass_pair_%s_%s_%s.png",
                      gsub("-", "_", row$Pair), row$Analysis, row$Window)
  ggplot2::ggsave(fig_file, p_pair, width = 10, height = 7, dpi = 300, bg = "white")
  cat(sprintf("[OK] Saved: %s\n", fig_file))
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
