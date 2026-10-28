# ==============================================================================
# Scherber_temporal_biomass_block_comparison.R
# ------------------------------------------------------------------------------
# Purpose: quantify how treating Block as a RANDOM effect vs a FIXED effect
# changes the biomass temporal model, with particular attention to the term
# that matters for the manuscript: the Date x Device_type interaction
# (trajectory parallelism between systems).
#
# Motivation. Block has only 2 levels (BLOCK1 = Maize, BLOCK2 = Meadow). A
# variance component estimated from 2 levels is poorly identified and tends to
# collapse to a singular fit, which mis-calibrates the standard errors of every
# other term. As a 2-level fixed contrast it is estimated without that problem.
# This script fits both specifications through identical code so the ANOVA
# tables and interaction p-values are directly comparable, and reports the
# Block variance in the random version so the collapse (if any) is visible.
#
# Models (both Tweedie, both with AR1 within device):
#   B_ranef : ... + (1 | Block/Site/Device) + ar1(times + 0 | Device)   [current]
#   B_fixef : ... + Block + (1 | Site/Device) + ar1(times + 0 | Device) [proposed]
# Each is validated with tmb_is_valid(); on failure it falls back to its own
# no-AR1 version, and the fallback is logged.
#
# No figures are produced here: this is a specification comparison, not a
# presentation output. Console log + one comparison table only.
# ==============================================================================

library(tidyverse)
library(glmmTMB)
library(car)
library(splines)

source("R/00_load_data.R")

dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console", showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/Scherber_temporal_biomass_block_comparison_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("BIOMASS TEMPORAL MODEL: BLOCK RANDOM vs BLOCK FIXED\n")
cat("Identical fixed structure and AR1; only the role of Block differs.\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# tmb_is_valid(): strict convergence check for a glmmTMB fit.
# (Reused verbatim from Scherber_temporal_abundance_FAIRD_ID.R so the validity
# criterion is identical across the project. A model can be returned without
# erroring yet be numerically degenerate: non-positive-definite Hessian, NA AIC.
# tryCatch does not catch this, so it is checked explicitly.)
# ==============================================================================
tmb_is_valid <- function(m) {
  if (is.null(m)) return(FALSE)
  pd_hess <- isTRUE(tryCatch(m$sdr$pdHess, error = function(e) FALSE))
  aic_ok  <- tryCatch(is.finite(AIC(m)), error = function(e) FALSE)
  pd_hess && aic_ok
}


# ==============================================================================
# SECTION 1: DATA PREPARATION
# (identical to Scherber_temporal_biomass.R so results are comparable)
# ==============================================================================

cat("=== DATA PREPARATION ===\n\n")

# FAIRD and ID daily biomass from etraps_daily (zero-capture days already
# filled as real zeros upstream -> Tweedie family is appropriate).
etraps_bio <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, biomass_mg) %>%
  dplyr::rename(Biomass_mg = biomass_mg)

# AMMOD daily biomass — rename to match structure.
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

# Restrict to overlapping date range (AMMOD only has Batch 1: 12 days).
date_range  <- range(ammod_bio$Date)
bio_overlap <- bio_combined %>%
  dplyr::filter(Date >= date_range[1] & Date <= date_range[2])

# Design factors.
# Block: Maize (Sites 1-2) vs Meadow (Sites 3-4), treated as a blocking factor.
bio_overlap$Block <- factor(
  ifelse(bio_overlap$Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")
)
# times: Date as a factor, required for the AR1 structure in glmmTMB.
bio_overlap$times <- factor(bio_overlap$Date)

cat(sprintf("Overlap period (%s to %s): %d observations\n",
            date_range[1], date_range[2], nrow(bio_overlap)))
cat(sprintf("Distinct days: %d\n", dplyr::n_distinct(bio_overlap$Date)))
cat("Blocks: BLOCK1 = Site1+Site2 (Maize), BLOCK2 = Site3+Site4 (Meadow)\n")
cat(sprintf("Block levels: %d | Site levels: %d | Device levels: %d\n",
            dplyr::n_distinct(bio_overlap$Block),
            dplyr::n_distinct(bio_overlap$Site),
            dplyr::n_distinct(bio_overlap$Device)))

cat("\nObservations per Device (overlap period):\n")
print(bio_overlap %>%
        dplyr::group_by(Block, Device_type, Device, Site) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop"))


# ==============================================================================
# SECTION 2: FIT ONE SPECIFICATION (helper)
# ------------------------------------------------------------------------------
# Fits the AR1 model for a given Block role; on validity failure, falls back to
# the no-AR1 version of the same specification. Returns the chosen fit, its
# label, and the validity flags, so the caller can report what actually fed the
# comparison.
# ==============================================================================

fit_spec <- function(block_role, data) {

  stopifnot(block_role %in% c("random", "fixed"))

  # Build the two candidate formulas (AR1 and no-AR1) for this Block role.
  if (block_role == "random") {
    f_ar1  <- Biomass_mg ~ bs(Date.num, 3) * Device_type +
      (1 | Block / Site / Device) + ar1(times + 0 | Device)
    f_noar <- Biomass_mg ~ bs(Date.num, 3) * Device_type +
      (1 | Block / Site / Device)
    tag <- "Block RANDOM (1|Block/Site/Device)"
  } else {
    f_ar1  <- Biomass_mg ~ bs(Date.num, 3) * Device_type + Block +
      (1 | Site / Device) + ar1(times + 0 | Device)
    f_noar <- Biomass_mg ~ bs(Date.num, 3) * Device_type + Block +
      (1 | Site / Device)
    tag <- "Block FIXED (+ Block, 1|Site/Device)"
  }

  cat(sprintf("\n--- Fitting: %s ---\n", tag))

  # Attempt the AR1 model first.
  m_ar1 <- tryCatch(
    glmmTMB::glmmTMB(f_ar1, family = glmmTMB::tweedie(), data = data),
    error = function(e) {
      cat("  [FAILED] AR1 model errored:\n    ", conditionMessage(e), "\n")
      NULL
    }
  )

  if (tmb_is_valid(m_ar1)) {
    cat("  [OK] AR1 model valid (positive-definite Hessian, finite AIC).\n")
    return(list(fit = m_ar1, label = paste(tag, "+ AR1"),
                block_role = block_role, ar1 = TRUE))
  }
  cat("  [INVALID] AR1 model not valid -> falling back to no-AR1 version.\n")

  # Fall back to the no-AR1 version of the same Block role.
  m_noar <- tryCatch(
    glmmTMB::glmmTMB(f_noar, family = glmmTMB::tweedie(), data = data),
    error = function(e) {
      cat("  [FAILED] no-AR1 model errored:\n    ", conditionMessage(e), "\n")
      NULL
    }
  )

  if (tmb_is_valid(m_noar)) {
    cat("  [OK] no-AR1 model valid.\n")
    return(list(fit = m_noar, label = paste(tag, "no AR1"),
                block_role = block_role, ar1 = FALSE))
  }

  cat("  [INVALID] no-AR1 model also not valid -> no usable fit for this role.\n")
  list(fit = NULL, label = paste(tag, "NO VALID FIT"),
       block_role = block_role, ar1 = NA)
}


# ==============================================================================
# SECTION 3: FIT BOTH SPECIFICATIONS
# ==============================================================================

cat("\n=== MODEL FITTING ===\n")

res_ranef <- fit_spec("random", bio_overlap)
res_fixef <- fit_spec("fixed",  bio_overlap)


# ==============================================================================
# SECTION 4: TYPE III ANOVA + KEY-TERM EXTRACTION
# ------------------------------------------------------------------------------
# For each valid fit, pull the Type III p-values for the interaction, the spline
# and Device_type. The interaction row is the manuscript-relevant one.
# ==============================================================================

interaction_term <- "bs(Date.num, 3):Device_type"
spline_term      <- "bs(Date.num, 3)"
device_term      <- "Device_type"

extract_terms <- function(res) {
  out <- list(label = res$label, ar1 = res$ar1,
              p_interaction = NA_real_, p_spline = NA_real_,
              p_device = NA_real_, p_block = NA_real_,
              block_var = NA_real_)

  if (is.null(res$fit)) return(out)

  av <- tryCatch(car::Anova(res$fit, type = "III"),
                 error = function(e) {
                   cat(sprintf("  [SKIP] Type III ANOVA failed for %s: %s\n",
                               res$label, conditionMessage(e)))
                   NULL
                 })

  if (!is.null(av)) {
    rn <- rownames(av)
    if (interaction_term %in% rn) out$p_interaction <- av[interaction_term, "Pr(>Chisq)"]
    if (spline_term      %in% rn) out$p_spline      <- av[spline_term,      "Pr(>Chisq)"]
    if (device_term      %in% rn) out$p_device      <- av[device_term,      "Pr(>Chisq)"]
    # Block p-value exists only in the fixed specification.
    if ("Block"          %in% rn) out$p_block       <- av["Block",          "Pr(>Chisq)"]
  }

  # Block variance exists only in the random specification. Extract it so a
  # collapse toward zero (the singular-fit symptom) is visible in the output.
  if (res$block_role == "random") {
    vc <- tryCatch(glmmTMB::VarCorr(res$fit)$cond, error = function(e) NULL)
    if (!is.null(vc) && "Block" %in% names(vc)) {
      out$block_var <- as.numeric(attr(vc$Block, "stddev"))^2
    }
  }

  out
}

cat("\n=== TYPE III ANOVA (both specifications) ===\n")

if (!is.null(res_ranef$fit)) {
  cat(sprintf("\n--- %s ---\n", res_ranef$label))
  print(car::Anova(res_ranef$fit, type = "III"))
  cat("\nVariance components (random specification):\n")
  print(glmmTMB::VarCorr(res_ranef$fit))
}

if (!is.null(res_fixef$fit)) {
  cat(sprintf("\n--- %s ---\n", res_fixef$label))
  print(car::Anova(res_fixef$fit, type = "III"))
  cat("\nFixed-effect coefficients incl. Block (fixed specification):\n")
  print(glmmTMB::fixef(res_fixef$fit)$cond)
}

terms_ranef <- extract_terms(res_ranef)
terms_fixef <- extract_terms(res_fixef)


# ==============================================================================
# SECTION 5: SIDE-BY-SIDE COMPARISON TABLE
# ==============================================================================

cat("\n=== SIDE-BY-SIDE COMPARISON ===\n\n")

fmt_p <- function(p) ifelse(is.na(p), "NA", sprintf("%.4f", p))

comparison <- data.frame(
  Quantity = c(
    "Model label",
    "AR1 retained",
    "p: Date x Device_type interaction",
    "p: spline bs(Date.num,3)",
    "p: Device_type",
    "p: Block (fixed only)",
    "Block variance (random only)",
    "AIC"
  ),
  Block_random = c(
    res_ranef$label,
    ifelse(is.na(terms_ranef$ar1), "NA", as.character(terms_ranef$ar1)),
    fmt_p(terms_ranef$p_interaction),
    fmt_p(terms_ranef$p_spline),
    fmt_p(terms_ranef$p_device),
    "— (not in model)",
    ifelse(is.na(terms_ranef$block_var), "NA",
           sprintf("%.4g", terms_ranef$block_var)),
    ifelse(is.null(res_ranef$fit), "NA", sprintf("%.2f", AIC(res_ranef$fit)))
  ),
  Block_fixed = c(
    res_fixef$label,
    ifelse(is.na(terms_fixef$ar1), "NA", as.character(terms_fixef$ar1)),
    fmt_p(terms_fixef$p_interaction),
    fmt_p(terms_fixef$p_spline),
    fmt_p(terms_fixef$p_device),
    fmt_p(terms_fixef$p_block),
    "— (not in model)",
    ifelse(is.null(res_fixef$fit), "NA", sprintf("%.2f", AIC(res_fixef$fit)))
  ),
  stringsAsFactors = FALSE
)

print(comparison, row.names = FALSE, right = FALSE)

readr::write_csv(comparison,
                 "outputs/tables/Scherber_biomass_block_comparison.csv")
cat("\n[OK] Saved: outputs/tables/Scherber_biomass_block_comparison.csv\n")


# ==============================================================================
# SECTION 6: READING GUIDE (printed, not interpreted for you)
# ==============================================================================

cat("\n=== HOW TO READ THIS ===\n")
cat(paste(
  "- The interaction p-value is the manuscript-relevant number. If it is",
  "  stable across the two specifications, the trajectory-parallelism",
  "  conclusion does not depend on the Block modelling choice.",
  "- Block variance (random spec): if it is ~0 or tiny relative to Site/Device,",
  "  the 2-level random effect has collapsed (singular fit) and the fixed",
  "  specification is the better-calibrated one.",
  "- Device_type p-value is the term most likely to move, because Block and",
  "  Device_type share mean-level variation under an unbalanced habitat layout",
  "  (AMMOD4/Meadow absent, Site3/Meadow reduced).",
  "- AIC is NOT comparable between the two specifications as a model-selection",
  "  criterion here, because a fixed and a random treatment of the same factor",
  "  are not nested in the usual sense; use it only as a rough sanity check.",
  sep = "\n"), "\n")

cat("\n==============================================================================\n")
cat("COMPARISON COMPLETE\n")
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================================\n")

while (sink.number() > 0) sink()
