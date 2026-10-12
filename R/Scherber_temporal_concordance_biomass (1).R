# ==============================================================================
# Scherber_temporal_concordance_biomass.R
# System-level temporal concordance of daily biomass (Stage 3, part B)
# ------------------------------------------------------------------------------
# Question: do the systems record the same temporal signal, and how does the
# agreement between systems compare with the consistency WITHIN each system?
# Between-system parallelism is tested in Scherber_temporal_biomass_5.R
# (constant ratio, log link; constant difference, identity link).
#
# Internal consistency uses the SAME GLMM as the between-system pair models
# (Tweedie, log link, bs(Date.num, 3) * Unit + Block, (1|Site/Device), AR1),
# with the two halves of a system (or its two Maize devices) as the units.
# For every comparison, internal and between systems, the script reports:
#   PARALLELISM (constant ratio): LRT of the Date x Unit interaction, at the
#     reported cascade level, moving down the cascade only if the additive
#     model is not valid there (same logic as pair_lrt() in version 5). For
#     the between-system rows this reproduces version 5 and serves as a
#     cross-check; for the internal rows it is the within-system reference,
#     i.e. whether two halves of one system are proportional to each other.
#     Unadjusted p-values: the internal tests are references, not a family.
#
# Unit of analysis: the daily SERIES OF A SYSTEM (or of one half of a system
# for internal consistency), not the single device, except where a system has
# only one device per block (see internal references below).
#
# Components, reported separately:
#   LOW frequency (shared shape): Pearson r between the two fitted unit curves
#     of the pair model, raw scale (mg), at the observed dates. r is invariant
#     to a scale factor and to a shift: it measures "same temporal shape up to
#     a linear transformation". Estimation only, no p-value: two smooth
#     single-peak curves correlate almost by construction, so the value is
#     read against the internal references and across pairs. 95% interval
#     from simulation of the fixed coefficients, N(beta_hat, vcov), of the
#     reported model (AR1 when valid). Reading rule: each draw adds estimation
#     noise to the curves, and noise lowers a correlation, so the simulated
#     distribution is shifted downwards and the lower limit is pessimistic. The
#     interval can support a claim of shared shape (lower limit still high),
#     not a claim of absent shared shape.
#   HIGH frequency (shared day-to-day fluctuation): Pearson residuals of each
#     observation relative to its own unit curve, from the model WITHOUT AR1
#     (fixed part plus Site/Device intercepts). In the pair models the AR1
#     term behaves as a device-day effect (rho close to 0), so its BLUPs are
#     the day-to-day deviations themselves and must not be removed.
#     Residuals are averaged within block, then across blocks, giving one
#     series per unit; r at lag 0 between the two series (lag 0 is justified:
#     all systems share the calendar-day boundary). Estimation only:
#     approximate 95% CI by Fisher z with n = number of dates. The interval
#     ignores the degrees of freedom used by the spline and any residual
#     autocorrelation, so it is indicative rather than exact.
#   TOTAL (anchor, descriptive): r between the daily unit means.
#
# Comparisons, fixed before running:
#   Between systems, PAIRED_BLOCK (block-days with both units):
#     AMMOD-FAIRD, AMMOD-ID, FAIRD-ID on A_12d (primary);
#     FAIRD-ID on B_22d full and window excluded (complementary).
#     AMMOD-ID Meadow = AMMOD3 vs ID4 on their common days.
#   Between systems, Maize only: the same five, used for disattenuation.
#   Internal references, on the same dates as the pairs they benchmark:
#     FAIRD split-half over both blocks, splits {1,3} vs {2,4} and
#     {1,4} vs {2,3}, summarised as the mean of the two splits;
#     Maize only: FAIRD1 vs FAIRD2, AMMOD1 vs AMMOD2, ID1 vs ID2. AMMOD and ID
#     have a single device with biomass in Meadow, so no other split exists.
#
# Disattenuation (HIGH component, Maize only): reliability of the 2-device
# system series from its 1-vs-1 reference by Spearman-Brown,
# rel = 2 r / (1 + r); r_dis = r_AB / sqrt(rel_A * rel_B); not computed when
# a reference r <= 0. Reading rules, set before running:
#   - values above 1 can arise because co-located systems share site-level
#     signal that the Site1 vs Site2 references do not share;
#   - values below 1 mix detection differences with differences in daily
#     recording coverage (FAIRD 07-19 h continuous; ID sessions within
#     07-09, 12-14 and 17-19 h; AMMOD 00-24 h).
#
# Sensitivity (single, declared): spline with 4 df instead of 3 for the
# between-system PAIRED_BLOCK pairs, both components.
#
# Windows: A_12d = AMMOD biomass window (Aug 23 to Sep 3), full ID series.
# B_22d = Aug 23 to Sep 13, full vs adverse window excluded for both units.
# `times` keeps all dates of the window as levels.
#
# Model cascade (as in the pair models): AR1 -> no AR1 -> no Block; two
# levels without Site/Device intercepts are appended as a last resort,
# because in the Maize-only references each unit is a single device and the
# intercepts are then not identifiable. The reported model must also have a
# finite fixed-effect covariance matrix, which the LOW interval needs.
#
# Part of the shared signal is seasonal and meteorological (one weather
# station, common to all sites); it cannot be partialled out.
#
# A parametric-bootstrap p-value for the HIGH component was implemented in a
# development version and removed here because of its run time (several
# hours); this version reports estimates and intervals only.
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGES: install only what is missing (skipped when everything is present).
# Base R only, so it runs before any library() call.
# ------------------------------------------------------------------------------
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing)
  }
}

# Required
install_if_missing(c("tidyverse", "glmmTMB", "car", "here", "RColorBrewer",
                     "scales", "MASS", "nlme"))

# Optional (extra variance-component output in Scherber_temporal_biomass_5.R);
# a failure here does not stop the script
tryCatch({
  install_if_missing(c("sjPlot", "ape", "remotes"))
  if (!requireNamespace("mixedup", quietly = TRUE)) {
    remotes::install_github("m-clark/mixedup")
  }
}, error = function(e) {
  message("Optional packages not installed: ", conditionMessage(e))
})

library(tidyverse)
library(glmmTMB)
library(splines)
library(ggplot2)

# ------------------------------------------------------------------------------
# PATHS: the only place where folders are defined; nothing below is hardcoded.
# PROJECT_DIR is the project root, i.e. the folder that contains
#   R/00_load_data.R, R/00_custom_colors.R and source_tables/ (input CSVs).
# By default it is located with here::here(), as in 00_load_data.R. To use
# another location, set an environment variable before running, e.g.
#   Sys.setenv(FAIRD_PROJECT_DIR = "C:/path/to/project")
# Outputs are written to OUTPUT_DIR, by default PROJECT_DIR/outputs; it can
# be changed the same way with FAIRD_OUTPUT_DIR.
# ------------------------------------------------------------------------------
PROJECT_DIR <- normalizePath(Sys.getenv("FAIRD_PROJECT_DIR", unset = here::here()),
                             winslash = "/", mustWork = TRUE)
OUTPUT_DIR  <- Sys.getenv("FAIRD_OUTPUT_DIR", unset = file.path(PROJECT_DIR, "outputs"))
R_DIR       <- file.path(PROJECT_DIR, "R")
FIG_DIR     <- file.path(OUTPUT_DIR, "figures")
TAB_DIR     <- file.path(OUTPUT_DIR, "tables")
CON_DIR     <- file.path(OUTPUT_DIR, "console")

for (f in c("00_load_data.R", "00_custom_colors.R")) {
  if (!file.exists(file.path(R_DIR, f))) {
    stop(sprintf("File not found: %s. Set FAIRD_PROJECT_DIR to the project root.",
                 file.path(R_DIR, f)))
  }
}
for (d in c(FIG_DIR, TAB_DIR, CON_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# 00_load_data.R locates its input files with here::i_am(), which resolves
# from the working directory; it is therefore sourced from PROJECT_DIR, and
# the original working directory is restored afterwards.
old_wd <- setwd(PROJECT_DIR)
source(file.path(R_DIR, "00_load_data.R"))
source(file.path(R_DIR, "00_custom_colors.R"))
setwd(old_wd)

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, "Scherber_temporal_concordance_biomass_output.txt"), split = TRUE)

cat("==============================================================================\n")
cat("SYSTEM-LEVEL TEMPORAL CONCORDANCE, DAILY BIOMASS\n")
cat("Components: LOW (shared shape), HIGH (shared day-to-day fluctuation), TOTAL\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
t_script <- Sys.time()


# ==============================================================================
# SECTION 0: PARAMETERS
# ==============================================================================

SPLINE_DF      <- 3     # primary spline df (settled)
SPLINE_DF_SENS <- 4     # single declared sensitivity
N_COEF_SIM     <- 1000  # coefficient draws for the LOW-frequency interval
set.seed(20261002)

cat(sprintf("Spline df: %d (sensitivity %d) | coefficient draws: %d\n\n",
            SPLINE_DF, SPLINE_DF_SENS, N_COEF_SIM))


# ==============================================================================
# SECTION 1: HELPERS
# ==============================================================================

# tmb_is_valid(): copied verbatim from Scherber_temporal_biomass_5.R.
tmb_is_valid <- function(m) {
  if (is.null(m)) return(FALSE)
  pd_hess <- isTRUE(tryCatch(m$sdr$pdHess, error = function(e) FALSE))
  aic_ok  <- tryCatch(is.finite(AIC(m)), error = function(e) FALSE)
  pd_hess && aic_ok
}

# nice_ceiling(): copied verbatim from Scherber_temporal_biomass_5.R.
nice_ceiling <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1 || !is.finite(x) || x <= 0) return(NA_real_)
  step <- 10^floor(log10(x)) / 2
  ceiling(x / step) * step
}

# fit_quiet(): Tweedie log-link fit without console output. Warnings are
# muffled; validity is decided by tmb_is_valid() wherever the fit is used.
fit_quiet <- function(formula, data, start = NULL) {
  withCallingHandlers(
    tryCatch(
      glmmTMB::glmmTMB(formula, family = glmmTMB::tweedie(), data = data,
                       start = start),
      error = function(e) NULL),
    warning = function(w) invokeRestart("muffleWarning"))
}

# tweedie_power(): estimated Tweedie power p of a fitted model.
tweedie_power <- function(m) {
  p <- tryCatch(unname(glmmTMB::family_params(m))[1], error = function(e) NA_real_)
  if (!isTRUE(is.finite(p))) {
    par <- m$fit$par
    nm  <- intersect(c("psi", "thetaf"), names(par))
    if (length(nm) > 0) p <- 1 + stats::plogis(par[[nm[1]]])
  }
  if (!isTRUE(is.finite(p))) stop("Tweedie power could not be extracted from the model.")
  p
}

fisher_ci <- function(r, n) {
  if (!is.finite(r) || n < 4) return(c(NA_real_, NA_real_))
  z <- atanh(r)
  tanh(z + c(-1, 1) * stats::qnorm(0.975) / sqrt(n - 3))
}


# ==============================================================================
# SECTION 2: DATA
# ==============================================================================

cat("=== DATA PREPARATION ===\n\n")

etraps_bio <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, biomass_mg) %>%
  dplyr::rename(Biomass_mg = biomass_mg)

ammod_bio <- data_list$ammod_biomass_for_stats %>%
  dplyr::mutate(Device_type = "AMMOD", Biomass_mg = Live_mass) %>%
  dplyr::select(Date, Device, Device_type, Site, Biomass_mg)

bio_all <- dplyr::bind_rows(etraps_bio, ammod_bio) %>%
  dplyr::mutate(
    Date        = as.Date(Date),
    Date.num    = as.numeric(Date),
    Device      = as.character(Device),
    Device_type = as.character(Device_type),
    Site        = as.character(Site),
    Block       = ifelse(Site %in% c("Site1", "Site2"), "BLOCK1", "BLOCK2")
  ) %>%
  dplyr::filter(!is.na(Biomass_mg))

dates_12 <- seq(min(ammod_bio$Date), max(ammod_bio$Date), by = "day")
dates_22 <- seq(min(etraps_bio$Date), max(etraps_bio$Date), by = "day")

adverse_dates <- data_list$effort_base %>%
  dplyr::filter(Field_Valid, !Field_Valid_Operational_ID) %>%
  dplyr::distinct(Date) %>%
  dplyr::arrange(Date) %>%
  dplyr::pull(Date) %>%
  as.Date()

windows <- list(
  A_12d      = list(dates = dates_12, drop = NULL,          label = "A_12d"),
  B_22d_full = list(dates = dates_22, drop = NULL,          label = "B_22d full"),
  B_22d_excl = list(dates = dates_22, drop = adverse_dates, label = "B_22d window excluded")
)

cat(sprintf("A_12d: %s to %s (%d dates)\n", min(dates_12), max(dates_12), length(dates_12)))
cat(sprintf("B_22d: %s to %s (%d dates)\n", min(dates_22), max(dates_22), length(dates_22)))
cat(sprintf("Adverse window (%d dates): %s\n\n", length(adverse_dates),
            paste(format(adverse_dates), collapse = ", ")))
cat("Observations per device (all valid days):\n")
print(bio_all %>% dplyr::count(Block, Device_type, Device, name = "n"), n = Inf)


# ==============================================================================
# SECTION 3: CORE FUNCTIONS
# ==============================================================================

# build_unit_data(): the two units of one comparison, window and exclusion
# applied, block-day pairing (both units present in the block on that date),
# `times` rebuilt with all dates of the window as levels.
build_unit_data <- function(cmp) {
  w   <- windows[[cmp$window]]
  map <- dplyr::bind_rows(
    tibble::tibble(Device = cmp$dev_A, Unit = cmp$lab_A),
    tibble::tibble(Device = cmp$dev_B, Unit = cmp$lab_B)
  )
  d <- bio_all %>%
    dplyr::inner_join(map, by = "Device") %>%
    dplyr::filter(Date %in% w$dates)
  if (!is.null(w$drop)) d <- d %>% dplyr::filter(!(Date %in% w$drop))
  if (cmp$version == "maize_only") d <- d %>% dplyr::filter(Block == "BLOCK1")
  d %>%
    dplyr::group_by(Block, Date) %>%
    dplyr::filter(dplyr::n_distinct(Unit) == 2) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Unit   = factor(Unit, levels = c(cmp$lab_A, cmp$lab_B)),
      Device = factor(Device),
      Site   = factor(Site),
      Block  = factor(Block),
      times  = factor(as.character(Date), levels = as.character(w$dates))
    ) %>%
    dplyr::arrange(Unit, Device, Date)
}

# make_unit_specs(): model cascade. `Unit` takes the place of Device_type.
make_unit_specs <- function(has_block, df) {
  sp  <- sprintf("Biomass_mg ~ bs(Date.num, %d) * Unit", df)
  sa  <- sprintf("Biomass_mg ~ bs(Date.num, %d) + Unit", df)   # additive (parallel)
  blk <- if (has_block) " + Block" else ""
  re  <- " + (1 | Site / Device)"
  ar  <- " + ar1(times + 0 | Device)"
  mk  <- function(rest, ar1) list(f = paste0(sp, rest), add = paste0(sa, rest), ar1 = ar1)
  lv  <- list()
  lv$AR1 <- mk(paste0(blk, re, ar), TRUE)
  if (has_block) lv$noAR1 <- mk(paste0(blk, re), FALSE)
  lv$simple   <- mk(re,               FALSE)
  lv$AR1_noRE <- mk(paste0(blk, ar),  TRUE)
  lv$fixed    <- mk(blk,              FALSE)
  lv
}

vcov_finite <- function(m) {
  isTRUE(tryCatch(all(is.finite(as.matrix(stats::vcov(m)$cond))),
                  error = function(e) FALSE))
}

first_valid <- function(specs, levels_to_try, dat, need_vcov = FALSE) {
  for (lvl in levels_to_try) {
    m <- fit_quiet(stats::as.formula(specs[[lvl]]$f), dat)
    if (tmb_is_valid(m)) {
      if (!need_vcov || vcov_finite(m)) return(list(model = m, level = lvl))
      cat(sprintf("  [INFO] level '%s' valid but its covariance matrix is not finite -> next level.\n", lvl))
    }
  }
  list(model = NULL, level = NA_character_)
}

# unit_lrt(): LRT of the Date x Unit interaction (constant ratio on the log
# link). Starts at the reported level; moves down the cascade only if the
# additive model is not valid there, refitting BOTH models at each step.
unit_lrt <- function(specs, rep_level, m_rep, dat) {
  out  <- list(p = NA_real_, level = NA_character_, chisq = NA_real_, df = NA_real_)
  if (is.null(m_rep)) return(out)
  lvls <- names(specs)
  for (lvl in lvls[match(rep_level, lvls):length(lvls)]) {
    m_full <- if (lvl == rep_level) m_rep else fit_quiet(stats::as.formula(specs[[lvl]]$f), dat)
    m_add  <- fit_quiet(stats::as.formula(specs[[lvl]]$add), dat)
    if (tmb_is_valid(m_full) && tmb_is_valid(m_add)) {
      ll_f <- stats::logLik(m_full)
      ll_a <- stats::logLik(m_add)
      stat <- 2 * (as.numeric(ll_f) - as.numeric(ll_a))
      ddf  <- attr(ll_f, "df") - attr(ll_a, "df")
      if (lvl != rep_level) {
        cat(sprintf("  [NOTE] LRT computed at level '%s', below the reported level '%s'.\n",
                    lvl, rep_level))
      }
      return(list(p = stats::pchisq(stat, df = ddf, lower.tail = FALSE),
                  level = lvl, chisq = stat, df = ddf))
    }
  }
  cat("  [SKIP] no cascade level with both models valid -> LRT not computed.\n")
  out
}

# low_frequency(): fitted unit curves at the observed dates (reference
# block), r between them, and its interval from coefficient simulation.
# The design matrix of the fit is reused, with the Block column set to the
# reference level; Block is additive and common to both units, so r does not
# depend on the block. The design is checked against predict().
low_frequency <- function(m, dat, n_sim) {
  X    <- as.matrix(glmmTMB::getME(m, "X"))
  beta <- glmmTMB::fixef(m)$cond
  if (nrow(X) != nrow(dat)) stop("Design matrix rows differ from data rows.")
  eta_chk <- stats::predict(m, newdata = dat, type = "link", re.form = NA)
  if (max(abs(drop(X %*% beta) - eta_chk)) > 1e-6) {
    stop("Design matrix does not reproduce predict(); LOW component not computed.")
  }
  grid <- dat %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    dplyr::group_by(Unit, Date) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Unit, Date)
  Xg       <- X[grid$row, , drop = FALSE]
  blk_cols <- grepl("^Block", colnames(Xg))
  Xg[, blk_cols] <- 0
  iA    <- grid$Unit == levels(dat$Unit)[1]
  iB    <- !iA
  dates <- sort(unique(dat$Date))
  if (sum(iA) != length(dates) || sum(iB) != length(dates) ||
      any(grid$Date[iA] != dates) || any(grid$Date[iB] != dates)) {
    stop("Units do not share the same dates; LOW component not computed.")
  }
  fit_hat <- exp(drop(Xg %*% beta))
  r_hat   <- stats::cor(fit_hat[iA], fit_hat[iB])

  V     <- as.matrix(stats::vcov(m)$cond)
  B     <- MASS::mvrnorm(n_sim, mu = beta, Sigma = V)
  E     <- exp(Xg %*% t(B))
  r_sim <- vapply(seq_len(n_sim), function(j) stats::cor(E[iA, j], E[iB, j]), numeric(1))

  # Display only: curves averaged over the blocks present in the data.
  w_blk <- if (any(blk_cols) && dplyr::n_distinct(dat$Block) == 2) {
    mean(c(1, exp(beta[blk_cols])))
  } else 1
  curves <- tibble::tibble(
    Unit  = grid$Unit,
    Date  = grid$Date,
    fit   = fit_hat * w_blk,
    lower = apply(E, 1, stats::quantile, probs = 0.025) * w_blk,
    upper = apply(E, 1, stats::quantile, probs = 0.975) * w_blk
  )
  list(r = r_hat,
       ci = unname(stats::quantile(r_sim, c(0.025, 0.975), na.rm = TRUE)),
       curves = curves)
}

# unit_residuals(): Pearson residuals from the no-AR1 model, averaged within
# block, then across blocks, one series per unit.
unit_residuals <- function(m_res, dat) {
  mu  <- stats::fitted(m_res)
  p   <- tweedie_power(m_res)
  phi <- tryCatch(stats::sigma(m_res), error = function(e) 1)
  dat %>%
    dplyr::mutate(res = (Biomass_mg - mu) / sqrt(phi * mu^p)) %>%
    dplyr::group_by(Unit, Block, Date) %>%
    dplyr::summarise(res = mean(res), .groups = "drop") %>%
    dplyr::group_by(Unit, Date) %>%
    dplyr::summarise(res = mean(res), n_blocks = dplyr::n(), .groups = "drop")
}

series_r <- function(series, value) {
  w <- series %>%
    dplyr::select(Unit, Date, dplyr::all_of(value)) %>%
    tidyr::pivot_wider(names_from = Unit, values_from = dplyr::all_of(value)) %>%
    tidyr::drop_na()
  list(r = stats::cor(w[[2]], w[[3]]), n = nrow(w))
}

# unit_means(): daily unit means, averaged within block, then across blocks.
unit_means <- function(dat) {
  dat %>%
    dplyr::group_by(Unit, Block, Date) %>%
    dplyr::summarise(mean_mg = mean(Biomass_mg), .groups = "drop") %>%
    dplyr::group_by(Unit, Date) %>%
    dplyr::summarise(mean_mg = mean(mean_mg), .groups = "drop")
}

# run_comparison(): one comparison, all components.
run_comparison <- function(cmp, df = SPLINE_DF, n_coef = N_COEF_SIM) {
  dat       <- build_unit_data(cmp)
  has_block <- dplyr::n_distinct(dat$Block) == 2
  specs     <- make_unit_specs(has_block, df)
  n_dates   <- dplyr::n_distinct(dat$Date)

  cat(sprintf("\n>>> %s | %s | %s | spline df %d\n", cmp$label, cmp$window, cmp$version, df))
  cat(sprintf("  %d observations, %d dates, blocks: %s\n", nrow(dat), n_dates,
              paste(levels(dat$Block), collapse = ", ")))
  print(dat %>% dplyr::count(Unit, Block, Device, name = "n"), n = Inf)

  rep_fit <- first_valid(specs, names(specs), dat, need_vcov = TRUE)
  res_fit <- first_valid(specs, names(specs)[!vapply(specs, `[[`, logical(1), "ar1")], dat)
  cat(sprintf("  Reported model level: %s | residual (no-AR1) model level: %s\n",
              rep_fit$level, res_fit$level))

  row <- tibble::tibble(
    id = cmp$id, type = cmp$type, comparison = cmp$label, window = cmp$window,
    version = cmp$version, role = cmp$role, split = cmp$split,
    spline_df = df, n_obs = nrow(dat), n_dates = n_dates,
    model_level = rep_fit$level, residual_level = res_fit$level,
    r_low = NA_real_, low_lo = NA_real_, low_hi = NA_real_,
    r_high = NA_real_, high_lo = NA_real_, high_hi = NA_real_,
    r_total = NA_real_,
    p_LRT = NA_real_, LRT_level = NA_character_, LRT_chisq = NA_real_, LRT_df = NA_real_
  )
  lrt <- unit_lrt(specs, rep_fit$level, rep_fit$model, dat)
  row$p_LRT     <- lrt$p
  row$LRT_level <- lrt$level
  row$LRT_chisq <- lrt$chisq
  row$LRT_df    <- lrt$df
  low <- NULL
  ser <- NULL

  if (!is.null(rep_fit$model)) {
    low <- tryCatch(low_frequency(rep_fit$model, dat, n_coef), error = function(e) {
      cat("  [SKIP] LOW component:", conditionMessage(e), "\n"); NULL
    })
    if (!is.null(low)) {
      row$r_low  <- low$r
      row$low_lo <- low$ci[1]
      row$low_hi <- low$ci[2]
    }
  } else {
    cat("  [ERROR] no valid reported model -> LOW component not computed.\n")
  }

  if (!is.null(res_fit$model)) {
    ser <- unit_residuals(res_fit$model, dat)
    hr  <- series_r(ser, "res")
    row$r_high  <- hr$r
    ci          <- fisher_ci(hr$r, hr$n)
    row$high_lo <- ci[1]
    row$high_hi <- ci[2]
  } else {
    cat("  [ERROR] no valid no-AR1 model -> HIGH component not computed.\n")
  }

  means <- unit_means(dat)
  row$r_total <- series_r(means, "mean_mg")$r

  cat(sprintf("  Parallelism (constant ratio): LRT p = %s (level %s)\n",
              ifelse(is.na(row$p_LRT), "not computed", sprintf("%.4f", row$p_LRT)),
              row$LRT_level))
  cat(sprintf("  LOW r = %.3f [%.3f, %.3f] | HIGH r = %.3f [%.3f, %.3f] | TOTAL r = %.3f\n",
              row$r_low, row$low_lo, row$low_hi, row$r_high, row$high_lo, row$high_hi,
              row$r_total))
  list(cmp = cmp, dat = dat, row = row, low = low, series = ser, means = means)
}


# ==============================================================================
# SECTION 4: COMPARISONS (fixed before running)
# ==============================================================================

sys_dev <- list(AMMOD = paste0("AMMOD", 1:4), FAIRD = paste0("FAIRD", 1:4),
                ID = paste0("ID", 1:4))

cmp_between <- function(a, b, window, version, role) {
  list(id = sprintf("%s-%s|%s|%s", a, b, window, version), type = "between",
       label = sprintf("%s vs %s", a, b), systems = c(a, b),
       lab_A = a, dev_A = sys_dev[[a]], lab_B = b, dev_B = sys_dev[[b]],
       window = window, version = version, role = role, split = NA_character_)
}

cmp_internal <- function(system, lab_A, dev_A, lab_B, dev_B, window, version, split) {
  list(id = sprintf("%s|%s|%s|%s", system, split, window, version), type = "internal",
       label = sprintf("%s: %s vs %s", system, lab_A, lab_B), systems = system,
       lab_A = lab_A, dev_A = dev_A, lab_B = lab_B, dev_B = dev_B,
       window = window, version = version, role = "reference", split = split)
}

comparisons <- list()
pairs_A <- list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))
for (p in pairs_A) {
  comparisons <- c(comparisons, list(
    cmp_between(p[1], p[2], "A_12d", "two_block",  "primary"),
    cmp_between(p[1], p[2], "A_12d", "maize_only", "benchmark")))
}
for (w in c("B_22d_full", "B_22d_excl")) {
  comparisons <- c(comparisons, list(
    cmp_between("FAIRD", "ID", w, "two_block",  "complementary"),
    cmp_between("FAIRD", "ID", w, "maize_only", "benchmark")))
}
for (w in c("A_12d", "B_22d_full", "B_22d_excl")) {
  comparisons <- c(comparisons, list(
    cmp_internal("FAIRD", "FAIRD_h13", c("FAIRD1", "FAIRD3"),
                 "FAIRD_h24", c("FAIRD2", "FAIRD4"), w, "two_block", "S1"),
    cmp_internal("FAIRD", "FAIRD_h14", c("FAIRD1", "FAIRD4"),
                 "FAIRD_h23", c("FAIRD2", "FAIRD3"), w, "two_block", "S2")))
  maize_systems <- if (w == "A_12d") c("FAIRD", "AMMOD", "ID") else c("FAIRD", "ID")
  for (s in maize_systems) {
    comparisons <- c(comparisons, list(
      cmp_internal(s, paste0(s, "1"), paste0(s, "1"), paste0(s, "2"), paste0(s, "2"),
                   w, "maize_only", "maize")))
  }
}

cat(sprintf("\n%d comparisons defined (%d between systems, %d internal references).\n",
            length(comparisons),
            sum(vapply(comparisons, function(x) x$type == "between", logical(1))),
            sum(vapply(comparisons, function(x) x$type == "internal", logical(1)))))


# ==============================================================================
# SECTION 5: RUN (primary spline df)
# ==============================================================================

cat("\n\n=== CONCORDANCE, SPLINE DF", SPLINE_DF, "===\n")

runs <- lapply(comparisons, function(cmp) run_comparison(cmp, df = SPLINE_DF))
names(runs) <- vapply(comparisons, `[[`, character(1), "id")
res_tbl <- dplyr::bind_rows(lapply(runs, `[[`, "row"))


# ==============================================================================
# SECTION 6: SPLIT-HALF SUMMARY, DISATTENUATION
# ==============================================================================

# FAIRD split-half over both blocks: mean of the two splits per window.
split_mean <- res_tbl %>%
  dplyr::filter(type == "internal", version == "two_block") %>%
  dplyr::group_by(window) %>%
  dplyr::summarise(r_low = mean(r_low), r_high = mean(r_high), r_total = mean(r_total),
                   n_splits = dplyr::n(), .groups = "drop")

# Disattenuation, HIGH component, Maize only.
ref_high <- res_tbl %>%
  dplyr::filter(type == "internal", version == "maize_only") %>%
  dplyr::mutate(system = sub(":.*$", "", comparison),
                rel = dplyr::if_else(r_high > 0, 2 * r_high / (1 + r_high), NA_real_)) %>%
  dplyr::select(system, window, r_ref_high = r_high, rel)

dis_tbl <- res_tbl %>%
  dplyr::filter(type == "between", version == "maize_only") %>%
  dplyr::mutate(sys_A = sub(" vs .*$", "", comparison),
                sys_B = sub("^.* vs ", "", comparison)) %>%
  dplyr::left_join(ref_high %>% dplyr::rename(sys_A = system, r_ref_A = r_ref_high, rel_A = rel),
                   by = c("sys_A", "window")) %>%
  dplyr::left_join(ref_high %>% dplyr::rename(sys_B = system, r_ref_B = r_ref_high, rel_B = rel),
                   by = c("sys_B", "window")) %>%
  dplyr::mutate(r_high_dis = dplyr::if_else(!is.na(rel_A) & !is.na(rel_B),
                                            r_high / sqrt(rel_A * rel_B), NA_real_)) %>%
  dplyr::select(comparison, window, r_high, r_ref_A, rel_A, r_ref_B, rel_B, r_high_dis)

res_tbl <- res_tbl %>%
  dplyr::left_join(dis_tbl %>% dplyr::select(comparison, window, rel_A, rel_B, r_high_dis) %>%
                     dplyr::mutate(version = "maize_only"),
                   by = c("comparison", "window", "version"))

fmt <- function(d) d %>% dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 3)))

cat("\n------------------------------------------------------------------------\n")
cat("RESULTS, SPLINE DF", SPLINE_DF, "\n")
cat("------------------------------------------------------------------------\n")
print(fmt(res_tbl %>%
            dplyr::select(comparison, window, version, role, split, n_dates,
                          model_level, residual_level, r_low, low_lo, low_hi,
                          r_high, high_lo, high_hi, r_total)),
      n = Inf, width = Inf)

cat("\n--- PARALLELISM (constant ratio, log link): within systems vs between systems ---\n")
print(fmt(res_tbl %>%
            dplyr::mutate(type = factor(type, levels = c("internal", "between"))) %>%
            dplyr::arrange(window, type, version, comparison) %>%
            dplyr::select(type, comparison, window, version, n_dates, model_level,
                          LRT_level, LRT_chisq, LRT_df, p_LRT)),
      n = Inf, width = Inf)

cat("\n--- FAIRD split-half over both blocks, mean of the two splits ---\n")
print(fmt(split_mean), n = Inf, width = Inf)

cat("\n--- Disattenuated HIGH-frequency r, Maize only (Spearman-Brown) ---\n")
print(fmt(dis_tbl), n = Inf, width = Inf)

cat("\n  NOTES:\n")
cat("  - LOW: r between fitted unit curves, raw scale, interval from coefficient\n")
cat("    simulation; lower limit pessimistic (supports shared shape, not its absence).\n")
cat("  - HIGH: r between unit series of Pearson residuals (no-AR1 model), lag 0;\n")
cat("    CI by Fisher z, indicative only (ignores spline df and autocorrelation).\n")
cat("  - Concordance components are estimates only. The only p-values are the\n")
cat("    LRT of the Date x Unit interaction (constant ratio), unadjusted; for the\n")
cat("    between-system rows they cross-check version 5, where Holm is applied.\n")
cat("  - TOTAL: r between daily unit means, descriptive.\n")
cat("  - maize_only rows: balanced 2-vs-2 devices; disattenuation only there.\n")
cat("  - An interval including 0 is not evidence of absent shared fluctuation.\n")


# ==============================================================================
# SECTION 7: SENSITIVITY, SPLINE DF 4 (between systems, two_block)
# ==============================================================================

cat("\n\n=== SENSITIVITY, SPLINE DF", SPLINE_DF_SENS, "===\n")

sens_cmp  <- Filter(function(x) x$type == "between" && x$version == "two_block", comparisons)
sens_runs <- lapply(sens_cmp, function(cmp) run_comparison(cmp, df = SPLINE_DF_SENS))
sens_tbl  <- dplyr::bind_rows(lapply(sens_runs, `[[`, "row"))

cat("\n--- Primary (df 3) vs sensitivity (df 4) ---\n")
print(fmt(dplyr::bind_rows(res_tbl %>% dplyr::filter(type == "between", version == "two_block"),
                           sens_tbl) %>%
            dplyr::select(comparison, window, spline_df, model_level, r_low, low_lo, low_hi,
                          r_high, high_lo, high_hi, r_total) %>%
            dplyr::arrange(window, comparison, spline_df)),
      n = Inf, width = Inf)


# ==============================================================================
# SECTION 8: FIGURES (linear axes only)
# ==============================================================================

cat("\n\n=== FIGURES ===\n")

fig_theme <- ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 9, color = "grey40"),
    plot.caption     = ggplot2::element_text(size = 8, color = "grey50"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top"
  )

fig_runs <- Filter(function(r) r$cmp$type == "between" && r$cmp$version == "two_block", runs)

for (r in fig_runs) {
  tag <- sprintf("%s_%s", gsub(" vs ", "_", r$cmp$label), r$cmp$window)
  row <- r$row

  # --- Shape figure: fitted unit curves and daily unit means ----------------
  if (!is.null(r$low)) {
    has_ammod <- "AMMOD" %in% r$cmp$systems
    other     <- setdiff(r$cmp$systems, "AMMOD")[1]
    k <- if (has_ammod) {
      signif(mean(r$means$mean_mg[r$means$Unit == "AMMOD"]) /
               mean(r$means$mean_mg[r$means$Unit == other]), 2)
    } else 1
    div    <- function(u, y) ifelse(u == "AMMOD", y / k, y)
    curves <- r$low$curves %>%
      dplyr::mutate(fit = div(Unit, fit), lower = div(Unit, lower), upper = div(Unit, upper))
    pts    <- r$means %>% dplyr::mutate(y = div(Unit, mean_mg))
    y_cap  <- nice_ceiling(max(c(pts$y, curves$fit), na.rm = TRUE) * 1.15)

    p_shape <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = curves,
                           ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = Unit),
                           alpha = 0.15) +
      ggplot2::geom_line(data = curves, ggplot2::aes(x = Date, y = fit, color = Unit),
                         linewidth = 1) +
      ggplot2::geom_point(data = pts, ggplot2::aes(x = Date, y = y, color = Unit),
                          size = 2, alpha = 0.7) +
      ggplot2::scale_color_manual(values = custom_colors) +
      ggplot2::scale_fill_manual(values = custom_colors) +
      ggplot2::scale_y_continuous(
        name     = if (has_ammod) sprintf("%s Biomass (mg)", other) else "Biomass (mg)",
        labels   = scales::label_comma(),
        sec.axis = if (has_ammod) ggplot2::sec_axis(~ . * k,
                     name = sprintf("AMMOD Biomass (mg), axis = left x %s", k),
                     labels = scales::label_comma()) else ggplot2::waiver()) +
      ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
      ggplot2::labs(
        title    = sprintf("Shared temporal shape: %s (%s)", r$cmp$label, windows[[r$cmp$window]]$label),
        subtitle = sprintf("LOW-frequency r = %.2f [%.2f, %.2f] (raw scale, coefficient simulation) | TOTAL r = %.2f | constant-ratio LRT p = %s | model level %s, spline df %d",
                           row$r_low, row$low_lo, row$low_hi, row$r_total,
                           ifelse(is.na(row$p_LRT), "n/a", sprintf("%.4f", row$p_LRT)),
                           row$model_level, SPLINE_DF),
        x = "Date", color = "System", fill = "System",
        caption  = "Curves: fixed effects averaged over blocks, 95% band from coefficient simulation. Points: daily system means (within block, then across blocks).") +
      fig_theme
    f_shape <- sprintf("%s/Scherber_concordance_shape_%s.png", FIG_DIR, tag)
    ggplot2::ggsave(f_shape, p_shape, width = 10, height = 6.5, dpi = 300, bg = "white")
    cat(sprintf("[OK] Saved: %s\n", f_shape))
  }

  # --- Fluctuation figure: unit series of Pearson residuals -----------------
  if (!is.null(r$series)) {
    p_fluct <- ggplot2::ggplot(r$series, ggplot2::aes(x = Date, y = res, color = Unit)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey60") +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = custom_colors) +
      ggplot2::labs(
        title    = sprintf("Shared day-to-day fluctuation: %s (%s)", r$cmp$label, windows[[r$cmp$window]]$label),
        subtitle = sprintf("HIGH-frequency r = %.2f [%.2f, %.2f] (Fisher z, indicative) | residual model level %s",
                           row$r_high, row$high_lo, row$high_hi, row$residual_level),
        x = "Date", y = "Pearson residual (system mean)", color = "System",
        caption  = "Residuals from each system curve (model without AR1), averaged within block, then across blocks.") +
      fig_theme
    f_fluct <- sprintf("%s/Scherber_concordance_fluct_%s.png", FIG_DIR, tag)
    ggplot2::ggsave(f_fluct, p_fluct, width = 10, height = 5, dpi = 300, bg = "white")
    cat(sprintf("[OK] Saved: %s\n", f_fluct))
  }
}

# --- Internal consistency figures: one per window and component -------------
# Same model and display as the between-system figures; each panel is one
# internal comparison, both units in the colour of their system, first unit
# solid line / filled points, second unit dashed line / open points. Ribbons
# are truncated at 1.3 times the panel maximum of curves and points so that
# wide bands do not flatten the curves (display only).
int_runs <- Filter(function(r) r$cmp$type == "internal", runs)
unit_role <- function(u, r) ifelse(as.character(u) == r$cmp$lab_A, "first unit", "second unit")
fmt_p <- function(x, d) ifelse(is.na(x), "n/a", sprintf(paste0("%.", d, "f"), x))

for (w in names(windows)) {
  rw <- Filter(function(r) r$cmp$window == w, int_runs)
  if (length(rw) == 0) next
  panel_shape <- function(r) sprintf("%s [%s] | constant-ratio LRT p = %s | LOW r = %s | TOTAL r = %s",
                                     r$cmp$label, r$cmp$version, fmt_p(r$row$p_LRT, 3),
                                     fmt_p(r$row$r_low, 2), fmt_p(r$row$r_total, 2))
  panel_fluct <- function(r) sprintf("%s [%s] | HIGH r = %s [%s, %s]",
                                     r$cmp$label, r$cmp$version, fmt_p(r$row$r_high, 2),
                                     fmt_p(r$row$high_lo, 2), fmt_p(r$row$high_hi, 2))
  curves_df <- dplyr::bind_rows(lapply(rw, function(r) {
    if (is.null(r$low)) return(NULL)
    r$low$curves %>% dplyr::mutate(System = r$cmp$systems[1], panel = panel_shape(r),
                                   role = unit_role(Unit, r))
  }))
  means_df <- dplyr::bind_rows(lapply(rw, function(r) {
    r$means %>% dplyr::mutate(System = r$cmp$systems[1], panel = panel_shape(r),
                              role = unit_role(Unit, r))
  }))
  series_df <- dplyr::bind_rows(lapply(rw, function(r) {
    if (is.null(r$series)) return(NULL)
    r$series %>% dplyr::mutate(System = r$cmp$systems[1], panel = panel_fluct(r),
                               role = unit_role(Unit, r))
  }))
  if (nrow(curves_df) > 0) {
    caps <- dplyr::bind_rows(
      curves_df %>% dplyr::select(panel, y = fit),
      means_df  %>% dplyr::select(panel, y = mean_mg)) %>%
      dplyr::group_by(panel) %>%
      dplyr::summarise(cap = 1.3 * max(y, na.rm = TRUE), .groups = "drop")
    curves_df <- curves_df %>%
      dplyr::left_join(caps, by = "panel") %>%
      dplyr::mutate(upper = pmin(upper, cap))
  }
  n_pan <- dplyr::n_distinct(means_df$panel)

  p_int_shape <- ggplot2::ggplot() +
    { if (nrow(curves_df) > 0) ggplot2::geom_ribbon(
        data = curves_df,
        ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = System, group = role),
        alpha = 0.12) } +
    { if (nrow(curves_df) > 0) ggplot2::geom_line(
        data = curves_df,
        ggplot2::aes(x = Date, y = fit, color = System, linetype = role), linewidth = 0.9) } +
    ggplot2::geom_point(data = means_df,
                        ggplot2::aes(x = Date, y = mean_mg, color = System, shape = role),
                        size = 1.8, alpha = 0.8) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors, guide = "none") +
    ggplot2::scale_linetype_manual(values = c("first unit" = "solid", "second unit" = "dashed")) +
    ggplot2::scale_shape_manual(values = c("first unit" = 16, "second unit" = 1)) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title    = sprintf("Internal consistency, temporal biomass: %s", windows[[w]]$label),
      subtitle = "Same GLMM as the between-system pairs, with the two units of one system as the factor",
      x = "Date", y = "Biomass (mg)", color = "System", linetype = NULL, shape = NULL,
      caption  = "FAIRD halves: devices {1,3} vs {2,4} and {1,4} vs {2,3}, one device per block each; AMMOD and ID: the two Maize devices. Curves at the observed dates.") +
    fig_theme +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 9, hjust = 0))
  f_is <- file.path(FIG_DIR, sprintf("Scherber_concordance_internal_shape_%s.png", w))
  ggplot2::ggsave(f_is, p_int_shape, width = 10, height = 1.2 + 2.6 * n_pan, dpi = 300, bg = "white")
  cat("[OK] Saved:", f_is, "\n")

  if (nrow(series_df) > 0) {
    p_int_fluct <- ggplot2::ggplot(series_df,
                                   ggplot2::aes(x = Date, y = res, color = System)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey60") +
      ggplot2::geom_line(ggplot2::aes(linetype = role), linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(shape = role), size = 1.8) +
      ggplot2::scale_color_manual(values = custom_colors) +
      ggplot2::scale_linetype_manual(values = c("first unit" = "solid", "second unit" = "dashed")) +
      ggplot2::scale_shape_manual(values = c("first unit" = 16, "second unit" = 1)) +
      ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
      ggplot2::labs(
        title    = sprintf("Internal consistency, day-to-day fluctuation: %s", windows[[w]]$label),
        subtitle = "Pearson residuals from each unit curve (model without AR1), averaged within block, then across blocks",
        x = "Date", y = "Pearson residual (unit mean)", color = "System",
        linetype = NULL, shape = NULL) +
      fig_theme +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 9, hjust = 0))
    f_if <- file.path(FIG_DIR, sprintf("Scherber_concordance_internal_fluct_%s.png", w))
    ggplot2::ggsave(f_if, p_int_fluct, width = 10, height = 1.2 + 2.3 * dplyr::n_distinct(series_df$panel),
                    dpi = 300, bg = "white")
    cat("[OK] Saved:", f_if, "\n")
  }
}

# --- Forest plot: all comparisons, both components --------------------------
forest_main <- dplyr::bind_rows(
  res_tbl %>% dplyr::filter(!(type == "internal" & version == "two_block")) %>%
    dplyr::transmute(comparison, window, version, type,
                     component = "LOW (shared shape)", r = r_low, lo = low_lo, hi = low_hi),
  res_tbl %>% dplyr::filter(!(type == "internal" & version == "two_block")) %>%
    dplyr::transmute(comparison, window, version, type,
                     component = "HIGH (day-to-day)", r = r_high, lo = high_lo, hi = high_hi),
  split_mean %>%
    dplyr::transmute(comparison = "FAIRD: split-half (mean of 2 splits)", window,
                     version = "two_block", type = "internal",
                     component = "LOW (shared shape)", r = r_low, lo = NA_real_, hi = NA_real_),
  split_mean %>%
    dplyr::transmute(comparison = "FAIRD: split-half (mean of 2 splits)", window,
                     version = "two_block", type = "internal",
                     component = "HIGH (day-to-day)", r = r_high, lo = NA_real_, hi = NA_real_),
  dis_tbl %>%
    dplyr::transmute(comparison = paste0(comparison, " (disattenuated)"), window,
                     version = "maize_only", type = "between, disattenuated",
                     component = "HIGH (day-to-day)", r = r_high_dis, lo = NA_real_, hi = NA_real_)
) %>%
  dplyr::filter(!is.na(r)) %>%
  dplyr::mutate(
    window    = factor(windows[as.character(window)] %>% vapply(`[[`, character(1), "label"),
                       levels = vapply(windows, `[[`, character(1), "label")),
    component = factor(component, levels = c("LOW (shared shape)", "HIGH (day-to-day)")),
    label     = sprintf("%s [%s]", comparison, version)
  )

x_max <- max(1.05, max(forest_main$r, na.rm = TRUE) + 0.05)
p_forest <- ggplot2::ggplot(forest_main, ggplot2::aes(x = r, y = label, shape = type)) +
  ggplot2::geom_vline(xintercept = 0, color = "grey60") +
  ggplot2::geom_vline(xintercept = 1, color = "grey80", linetype = "dashed") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), width = 0.25,
                         orientation = "y", na.rm = TRUE) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_shape_manual(values = c("between" = 16, "internal" = 1,
                                         "between, disattenuated" = 17)) +
  ggplot2::facet_grid(window ~ component, scales = "free_y", space = "free_y") +
  ggplot2::coord_cartesian(xlim = c(-1, x_max)) +
  ggplot2::labs(
    title    = "System-level temporal concordance of daily biomass",
    subtitle = sprintf("Pearson r with 95%% intervals (LOW: coefficient simulation; HIGH: Fisher z) | spline df %d", SPLINE_DF),
    x = "r", y = NULL, shape = NULL,
    caption  = "Internal references benchmark the pairs on the same dates. Disattenuated values: Maize only, Spearman-Brown; values above 1 can reflect shared site-level signal.") +
  fig_theme
f_forest <- file.path(FIG_DIR, "Scherber_concordance_forest.png")
ggplot2::ggsave(f_forest, p_forest, width = 12, height = 10, dpi = 300, bg = "white")
cat(sprintf("[OK] Saved: %s\n", f_forest))


# ==============================================================================
# SECTION 9: EXPORT
# ==============================================================================

cat("\n=== EXPORT ===\n\n")
out_tbl <- dplyr::bind_rows(res_tbl, sens_tbl)
readr::write_csv(out_tbl, file.path(TAB_DIR, "Scherber_concordance_biomass.csv"))
readr::write_csv(dis_tbl, file.path(TAB_DIR, "Scherber_concordance_disattenuation.csv"))
cat("[OK] Saved:", file.path(TAB_DIR, "Scherber_concordance_biomass.csv"), "\n")
cat("[OK] Saved:", file.path(TAB_DIR, "Scherber_concordance_disattenuation.csv"), "\n")

cat("\n==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Finished: %s (%.1f min)\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(Sys.time(), t_script, units = "mins"))))
cat("==============================================================================\n")

while (sink.number() > 0) sink()
