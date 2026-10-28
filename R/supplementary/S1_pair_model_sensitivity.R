# ==============================================================================
# supplementary/S1_pair_model_sensitivity.R
# Supplementary Material S1: sensitivity of the biomass trajectory tests
# ------------------------------------------------------------------------------
# Complements 3-1a_quantity_trajectories.R (same data, same functions in
# R/00_trajectory_functions.R). The primary pair models are refitted here as
# the reference for every sensitivity.
#
# Parts:
#   1. Block as a random effect instead of fixed (three systems): does the
#      interaction conclusion depend on this choice? A two-level random
#      effect is expected to collapse to a singular fit.
#   2. Gaussian cross-check of the three-system model (nlme::lme, random
#      slope over time, corAR1).
#   3. Pairing criterion: PAIRED_SITE (co-located site-days only) and FULL
#      (all observations of the two systems), 12-day window. A variant whose
#      data are identical to PAIRED_BLOCK is not refitted.
#   4. Post hoc: AMMOD vs FAIRD without the first three days (the ratio of
#      the two systems is highest on those days at every site); chosen after
#      seeing the data and reported as such.
#   5. Identity link: the primary and complementary pair runs refitted with
#      an identity link, so the interaction tests a constant DIFFERENCE on the
#      scale of the figures instead of a constant RATIO. AMMOD divided by the
#      data-derived factor k of each pair. Cascade: Tweedie identity (AR1 ->
#      no AR1 -> no Block), then Gaussian identity with dispersion by system.
#      Starting values for the fixed effects from OLS on the same formula
#      (the glmmTMB default of zero gives a zero mean under the identity
#      link, where the Tweedie likelihood is not defined).
# ==============================================================================

required_packages <- c("tidyverse", "here", "glmmTMB", "car", "splines", "nlme", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)
library(glmmTMB)
library(splines)

source(here::here("R", "00_load_data.R"))
source(here::here("R", "00_custom_colors.R"))
source(here::here("R", "00_trajectory_functions.R"))

SCRIPT  <- "S1_pair_model_sensitivity"
FIG_DIR <- here::here("outputs", "figures")
TAB_DIR <- here::here("outputs", "tables")
CON_DIR <- here::here("outputs", "console")
for (d in c(FIG_DIR, TAB_DIR, CON_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(file.path(CON_DIR, paste0(SCRIPT, "_output.txt")), split = TRUE)

cat("==============================================================================\n")
cat("SUPPLEMENTARY S1: SENSITIVITY OF THE BIOMASS TRAJECTORY TESTS\n")
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
# SECTION 2: BLOCK FIXED VS BLOCK RANDOM (three systems)
# ==============================================================================

cat("\n\n=== BLOCK FIXED VS BLOCK RANDOM ===\n")

B1 <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + Block +
    (1 | Site / Device) + ar1(times + 0 | Device),
  bio_overlap, "Model B1 (as in 3-1a): Block fixed + (1|Site/Device) + AR1")
B_ranef <- fit_tmb(
  Biomass_mg ~ bs(Date.num, 3) * Device_type +
    (1 | Block / Site / Device) + ar1(times + 0 | Device),
  bio_overlap, "Model B_ranef: Block RANDOM + (1|Site/Device) + AR1")

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
  Specification  = c("Block FIXED (reported)", "Block RANDOM (check)"),
  p_interaction  = c(p_fixed, p_random),
  Block_variance = c(NA_real_, block_var),
  stringsAsFactors = FALSE)
print(robustness, row.names = FALSE)
cat("\n  A Block variance at or near zero indicates that the two-level random\n")
cat("  effect is a singular fit; the fixed specification is reported.\n")

readr::write_csv(robustness, file.path(TAB_DIR, paste0(SCRIPT, "_block_robustness.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_block_robustness.csv")), "\n")


# ==============================================================================
# SECTION 3: GAUSSIAN CROSS-CHECK (nlme)
# ==============================================================================

cat("\n\n=== GAUSSIAN CROSS-CHECK: nlme::lme, random slope + corAR1 ===\n")
cat("Formula: Biomass_mg ~ bs(Date.num,3) * Device_type + Block,",
    "random = ~Date.num | Site/Device, correlation = corAR1()\n")

BX <- tryCatch(
  nlme::lme(
    Biomass_mg ~ bs(Date.num, 3) * Device_type + Block,
    random      = ~ Date.num | Site / Device,
    correlation = nlme::corAR1(),
    data        = bio_overlap,
    control     = list(opt = "optim")),
  error = function(e) {
    cat("  [FAILED] nlme::lme did not converge:\n    ", conditionMessage(e), "\n")
    NULL
  })

if (!is.null(BX)) {
  cat("  [OK] converged.\n\n")
  print(summary(BX))
  cat("\n--- Type II ANOVA (car::Anova) ---\n")
  print(car::Anova(BX))
}


# ==============================================================================
# SECTION 4: PAIRING CRITERION AND POST HOC CHECK (log link)
# ==============================================================================

cat("\n\n=== PAIR MODELS: PRIMARY (reference) AND SENSITIVITIES ===\n")

pair_types <- list(c("AMMOD", "FAIRD"), c("AMMOD", "ID"), c("FAIRD", "ID"))
pair_runs  <- list()

for (tp in pair_types) {
  primary <- run_pair(tp, bio_overlap, time_levels_12, "Block", "A_12d", "full", "primary")
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

n_drop_initial <- 3
initial_dates  <- sort(unique(af_primary$data$Date))[seq_len(n_drop_initial)]
cat(sprintf("\nPost hoc sensitivity: dates removed = %s\n",
            paste(format(initial_dates), collapse = ", ")))
r <- run_pair(c("AMMOD", "FAIRD"), bio_overlap, time_levels_12, "Block", "A_12d",
              sprintf("drop_first%d", n_drop_initial), "sensitivity_posthoc",
              drop_dates = initial_dates)
pair_runs[[r$label]] <- r

r <- run_pair(c("FAIRD", "ID"), bio_22, time_levels_22, "Block",
              "B_22d", "full", "complementary")
pair_runs[[r$label]] <- r
r <- run_pair(c("FAIRD", "ID"), bio_22, time_levels_22, "Block",
              "B_22d", "window_excluded", "complementary", drop_dates = adverse_dates)
pair_runs[[r$label]] <- r

pair_tbl <- dplyr::bind_rows(lapply(pair_runs, `[[`, "row"))
is_family <- pair_tbl$Analysis == "A_12d" & pair_tbl$Role == "primary"
n_family  <- sum(is_family)
pair_tbl$p_LRT_Holm  <- NA_real_
pair_tbl$p_Wald_Holm <- NA_real_
pair_tbl$p_LRT_Holm[is_family] <- stats::p.adjust(
  pair_tbl$p_interaction_LRT[is_family], method = "holm", n = n_family)
pair_tbl$p_Wald_Holm[is_family] <- stats::p.adjust(
  pair_tbl$p_interaction_Wald[is_family], method = "holm", n = n_family)

cat("\n------------------------------------------------------------------------\n")
cat("DATE x DEVICE_TYPE INTERACTION, LOG LINK: PRIMARY AND SENSITIVITY VARIANTS\n")
cat("------------------------------------------------------------------------\n")
print(pair_tbl %>%
        dplyr::select(Pair, Analysis, Window, Pairing, Role, n_obs, n_days,
                      n_devices, Model, df_interaction, p_interaction_Wald,
                      p_interaction_LRT, LRT_level, p_Wald_at_LRT_level,
                      p_LRT_Holm, p_Wald_Holm) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with("p_"), ~ round(.x, 4))),
      row.names = FALSE)
cat("\n  NOTES:\n")
cat("  - primary and complementary rows reproduce 3-1a (reference).\n")
cat("  - AMMOD-ID PAIRED_SITE is Maize only (Sites 1-2, Block not estimable).\n")
cat("  - drop_first rows are a POST HOC check, not a model choice.\n")

readr::write_csv(pair_tbl, file.path(TAB_DIR, paste0(SCRIPT, "_pairing.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_pairing.csv")), "\n")


# ==============================================================================
# SECTION 5: PAIR MODELS ON THE RAW SCALE (IDENTITY LINK)
# ==============================================================================

cat("\n\n=== PAIR MODELS ON THE RAW SCALE (IDENTITY LINK) ===\n")

identity_families <- list(
  tweedie_identity = list(
    family = glmmTMB::tweedie(link = "identity"),
    disp   = ~ 1,
    name   = "Tweedie identity"
  ),
  gaussian_identity = list(
    family = stats::gaussian(link = "identity"),
    disp   = ~ Device_type,
    name   = "Gaussian identity, dispersion by Device_type"
  )
)

# make_identity_specs(): same cascade levels as make_pair_specs(), on the
# scaled response, with the fixed-only formulas needed for the OLS start.
make_identity_specs <- function(has_block) {
  re  <- " + (1 | Site / Device)"
  ar  <- " + ar1(times + 0 | Device)"
  mk <- function(use_block, use_ar, name) {
    b    <- if (use_block) " + Block" else ""
    rand <- paste0(re, if (use_ar) ar else "")
    list(
      full       = paste0("Biomass_scaled ~ bs(Date.num, 3) * Device_type", b, rand),
      add        = paste0("Biomass_scaled ~ bs(Date.num, 3) + Device_type", b, rand),
      fixed_full = paste0("Biomass_scaled ~ bs(Date.num, 3) * Device_type", b),
      fixed_add  = paste0("Biomass_scaled ~ bs(Date.num, 3) + Device_type", b),
      name       = name
    )
  }
  if (has_block) {
    list(
      AR1    = mk(TRUE,  TRUE,  "Block fixed + (1|Site/Device) + AR1"),
      noAR1  = mk(TRUE,  FALSE, "Block fixed + (1|Site/Device), no AR1"),
      simple = mk(FALSE, FALSE, "(1|Site/Device), no Block, no AR1")
    )
  } else {
    list(
      AR1    = mk(FALSE, TRUE,  "(1|Site/Device) + AR1, single Block"),
      simple = mk(FALSE, FALSE, "(1|Site/Device), single Block, no AR1")
    )
  }
}

# fit_tmb_family(): fit one model with a given family and OLS starting values.
# Warnings are collected and counted rather than printed one by one; validity
# is decided by tmb_is_valid(), as everywhere else in the script.
fit_tmb_family <- function(formula_str, fixed_str, dat, fam, label) {
  cat(sprintf("\n--- %s ---\n", label))
  beta0 <- tryCatch(
    stats::coef(stats::lm(stats::as.formula(fixed_str), data = dat)),
    error = function(e) NULL)
  start <- if (!is.null(beta0) && all(is.finite(beta0))) {
    list(beta = unname(beta0))
  } else NULL
  warn_msgs <- character(0)
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB::glmmTMB(stats::as.formula(formula_str), family = fam$family,
                       dispformula = fam$disp, data = dat, start = start),
      error = function(e) {
        cat("  [FAILED] model errored:\n    ", conditionMessage(e), "\n")
        NULL
      }),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(warn_msgs) > 0) {
    cat(sprintf("  [NOTE] %d warning(s) during fitting, first: %s\n",
                length(warn_msgs), warn_msgs[1]))
  }
  if (tmb_is_valid(m)) {
    cat("  [VALID] positive-definite Hessian, finite AIC.\n")
  } else if (!is.null(m)) {
    cat("  [INVALID] fit returned but numerically degenerate (pdHess/AIC check failed).\n")
  }
  m
}

fit_identity_cascade <- function(dat, label, specs) {
  cat(sprintf("\n##### IDENTITY | %s #####\n", label))
  fam_ids <- names(identity_families)
  for (fam_id in fam_ids) {
    fam <- identity_families[[fam_id]]
    for (lvl in names(specs)) {
      spec <- specs[[lvl]]
      m <- fit_tmb_family(spec$full, spec$fixed_full, dat, fam,
                          sprintf("%s | %s | %s", label, fam$name, spec$name))
      if (tmb_is_valid(m)) {
        if (fam_id != fam_ids[1] || lvl != names(specs)[1]) {
          cat(sprintf("  [WARNING] higher cascade levels invalid -> using '%s', level '%s'.\n",
                      fam$name, lvl))
        }
        return(list(model = m, family_id = fam_id, level = lvl,
                    spec = spec, specs = specs))
      }
    }
  }
  cat("  [ERROR] no valid fit anywhere in the identity-link cascade.\n")
  list(model = NULL, family_id = NA_character_, level = NA_character_,
       spec = NULL, specs = specs)
}

# identity_lrt(): LRT of the interaction (constant difference vs not), full
# and additive models with the same family, random and AR1 structure. Moves
# down the levels of the reported family and then, if needed, to the next
# family, refitting BOTH models at each step (same logic as pair_lrt()).
identity_lrt <- function(res, dat, label) {
  empty <- list(p_lrt = NA_real_, p_wald = NA_real_,
                level = NA_character_, family = NA_character_)
  if (is.null(res$model)) return(empty)
  fam_ids <- names(identity_families)
  lvls    <- names(res$specs)
  for (fi in match(res$family_id, fam_ids):length(fam_ids)) {
    fam_id    <- fam_ids[fi]
    fam       <- identity_families[[fam_id]]
    lvl_start <- if (fam_id == res$family_id) match(res$level, lvls) else 1
    for (lvl in lvls[lvl_start:length(lvls)]) {
      spec   <- res$specs[[lvl]]
      m_full <- if (fam_id == res$family_id && lvl == res$level) res$model else
        fit_tmb_family(spec$full, spec$fixed_full, dat, fam,
                       sprintf("%s | LRT interaction | %s | %s", label, fam$name, spec$name))
      m_add <- fit_tmb_family(spec$add, spec$fixed_add, dat, fam,
                              sprintf("%s | LRT additive | %s | %s", label, fam$name, spec$name))
      if (tmb_is_valid(m_full) && tmb_is_valid(m_add)) {
        ll_full <- stats::logLik(m_full)
        ll_add  <- stats::logLik(m_add)
        lr_stat <- 2 * (as.numeric(ll_full) - as.numeric(ll_add))
        lr_df   <- attr(ll_full, "df") - attr(ll_add, "df")
        p_lrt   <- stats::pchisq(lr_stat, df = lr_df, lower.tail = FALSE)
        cat(sprintf("\n--- LRT (identity link): interaction vs additive, %s, %s, level '%s' ---\n",
                    label, fam$name, lvl))
        if (fam_id != res$family_id || lvl != res$level) {
          cat("  [NOTE] LRT computed at a LOWER cascade step than the reported model,\n")
          cat("         because the additive model was not valid there.\n")
        }
        print(data.frame(
          Model  = c("additive (constant difference)", "interaction"),
          df     = c(attr(ll_add, "df"), attr(ll_full, "df")),
          logLik = round(c(as.numeric(ll_add), as.numeric(ll_full)), 3),
          AIC    = round(c(AIC(m_add), AIC(m_full)), 2),
          Chisq  = c(NA, round(lr_stat, 3)),
          Chi_df = c(NA, lr_df),
          p      = c(NA, p_lrt)
        ), row.names = FALSE)
        p_wald <- tryCatch({
          av <- car::Anova(m_full, type = "III")
          av[interaction_term, "Pr(>Chisq)"]
        }, error = function(e) NA_real_)
        cat(sprintf("  Same step: Wald p = %.4f | LRT p = %.4f\n", p_wald, p_lrt))
        return(list(p_lrt = p_lrt, p_wald = p_wald, level = lvl, family = fam$name))
      }
      cat(sprintf("  [INFO] LRT not possible at '%s', level '%s' (full: %s, additive: %s).\n",
                  fam$name, lvl,
                  ifelse(tmb_is_valid(m_full), "valid", "INVALID"),
                  ifelse(tmb_is_valid(m_add),  "valid", "INVALID")))
    }
  }
  cat("  [SKIP] no cascade step with both models valid -> LRT not computed.\n")
  empty
}

# fitted_at_dates(): fixed-effect fitted curves of both systems at the
# observed dates, reference block, response scale. One column per system.
fitted_at_dates <- function(m, dat) {
  grid <- expand.grid(
    Date.num    = sort(unique(dat$Date.num)),
    Device_type = levels(dat$Device_type)
  ) %>%
    dplyr::mutate(
      Site   = NA,
      Device = NA,
      times  = NA,
      Block  = factor(levels(dat$Block)[1], levels = levels(dat$Block))
    )
  grid$fit <- predict(m, newdata = grid, type = "response", re.form = NA)
  grid %>%
    dplyr::mutate(Device_type = as.character(Device_type)) %>%
    tidyr::pivot_wider(id_cols = Date.num, names_from = Device_type,
                       values_from = fit)
}

# run_identity(): identity-link version of one log-link pair run, on the same
# data object, so that both scales are compared on identical observations.
run_identity <- function(r) {
  has_ammod <- "AMMOD" %in% r$types
  other     <- setdiff(r$types, "AMMOD")[1]
  k_pair    <- if (has_ammod) signif(system_ratio(r$data, "AMMOD", other)$overall, 2) else 1
  dat <- r$data %>%
    dplyr::mutate(Biomass_scaled = dplyr::if_else(Device_type == "AMMOD",
                                                  Biomass_mg / k_pair, Biomass_mg))
  cat(sprintf("\n\n>>> IDENTITY PAIR %s\n", r$label))
  cat(sprintf("  %d observations, %d days | AMMOD scale factor k = %s\n",
              nrow(dat), dplyr::n_distinct(dat$Date), k_pair))
  has_block <- dplyr::n_distinct(dat$Block) == 2
  res <- fit_identity_cascade(dat, r$label, make_identity_specs(has_block))

  row <- data.frame(
    Pair = paste(r$types, collapse = "-"), Analysis = r$analysis,
    Window = r$window, Pairing = r$row$Pairing, Role = r$role,
    k_AMMOD = k_pair, n_obs = nrow(dat), n_days = dplyr::n_distinct(dat$Date),
    Family = NA_character_, Model = NA_character_,
    df_interaction = NA_real_, p_interaction_Wald = NA_real_,
    p_interaction_LRT = NA_real_, LRT_family = NA_character_,
    LRT_level = NA_character_, p_Wald_at_LRT_level = NA_real_,
    diff_min = NA_real_, diff_max = NA_real_,
    ratio_log_min = NA_real_, ratio_log_max = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(res$model)) {
    m <- res$model
    row$Family <- identity_families[[res$family_id]]$name
    row$Model  <- res$spec$name
    av <- tryCatch(car::Anova(m, type = "III"), error = function(e) {
      cat("  [SKIP] Type III ANOVA failed:", conditionMessage(e), "\n"); NULL
    })
    if (!is.null(av)) {
      cat(sprintf("\n--- Type III ANOVA (identity link), %s ---\n", r$label)); print(av)
      if (interaction_term %in% rownames(av)) {
        row$df_interaction     <- av[interaction_term, "Df"]
        row$p_interaction_Wald <- av[interaction_term, "Pr(>Chisq)"]
      }
    }
    lrt <- identity_lrt(res, dat, r$label)
    row$p_interaction_LRT   <- lrt$p_lrt
    row$LRT_family          <- lrt$family
    row$LRT_level           <- lrt$level
    row$p_Wald_at_LRT_level <- lrt$p_wald

    # Fitted difference (identity model), first system minus second, mg on
    # the figure scale (AMMOD divided by k).
    fd <- fitted_at_dates(m, dat)
    d  <- fd[[r$types[1]]] - fd[[r$types[2]]]
    row$diff_min <- min(d)
    row$diff_max <- max(d)
    cat(sprintf("\n--- Fitted curves at observed dates (identity link, %s minus %s) ---\n",
                r$types[1], r$types[2]))
    print(fd %>%
            dplyr::mutate(Date = as.Date(Date.num, origin = "1970-01-01"),
                          difference = d) %>%
            dplyr::select(-Date.num) %>%
            dplyr::relocate(Date) %>%
            dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 1))),
          n = Inf)
  }
  # Fitted ratio (log-link pair model), first system over second, raw mg.
  if (!is.null(r$res$model)) {
    fr <- fitted_at_dates(r$res$model, r$data)
    q  <- fr[[r$types[1]]] / fr[[r$types[2]]]
    row$ratio_log_min <- min(q)
    row$ratio_log_max <- max(q)
  }
  list(res = res, data = dat, row = row, types = r$types, label = r$label,
       analysis = r$analysis, window = r$window, role = r$role, k = k_pair)
}

identity_runs <- lapply(
  Filter(function(r) r$role %in% c("primary", "complementary"), pair_runs),
  run_identity)

# --- Results table, Holm over the same primary family as the log-link models --
id_tbl <- dplyr::bind_rows(lapply(identity_runs, `[[`, "row"))
is_family_id <- id_tbl$Analysis == "A_12d" & id_tbl$Role == "primary"
if (anyNA(id_tbl$p_interaction_LRT[is_family_id])) {
  cat("\n[WARNING] at least one primary pair has no identity LRT p-value; Holm keeps\n")
  cat("          n = 3 comparisons so the correction is not weakened.\n")
}
id_tbl$p_LRT_Holm <- NA_real_
id_tbl$p_LRT_Holm[is_family_id] <- stats::p.adjust(
  id_tbl$p_interaction_LRT[is_family_id], method = "holm", n = sum(is_family_id))

# Side by side with the log-link results (same data objects).
scale_cmp <- id_tbl %>%
  dplyr::left_join(
    pair_tbl %>%
      dplyr::select(Pair, Analysis, Window, Pairing, Role,
                    p_log_LRT = p_interaction_LRT, p_log_LRT_Holm = p_LRT_Holm),
    by = c("Pair", "Analysis", "Window", "Pairing", "Role")) %>%
  dplyr::rename(p_identity_LRT = p_interaction_LRT,
                p_identity_LRT_Holm = p_LRT_Holm)

cat("\n------------------------------------------------------------------------\n")
cat("DATE x DEVICE_TYPE INTERACTION: LOG LINK (constant ratio) vs IDENTITY LINK\n")
cat("(constant difference), same data per row\n")
cat("------------------------------------------------------------------------\n")
print(scale_cmp %>%
        dplyr::select(Pair, Analysis, Window, Role, n_obs, n_days, k_AMMOD,
                      p_log_LRT, p_log_LRT_Holm, ratio_log_min, ratio_log_max,
                      Family, Model, p_interaction_Wald, p_identity_LRT,
                      p_identity_LRT_Holm, LRT_family, LRT_level,
                      diff_min, diff_max) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with("p_"), ~ round(.x, 4)),
                      dplyr::across(dplyr::starts_with("ratio_"), ~ round(.x, 2)),
                      dplyr::across(dplyr::starts_with("diff_"), ~ round(.x, 1))),
      row.names = FALSE)
cat("\n  NOTES:\n")
cat("  - Log link: the interaction tests a constant RATIO between system curves.\n")
cat("  - Identity link: the interaction tests a constant DIFFERENCE, on the scale\n")
cat("    of the pair figures (AMMOD divided by k_AMMOD).\n")
cat("  - ratio_log_*: fitted ratio, first / second system, log model, raw mg.\n")
cat("  - diff_*: fitted difference, first minus second system, identity model,\n")
cat("    mg on the figure scale. Both at the observed dates, reference block.\n")
cat("  - Primary family = the three A_12d PAIRED_BLOCK rows; B_22d rows are\n")
cat("    complementary and outside the Holm family.\n")
cat("  - A non-significant interaction is not evidence of parallelism.\n")

readr::write_csv(scale_cmp, file.path(TAB_DIR, paste0(SCRIPT, "_identity_link.csv")))
cat("[OK] Saved:", file.path(TAB_DIR, paste0(SCRIPT, "_identity_link.csv")), "\n")



# ==============================================================================
# SECTION 6: FIGURES
# ==============================================================================

cat("\n\n=== FIGURES ===\n")

# --- One figure per identity-link pair model ------------------------------------
# Same layout as the log-link pair figures. The model is fitted on AMMOD/k, so
# AMMOD predictions are multiplied back by k here and make_dual_plot() divides
# them again, which puts curves and points on the same axis. Under the identity
# link the CI is symmetric and its lower bound can fall below zero; the view
# starts at zero (coord_cartesian), nothing is dropped from the ribbons.
for (r in identity_runs) {
  if (is.null(r$res$model)) {
    cat(sprintf("\n[SKIP] Identity pair figure %s: no valid model.\n", r$label))
    next
  }
  has_ammod <- "AMMOD" %in% r$types
  other     <- setdiff(r$types, "AMMOD")[1]
  row       <- r$row
  ref_block <- levels(r$data$Block)[1]
  pred_id   <- predict_fixed(r$res$model, r$data, ref_block) %>%
    dplyr::mutate(
      fit   = dplyr::if_else(Device_type == "AMMOD", fit   * r$k, fit),
      lower = dplyr::if_else(Device_type == "AMMOD", lower * r$k, lower),
      upper = dplyr::if_else(Device_type == "AMMOD", upper * r$k, upper)
    )

  window_note <- if (r$window == "window_excluded") {
    sprintf(" | excluded: %s to %s", format(min(adverse_dates)), format(max(adverse_dates)))
  } else ""
  axis_note <- if (has_ammod) sprintf(" | AMMOD right axis, scale factor %s (data-derived)", r$k) else ""

  p_id <- make_dual_plot(
    pred      = pred_id,
    obs       = r$data,
    k         = r$k,
    title     = sprintf("Biomass trajectories: %s (identity link, %s, %s)",
                        sub("-", " vs ", row$Pair), row$Analysis, row$Window),
    subtitle  = sprintf(
      "%s | %s | %s | predictions at %s | constant-difference test: p = %.4f (Wald), %.4f (LRT, %s, level %s)%s%s",
      row$Family, row$Model, row$Pairing, ref_block,
      row$p_interaction_Wald, row$p_interaction_LRT, row$LRT_family,
      row$LRT_level, window_note, axis_note
    ),
    left_name = if (has_ammod) sprintf("%s Biomass (mg)", other) else "Biomass (mg)",
    dual      = has_ammod
  )
  print(p_id)
  fig_file <- sprintf("%s/%s_identity_%s_%s_%s.png", FIG_DIR, SCRIPT,
                      gsub("-", "_", row$Pair), row$Analysis, row$Window)
  ggplot2::ggsave(fig_file, p_id, width = 10, height = 7, dpi = 300, bg = "white")
  cat(sprintf("[OK] Saved: %s\n", fig_file))
}

cat("\n==============================================================================\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

while (sink.number() > 0) sink()
