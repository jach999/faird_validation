# ==============================================================================
# 00_trajectory_functions.R
# Shared functions for the biomass trajectory models (Tweedie GLMM, log link)
# ------------------------------------------------------------------------------
# Sourced by 3-1a_quantity_trajectories.R and
# supplementary/S1_pair_model_sensitivity.R. Defines no data and runs nothing.
#
# Model family: glmmTMB, Tweedie, log link. Fixed part bs(Date.num, 3) *
# Device_type (+ Block when two blocks are present); random part
# (1 | Site / Device); temporal autocorrelation ar1(times + 0 | Device).
# A fit counts as valid only if its Hessian is positive definite and its AIC
# finite (tmb_is_valid); every cascade moves to the next level otherwise.
# ==============================================================================

# Term of the Date x system interaction in car::Anova() output.
interaction_term <- "bs(Date.num, 3):Device_type"

# Labels of the pairing criteria for two-system comparisons.
pairing_label <- c(Block = "PAIRED_BLOCK", Site = "PAIRED_SITE", none = "FULL")

# tmb_is_valid(): TRUE only if the fit exists, its sdreport Hessian is
# positive definite and its AIC is finite. A glmmTMB fit can be returned
# without an error and still be numerically degenerate.
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

# get_p(): p-value of one term in the Type III Wald table of a valid fit.
get_p <- function(m, term) {
  if (!tmb_is_valid(m)) return(NA_real_)
  av <- tryCatch(car::Anova(m, type = "III"), error = function(e) NULL)
  if (is.null(av) || !(term %in% rownames(av))) return(NA_real_)
  av[term, "Pr(>Chisq)"]
}

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

# build_internal_data(): two units of one system, PAIRED_BLOCK pairing (both
# units present in the block on that date), optional Maize-only subset and
# date exclusion, `times` rebuilt with all levels of the window.
build_internal_data <- function(dat, units, time_levels, maize_only = FALSE,
                                drop_dates = NULL) {
  map <- tibble::tibble(Device = unlist(units, use.names = FALSE),
                        Unit   = rep(names(units), lengths(units)))
  out <- dat %>%
    dplyr::mutate(Device = as.character(Device)) %>%
    dplyr::inner_join(map, by = "Device") %>%
    dplyr::filter(!is.na(Biomass_mg))
  if (!is.null(drop_dates)) out <- out %>% dplyr::filter(!(Date %in% drop_dates))
  if (maize_only) out <- out %>% dplyr::filter(Block == "BLOCK1")
  out <- out %>%
    dplyr::group_by(Block, Date) %>%
    dplyr::filter(dplyr::n_distinct(Unit) == 2) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      System      = as.character(Device_type),
      Device_type = factor(Unit, levels = names(units)),
      Device      = factor(Device)
    )
  out <- droplevels(out)
  out$times <- factor(as.character(out$Date), levels = time_levels)
  out
}

# make_noRE_specs(): last-resort cascade levels for single-device units
# (Maize-only internal comparisons, single Block).
make_noRE_specs <- function() {
  mk <- function(rhs, name) list(
    full = paste0("Biomass_mg ~ bs(Date.num, 3) * Device_type", rhs),
    add  = paste0("Biomass_mg ~ bs(Date.num, 3) + Device_type", rhs),
    name = name
  )
  list(
    AR1_noRE = mk(" + ar1(times + 0 | Device)", "AR1, no Site/Device intercepts, single Block"),
    fixed    = mk("",                           "fixed effects only, single Block")
  )
}

run_internal <- function(system, units, dat, time_levels, analysis, window_label,
                         version, split, maize_only = FALSE, drop_dates = NULL) {
  comparison <- paste(names(units), collapse = " vs ")
  label <- sprintf("INTERNAL %s: %s | %s | %s | %s", system, comparison, analysis,
                   window_label, version)
  idat <- build_internal_data(dat, units, time_levels, maize_only, drop_dates)
  cat(sprintf("\n\n>>> %s\n", label))
  cat(sprintf("  %d observations, %d days (%s to %s)\n", nrow(idat),
              dplyr::n_distinct(idat$Date), min(idat$Date), max(idat$Date)))
  cat("  Observations per unit and device (Device_type holds the unit label):\n")
  print(idat %>% dplyr::count(Block, Device_type, Device, Site, name = "n"))
  has_block <- dplyr::n_distinct(idat$Block) == 2
  if (!has_block) cat("  [NOTE] single Block level -> Block left out of the model.\n")
  specs <- make_pair_specs(has_block)
  if (maize_only) specs <- c(specs, make_noRE_specs())
  res <- fit_pair_cascade(idat, label, specs)
  row <- summarise_pair(res, idat, label) %>%
    dplyr::mutate(System = system, Comparison = comparison, Analysis = analysis,
                  Window = window_label, Version = version, Split = split, .before = 1)
  list(res = res, data = idat, row = row, units = names(units), system = system,
       label = label, analysis = analysis, window = window_label, version = version)
}

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

# make_internal_plot(): internal-consistency figure, two units of one system
# (fixed-effect curves at the reference block, 95% CI, observed unit-days),
# single linear axis because both units belong to one
# system. Both units in the colour of their system; first unit solid line and
# filled points, second unit dashed line and open points.
make_internal_plot <- function(pred, obs, units, system, title, subtitle) {
  pred_s <- pred %>%
    dplyr::mutate(Unit = factor(as.character(Device_type), levels = units), System = system)
  obs_s <- obs %>%
    dplyr::mutate(Unit = factor(as.character(Device_type), levels = units))
  y_cap <- nice_ceiling(stats::quantile(obs_s$Biomass_mg, 0.98, na.rm = TRUE))
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_s,
      ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = System, group = Unit),
      alpha = 0.12) +
    ggplot2::geom_line(
      data = pred_s,
      ggplot2::aes(x = Date, y = fit, color = System, linetype = Unit), linewidth = 1) +
    ggplot2::geom_point(
      data = obs_s,
      ggplot2::aes(x = Date, y = Biomass_mg, color = System, shape = Unit),
      alpha = 0.6, size = 1.8) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_fill_manual(values = custom_colors, guide = "none") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
    ggplot2::scale_shape_manual(values = c(16, 1)) +
    ggplot2::scale_y_continuous(name = "Biomass (mg)", labels = scales::label_comma()) +
    ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
    ggplot2::labs(
      title = title, subtitle = subtitle, x = "Date",
      color = "System", linetype = "Unit", shape = "Unit",
      caption = sprintf("View capped at %s (rounded 98th percentile of plotted points).",
                        scales::label_comma()(y_cap))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 9, color = "grey40"),
      plot.caption     = ggplot2::element_text(size = 8, color = "grey50"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top")
}
