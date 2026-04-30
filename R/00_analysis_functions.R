# ==============================================================================
# ANALYSIS FUNCTIONS - FAIR-Device Validation Project
# ==============================================================================
# Description: Core analysis functions adapted to project data structure
# Author: [Your name]
# Date: 2025-10-10
# ==============================================================================

library(tidyverse)
library(patchwork)
library(vegan)

# ==============================================================================
# HYPOTHESIS 1: INTRA-DEVICE CONSISTENCY
# ==============================================================================

# --- Function 1.1: Device Consistency Analysis ---
analyze_device_consistency <- function(daily_data, device1, device2,
                                       metric = "biomass_mg") {
  #' Analyze consistency between two replicate devices
  #'
  #' @param daily_data Dataframe with daily aggregated data (from data_list$etraps_daily)
  #' @param device1 Name of first device (e.g., "FAIRD1")
  #' @param device2 Name of second device (e.g., "FAIRD2")
  #' @param metric Column name to analyze ("biomass_mg" or "abundance")
  #' @return List with correlation results and plot

  df <- daily_data %>%
    dplyr::filter(Device %in% c(device1, device2)) %>%
    dplyr::select(Date, Device, !!sym(metric)) %>%
    tidyr::pivot_wider(names_from = Device, values_from = !!sym(metric)) %>%
    tidyr::drop_na()

  if (nrow(df) < 3) {
    warning(sprintf("Insufficient data for %s vs %s", device1, device2))
    return(NULL)
  }

  cor_pearson  <- cor.test(df[[device1]], df[[device2]], method = "pearson")
  cor_spearman <- cor.test(df[[device1]], df[[device2]], method = "spearman")
  model        <- lm(df[[device2]] ~ df[[device1]])
  r_squared    <- summary(model)$r.squared

  p <- ggplot2::ggplot(df, ggplot2::aes(x = !!sym(device1), y = !!sym(device2))) +
    ggplot2::geom_point(size = 3, alpha = 0.6, color = "steelblue") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 1) +
    ggplot2::labs(
      title    = sprintf("Consistency Analysis: %s vs %s", device1, device2),
      subtitle = sprintf("Metric: %s | n = %d days | R² = %.3f | r = %.3f",
                         metric, nrow(df), r_squared,
                         as.numeric(cor_pearson$estimate)),
      x = sprintf("%s %s", device1, metric),
      y = sprintf("%s %s", device2, metric)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(face = "bold"),
      panel.grid.minor  = ggplot2::element_blank()
    )

  list(
    devices      = c(device1, device2),
    metric       = metric,
    n_days       = nrow(df),
    pearson_r    = as.numeric(cor_pearson$estimate),
    pearson_p    = cor_pearson$p.value,
    spearman_rho = as.numeric(cor_spearman$estimate),
    spearman_p   = cor_spearman$p.value,
    r_squared    = r_squared,
    model        = model,
    data         = df,
    plot         = p
  )
}

# --- Function 1.2: Coefficient of Variation ---
calculate_CV <- function(daily_data, device1, device2, metric = "biomass_mg") {
  #' Calculate Coefficient of Variation between replicate devices
  #'
  #' @param daily_data Dataframe with daily aggregated data
  #' @param device1 Name of first device
  #' @param device2 Name of second device
  #' @param metric Column name to analyze
  #' @return List with CV statistics and plot

  df <- daily_data %>%
    dplyr::filter(Device %in% c(device1, device2)) %>%
    dplyr::select(Date, Device, !!sym(metric)) %>%
    tidyr::pivot_wider(names_from = Device, values_from = !!sym(metric))

  df_cv <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Mean = mean(c(!!sym(device1), !!sym(device2)), na.rm = TRUE),
      SD   = sd(c(!!sym(device1), !!sym(device2)), na.rm = TRUE),
      CV   = (SD / Mean) * 100
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.finite(CV))

  cv_stats <- df_cv %>%
    dplyr::summarise(
      n_days       = dplyr::n(),
      CV_mean      = mean(CV, na.rm = TRUE),
      CV_median    = median(CV, na.rm = TRUE),
      CV_sd        = sd(CV, na.rm = TRUE),
      CV_min       = min(CV, na.rm = TRUE),
      CV_max       = max(CV, na.rm = TRUE),
      prop_below_20 = sum(CV < 20, na.rm = TRUE) / dplyr::n() * 100
    )

  p <- ggplot2::ggplot(df_cv, ggplot2::aes(x = Date, y = CV)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::geom_hline(yintercept = mean(df_cv$CV, na.rm = TRUE),
                        linetype = "solid", color = "darkgreen", linewidth = 1) +
    ggplot2::annotate("text", x = min(df_cv$Date), y = 22,
                      label = "CV = 20% threshold", hjust = 0, color = "red", size = 3.5) +
    ggplot2::annotate("text", x = min(df_cv$Date),
                      y = mean(df_cv$CV, na.rm = TRUE) + 2,
                      label = sprintf("Mean CV = %.1f%%", mean(df_cv$CV, na.rm = TRUE)),
                      hjust = 0, color = "darkgreen", size = 3.5) +
    ggplot2::labs(
      title    = sprintf("Coefficient of Variation: %s vs %s", device1, device2),
      subtitle = sprintf("Metric: %s", metric),
      x = "Date",
      y = "CV (%)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  list(
    devices  = c(device1, device2),
    metric   = metric,
    cv_data  = df_cv,
    cv_stats = cv_stats,
    plot     = p
  )
}

# --- Function 1.3: Robust device consistency (log10 transformation + diagnostics) ---
analyze_device_consistency_robust <- function(daily_data, device1, device2,
                                              metric        = "biomass_mg",
                                              log_transform = TRUE) {

  df <- daily_data %>%
    dplyr::filter(Device %in% c(device1, device2)) %>%
    dplyr::select(Date, Device, !!sym(metric)) %>%
    tidyr::pivot_wider(names_from = Device, values_from = !!sym(metric)) %>%
    tidyr::drop_na()

  # Apply log10(x + 1) transformation
  if (log_transform) {
    df <- df %>%
      dplyr::mutate(
        dplyr::across(c(!!sym(device1), !!sym(device2)),
                      ~log10(.x + 1),
                      .names = "{.col}_log")
      )
    col1 <- paste0(device1, "_log")
    col2 <- paste0(device2, "_log")
  } else {
    col1 <- device1
    col2 <- device2
  }

  cor_pearson  <- cor.test(df[[col1]],    df[[col2]],    method = "pearson")
  cor_spearman <- cor.test(df[[device1]], df[[device2]], method = "spearman")

  model       <- lm(df[[col2]] ~ df[[col1]])
  shapiro_test <- shapiro.test(residuals(model))

  list(
    transformation       = ifelse(log_transform, "log10(x+1)", "none"),
    pearson_r            = cor_pearson$estimate,
    pearson_p            = cor_pearson$p.value,
    spearman_rho         = cor_spearman$estimate,
    spearman_p           = cor_spearman$p.value,
    residuals_normality_p = shapiro_test$p.value,
    residuals_normal      = shapiro_test$p.value > 0.05
  )
}

# ==============================================================================
# HYPOTHESIS 2: VALIDATION AGAINST AMMOD
# ==============================================================================

# --- Function 2.1: FAIRD-AMMOD Correlation ---
analyze_faird_ammod_correlation <- function(etraps_daily, ammod_biomass,
                                            site_name) {
  #' Correlate FAIRD biomass with AMMOD biomass by site
  #'
  #' @param etraps_daily Daily e-traps data (from data_list$etraps_daily)
  #' @param ammod_biomass AMMOD biomass data (from data_list$ammod_biomass)
  #' @param site_name Site identifier (e.g., "Site1", "Site2")
  #' @return List with correlation results and plots

  faird_df <- etraps_daily %>%
    dplyr::filter(Site == site_name, Device_type == "FAIRD") %>%
    dplyr::group_by(Date, Site) %>%
    dplyr::summarise(faird_biomass_mg = sum(biomass_mg, na.rm = TRUE), .groups = "drop")

  # Live_mass column holds gravimetric biomass in mg
  ammod_df <- ammod_biomass %>%
    dplyr::filter(Site == site_name) %>%
    dplyr::select(Date, Site, Live_mass) %>%
    dplyr::rename(ammod_biomass_mg = Live_mass)

  df <- faird_df %>%
    dplyr::inner_join(ammod_df, by = c("Date", "Site")) %>%
    tidyr::drop_na()

  if (nrow(df) < 3) {
    warning(sprintf("Insufficient data for site %s", site_name))
    return(NULL)
  }

  cor_pearson  <- cor.test(df$ammod_biomass_mg, df$faird_biomass_mg, method = "pearson")
  cor_spearman <- cor.test(df$ammod_biomass_mg, df$faird_biomass_mg, method = "spearman")
  model         <- lm(faird_biomass_mg ~ ammod_biomass_mg, data = df)
  model_summary <- summary(model)
  shapiro_test  <- shapiro.test(residuals(model))

  p_scatter <- ggplot2::ggplot(df, ggplot2::aes(x = ammod_biomass_mg, y = faird_biomass_mg)) +
    ggplot2::geom_point(size = 3, alpha = 0.6, color = "darkgreen") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
    ggplot2::labs(
      title    = sprintf("FAIRD vs AMMOD Biomass - %s", site_name),
      subtitle = sprintf("n = %d days | R² = %.3f | r = %.3f",
                         nrow(df), model_summary$r.squared,
                         as.numeric(cor_pearson$estimate)),
      x = "AMMOD Biomass (mg live weight)",
      y = "FAIRD Biomass (mg live weight)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  p_qq <- ggplot2::ggplot(data.frame(residuals = residuals(model)),
                           ggplot2::aes(sample = residuals)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = "red") +
    ggplot2::labs(
      title = sprintf("Q-Q Plot - %s", site_name),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  p_resid <- ggplot2::ggplot(data.frame(fitted    = fitted(model),
                                        residuals = residuals(model)),
                              ggplot2::aes(x = fitted, y = residuals)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(se = FALSE, color = "blue") +
    ggplot2::labs(
      title = sprintf("Residuals vs Fitted - %s", site_name),
      x = "Fitted values",
      y = "Residuals"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  list(
    site         = site_name,
    n_days       = nrow(df),
    pearson_r    = as.numeric(cor_pearson$estimate),
    pearson_p    = cor_pearson$p.value,
    spearman_rho = as.numeric(cor_spearman$estimate),
    spearman_p   = cor_spearman$p.value,
    r_squared    = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared,
    model         = model,
    shapiro_p     = shapiro_test$p.value,
    data          = df,
    plot_scatter  = p_scatter,
    plot_qq       = p_qq,
    plot_residuals = p_resid
  )
}

# --- Function 2.2: Temporal Synchronization Analysis ---
temporal_synchronization <- function(etraps_daily, ammod_biomass, site_name) {
  #' Analyze temporal synchronization between FAIRD and AMMOD
  #'
  #' @param etraps_daily Daily e-traps data
  #' @param ammod_biomass AMMOD biomass data
  #' @param site_name Site identifier
  #' @return List with synchronization metrics and plots

  faird_df <- etraps_daily %>%
    dplyr::filter(Site == site_name, Device_type == "FAIRD") %>%
    dplyr::group_by(Date, Site) %>%
    dplyr::summarise(faird_biomass_mg = sum(biomass_mg, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(Date)

  ammod_df <- ammod_biomass %>%
    dplyr::filter(Site == site_name) %>%
    dplyr::select(Date, Site, Live_mass) %>%
    dplyr::rename(ammod_biomass_mg = Live_mass) %>%
    dplyr::arrange(Date)

  df <- faird_df %>%
    dplyr::inner_join(ammod_df, by = c("Date", "Site")) %>%
    tidyr::drop_na() %>%
    dplyr::arrange(Date)

  if (nrow(df) < 5) {
    warning(sprintf("Insufficient data for temporal analysis at %s", site_name))
    return(NULL)
  }

  ts_faird <- ts(df$faird_biomass_mg, start = 1, frequency = 1)
  ts_ammod <- ts(df$ammod_biomass_mg, start = 1, frequency = 1)

  ccf_result <- ccf(ts_faird, ts_ammod, lag.max = 5, plot = FALSE)
  max_lag    <- ccf_result$lag[which.max(abs(ccf_result$acf))]
  max_corr   <- max(abs(ccf_result$acf))

  df_trends <- df %>%
    dplyr::mutate(
      faird_change    = c(NA, diff(faird_biomass_mg)),
      ammod_change    = c(NA, diff(ammod_biomass_mg)),
      faird_direction = sign(faird_change),
      ammod_direction = sign(ammod_change),
      concordant      = (faird_direction == ammod_direction) & (faird_direction != 0)
    ) %>%
    dplyr::filter(!is.na(faird_change))

  concordance_pct <- sum(df_trends$concordant, na.rm = TRUE) / nrow(df_trends) * 100

  if (nrow(df_trends) > 2) {
    kendall_test <- cor.test(df_trends$faird_change, df_trends$ammod_change,
                             method = "kendall", use = "complete.obs")
  } else {
    kendall_test <- list(estimate = NA, p.value = NA)
  }

  df_long <- df %>%
    tidyr::pivot_longer(cols = c(faird_biomass_mg, ammod_biomass_mg),
                        names_to = "Device", values_to = "Biomass") %>%
    dplyr::mutate(Device = dplyr::recode(Device,
                                         "faird_biomass_mg" = "FAIRD",
                                         "ammod_biomass_mg" = "AMMOD"))

  p_timeseries <- ggplot2::ggplot(df_long,
                                   ggplot2::aes(x = Date, y = Biomass, color = Device)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = c("FAIRD" = "darkgreen", "AMMOD" = "darkorange")) +
    ggplot2::labs(
      title    = sprintf("Temporal Synchronization - %s", site_name),
      subtitle = sprintf("Concordance: %.1f%% | Kendall τ: %.3f",
                         concordance_pct,
                         ifelse(is.na(kendall_test$estimate), 0,
                                kendall_test$estimate)),
      x = "Date",
      y = "Daily Biomass (mg live weight)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )

  ccf_df <- data.frame(
    Lag         = ccf_result$lag,
    Correlation = ccf_result$acf
  )

  p_ccf <- ggplot2::ggplot(ccf_df, ggplot2::aes(x = Lag, y = Correlation)) +
    ggplot2::geom_hline(yintercept = 0, color = "gray50") +
    ggplot2::geom_segment(ggplot2::aes(xend = Lag, yend = 0),
                          color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "red") +
    ggplot2::labs(
      title    = sprintf("Cross-Correlation Function - %s", site_name),
      subtitle = sprintf("Max correlation: %.3f at lag %d", max_corr, max_lag),
      x = "Lag (days)",
      y = "Correlation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  list(
    site            = site_name,
    n_days          = nrow(df),
    max_ccf         = max_corr,
    optimal_lag     = max_lag,
    concordance_pct = concordance_pct,
    kendall_tau     = ifelse(is.na(kendall_test$estimate), NA,
                             as.numeric(kendall_test$estimate)),
    kendall_p       = kendall_test$p.value,
    data            = df,
    trends_data     = df_trends,
    plot_timeseries = p_timeseries,
    plot_ccf        = p_ccf
  )
}

# ==============================================================================
# HYPOTHESIS 3: COMPARISON WITH INSECT DETECT
# ==============================================================================

# --- Function 3.1: Taxonomic Richness Comparison ---
compare_taxonomic_richness <- function(etraps_data) {
  #' Compare taxonomic richness between FAIRD and ID
  #'
  #' @param etraps_data Combined e-traps individual-level data
  #' @return List with richness comparison and plots

  richness <- etraps_data %>%
    dplyr::group_by(Device_type) %>%
    dplyr::summarise(
      N_individuals = dplyr::n(),
      N_orders      = dplyr::n_distinct(Order,  na.rm = TRUE),
      N_families    = dplyr::n_distinct(Family, na.rm = TRUE),
      N_genera      = dplyr::n_distinct(Genus,  na.rm = TRUE)
    )

  composition <- etraps_data %>%
    dplyr::group_by(Device_type, Order) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Device_type) %>%
    dplyr::mutate(Proportion = Count / sum(Count) * 100)

  contingency_table <- etraps_data %>%
    dplyr::count(Device_type, Order) %>%
    tidyr::pivot_wider(names_from = Device_type, values_from = n, values_fill = 0)

  chi_test <- chisq.test(contingency_table[, -1])

  p_composition <- ggplot2::ggplot(composition,
                                    ggplot2::aes(x = reorder(Order, -Proportion),
                                                 y = Proportion,
                                                 fill = Device_type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::scale_fill_manual(values = c("FAIRD" = "darkgreen", "ID" = "purple")) +
    ggplot2::labs(
      title = "Taxonomic Composition: FAIRD vs Insect Detect",
      x     = "Order",
      y     = "Proportion (%)",
      fill  = "Device"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title  = ggplot2::element_text(face = "bold")
    )

  comm_matrix <- etraps_data %>%
    dplyr::count(Device_type, Order) %>%
    tidyr::pivot_wider(names_from = Order, values_from = n, values_fill = 0) %>%
    tibble::column_to_rownames("Device_type")

  diversity_indices <- data.frame(
    Device   = rownames(comm_matrix),
    Shannon  = vegan::diversity(comm_matrix, index = "shannon"),
    Simpson  = vegan::diversity(comm_matrix, index = "simpson"),
    Richness = vegan::specnumber(comm_matrix)
  )

  list(
    richness          = richness,
    composition       = composition,
    chi_squared       = chi_test,
    diversity_indices = diversity_indices,
    plot              = p_composition
  )
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# --- Create summary table ---
create_summary_table <- function(results_list, type = "correlation") {
  #' Create formatted summary table from analysis results
  #'
  #' @param results_list List of analysis results
  #' @param type Type of table ("correlation", "cv", "temporal")
  #' @return Tibble with formatted results

  if (type == "correlation") {
    purrr::map_dfr(results_list, function(x) {
      if (is.null(x)) return(NULL)
      tibble::tibble(
        Comparison   = paste(x$devices, collapse = " vs "),
        Site         = x$site %||% NA,
        n_days       = x$n_days,
        Pearson_r    = round(x$pearson_r, 3),
        Pearson_p    = format.pval(x$pearson_p, digits = 3),
        Spearman_rho = round(x$spearman_rho, 3),
        Spearman_p   = format.pval(x$spearman_p, digits = 3),
        R_squared    = round(x$r_squared, 3)
      )
    })
  } else if (type == "cv") {
    purrr::map_dfr(results_list, function(x) {
      if (is.null(x)) return(NULL)
      tibble::tibble(
        Comparison    = paste(x$devices, collapse = " vs "),
        CV_mean       = round(x$cv_stats$CV_mean, 2),
        CV_median     = round(x$cv_stats$CV_median, 2),
        CV_SD         = round(x$cv_stats$CV_sd, 2),
        prop_below_20 = round(x$cv_stats$prop_below_20, 1)
      )
    })
  } else if (type == "temporal") {
    purrr::map_dfr(results_list, function(x) {
      if (is.null(x)) return(NULL)
      tibble::tibble(
        Site            = x$site,
        n_days          = x$n_days,
        Max_CCF         = round(x$max_ccf, 3),
        Optimal_lag     = x$optimal_lag,
        Concordance_pct = round(x$concordance_pct, 1),
        Kendall_tau     = round(x$kendall_tau, 3),
        Kendall_p       = format.pval(x$kendall_p, digits = 3)
      )
    })
  }
}

# --- Taxonomic classification fallback ---
# Most specific level available, from Superfamily upward
get_best_classification <- function(row) {
  if (row$Superfamily != "#N/C") {
    return(row$Superfamily)
  } else if (row$Infraorder != "#N/C") {
    return(paste0(row$Infraorder, " (Infraorder)"))
  } else if (row$Suborder != "#N/C") {
    return(paste0(row$Suborder, " (Suborder)"))
  } else if (row$Order != "#N/C") {
    return(paste0(row$Order, " (Order)"))
  } else if (row$Class != "#N/C") {
    return(paste0(row$Class, " (Class)"))
  } else {
    return("Unclassified")
  }
}

# --- Add significance stars to model output ---
add_significance <- function(model) {
  model <- as.data.frame(model)
  model %>%
    dplyr::mutate(significance = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ))
}

# --- Reorder and rename ANOVA output columns ---
reorder_rename_anova <- function(anova_result) {
  anova_df <- as.data.frame(anova_result)
  data.frame(
    term      = rownames(anova_df),
    df        = anova_df$Df,
    sumsq     = anova_df$`Sum Sq`,
    meansq    = anova_df$`Sum Sq` / anova_df$Df,
    statistic = anova_df$`F value`,
    p.value   = anova_df$`Pr(>F)`
  )
}

cat("\n=== ANALYSIS FUNCTIONS LOADED ===\n")
cat("Available functions:\n")
cat("  - analyze_device_consistency()\n")
cat("  - calculate_CV()\n")
cat("  - analyze_device_consistency_robust()\n")
cat("  - analyze_faird_ammod_correlation()\n")
cat("  - temporal_synchronization()\n")
cat("  - compare_taxonomic_richness()\n")
cat("  - create_summary_table()\n")
cat("  - get_best_classification()\n")
cat("  - add_significance()\n")
cat("  - reorder_rename_anova()\n")
