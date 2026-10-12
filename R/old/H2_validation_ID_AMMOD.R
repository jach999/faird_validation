# ==============================================================================
# H2: ID-AMMOD VALIDATION ANALYSIS (EXPLORATORY)
# ==============================================================================
# Purpose: Exploratory comparison of Insect Detect (ID) vs AMMOD ground truth
#          across the operational overlap window.
#
# IMPORTANT STATISTICAL CAVEAT:
#   The overlap between the AMMOD sampling window (Aug 23 - Sep 3, n=12 days)
#   and the ID operational period (n=17 days, excluding Aug 28 - Sep 1 due to
#   battery failures) yields only n=7 paired days per site:
#     Aug 23-27 (5 days) + Sep 2-3 (2 days)
#   With n=7, statistical power is critically low:
#     - CCF confidence bound: ±1.96/sqrt(7) = ±0.741
#     - Regression significance requires very large effect sizes
#   Results are therefore presented as EXPLORATORY only.
#   No causal or confirmatory interpretations should be drawn.
#
# Sites included:
#   Site1 (Maize Edge): ID1 vs AMMOD1
#   Site2 (Maize Core): ID2 vs AMMOD2
#   Site3: ID3 failed (no data)
#   Site4: AMMOD4 biomass unavailable
#   Meadow cross (AMMOD3+ID4): n=8 AMMOD days only — insufficient, excluded
#
# Methodology: Parallel to H2_validation.R (FAIRD-AMMOD) for direct comparison.
#   All variable names carry _id suffix for cross-script traceability.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(car)
library(broom)
library(ggrepel)
library(lmtest)

# ==============================================================================
# CONSOLE OUTPUT CAPTURE
# ==============================================================================

script_name <- sub("\\.R$", "", basename(sys.frame(1)$ofile))

output_file <- file.path(
  "outputs",
  paste0(script_name, "_output_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

sink(output_file, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Script:", paste0(script_name, ".R"), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("========================================\n")
cat("H2: ID-AMMOD VALIDATION (EXPLORATORY)\n")
cat("========================================\n\n")

cat("*** EXPLORATORY ANALYSIS — n=7 paired days per site ***\n")
cat("*** Statistical results are indicative only.         ***\n\n")

# ==============================================================================
# CUSTOM COLORS
# ==============================================================================

source("custom_colors.R")

cat("[OK] Colors: ID=", custom_colors["ID"],
    ", AMMOD=", custom_colors["AMMOD"], "\n\n")

# ==============================================================================
# SECTION 0: DATA PREPARATION
# ==============================================================================

cat("=== SECTION 0: DATA PREPARATION ===\n\n")

load("data_cleaned.RData")

etraps_daily       <- data_list$etraps_daily
ammod_biomass      <- data_list$ammod_biomass_for_stats
etraps_daily_overlap  <- data_list$etraps_daily_overlap
ammod_biomass_overlap <- data_list$ammod_biomass_overlap

cat("E-traps daily:", nrow(etraps_daily), "obs\n")
cat("AMMOD biomass:", nrow(ammod_biomass), "obs\n\n")

# Operational exclusion window (battery failure period)
ID_EXCLUDE_START <- as.Date("2023-08-28")
ID_EXCLUDE_END   <- as.Date("2023-09-01")

cat(sprintf("ID operational exclusion: %s to %s (%d days removed)\n\n",
            ID_EXCLUDE_START, ID_EXCLUDE_END,
            as.integer(ID_EXCLUDE_END - ID_EXCLUDE_START) + 1))

# Build paired ID-AMMOD dataset
# Sites: Site1 (ID1 vs AMMOD1) and Site2 (ID2 vs AMMOD2)
id_ammod_paired_id <- etraps_daily_overlap %>%
  dplyr::filter(
    Device_type == "ID",
    Site %in% c("Site1", "Site2"),
    # Exclude ID3 (device failure) — already filtered by Site
    # Apply operational exclusion
    !(Date >= ID_EXCLUDE_START & Date <= ID_EXCLUDE_END)
  ) %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    id_biomass_mg  = sum(biomass_mg,  na.rm = TRUE),
    id_abundance   = sum(abundance,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::inner_join(
    ammod_biomass_overlap %>%
      dplyr::filter(Site %in% c("Site1", "Site2")) %>%
      dplyr::select(Date, Site, ammod_live_mass = Live_mass),
    by = c("Date", "Site")
  )

cat("Paired ID-AMMOD observations:\n")
id_ammod_paired_id %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    n_days     = dplyr::n(),
    date_min   = min(Date),
    date_max   = max(Date),
    .groups = "drop"
  ) %>%
  print()
cat("\n")

cat(sprintf("Total paired obs: %d  |  Sites: %s\n\n",
            nrow(id_ammod_paired_id),
            paste(sort(unique(id_ammod_paired_id$Site)), collapse = ", ")))

# ==============================================================================
# ANALYSIS FUNCTIONS
# ==============================================================================

# Regression + diagnostics on log10(x+1) scale
analyze_id_ammod_site_id <- function(paired_data_id, site_name_id) {
  
  site_data_raw_id <- paired_data_id %>%
    dplyr::filter(Site == site_name_id) %>%
    dplyr::select(Date,
                  ID_biomass   = id_biomass_mg,
                  AMMOD_biomass = ammod_live_mass)
  
  if(nrow(site_data_raw_id) < 3) {
    cat(sprintf("   Skipping %s: insufficient data (n=%d)\n",
                site_name_id, nrow(site_data_raw_id)))
    return(NULL)
  }
  
  # Log10(x+1) transformation
  site_data_log_id <- site_data_raw_id %>%
    dplyr::mutate(
      ID_biomass_log    = log10(ID_biomass    + 1),
      AMMOD_biomass_log = log10(AMMOD_biomass + 1)
    )
  
  # OLS regression on log data
  model_id     <- lm(ID_biomass_log ~ AMMOD_biomass_log, data = site_data_log_id)
  shapiro_id   <- shapiro.test(residuals(model_id))
  bp_id        <- car::ncvTest(model_id)
  dw_id        <- lmtest::dwtest(model_id)
  
  model_coef_id  <- coef(summary(model_id))
  slope_id       <- model_coef_id["AMMOD_biomass_log", "Estimate"]
  slope_se_id    <- model_coef_id["AMMOD_biomass_log", "Std. Error"]
  slope_ci_id    <- confint(model_id)["AMMOD_biomass_log", ]
  intercept_id   <- model_coef_id["(Intercept)", "Estimate"]
  
  # Pearson on log data
  pearson_id <- stats::cor.test(
    site_data_log_id$ID_biomass_log,
    site_data_log_id$AMMOD_biomass_log,
    method = "pearson"
  )
  
  list(
    site          = site_name_id,
    n             = nrow(site_data_raw_id),
    data_raw      = site_data_raw_id,
    data_log      = site_data_log_id,
    model         = model_id,
    r_squared     = summary(model_id)$r.squared,
    adj_r_squared = summary(model_id)$adj.r.squared,
    pearson_r     = pearson_id$estimate,
    pearson_p     = pearson_id$p.value,
    shapiro_p     = shapiro_id$p.value,
    bp_p          = bp_id$p,
    dw_statistic  = as.numeric(dw_id$statistic),
    dw_p          = dw_id$p.value,
    slope         = slope_id,
    slope_se      = slope_se_id,
    slope_ci_low  = slope_ci_id[1],
    slope_ci_high = slope_ci_id[2],
    intercept     = intercept_id
  )
}

# Directional concordance + Kendall tau
analyze_temporal_sync_id <- function(site_data_raw_id) {
  if(nrow(site_data_raw_id) < 3) return(NULL)
  
  site_data_id <- site_data_raw_id %>%
    dplyr::arrange(Date) %>%
    dplyr::mutate(
      ID_change    = ID_biomass    - dplyr::lag(ID_biomass),
      AMMOD_change = AMMOD_biomass - dplyr::lag(AMMOD_biomass),
      concordant   = sign(ID_change) == sign(AMMOD_change)
    ) %>%
    dplyr::filter(!is.na(concordant))
  
  if(nrow(site_data_id) == 0) return(NULL)
  
  concordance_pct_id <- 100 * sum(site_data_id$concordant, na.rm = TRUE) /
    nrow(site_data_id)
  
  kendall_id <- stats::cor.test(
    site_data_id$ID_change,
    site_data_id$AMMOD_change,
    method = "kendall",
    exact  = FALSE
  )
  
  list(
    n_changes       = nrow(site_data_id),
    concordance_pct = concordance_pct_id,
    kendall_tau     = kendall_id$estimate,
    kendall_p       = kendall_id$p.value
  )
}

# Significance stars
get_sig_stars_id_h2 <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# ==============================================================================
# SECTION 1: CORRELATIONS (LOG-TRANSFORMED)
# ==============================================================================

cat("\n=== SECTION 1: CORRELATIONS (Log-transformed) ===\n\n")
cat("*** EXPLORATORY: n=7 per site — interpret with caution ***\n\n")

sites_id <- sort(unique(id_ammod_paired_id$Site))  # Site1, Site2
cat(sprintf("Sites: %s\n\n", paste(sites_id, collapse = ", ")))

site_results_id <- list()

for(site_id in sites_id) {
  cat(sprintf("--- %s ---\n", site_id))
  result_id <- analyze_id_ammod_site_id(id_ammod_paired_id, site_id)
  
  if(!is.null(result_id)) {
    site_results_id[[site_id]] <- result_id
    
    cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f\n",
                result_id$n, result_id$r_squared,
                result_id$adj_r_squared, result_id$pearson_p))
    
    cat(sprintf("   Slope: %.3f +/- %.3f SE  [95%% CI: %.3f - %.3f]  (Intercept: %.3f)\n",
                result_id$slope, result_id$slope_se,
                result_id$slope_ci_low, result_id$slope_ci_high,
                result_id$intercept))
    
    cat(sprintf("   Diagnostics: Shapiro p=%.3f | Breusch-Pagan p=%.3f\n",
                result_id$shapiro_p, result_id$bp_p))
    
    dw_flag_id   <- ifelse(result_id$dw_p >= 0.05, "[OK]", "[WARNING]")
    dw_interp_id <- ifelse(result_id$dw_p >= 0.05, "no significant autocorrelation",
                           ifelse(result_id$dw_statistic < 2,
                                  "positive autocorrelation",
                                  "negative autocorrelation"))
    cat(sprintf("   Durbin-Watson: DW=%.3f | p=%.4f -> %s %s\n\n",
                result_id$dw_statistic, result_id$dw_p, dw_interp_id, dw_flag_id))
  }
}

# ==============================================================================
# SECTION 2: TEMPORAL SYNC
# ==============================================================================

cat("\n=== SECTION 2: TEMPORAL SYNC ===\n\n")

for(site_id in sites_id) {
  if(!is.null(site_results_id[[site_id]])) {
    
    sync_result_id <- analyze_temporal_sync_id(site_results_id[[site_id]]$data_raw)
    
    if(!is.null(sync_result_id)) {
      site_results_id[[site_id]]$sync <- sync_result_id
      cat(sprintf("%s: %.1f%% concordance | Kendall tau = %.3f (p = %.4f)\n",
                  site_id,
                  sync_result_id$concordance_pct,
                  as.numeric(sync_result_id$kendall_tau),
                  sync_result_id$kendall_p))
    }
  }
}

# ==============================================================================
# SECTION 3: RESULTS TABLE
# ==============================================================================

cat("\n=== SECTION 3: RESULTS TABLE ===\n\n")
cat("*** All statistics exploratory — n=7 per site ***\n\n")

site_environments_id <- tibble::tibble(
  Site        = c("Site1", "Site2"),
  Environment = c("Maize Edge", "Maize Core")
)

results_table_id <- tibble::tibble(Site = sites_id) %>%
  dplyr::left_join(site_environments_id, by = "Site") %>%
  dplyr::mutate(
    n_days = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]])) site_results_id[[s]]$n else NA_integer_),
    
    R2_log = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]])) site_results_id[[s]]$r_squared else NA_real_),
    
    Pearson_p = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]])) site_results_id[[s]]$pearson_p else NA_real_),
    
    R2_label = paste0(
      formatC(R2_log, digits = 3, format = "f"),
      " ", sapply(Pearson_p, get_sig_stars_id_h2)
    ),
    
    Slope = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]])) site_results_id[[s]]$slope else NA_real_),
    
    Slope_CI_low = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]])) site_results_id[[s]]$slope_ci_low else NA_real_),
    
    Slope_CI_high = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]])) site_results_id[[s]]$slope_ci_high else NA_real_),
    
    Slope_label = sprintf("%.2f [%.2f, %.2f]",
                          Slope, Slope_CI_low, Slope_CI_high),
    
    Concordance_pct = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]]$sync))
        site_results_id[[s]]$sync$concordance_pct else NA_real_),
    
    Kendall_tau = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]]$sync))
        as.numeric(site_results_id[[s]]$sync$kendall_tau) else NA_real_),
    
    Kendall_p = sapply(Site, function(s)
      if(!is.null(site_results_id[[s]]$sync))
        site_results_id[[s]]$sync$kendall_p else NA_real_),
    
    Kendall_label = paste0(
      formatC(Kendall_tau, digits = 3, format = "f"),
      " ", sapply(Kendall_p, get_sig_stars_id_h2)
    ),
    
    Concordance_label = sprintf("%.1f%%", Concordance_pct)
  )

# Display version
table_display_id <- results_table_id %>%
  dplyr::mutate(Devices = "ID vs AMMOD") %>%
  dplyr::select(
    Site,
    Habitat          = Environment,
    Devices,
    n                = n_days,
    `R² (log)`       = R2_label,
    `Slope [95% CI]` = Slope_label,
    `Concordance`    = Concordance_label,
    `Kendall τ`      = Kendall_label
  )

print(table_display_id)

# Export full numeric version
table_export_id <- results_table_id %>%
  dplyr::mutate(Devices = "ID vs AMMOD") %>%
  dplyr::select(
    Site, Environment, Devices,
    n = n_days,
    R2_log, Pearson_p,
    Slope, Slope_CI_low, Slope_CI_high,
    Concordance_pct,
    Kendall_tau, Kendall_p
  ) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

readr::write_csv(table_export_id, "outputs/H2_ID_Validation_Results_Table.csv")
cat("\n[OK] Table saved: outputs/H2_ID_Validation_Results_Table.csv\n\n")

# ==============================================================================
# SECTION 4: CONCORDANCE SCATTER PLOTS (LOG SCALE)
# ==============================================================================

cat("\n=== SECTION 4: CONCORDANCE SCATTER PLOTS (Log Scale) ===\n\n")

plot_list_id <- list()

for(site_id in sites_id) {
  if(!is.null(site_results_id[[site_id]])) {
    result_id  <- site_results_id[[site_id]]
    site_color_id <- custom_colors[site_id]
    
    # Stats annotation
    concordance_pct_id <- if(!is.null(result_id$sync))
      result_id$sync$concordance_pct else NA_real_
    kendall_tau_id <- if(!is.null(result_id$sync))
      as.numeric(result_id$sync$kendall_tau) else NA_real_
    kendall_p_id <- if(!is.null(result_id$sync))
      result_id$sync$kendall_p else NA_real_
    
    if(!is.na(concordance_pct_id) && !is.na(kendall_tau_id)) {
      stats_label_id <- sprintf(
        "R\u00b2 = %.3f\nConcordance: %.1f%%\nKendall \u03c4 = %.3f (p = %.3f)\n[Exploratory, n=%d]",
        result_id$r_squared, concordance_pct_id,
        kendall_tau_id, kendall_p_id, result_id$n
      )
    } else {
      stats_label_id <- sprintf("R\u00b2 = %.3f\n[Exploratory, n=%d]",
                                result_id$r_squared, result_id$n)
    }
    
    # Date labels
    plot_data_id <- result_id$data_log %>%
      dplyr::mutate(Date_label = format(Date, "%b %d"))
    
    env_label_id <- results_table_id$Environment[results_table_id$Site == site_id]
    
    p_id <- ggplot(plot_data_id,
                   aes(x = AMMOD_biomass_log, y = ID_biomass_log)) +
      geom_point(size = 3, alpha = 0.7, color = "gray30") +
      
      ggrepel::geom_text_repel(
        aes(label = Date_label),
        size               = 2.5,
        max.overlaps       = 20,
        box.padding        = 0.35,
        point.padding      = 0.3,
        segment.size       = 0.2,
        segment.alpha      = 0.6,
        min.segment.length = 0,
        force              = 2,
        color              = "gray30"
      ) +
      
      geom_smooth(
        method    = "lm",
        se        = TRUE,
        color     = "gray30",
        fill      = site_color_id,
        alpha     = 0.3,
        linewidth = 1
      ) +
      
      labs(
        title   = sprintf("%s (%s) — Exploratory", site_id, env_label_id),
        caption = "Caution: n=7. CI bound ±0.741. Results not confirmatory.",
        x = "log\u2081\u2080(AMMOD Live Mass + 1)",
        y = "log\u2081\u2080(ID Biomass + 1)"
      ) +
      
      annotate(
        "text",
        x = -Inf, y = Inf,
        label    = stats_label_id,
        hjust    = -0.1, vjust = 1.2,
        size     = 3.5,
        fontface = "italic"
      ) +
      
      theme_minimal(base_size = 11) +
      theme(
        plot.title   = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, face = "italic",
                                    color = "gray50", size = 8)
      )
    
    plot_list_id[[site_id]] <- p_id
    
    # Save individual
    ggsave(sprintf("outputs/H2_ID_Concordance_%s.png", site_id),
           plot = p_id, width = 6, height = 5.5, dpi = 300, bg = "white")
    ggsave(sprintf("outputs/H2_ID_Concordance_%s.pdf", site_id),
           plot = p_id, width = 6, height = 5.5, device = cairo_pdf)
    cat(sprintf("[OK] Saved: outputs/H2_ID_Concordance_%s.png / .pdf\n", site_id))
    
    print(p_id)
  }
}

# Combined
combined_id <- patchwork::wrap_plots(plot_list_id, ncol = length(plot_list_id))
print(combined_id)
ggsave("outputs/H2_ID_Concordance_Combined.png", combined_id,
       width = 6 * length(plot_list_id), height = 5.5, dpi = 300, bg = "white")
ggsave("outputs/H2_ID_Concordance_Combined.pdf", combined_id,
       width = 6 * length(plot_list_id), height = 5.5, device = cairo_pdf)
cat("[OK] Combined concordance plot saved\n\n")

# ==============================================================================
# SECTION 5: TEMPORAL DYNAMICS (PSEUDO-LOG)
# ==============================================================================

cat("\n=== SECTION 5: TEMPORAL DYNAMICS ===\n\n")

temporal_plots_id <- list()

for(site_id in sites_id) {
  if(!is.null(site_results_id[[site_id]])) {
    result_id <- site_results_id[[site_id]]
    sync_id   <- result_id$sync
    
    # Subtitle with stats
    if(!is.null(sync_id)) {
      subtitle_id <- sprintf(
        "R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f) | [Exploratory, n=%d]",
        result_id$r_squared,
        sync_id$concordance_pct,
        as.numeric(sync_id$kendall_tau),
        sync_id$kendall_p,
        result_id$n
      )
    } else {
      subtitle_id <- sprintf("[Exploratory, n=%d]", result_id$n)
    }
    
    env_label_id <- results_table_id$Environment[results_table_id$Site == site_id]
    
    site_temporal_id <- result_id$data_raw %>%
      dplyr::select(Date, ID = ID_biomass, AMMOD = AMMOD_biomass) %>%
      tidyr::pivot_longer(
        cols      = c(ID, AMMOD),
        names_to  = "Method",
        values_to = "Biomass"
      )
    
    cat(sprintf("%s: ID=%s, AMMOD=%s\n",
                site_id, custom_colors["ID"], custom_colors["AMMOD"]))
    
    p_id <- ggplot(site_temporal_id,
                   aes(x = Date, y = Biomass,
                       color = Method, group = Method)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5, alpha = 0.8) +
      # Gap shading: ID battery failure period (Aug 28 - Sep 1)
      annotate(
        "rect",
        xmin  = as.Date("2023-08-27") + 0.5,
        xmax  = as.Date("2023-09-01") + 0.5,
        ymin  = -Inf,
        ymax  = Inf,
        fill  = "gray80",
        alpha = 0.4
      ) +
      annotate(
        "text",
        x     = as.Date("2023-08-30"),
        y     = Inf,
        label = "ID offline\n(battery failure)",
        vjust = 1.3,
        size  = 3,
        color = "gray40",
        fontface = "italic"
      ) +
      scale_color_manual(values = custom_colors, name = "Device Type") +
      
      # Pseudo-log for zeros
      scale_y_continuous(
        trans  = scales::pseudo_log_trans(base = 10, sigma = 1),
        labels = scales::comma,
        breaks = c(0, 10, 100, 1000, 5000)
      ) +
      
      labs(
        title    = sprintf("Temporal Concordance: ID vs AMMOD (%s — %s)",
                           site_id, env_label_id),
        subtitle = subtitle_id,
        x        = "Date",
        y        = "Biomass (mg)"
      ) +
      
      theme_minimal(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1)
      )
    
    temporal_plots_id[[site_id]] <- p_id
    
    # Save individual
    ggsave(sprintf("outputs/H2_ID_Temporal_%s.png", site_id),
           plot = p_id, width = 8, height = 5, dpi = 300, bg = "white")
    ggsave(sprintf("outputs/H2_ID_Temporal_%s.pdf", site_id),
           plot = p_id, width = 8, height = 5, device = cairo_pdf)
    cat(sprintf("[OK] Saved: outputs/H2_ID_Temporal_%s.png / .pdf\n", site_id))
    
    print(p_id)
  }
}

# Combined temporal
combined_temporal_id <- patchwork::wrap_plots(
  temporal_plots_id, ncol = 1, guides = "collect"
) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 10)
  )

print(combined_temporal_id)
ggsave("outputs/H2_ID_Temporal_Combined.png", combined_temporal_id,
       width = 10, height = 4 * length(temporal_plots_id),
       dpi = 300, bg = "white")
ggsave("outputs/H2_ID_Temporal_Combined.pdf", combined_temporal_id,
       width = 10, height = 4 * length(temporal_plots_id),
       device = cairo_pdf)
cat("[OK] Combined temporal plot saved\n\n")

# ==============================================================================
# SECTION 6: CROSS-CORRELATION (EXPLORATORY)
# ==============================================================================

cat("\n=== SECTION 6: CROSS-CORRELATION (Exploratory) ===\n\n")
cat("CI bound: +/-1.96/sqrt(7) = +/-0.741\n")
cat("Note: With n=7, very few lags will be estimable.\n")
cat("      Lags >=3 have <4 pairs and are unreliable.\n\n")

for(site_id in sites_id) {
  if(!is.null(site_results_id[[site_id]])) {
    cat(sprintf("--- %s ---\n", site_id))
    
    ts_raw_id <- site_results_id[[site_id]]$data_raw %>%
      dplyr::arrange(Date)
    
    n_id       <- nrow(ts_raw_id)
    ci_id      <- 1.96 / sqrt(n_id)
    
    cat(sprintf("   n = %d days | CI bound: +/-%.3f\n", n_id, ci_id))
    
    # Log-transform for CCF (stabilizes variance)
    id_ts_id   <- log10(ts_raw_id$ID_biomass    + 1)
    ammod_ts_id <- log10(ts_raw_id$AMMOD_biomass + 1)
    
    # Cross-correlation: ID vs AMMOD
    ccf_result_id <- stats::ccf(
      id_ts_id, ammod_ts_id,
      lag.max = min(3, floor(n_id / 2) - 1),  # Conservative max lag for n=7
      plot    = FALSE,
      na.action = na.pass
    )
    
    # Lag 0 (same-day synchrony)
    lag0_idx_id <- which(ccf_result_id$lag == 0)
    ccf_lag0_id <- ccf_result_id$acf[lag0_idx_id]
    
    cat(sprintf("   >>> LAG 0 (same-day synchrony): CCF = %.3f %s\n",
                ccf_lag0_id,
                ifelse(abs(ccf_lag0_id) > ci_id, "SIGNIFICANT [*]", "(ns)")))
    
    # All lags
    cat("   Detailed CCF (ID -> AMMOD):\n")
    for(i in seq_along(ccf_result_id$lag)) {
      lag_val  <- ccf_result_id$lag[i]
      ccf_val  <- ccf_result_id$acf[i]
      sig_flag <- ifelse(abs(ccf_val) > ci_id, "[*]", "")
      cat(sprintf("     Lag %+d: CCF = %+.3f %s\n", lag_val, ccf_val, sig_flag))
    }
    cat("\n")
  }
}

cat("=== CROSS-CORRELATION INTERPRETATION ===\n\n")
cat("A consistent ID lead of +1 day (as observed in FAIRD) would suggest\n")
cat("the same date-labeling offset in AMMOD data (rotation date vs. capture date).\n")
cat("An absent or inconsistent lead would suggest the pattern is site-specific.\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("========================================\n")
cat("H2 ID VALIDATION — SUMMARY\n")
cat("========================================\n\n")

cat("STATISTICAL CAVEAT:\n")
cat("  n=7 paired days per site (AMMOD window Aug 23-Sep 3 intersected\n")
cat("  with ID operational period, excluding Aug 28-Sep 1 battery failure)\n")
cat("  CCF CI bound: +/-0.741 — only very large effects detectable\n\n")

cat("RESULTS:\n")
for(site_id in sites_id) {
  if(!is.null(site_results_id[[site_id]])) {
    r_id    <- site_results_id[[site_id]]
    env_id  <- results_table_id$Environment[results_table_id$Site == site_id]
    conc_id <- if(!is.null(r_id$sync)) r_id$sync$concordance_pct else NA
    tau_id  <- if(!is.null(r_id$sync)) as.numeric(r_id$sync$kendall_tau) else NA
    tau_p_id <- if(!is.null(r_id$sync)) r_id$sync$kendall_p else NA
    
    cat(sprintf("  %s (%s): R2=%.3f%s | Concordance=%.1f%% | Tau=%.3f%s\n",
                site_id, env_id,
                r_id$r_squared,
                get_sig_stars_id_h2(r_id$pearson_p),
                conc_id,
                tau_id,
                get_sig_stars_id_h2(tau_p_id)))
  }
}

cat("\nFILES GENERATED:\n")
for(site_id in sites_id) {
  cat(sprintf("  [OK] outputs/H2_ID_Concordance_%s.png / .pdf\n", site_id))
  cat(sprintf("  [OK] outputs/H2_ID_Temporal_%s.png / .pdf\n", site_id))
}
cat("  [OK] outputs/H2_ID_Concordance_Combined.png / .pdf\n")
cat("  [OK] outputs/H2_ID_Temporal_Combined.png / .pdf\n")
cat("  [OK] outputs/H2_ID_Validation_Results_Table.csv\n\n")

cat("========================================\n")
cat("H2 ID COMPLETE\n")
cat("========================================\n\n")

# ==============================================================================
# CLOSE LOG
# ==============================================================================
sink()
cat(">> Output saved to:", output_file, "\n")