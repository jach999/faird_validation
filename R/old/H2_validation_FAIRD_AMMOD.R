# ==============================================================================
# H2: FAIRD-AMMOD VALIDATION ANALYSIS
# ==============================================================================
# Updated: 2025-10-28 (JCF/Gemini) - Applied log10(x+1) transformation
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(car)
library(broom)
library(RColorBrewer)
library(dunn.test)
library(ggrepel)  # For date labels without overlap
library(lmtest)   # For Durbin-Watson test

# ============================================================
# CONSOLE OUTPUT CAPTURE
# ============================================================
script_name <- sub("\\.R$", "", basename(sys.frame(1)$ofile))

output_file <- file.path("outputs",
                         paste0(script_name, "_output_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

# Close any previously open sinks to prevent output accumulation
while (sink.number() > 0) sink()
sink(output_file, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Script:", paste0(script_name, ".R"), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("\n========================================\n")
cat("H2: FAIRD-AMMOD VALIDATION ANALYSIS\n")
cat("========================================\n\n")

# ==============================================================================
# CUSTOM COLORS
# ==============================================================================

# Load custom color palette
source("custom_colors.R")

cat("[OK] Colors: FAIRD=", custom_colors["FAIRD"], 
    ", ID=", custom_colors["ID"], 
    ", AMMOD=", custom_colors["AMMOD"], "\n\n")

# Load data

# ==============================================================================
# EFFECT SIZE FUNCTION FOR KRUSKAL-WALLIS
# ==============================================================================

calculate_epsilon_squared <- function(kw_result, n, k) {
  # Calculate epsilon-squared effect size for Kruskal-Wallis test
  # Arguments:
  #   kw_result: result from kruskal.test()
  #   n: total number of observations
  #   k: number of groups
  # Returns:
  #   list with epsilon_squared, interpretation, variance_explained_pct
  # Interpretation thresholds (Cohen, 1988):
  #   Small: 0.01, Medium: 0.06, Large: 0.14
  
  H <- as.numeric(kw_result$statistic)
  epsilon_sq <- (H - k + 1) / (n - k)
  
  interpretation <- dplyr::case_when(
    epsilon_sq < 0.01 ~ "negligible",
    epsilon_sq < 0.06 ~ "small",
    epsilon_sq < 0.14 ~ "medium",
    TRUE ~ "large"
  )
  
  list(
    epsilon_squared = epsilon_sq,
    interpretation = interpretation,
    variance_explained_pct = epsilon_sq * 100
  )
}

load("data_cleaned.RData")

# ==============================================================================
# SECTION 0: DATA PREPARATION
# ==============================================================================

cat("=== LOADING DATA ===\n\n")

etraps_daily <- data_list$etraps_daily
# *** USE PRORATED AMMOD3 DATA ***
ammod_biomass <- data_list$ammod_biomass_for_stats

etraps_daily_overlap <- data_list$etraps_daily_overlap
# *** USE PRORATED AMMOD3 DATA ***
ammod_biomass_overlap <- data_list$ammod_biomass_overlap 

cat("E-traps daily:", nrow(etraps_daily), "obs\n")
cat("AMMOD biomass:", nrow(ammod_biomass), "obs\n")

# Create paired dataset using dplyr:: to avoid conflicts
faird_ammod_paired <- etraps_daily_overlap %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    faird_biomass_mg = sum(biomass_mg, na.rm = TRUE),
    faird_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::inner_join(
    ammod_biomass_overlap %>%
      dplyr::select(Date, Site, ammod_live_mass = Live_mass),
    by = c("Date", "Site")
  )

cat("Paired FAIRD-AMMOD:", nrow(faird_ammod_paired), "obs\n\n")

# ==============================================================================
# SECTION 1-3: ANALYSES (WITH LOG TRANSFORMATION)
# ==============================================================================

analyze_faird_ammod_site <- function(paired_data, site_name) {
  
  # 1. Get raw data
  site_data_raw <- paired_data %>%
    dplyr::filter(Site == site_name) %>%
    dplyr::select(Date, FAIRD_biomass = faird_biomass_mg, AMMOD_biomass = ammod_live_mass)
  
  if(nrow(site_data_raw) < 3) {
    cat(sprintf("   Skipping %s: Insufficient data (n=%d)\n", site_name, nrow(site_data_raw)))
    return(NULL)
  }
  
  # 2. +++ APPLY LOG10(X+1) TRANSFORMATION +++
  site_data_log <- site_data_raw %>%
    dplyr::mutate(
      FAIRD_biomass_log = log10(FAIRD_biomass + 1),
      AMMOD_biomass_log = log10(AMMOD_biomass + 1)
    )
  
  # 3. Run model and diagnostics on LOG DATA
  model        <- lm(FAIRD_biomass_log ~ AMMOD_biomass_log, data = site_data_log)
  shapiro_test <- shapiro.test(residuals(model))
  bp_test      <- car::ncvTest(model)
  dw_test      <- lmtest::dwtest(model)  # Durbin-Watson: tests temporal autocorrelation in residuals
  
  # Extract slope, SE, and 95% CI from OLS regression (log scale)
  model_coef <- coef(summary(model))
  slope      <- model_coef["AMMOD_biomass_log", "Estimate"]
  slope_se   <- model_coef["AMMOD_biomass_log", "Std. Error"]
  slope_ci   <- confint(model)["AMMOD_biomass_log", ]  # 2.5% and 97.5%
  intercept  <- model_coef["(Intercept)", "Estimate"]
  
  # 4. Run Pearson on LOG DATA
  pearson_test <- cor.test(site_data_log$FAIRD_biomass_log,
                           site_data_log$AMMOD_biomass_log,
                           method = "pearson")
  
  # 5. Run Spearman on RAW DATA (non-parametric verification)
  spearman_test <- cor.test(site_data_raw$FAIRD_biomass,
                            site_data_raw$AMMOD_biomass,
                            method = "spearman", exact = FALSE)
  
  list(
    site          = site_name,
    n             = nrow(site_data_raw),
    data_raw      = site_data_raw,
    data_log      = site_data_log,
    model         = model,
    r_squared     = summary(model)$r.squared,
    adj_r_squared = summary(model)$adj.r.squared,
    pearson_r     = pearson_test$estimate,
    pearson_p     = pearson_test$p.value,
    spearman_rho  = spearman_test$estimate,
    spearman_p    = spearman_test$p.value,
    shapiro_p     = shapiro_test$p.value,
    bp_p          = bp_test$p,
    dw_statistic  = as.numeric(dw_test$statistic),  # DW ~ 2 = no autocorrelation
    dw_p          = dw_test$p.value,                # p < 0.05 = significant autocorrelation
    residuals_normal = shapiro_test$p.value > 0.05,
    slope         = slope,
    slope_se      = slope_se,
    slope_ci_low  = slope_ci[1],
    slope_ci_high = slope_ci[2],
    intercept     = intercept
  )
}

analyze_temporal_sync <- function(site_data_raw) {
  # This function uses changes (non-parametric), so raw data is fine
  if(nrow(site_data_raw) < 3) return(NULL)
  
  site_data <- site_data_raw %>%
    dplyr::arrange(Date) %>%
    dplyr::mutate(
      FAIRD_change = FAIRD_biomass - lag(FAIRD_biomass),
      AMMOD_change = AMMOD_biomass - lag(AMMOD_biomass),
      concordant = sign(FAIRD_change) == sign(AMMOD_change)
    ) %>%
    dplyr::filter(!is.na(concordant))
  
  if(nrow(site_data) == 0) return(NULL)
  
  concordance_pct <- 100 * sum(site_data$concordant, na.rm = TRUE) / nrow(site_data)
  kendall_test <- cor.test(site_data$FAIRD_change, site_data$AMMOD_change, method = "kendall", exact = FALSE)
  
  list(n_changes = nrow(site_data), concordance_pct = concordance_pct,
       kendall_tau = kendall_test$estimate, kendall_p = kendall_test$p.value)
}

cat("\n=== SECTION 1: CORRELATIONS (Log-transformed) ===\n\n")

sites <- unique(faird_ammod_paired$Site)
cat(sprintf("Found %d valid sites for FAIRD-AMMOD correlation: %s\n\n",
            length(sites), paste(sort(sites), collapse = ", ")))

site_results <- list()

for(site in sites) {
  cat(sprintf("--- %s ---\n", site))
  result <- analyze_faird_ammod_site(faird_ammod_paired, site)
  
  if(!is.null(result)) {
    site_results[[site]] <- result
    
    cat(sprintf("n=%d | R2=%.3f (Adj.R2=%.3f) | p=%.4f\n",
                result$n, result$r_squared, result$adj_r_squared, result$pearson_p))
    
    cat(sprintf("   Slope: %.3f +/- %.3f SE  [95%% CI: %.3f - %.3f]  (Intercept: %.3f)\n",
                result$slope, result$slope_se,
                result$slope_ci_low, result$slope_ci_high,
                result$intercept))
    
    cat(sprintf("   Diagnostics: Shapiro p=%.3f | Breusch-Pagan p=%.3f\n",
                result$shapiro_p, result$bp_p))
    
    dw_flag   <- ifelse(result$dw_p >= 0.05, "[OK]", "[WARNING]")
    dw_interp <- ifelse(result$dw_p >= 0.05, "no significant autocorrelation",
                        ifelse(result$dw_statistic < 2, "positive autocorrelation",
                               "negative autocorrelation"))
    cat(sprintf("   Durbin-Watson: DW=%.3f | p=%.4f -> %s %s\n\n",
                result$dw_statistic, result$dw_p, dw_interp, dw_flag))
  }
}

cat("\n=== SECTION 2: TEMPORAL SYNC ===\n\n")

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    # Pass raw data (data_raw) to synchronization function
    sync_result <- analyze_temporal_sync(site_results[[site]]$data_raw)
    if(!is.null(sync_result)) {
      site_results[[site]]$sync <- sync_result
      cat(sprintf("%s: %.1f%% concordance\n", site, sync_result$concordance_pct))
    }
  }
}

cat("\n=== SECTION 3.1: VALIDATION RESULTS TABLE ===\n\n")

site_environments <- tibble(
  Site = c("Site1", "Site2", "Site3", "Site4"),
  Environment = c("Maize", "Maize", "Meadow", "Meadow") 
)

results_table <- tibble(
  Site = sites
) %>%
  dplyr::left_join(site_environments, by = "Site") %>%
  dplyr::mutate(
    n_days = sapply(Site, function(s) if(!is.null(site_results[[s]])) site_results[[s]]$n else NA),
    R2_log = sapply(Site, function(s) if(!is.null(site_results[[s]])) site_results[[s]]$r_squared else NA),
    Pearson_p_log = sapply(Site, function(s) if(!is.null(site_results[[s]])) site_results[[s]]$pearson_p else NA),
    Spearman_rho = sapply(Site, function(s) if(!is.null(site_results[[s]])) site_results[[s]]$spearman_rho else NA),
    Spearman_p = sapply(Site, function(s) if(!is.null(site_results[[s]])) site_results[[s]]$spearman_p else NA)
  )

print(results_table)
write_csv(results_table, "outputs/H2_FAIRD_AMMOD_Validation_Log.csv")
cat("\n[OK] Table saved\n\n")


cat("\n=== SECTION 3.2: SUMMARY TABLE ===\n\n")

site_environments <- tibble::tibble(
  Site        = c("Site1", "Site2", "Site3", "Site4"),
  Environment = c("Maize", "Maize", "Meadow", "Meadow")
)

get_sig_stars_h2 <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Cross-correlation values from H2_validation_cross_corr.R output
# CCF lag-0 confidence bounds: ±0.566 (n=12), ±0.693 (n=8)
ccf_data <- tibble::tibble(
  Site        = c("Site1", "Site2", "Site3"),
  CCF_lag0    = c(0.382,  0.696, -0.010),
  CCF_lag0_p  = c(0.566,  0.566,  0.693),   # CI bound (not a p-value)
  FAIRD_lead  = c("+1 day", "+1 day", "+2 days")
) %>%
  dplyr::mutate(
    # Significant if |CCF| exceeds CI bound (±1.96/sqrt(n))
    CCF_sig = dplyr::if_else(abs(CCF_lag0) > CCF_lag0_p, "*", "ns"),
    CCF_label = sprintf("%.3f %s", CCF_lag0, CCF_sig)
  )

# Build main table
results_table <- tibble::tibble(Site = sites) %>%
  dplyr::left_join(site_environments, by = "Site") %>%
  dplyr::mutate(
    n_days = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$n else NA_integer_),
    
    R2_log = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$r_squared else NA_real_),
    
    Pearson_p = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$pearson_p else NA_real_),
    
    R2_label = paste0(
      formatC(R2_log, digits = 3, format = "f"),
      " ", sapply(Pearson_p, get_sig_stars_h2)
    ),
    
    Slope = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$slope else NA_real_),
    
    Slope_CI_low = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$slope_ci_low else NA_real_),
    
    Slope_CI_high = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$slope_ci_high else NA_real_),
    
    Slope_label = sprintf("%.2f [%.2f, %.2f]",
                          Slope, Slope_CI_low, Slope_CI_high),
    
    Concordance_pct = sapply(Site, function(s)
      if(!is.null(site_results[[s]]$sync))
        site_results[[s]]$sync$concordance_pct else NA_real_),
    
    Kendall_tau = sapply(Site, function(s)
      if(!is.null(site_results[[s]]$sync))
        as.numeric(site_results[[s]]$sync$kendall_tau) else NA_real_),
    
    Kendall_p = sapply(Site, function(s)
      if(!is.null(site_results[[s]]$sync))
        site_results[[s]]$sync$kendall_p else NA_real_),
    
    Kendall_label = paste0(
      formatC(Kendall_tau, digits = 3, format = "f"),
      " ", sapply(Kendall_p, get_sig_stars_h2)
      ),
    
    Spearman_rho = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$spearman_rho else NA_real_),

    Spearman_p = sapply(Site, function(s)
      if(!is.null(site_results[[s]])) site_results[[s]]$spearman_p else NA_real_),
    
    Spearman_label = paste0(
      formatC(Spearman_rho, digits = 3, format = "f"),
      " ", sapply(Spearman_p, get_sig_stars_h2)
     ),
    
    Concordance_label = sprintf("%.1f%%", Concordance_pct)
  ) %>%
  dplyr::left_join(
    ccf_data %>% dplyr::select(Site, CCF_label, FAIRD_lead),
    by = "Site"
  )

# Display version
table_display <- results_table %>%
  dplyr::select(
    Site,
    Habitat        = Environment,
    n              = n_days,
    `R² (log)`     = R2_label,
    `Slope [95% CI]`  = Slope_label,
    `Concordance`  = Concordance_label,
    `Spearman ρ` = Spearman_label,
    `Kendall τ`    = Kendall_label,
    `CCF lag 0`    = CCF_label,
    `FAIRD lead`   = FAIRD_lead
  )

print(table_display)

# Export — full numeric version
table_export <- results_table %>%
  dplyr::select(
    Site, Environment, n = n_days,
    R2_log, Pearson_p,
    Slope, Slope_CI_low, Slope_CI_high,
    Concordance_pct,
    Spearman_rho, Spearman_p,
    Kendall_tau, Kendall_p
  ) %>%
  dplyr::left_join(
    ccf_data %>% dplyr::select(Site, CCF_lag0, FAIRD_lead),
    by = "Site"
  ) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

# Publication version (with significance stars)
readr::write_csv(table_display,
                 "outputs/H2_FAIRD_AMMOD_Publication_Table.csv")
cat("[OK] Table saved: outputs/H2_FAIRD_AMMOD_Publication_Table.csv\n\n")

# Numeric backup version (for PCV verification and supplementary)
readr::write_csv(table_export,
                 "outputs/H2_FAIRD_AMMOD_Validation_Results_Table.csv")
cat("[OK] Table saved: outputs/H2_FAIRD_AMMOD_Validation_Results_Table.csv\n\n")

# ==============================================================================
# SECTION 4: CORRELATION PLOTS (ON LOG SCALE)
# ==============================================================================

cat("\n=== SECTION 4: CORRELATION PLOTS (Log Scale) ===\n\n")

plot_list <- list()

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    result <- site_results[[site]]
    site_color <- custom_colors[site]
    
    # +++ PREPARE DATE LABELS +++
    # Format dates as "Aug 23" for readability
    plot_data <- result$data_log %>%
      dplyr::mutate(
        Date_label = format(Date, "%b %d"),
        # AMMOD3 biomass data unavailable due probes processing change (ATL) dates (Aug 31 - Sep 3)
        Is_problematic = (site == "Site3" & 
                            Date >= as.Date("2023-08-31") & 
                            Date <= as.Date("2023-09-03"))
      )
    
    # +++ SCATTER PLOT WITH DATE LABELS +++
    p <- ggplot(plot_data, aes(x = AMMOD_biomass_log, y = FAIRD_biomass_log)) +
      # Points - color problematic dates in red for Site3
      geom_point(aes(color = Is_problematic), size = 3, alpha = 0.7) + 
      scale_color_manual(values = c("FALSE" = "gray30", "TRUE" = "red"), 
                         guide = "none") +
      
      # Add date labels with ggrepel to avoid overlap
      ggrepel::geom_text_repel(
        aes(label = Date_label, color = Is_problematic),
        size = 2.5,
        max.overlaps = 20,
        box.padding = 0.35,
        point.padding = 0.3,
        segment.size = 0.2,
        segment.alpha = 0.6,
        min.segment.length = 0,
        force = 2,
        show.legend = FALSE
      ) +
      
      # Regression line
      geom_smooth(
        method = "lm", 
        se = TRUE, 
        color = "gray30",       
        fill = site_color,      
        alpha = 0.3, 
        linewidth = 1
      ) +
      
      # Labels and annotations
      labs(
        title = sprintf("%s (%s)", site, 
                        results_table$Environment[results_table$Site == site]),
        x = "log10(AMMOD Live Mass + 1)", 
        y = "log10(FAIRD Biomass + 1)",
      ) +
      
      # R² and p-value annotation
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("R² = %.3f\np = %.3f", 
                               result$r_squared, result$pearson_p),
               hjust = -0.1, vjust = 1.2, size = 3.5, fontface = "italic") +
      
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, face = "italic", 
                                    color = "red", size = 9)
      )
    
    plot_list[[site]] <- p
    
    # +++ NEW: SAVE INDIVIDUAL PLOTS +++
    ind_filename <- sprintf("outputs/H2_Correlation_%s", site)
    ggsave(paste0(ind_filename, ".png"), plot = p, width = 6, height = 6, dpi = 300, bg = "white")
    ggsave(paste0(ind_filename, ".pdf"), plot = p, width = 6, height = 6, device = cairo_pdf)
    cat(sprintf("   [OK] Saved individual plot: %s.png\n", ind_filename))
  }
}

# --- Save Combined Plot ---
combined <- wrap_plots(plot_list, ncol = length(plot_list))
print(combined)
ggsave("outputs/H2_Validation_Log_with_dates.png", combined, width = 5 * length(plot_list), height = 5, dpi = 300, bg = "white")
ggsave("outputs/H2_Validation_Log_with_dates.pdf", combined, width = 5 * length(plot_list), height = 5, device = cairo_pdf)

cat("\n[OK] Combined correlation plot saved\n\n")

# ==============================================================================
# SECTION 5: TEMPORAL DYNAMICS (WITH PSEUDO-LOG TRANSFORMATION)
# ==============================================================================

cat("\n=== SECTION 5: TEMPORAL DYNAMICS ===\n\n")

temporal_plots <- list()

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    result <- site_results[[site]]
    sync <- result$sync
    
    # Prepare data (use data_raw)
    site_temporal <- result$data_raw %>%
      dplyr::select(Date, FAIRD = FAIRD_biomass, AMMOD = AMMOD_biomass) %>%
      tidyr::pivot_longer(cols = c(FAIRD, AMMOD), names_to = "Method", values_to = "Biomass")
    
    faird_color <- custom_colors["FAIRD"]
    ammod_color <- custom_colors["AMMOD"]
    
    cat(sprintf("%s: FAIRD=%s, AMMOD=%s\n", site, faird_color, ammod_color))
    
    # +++ NEW: PREPARE SUBTITLE WITH STATS (R2, Concordance, Tau, p) +++
    subtitle_text <- ""
    if(!is.null(sync)) {
      subtitle_text <- sprintf("R\u00b2 = %.3f | Concordance: %.1f%% | Kendall \u03c4 = %.3f (p = %.3f)",
                               result$r_squared,
                               sync$concordance_pct,
                               sync$kendall_tau,
                               sync$kendall_p)
    }
    
    # Create plot with colors
    p <- ggplot(site_temporal, aes(x = Date, y = Biomass, color = Method, group = Method)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5, alpha = 0.8) +
      scale_color_manual(
        values = custom_colors,
        name = "Device Type"
      ) +
      
      # USE PSEUDO-LOG TO HANDLE ZEROS
      scale_y_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 1),
        labels = scales::comma,
        breaks = c(0, 10, 100, 1000, 5000)
      ) +
      labs(
        title = sprintf("Temporal Concordance: FAIRD vs AMMOD (%s)", site),
        subtitle = subtitle_text,
        x = "Date",
        y = "Biomass (mg)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 11, face = "bold"),
        axis.title.y = element_text(size = 11, face = "bold")
      )
    
    temporal_plots[[site]] <- p
    
    # +++ NEW: SAVE INDIVIDUAL PLOTS +++
    ind_filename <- sprintf("outputs/H2_Temporal_%s", site)
    ggsave(paste0(ind_filename, ".png"), plot = p, width  = 7, height = 6, dpi = 300, bg = "white")
    ggsave(paste0(ind_filename, ".pdf"), plot = p, width = 8, height = 5, device = cairo_pdf)
    cat(sprintf("   [OK] Saved individual temporal plot: %s.png\n", ind_filename))
  }
}

# --- Save Combined Temporal Plot ---
combined_temporal <- wrap_plots(temporal_plots, ncol = 1, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10))

print(combined_temporal)

ggsave("outputs/H2_FAIRD_AMMOD_Temporal_Log.png", combined_temporal, width = 10, height = 4 * length(temporal_plots), dpi = 300, bg = "white")
ggsave("outputs/H2_FAIRD_AMMOD_Temporal_Log.pdf", combined_temporal, width = 10, height = 4 * length(temporal_plots), device = cairo_pdf)

cat("\n[OK] Combined temporal dynamics saved (with pseudo-log)\n\n")

# ==============================================================================
# SECTION 5B: COMBINED FIGURES PER SITE (A: scatter | B: temporal concordance)
# ==============================================================================

cat("\n=== SECTION 5B: COMBINED FIGURES PER SITE ===\n\n")

for (site in sites) {
  if (!is.null(plot_list[[site]]) && !is.null(temporal_plots[[site]])) {
    
    scatter_A <- plot_list[[site]] +
      labs(title = sprintf("A) %s (%s)",
                           site,
                           results_table$Environment[results_table$Site == site]))
    
    temporal_B <- temporal_plots[[site]] +
      labs(title = sprintf("B) Temporal Concordance: FAIRD vs AMMOD (%s)", site))
    
    combined_site <- scatter_A + temporal_B +
      patchwork::plot_layout(ncol = 2, widths = c(1, 1.4))
    
    print(combined_site)
    
    out_base <- sprintf("outputs/H2_Combined_%s", site)
    
    ggsave(paste0(out_base, ".png"),
           plot   = combined_site,
           width  = 14, height = 6, dpi = 300, bg = "white")
    
    ggsave(paste0(out_base, ".pdf"),
           plot   = combined_site,
           width  = 14, height = 6, device = cairo_pdf)
    
    cat(sprintf("   [OK] Combined figure saved: %s.png\n", out_base))
  }
}

cat("\n[OK] All per-site combined figures saved\n\n")


# ==============================================================================
# SECTION 6: DEVICE COMPARISON
# ==============================================================================

cat("\n=== SECTION 6: DEVICE COMPARISON ===\n\n")

LETTER_SIZE <- 3.5
LETTER_FONTFACE <- "plain"
LETTER_COLOR <- "black"
Y_OFFSET <- 5 
LETTER_VJUST <- -0.5  # Vertical adjustment for letter position

cat("CUSTOMIZATION SETTINGS:\n")
cat(sprintf("  Letter size: %.1f\n", LETTER_SIZE))
cat(sprintf("  Letter style: %s\n", LETTER_FONTFACE))
cat(sprintf("  Letter color: %s\n", LETTER_COLOR))
cat(sprintf("  Vertical offset: %.1f\n\n", Y_OFFSET))

etraps_daily <- data_list$etraps_daily

device_data <- etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::mutate(Group = paste(Device_type, Site, sep = ":"))

device_summary <- device_data %>%
  dplyr::group_by(Device_type, Site, Group) %>%
  dplyr::summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean_abundance))

letter_assignments <- tribble(
  ~Group,           ~Letters,  ~Rank,
  "ID:Site2",       "d",       1,
  "ID:Site4",       "c",       2,
  "ID:Site1",       "bc",      3,
  "FAIRD:Site4",    "ab",      4,
  "FAIRD:Site1",    "a",       5,
  "FAIRD:Site3",    "a",       6,
  "FAIRD:Site2",    "a",       7
)

device_summary <- device_summary %>%
  dplyr::left_join(
    letter_assignments %>% dplyr::select(Group, Letters),
    by = "Group"
  ) %>%
  dplyr::mutate(
    letter_y = mean_abundance + sd_abundance + Y_OFFSET
  )

cat("Groups with letter assignments:\n")
print(device_summary %>% 
        dplyr::select(Group, mean_abundance, Letters) %>%
        dplyr::arrange(desc(mean_abundance)))
cat("\n")

cat("=== CREATING BOXPLOT ===\n\n")

p_boxplot_letters <- ggplot(device_data, 
                            aes(x = Site, y = abundance, fill = Device_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(
    data = device_summary,
    aes(x = Site, y = letter_y, label = Letters, group = Device_type),
    position = position_dodge(width = 0.8),
    size = LETTER_SIZE,
    fontface = LETTER_FONTFACE,
    color = LETTER_COLOR,
    vjust = LETTER_VJUST
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Distribution of Daily Abundance by Device and Site",
    subtitle = "Letters indicate significance groups (Dunn's test with Bonferroni correction, ± = 0.05).\nGroups sharing letters are not significantly different.",
    x = "Site",
    y = "Daily Abundance (individuals)",
    fill = "Device Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(p_boxplot_letters)
ggsave("outputs/H2_Section6_Boxplot_WithLetters.png", p_boxplot_letters, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("outputs/H2_Section6_Boxplot_WithLetters.pdf", p_boxplot_letters, width = 10, height = 6, device = "pdf")
cat("[OK] Boxplot saved\n\n")

cat("=== CREATING BAR PLOT ===\n\n")

p_bars_letters <- ggplot(device_summary, 
                         aes(x = Site, y = mean_abundance, fill = Device_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_abundance - sd_abundance, 
        ymax = mean_abundance + sd_abundance),
    position = position_dodge(width = 0.8), 
    width = 0.3
  ) +
  geom_text(
    aes(label = Letters, y = letter_y),
    position = position_dodge(width = 0.8),
    size = LETTER_SIZE,
    fontface = LETTER_FONTFACE,
    color = LETTER_COLOR
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Average Daily Abundance by Device Type and Site",
    subtitle = "Error bars: ± 1 SD. Letters indicate significance groups (Dunn's test, Bonferroni, ± = 0.05).",
    x = "Site",
    y = "Mean Daily Abundance (individuals)",
    fill = "Device Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(p_bars_letters)
ggsave("outputs/H2_FAIRD_ID_Devices_WithLetters.png", p_bars_letters, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("outputs/H2_FAIRD_ID_Devices_WithLetters.pdf", p_bars_letters, width = 10, height = 6, device = "pdf")
cat("[OK] Bar plot saved\n\n")

letter_table <- device_summary %>%
  dplyr::select(Device_type, Site, Group, n, mean_abundance, 
                sd_abundance, Letters) %>%
  dplyr::arrange(Device_type, Site)

write_csv(letter_table, "outputs/H2_Section6_Significance_Letters.csv")

# SECTION 6B: THREE-WAY DEVICE COMPARISON (Sites with AMMOD available)
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 6B: THREE-WAY DEVICE COMPARISON\n")
cat("========================================\n\n")

cat("Preparing combined dataset with FAIRD, ID, and AMMOD...\n\n")

# Prepare e-traps data (FAIRD and ID)
etraps_for_comparison <- etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Site, Device_type, abundance) %>%
  dplyr::rename(daily_value = abundance)

# Prepare AMMOD data (aggregate by date and site)
ammod_for_comparison <- ammod_biomass %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    daily_value = sum(Live_mass, na.rm = TRUE), # Using biomass as proxy for abundance
    .groups = "drop"
  ) %>%
  dplyr::mutate(Device_type = "AMMOD")

# Combine all three devices
combined_devices <- dplyr::bind_rows(
  etraps_for_comparison,
  ammod_for_comparison
)

# Identify sites with AMMOD data
sites_with_ammod <- combined_devices %>%
  dplyr::filter(Device_type == "AMMOD") %>%
  dplyr::pull(Site) %>%
  unique() %>%
  sort()

cat("Sites with AMMOD available:", paste(sites_with_ammod, collapse = ", "), "\n\n")

# Filter to only sites with AMMOD
comparison_data <- combined_devices %>%
  dplyr::filter(Site %in% sites_with_ammod)

# Summary statistics
comparison_summary <- comparison_data %>%
  dplyr::group_by(Device_type) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(daily_value, na.rm = TRUE),
    median = median(daily_value, na.rm = TRUE),
    sd = sd(daily_value, na.rm = TRUE),
    min = min(daily_value, na.rm = TRUE),
    max = max(daily_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean))

cat("=== DESCRIPTIVE STATISTICS (Sites with AMMOD) ===\n\n")
print(comparison_summary, n = Inf)
cat("\n")

# NOTE: Inferential comparison (Kruskal-Wallis + Dunn post-hoc) is performed 
# in Section 6C using biomass-converted units (mg/day) for all three devices, 
# enabling statistically valid comparison on a homogeneous scale.

# ==============================================================================
# SECTION 6C: ARITHMETIC vs GEOMETRIC MEANS OF DAILY BIOMASS (mg/day)
# ==============================================================================
# Purpose: Compare arithmetic and geometric means across devices in the same
#          units (mg/day) to contextualise Mixed Model LS means (back-transformed
#          = geometric means) reported in the SAS PROC MIXED output.
# Scope: Same as Section 6B — sites with AMMOD available (Site1, Site2, Site3).
#        FAIRD: n=66 (3 sites × 22 days)
#        ID:    n=44 (2 sites × 22 days; Site3 excluded, ID3 failure)
#        AMMOD: n=32 (Site1 n=12, Site2 n=12, Site3 n=8)
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 6C: ARITHMETIC vs GEOMETRIC MEANS (mg/day)\n")
cat("========================================\n\n")

# --- FAIRD and ID: biomass_mg from etraps_daily (Sites 1-3) ---
etraps_biomass_6c <- data_list$etraps_daily %>%
  dplyr::filter(
    Device_type %in% c("FAIRD", "ID"),
    Site %in% sites_with_ammod  # reuses sites_with_ammod from Section 6B
  ) %>%
  dplyr::select(Date, Site, Device_type, biomass_mg)

# --- AMMOD: Live_mass from ammod_biomass (Sites 1-3) ---
ammod_biomass_6c <- ammod_biomass %>%
  dplyr::filter(
    Site %in% sites_with_ammod,
    !is.na(Live_mass)
  ) %>%
  dplyr::select(Date, Site, Device_type, biomass_mg = Live_mass)

# --- Combine ---
combined_biomass_6c <- dplyr::bind_rows(etraps_biomass_6c, ammod_biomass_6c)

# Calculate totals from the intermediate table
device_totals_6c <- combined_biomass_6c %>%
  group_by(Device_type) %>%
  summarise(
    Total_Biomass_mg = sum(biomass_mg, na.rm = TRUE),
    .groups = 'drop'
  )
cat("Biomass Totals per Device Type\n\n")
print(device_totals_6c, n = Inf)
cat("\n")

# --- Calculate arithmetic and geometric means ---
# Geometric mean uses log(x + 1) to handle zeros, then back-transforms
biomass_means_6c <- combined_biomass_6c %>%
  dplyr::group_by(Device_type) %>%
  dplyr::summarise(
    n                  = dplyr::n(),
    n_nonzero          = sum(biomass_mg > 0),
    arithmetic_mean_mg = mean(biomass_mg, na.rm = TRUE),
    arithmetic_sd_mg   = sd(biomass_mg, na.rm = TRUE),
    geometric_mean_mg  = exp(mean(log(biomass_mg + 1), na.rm = TRUE)) - 1,
    median_mg          = median(biomass_mg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ratio_arith_geom = round(arithmetic_mean_mg / geometric_mean_mg, 2)
  ) %>%
  dplyr::arrange(desc(arithmetic_mean_mg))

cat("=== BIOMASS DESCRIPTIVE STATISTICS (mg/day) ===\n")
cat("Scope: Sites with AMMOD available (Site1, Site2, Site3)\n\n")
print(biomass_means_6c, n = Inf)
cat("\n")

cat("NOTE: Geometric mean = exp(mean(log(x+1))) - 1, handles zero-capture days.\n")
cat("      Mixed Model LS means (SAS PROC MIXED, back-transformed from log10)\n")
cat("      correspond to geometric means and are therefore lower than\n")
cat("      arithmetic means when biomass distributions are right-skewed.\n")
cat("      ratio_arith_geom: how many times larger the arithmetic mean is\n")
cat("      relative to the geometric mean for each device.\n\n")

# --- Save results ---
readr::write_csv(biomass_means_6c, "outputs/H2_Section6C_Biomass_Arithmetic_vs_Geometric.csv")
cat("[OK] Biomass means table saved to: outputs/H2_Section6C_Biomass_Arithmetic_vs_Geometric.csv\n\n")

# --- Kruskal-Wallis on homogeneous units (biomass_mg) ---
cat("=== KRUSKAL-WALLIS TEST (Three Devices, homogeneous units: mg/day) ===\n\n")

kw_biomass_6c <- kruskal.test(biomass_mg ~ Device_type, data = combined_biomass_6c)

n_6c <- nrow(combined_biomass_6c)
k_6c <- length(unique(combined_biomass_6c$Device_type))
effect_size_6c <- calculate_epsilon_squared(kw_biomass_6c, n_6c, k_6c)

cat(sprintf("  • Kruskal-Wallis: χ² = %.3f, df = %d, p = %.4f %s\n",
            kw_biomass_6c$statistic,
            kw_biomass_6c$parameter,
            kw_biomass_6c$p.value,
            ifelse(kw_biomass_6c$p.value < 0.001, "***",
                   ifelse(kw_biomass_6c$p.value < 0.01, "**",
                          ifelse(kw_biomass_6c$p.value < 0.05, "*", "(NS)")))))
cat(sprintf("  • Effect size (ε²): %.3f (%s)\n",
            effect_size_6c$epsilon_squared,
            effect_size_6c$interpretation))
cat(sprintf("    - Variance explained: %.1f%%\n\n",
            effect_size_6c$variance_explained_pct))

# --- Dunn post-hoc ---
cat("=== DUNN'S POST-HOC TEST (Bonferroni correction) ===\n\n")

dunn_6c <- dunn.test::dunn.test(
  x = combined_biomass_6c$biomass_mg,
  g = combined_biomass_6c$Device_type,
  method = "bonferroni",
  kw = FALSE,
  label = TRUE,
  wrap = FALSE,
  table = FALSE,
  list = TRUE,
  altp = FALSE
)

dunn_table_6c <- tibble::tibble(
  Comparison  = dunn_6c$comparisons,
  Z_statistic = round(as.numeric(dunn_6c$Z), 3),
  p_adjusted  = round(as.numeric(dunn_6c$P.adjusted), 4),
  Significant = as.numeric(dunn_6c$P.adjusted) < 0.05
)

print(dunn_table_6c, n = Inf)
cat("\n")

readr::write_csv(dunn_table_6c, "outputs/H2_Section6C_Dunn_Biomass.csv")
cat("[OK] Dunn results saved to: outputs/H2_Section6C_Dunn_Biomass.csv\n\n")

# ==============================================================================
# SECTION 7: WEEKLY PATTERNS - 3 WEEKS OF 7 DAYS EACH
# ==============================================================================

cat("\n=== SECTION 7: WEEKLY PATTERNS ===\n\n")

etraps_weekly <- etraps_daily %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    days_since_start = as.numeric(Date - min(Date)),
    Week = pmin(floor(days_since_start / 7) + 1, 3)
  ) %>%
  dplyr::filter(Week <= 3) %>%
  dplyr::group_by(Week, Device_type, Site) %>%
  dplyr::summarise(w = mean(abundance, na.rm = TRUE), n = n(), .groups = "drop")

cat("Weekly patterns summary:\n")
cat(sprintf("Study period: 3 weeks (%d days)\n", 21))
cat(sprintf("Weeks in data: %s\n", paste(sort(unique(etraps_weekly$Week)), collapse = ", ")))
print(etraps_weekly)

plot_weekly <- ggplot(etraps_weekly, aes(x = factor(Week), y = w, fill = Device_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Site, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = c("1", "2", "3")) +
  labs(
    title = "Weekly Average Abundance by Device Type and Site",
    x = "Week",
    y = "Mean Daily Abundance",
    fill = "Device Type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(plot_weekly)
ggsave("outputs/H2_FAIRD_ID_Weekly.png", plot_weekly, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("outputs/H2_FAIRD_ID_Weekly.pdf", plot_weekly, width = 12, height = 8, device = cairo_pdf)
cat("\n[OK] Weekly patterns saved (3 weeks of 7 days each)\n\n")


# ==============================================================================
# H2 SECTION 8: ENVIRONMENT COMPARISON (IMPROVED VERSION)
# ==============================================================================

cat("\n=== SECTION 8: ENVIRONMENT COMPARISON (IMPROVED) ===\n\n")

environment_data <- etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::mutate(
    Environment = case_when(
      Site %in% c("Site1", "Site2") ~ "Maize",
      Site == "Site4" ~ "Meadow",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(Environment %in% c("Maize", "Meadow"))

env_summary <- environment_data %>%
  dplyr::group_by(Device_type, Environment) %>%
  dplyr::summarise(
    n = n(),
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    median_abundance = median(abundance, na.rm = TRUE),
    q75 = quantile(abundance, 0.75, na.rm = TRUE),
    max_abundance = max(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Environment, desc(mean_abundance))

cat("Summary by environment:\n")
print(env_summary)
cat("\n")

maize_data <- environment_data %>% dplyr::filter(Environment == "Maize")
meadow_data <- environment_data %>% dplyr::filter(Environment == "Meadow")

kw_maize <- kruskal.test(abundance ~ Device_type, data = maize_data)
kw_meadow <- kruskal.test(abundance ~ Device_type, data = meadow_data)

mw_maize <- wilcox.test(abundance ~ Device_type, data = maize_data, exact = FALSE)
mw_meadow <- wilcox.test(abundance ~ Device_type, data = meadow_data, exact = FALSE)

cat("Statistical tests:\n")
cat(sprintf("  Maize:  p = %.4f %s\n", mw_maize$p.value, 
            ifelse(mw_maize$p.value < 0.05, "***", "(NS)")))
cat(sprintf("  Meadow: p = %.4f %s\n\n", mw_meadow$p.value,
            ifelse(mw_meadow$p.value < 0.05, "***", "(NS)")))

calculate_whisker_top <- function(data) {
  if (length(data) == 0 || all(is.na(data))) return(NA)
  q1 <- quantile(data, 0.25, na.rm = TRUE)
  q3 <- quantile(data, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  upper_whisker <- min(max(data[data <= q3 + 1.5 * iqr], na.rm = TRUE), 
                       max(data, na.rm = TRUE))
  return(upper_whisker)
}

env_summary <- env_summary %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    whisker_top = calculate_whisker_top(
      environment_data$abundance[
        environment_data$Device_type == Device_type & 
          environment_data$Environment == Environment
      ]
    )
  ) %>%
  dplyr::ungroup()

env_summary <- env_summary %>%
  dplyr::mutate(
    Letters = case_when(
      Environment == "Maize" & Device_type == "FAIRD" ~ ifelse(mw_maize$p.value < 0.05, "a", "a"),
      Environment == "Maize" & Device_type == "ID" ~ ifelse(mw_maize$p.value < 0.05, "b", "a"),
      Environment == "Meadow" & Device_type == "FAIRD" ~ ifelse(mw_meadow$p.value < 0.05, "a", "a"),
      Environment == "Meadow" & Device_type == "ID" ~ ifelse(mw_meadow$p.value < 0.05, "b", "a"),
      TRUE ~ "?"
    ),
    letter_y = whisker_top + 8
  )

cat("Letter positions:\n")
print(env_summary %>% 
        dplyr::select(Environment, Device_type, whisker_top, letter_y, Letters))
cat("\n")

cat("Creating integrated plot with shared Y-axis...\n\n")

environment_data <- environment_data %>%
  dplyr::mutate(
    Group = paste(Environment, Device_type, sep = "\n"),
    Group = factor(Group, levels = c("Maize\nFAIRD", "Maize\nID", 
                                     "Meadow\nFAIRD", "Meadow\nID"))
  )

env_summary <- env_summary %>%
  dplyr::mutate(
    Group = paste(Environment, Device_type, sep = "\n"),
    Group = factor(Group, levels = c("Maize\nFAIRD", "Maize\nID", 
                                     "Meadow\nFAIRD", "Meadow\nID"))
  )

plot_integrated <- ggplot(environment_data, 
                          aes(x = Group, y = abundance, fill = Device_type)) +
  geom_boxplot(width = 0.6) +
  geom_text(
    data = env_summary,
    aes(x = Group, y = letter_y, label = Letters),
    size = LETTER_SIZE,
    fontface = LETTER_FONTFACE,
    color = LETTER_COLOR,
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_fill_manual(values = custom_colors) +
  annotate("text", x = 1.5, y = Inf, label = "Maize (Site1+Site2)",
           vjust = 1.5, fontface = "bold", size = 4) +
  annotate("text", x = 3.5, y = Inf, label = "Meadow (Site4)",
           vjust = 1.5, fontface = "bold", size = 4) +
  labs(
    title = "Device Performance by Environment Type",
    subtitle = sprintf("Letters indicate significance within each environment (Mann-Whitney U, ± = 0.05).\nMaize: p = %.3f; Meadow: p = %.3f",
                       mw_maize$p.value, mw_meadow$p.value),
    x = NULL,
    y = "Daily Abundance (individuals)",
    fill = "Device Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10)
  )

print(plot_integrated)

ggsave("outputs/H2_Section8_Environment_Comparison.png",
       plot_integrated, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("outputs/H2_Section8_Environment_Comparison.pdf",
       plot_integrated, width = 10, height = 6, device = "pdf")

cat("[OK] Integrated environment plot saved\n\n")

plot_faceted <- ggplot(environment_data %>% 
                         dplyr::mutate(Group_simple = Device_type), 
                       aes(x = Group_simple, y = abundance, fill = Device_type)) +
  geom_boxplot(width = 0.6) +
  geom_text(
    data = env_summary %>% dplyr::mutate(Group_simple = Device_type),
    aes(x = Group_simple, y = letter_y, label = Letters),
    size = LETTER_SIZE,
    fontface = LETTER_FONTFACE,
    color = LETTER_COLOR,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Environment, scales = "free_x") +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Device Performance by Environment Type",
    subtitle = "Letters indicate significance (Mann-Whitney U, ± = 0.05)",
    x = "Device Type",
    y = "Daily Abundance (individuals)",
    fill = "Device"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank()
  )

print(plot_faceted)

ggsave("outputs/H2_Section8_Environment_Faceted.png",
       plot_faceted, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("outputs/H2_Section8_Environment_Faceted.pdf",
       plot_faceted, width = 10, height = 6, device = "pdf")

cat("[OK] Faceted environment plot saved (alternative version)\n\n")

results_env <- tibble(
  Environment = c("Maize", "Meadow"),
  Sites = c("Site1+Site2", "Site4"),
  FAIRD_n = c(
    sum(environment_data$Device_type == "FAIRD" & environment_data$Environment == "Maize"),
    sum(environment_data$Device_type == "FAIRD" & environment_data$Environment == "Meadow")
  ),
  ID_n = c(
    sum(environment_data$Device_type == "ID" & environment_data$Environment == "Maize"),
    sum(environment_data$Device_type == "ID" & environment_data$Environment == "Meadow")
  ),
  FAIRD_mean = env_summary$mean_abundance[env_summary$Device_type == "FAIRD" & env_summary$Environment == "Maize"],
  FAIRD_sd = env_summary$sd_abundance[env_summary$Device_type == "FAIRD" & env_summary$Environment == "Maize"],
  ID_mean = env_summary$mean_abundance[env_summary$Device_type == "ID" & env_summary$Environment == "Maize"],
  ID_sd = env_summary$sd_abundance[env_summary$Device_type == "ID" & env_summary$Environment == "Maize"],
  p_value = c(mw_maize$p.value, mw_meadow$p.value),
  Significant = c(mw_maize$p.value < 0.05, mw_meadow$p.value < 0.05)
)

write_csv(results_env, "outputs/H2_Section8_Environment_Results.csv")

cat("\n================================================================\n")
cat("|                    SUMMARY                                   |\n")
cat("================================================================\n\n")

cat("FILES GENERATED:\n")
cat("  [OK] outputs/H2_Section8_Environment_Comparison.png/pdf (integrated)\n")
cat("  [OK] outputs/H2_Section8_Environment_Faceted.png/pdf (alternative)\n")
cat("  [OK] outputs/H2_Section8_Environment_Results.csv\n\n")

cat("========================================\n")
cat("SECTION 8 COMPLETE (IMPROVED)\n")

# SECTION 8B: GLOBAL HABITAT COMPARISON
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 8B: GLOBAL HABITAT COMPARISON\n")
cat("========================================\n\n")

cat("Preparing habitat comparison dataset...\n\n")

# Create habitat classification for all devices including AMMOD
habitat_data_full <- dplyr::bind_rows(
  # E-traps (FAIRD and ID)
  etraps_daily %>%
    dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
    dplyr::mutate(
      Environment = case_when(
        Site %in% c("Site1", "Site2") ~ "Maize",
        Site == "Site4" ~ "Meadow",
        TRUE ~ "Other"
      ),
      value = abundance
    ) %>%
    dplyr::select(Date, Site, Device_type, Environment, value),
  
  # AMMOD
  ammod_biomass %>%
    dplyr::mutate(
      Environment = case_when(
        Site %in% c("Site1", "Site2") ~ "Maize",
        Site == "Site4" ~ "Meadow",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::group_by(Date, Site, Device_type, Environment) %>%
    dplyr::summarise(value = sum(Live_mass, na.rm = TRUE), .groups = "drop")
) %>%
  dplyr::filter(Environment %in% c("Maize", "Meadow"))

# Create combined factor for habitat × device
habitat_data_full <- habitat_data_full %>%
  dplyr::mutate(
    Habitat_Device = paste(Environment, Device_type, sep = "_"),
    Habitat_Device = factor(Habitat_Device)
  )

# Summary by habitat and device
habitat_summary_full <- habitat_data_full %>%
  dplyr::group_by(Environment, Device_type) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Environment, desc(mean))

cat("=== DESCRIPTIVE STATISTICS BY HABITAT ===\n\n")
print(habitat_summary_full, n = Inf)
cat("\n")

# ==============================================================================
# SECTION 8C: HABITAT × DEVICE COMPARISON (homogeneous units: mg/day)
# ==============================================================================
# Purpose: Corrected version of Section 8B using biomass_mg for all devices,
#          enabling statistically valid Kruskal-Wallis comparison.
#          Section 8B descriptive stats are retained for natural-unit context.
# Scope:
#   Maize_AMMOD: Site1 + Site2, n=24 (12 days each)
#   Maize_FAIRD: Site1 + Site2, n=44 (22 days each)
#   Maize_ID:    Site1 + Site2, n=44 (22 days each)
#   Meadow_FAIRD: Site4, n=22 (22 days)
#   Meadow_ID:    Site4, n=22 (22 days)
#   NOTE: Meadow_AMMOD absent (AMMOD3 n=8 days only, AMMOD4 unavailable)
# ==============================================================================

cat("\n========================================\n")
cat("SECTION 8C: HABITAT × DEVICE COMPARISON (homogeneous units: mg/day)\n")
cat("========================================\n\n")

# --- Build homogeneous dataset ---
habitat_biomass_8c <- dplyr::bind_rows(
  # FAIRD and ID: biomass_mg (individuals weighted by allometric mass)
  data_list$etraps_daily %>%
    dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
    dplyr::mutate(
      Environment = dplyr::case_when(
        Site %in% c("Site1", "Site2") ~ "Maize",
        Site == "Site4"               ~ "Meadow",
        TRUE                          ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Environment)) %>%
    dplyr::select(Date, Site, Device_type, Environment, biomass_mg),
  
  # AMMOD: Live_mass (gravimetric, mg/day) — Maize only (Site1-2)
  ammod_biomass %>%
    dplyr::filter(!is.na(Live_mass)) %>%
    dplyr::mutate(
      Environment = dplyr::case_when(
        Site %in% c("Site1", "Site2") ~ "Maize",
        TRUE                          ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Environment)) %>%
    dplyr::group_by(Date, Site, Device_type, Environment) %>%
    dplyr::summarise(biomass_mg = sum(Live_mass, na.rm = TRUE), .groups = "drop")
) %>%
  dplyr::mutate(
    Habitat_Device = paste(Environment, Device_type, sep = "_"),
    Habitat_Device = factor(Habitat_Device)
  )

# --- Descriptive statistics ---
habitat_summary_8c <- habitat_biomass_8c %>%
  dplyr::group_by(Environment, Device_type, Habitat_Device) %>%
  dplyr::summarise(
    n          = dplyr::n(),
    n_nonzero  = sum(biomass_mg > 0),
    mean_mg    = mean(biomass_mg, na.rm = TRUE),
    median_mg  = median(biomass_mg, na.rm = TRUE),
    sd_mg      = sd(biomass_mg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Environment, dplyr::desc(mean_mg))

cat("=== DESCRIPTIVE STATISTICS BY HABITAT (mg/day) ===\n\n")
print(habitat_summary_8c, n = Inf)
cat("\n")

# --- Kruskal-Wallis ---
cat("=== KRUSKAL-WALLIS TEST (Habitat × Device, homogeneous units) ===\n\n")

kw_8c <- kruskal.test(biomass_mg ~ Habitat_Device, data = habitat_biomass_8c)

n_8c <- nrow(habitat_biomass_8c)
k_8c <- length(unique(habitat_biomass_8c$Habitat_Device))
effect_8c <- calculate_epsilon_squared(kw_8c, n_8c, k_8c)

cat(sprintf("  • Kruskal-Wallis: χ² = %.3f, df = %d, p = %.4f %s\n",
            kw_8c$statistic,
            kw_8c$parameter,
            kw_8c$p.value,
            ifelse(kw_8c$p.value < 0.001, "***",
                   ifelse(kw_8c$p.value < 0.01, "**",
                          ifelse(kw_8c$p.value < 0.05, "*", "(NS)")))))
cat(sprintf("  • Effect size (ε²): %.3f (%s)\n",
            effect_8c$epsilon_squared,
            effect_8c$interpretation))
cat(sprintf("    - Variance explained: %.1f%%\n\n",
            effect_8c$variance_explained_pct))

# --- Dunn post-hoc ---
cat("=== DUNN'S POST-HOC TEST (Bonferroni correction) ===\n\n")

dunn_8c <- dunn.test::dunn.test(
  x      = habitat_biomass_8c$biomass_mg,
  g      = habitat_biomass_8c$Habitat_Device,
  method = "bonferroni",
  kw     = FALSE,
  label  = TRUE,
  wrap   = FALSE,
  table  = FALSE,
  list   = TRUE,
  altp   = FALSE
)

dunn_table_8c <- tibble::tibble(
  Comparison  = dunn_8c$comparisons,
  Z_statistic = round(as.numeric(dunn_8c$Z), 3),
  p_adjusted  = round(as.numeric(dunn_8c$P.adjusted), 4),
  Significant = as.numeric(dunn_8c$P.adjusted) < 0.05
) %>%
  dplyr::arrange(p_adjusted)

cat("=== FULL PAIRWISE COMPARISONS TABLE (sorted by p-value) ===\n\n")
print(dunn_table_8c, n = Inf)
cat("\n")

readr::write_csv(dunn_table_8c, "outputs/H2_Section8C_Dunn_Habitat_Biomass.csv")
cat("[OK] Dunn results saved to: outputs/H2_Section8C_Dunn_Habitat_Biomass.csv\n\n")

# --- Save KW summary ---
kw_summary_8c <- tibble::tibble(
  Test                   = "Kruskal-Wallis",
  Statistic              = round(kw_8c$statistic, 3),
  df                     = kw_8c$parameter,
  p_value                = kw_8c$p.value,
  epsilon_squared        = round(effect_8c$epsilon_squared, 3),
  effect_interpretation  = effect_8c$interpretation,
  variance_explained_pct = round(effect_8c$variance_explained_pct, 1),
  Significant            = kw_8c$p.value < 0.05
)

readr::write_csv(kw_summary_8c, "outputs/H2_Section8C_KruskalWallis_Habitat_Biomass.csv")
cat("[OK] KW results saved to: outputs/H2_Section8C_KruskalWallis_Habitat_Biomass.csv\n\n")

# ==============================================================================
cat("========================================\n\n")

cat("\n========================================\n")
cat("H2 COMPLETE\n")
cat("========================================\n\n")#

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")