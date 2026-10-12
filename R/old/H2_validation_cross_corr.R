# ==============================================================================
# H2: TEMPORAL CROSS-CORRELATION ANALYSIS (MULTIPLE TIME SERIES)
# ==============================================================================
# Created: 2025-12-09
# Purpose: Generate cross-correlation plots between FAIRD and AMMOD time series
#          Based on Crawley's "The R Book" section on Multiple Time Series
# Reference: Crawley (2013), Chapter 24.8 - Multiple time series
# ==============================================================================

library(tidyverse)
library(patchwork)

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
cat("H2: TEMPORAL CROSS-CORRELATION ANALYSIS\n")
cat("========================================\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("=== LOADING DATA ===\n\n")

# Load custom colors
source("custom_colors.R")

# Load cleaned data
load("data_cleaned.RData")

etraps_daily_overlap <- data_list$etraps_daily_overlap
ammod_biomass_overlap <- data_list$ammod_biomass_overlap

# ==============================================================================
# DATA PREPARATION FOR TIME SERIES ANALYSIS
# ==============================================================================

cat("=== PREPARING TIME SERIES DATA ===\n\n")

# Create paired dataset for time series analysis
# FAIRD: inferred biomass (mg live mass, allometric inference — consistent with
#         H2_validation_FAIRD_AMMOD.R and H2_validation_aggregated_levels.R)
# AMMOD: gravimetric live mass (mg, converted from dry mass)
# Both series are on the same biomass scale prior to CCF computation.
# NOTE: raw abundance (counts) is NOT used for FAIRD — using biomass_mg is
# consistent with the OLS validation and avoids mixing measurement scales.

ts_data <- etraps_daily_overlap %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    faird_value = sum(biomass_mg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::inner_join(
    ammod_biomass_overlap %>%
      dplyr::select(Date, Site, ammod_value = Live_mass),
    by = c("Date", "Site")
  ) %>%
  dplyr::arrange(Site, Date)

# Check data availability by site
sites_available <- ts_data %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    n_days = n(),
    date_range = paste(min(Date), "to", max(Date)),
    faird_min = min(faird_value),
    faird_max = max(faird_value),
    faird_mean = round(mean(faird_value), 1),
    ammod_min = round(min(ammod_value), 2),
    ammod_max = round(max(ammod_value), 2),
    ammod_mean = round(mean(ammod_value), 2),
    .groups = "drop"
  )

cat("Data availability by site:\n\n")
print(sites_available, n = Inf)
cat("\n")

cat("Data ranges:\n")
cat(sprintf("  Total site-days: %d\n", nrow(ts_data)))
cat(sprintf("  Sites included: %s\n", paste(unique(ts_data$Site), collapse = ", ")))
cat(sprintf("  Date range: %s to %s\n", min(ts_data$Date), max(ts_data$Date)))
cat("\n")

# ==============================================================================
# FUNCTION: PLOT TIME SERIES (Crawley-style stacked plot)
# ==============================================================================

plot_timeseries_stacked <- function(site_data, site_name) {
  # Create stacked time series plot similar to Crawley's plot.ts(cbind(x,y))
  # Arguments:
  #   site_data: data frame with Date, faird_value, ammod_value
  #   site_name: character string for plot title
  
  if(nrow(site_data) < 3) {
    cat(sprintf("   Skipping %s: Insufficient data (n=%d)\n", site_name, nrow(site_data)))
    return(NULL)
  }
  
  # Prepare data for plotting
  plot_data <- site_data %>%
    dplyr::select(Date, faird_value, ammod_value) %>%
    tidyr::pivot_longer(
      cols = c(faird_value, ammod_value),
      names_to = "Device",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Device = dplyr::case_when(
        Device == "faird_value" ~ "FAIRD (biomass mg)",
        Device == "ammod_value" ~ "AMMOD (biomass mg)"
      ),
      Device = factor(Device, levels = c("FAIRD (biomass mg)", "AMMOD (biomass mg)"))
    )
  
  # Define colors specifically for this plot (avoid warning)
  plot_colors <- c(
    "FAIRD (biomass mg)" = unname(custom_colors["FAIRD"]),
    "AMMOD (biomass mg)" = unname(custom_colors["AMMOD"])
  )
  
  # Create stacked time series plot
  p <- ggplot(plot_data, aes(x = Date, y = Value, color = Device)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~Device, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = plot_colors) +
    labs(
      title = paste("Time Series Comparison:", site_name),
      subtitle = sprintf("n = %d days", nrow(site_data)),
      x = "Date",
      y = "Value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  return(p)
}

# ==============================================================================
# FUNCTION: CROSS-CORRELATION ANALYSIS (Crawley-style ACF plot)
# ==============================================================================

analyze_crosscorrelation <- function(site_data, site_name) {
  # Perform cross-correlation analysis using ACF (NOT PACF)
  # ACF measures simple correlation between series
  # PACF would partial out autocorrelation effects (too conservative here)
  # 
  # CRITICAL: Standardize series to have mean=0, sd=1 before analysis
  # This prevents scale differences from affecting correlation estimates
  # Arguments:
  #   site_data: data frame with Date, faird_value, ammod_value
  #   site_name: character string for identification
  # Returns:
  #   List with ACF results and plot
  
  if(nrow(site_data) < 3) {
    cat(sprintf("   Skipping %s: Insufficient data (n=%d)\n", site_name, nrow(site_data)))
    return(NULL)
  }
  
  # Extract time series as vectors (ordered by Date)
  site_data <- site_data %>% dplyr::arrange(Date)
  faird_ts_raw <- site_data$faird_value
  ammod_ts_raw <- site_data$ammod_value
  
  # STANDARDIZE series (mean=0, sd=1) to make scales comparable
  # This is ESSENTIAL when series have very different scales
  # FAIRD: counts (0-100), AMMOD: biomass mg (500-7000)
  faird_ts <- scale(faird_ts_raw)[,1]  # Extract vector from matrix
  ammod_ts <- scale(ammod_ts_raw)[,1]
  
  # Convert to time series objects
  faird_ts_obj <- ts(faird_ts)
  ammod_ts_obj <- ts(ammod_ts)
  
  # Calculate ACF (NOT PACF) for multivariate series
  # This gives simple cross-correlation without partialing effects
  # Following Crawley (2013) The R Book, Chapter 24.8
  acf_results <- acf(cbind(faird_ts_obj, ammod_ts_obj), 
                     plot = FALSE, 
                     lag.max = min(10, floor(nrow(site_data)/3)))
  
  # Extract results for interpretation
  # acf_results$acf is a 3D array: [lag, series, series]
  # [lag, 1, 1] = FAIRD autocorrelation
  # [lag, 2, 2] = AMMOD autocorrelation  
  # [lag, 1, 2] = FAIRD & AMMOD cross-correlation
  # [lag, 2, 1] = AMMOD & FAIRD cross-correlation
  
  # NOTE: ACF includes lag 0, PACF does not
  # So acf[1,,] = lag 0, acf[2,,] = lag 1, etc.
  
  # Find significant lags (beyond confidence bounds)
  ci_level <- 0.05
  ci_bound <- qnorm(1 - ci_level/2) / sqrt(nrow(site_data))
  
  # Extract cross-correlation values
  # Skip lag 0 (index 1) for finding "lead-lag" relationships
  # But we'll report lag 0 separately as it's most important
  faird_ammod_ccf <- acf_results$acf[, 1, 2]  # FAIRD & AMMOD (all lags including 0)
  ammod_faird_ccf <- acf_results$acf[, 2, 1]  # AMMOD & FAIRD (all lags including 0)
  
  # Build explicit lag table with CCF values and CI for all lags
  n_lags     <- length(faird_ammod_ccf)
  lag_indices <- 0:(n_lags - 1)  # Actual lag values (0, 1, 2, ...)
  
  lag_table <- tibble::tibble(
    lag      = lag_indices,
    ccf_FA   = faird_ammod_ccf,
    ccf_AF   = ammod_faird_ccf,
    ci_bound = ci_bound,
    sig_FA   = abs(faird_ammod_ccf) > ci_bound,
    sig_AF   = abs(ammod_faird_ccf) > ci_bound
  )
  # Identify significant lags (excluding lag 0)
  sig_lags_fa <- which(abs(faird_ammod_ccf[-1]) > ci_bound)  # Remove lag 0, then find sig
  sig_lags_af <- which(abs(ammod_faird_ccf[-1]) > ci_bound)
  
  list(
    site      = site_name,
    n         = nrow(site_data),
    acf_obj   = acf_results,
    ci_bound  = ci_bound,
    ccf_lag0  = faird_ammod_ccf[1],
    lag_table = lag_table,
    sig_lags_faird_to_ammod = sig_lags_fa,
    sig_lags_ammod_to_faird = sig_lags_af,
    max_ccf_fa = ifelse(length(sig_lags_fa) > 0,
                        max(abs(faird_ammod_ccf[-1][sig_lags_fa])), 0),
    max_ccf_af = ifelse(length(sig_lags_af) > 0,
                        max(abs(ammod_faird_ccf[-1][sig_lags_af])), 0)
  )
}

# ==============================================================================
# FUNCTION: PLOT CROSS-CORRELATION (4-panel Crawley-style)
# ==============================================================================

plot_crosscorrelation <- function(ccf_result, site_name) {
  # Create 4-panel cross-correlation plot similar to Crawley's style
  # Arguments:
  #   ccf_result: output from analyze_crosscorrelation()
  #   site_name: character string for plot title
  
  if(is.null(ccf_result)) return(NULL)
  
  # Create temporary PNG file for base R plotting
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, width = 1200, height = 900, res = 120)
  
  # Set up 2x2 panel layout
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  # Get PACF object
  pacf_obj <- ccf_result$pacf_obj
  
  # Panel titles
  titles <- c(
    "FAIRD", 
    "FAIRD & AMMOD",
    "AMMOD & FAIRD",
    "AMMOD"
  )
  
  # Colors
  bar_color <- "red"
  
  # Plot 4 panels
  for(i in 1:4) {
    if(i == 1) {
      # FAIRD autocorrelation [lag, 1, 1]
      acf_vals <- pacf_obj$acf[, 1, 1]
    } else if(i == 2) {
      # FAIRD -> AMMOD cross-correlation [lag, 1, 2]
      acf_vals <- pacf_obj$acf[, 1, 2]
    } else if(i == 3) {
      # AMMOD -> FAIRD cross-correlation [lag, 2, 1]
      acf_vals <- pacf_obj$acf[, 2, 1]
    } else {
      # AMMOD autocorrelation [lag, 2, 2]
      acf_vals <- pacf_obj$acf[, 2, 2]
    }
    
    # Plot
    lags <- 1:length(acf_vals)
    plot(lags, acf_vals, 
         type = "h", 
         col = bar_color, 
         lwd = 2,
         ylim = c(-1, 1),
         xlab = "Lag", 
         ylab = "Partial ACF",
         main = titles[i],
         las = 1)
    
    # Add confidence bounds (dashed horizontal lines)
    abline(h = c(-ccf_result$ci_bound, ccf_result$ci_bound), 
           lty = 2, 
           col = "blue")
    abline(h = 0, col = "black")
  }
  
  # Add overall title
  mtext(paste("Cross-Correlation Analysis:", site_name, 
              sprintf("(n = %d days)", ccf_result$n)),
        side = 3, line = -2, outer = TRUE, cex = 1.2, font = 2)
  
  dev.off()
  
  # Read back the image
  img <- png::readPNG(temp_file)
  unlink(temp_file)
  
  # Return as grid object for display
  grid::grid.raster(img)
  
  return(img)
}

# ==============================================================================
# SECTION 1: TIME SERIES PLOTS (STACKED)
# ==============================================================================

cat("\n=== SECTION 1: STACKED TIME SERIES PLOTS ===\n\n")

sites <- unique(ts_data$Site)
ts_plots <- list()

for(site in sites) {
  cat(sprintf("Generating time series plot for %s...\n", site))
  
  site_data <- ts_data %>%
    dplyr::filter(Site == site)
  
  p <- plot_timeseries_stacked(site_data, site)
  
  if(!is.null(p)) {
    ts_plots[[site]] <- p
    
    # Save individual plot
    ggsave(
      filename = sprintf("outputs/H2_CrossCorr_TimeSeries_%s.png", site),
      plot = p,
      width = 10, 
      height = 6, 
      dpi = 300,
      bg = "white"
    )
    
    ggsave(
      filename = sprintf("outputs/H2_CrossCorr_TimeSeries_%s.pdf", site),
      plot = p,
      width = 10, 
      height = 6,
      device = "pdf"
    )
    
    cat(sprintf("   [OK] Saved: outputs/H2_CrossCorr_TimeSeries_%s.png/pdf\n", site))
    
    # Display plot in RStudio
    print(p)
  }
}

cat("\n")

# Create combined plot with all sites
if(length(ts_plots) > 0) {
  cat("Creating combined time series plot...\n")
  
  combined_ts <- patchwork::wrap_plots(ts_plots, ncol = 2) +
    patchwork::plot_annotation(
      title = "FAIRD-AMMOD Time Series Comparison (All Sites)",
      theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    )
  
  ggsave(
    "outputs/H2_CrossCorr_TimeSeries_Combined.png",
    combined_ts,
    width = 16,
    height = 12,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    "outputs/H2_CrossCorr_TimeSeries_Combined.pdf",
    combined_ts,
    width = 16,
    height = 12,
    device = "pdf"
  )
  
  cat("   [OK] Saved: outputs/H2_CrossCorr_TimeSeries_Combined.png/pdf\n\n")
  
  # Display combined plot in RStudio
  print(combined_ts)
}

# ==============================================================================
# SECTION 2: CROSS-CORRELATION ANALYSIS
# ==============================================================================

cat("\n=== SECTION 2: CROSS-CORRELATION ANALYSIS ===\n\n")

ccf_results <- list()

for(site in sites) {
  cat(sprintf("--- %s ---\n", site))
  
  site_data <- ts_data %>%
    dplyr::filter(Site == site)
  
  ccf_res <- analyze_crosscorrelation(site_data, site)
  
  if(!is.null(ccf_res)) {
    ccf_results[[site]] <- ccf_res
    
    # Print summary
    cat(sprintf("   n = %d days\n", ccf_res$n))
    cat(sprintf("   Confidence bound: ±%.3f\n", ccf_res$ci_bound))
    
    # Display lag 0 (same-day synchrony) - most important for validation
    # With ACF, this is stored directly in ccf_lag0 field
    cat(sprintf("   >>> LAG 0 (same-day synchrony): CCF = %.3f %s\n", 
                ccf_res$ccf_lag0,
                ifelse(abs(ccf_res$ccf_lag0) > ccf_res$ci_bound, "★ SIGNIFICANT", "(ns)")))
    
    if(length(ccf_res$sig_lags_faird_to_ammod) > 0) {
      cat(sprintf("   Significant lags beyond lag-0 (FAIRD & AMMOD): %s\n",
                  paste(ccf_res$sig_lags_faird_to_ammod, collapse = ", ")))
      cat(sprintf("   Max |CCF| (FAIRD & AMMOD, excluding lag 0): %.3f\n", ccf_res$max_ccf_fa))
    } else {
      cat("   No significant lags beyond lag-0 (FAIRD & AMMOD)\n")
    }
    
    if(length(ccf_res$sig_lags_ammod_to_faird) > 0) {
      cat(sprintf("   Significant lags beyond lag-0 (AMMOD & FAIRD): %s\n",
                  paste(ccf_res$sig_lags_ammod_to_faird, collapse = ", ")))
      cat(sprintf("   Max |CCF| (AMMOD & FAIRD, excluding lag 0): %.3f\n", ccf_res$max_ccf_af))
    } else {
      cat("   No significant lags beyond lag-0 (AMMOD & FAIRD)\n")
    }
    
    cat("\n")
  }
}

# ==============================================================================
# SECTION 3: CROSS-CORRELATION PLOTS (4-PANEL CRAWLEY-STYLE)
# ==============================================================================

cat("\n=== SECTION 3: GENERATING CROSS-CORRELATION PLOTS ===\n\n")

for(site in sites) {
  if(!is.null(ccf_results[[site]])) {
    cat(sprintf("Generating cross-correlation plot for %s...\n", site))
    
    # First, display in RStudio (normal graphics device)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    
    ccf_obj <- ccf_results[[site]]
    acf_obj <- ccf_obj$acf_obj  # Now using ACF, not PACF
    
    titles <- c("FAIRD", "FAIRD & AMMOD", "AMMOD & FAIRD", "AMMOD")
    
    for(i in 1:4) {
      if(i == 1) {
        # FAIRD autocorrelation
        acf_vals <- acf_obj$acf[, 1, 1]
      } else if(i == 2) {
        # FAIRD & AMMOD cross-correlation
        acf_vals <- acf_obj$acf[, 1, 2]
      } else if(i == 3) {
        # AMMOD & FAIRD cross-correlation
        acf_vals <- acf_obj$acf[, 2, 1]
      } else {
        # AMMOD autocorrelation
        acf_vals <- acf_obj$acf[, 2, 2]
      }
      
      # ACF includes lag 0 at index 1, so lags are 0, 1, 2, 3, ...
      lags <- 0:(length(acf_vals)-1)
      plot(lags, acf_vals, 
           type = "h", 
           col = "red", 
           lwd = 2,
           ylim = c(-1, 1),
           xlab = "Lag", 
           ylab = "ACF",  # Changed from "Partial ACF"
           main = titles[i],
           las = 1)
      
      abline(h = c(-ccf_obj$ci_bound, ccf_obj$ci_bound), 
             lty = 2, 
             col = "blue")
      abline(h = 0, col = "black")
    }
    
    mtext(paste("Cross-Correlation Analysis:", site, 
                sprintf("(n = %d days)", ccf_obj$n)),
          side = 3, line = -2, outer = TRUE, cex = 1.2, font = 2)
    
    # Now save to file
    png(
      filename = sprintf("outputs/H2_CrossCorr_ACF_%s.png", site),  # Changed from PACF
      width = 1200, 
      height = 900, 
      res = 120
    )
    
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    
    for(i in 1:4) {
      if(i == 1) {
        acf_vals <- acf_obj$acf[, 1, 1]
      } else if(i == 2) {
        acf_vals <- acf_obj$acf[, 1, 2]
      } else if(i == 3) {
        acf_vals <- acf_obj$acf[, 2, 1]
      } else {
        acf_vals <- acf_obj$acf[, 2, 2]
      }
      
      lags <- 0:(length(acf_vals)-1)
      plot(lags, acf_vals, 
           type = "h", 
           col = "red", 
           lwd = 2,
           ylim = c(-1, 1),
           xlab = "Lag", 
           ylab = "ACF",
           main = titles[i],
           las = 1)
      
      abline(h = c(-ccf_obj$ci_bound, ccf_obj$ci_bound), 
             lty = 2, 
             col = "blue")
      abline(h = 0, col = "black")
    }
    
    mtext(paste("Cross-Correlation Analysis:", site, 
                sprintf("(n = %d days)", ccf_obj$n)),
          side = 3, line = -2, outer = TRUE, cex = 1.2, font = 2)
    
    dev.off()
    
    cat(sprintf("   [OK] Saved: outputs/H2_CrossCorr_ACF_%s.png\n", site))
  }
}

# Reset graphics parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

cat("\n")

# ==============================================================================
# SECTION 3B: COMBINED PUBLICATION FIGURES (A: time series | B: CCF FAIRD→AMMOD)
# ==============================================================================

cat("\n=== SECTION 3B: COMBINED PUBLICATION FIGURES ===\n\n")

for (site in sites) {
  if (is.null(ts_plots[[site]]) || is.null(ccf_results[[site]])) next
  
  # --- Build panel B: ggplot CCF (FAIRD→AMMOD only, lag_table already computed) ---
  lt      <- ccf_results[[site]]$lag_table
  ci      <- ccf_results[[site]]$ci_bound
  n_days  <- ccf_results[[site]]$n
  lag0_val <- ccf_results[[site]]$ccf_lag0
  lag0_sig <- abs(lag0_val) > ci
  
  plot_ccf <- ggplot(lt, aes(x = lag, y = ccf_FA)) +
    geom_col(
      aes(fill = sig_FA),
      width = 0.6,
      color = "gray30",
      linewidth = 0.3
    ) +
    scale_fill_manual(
      values = c("TRUE"  = unname(custom_colors["FAIRD"]),
                 "FALSE" = "gray70"),
      labels = c("TRUE" = "Significant", "FALSE" = "ns"),
      name   = NULL
    ) +
    geom_hline(yintercept =  ci, linetype = "dashed",
               color = "dodgerblue", linewidth = 0.8) +
    geom_hline(yintercept = -ci, linetype = "dashed",
               color = "dodgerblue", linewidth = 0.8) +
    geom_hline(yintercept = 0,   color = "black",
               linewidth = 0.4) +
    scale_x_continuous(breaks = lt$lag) +
    ylim(-1, 1) +
    labs(
      title    = sprintf("B) Cross-Correlation: FAIRD \u2192 AMMOD (%s)", site),
      subtitle = sprintf("Lag-0 CCF = %.3f%s | n = %d days | CI = \u00b1%.3f",
                         lag0_val,
                         ifelse(lag0_sig, " *", " ns"),
                         n_days, ci),
      x = "Lag (days)",
      y = "Cross-correlation (ACF)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "gray80", fill = NA,
                                      linewidth = 0.5)
    )
  
  # --- Relabel panel A ---
  plot_ts_A <- ts_plots[[site]] +
    labs(title = sprintf("A) Time Series: FAIRD vs AMMOD (%s)", site))
  
  # --- Combine A + B ---
  combined_ccf <- plot_ts_A + plot_ccf +
    patchwork::plot_layout(ncol = 2, widths = c(1.4, 1))
  
  print(combined_ccf)
  
  out_base <- sprintf("outputs/H2_CrossCorr_Combined_%s", site)
  
  ggsave(paste0(out_base, ".png"),
         plot   = combined_ccf,
         width  = 14, height = 6, dpi = 300, bg = "white")
  
  ggsave(paste0(out_base, ".pdf"),
         plot   = combined_ccf,
         width  = 14, height = 6, device = cairo_pdf)
  
  cat(sprintf("   [OK] Combined figure saved: %s.png/.pdf\n", out_base))
}

cat("\n[OK] All combined CCF figures saved\n\n")


# ==============================================================================
# SECTION 4: AGGREGATED LEVEL CCF ANALYSIS (INTRA-HABITAT AND OVERALL)
# ==============================================================================
# Standard approach (ArcGIS Pro documentation; Zuur et al. 2007):
#   Compute CCF independently per site, then average the CCF values across
#   sites for each lag. This preserves the integrity of each temporal series
#   and avoids the artificial autocorrelation introduced by stacking.
#
# No formal significance test is reported for the averaged CCF — consistent
# with the standard: p-values and CIs are only valid for individual pairwise
# series, not for averages across multiple series.
#
# Lag availability:
#   Sites 1 & 2: lag.max = 4 (n=12 days each)
#   Site 3     : lag.max = 2 (n=8 days)
#   → Intra-Habitat Maize average: lags 0-4 (Sites 1+2 only)
#   → Overall average            : lags 0-2 (all three sites; limited by Site3)
# ==============================================================================

cat("\n=== SECTION 4: AGGREGATED LEVEL CCF ANALYSIS (averaged across sites) ===\n\n")
cat("Method: mean CCF per lag across sites (standard approach).\n")
cat("No formal significance test for averaged CCF values.\n\n")

# ------------------------------------------------------------------------------
# SECTION 4A: EXTRACT PER-SITE LAG TABLES FROM ccf_results
# ------------------------------------------------------------------------------

cat("--- 4A: Per-site lag tables ---\n\n")

# Build a single long data frame: one row per (site, lag)
site_lag_long <- purrr::map_dfr(names(ccf_results), function(s) {
  ccf_results[[s]]$lag_table %>%
    dplyr::mutate(
      Site    = s,
      Habitat = dplyr::case_when(
        s %in% c("Site1", "Site2") ~ "Maize",
        s == "Site3"               ~ "Meadow",
        TRUE                       ~ NA_character_
      )
    )
})

cat("Per-site CCF values (FAIRD -> AMMOD, lag 0 only):\n")
site_lag_long %>%
  dplyr::filter(lag == 0) %>%
  dplyr::select(Site, Habitat, lag, ccf_FA) %>%
  print()
cat("\n")

# ------------------------------------------------------------------------------
# SECTION 4B: COMPUTE AVERAGED CCF PER LEVEL (Fisher's z-transformation)
# ------------------------------------------------------------------------------
# Fisher (1921): raw Pearson r values cannot be arithmetically averaged because
# their sampling distribution is skewed near ±1. The correct procedure is:
#   1. Transform each r_s to z_s = arctanh(r_s)  [variance-stabilizing]
#   2. Average z_s across sites
#   3. Back-transform: r_avg = tanh(mean_z)
#
# Reference: Fisher (1921); applied to multi-site CCF aggregation per Gemini
# Deep Research synthesis (2026-03-21) and time series windowing literature
# (Windowed Cross-Lagged Correlation, WCLC).
#
# Note: For our values (r = 0.382, 0.696, -0.010) the numerical difference
# vs. raw averaging is small but the method is statistically correct.
# ------------------------------------------------------------------------------

cat("--- 4B: Averaged CCF computation (Fisher's z-transformation) ---\n\n")

# Helper: Fisher z-transform and back-transform
fisher_z    <- function(r) atanh(r)          # arctanh: r -> z
fisher_back <- function(z) tanh(z)           # tanh:    z -> r

# Add Fisher z-transformed columns to site_lag_long
site_lag_long <- site_lag_long %>%
  dplyr::mutate(
    z_FA = fisher_z(ccf_FA),
    z_AF = fisher_z(ccf_AF)
  )

# ------------------------------------------------------------------------------
# Intra-Habitat Level (Maize): Fisher's z average of Site1 and Site2, lags 0-4
# NOTE: For Intra-Habitat Maize the SUMMED approach (Section 4B) is the primary
#       method. This averaging is retained for methodological comparison and
#       for the Overall Level where summing across habitats is not appropriate.
# ------------------------------------------------------------------------------
maize_avg <- site_lag_long %>%
  dplyr::filter(Habitat == "Maize") %>%
  dplyr::group_by(lag) %>%
  dplyr::summarise(
    mean_z_FA   = mean(z_FA),
    mean_z_AF   = mean(z_AF),
    mean_ccf_FA = fisher_back(mean(z_FA)),   # back-transformed
    mean_ccf_AF = fisher_back(mean(z_AF)),
    n_sites     = dplyr::n(),
    .groups     = "drop"
  ) %>%
  dplyr::arrange(lag)

# ------------------------------------------------------------------------------
# Overall Level: Fisher's z average of Sites 1+2+3, lags 0-2
# Reported as a descriptive metric only (no formal CI); not interpretable as
# a primary validation result given the inconclusive Site 3 (n=8, CCF≈0).
# ------------------------------------------------------------------------------
max_lag_overall <- min(purrr::map_dbl(names(ccf_results),
                                      ~max(ccf_results[[.x]]$lag_table$lag)))

overall_avg <- site_lag_long %>%
  dplyr::filter(lag <= max_lag_overall) %>%
  dplyr::group_by(lag) %>%
  dplyr::summarise(
    mean_z_FA   = mean(z_FA),
    mean_z_AF   = mean(z_AF),
    mean_ccf_FA = fisher_back(mean(z_FA)),
    mean_ccf_AF = fisher_back(mean(z_AF)),
    n_sites     = dplyr::n(),
    .groups     = "drop"
  ) %>%
  dplyr::arrange(lag)

# ------------------------------------------------------------------------------
# Print summaries with both raw-average and Fisher's z for transparency
# ------------------------------------------------------------------------------

cat("Intra-Habitat Level (Maize) — averaged CCF via Fisher's z (Sites 1+2):\n\n")
cat(sprintf("  Individual site values:\n"))
cat(sprintf("    Site1: r = %.3f  ->  z = %.3f\n",
            ccf_results[["Site1"]]$ccf_lag0,
            fisher_z(ccf_results[["Site1"]]$ccf_lag0)))
cat(sprintf("    Site2: r = %.3f  ->  z = %.3f\n",
            ccf_results[["Site2"]]$ccf_lag0,
            fisher_z(ccf_results[["Site2"]]$ccf_lag0)))
cat(sprintf("    Mean z = %.3f  ->  r_avg (Fisher) = %.3f\n",
            maize_avg$mean_z_FA[1], maize_avg$mean_ccf_FA[1]))
cat(sprintf("    Raw arithmetic mean r = %.3f  (bias vs Fisher: %.4f)\n\n",
            mean(c(ccf_results[["Site1"]]$ccf_lag0,
                   ccf_results[["Site2"]]$ccf_lag0)),
            mean(c(ccf_results[["Site1"]]$ccf_lag0,
                   ccf_results[["Site2"]]$ccf_lag0)) - maize_avg$mean_ccf_FA[1]))

cat(sprintf("  %-5s  %-12s  %-12s  %-12s  %s\n",
            "Lag", "Mean z(FA)", "r_avg(FA)", "r_avg(AF)", "n_sites"))
for (i in seq_len(nrow(maize_avg))) {
  cat(sprintf("  %-5d  %+.3f        %+.3f        %+.3f        %d\n",
              maize_avg$lag[i],
              maize_avg$mean_z_FA[i],
              maize_avg$mean_ccf_FA[i],
              maize_avg$mean_ccf_AF[i],
              maize_avg$n_sites[i]))
}
cat("  Note: no formal CI — descriptive metric (averaged series).\n\n")

cat(sprintf("Overall Level — averaged CCF via Fisher's z (Sites 1+2+3, lags 0-%d):\n\n",
            max_lag_overall))
cat(sprintf("  Individual site values (lag 0):\n"))
for (s in names(ccf_results)) {
  r_s <- ccf_results[[s]]$ccf_lag0
  cat(sprintf("    %s: r = %.3f  ->  z = %.3f\n", s, r_s, fisher_z(r_s)))
}
cat(sprintf("    Mean z = %.3f  ->  r_avg (Fisher) = %.3f\n",
            overall_avg$mean_z_FA[1], overall_avg$mean_ccf_FA[1]))
cat(sprintf("    Raw arithmetic mean r = %.3f  (bias vs Fisher: %.4f)\n\n",
            mean(c(ccf_results[["Site1"]]$ccf_lag0,
                   ccf_results[["Site2"]]$ccf_lag0,
                   ccf_results[["Site3"]]$ccf_lag0)),
            mean(c(ccf_results[["Site1"]]$ccf_lag0,
                   ccf_results[["Site2"]]$ccf_lag0,
                   ccf_results[["Site3"]]$ccf_lag0)) - overall_avg$mean_ccf_FA[1]))

cat(sprintf("  %-5s  %-12s  %-12s  %-12s  %s\n",
            "Lag", "Mean z(FA)", "r_avg(FA)", "r_avg(AF)", "n_sites"))
for (i in seq_len(nrow(overall_avg))) {
  cat(sprintf("  %-5d  %+.3f        %+.3f        %+.3f        %d\n",
              overall_avg$lag[i],
              overall_avg$mean_z_FA[i],
              overall_avg$mean_ccf_FA[i],
              overall_avg$mean_ccf_AF[i],
              overall_avg$n_sites[i]))
}
cat("  Note: reported as descriptive metric only — not interpretable as primary\n")
cat("  validation result given inconclusive Site 3 (n=8, CCF lag-0 ≈ 0).\n\n")

# ------------------------------------------------------------------------------
# SECTION 4C: PUBLICATION FIGURES — bars (mean) + points (per-site values)
# ------------------------------------------------------------------------------

cat("--- 4C: Aggregated level CCF figures ---\n\n")

# Helper to build the CCF bar+dot figure for one aggregated level
plot_avg_ccf <- function(avg_df, site_lag_df, sites_included,
                         level_label, file_label) {
  
  # Per-site points: filter to relevant sites and lags
  max_lag <- max(avg_df$lag)
  site_pts <- site_lag_df %>%
    dplyr::filter(Site %in% sites_included, lag <= max_lag) %>%
    dplyr::select(Site, lag, ccf_FA)
  
  # Mean lag-0 for subtitle
  lag0_mean <- avg_df$mean_ccf_FA[avg_df$lag == 0]
  
  p <- ggplot() +
    # Mean CCF bars
    geom_col(
      data  = avg_df,
      aes(x = lag, y = mean_ccf_FA),
      fill      = "gray60",
      color     = "gray30",
      width     = 0.5,
      linewidth = 0.3,
      alpha     = 0.85
    ) +
    # Per-site individual points
    geom_point(
      data = site_pts,
      aes(x = lag, y = ccf_FA, color = Site),
      size   = 3,
      shape  = 19,
      stroke = 0.5
    ) +
    scale_color_manual(values = custom_colors, name = "Site") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_x_continuous(breaks = avg_df$lag) +
    ylim(-1, 1) +
    labs(
      title    = sprintf("Cross-Correlation: FAIRD \u2192 AMMOD (%s)", level_label),
      subtitle = sprintf(
        paste0("Bars = mean CCF across sites | Points = individual site values\n",
               "Mean lag-0 CCF = %.3f | No formal CI (averaged series)"),
        lag0_mean
      ),
      x = "Lag (days)",
      y = "Cross-correlation (ACF)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, color = "gray40"),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "gray80", fill = NA,
                                      linewidth = 0.5)
    )
  
  print(p)
  
  ggsave(sprintf("outputs/H2_CrossCorr_Avg_%s.png", file_label),
         plot = p, width = 8, height = 5, dpi = 300, bg = "white")
  ggsave(sprintf("outputs/H2_CrossCorr_Avg_%s.pdf", file_label),
         plot = p, width = 8, height = 5, device = cairo_pdf)
  
  cat(sprintf("   [OK] Saved: outputs/H2_CrossCorr_Avg_%s.png/.pdf\n", file_label))
  
  invisible(p)
}

# Intra-Habitat Maize
plot_avg_ccf(
  avg_df         = maize_avg,
  site_lag_df    = site_lag_long,
  sites_included = c("Site1", "Site2"),
  level_label    = "Intra-Habitat Level (Maize)",
  file_label     = "Maize"
)

# Overall
plot_avg_ccf(
  avg_df         = overall_avg,
  site_lag_df    = site_lag_long,
  sites_included = c("Site1", "Site2", "Site3"),
  level_label    = "Overall Level",
  file_label     = "Overall"
)

cat("\n[OK] All aggregated level CCF figures saved\n\n")

# Store for downstream summary
agg_ccf_results <- list(
  "Intra-Habitat Level (Maize)" = list(
    mean_ccf_lag0       = maize_avg$mean_ccf_FA[1],   # Fisher's z back-transformed
    mean_ccf_lag0_raw   = mean(c(ccf_results[["Site1"]]$ccf_lag0,
                                 ccf_results[["Site2"]]$ccf_lag0)),
    avg_table           = maize_avg,
    sites               = c("Site1", "Site2"),
    method              = "Fisher's z-transformation"
  ),
  "Overall Level (All sites)" = list(
    mean_ccf_lag0       = overall_avg$mean_ccf_FA[1],
    mean_ccf_lag0_raw   = mean(c(ccf_results[["Site1"]]$ccf_lag0,
                                 ccf_results[["Site2"]]$ccf_lag0,
                                 ccf_results[["Site3"]]$ccf_lag0)),
    avg_table           = overall_avg,
    sites               = c("Site1", "Site2", "Site3"),
    method              = "Fisher's z-transformation (descriptive only)"
  )
)


# ==============================================================================
# SECTION 4B: SUMMED BIOMASS CCF — INTRA-HABITAT LEVEL (MAIZE)
# ==============================================================================
# Alternative approach: sum daily FAIRD inferred biomass and AMMOD gravimetric biomass across
# Sites 1 and 2 for each day, yielding a single genuine paired time series
# (n=12 days). This avoids both stacking artifacts and the absence of formal
# CI in the averaging approach.
#
# Interpretation: represents total daily inferred biomass of the Maize habitat
# rather than an average or pooled series. CCF CI is formally valid (single pair).
#
# Note: Site3 is excluded. An Overall Level equivalent would reduce to n=8
# (overlap period limited by Site3), which is too small to be meaningful.
# ==============================================================================

cat("\n=== SECTION 4B: SUMMED BIOMASS CCF (Maize, Sites 1+2) ===\n\n")
cat("Method: daily FAIRD inferred biomass and AMMOD gravimetric biomass summed across Sites 1+2.\n")
cat("Result: single genuine time series (n=12 days); formal CI valid.\n\n")

# ------------------------------------------------------------------------------
# 4B-1: BUILD SUMMED SERIES
# ------------------------------------------------------------------------------

ts_maize_summed <- ts_data %>%
  dplyr::filter(Site %in% c("Site1", "Site2")) %>%
  dplyr::group_by(Date) %>%
  dplyr::summarise(
    faird_value = sum(faird_value, na.rm = TRUE),
    ammod_value = sum(ammod_value, na.rm = TRUE),
    n_sites     = dplyr::n(),
    .groups     = "drop"
  ) %>%
  dplyr::filter(n_sites == 2) %>%   # retain only days with both sites present
  dplyr::arrange(Date)

cat(sprintf("Summed series: n = %d days (days with both sites present)\n",
            nrow(ts_maize_summed)))
cat(sprintf("Date range: %s to %s\n\n",
            min(ts_maize_summed$Date), max(ts_maize_summed$Date)))

# Per-day check
cat("Daily summed values (first 6 rows):\n")
print(head(ts_maize_summed %>% dplyr::select(Date, faird_value, ammod_value)), n = 6)
cat("\n")

# ------------------------------------------------------------------------------
# 4B-2: CCF ON SUMMED SERIES (reuses analyze_crosscorrelation())
# ------------------------------------------------------------------------------

cat("--- CCF on summed Maize series ---\n")

ccf_maize_summed <- analyze_crosscorrelation(ts_maize_summed, "Maize_Summed")

if (!is.null(ccf_maize_summed)) {
  
  cat(sprintf("   n = %d days\n", ccf_maize_summed$n))
  cat(sprintf("   Confidence bound: +/-%.3f\n", ccf_maize_summed$ci_bound))
  cat(sprintf("   >>> LAG 0 (same-day synchrony): CCF = %.3f %s\n",
              ccf_maize_summed$ccf_lag0,
              ifelse(abs(ccf_maize_summed$ccf_lag0) > ccf_maize_summed$ci_bound,
                     "* SIGNIFICANT", "(ns)")))
  
  if (length(ccf_maize_summed$sig_lags_faird_to_ammod) > 0) {
    cat(sprintf("   Significant lags beyond lag-0 (FAIRD & AMMOD): %s\n",
                paste(ccf_maize_summed$sig_lags_faird_to_ammod, collapse = ", ")))
  } else {
    cat("   No significant lags beyond lag-0 (FAIRD & AMMOD)\n")
  }
  
  lt  <- ccf_maize_summed$lag_table
  ci  <- ccf_maize_summed$ci_bound
  cat(sprintf("\n   Detailed CCF by lag (n=%d, CI = +/-%.3f):\n",
              ccf_maize_summed$n, ci))
  cat(sprintf("     %-5s  %-16s  %-16s\n", "Lag", "CCF(FAIRD->AMMOD)",
              "CCF(AMMOD->FAIRD)"))
  for (i in seq_len(nrow(lt))) {
    cat(sprintf("     %-5d  %+.3f %-7s    %+.3f %-7s\n",
                lt$lag[i],
                lt$ccf_FA[i], ifelse(lt$sig_FA[i], "[*]", ""),
                lt$ccf_AF[i], ifelse(lt$sig_AF[i], "[*]", "")))
  }
  cat("   [*] = |CCF| exceeds CI bound (significant at alpha=0.05)\n\n")
}

# ------------------------------------------------------------------------------
# 4B-3: TIME SERIES PLOT (summed series)
# ------------------------------------------------------------------------------

cat("--- Generating time series plot (summed series) ---\n")

ts_summed_plot_data <- ts_maize_summed %>%
  tidyr::pivot_longer(
    cols      = c(faird_value, ammod_value),
    names_to  = "Device",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Device = dplyr::case_when(
      Device == "faird_value" ~ "FAIRD (biomass mg, summed)",
      Device == "ammod_value" ~ "AMMOD (biomass mg, summed)"
    ),
    Device = factor(Device,
                    levels = c("FAIRD (biomass mg, summed)",
                               "AMMOD (biomass mg, summed)"))
  )

p_summed_ts <- ggplot(ts_summed_plot_data,
                      aes(x = Date, y = Value, color = Device)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~Device, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_color_manual(
    values = c(
      "FAIRD (biomass mg, summed)"    = unname(custom_colors["FAIRD"]),
      "AMMOD (biomass mg, summed)"   = unname(custom_colors["AMMOD"])
    )
  ) +
  labs(
    title    = "A) Time Series: Summed FAIRD + AMMOD (Maize, Sites 1+2)",
    subtitle = sprintf("n = %d days | Daily values summed across Sites 1 and 2",
                       nrow(ts_maize_summed)),
    x = "Date", y = "Value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13)
  )

print(p_summed_ts)

# ------------------------------------------------------------------------------
# 4B-4: COMBINED PUBLICATION FIGURE (time series + CCF bar)
# ------------------------------------------------------------------------------

if (!is.null(ccf_maize_summed)) {
  
  lt       <- ccf_maize_summed$lag_table
  ci       <- ccf_maize_summed$ci_bound
  lag0_val <- ccf_maize_summed$ccf_lag0
  lag0_sig <- abs(lag0_val) > ci
  
  p_summed_ccf <- ggplot(lt, aes(x = lag, y = ccf_FA)) +
    geom_col(
      aes(fill = sig_FA),
      width     = 0.6,
      color     = "gray30",
      linewidth = 0.3
    ) +
    scale_fill_manual(
      values = c("TRUE"  = unname(custom_colors["FAIRD"]),
                 "FALSE" = "gray70"),
      labels = c("TRUE" = "Significant", "FALSE" = "ns"),
      name   = NULL
    ) +
    geom_hline(yintercept =  ci, linetype = "dashed",
               color = "dodgerblue", linewidth = 0.8) +
    geom_hline(yintercept = -ci, linetype = "dashed",
               color = "dodgerblue", linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_x_continuous(breaks = lt$lag) +
    ylim(-1, 1) +
    labs(
      title    = "B) Cross-Correlation: FAIRD \u2192 AMMOD (Maize summed, Sites 1+2)",
      subtitle = sprintf(
        "Lag-0 CCF = %.3f%s | n = %d days | CI = \u00b1%.3f\n(summed series; formal CI valid)",
        lag0_val,
        ifelse(lag0_sig, " *", " ns"),
        ccf_maize_summed$n, ci
      ),
      x = "Lag (days)",
      y = "Cross-correlation (ACF)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, color = "gray40"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "gray80", fill = NA,
                                      linewidth = 0.5)
    )
  
  combined_summed <- p_summed_ts + p_summed_ccf +
    patchwork::plot_layout(ncol = 2, widths = c(1.4, 1))
  
  print(combined_summed)
  
  ggsave("outputs/H2_CrossCorr_Summed_Maize.png",
         plot = combined_summed, width = 14, height = 6,
         dpi = 300, bg = "white")
  ggsave("outputs/H2_CrossCorr_Summed_Maize.pdf",
         plot = combined_summed, width = 14, height = 6,
         device = cairo_pdf)
  
  cat("   [OK] Saved: outputs/H2_CrossCorr_Summed_Maize.png/.pdf\n\n")
}

# ------------------------------------------------------------------------------
# 4B-5: COMPARISON TABLE — three approaches for Intra-Habitat Level (Maize)
# ------------------------------------------------------------------------------

cat("--- Comparison: three CCF approaches for Intra-Habitat Level (Maize) ---\n\n")

cat(sprintf("  %-35s  %-6s  %-10s  %-10s  %-20s\n",
            "Approach", "n", "lag-0 CCF", "lag+1 CCF", "Notes"))
cat(sprintf("  %-35s  %-6s  %-10s  %-10s  %-20s\n",
            "-----------------------------------", "------",
            "----------", "----------", "--------------------"))

# Site-level (already computed)
for (s in c("Site1", "Site2")) {
  lt    <- ccf_results[[s]]$lag_table
  lag1  <- if (nrow(lt) >= 2) sprintf("%+.3f", lt$ccf_FA[2]) else "—"
  sig0  <- ifelse(abs(ccf_results[[s]]$ccf_lag0) > ccf_results[[s]]$ci_bound, "*", "ns")
  cat(sprintf("  %-35s  %-6d  %+.3f %-3s   %-10s  CI valid (pairwise)\n",
              paste0("Site Level (", s, ")"),
              ccf_results[[s]]$n,
              ccf_results[[s]]$ccf_lag0, sig0, lag1))
}

# Averaged approach
lag1_avg <- if (nrow(maize_avg) >= 2) sprintf("%+.3f", maize_avg$mean_ccf_FA[2]) else "—"
cat(sprintf("  %-35s  %-6s  %+.3f %-3s   %-10s  No formal CI\n",
            "Intra-Habitat Avg (Sites 1+2)",
            "2x12",
            maize_avg$mean_ccf_FA[1], "—", lag1_avg))

# Summed approach
if (!is.null(ccf_maize_summed)) {
  lt_s   <- ccf_maize_summed$lag_table
  lag1_s <- if (nrow(lt_s) >= 2) sprintf("%+.3f", lt_s$ccf_FA[2]) else "—"
  sig0_s <- ifelse(abs(ccf_maize_summed$ccf_lag0) > ccf_maize_summed$ci_bound, "*", "ns")
  cat(sprintf("  %-35s  %-6d  %+.3f %-3s   %-10s  CI valid (genuine series)\n",
              "Intra-Habitat Summed (Sites 1+2)",
              ccf_maize_summed$n,
              ccf_maize_summed$ccf_lag0, sig0_s, lag1_s))
}

cat("\n[OK] Section 4B complete\n\n")


# ==============================================================================
# SECTION 5: SUMMARY TABLE
# ==============================================================================

cat("\n=== SECTION 4: CROSS-CORRELATION SUMMARY ===\n\n")

if(length(ccf_results) > 0) {
  
  summary_table <- tibble::tibble(
    Site = names(ccf_results)
  ) %>%
    dplyr::mutate(
      n_days = purrr::map_dbl(Site, ~ccf_results[[.x]]$n),
      ci_bound = purrr::map_dbl(Site, ~ccf_results[[.x]]$ci_bound),
      n_sig_lags_FA = purrr::map_dbl(Site, ~length(ccf_results[[.x]]$sig_lags_faird_to_ammod)),
      n_sig_lags_AF = purrr::map_dbl(Site, ~length(ccf_results[[.x]]$sig_lags_ammod_to_faird)),
      max_ccf_FA = purrr::map_dbl(Site, ~ccf_results[[.x]]$max_ccf_fa),
      max_ccf_AF = purrr::map_dbl(Site, ~ccf_results[[.x]]$max_ccf_af),
      ccf_lag0 = purrr::map_dbl(Site, ~ccf_results[[.x]]$ccf_lag0)  # Now using direct field
    ) %>%
    dplyr::arrange(Site)
  
  cat("Cross-correlation summary by site:\n\n")
  print(summary_table, n = Inf)
  cat("\n")
  
  # Additional interpretation table
  cat("Interpretation key:\n")
  cat("  • n_days: Number of days analyzed\n")
  cat("  • ci_bound: ±95% confidence bound for significance\n")
  cat("  • n_sig_lags_FA/AF: Number of significant cross-correlation lags (excluding lag 0)\n")
  cat("  • max_ccf_FA/AF: Maximum absolute cross-correlation value (excluding lag 0)\n")
  cat("  • ccf_lag0: Same-day synchrony (most important for validation)\n\n")
  
  # Save summary
  write_csv(summary_table, "outputs/H2_CrossCorr_Summary.csv")
  cat("[OK] Summary saved: outputs/H2_CrossCorr_Summary.csv\n\n")
  
  # Interpretation
  cat("=== INTERPRETATION NOTES ===\n\n")
  cat("1. ACF panels show autocorrelation structure of each series\n")
  cat("   - Lag 0 always = 1.0 (perfect correlation with itself)\n")
  cat("   - Significant lags beyond 0 indicate temporal dependencies\n")
  cat("   - Values beyond dashed blue lines are significant (α=0.05)\n\n")
  
  cat("2. Cross-correlation panels show relationships between FAIRD and AMMOD\n")
  cat("   - 'FAIRD & AMMOD': How FAIRD values correlate with AMMOD at each lag\n")
  cat("   - 'AMMOD & FAIRD': How AMMOD values correlate with FAIRD at each lag\n")
  cat("   - Lag 0: Same-day correlation (most important for validation)\n")
  cat("   - Lag 1+: Temporal lead-lag relationships\n\n")
  
  cat("3. Confidence bounds (dashed blue lines):\n")
  cat("   - Calculated as ±1.96/√n (95% confidence)\n")
  cat("   - Values outside bounds indicate significant correlation\n\n")
  
  cat("4. Reference: Crawley (2013) The R Book, Chapter 24.8\n\n")
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n================================================================\n")
cat("|                  ANALYSIS COMPLETE                           |\n")
cat("================================================================\n\n")

cat("FILES GENERATED:\n\n")

cat("Time Series Plots:\n")
for(site in sites) {
  if(!is.null(ts_plots[[site]])) {
    cat(sprintf("  [OK] outputs/H2_CrossCorr_TimeSeries_%s.png/pdf\n", site))
  }
}
cat("  [OK] outputs/H2_CrossCorr_TimeSeries_Combined.png/pdf\n")
if (exists("agg_file_labels")) {
  for (lbl in agg_file_labels) {
    cat(sprintf("  [OK] outputs/H2_CrossCorr_TimeSeries_%s.png/pdf\n", lbl))
  }
}
cat("\n")

cat("Aggregated Level CCF Figures:\n")
cat("  [OK] outputs/H2_CrossCorr_Avg_Maize.png/.pdf\n")
cat("  [OK] outputs/H2_CrossCorr_Avg_Overall.png/.pdf\n")
cat("  [OK] outputs/H2_CrossCorr_Summed_Maize.png/.pdf\n\n")

cat("Cross-Correlation Plots:\n")
for(site in sites) {
  if(!is.null(ccf_results[[site]])) {
    cat(sprintf("  [OK] outputs/H2_CrossCorr_ACF_%s.png\n", site))
    cat(sprintf("  [OK] outputs/H2_CrossCorr_Combined_%s.png/pdf\n", site))
  }
}
if (exists("agg_file_labels")) {
  for (lbl in agg_file_labels) {
    cat(sprintf("  [OK] outputs/H2_CrossCorr_Combined_%s.png/pdf\n", lbl))
  }
}
cat("\n")

cat("\nSummary Table:\n")
cat("  [OK] outputs/H2_CrossCorr_Summary.csv\n\n")

# ==============================================================================
# EXECUTIVE SUMMARY - QUICK RESULTS INTERPRETATION
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("|              EXECUTIVE SUMMARY - KEY FINDINGS                |\n")
cat("================================================================\n\n")

if(length(ccf_results) > 0) {
  
  # Extract lag-0 CCF for all sites (most important for validation)
  lag0_results <- tibble::tibble(
    Site = names(ccf_results)
  ) %>%
    dplyr::mutate(
      n_days = purrr::map_dbl(Site, ~ccf_results[[.x]]$n),
      ci_bound = purrr::map_dbl(Site, ~ccf_results[[.x]]$ci_bound),
      ccf_lag0 = purrr::map_dbl(Site, ~ccf_results[[.x]]$ccf_lag0),  # Direct field access
      significant = abs(ccf_lag0) > ci_bound
    ) %>%
    dplyr::arrange(desc(ccf_lag0))
  
  cat("1. SAME-DAY TEMPORAL SYNCHRONY (Lag 0 - Most Important for Validation):\n")
  cat("   ────────────────────────────────────────────────────────────\n")
  for(i in 1:nrow(lag0_results)) {
    row <- lag0_results[i, ]
    sig_marker <- ifelse(row$significant, "★ SIGNIFICANT", "(ns)")
    interpretation <- dplyr::case_when(
      abs(row$ccf_lag0) >= 0.7 ~ "STRONG",
      abs(row$ccf_lag0) >= 0.5 ~ "MODERATE", 
      abs(row$ccf_lag0) >= 0.3 ~ "WEAK",
      TRUE ~ "VERY WEAK"
    )
    
    cat(sprintf("   %s: CCF = %+.3f %s [%s synchrony]\n",
                row$Site, row$ccf_lag0, sig_marker, interpretation))
  }
  cat("\n")
  
  # Overall assessment
  n_significant <- sum(lag0_results$significant)
  mean_ccf <- mean(abs(lag0_results$ccf_lag0))
  
  cat("2. OVERALL VALIDATION ASSESSMENT:\n")
  cat("   ────────────────────────────────────────────────────────────\n")
  cat(sprintf("   Sites with significant lag-0 CCF: %d/%d (%.1f%%)\n",
              n_significant, nrow(lag0_results),
              100 * n_significant / nrow(lag0_results)))
  cat(sprintf("   Mean |CCF| at lag 0: %.3f\n", mean_ccf))
  
  overall_quality <- dplyr::case_when(
    n_significant >= 2 & mean_ccf >= 0.5 ~ "EXCELLENT",
    n_significant >= 2 & mean_ccf >= 0.3 ~ "GOOD",
    n_significant >= 1 & mean_ccf >= 0.3 ~ "MODERATE",
    TRUE ~ "WEAK"
  )
  
  cat(sprintf("   >>> VALIDATION QUALITY: %s\n", overall_quality))
  cat("\n")
  
  # Interpretation guidance
  cat("3. INTERPRETATION:\n")
  cat("   ────────────────────────────────────────────────────────────\n")
  
  if(overall_quality == "EXCELLENT") {
    cat("   ✓ Strong temporal synchrony detected between FAIRD and AMMOD\n")
    cat("   ✓ Devices track same-day fluctuations consistently\n")
    cat("   ✓ Validates FAIRD's capability for dynamic monitoring\n")
    cat("   ✓ SUCCESSFUL temporal validation\n")
  } else if(overall_quality == "GOOD") {
    cat("   ✓ Moderate to strong synchrony in most sites\n")
    cat("   ✓ Devices show coordinated responses to daily dynamics\n")
    cat("   ⚠ Some site-specific variation present\n")
    cat("   ✓ SUCCESSFUL with site considerations\n")
  } else if(overall_quality == "MODERATE") {
    cat("   ⚠ Limited temporal synchrony detected\n")
    cat("   ⚠ Synchrony may be site-dependent\n")
    cat("   ⚠ Methodological differences should be discussed\n")
    cat("   → Different temporal scales or community components should be considered?\n")
  } else {
    cat("   ⚠ Weak temporal synchrony across sites\n")
    cat("   → Recommendation: Honest reporting of limitations\n")
    cat("   → Discuss: Devices may track different temporal patterns\n")
    cat("   → Consider: Small sample size (n=%d) or high variability\n", 
        min(lag0_results$n_days))
  }
  cat("\n")
  
  # Lead-lag analysis
  has_leadlag <- FALSE
  leadlag_summary <- list()
  
  for(site in names(ccf_results)) {
    n_fa <- length(ccf_results[[site]]$sig_lags_faird_to_ammod)
    n_af <- length(ccf_results[[site]]$sig_lags_ammod_to_faird)
    
    # Check if there are significant lags beyond lag 0
    sig_lags_fa <- ccf_results[[site]]$sig_lags_faird_to_ammod
    sig_lags_af <- ccf_results[[site]]$sig_lags_ammod_to_faird
    
    # Remove lag 0 if present (we already reported it)
    sig_lags_fa_nonzero <- sig_lags_fa   # lag 0 already excluded upstream via [-1]
    sig_lags_af_nonzero <- sig_lags_af
    
    if(length(sig_lags_fa_nonzero) > 0 || length(sig_lags_af_nonzero) > 0) {
      has_leadlag <- TRUE
      leadlag_summary[[site]] <- list(
        fa = sig_lags_fa_nonzero,
        af = sig_lags_af_nonzero
      )
    }
  }
  
  cat("4. LEAD-LAG RELATIONSHIPS:\n")
  cat("   ────────────────────────────────────────────────────────────\n")
  
  if(has_leadlag) {
    cat("   Significant temporal lags detected (beyond same-day):\n\n")
    
    for(site in names(leadlag_summary)) {
      info <- leadlag_summary[[site]]
      
      if(length(info$fa) > 0) {
        actual_lags <- info$fa   # index already equals actual lag number
        cat(sprintf("   %s: FAIRD leads by %s day(s)\n",
                    site, paste(actual_lags, collapse = ", ")))
        cat("      → FAIRD detects changes before they appear in AMMOD\n")
      }
      
      if(length(info$af) > 0) {
        actual_lags <- info$af   # index already equals actual lag number
        cat(sprintf("   %s: AMMOD leads by %s day(s)\n",
                    site, paste(actual_lags, collapse = ", ")))
        cat("      → AMMOD accumulates biomass before FAIRD detection\n")
      }
      cat("\n")
    }
    
    cat("   Interpretation:\n")
    cat("   • Lead-lag patterns reflect methodological differences\n")
    cat("   • Continuous (FAIRD) vs. 24h integrated (AMMOD) sampling\n")
    cat("   • NOT a validation problem - additional information!\n")
  } else {
    cat("   No significant lead-lag relationships detected\n")
    cat("   → Devices show simultaneous responses (lag 0 only)\n")
    cat("   → Indicates similar temporal resolution\n")
  }
  cat("\n")
  
  # Autocorrelation structure
  cat("5. AUTOCORRELATION STRUCTURE:\n")
  cat("   ────────────────────────────────────────────────────────────\n")
  
  auto_consistent <- TRUE
  for(site in names(ccf_results)) {
    acf_obj <- ccf_results[[site]]$acf_obj  # Changed from pacf_obj to acf_obj
    
    # Check lag 1 autocorrelation for both devices
    # ACF includes lag 0 at index 1, so lag 1 is at index 2
    faird_auto_lag1 <- acf_obj$acf[2, 1, 1]  # index 2 = lag 1
    ammod_auto_lag1 <- acf_obj$acf[2, 2, 2]  # index 2 = lag 1
    ci <- ccf_results[[site]]$ci_bound
    
    faird_sig <- abs(faird_auto_lag1) > ci
    ammod_sig <- abs(ammod_auto_lag1) > ci
    
    cat(sprintf("   %s:\n", site))
    cat(sprintf("      FAIRD lag-1 autocorr: %+.3f %s\n",
                faird_auto_lag1, ifelse(faird_sig, "★", "")))
    cat(sprintf("      AMMOD lag-1 autocorr: %+.3f %s\n",
                ammod_auto_lag1, ifelse(ammod_sig, "★", "")))
    
    # Check if similar
    if((faird_sig && ammod_sig && sign(faird_auto_lag1) == sign(ammod_auto_lag1)) ||
       (!faird_sig && !ammod_sig)) {
      cat("      → Similar temporal structure ✓\n")
    } else {
      cat("      → Different temporal dynamics\n")
      auto_consistent <- FALSE
    }
    cat("\n")
  }
  
  if(auto_consistent) {
    cat("   ✓ Devices show consistent autocorrelation patterns\n")
    cat("   → Both respond similarly to environmental persistence\n")
  } else {
    cat("   ⚠ Devices show different autocorrelation patterns\n")
    cat("   → May reflect different sensitivities or time scales\n")
  }
  cat("\n")
  
  # Sample size considerations
  cat("6. STATISTICAL CONSIDERATIONS:\n")
  cat("   ────────────────────────────────────────────────────────────\n")
  cat(sprintf("   Sample size: n = %d-%d days per site\n",
              min(lag0_results$n_days), max(lag0_results$n_days)))
  cat(sprintf("   Confidence bounds: ±%.3f to ±%.3f\n",
              min(lag0_results$ci_bound), max(lag0_results$ci_bound)))
  cat("\n")
  
  # Detailed CCF values by lag per site
  cat("   Detailed CCF by lag:\n\n")
  for(site in names(ccf_results)) {
    lt <- ccf_results[[site]]$lag_table
    ci <- ccf_results[[site]]$ci_bound
    cat(sprintf("   %s (n=%d, CI = +/-%.3f):\n",
                site, ccf_results[[site]]$n, ci))
    cat(sprintf("     %-5s  %-16s  %-16s\n",
                "Lag", "CCF(FAIRD->AMMOD)", "CCF(AMMOD->FAIRD)"))
    for(i in 1:nrow(lt)) {
      cat(sprintf("     %-5d  %+.3f %-7s    %+.3f %-7s\n",
                  lt$lag[i],
                  lt$ccf_FA[i], ifelse(lt$sig_FA[i], "[*]", ""),
                  lt$ccf_AF[i], ifelse(lt$sig_AF[i], "[*]", "")))
    }
    cat("\n")
  }
  cat("   [*] = |CCF| exceeds CI bound (significant at alpha=0.05)\n\n")
  
  if(min(lag0_results$n_days) < 10) {
    cat("   ⚠ Small sample size limits statistical power\n")
    cat("   → Focus on effect sizes (CCF magnitude) over p-values\n")
    cat("   → Consistent patterns across sites = stronger evidence\n")
  } else {
    cat("   ✓ Adequate sample size for cross-correlation analysis\n")
  }
  cat("\n")
  
} else {
  cat("   No cross-correlation results to display\n")
  cat("   → Check data availability and filtering criteria\n\n")
}

cat("========================================\n")
cat("CROSS-CORRELATION ANALYSIS COMPLETE\n")
cat("========================================\n\n")

# ---- Aggregated levels quick summary (for manuscript §3.3.4) ----
if (exists("agg_ccf_results") && length(agg_ccf_results) > 0) {
  cat("--- AGGREGATED LEVEL CCF SUMMARY (Fisher's z-transformation) ---\n\n")
  cat("  PRIMARY RESULT — Intra-Habitat Level (Maize, summed series):\n")
  if (!is.null(ccf_maize_summed)) {
    cat(sprintf("    CCF lag-0 = %.3f %s (n=12, CI = +/-%.3f; formal CI valid)\n",
                ccf_maize_summed$ccf_lag0,
                ifelse(abs(ccf_maize_summed$ccf_lag0) > ccf_maize_summed$ci_bound,
                       "*", "ns"),
                ccf_maize_summed$ci_bound))
  }
  cat("\n  DESCRIPTIVE — Overall Level (Fisher's z average, Sites 1+2+3):\n")
  if (exists("overall_avg")) {
    cat(sprintf("    r_avg (Fisher's z) lag-0 = %.3f (no formal CI)\n",
                overall_avg$mean_ccf_FA[1]))
    cat(sprintf("    Raw arithmetic mean lag-0 = %.3f (shown for comparison)\n",
                mean(c(ccf_results[["Site1"]]$ccf_lag0,
                       ccf_results[["Site2"]]$ccf_lag0,
                       ccf_results[["Site3"]]$ccf_lag0))))
    cat("    Note: not interpretable as primary validation result\n")
    cat("    (Site 3: n=8, CCF lag-0 ≈ 0 — inconclusive power)\n")
  }
  cat("\n")
}

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")