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
# Note: Using ABUNDANCE for FAIRD (counts) and LIVE_MASS for AMMOD (biomass)
# These are the actual measured quantities, not transformed values

ts_data <- etraps_daily_overlap %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::group_by(Date, Site) %>%
  dplyr::summarise(
    faird_value = sum(abundance, na.rm = TRUE),
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
        Device == "faird_value" ~ "FAIRD (abundance)",
        Device == "ammod_value" ~ "AMMOD (biomass mg)"
      ),
      Device = factor(Device, levels = c("FAIRD (abundance)", "AMMOD (biomass mg)"))
    )
  
  # Define colors specifically for this plot (avoid warning)
  plot_colors <- c(
    "FAIRD (abundance)" = unname(custom_colors["FAIRD"]),
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
# SECTION 4: SUMMARY TABLE
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
cat("  [OK] outputs/H2_CrossCorr_TimeSeries_Combined.png/pdf\n\n")

cat("Cross-Correlation Plots:\n")
for(site in sites) {
  if(!is.null(ccf_results[[site]])) {
    cat(sprintf("  [OK] outputs/H2_CrossCorr_ACF_%s.png\n", site))
  }
}

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

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")