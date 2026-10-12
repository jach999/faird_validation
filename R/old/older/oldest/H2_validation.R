# ==============================================================================
# H2: FAIRD-AMMOD VALIDATION ANALYSIS - WORKING VERSION
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(car)
library(broom)
library(RColorBrewer)

cat("\n========================================\n")
cat("H2: FAIRD-AMMOD VALIDATION ANALYSIS\n")
cat("========================================\n\n")

# ==============================================================================
# CUSTOM COLORS
# ==============================================================================

# Load custom color palette
source("custom_colors.R")

cat("✅ Colors: FAIRD=", custom_colors["FAIRD"], 
    ", ID=", custom_colors["ID"], 
    ", AMMOD=", custom_colors["AMMOD"], "\n\n")

# Load data
load("data_cleaned.RData")

# ==============================================================================
# SECTION 0: DATA PREPARATION
# ==============================================================================

cat("=== LOADING DATA ===\n\n")

etraps_daily <- data_list$etraps_daily
ammod_biomass <- data_list$ammod_biomass_for_stats
etraps_daily_overlap <- data_list$etraps_daily_overlap
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
# SECTION 1-3: ANALYSES
# ==============================================================================

analyze_faird_ammod_site <- function(paired_data, site_name) {
  site_data <- paired_data %>%
    dplyr::filter(Site == site_name) %>%
    dplyr::select(Date, FAIRD_biomass = faird_biomass_mg, AMMOD_biomass = ammod_live_mass)
  
  if(nrow(site_data) < 3) return(NULL)
  
  model <- lm(FAIRD_biomass ~ AMMOD_biomass, data = site_data)
  shapiro_test <- shapiro.test(residuals(model))
  bp_test <- car::ncvTest(model)
  pearson_test <- cor.test(site_data$FAIRD_biomass, site_data$AMMOD_biomass, method = "pearson")
  spearman_test <- cor.test(site_data$FAIRD_biomass, site_data$AMMOD_biomass, method = "spearman", exact = FALSE)
  
  list(
    site = site_name, n = nrow(site_data), data = site_data, model = model,
    r_squared = summary(model)$r.squared,
    pearson_r = pearson_test$estimate, pearson_p = pearson_test$p.value,
    spearman_rho = spearman_test$estimate, spearman_p = spearman_test$p.value,
    shapiro_p = shapiro_test$p.value, bp_p = bp_test$p,
    residuals_normal = shapiro_test$p.value > 0.05
  )
}

analyze_temporal_sync <- function(site_data) {
  if(nrow(site_data) < 3) return(NULL)
  
  site_data <- site_data %>%
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

cat("\n=== SECTION 1: CORRELATIONS ===\n\n")

sites <- c("Site1", "Site2", "Site3")
site_results <- list()

for(site in sites) {
  cat(sprintf("--- %s ---\n", site))
  result <- analyze_faird_ammod_site(faird_ammod_paired, site)
  
  if(!is.null(result)) {
    site_results[[site]] <- result
    cat(sprintf("n=%d, R²=%.3f, p=%.4f\n\n", result$n, result$r_squared, result$pearson_p))
  }
}

cat("\n=== SECTION 2: TEMPORAL SYNC ===\n\n")

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    sync_result <- analyze_temporal_sync(site_results[[site]]$data)
    if(!is.null(sync_result)) {
      site_results[[site]]$sync <- sync_result
      cat(sprintf("%s: %.1f%% concordance\n", site, sync_result$concordance_pct))
    }
  }
}

cat("\n=== SECTION 3: TABLE ===\n\n")

results_table <- tibble(
  Site = sites,
  Environment = c("Maize-Edge", "Maize-Core", "Meadow"),
  n_days = sapply(site_results, function(x) if(!is.null(x)) x$n else NA),
  R2 = sapply(site_results, function(x) if(!is.null(x)) x$r_squared else NA),
  Pearson_p = sapply(site_results, function(x) if(!is.null(x)) x$pearson_p else NA)
)

print(results_table)
write_csv(results_table, "Table2_H2_Validation.csv")
cat("\n✅ Table saved\n\n")

# ==============================================================================
# SECTION 4: CORRELATION PLOTS
# ==============================================================================

cat("\n=== SECTION 4: CORRELATION PLOTS ===\n\n")

plot_list <- list()

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    result <- site_results[[site]]
    site_color <- custom_colors[site]
    
    p <- ggplot(result$data, aes(x = AMMOD_biomass, y = FAIRD_biomass)) +
      geom_point(size = 3, alpha = 0.7, color = "gray30") +  # Dark gray points
      geom_smooth(
        method = "lm", 
        se = TRUE, 
        color = "gray30",        # Dark gray line
        fill = site_color,       # Site-colored confidence band
        alpha = 0.3, 
        linewidth = 1
      ) +
      labs(
        title = sprintf("%s (%s)", site, results_table$Environment[results_table$Site == site]),
        x = "AMMOD Live Mass (mg)", 
        y = "FAIRD Biomass (mg)"
      ) +
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("R² = %.3f\np = %.3f", result$r_squared, result$pearson_p),
               hjust = -0.1, vjust = 1.2, size = 3.5, fontface = "italic") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    
    plot_list[[site]] <- p
  }
}

combined <- wrap_plots(plot_list, ncol = 3)
print(combined)
ggsave("Figure3_H2_Validation.png", combined, width = 14, height = 5, dpi = 300, bg = "white")
ggsave("Figure3_H2_Validation.pdf", combined, width = 14, height = 5, device = cairo_pdf)
cat("✅ Correlation plots saved\n\n")

# ==============================================================================
# SECTION 5: TEMPORAL DYNAMICS - WITH COLORS AND SHARED LABELS
# ==============================================================================

cat("\n=== SECTION 5: TEMPORAL DYNAMICS ===\n\n")

temporal_plots <- list()
site_number <- 0

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    site_number <- site_number + 1
    
    # Prepare data for this site
    site_temporal <- site_results[[site]]$data %>%
      dplyr::select(Date, FAIRD = FAIRD_biomass, AMMOD = AMMOD_biomass) %>%
      tidyr::pivot_longer(cols = c(FAIRD, AMMOD), names_to = "Method", values_to = "Biomass")
    
    # Get colors - FAIRD always blue, AMMOD always darkorange
    faird_color <- custom_colors["FAIRD"]    # blue
    ammod_color <- custom_colors["AMMOD"]    # darkorange
    
    cat(sprintf("%s: FAIRD=%s, AMMOD=%s\n", site, faird_color, ammod_color))
    
    # Determine labels based on position
    x_label <- if(site_number == 3) "Date" else NULL
    y_label <- if(site_number == 2) "Biomass (mg) [log₁₀ scale]" else NULL
    
    # Create plot with colors
    p <- ggplot(site_temporal, aes(x = Date, y = Biomass, color = Method, group = Method)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5, alpha = 0.8) +
      scale_color_manual(
        values = custom_colors,
        name = "Device Type"  # Change legend title
      ) +
      scale_y_log10(
        labels = scales::comma,
        breaks = scales::trans_breaks("log10", function(x) 10^x)
      ) +
      annotation_logticks(sides = "l") +
      labs(
        title = site,
        x = x_label,
        y = y_label
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "none",  # Remove individual legends
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 11, face = "bold"),
        axis.title.y = element_text(size = 11, face = "bold")
      )
    
    temporal_plots[[site]] <- p
  }
}

# Combine plots with shared legend at bottom
combined_temporal <- wrap_plots(temporal_plots, ncol = 1, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10))

print(combined_temporal)

# Save figures
ggsave("Figure_S2_Temporal.png", combined_temporal, width = 10, height = 12, dpi = 300, bg = "white")
ggsave("Figure_S2_Temporal.pdf", combined_temporal, width = 10, height = 12, device = cairo_pdf)

cat("\n✅ Temporal dynamics saved: FAIRD (blue) + AMMOD (darkorange)\n\n")

# ==============================================================================
# SECTION 6: DEVICE COMPARISON - EXACTLY like analysis_general.R
# ==============================================================================

cat("\n=== SECTION 6: DEVICE COMPARISON ===\n\n")

# Prepare data exactly like analysis_general.R does
device_summary <- etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::group_by(Device_type, Site) %>%
  dplyr::summarise(
    w = mean(abundance, na.rm = TRUE),
    sd = sd(abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(device_summary)

# Plot EXACTLY like analysis_general.R
# Key: fill = Device_type (the column), and pass custom_colors directly
plot_device <- ggplot(device_summary, aes(x = Site, y = w, fill = Device_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = w - sd, ymax = w + sd),
                position = position_dodge(width = 0.8), width = 0.3) +
  # CRITICAL: Pass entire custom_colors vector, not subset
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Average Daily Abundance by Device Type and Site",
    x = "Site",
    y = "Mean Daily Abundance (individuals)",
    fill = "Device Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

print(plot_device)
ggsave("Figure_S3_Devices.png", plot_device, width = 10, height = 6, dpi = 300, bg = "white")
cat("\n✅ Device comparison saved\n\n")

# ==============================================================================
# SECTION 7: WEEKLY PATTERNS - EXACTLY like analysis_general.R
# ==============================================================================

cat("\n=== SECTION 7: WEEKLY PATTERNS ===\n\n")

etraps_weekly <- etraps_daily %>%
  dplyr::mutate(Week = lubridate::isoweek(Date)) %>%
  dplyr::group_by(Week, Device_type, Site) %>%
  dplyr::summarise(w = mean(abundance, na.rm = TRUE), n = n(), .groups = "drop")

print(head(etraps_weekly, 10))

# Plot EXACTLY like analysis_general.R
plot_weekly <- ggplot(etraps_weekly, aes(x = factor(Week), y = w, fill = Device_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Site, ncol = 2, scales = "free_y") +
  # CRITICAL: Pass entire custom_colors vector
  scale_fill_manual(values = custom_colors) +
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
ggsave("Figure_S4_Weekly.png", plot_weekly, width = 12, height = 8, dpi = 300, bg = "white")
cat("\n✅ Weekly saved\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("H2 COMPLETE - All figures with colors!\n")
cat("========================================\n\n")