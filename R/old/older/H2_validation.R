# ==============================================================================
# H2: FAIRD-AMMOD VALIDATION ANALYSIS - WORKING VERSION
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
  model <- lm(FAIRD_biomass_log ~ AMMOD_biomass_log, data = site_data_log)
  shapiro_test <- shapiro.test(residuals(model))
  bp_test <- car::ncvTest(model)
  
  # 4. Run Pearson on LOG DATA
  pearson_test <- cor.test(site_data_log$FAIRD_biomass_log, 
                           site_data_log$AMMOD_biomass_log, 
                           method = "pearson")
  
  # 5. Run Spearman on RAW DATA (non-parametric verification)
  spearman_test <- cor.test(site_data_raw$FAIRD_biomass, 
                            site_data_raw$AMMOD_biomass, 
                            method = "spearman", exact = FALSE)
  
  list(
    site = site_name, 
    n = nrow(site_data_raw), 
    data_raw = site_data_raw, 
    data_log = site_data_log, # Return log data for plotting
    model = model,
    r_squared = summary(model)$r.squared,
    pearson_r = pearson_test$estimate, 
    pearson_p = pearson_test$p.value,
    spearman_rho = spearman_test$estimate, 
    spearman_p = spearman_test$p.value,
    shapiro_p = shapiro_test$p.value, 
    bp_p = bp_test$p,
    residuals_normal = shapiro_test$p.value > 0.05
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
    cat(sprintf("n=%d, RÃ‚Â²(log)=%.3f, p=%.4f\n", result$n, result$r_squared, result$pearson_p))
    cat(sprintf("   Diagnostics: Shapiro p=%.3f, Breusch-Pagan p=%.3f\n\n", 
                result$shapiro_p, result$bp_p))
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

cat("\n=== SECTION 3: TABLE ===\n\n")

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
write_csv(results_table, "outputs/Table2_H2_Validation_Log.csv")
cat("\n[OK] Table saved\n\n")

# ==============================================================================
# SECTION 4: CORRELATION PLOTS (ON LOG SCALE)
# ==============================================================================

cat("\n=== SECTION 4: CORRELATION PLOTS (Log Scale) ===\n\n")

plot_list <- list()

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    result <- site_results[[site]]
    site_color <- custom_colors[site]
    
    # +++ GRAFICAR USANDO data_log +++
    p <- ggplot(result$data_log, aes(x = AMMOD_biomass_log, y = FAIRD_biomass_log)) +
      geom_point(size = 3, alpha = 0.7, color = "gray30") + 
      geom_smooth(
        method = "lm", 
        se = TRUE, 
        color = "gray30",       
        fill = site_color,      
        alpha = 0.3, 
        linewidth = 1
      ) +
      labs(
        title = sprintf("%s (%s)", site, results_table$Environment[results_table$Site == site]),
        # +++ ETIQUETAS LOG +++
        x = "log₁₀(AMMOD Live Mass + 1)", 
        y = "log₁₀(FAIRD Biomass + 1)"
      ) +
      # +++ USE LOG MODEL RESULTS +++
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("R² = %.3f\np = %.3f", result$r_squared, result$pearson_p),
               hjust = -0.1, vjust = 1.2, size = 3.5, fontface = "italic") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    
    plot_list[[site]] <- p
  }
}

combined <- wrap_plots(plot_list, ncol = length(plot_list))
print(combined)
ggsave("outputs/Figure3_H2_Validation_Log.png", combined, width = 5 * length(plot_list), height = 5, dpi = 300, bg = "white")
ggsave("outputs/Figure3_H2_Validation_Log.pdf", combined, width = 5 * length(plot_list), height = 5, device = cairo_pdf)

cat("[OK] Correlation plots saved\n\n")

# ==============================================================================
# SECTION 5: TEMPORAL DYNAMICS (WITH PSEUDO-LOG TRANSFORMATION)
# ==============================================================================

cat("\n=== SECTION 5: TEMPORAL DYNAMICS ===\n\n")

temporal_plots <- list()
site_number <- 0

for(site in sites) {
  if(!is.null(site_results[[site]])) {
    site_number <- site_number + 1
    
    # Prepare data (use data_raw)
    site_temporal <- site_results[[site]]$data_raw %>%
      dplyr::select(Date, FAIRD = FAIRD_biomass, AMMOD = AMMOD_biomass) %>%
      tidyr::pivot_longer(cols = c(FAIRD, AMMOD), names_to = "Method", values_to = "Biomass")
    
    faird_color <- custom_colors["FAIRD"]
    ammod_color <- custom_colors["AMMOD"]
    
    cat(sprintf("%s: FAIRD=%s, AMMOD=%s\n", site, faird_color, ammod_color))
    
    x_label <- if(site_number == length(sites)) "Date" else NULL
    y_label <- if(site_number == round(length(sites) / 2)) "Biomass (mg)" else NULL
    
    # Create plot with colors
    p <- ggplot(site_temporal, aes(x = Date, y = Biomass, color = Method, group = Method)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5, alpha = 0.8) +
      scale_color_manual(
        values = custom_colors,
        name = "Device Type"
      ) +
      
      # +++ CORRECTION: USE PSEUDO-LOG TO HANDLE ZEROS +++
      scale_y_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 1),
        labels = scales::comma,
        breaks = c(0, 10, 100, 1000, 5000) # Breaks manuales
      ) +
      labs(
        title = site,
        x = x_label,
        y = y_label
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "none",
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

ggsave("Figure_S2_Temporal_Log.png", combined_temporal, width = 10, height = 4 * length(temporal_plots), dpi = 300, bg = "white")
ggsave("Figure_S2_Temporal_Log.pdf", combined_temporal, width = 10, height = 4 * length(temporal_plots), device = cairo_pdf)

cat("\n[OK] Temporal dynamics saved (with pseudo-log)\n\n")

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
    subtitle = "Letters indicate significance groups (Dunn's test with Bonferroni correction, ÃŽÂ± = 0.05).\nGroups sharing letters are not significantly different.",
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

cat("=== CREATING BAR PLOT (Figure S3 style) ===\n\n")

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
    subtitle = "Error bars: Ã‚Â± 1 SD. Letters indicate significance groups (Dunn's test, Bonferroni, ÃŽÂ± = 0.05).",
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
ggsave("outputs/Figure_S3_Devices_WithLetters.png", p_bars_letters, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("outputs/Figure_S3_Devices_WithLetters.pdf", p_bars_letters, width = 10, height = 6, device = "pdf")
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

# ==============================================================================
# KRUSKAL-WALLIS TEST (Three-way comparison)
# ==============================================================================

cat("=== KRUSKAL-WALLIS TEST (Three Devices) ===\n\n")

kw_devices <- kruskal.test(daily_value ~ Device_type, data = comparison_data)

# Calculate effect size
n_devices <- nrow(comparison_data)
k_devices <- length(unique(comparison_data$Device_type))
effect_size_devices <- calculate_epsilon_squared(kw_devices, n_devices, k_devices)

cat("DEVICE COMPARISON (Sites with AMMOD available):\n")
cat(sprintf("  • Kruskal-Wallis: χ² = %.3f, df = %d, p = %.4f %s\n",
            kw_devices$statistic,
            kw_devices$parameter,
            kw_devices$p.value,
            ifelse(kw_devices$p.value < 0.05, "***", "(NS)")))
cat(sprintf("  • Effect size (ε²): %.3f (%s)\n",
            effect_size_devices$epsilon_squared,
            effect_size_devices$interpretation))
cat(sprintf("    - Variance explained: %.1f%%\n",
            effect_size_devices$variance_explained_pct))
cat("\n")

if(kw_devices$p.value < 0.05) {
  cat("Significant difference detected. Proceeding with post-hoc tests...\n\n")
  
  # ==============================================================================
  # DUNN'S POST-HOC TEST with Bonferroni correction
  # ==============================================================================
  
  cat("=== DUNN'S POST-HOC TEST (Bonferroni correction) ===\n\n")
  
  # Run Dunn's test
  dunn_result <- dunn.test::dunn.test(
    x = comparison_data$daily_value,
    g = comparison_data$Device_type,
    method = "bonferroni",
    kw = FALSE, # Don't repeat KW test
    label = TRUE,
    wrap = FALSE,
    table = FALSE,
    list = TRUE,
    altp = FALSE  # Changed to FALSE to get numeric p-values
  )
  
  # Extract and format results
  comparisons <- dunn_result$comparisons
  z_values <- as.numeric(dunn_result$Z)
  p_values <- as.numeric(dunn_result$P.adjusted)
  
  cat("  • Dunn's post-hoc (Bonferroni):\n")
  for(i in seq_along(comparisons)) {
    sig_marker <- ifelse(p_values[i] < 0.001, "***",
                         ifelse(p_values[i] < 0.01, "**",
                                ifelse(p_values[i] < 0.05, "*", "(NS)")))
    
    cat(sprintf("    - %s: Z = %.3f, p = %.4f %s\n",
                comparisons[i],
                z_values[i],
                p_values[i],
                sig_marker))
  }
  cat("\n")
  
  # Create summary table
  dunn_table <- tibble::tibble(
    Comparison = comparisons,
    Z_statistic = round(z_values, 3),
    p_adjusted = round(p_values, 4),
    Significant = p_values < 0.05
  )
  
  cat("=== PAIRWISE COMPARISONS TABLE ===\n\n")
  print(dunn_table, n = Inf)
  cat("\n")
  
  # Save results
  write_csv(dunn_table, "outputs/H2_Section6B_Dunn_PostHoc.csv")
  cat("[OK] Dunn's test results saved to: outputs/H2_Section6B_Dunn_PostHoc.csv\n\n")
  
} else {
  cat("No significant difference found. Post-hoc tests not necessary.\n\n")
}

# Save KW results
kw_summary <- tibble::tibble(
  Test = "Kruskal-Wallis",
  Statistic = kw_devices$statistic,
  df = kw_devices$parameter,
  p_value = kw_devices$p.value,
  epsilon_squared = effect_size_devices$epsilon_squared,
  effect_interpretation = effect_size_devices$interpretation,
  variance_explained_pct = effect_size_devices$variance_explained_pct,
  Significant = kw_devices$p.value < 0.05
)

write_csv(kw_summary, "outputs/H2_Section6B_KruskalWallis.csv")
cat("[OK] Kruskal-Wallis results saved to: outputs/H2_Section6B_KruskalWallis.csv\n\n")

# ==============================================================================

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
ggsave("outputs/Figure_S4_Weekly.png", plot_weekly, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("outputs/Figure_S4_Weekly.pdf", plot_weekly, width = 12, height = 8, device = cairo_pdf)
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
    subtitle = sprintf("Letters indicate significance within each environment (Mann-Whitney U, ÃŽÂ± = 0.05).\nMaize: p = %.3f; Meadow: p = %.3f",
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
    subtitle = "Letters indicate significance (Mann-Whitney U, ÃŽÂ± = 0.05)",
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
# KRUSKAL-WALLIS TEST (Global habitat comparison)
# ==============================================================================

cat("=== KRUSKAL-WALLIS TEST (Habitat × Device) ===\n\n")

kw_habitat <- kruskal.test(value ~ Habitat_Device, data = habitat_data_full)

# Calculate effect size
n_habitat <- nrow(habitat_data_full)
k_habitat <- length(unique(habitat_data_full$Habitat_Device))
effect_size_habitat <- calculate_epsilon_squared(kw_habitat, n_habitat, k_habitat)

cat("HABITAT COMPARISON (All devices and environments):\n")
cat(sprintf("  • Kruskal-Wallis: χ² = %.3f, df = %d, p = %.4f %s\n",
            kw_habitat$statistic,
            kw_habitat$parameter,
            kw_habitat$p.value,
            ifelse(kw_habitat$p.value < 0.05, "***", "(NS)")))
cat(sprintf("  • Effect size (ε²): %.3f (%s)\n",
            effect_size_habitat$epsilon_squared,
            effect_size_habitat$interpretation))
cat(sprintf("    - Variance explained: %.1f%%\n",
            effect_size_habitat$variance_explained_pct))
cat("\n")

if(kw_habitat$p.value < 0.05) {
  cat("Significant difference detected. Proceeding with post-hoc tests...\n\n")
  
  cat("=== DUNN'S POST-HOC TEST (Bonferroni correction) ===\n\n")
  
  # Run Dunn's test
  dunn_habitat <- dunn.test::dunn.test(
    x = habitat_data_full$value,
    g = habitat_data_full$Habitat_Device,
    method = "bonferroni",
    kw = FALSE,
    label = TRUE,
    wrap = FALSE,
    table = FALSE,
    list = TRUE,
    altp = FALSE  # Changed to FALSE to get numeric p-values
  )
  
  # Extract results
  comparisons_hab <- dunn_habitat$comparisons
  z_values_hab <- as.numeric(dunn_habitat$Z)
  p_values_hab <- as.numeric(dunn_habitat$P.adjusted)
  
  cat("  • Significant pairwise comparisons (p < 0.05):\n")
  
  significant_count <- 0
  for(i in seq_along(comparisons_hab)) {
    if(p_values_hab[i] < 0.05) {
      significant_count <- significant_count + 1
      sig_marker <- ifelse(p_values_hab[i] < 0.001, "***",
                           ifelse(p_values_hab[i] < 0.01, "**", "*"))
      
      cat(sprintf("    - %s: Z = %.3f, p = %.4f %s\n",
                  comparisons_hab[i],
                  z_values_hab[i],
                  p_values_hab[i],
                  sig_marker))
    }
  }
  
  if(significant_count == 0) {
    cat("    (None - all p-values > 0.05 after Bonferroni correction)\n")
  }
  cat("\n")
  
  # Full table (including non-significant)
  dunn_habitat_table <- tibble::tibble(
    Comparison = comparisons_hab,
    Z_statistic = round(z_values_hab, 3),
    p_adjusted = round(p_values_hab, 4),
    Significant = p_values_hab < 0.05
  ) %>%
    dplyr::arrange(p_adjusted)
  
  cat("=== FULL PAIRWISE COMPARISONS TABLE (sorted by p-value) ===\n\n")
  print(dunn_habitat_table, n = Inf)
  cat("\n")
  
  # Save results
  write_csv(dunn_habitat_table, "outputs/H2_Section8B_Dunn_Habitat.csv")
  cat("[OK] Dunn's test results saved to: outputs/H2_Section8B_Dunn_Habitat.csv\n\n")
  
} else {
  cat("No significant difference found. Post-hoc tests not necessary.\n\n")
}

# Save KW results
kw_habitat_summary <- tibble::tibble(
  Test = "Kruskal-Wallis (Habitat)",
  Statistic = kw_habitat$statistic,
  df = kw_habitat$parameter,
  p_value = kw_habitat$p.value,
  epsilon_squared = effect_size_habitat$epsilon_squared,
  effect_interpretation = effect_size_habitat$interpretation,
  variance_explained_pct = effect_size_habitat$variance_explained_pct,
  Significant = kw_habitat$p.value < 0.05
)

write_csv(kw_habitat_summary, "outputs/H2_Section8B_KruskalWallis_Habitat.csv")
cat("[OK] Kruskal-Wallis results saved to: outputs/H2_Section8B_KruskalWallis_Habitat.csv\n\n")

# ==============================================================================
cat("========================================\n\n")

cat("\n========================================\n")
cat("H2 COMPLETE\n")
cat("========================================\n\n")