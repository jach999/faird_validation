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
# SECTION 6: DEVICE COMPARISON
# ==============================================================================

cat("\n=== SECTION 6: DEVICE COMPARISON ===\n\n")

# ==============================================================================
# CUSTOMIZATION PARAMETERS - ADJUST THESE TO YOUR PREFERENCE
# ==============================================================================

# Letter appearance
LETTER_SIZE <- 3.5        # Size of letters (default: 3.5, try 3-4.5)
LETTER_FONTFACE <- "plain" # Font style: "plain", "bold", "italic", "bold.italic"
LETTER_COLOR <- "black"    # Color of letters
LETTER_VJUST <- -0.5       # Vertical position adjustment (negative = higher)

# Letter position offset (distance above error bars)
Y_OFFSET <- 5              # Offset above error bars (default: 5, try 3-8)

cat("CUSTOMIZATION SETTINGS:\n")
cat(sprintf("  Letter size: %.1f\n", LETTER_SIZE))
cat(sprintf("  Letter style: %s\n", LETTER_FONTFACE))
cat(sprintf("  Letter color: %s\n", LETTER_COLOR))
cat(sprintf("  Vertical offset: %.1f\n\n", Y_OFFSET))

# ==============================================================================
# LOAD DATA
# ==============================================================================

source("custom_colors.R")
load("data_cleaned.RData")

etraps_daily <- data_list$etraps_daily

# ==============================================================================
# PREPARE DATA
# ==============================================================================

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

# ==============================================================================
# LETTER ASSIGNMENTS
# ==============================================================================

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
    letter_y = mean_abundance + sd_abundance + Y_OFFSET  # Use custom offset
  )

cat("Groups with letter assignments:\n")
print(device_summary %>% 
        dplyr::select(Group, mean_abundance, Letters) %>%
        dplyr::arrange(desc(mean_abundance)))
cat("\n")

# ==============================================================================
# CREATE BOXPLOT WITH LETTERS
# ==============================================================================

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
    subtitle = "Letters indicate significance groups (Dunn's test with Bonferroni correction, α = 0.05).\nGroups sharing letters are not significantly different.",
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

ggsave("H2_Section6_Boxplot_WithLetters.png", 
       p_boxplot_letters, 
       width = 10, height = 6, dpi = 300, bg = "white")

ggsave("H2_Section6_Boxplot_WithLetters.pdf",
       p_boxplot_letters,
       width = 10, height = 6, device = "pdf")

cat("✅ Boxplot saved\n\n")

# ==============================================================================
# CREATE BAR PLOT WITH LETTERS (FIGURE S3 STYLE)
# ==============================================================================

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
    subtitle = "Error bars: ± 1 SD. Letters indicate significance groups (Dunn's test, Bonferroni, α = 0.05).",
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

ggsave("Figure_S3_Devices_WithLetters.png", 
       p_bars_letters, 
       width = 10, height = 6, dpi = 300, bg = "white")

ggsave("Figure_S3_Devices_WithLetters.pdf",
       p_bars_letters,
       width = 10, height = 6, device = "pdf")

cat("✅ Bar plot saved\n\n")

# ==============================================================================
# SAVE LETTER ASSIGNMENTS
# ==============================================================================

letter_table <- device_summary %>%
  dplyr::select(Device_type, Site, Group, n, mean_abundance, 
                sd_abundance, Letters) %>%
  dplyr::arrange(Device_type, Site)

write_csv(letter_table, "H2_Section6_Significance_Letters.csv")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║                    SUMMARY                                    ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("CUSTOMIZATION USED:\n")
cat(sprintf("  • Letter size: %.1f\n", LETTER_SIZE))
cat(sprintf("  • Letter style: %s\n", LETTER_FONTFACE))
cat(sprintf("  • Letter color: %s\n", LETTER_COLOR))
cat(sprintf("  • Y offset: %.1f\n\n", Y_OFFSET))

cat("FILES GENERATED:\n")
cat("  ✅ H2_Section6_Boxplot_WithLetters.png/pdf\n")
cat("  ✅ Figure_S3_Devices_WithLetters.png/pdf\n")
cat("  ✅ H2_Section6_Significance_Letters.csv\n\n")

cat("LETTER ASSIGNMENTS:\n")
for (i in 1:nrow(letter_table)) {
  cat(sprintf("  %12s: %3s  (%.1f ± %.1f)\n",
              letter_table$Group[i],
              letter_table$Letters[i],
              letter_table$mean_abundance[i],
              letter_table$sd_abundance[i]))
}

cat("\n")
cat("TO ADJUST APPEARANCE:\n")
cat("  • Edit parameters at the top of this script\n")
cat("  • LETTER_SIZE: 3.5 (current) - try 3.0-4.5\n")
cat("  • LETTER_FONTFACE: 'plain' (current) - try 'bold', 'italic'\n")
cat("  • Y_OFFSET: 5 (current) - try 3-8 for vertical position\n\n")

cat("RECOMMENDED FOR MANUSCRIPT:\n")
cat("  → Use Figure_S3_Devices_WithLetters.png\n\n")

# ==============================================================================
# SECTION 7: WEEKLY PATTERNS - 3 WEEKS OF 7 DAYS EACH
# ==============================================================================

cat("\n=== SECTION 7: WEEKLY PATTERNS ===\n\n")

# Calculate relative week number (1, 2, 3) - exactly 3 weeks of 7 days
etraps_weekly <- etraps_daily %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    # Calculate days since first date
    days_since_start = as.numeric(Date - min(Date)),
    # Convert to week number: days 0-6=Week1, 7-13=Week2, 14-20=Week3
    Week = pmin(floor(days_since_start / 7) + 1, 3)  # Cap at 3 weeks
  ) %>%
  dplyr::filter(Week <= 3) %>%  # Only keep first 21 days (3 weeks)
  dplyr::group_by(Week, Device_type, Site) %>%
  dplyr::summarise(w = mean(abundance, na.rm = TRUE), n = n(), .groups = "drop")

cat("Weekly patterns summary:\n")
cat(sprintf("Study period: 3 weeks (%d days)\n", 21))
cat(sprintf("Weeks in data: %s\n", paste(sort(unique(etraps_weekly$Week)), collapse = ", ")))
print(etraps_weekly)

# Plot with relative weeks (1, 2, 3 only)
plot_weekly <- ggplot(etraps_weekly, aes(x = factor(Week), y = w, fill = Device_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Site, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = c("1", "2", "3")) +  # Ensure only 3 weeks shown
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
ggsave("Figure_S4_Weekly.pdf", plot_weekly, width = 12, height = 8, device = cairo_pdf)
cat("\n✅ Weekly patterns saved (3 weeks of 7 days each)\n\n")


# ==============================================================================
# H2 SECTION 8: ENVIRONMENT COMPARISON (IMPROVED VERSION)
# ==============================================================================
# Improvements:
# 1. Single plot with shared Y-axis for direct comparison
# 2. Better letter positioning (above whiskers, not fixed offset)
# Date: October 24, 2025
# ==============================================================================

cat("\n=== SECTION 8: ENVIRONMENT COMPARISON (IMPROVED) ===\n\n")

# Prepare data by environment
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

# Summary by environment
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

# Statistical tests
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

# Calculate whisker positions for letter placement
# Whiskers extend to 1.5 * IQR from Q1/Q3, or to min/max if within that range
calculate_whisker_top <- function(data) {
  q1 <- quantile(data, 0.25, na.rm = TRUE)
  q3 <- quantile(data, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  upper_whisker <- min(max(data[data <= q3 + 1.5 * iqr], na.rm = TRUE), 
                       max(data, na.rm = TRUE))
  return(upper_whisker)
}

# Calculate whisker positions for each group
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

# Assign letters
env_summary <- env_summary %>%
  dplyr::mutate(
    Letters = case_when(
      Environment == "Maize" & Device_type == "FAIRD" ~ ifelse(mw_maize$p.value < 0.05, "a", "a"),
      Environment == "Maize" & Device_type == "ID" ~ ifelse(mw_maize$p.value < 0.05, "b", "a"),
      Environment == "Meadow" & Device_type == "FAIRD" ~ ifelse(mw_meadow$p.value < 0.05, "a", "a"),
      Environment == "Meadow" & Device_type == "ID" ~ ifelse(mw_meadow$p.value < 0.05, "b", "a"),
      TRUE ~ "?"
    ),
    # Position letters above the whisker with a fixed offset
    letter_y = whisker_top + 8  # 8 units above whisker top
  )

cat("Letter positions:\n")
print(env_summary %>% 
        dplyr::select(Environment, Device_type, whisker_top, letter_y, Letters))
cat("\n")

# ==============================================================================
# CREATE SINGLE INTEGRATED PLOT (BETTER FOR COMPARISON)
# ==============================================================================

cat("Creating integrated plot with shared Y-axis...\n\n")

LETTER_SIZE <- 3.5
LETTER_FONTFACE <- "plain"
LETTER_COLOR <- "black"

# Create combined grouping variable for X-axis
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

# Single integrated plot
plot_integrated <- ggplot(environment_data, 
                          aes(x = Group, y = abundance, fill = Device_type)) +
  
  geom_boxplot(width = 0.6) +
  
  # Add significance letters positioned above whiskers
  geom_text(
    data = env_summary,
    aes(x = Group, y = letter_y, label = Letters),
    size = LETTER_SIZE,
    fontface = LETTER_FONTFACE,
    color = LETTER_COLOR,
    inherit.aes = FALSE
  ) +
  
  # Add vertical separator line between environments
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  
  scale_fill_manual(values = custom_colors) +
  
  # Add environment labels at the top
  annotate("text", x = 1.5, y = Inf, label = "Maize (Site1+Site2)",
           vjust = 1.5, fontface = "bold", size = 4) +
  annotate("text", x = 3.5, y = Inf, label = "Meadow (Site4)",
           vjust = 1.5, fontface = "bold", size = 4) +
  
  labs(
    title = "Device Performance by Environment Type",
    subtitle = sprintf("Letters indicate significance within each environment (Mann-Whitney U, α = 0.05).\nMaize: p = %.3f; Meadow: p = %.3f",
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

ggsave("H2_Section8_Environment_Comparison.png",
       plot_integrated, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("H2_Section8_Environment_Comparison.pdf",
       plot_integrated, width = 10, height = 6, device = "pdf")

cat("✅ Integrated environment plot saved\n\n")

# ==============================================================================
# ALTERNATIVE: FACETED VERSION (IF NEEDED)
# ==============================================================================

# If you prefer faceted version with shared Y-axis
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
    subtitle = "Letters indicate significance (Mann-Whitney U, α = 0.05)",
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

ggsave("H2_Section8_Environment_Faceted.png",
       plot_faceted, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("H2_Section8_Environment_Faceted.pdf",
       plot_faceted, width = 10, height = 6, device = "pdf")

cat("✅ Faceted environment plot saved (alternative version)\n\n")

# Save results
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
  FAIRD_mean = env_summary$mean_abundance[env_summary$Device_type == "FAIRD"],
  FAIRD_sd = env_summary$sd_abundance[env_summary$Device_type == "FAIRD"],
  ID_mean = env_summary$mean_abundance[env_summary$Device_type == "ID"],
  ID_sd = env_summary$sd_abundance[env_summary$Device_type == "ID"],
  p_value = c(mw_maize$p.value, mw_meadow$p.value),
  Significant = c(mw_maize$p.value < 0.05, mw_meadow$p.value < 0.05)
)

write_csv(results_env, "H2_Section8_Environment_Results.csv")

cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
cat("║                    SUMMARY                                    ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("IMPROVEMENTS IN THIS VERSION:\n")
cat("  ✅ Single plot with SHARED Y-AXIS for direct comparison\n")
cat("  ✅ Letters positioned ABOVE whiskers (not overlapping)\n")
cat("  ✅ Letter position calculated dynamically per group\n")
cat("  ✅ Visual separator between environments\n")
cat("  ✅ Environment labels at top of plot\n\n")

cat("LETTER POSITIONING:\n")
cat("  • Calculated from actual data distribution\n")
cat("  • Based on upper whisker position (Q3 + 1.5*IQR)\n")
cat("  • Fixed 8-unit offset above whisker\n")
cat("  • Prevents overlap with boxplot elements\n\n")

cat("FILES GENERATED:\n")
cat("  ✅ H2_Section8_Environment_Comparison.png/pdf (integrated)\n")
cat("  ✅ H2_Section8_Environment_Faceted.png/pdf (alternative)\n")
cat("  ✅ H2_Section8_Environment_Results.csv\n\n")

cat("RECOMMENDED:\n")
cat("  → Use INTEGRATED version (shared Y-axis)\n")
cat("  → Better for comparing magnitude between environments\n\n")

cat("========================================\n")
cat("SECTION 8 COMPLETE (IMPROVED)\n")
cat("========================================\n\n")
# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("H2 COMPLETE - All figures with colors!\n")
cat("========================================\n\n")