# ==============================================================================
# H3: TAXONOMIC BREADTH COMPARISON
# ==============================================================================
# Description: Complete taxonomic comparison between FAIRD, ID, and AMMOD
# Date: October 22, 2025
# ==============================================================================

# --- Setup ---
library(tidyverse)
library(vegan)      # For diversity indices
library(ggpubr)     # For plot annotations
library(patchwork)  # For combining plots

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

# Load data
cat("\n========================================\n")
cat("H3: TAXONOMIC BREADTH COMPARISON\n")
cat("========================================\n\n")

source("load_data.R")
load("data_cleaned.RData")

# Load custom colors
source("custom_colors.R")

cat("✅ Colors: FAIRD=", custom_colors["FAIRD"], 
    ", ID=", custom_colors["ID"], 
    ", AMMOD=", custom_colors["AMMOD"], "\n\n")

etraps_data <- data_list$etraps
ammod_taxonomy <- data_list$ammod_taxonomy

# ==============================================================================
# [NEW] OTU/SPECIES COUNT VERIFICATION
# ==============================================================================
cat("\n=== VERIFICACIÓN DE RIQUEZA TAXONÓMICA (OTUs/Especies) ===\n")

# 1. Definir fechas válidas (Batch 1 + Batch 2)
# Excluyendo cambio de botella (4 Sep) y fallo mecánico AMMOD4 (28 Aug - 3 Sep)
valid_dates_batch1 <- seq(as.Date("2023-08-23"), as.Date("2023-09-03"), by="day")
valid_dates_batch2 <- seq(as.Date("2023-09-05"), as.Date("2023-09-13"), by="day")
all_valid_dates <- c(valid_dates_batch1, valid_dates_batch2)

ammod4_failure_dates <- seq(as.Date("2023-08-28"), as.Date("2023-09-03"), by="day")

# 2. Filtrar el dataset
analyzed_taxa <- data_list$ammod_taxonomy %>%
  dplyr::filter(Date %in% all_valid_dates) %>%                       # Solo fechas válidas
  dplyr::filter(!(Device == "AMMOD4" & Date %in% ammod4_failure_dates)) # Excluir fallo AMMOD4

# 3. Contar Especies Únicas (Esto es tu "N OTUs")
# Usamos na.rm = TRUE por si hay alguna fila vacía
n_total_species <- length(unique(na.omit(analyzed_taxa$Species)))

cat("Total de Especies únicas (OTUs) en el periodo analizado:", n_total_species, "\n")
cat("============================================================\n\n")

# ==============================================================================
# SECTION 1: TAXONOMIC RICHNESS SUMMARY
# ==============================================================================

cat("\n=== SECTION 1: TAXONOMIC RICHNESS ===\n\n")

# Calculate richness for each device type
richness_summary <- bind_rows(
  
  # FAIRD richness
  etraps_data %>%
    dplyr::filter(Device_type == "FAIRD") %>%
    dplyr::summarise(
      Device_type = "FAIRD",
      N_individuals = n(),
      N_orders = n_distinct(Order, na.rm = TRUE),
      N_families = n_distinct(Family, na.rm = TRUE),
      N_genera = n_distinct(Genus, na.rm = TRUE)
    ),
  
  # ID richness
  etraps_data %>%
    dplyr::filter(Device_type == "ID") %>%
    dplyr::summarise(
      Device_type = "ID",
      N_individuals = n(),
      N_orders = n_distinct(Order, na.rm = TRUE),
      N_families = n_distinct(Family, na.rm = TRUE),
      N_genera = n_distinct(Genus, na.rm = TRUE)
    ),
  
  # AMMOD richness (from metabarcoding)
  ammod_taxonomy %>%
    dplyr::summarise(
      Device_type = "AMMOD",
      N_individuals = NA,  # Not applicable (metabarcoding)
      N_orders = n_distinct(Order, na.rm = TRUE),
      N_families = n_distinct(Family, na.rm = TRUE),
      N_genera = n_distinct(Genus, na.rm = TRUE)
    )
) %>%
  dplyr::mutate(
    # Calculate advantage of FAIRD over ID
    Genera_advantage_vs_ID = ifelse(Device_type == "FAIRD", 
                                    round((N_genera / 5 - 1) * 100, 0), 
                                    NA)
  )

print(richness_summary)

cat("\n# KEY FINDING:\n")
cat(sprintf("   FAIRD captures %d%% MORE genera than Insect Detect\n",
            richness_summary$Genera_advantage_vs_ID[1]))
cat(sprintf("   FAIRD: %d genera vs ID: %d genera\n\n",
            richness_summary$N_genera[1],
            richness_summary$N_genera[2]))

# ==============================================================================
# SECTION 2: TAXONOMIC COMPOSITION ANALYSIS
# ==============================================================================

cat("\n=== SECTION 2: TAXONOMIC COMPOSITION ===\n\n")

# Order-level composition for e-traps
composition_etraps <- etraps_data %>%
  dplyr::group_by(Device_type, Order) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  dplyr::group_by(Device_type) %>%
  dplyr::mutate(
    Total = sum(Count),
    Percentage = (Count / Total) * 100
  ) %>%
  dplyr::ungroup()

cat("--- E-traps Composition (by Order) ---\n")

# FAIRD composition
cat("\nFAIRD:\n")
composition_etraps %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::arrange(desc(Percentage)) %>%
  dplyr::select(Order, Count, Percentage) %>%
  print(n = 15)

# ID composition
cat("\nInsect Detect:\n")
composition_etraps %>%
  dplyr::filter(Device_type == "ID") %>%
  dplyr::arrange(desc(Percentage)) %>%
  dplyr::select(Order, Count, Percentage) %>%
  print(n = 15)

# Calculate ID specialization
id_diptera_pct <- composition_etraps %>%
  dplyr::filter(Device_type == "ID", Order == "Diptera") %>%
  dplyr::pull(Percentage)

cat(sprintf("\n# ID Specialization: %.1f%% Diptera\n", id_diptera_pct))
cat("   â†’ ID functions as Diptera specialist\n")
cat("   â†’ FAIRD functions as generalist sampler\n\n")

# ==============================================================================
# SECTION 2B: DIPTERA FAMILY COMPOSITION (FAIRD vs ID)
# ==============================================================================

cat("\n=== SECTION 2B: DIPTERA FAMILY COMPOSITION ===\n\n")
cat("Objective: Verify if ID is specialized in Syrphidae (hover flies)\n")
cat("           as claimed by the developer\n\n")

# Filter only Diptera records
diptera_data <- etraps_data %>%
  dplyr::filter(Order == "Diptera")

# Calculate family composition within Diptera for each device
diptera_families <- diptera_data %>%
  dplyr::group_by(Device_type, Family) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  dplyr::group_by(Device_type) %>%
  dplyr::mutate(
    Total = sum(Count),
    Percentage = (Count / Total) * 100
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Device_type, desc(Percentage))

# Display FAIRD Diptera families
cat("--- FAIRD Diptera Family Composition ---\n")
diptera_families %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::select(Family, Count, Percentage) %>%
  print(n = 20)

# Display ID Diptera families
cat("\n--- Insect Detect Diptera Family Composition ---\n")
diptera_families %>%
  dplyr::filter(Device_type == "ID") %>%
  dplyr::select(Family, Count, Percentage) %>%
  print(n = 20)

# Calculate Syrphidae specialization for ID
id_syrphidae_pct <- diptera_families %>%
  dplyr::filter(Device_type == "ID", Family == "Syrphidae") %>%
  dplyr::pull(Percentage)

faird_syrphidae_pct <- diptera_families %>%
  dplyr::filter(Device_type == "FAIRD", Family == "Syrphidae") %>%
  dplyr::pull(Percentage)

# Handle case where Syrphidae might not be present
if(length(id_syrphidae_pct) == 0) id_syrphidae_pct <- 0
if(length(faird_syrphidae_pct) == 0) faird_syrphidae_pct <- 0

cat("\n📊 SYRPHIDAE SPECIALIZATION ANALYSIS:\n")
cat(sprintf("   ID: %.1f%% of Diptera are Syrphidae\n", id_syrphidae_pct))
cat(sprintf("   FAIRD: %.1f%% of Diptera are Syrphidae\n", faird_syrphidae_pct))

if(id_syrphidae_pct > faird_syrphidae_pct) {
  cat(sprintf("   ✅ CONFIRMED: ID shows %.1fx higher Syrphidae proportion\n", 
              id_syrphidae_pct / faird_syrphidae_pct))
  cat("   → Developer claim VALIDATED\n")
} else if(id_syrphidae_pct > 0) {
  cat("   ⚠️  ID does capture Syrphidae but not predominantly\n")
  cat("   → Mixed evidence for developer claim\n")
} else {
  cat("   ❌ No Syrphidae detected in ID captures\n")
  cat("   → Developer claim not supported by data\n")
}

# Calculate family richness within Diptera
diptera_richness <- diptera_data %>%
  dplyr::group_by(Device_type) %>%
  dplyr::summarise(
    N_families = n_distinct(Family, na.rm = TRUE),
    N_genera = n_distinct(Genus, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Diptera Taxonomic Richness ---\n")
print(diptera_richness)

cat("\n💡 INTERPRETATION:\n")
cat("   If ID shows high Syrphidae proportion → Specialized design working\n")
cat("   If ID shows diverse families → Functions as general Diptera trap\n\n")

# ==============================================================================
# SECTION 2C: DIPTERA SUPERFAMILY COMPOSITION (FAIRD vs ID)
# ==============================================================================

cat("\n=== SECTION 2C: DIPTERA SUPERFAMILY COMPOSITION ===\n\n")
cat("Objective: Analyze at Superfamily level (more appropriate taxonomically)\n")
cat("           Syrphoidea encompasses all hover flies (syrphids)\n\n")

# Calculate superfamily composition within Diptera for each device
diptera_superfamilies <- diptera_data %>%
  dplyr::group_by(Device_type, Superfamily) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  dplyr::group_by(Device_type) %>%
  dplyr::mutate(
    Total = sum(Count),
    Percentage = (Count / Total) * 100
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Device_type, desc(Percentage))

# Display FAIRD Diptera superfamilies
cat("--- FAIRD Diptera Superfamily Composition ---\n")
diptera_superfamilies %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::select(Superfamily, Count, Percentage) %>%
  print(n = 20)

# Display ID Diptera superfamilies
cat("\n--- Insect Detect Diptera Superfamily Composition ---\n")
diptera_superfamilies %>%
  dplyr::filter(Device_type == "ID") %>%
  dplyr::select(Superfamily, Count, Percentage) %>%
  print(n = 20)

# Calculate Syrphoidea specialization for ID
id_syrphoidea_pct <- diptera_superfamilies %>%
  dplyr::filter(Device_type == "ID", Superfamily == "Syrphoidea") %>%
  dplyr::pull(Percentage)

faird_syrphoidea_pct <- diptera_superfamilies %>%
  dplyr::filter(Device_type == "FAIRD", Superfamily == "Syrphoidea") %>%
  dplyr::pull(Percentage)

# Handle case where Syrphoidea might not be present
if(length(id_syrphoidea_pct) == 0) id_syrphoidea_pct <- 0
if(length(faird_syrphoidea_pct) == 0) faird_syrphoidea_pct <- 0

cat("\n📊 SYRPHOIDEA SPECIALIZATION ANALYSIS:\n")
cat(sprintf("   ID: %.1f%% of Diptera are Syrphoidea\n", id_syrphoidea_pct))
cat(sprintf("   FAIRD: %.1f%% of Diptera are Syrphoidea\n", faird_syrphoidea_pct))

if(id_syrphoidea_pct > faird_syrphoidea_pct) {
  cat(sprintf("   ✅ CONFIRMED: ID shows %.1fx higher Syrphoidea proportion\n", 
              id_syrphoidea_pct / faird_syrphoidea_pct))
  cat("   → Developer claim VALIDATED at Superfamily level\n")
} else if(id_syrphoidea_pct > 0) {
  cat("   ⚠️  ID does capture Syrphoidea but not predominantly\n")
  cat("   → Mixed evidence for developer claim\n")
} else {
  cat("   ❌ No Syrphoidea detected in ID captures\n")
  cat("   → Developer claim not supported by data\n")
}

# Calculate superfamily richness within Diptera
diptera_superfamily_richness <- diptera_data %>%
  dplyr::group_by(Device_type) %>%
  dplyr::summarise(
    N_superfamilies = n_distinct(Superfamily, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Diptera Superfamily Richness ---\n")
print(diptera_superfamily_richness)

cat("\n💡 TAXONOMIC NOTE:\n")
cat("   Superfamily Syrphoidea includes Family Syrphidae (hover flies)\n")
cat("   This level provides better taxonomic resolution than Family alone\n")
cat("   Syrphoidea = all syrphids (hover flies, pollinating Diptera)\n\n")

# ==============================================================================
# SECTION 3: GENUS-LEVEL DIVERSITY
# ==============================================================================

cat("\n=== SECTION 3: GENUS-LEVEL DIVERSITY ===\n\n")

# Get genus lists for each device type
genera_faird <- etraps_data %>%
  dplyr::filter(Device_type == "FAIRD", !is.na(Genus)) %>%
  dplyr::distinct(Genus) %>%
  dplyr::pull(Genus)

genera_id <- etraps_data %>%
  dplyr::filter(Device_type == "ID", !is.na(Genus)) %>%
  dplyr::distinct(Genus) %>%
  dplyr::pull(Genus)

genera_ammod <- ammod_taxonomy %>%
  dplyr::filter(!is.na(Genus)) %>%
  dplyr::distinct(Genus) %>%
  dplyr::pull(Genus)

# Calculate overlaps
overlap_faird_id <- length(intersect(genera_faird, genera_id))
overlap_faird_ammod <- length(intersect(genera_faird, genera_ammod))
overlap_id_ammod <- length(intersect(genera_id, genera_ammod))

unique_faird <- length(setdiff(genera_faird, union(genera_id, genera_ammod)))
unique_id <- length(setdiff(genera_id, union(genera_faird, genera_ammod)))
unique_ammod <- length(setdiff(genera_ammod, union(genera_faird, genera_id)))

cat("Genus-level overlaps:\n")
cat(sprintf("  FAIRD âˆ© ID: %d genera\n", overlap_faird_id))
cat(sprintf("  FAIRD âˆ© AMMOD: %d genera\n", overlap_faird_ammod))
cat(sprintf("  ID âˆ© AMMOD: %d genera\n\n", overlap_id_ammod))

cat("Unique genera:\n")
cat(sprintf("  FAIRD only: %d genera\n", unique_faird))
cat(sprintf("  ID only: %d genera\n", unique_id))
cat(sprintf("  AMMOD only: %d genera\n\n", unique_ammod))

# ==============================================================================
# SECTION 4: CREATE RESULTS TABLE
# ==============================================================================

cat("\n=== SECTION 4: RESULTS TABLE (TABLE 3) ===\n\n")

results_table <- richness_summary %>%
  dplyr::select(Device_type, N_individuals, N_orders, N_families, N_genera) %>%
  dplyr::mutate(
    Comparison_vs_ID = case_when(
      Device_type == "FAIRD" ~ sprintf("+%d%%", 
                                       round((N_genera / 5 - 1) * 100, 0)),
      Device_type == "ID" ~ "Baseline",
      Device_type == "AMMOD" ~ "Reference"
    ),
    Notes = case_when(
      Device_type == "FAIRD" ~ "Semi-automated, generalist",
      Device_type == "ID" ~ ">95% Diptera specialist",
      Device_type == "AMMOD" ~ "Metabarcoding (qualitative)"
    )
  )

print(results_table)

# Save table
write_csv(results_table, "outputs/Table3_H3_Taxonomic_Results.csv")
cat("\nTable saved: outputs/Table3_H3_Taxonomic_Results.csv\n\n")

# ==============================================================================
# SECTION 5: VISUALIZATION - RICHNESS BAR PLOT
# ==============================================================================

cat("\n=== SECTION 5: RICHNESS BAR PLOT ===\n\n")

# Prepare data for plotting (E-traps only for fair comparison)
richness_plot_data <- richness_summary %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  tidyr::pivot_longer(
    cols = c(N_orders, N_families, N_genera),
    names_to = "Taxonomic_level",
    values_to = "Count"
  ) %>%
  dplyr::mutate(
    Taxonomic_level = factor(Taxonomic_level,
                             levels = c("N_orders", "N_families", "N_genera"),
                             labels = c("Orders", "Families", "Genera"))
  )

plot_richness <- ggplot(richness_plot_data, 
                        aes(x = Taxonomic_level, y = Count, 
                            fill = Device_type)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(
    values = custom_colors[c("FAIRD", "ID")],
    labels = c("FAIRD", "Insect Detect")
  ) +
  labs(
    title = "Taxonomic Richness Comparison",
    subtitle = "FAIRD captures 980% more genera than Insect Detect",
    x = "Taxonomic Level",
    y = "Number of Taxa",
    fill = "Device Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, max(richness_plot_data$Count) * 1.15)

print(plot_richness)

# Save plot
ggsave("outputs/H3_Richness.png", 
       plot = plot_richness, 
       width = 8, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H3_Richness.pdf",
       plot = plot_richness,
       width = 8, height = 6, device = cairo_pdf)

cat("Richness plot saved:\n")
cat("   - outputs/H3_Richness.png\n")
cat("   - H3_Richness.pdf\n\n")

# ==============================================================================
# SECTION 6: VISUALIZATION - COMPOSITION STACKED BAR
# ==============================================================================

cat("\n=== SECTION 6: COMPOSITION STACKED BAR ===\n\n")

# Prepare data - top orders only for clarity
top_orders <- composition_etraps %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(Total_count = sum(Count), .groups = "drop") %>%
  dplyr::arrange(desc(Total_count)) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::pull(Order)

composition_plot_data <- composition_etraps %>%
  dplyr::mutate(
    Order_grouped = ifelse(Order %in% top_orders, Order, "Other")
  ) %>%
  dplyr::group_by(Device_type, Order_grouped) %>%
  dplyr::summarise(Percentage = sum(Percentage), .groups = "drop") %>%
  dplyr::mutate(
    Order_grouped = factor(Order_grouped, 
                           levels = c(top_orders, "Other"))
  )

# Use order_colors from custom_colors.R
# Note: custom_colors.R already defines order_colors with the main orders
# We'll use those colors and add any missing ones if needed

plot_composition <- ggplot(composition_plot_data, 
                           aes(x = Device_type, y = Percentage, 
                               fill = Order_grouped)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(Percentage > 5, 
                               sprintf("%.0f%%", Percentage), "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = order_colors, name = "Order") +
  scale_x_discrete(labels = c("FAIRD", "Insect Detect")) +
  labs(
    title = "Taxonomic Composition Comparison",
    subtitle = "ID specialized (>95% Diptera), FAIRD generalist sampler",
    x = "Device Type",
    y = "Percentage of Captures (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

print(plot_composition)

# Save plot
ggsave("outputs/H3_Composition.png", 
       plot = plot_composition, 
       width = 10, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H3_Composition.pdf",
       plot = plot_composition,
       width = 10, height = 6, device = cairo_pdf)

cat("Composition plot saved:\n")
cat("   - outputs/H3_Composition.png\n")
cat("   - H3_Composition.pdf\n\n")

# ==============================================================================
# SECTION 6B: VISUALIZATION - DIPTERA FAMILY COMPOSITION
# ==============================================================================

cat("\n=== SECTION 6B: DIPTERA FAMILY COMPOSITION PLOT ===\n\n")

# Prepare data for visualization - top families within Diptera
top_diptera_families <- diptera_families %>%
  dplyr::group_by(Family) %>%
  dplyr::summarise(Total_count = sum(Count), .groups = "drop") %>%
  dplyr::arrange(desc(Total_count)) %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::pull(Family)

diptera_plot_data <- diptera_families %>%
  dplyr::mutate(
    Family_grouped = ifelse(Family %in% top_diptera_families, Family, "Other")
  ) %>%
  dplyr::group_by(Device_type, Family_grouped) %>%
  dplyr::summarise(
    Count = sum(Count),
    Percentage = sum(Percentage), 
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Family_grouped = factor(Family_grouped, 
                            levels = c(top_diptera_families, "Other"))
  )

# Create color palette for families
family_colors <- c(
  "Syrphidae" = "#E31A1C",      # Red - highlight Syrphidae
  "Muscidae" = "#377EB8",       # Blue
  "Calliphoridae" = "#4DAF4A",  # Green
  "Tachinidae" = "#984EA3",     # Purple
  "Anthomyiidae" = "#FF7F00",   # Orange
  "Drosophilidae" = "#FFFF33",  # Yellow
  "Other" = "#999999"           # Gray
)

# Adjust colors based on actual families present
actual_families <- unique(diptera_plot_data$Family_grouped)
available_colors <- family_colors[names(family_colors) %in% actual_families]

# Add any missing families with default colors
if(length(available_colors) < length(actual_families)) {
  missing_families <- setdiff(actual_families, names(available_colors))
  additional_colors <- RColorBrewer::brewer.pal(
    min(length(missing_families), 8), "Set2"
  )[1:length(missing_families)]
  names(additional_colors) <- missing_families
  available_colors <- c(available_colors, additional_colors)
}

plot_diptera_families <- ggplot(diptera_plot_data, 
                                aes(x = Device_type, y = Percentage, 
                                    fill = Family_grouped)) +
  geom_col(position = "stack", width = 0.5) +
  geom_text(aes(label = ifelse(Percentage > 8, 
                               sprintf("%.0f%%", Percentage), "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = available_colors, name = "Family") +
  scale_x_discrete(labels = c("FAIRD", "Insect Detect")) +
  labs(
    title = "Diptera Family Composition: FAIRD vs Insect Detect",
    subtitle = "Testing Syrphidae specialization claim for ID",
    x = "Device Type",
    y = "Percentage within Diptera (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

print(plot_diptera_families)

# Save plot
ggsave("outputs/H3_Diptera_Families.png", 
       plot = plot_diptera_families, 
       width = 10, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H3_Diptera_Families.pdf",
       plot = plot_diptera_families,
       width = 10, height = 6, device = cairo_pdf)

cat("Diptera family composition plot saved:\n")
cat("   - outputs/H3_Diptera_Families.png\n")
cat("   - H3_Diptera_Families.pdf\n\n")

# =============================================================================
# SECTION 6C: TAXONOMIC COMPOSITION (DIPTERA) - CORRECTED VERSION
# =============================================================================
# Objective: Visualize Diptera composition by device, with hierarchical
#            stacking order (Family on top, 'Other' on bottom).
#            Legend is reversed to match the stack.
#            Groups <= 1% are lumped, counts in subtitle.
#            Uses colors defined in 'custom_colors.R'.
# =============================================================================

cat("\n=== SECTION 6C: DIPTERA TAXONOMIC COMPOSITION (Corrected) ===\n\n")

# 1. Filter for Diptera data from the main dataset
diptera_data <- etraps_data %>%
  dplyr::filter(Order == "Diptera")

# 2. Apply the function to create the base plot label
# ✅ FIX: Use rowwise() to apply get_best_classification() to each row
diptera_plot_data <- diptera_data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    plot_label = get_best_classification(dplyr::cur_data())
  ) %>%
  dplyr::ungroup()

# 3. Identify global categories that are <= 1% of the total
cat("  Identifying global taxonomic groups <= 1% of total observations...\n")
total_diptera_obs <- nrow(diptera_plot_data)
labels_to_lump <- diptera_plot_data %>%
  dplyr::count(plot_label) %>%
  dplyr::mutate(prop = n / total_diptera_obs) %>%
  dplyr::filter(prop <= 0.01) %>%
  dplyr::pull(plot_label)
cat(sprintf("  Identified %d globally rare taxa (<= 1%% total).\n", length(labels_to_lump)))


# 4. Count unique rare taxa *within* each device (for subtitle)
cat("  Counting unique rare taxa *within* each device (for subtitle)...\n")
all_devices <- unique(diptera_plot_data$Device_type)
per_device_lump_counts <- diptera_plot_data %>%
  dplyr::filter(plot_label %in% labels_to_lump) %>%
  dplyr::group_by(Device_type) %>%
  dplyr::summarise(n_lumped = n_distinct(plot_label), .groups = "drop") %>%
  dplyr::right_join(tibble(Device_type = all_devices), by = "Device_type") %>%
  tidyr::replace_na(list(n_lumped = 0))
cat("  Per-device rare taxa counts:\n")
print(per_device_lump_counts %>% dplyr::select(Device_type, n_lumped))


# 5. Create the dynamic subtitle string
subtitle_string_parts <- per_device_lump_counts %>%
  dplyr::mutate(label_part = sprintf("%s (%d taxa)", Device_type, n_lumped)) %>%
  dplyr::pull(label_part)
plot_subtitle <- paste("Groups <= 1% lumped into 'Other'. Counts:", 
                       paste(subtitle_string_parts, collapse = ", "))
cat(sprintf("  Generated subtitle: %s\n", plot_subtitle))


# 6. Apply lumping to a *single* 'Other' category
cat("  Collapsing all rare groups into a single 'Other' category...\n")
diptera_plot_data <- diptera_plot_data %>%
  dplyr::mutate(
    plot_label_lumped = forcats::fct_collapse(
      factor(plot_label),
      "Other" = labels_to_lump
    )
  )

# --- (REVISED) STEP 6.5: ROBUSTLY REORDER FACTOR LEVELS ---
cat("  Robustly reordering factor levels for hierarchical stacking...\n")

# Get all current levels
current_levels <- levels(diptera_plot_data$plot_label_lumped)

# Create a data frame of levels and assign a numeric rank for sorting
# Rank 1 = bottom of bar, Rank 7 = top of bar
level_ranks <- tibble(label = current_levels) %>%
  dplyr::mutate(
    rank_order = dplyr::case_when(
      label == "Other" ~ 1,
      stringr::str_detect(label, stringr::fixed("(Class)")) ~ 2,
      stringr::str_detect(label, stringr::fixed("(Order)")) ~ 3,
      stringr::str_detect(label, stringr::fixed("(Suborder)")) ~ 4,
      stringr::str_detect(label, stringr::fixed("(Infraorder)")) ~ 5,
      stringr::str_detect(label, stringr::fixed("(Superfamily)")) ~ 6,
      TRUE ~ 7 # All others (Families, Genera, etc.) go on top
    )
  ) %>%
  # Sort by the new rank, then alphabetically within that rank
  dplyr::arrange(rank_order, label)

# Pull the sorted labels as the new factor order
new_level_order <- level_ranks$label

# Apply the new, robust order to the factor
diptera_plot_data <- diptera_plot_data %>%
  dplyr::mutate(plot_label_lumped = forcats::fct_relevel(plot_label_lumped, new_level_order))

cat("  Factor levels reordered successfully.\n")
# -------------------------------------------------------------------------


# 7. Create the custom color map
cat("\n  Applying custom color map.\n")

# Get all unique final levels (now in the correct order)
plot_levels <- levels(diptera_plot_data$plot_label_lumped)

other_color <- order_colors["Other"]
main_levels <- setdiff(plot_levels, "Other")
main_level_count <- length(main_levels)

if (main_level_count > 0) {
  main_colors <- colorRampPalette(RColorBrewer::brewer.pal(min(main_level_count, 11), "Set3"))(main_level_count)
  diptera_plot_colors <- setNames(main_colors, main_levels)
} else {
  diptera_plot_colors <- c() 
}

diptera_plot_colors["Other"] <- other_color


# 8. Generate the composition plot (with corrected stacking)
p_diptera_composition <- ggplot(diptera_plot_data, 
                                aes(x = Device_type, fill = plot_label_lumped)) +
  
  # --- (KEY CHANGE HERE) ---
  # Use position_fill(reverse = TRUE) to stack in the
  # *same* order as the factor levels (Rank 1 at bottom, Rank 7 at top).
  geom_bar(position = position_fill(reverse = TRUE)) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  
  # scale_fill_manual respects the factor order
  scale_fill_manual(values = diptera_plot_colors) +
  
  # Reverse the legend order to match the stacking order
  # (Top of stack = Top of legend)
  guides(fill = guide_legend(reverse = TRUE)) +
  
  labs(
    title = "Taxonomic Composition of Diptera up to Family-level by Device",
    subtitle = plot_subtitle, 
    x = "Device Type",
    y = "Proportion",
    fill = "Taxonomic Group" 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    legend.position = "right"
  )

print(p_diptera_composition)

cat("\n✅ Diptera composition plot (hierarchically stacked) generated.\n")


# 9. Save plot
ggsave("outputs/H3_Diptera_Composition.png", 
       plot = p_diptera_composition, 
       width = 10, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H3_Diptera_Composition.pdf",
       plot = p_diptera_composition,
       width = 10, height = 6, device = cairo_pdf)

cat("\n✅ Diptera composition plot saved:\n")
cat("   - outputs/H3_Diptera_Composition.png\n")
cat("   - H3_Diptera_Composition.pdf\n\n")

# =============================================================================
# SECTION 6D: DIPTERA SUPERFAMILY COMPOSITION (Superfamily-level only)
# =============================================================================
# Objective: Visualize Diptera composition by device at SUPERFAMILY level
#            This avoids mixing taxonomic levels and focuses on Syrphoidea
#            Uses the same hierarchical approach as 6C but stops at Superfamily
# =============================================================================

cat("\n=== SECTION 6D: DIPTERA SUPERFAMILY COMPOSITION (Superfamily-level) ===\n\n")

# 1. Filter for Diptera data from the main dataset
diptera_data_sf <- etraps_data %>%
  filter(Order == "Diptera")

# 2. Define function to get best classification UP TO Superfamily level
get_best_classification_superfamily <- function(Superfamily, Infraorder, 
                                                Suborder, Order, Class) {
  # Vectorized function - stops at Superfamily as most detailed level
  case_when(
    # Level 1: Use Superfamily if available (MOST DETAILED)
    !is.na(Superfamily) & Superfamily != "#N/C" & Superfamily != "" ~ Superfamily,
    
    # Level 2: Use Infraorder if Superfamily not available
    !is.na(Infraorder) & Infraorder != "#N/C" & Infraorder != "" ~ 
      paste0(Infraorder, " (Infraorder)"),
    
    # Level 3: Use Suborder if Infraorder not available
    !is.na(Suborder) & Suborder != "#N/C" & Suborder != "" ~ 
      paste0(Suborder, " (Suborder)"),
    
    # Level 4: Use Order if Suborder not available
    !is.na(Order) & Order != "#N/C" & Order != "" ~ 
      paste0(Order, " (Order)"),
    
    # Level 5: Use Class as last resort
    !is.na(Class) & Class != "#N/C" & Class != "" ~ 
      paste0(Class, " (Class)"),
    
    # Default: If nothing is available
    TRUE ~ "Unclassified"
  )
}

# 3. Apply the function to create plot labels (Superfamily-level)
diptera_plot_data_sf <- diptera_data_sf %>%
  mutate(
    plot_label = get_best_classification_superfamily(
      Superfamily, 
      Infraorder, 
      Suborder, 
      Order, 
      Class
    )
  )

cat("  Classification levels used:\n")
level_counts <- diptera_plot_data_sf %>%
  mutate(
    level = case_when(
      str_detect(plot_label, "Infraorder") ~ "Infraorder",
      str_detect(plot_label, "Suborder") ~ "Suborder",
      str_detect(plot_label, "Order") ~ "Order",
      str_detect(plot_label, "Class") ~ "Class",
      plot_label == "Unclassified" ~ "Unclassified",
      TRUE ~ "Superfamily"
    )
  ) %>%
  count(level) %>%
  arrange(match(level, c("Superfamily", "Infraorder", "Suborder", 
                         "Order", "Class", "Unclassified")))
print(level_counts)

# 4. Identify global categories that are <= 1% of the total
cat("\n  Identifying global taxonomic groups <= 1% of total observations...\n")
total_diptera_obs_sf <- nrow(diptera_plot_data_sf)
labels_to_lump_sf <- diptera_plot_data_sf %>%
  count(plot_label) %>%
  mutate(prop = n / total_diptera_obs_sf) %>%
  filter(prop <= 0.01) %>%
  pull(plot_label)
cat(sprintf("  Identified %d globally rare taxa (<= 1%% total).\n", 
            length(labels_to_lump_sf)))

# 5. Count unique rare taxa within each device (for subtitle)
cat("  Counting unique rare taxa within each device...\n")
all_devices_sf <- unique(diptera_plot_data_sf$Device_type)
per_device_lump_counts_sf <- diptera_plot_data_sf %>%
  filter(plot_label %in% labels_to_lump_sf) %>%
  group_by(Device_type) %>%
  summarise(n_lumped = n_distinct(plot_label)) %>%
  right_join(tibble(Device_type = all_devices_sf), by = "Device_type") %>%
  replace_na(list(n_lumped = 0))

cat("  Per-device rare taxa counts:\n")
# FIX: Use dplyr:: prefix to avoid conflicts
print(per_device_lump_counts_sf %>% dplyr::select(Device_type, n_lumped))

# 6. Create dynamic subtitle
subtitle_string_parts_sf <- per_device_lump_counts_sf %>%
  mutate(label_part = sprintf("%s (%d taxa)", Device_type, n_lumped)) %>%
  pull(label_part)
plot_subtitle_sf <- paste("Groups <= 1% lumped into 'Other'. Counts:", 
                          paste(subtitle_string_parts_sf, collapse = ", "))
cat(sprintf("  Generated subtitle: %s\n", plot_subtitle_sf))

# 7. Apply lumping to a single 'Other' category
cat("  Collapsing all rare groups into 'Other'...\n")
diptera_plot_data_sf <- diptera_plot_data_sf %>%
  mutate(
    plot_label_lumped = fct_collapse(
      factor(plot_label),
      "Other" = labels_to_lump_sf
    )
  )

# 8. Robustly reorder factor levels for hierarchical stacking
cat("  Reordering factor levels for hierarchical stacking...\n")

current_levels_sf <- levels(diptera_plot_data_sf$plot_label_lumped)

level_ranks_sf <- tibble(label = current_levels_sf) %>%
  mutate(
    rank_order = case_when(
      label == "Other" ~ 1,                              # Bottom
      str_detect(label, fixed("(Class)")) ~ 2,
      str_detect(label, fixed("(Order)")) ~ 3,
      str_detect(label, fixed("(Suborder)")) ~ 4,
      str_detect(label, fixed("(Infraorder)")) ~ 5,
      TRUE ~ 6  # Superfamilies at top (most detailed level)
    )
  ) %>%
  arrange(rank_order, label)

new_level_order_sf <- level_ranks_sf$label

diptera_plot_data_sf <- diptera_plot_data_sf %>%
  mutate(plot_label_lumped = fct_relevel(plot_label_lumped, new_level_order_sf))

cat("  Factor levels reordered successfully.\n")

# 9. Create custom color map with Syrphoidea highlighted
cat("\n  Applying custom color map (Syrphoidea highlighted in red).\n")

plot_levels_sf <- levels(diptera_plot_data_sf$plot_label_lumped)
other_color_sf <- order_colors["Other"]
main_levels_sf <- setdiff(plot_levels_sf, "Other")
main_level_count_sf <- length(main_levels_sf)

# Create base colors
if (main_level_count_sf > 0) {
  main_colors_sf <- colorRampPalette(
    brewer.pal(min(main_level_count_sf, 11), "Set3")
  )(main_level_count_sf)
  diptera_plot_colors_sf <- setNames(main_colors_sf, main_levels_sf)
} else {
  diptera_plot_colors_sf <- c() 
}

# Highlight Syrphoidea in red if present
if ("Syrphoidea" %in% names(diptera_plot_colors_sf)) {
  diptera_plot_colors_sf["Syrphoidea"] <- "#E31A1C"  # Red for Syrphoidea
  cat("  ✅ Syrphoidea highlighted in RED\n")
}

# Add Other color
diptera_plot_colors_sf["Other"] <- other_color_sf

# 10. Generate the composition plot
p_diptera_superfamily <- ggplot(diptera_plot_data_sf, 
                                aes(x = Device_type, fill = plot_label_lumped)) +
  
  # Stack in same order as factor levels (reverse = TRUE for correct order)
  geom_bar(position = position_fill(reverse = TRUE)) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  
  # Use custom colors
  scale_fill_manual(values = diptera_plot_colors_sf) +
  
  # Reverse legend to match stacking order (top of stack = top of legend)
  guides(fill = guide_legend(reverse = TRUE)) +
  
  labs(
    title = "Diptera Taxonomic Composition up to Superfamily-level by Device",
    subtitle = plot_subtitle_sf,
    x = "Device Type",
    y = "Proportion",
    fill = "Taxonomic Group" 
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

print(p_diptera_superfamily)

cat("\n✅ Diptera superfamily-level composition plot generated.\n")

# 11. Save plot
ggsave("outputs/H3_Diptera_Superfamily_Only.png", 
       plot = p_diptera_superfamily, 
       width = 10, height = 6, dpi = 300, bg = "white")

ggsave("outputs/H3_Diptera_Superfamily_Only.pdf",
       plot = p_diptera_superfamily,
       width = 10, height = 6, device = cairo_pdf)

cat("✅ Diptera superfamily-level composition plot saved:\n")
cat("   - outputs/H3_Diptera_Superfamily_Only.png\n")
cat("   - H3_Diptera_Superfamily_Only.pdf\n\n")

cat("💡 INTERPRETATION NOTE:\n")
cat("   This plot uses SUPERFAMILY as the most detailed level\n")
cat("   → Avoids mixing Family with Superfamily labels\n")
cat("   → Focuses on Syrphoidea (hover flies) vs other Diptera groups\n")
cat("   → 'Brachycera (Suborder)' represents specimens not identified to Superfamily\n")
cat("   → Based on ID AI categories: 'hfly_*' = Syrphoidea, 'fly*' = Non-Syrphoidea\n\n")

# 12. Compare with Section 2C results
cat("📊 COMPARISON WITH SECTION 2C:\n")
cat("   Section 2C: Analyzes ONLY specimens with Superfamily identified\n")
cat("   Section 6D: Shows ALL specimens, using best available level up to Superfamily\n")
cat("   → Both approaches valid for different purposes\n")
cat("   → Section 2C gives exact % Syrphoidea among identified specimens\n")
cat("   → Section 6D shows complete composition including unidentified specimens\n\n")

# ==============================================================================
# SECTION 7: DEVICE RELIABILITY COMPARISON
# ==============================================================================

cat("\n=== SECTION 7: DEVICE RELIABILITY ===\n\n")

# Check device failures
device_summary <- etraps_data %>%
  dplyr::group_by(Device_type, Device) %>%
  dplyr::summarise(
    N_detections = n(),
    N_days_active = n_distinct(Date),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Device_type, Device)

print(device_summary)

cat("\n# Reliability Assessment:\n")

# ID3 failure
id3_records <- device_summary %>%
  dplyr::filter(Device == "ID3") %>%
  nrow()

if(id3_records == 0) {
  cat("   âŒ ID3: COMPLETE FAILURE (0 detections)\n")
} else {
  cat("   âš ï¸ ID3: Minimal data (device issues)\n")
}

# FAIRD robustness
faird_min <- device_summary %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::pull(N_detections) %>%
  min()

cat(sprintf(" - FAIRD: All devices functional (min: %d detections)\n\n", 
            faird_min))

# ==============================================================================
# SECTION 8: DIVERSITY INDICES (OPTIONAL) - CORREGIDO con pick()
# ==============================================================================

cat("\n=== SECTION 8: DIVERSITY INDICES (OPTIONAL) ===\n\n")

# Calculate for e-traps only
diversity_etraps <- etraps_data %>%
  dplyr::group_by(Device_type, Order) %>%
  dplyr::summarise(Abundance = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0) %>%
  dplyr::rowwise() %>% # Necesario para calcular diversidad fila por fila (por Device_type)
  dplyr::mutate(
    # USANDO pick() para seleccionar todas las columnas de ORDEN excepto Device_type
    Abundance_Matrix = dplyr::pick(dplyr::everything(), -Device_type),
    Shannon = vegan::diversity(Abundance_Matrix, index = "shannon"),
    Simpson = vegan::diversity(Abundance_Matrix, index = "simpson")
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Device_type, Shannon, Simpson)

cat("Diversity indices (Order level):\n")
print(diversity_etraps)

cat("\n")

# ==============================================================================
# SECTION 9: TAXONOMIC RESOLUTION SUMMARY BY ORDER AND DEVICE
# ==============================================================================
# Purpose: Calculate the percentage of individuals identified at each taxonomic
#          level (Genus, Family, Superfamily, Order-or-higher) for each device
#          type and order. Verifies/corrects PCV14 Table (FAIRD Genus ID rates).
#
# NOTE: "#N/C" values in e-trap CSVs are NOT converted to NA on load (no na=
#       argument in read_csv for e-trap files). They must be treated as invalid
#       here explicitly alongside standard NA checks.
# ==============================================================================

cat("\n=== SECTION 9: TAXONOMIC RESOLUTION BY ORDER AND DEVICE ===\n\n")

# --- Helper: classify each individual by its finest resolved taxonomic level ---

# Valid value = not NA AND not the "#N/C" placeholder
is_valid <- function(x) !is.na(x) & x != "#N/C" & x != ""

resolution_data <- etraps_data %>%
  dplyr::mutate(
    finest_level = dplyr::case_when(
      is_valid(Genus)       ~ "Genus",
      is_valid(Family)      ~ "Family",
      is_valid(Superfamily) ~ "Superfamily",
      is_valid(Infraorder)  ~ "Infraorder",
      is_valid(Suborder)    ~ "Suborder",
      is_valid(Order)       ~ "Order",
      TRUE                  ~ "Class_or_unresolved"
    ),
    finest_level = factor(finest_level,
                          levels = c("Genus", "Family", "Superfamily",
                                     "Infraorder", "Suborder", "Order",
                                     "Class_or_unresolved"))
  )

# ==============================================================================
# TABLE A: Resolution by Device x Order (mirrors PCV14 table)
# ==============================================================================

cat("--- TABLE A: Resolution by Device x Order ---\n\n")

resolution_by_order <- resolution_data %>%
  dplyr::group_by(Device_type, Order, finest_level) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Device_type, Order) %>%
  dplyr::mutate(
    total = sum(n),
    pct   = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup()

# Pivot wider for a compact table: one row per Device x Order
resolution_table_order <- resolution_by_order %>%
  dplyr::select(Device_type, Order, finest_level, n, pct, total) %>%
  tidyr::pivot_wider(
    names_from  = finest_level,
    values_from = c(n, pct),
    values_fill = list(n = 0, pct = 0.0)
  ) %>%
  # Keep only columns that actually exist (avoids errors if some levels absent)
  dplyr::select(
    Device_type, Order, total,
    dplyr::any_of(c(
      "n_Genus",            "pct_Genus",
      "n_Family",           "pct_Family",
      "n_Superfamily",      "pct_Superfamily",
      "n_Infraorder",       "pct_Infraorder",
      "n_Suborder",         "pct_Suborder",
      "n_Order",            "pct_Order",
      "n_Class_or_unresolved", "pct_Class_or_unresolved"
    ))
  ) %>%
  dplyr::arrange(Device_type, Order)

print(resolution_table_order, n = 50)

# ==============================================================================
# TABLE B: Overall resolution per Device (all orders combined)
# ==============================================================================

cat("\n--- TABLE B: Overall Resolution per Device (all Orders combined) ---\n\n")

resolution_overall <- resolution_data %>%
  dplyr::group_by(Device_type, finest_level) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Device_type) %>%
  dplyr::mutate(
    total = sum(n),
    pct   = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Device_type, finest_level)

print(resolution_overall)

# ==============================================================================
# TABLE C: Condensed summary matching PCV14 format
#          (Genus ID %, Family only %, Superfamily/lower %) per Order per Device
# ==============================================================================

cat("\n--- TABLE C: Condensed summary (PCV14-style) ---\n")
cat("    [Genus ID %] | [Family only %] | [Superfamily or lower %]\n\n")

resolution_condensed <- resolution_data %>%
  dplyr::mutate(
    resolution_class = dplyr::case_when(
      finest_level == "Genus"                 ~ "Genus_ID",
      finest_level == "Family"                ~ "Family_only",
      TRUE                                    ~ "Superfamily_or_lower"
    )
  ) %>%
  dplyr::group_by(Device_type, Order, resolution_class) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Device_type, Order) %>%
  dplyr::mutate(
    total = sum(n),
    pct   = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    names_from  = resolution_class,
    values_from = c(n, pct),
    values_fill = list(n = 0, pct = 0.0)
  ) %>%
  dplyr::arrange(Device_type, dplyr::desc(total)) %>%
  dplyr::select(
    Device_type, Order, total,
    dplyr::any_of(c(
      "pct_Genus_ID", "n_Genus_ID",
      "pct_Family_only", "n_Family_only",
      "pct_Superfamily_or_lower", "n_Superfamily_or_lower"
    ))
  )

print(resolution_condensed, n = 50)

# ==============================================================================
# TABLE D: OVERALL condensed per Device (PCV14 "Overall" row)
# ==============================================================================

cat("\n--- TABLE D: Overall condensed per Device (PCV14 'Overall' row) ---\n\n")

resolution_overall_condensed <- resolution_data %>%
  dplyr::mutate(
    resolution_class = dplyr::case_when(
      finest_level == "Genus"  ~ "Genus_ID",
      finest_level == "Family" ~ "Family_only",
      TRUE                     ~ "Superfamily_or_lower"
    )
  ) %>%
  dplyr::group_by(Device_type, resolution_class) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Device_type) %>%
  dplyr::mutate(
    total = sum(n),
    pct   = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    names_from  = resolution_class,
    values_from = c(n, pct),
    values_fill = list(n = 0, pct = 0.0)
  )

print(resolution_overall_condensed)

cat("\n--- Sanity check: pct columns should sum to ~100% per row ---\n")
resolution_overall_condensed %>%
  dplyr::mutate(
    pct_sum = rowSums(dplyr::across(dplyr::starts_with("pct_")), na.rm = TRUE)
  ) %>%
  dplyr::select(Device_type, total, pct_sum) %>%
  print()

# ==============================================================================
# SAVE TO CSV
# ==============================================================================

readr::write_csv(resolution_table_order,
                 "outputs/Table_H3_Taxonomic_Resolution_by_Order.csv")

readr::write_csv(resolution_condensed,
                 "outputs/Table_H3_Taxonomic_Resolution_Condensed.csv")

cat("\nFiles saved:\n")
cat("  - outputs/Table_H3_Taxonomic_Resolution_by_Order.csv\n")
cat("  - outputs/Table_H3_Taxonomic_Resolution_Condensed.csv\n\n")

# ==============================================================================
# VERIFICATION: Diptera dominance by Site and Habitat (Ambient)
# ==============================================================================

cat("\n=== VERIFICATION: DIPTERA DOMINANCE BY HABITAT AND SITE ===\n\n")

# --- 1. FAIRD composition by Ambient (Maize vs Meadow) ---

cat("--- TABLE 1: FAIRD Order composition by Habitat (Ambient) ---\n\n")

faird_composition_ambient <- etraps_data %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::group_by(Ambient, Order) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Ambient) %>%
  dplyr::mutate(
    total   = sum(n),
    pct     = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Ambient, dplyr::desc(pct))

print(faird_composition_ambient, n = 30)

# Extract Diptera % by Ambient for direct comparison
cat("\n--- KEY: Diptera % by Habitat ---\n")
faird_composition_ambient %>%
  dplyr::filter(Order == "Diptera") %>%
  dplyr::select(Ambient, n, total, pct) %>%
  dplyr::rename(Diptera_n = n, Total_FAIRD_n = total, Diptera_pct = pct) %>%
  print()

# --- 2. FAIRD composition by Site ---

cat("\n--- TABLE 2: FAIRD Order composition by Site ---\n\n")

faird_composition_site <- etraps_data %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::group_by(Site, Ambient, Order) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    total   = sum(n),
    pct     = round(n / total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Ambient, Site, dplyr::desc(pct))

print(faird_composition_site, n = 40)

# Extract Diptera % by Site for direct comparison
cat("\n--- KEY: Diptera % by Site ---\n")
faird_composition_site %>%
  dplyr::filter(Order == "Diptera") %>%
  dplyr::select(Site, Ambient, n, total, pct) %>%
  dplyr::rename(Diptera_n = n, Total_FAIRD_n = total, Diptera_pct = pct) %>%
  print()

# --- 3. Overall reference (all sites) ---

cat("\n--- REFERENCE: Diptera % across ALL FAIRD sites (PCV14 value) ---\n")
etraps_data %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    total = sum(n),
    pct   = round(n / total * 100, 1)
  ) %>%
  dplyr::filter(Order == "Diptera") %>%
  dplyr::select(Order, n, total, pct) %>%
  print()

cat("\n--- INTERPRETATION GUIDE ---\n")
cat("PCV14 claim: '87.4% of FAIRD catches' in 'maize sites'\n")
cat("  -> If Maize Diptera pct ≈ 87.4%: claim is CORRECT (value is maize-specific)\n")
cat("  -> If Maize Diptera pct differs from 87.4%: claim is INCORRECT\n")
cat("     (87.4% is the OVERALL figure, not maize-specific)\n\n")

cat("=== VERIFICATION COMPLETE ===\n\n")

cat("=== SECTION 9 COMPLETE ===\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("H3 ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n")
cat("- outputs/Table3_H3_Taxonomic_Results.csv\n")
cat("- outputs/H3_Richness.png (Bar plot)\n")
cat("- outputs/H3_Richness.pdf\n")
cat("- outputs/H3_Composition.png (Stacked bar)\n")
cat("- outputs/H3_Composition.pdf\n")
cat("- outputs/H3_Diptera_Families.png (Diptera families)\n")
cat("- outputs/H3_Diptera_Families.pdf\n")
cat("- outputs/H3_Diptera_Superfamilies.png (Diptera superfamilies)\n")
cat("- outputs/H3_Diptera_Superfamilies.pdf\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("- FAIRD: %d genera (GENERALIST sampler)\n", 
            richness_summary$N_genera[1]))
cat(sprintf("- ID: %d genera (DIPTERA specialist, %.1f%% Diptera)\n", 
            richness_summary$N_genera[2], id_diptera_pct))
cat(sprintf("- FAIRD captures %d%% MORE genera than ID\n",
            richness_summary$Genera_advantage_vs_ID[1]))

# Add Syrphidae findings if data is available (Family level)
if(id_syrphidae_pct > 0 | faird_syrphidae_pct > 0) {
  cat(sprintf("- Syrphidae (Family) in ID: %.1f%% (vs FAIRD: %.1f%%)\n",
              id_syrphidae_pct, faird_syrphidae_pct))
  if(id_syrphidae_pct > faird_syrphidae_pct) {
    cat("- ID Syrphidae specialization CONFIRMED at Family level\n")
  } else {
    cat("- No strong Syrphidae specialization at Family level\n")
  }
}

# Add Syrphoidea findings if data is available (Superfamily level - more appropriate)
if(id_syrphoidea_pct > 0 | faird_syrphoidea_pct > 0) {
  cat(sprintf("- Syrphoidea (Superfamily) in ID: %.1f%% (vs FAIRD: %.1f%%)\n",
              id_syrphoidea_pct, faird_syrphoidea_pct))
  if(id_syrphoidea_pct > faird_syrphoidea_pct) {
    cat("- ID hover fly specialization CONFIRMED at Superfamily level\n")
    cat("- Taxonomically appropriate analysis\n")
  } else {
    cat("- No strong Syrphoidea specialization at Superfamily level\n")
  }
}

cat("- ID3 complete failure vs FAIRD robust performance\n")
cat("- Trade-off: Automation efficiency vs taxonomic breadth\n\n")

# ==============================================================================
cat("========================================\n\n")

cat("\n========================================\n")
cat("H3 COMPLETE\n")
cat("========================================\n\n")

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")
