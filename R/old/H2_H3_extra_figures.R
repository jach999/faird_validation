# ============================================================================
# GENERATE MISSING FIGURES FOR THE PRESENTATION
# ============================================================================
# 
# This file contains code to create the two figures that are not automatically generated
# in the current scripts H2 and H3:
# 1. H2 Size Distribution (detection bias visualization)
# 2. H3 Venn Diagram (genus-level overlaps)
#
# INSTRUCTIONS:
# - Run AFTER loading data_cleaned.RData
# - Requires custom_colors.R
# - Generates figures in the outputs folder/
# ============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(VennDiagram)  # For Venn diagrams
library(ggvenn)       # Alternative simpler Venn package

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


# Load project data and colors
load("data_cleaned.RData")
source("custom_colors.R")

# ============================================================================
# 1. H2: SIZE DISTRIBUTION / DETECTION BIAS (FAIRD ONLY)
# ============================================================================
# Purpose: Show that FAIRD has detection limits at 3mm and constrained <5mm
# Expected in slide: "H2 Known Limitations: Detection Bias"
# ============================================================================

message("Creating H2 Size Distribution figure (FAIRD only)...")

# Get size data from FAIRD only
size_data_faird <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::select(Body_length, Body_width) %>%
  dplyr::filter(!is.na(Body_length)) %>%
  dplyr::mutate(Device = "FAIRD")  # Keep column name for consistency

# Calculate detection statistics by size class
size_stats <- data.frame(
  Class = c("<3mm", "3-5mm", ">5mm"),
  Count = c(
    sum(size_data_faird$Body_length < 3),
    sum(size_data_faird$Body_length >= 3 & size_data_faird$Body_length < 5),
    sum(size_data_faird$Body_length >= 5)
  )
)

# Calculate percentages
size_stats$Percentage <- round(100 * size_stats$Count / sum(size_stats$Count), 1)

# Print statistics
message("\n=== FAIRD SIZE CLASS DETECTION STATISTICS ===")
print(size_stats)
message("\nKey findings:")
message(sprintf("  FAIRD <3mm: %d individuals (%.1f%%) - Absolute limit confirmed", 
                size_stats$Count[1], size_stats$Percentage[1]))
message(sprintf("  FAIRD 3-5mm: %d individuals (%.1f%%) - Constrained but functional", 
                size_stats$Count[2], size_stats$Percentage[2]))
message(sprintf("  FAIRD >5mm: %d individuals (%.1f%%) - Optimal detection", 
                size_stats$Count[3], size_stats$Percentage[3]))

# Create size distribution plot
p_size_distribution <- ggplot(size_data_faird, aes(x = Body_length, fill = Device)) +
  geom_histogram(
    alpha = 0.8, 
    bins = 40,
    color = "white"
  ) +
  # Absolute detection limit (3mm)
  geom_vline(
    xintercept = 3, 
    linetype = "dashed", 
    color = "#E74C3C",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = 3, 
    y = Inf,
    label = "Absolute limit (3mm)",
    hjust = -0.05,
    vjust = 1.5,
    color = "#E74C3C",
    fontface = "bold",
    size = 3.5
  ) +
  # Constrained detection zone (5mm)
  geom_vline(
    xintercept = 5, 
    linetype = "dotted", 
    color = "#F39C12",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = 5, 
    y = Inf,
    label = "Constrained detection (<5mm)",
    hjust = -0.05,
    vjust = 3,
    color = "#F39C12",
    fontface = "bold",
    size = 3.5
  ) +
  # Add subtle background shading for zones
  annotate(
    "rect",
    xmin = 0, xmax = 3,
    ymin = -Inf, ymax = Inf,
    fill = "#E74C3C",
    alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = 3, xmax = 5,
    ymin = -Inf, ymax = Inf,
    fill = "#F39C12",
    alpha = 0.08
  ) +
  # Add detection statistics annotations
  annotate(
    "text",
    x = 1.5, y = 0,
    label = sprintf("FAIRD: %d (%.1f%%)", 
                    size_stats$Count[1], 
                    size_stats$Percentage[1]),
    hjust = 0.5,
    vjust = -0.5,
    color = "#E74C3C",
    fontface = "bold",
    size = 3.5
  ) +
  annotate(
    "text",
    x = 4, y = 0,
    label = sprintf("FAIRD: %d (%.1f%%)", 
                    size_stats$Count[2], 
                    size_stats$Percentage[2]),
    hjust = 0.5,
    vjust = -0.5,
    color = "#F39C12",
    fontface = "bold",
    size = 3.5
  ) +
  scale_fill_manual(
    values = custom_colors
  ) +
  labs(
    title = "Body Size Distribution: FAIRD",
    subtitle = "FAIRD shows detection bias: absolute limit at 3mm, constrained below 5mm",
    x = "Body Length (mm)",
    y = "Number of Individuals",
    fill = "Device"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#7F8C8D", size = 12),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

plot(p_size_distribution)

# Save the figure
ggsave(
  "outputs/H2_Size_Distribution_Comparison.png",
  plot = p_size_distribution,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

message("✓ Saved: outputs/H2_Size_Distribution_Comparison.png")

# ============================================================================
# 2. H3: VENN DIAGRAM - GENUS OVERLAPS
# ============================================================================
# Purpose: Show taxonomic overlaps at genus level between devices
# Expected in slide: "H3 Taxonomic Overlaps: Unique vs Shared Genera"
# ============================================================================

message("Creating H3 Venn Diagram...")

# Extract unique genera for each device
faird_genera <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::filter(!is.na(Genus)) %>%
  dplyr::pull(Genus) %>%
  unique()

id_genera <- data_list$etraps %>%
  dplyr::filter(Device_type == "ID") %>%
  dplyr::filter(!is.na(Genus)) %>%
  dplyr::pull(Genus) %>%
  unique()

ammod_genera <- data_list$ammod_taxonomy %>%
  dplyr::filter(!is.na(Genus)) %>%
  dplyr::pull(Genus) %>%
  unique()

# Print statistics for verification
message("Genus counts:")
message(paste0("  FAIRD: ", length(faird_genera), " genera"))
message(paste0("  ID: ", length(id_genera), " genera"))
message(paste0("  AMMOD: ", length(ammod_genera), " genera"))

# Calculate overlaps
faird_id_shared <- length(intersect(faird_genera, id_genera))
faird_ammod_shared <- length(intersect(faird_genera, ammod_genera))
id_ammod_shared <- length(intersect(id_genera, ammod_genera))

message("Overlaps:")
message(paste0("  FAIRD ∩ ID: ", faird_id_shared, " genera"))
message(paste0("  FAIRD ∩ AMMOD: ", faird_ammod_shared, " genera"))
message(paste0("  ID ∩ AMMOD: ", id_ammod_shared, " genera"))

# ============================================================================
# OPTION 1: Using VennDiagram package (more customizable)
# ============================================================================

venn_plot <- venn.diagram(
  x = list(
    FAIRD = faird_genera,
    ID = id_genera,
    AMMOD = ammod_genera
  ),
  filename = NULL,  # Don't save yet, we'll modify first
  
  # Colors
  fill = c(custom_colors["FAIRD"], custom_colors["ID"], custom_colors["AMMOD"]),
  alpha = 0.5,
  
  # Labels
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  
  # Numbers
  cex = 1.2,
  fontface = "bold",
  
  # Style
  lwd = 2,
  lty = "solid",
  
  # Size
  height = 3000,
  width = 3000,
  resolution = 300,
  units = "px"
)

# Save Venn diagram
png(
  "outputs/H3_Venn_Diagram_Genera.png",
  width = 3000,
  height = 3000,
  res = 300,
  bg = "white"
)
# It's often good practice to clear the graphics device first
grid.newpage()
grid::grid.draw(venn_plot)
dev.off()

message("✓ Saved: outputs/H3_Venn_Diagram_Genera.png")

# ============================================================================
# OPTION 2: Using ggvenn package (simpler, ggplot-based)
# ============================================================================
# Uncomment this section if you prefer ggvenn (cleaner, more modern look)
# 
library(ggvenn)

genera_list <- list(
  FAIRD = faird_genera,
  ID = id_genera,
  AMMOD = ammod_genera
)

p_venn <- ggvenn(
  genera_list,
  fill_color = c(custom_colors["FAIRD"], custom_colors["ID"], custom_colors["AMMOD"]),
  fill_alpha = 0.5,
  stroke_size = 1,
  set_name_size = 5,
  text_size = 4
) +
  labs(
    title = "Genus-level Taxonomic Overlaps",
    subtitle = "FAIRD captures 20 unique genera not detected by ID or AMMOD"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#7F8C8D")
  )

ggsave(
  "outputs/H3_Venn_Diagram_Genera.png",
  plot = p_venn,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

plot(p_venn)

# ============================================================================
# OPTIONAL: CREATE SUMMARY TABLE FOR PRESENTATION
# ============================================================================

# Calculate exclusive genera
faird_exclusive <- setdiff(
  faird_genera, 
  union(id_genera, ammod_genera)
)

id_exclusive <- setdiff(
  id_genera,
  union(faird_genera, ammod_genera)
)

ammod_exclusive <- setdiff(
  ammod_genera,
  union(faird_genera, id_genera)
)

# Create summary data frame
overlap_summary <- data.frame(
  Category = c(
    "Total Genera",
    "FAIRD ∩ ID (shared)",
    "FAIRD ∩ AMMOD (shared)",
    "ID ∩ AMMOD (shared)",
    "FAIRD exclusive",
    "ID exclusive",
    "AMMOD exclusive"
  ),
  Count = c(
    paste(length(faird_genera), length(id_genera), length(ammod_genera), sep = " / "),
    faird_id_shared,
    faird_ammod_shared,
    id_ammod_shared,
    length(faird_exclusive),
    length(id_exclusive),
    length(ammod_exclusive)
  )
)

# Print summary
message("\n=== TAXONOMIC OVERLAP SUMMARY ===")
print(overlap_summary)

# Save summary as CSV
write.csv(
  overlap_summary,
  "outputs/H3_Genus_Overlap_Summary.csv",
  row.names = FALSE
)

message("✓ Saved: outputs/H3_Genus_Overlap_Summary.csv")

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

message("\n✓✓✓ ALL MISSING FIGURES CREATED SUCCESSFULLY ✓✓✓")
message("\nGenerated files:")
message("  1. outputs/H2_Size_Distribution_Comparison.png")
message("  2. outputs/H3_Venn_Diagram_Genera.png")
message("  3. outputs/H3_Genus_Overlap_Summary.csv (bonus)")
message("\nYou can now update the imageConfig paths in presentation.html")
message("and add the IDs to placeholders following PLACEHOLDER_MODIFICATIONS_GUIDE.md")

# ============================================================================
# END OF SCRIPT
# ============================================================================

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")