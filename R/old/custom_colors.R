# ==============================================================================
# custom_colors.R
# Project-wide color palette for ggplot2 figures
#
# IMPORTANT — before using these colors in any ggplot2 script, read:
#   custom_colors_guide.md
#
# Key rule: always pass `custom_colors` directly to scale_color_manual() /
# scale_fill_manual(). Never subset it to build intermediate named vectors.
# ==============================================================================

# Load the RColorBrewer library, needed for the palette
library(RColorBrewer)

# Color mapping
custom_colors <- c(
  "Site1" = "gold",
  "Site2" = "darkorange",
  "Site3" = "greenyellow",
  "Site4" = "darkgreen",
  "Maize" = "gold",
  "Meadow" = "darkgreen",
  "FAIRD1" = "gold",
  "FAIRD2" = "darkorange",
  "FAIRD3" = "greenyellow",
  "FAIRD4" = "darkgreen",
  "ID1" = "gold",
  "ID2" = "darkorange",
  "ID3" = "greenyellow",
  "ID4" = "darkgreen",
  "AMMOD1" = "gold",
  "AMMOD2" = "darkorange",
  "AMMOD3" = "greenyellow",
  "AMMOD4" = "darkgreen",
  "FAIRD" = "blue",
  "ID" = "magenta",
  "AMMOD" = "red",
  "FAIRD:Maize" = "gold",
  "FAIRD:Meadow" = "yellowgreen",
  "ID:Maize" = "gold",
  "ID:Meadow" = "yellowgreen",
  "AMMOD:Maize" = "gold",
  "AMMOD:Meadow" = "yellowgreen",
  "AIRTEMP" = "red",
  "WINDSPEED" = "cyan",
  "PRECIPITATION" = "lightskyblue"
)

# Color map according device selection
device_colors <- c(
  "1" = "blue",       # FAIR-Device
  "2" = "magenta",    # Insect Detect
  "3" = "limegreen"   # both
)

# --- Taxonomic Color Mapping ---
# Create a color palette with 11 distinct colors from the "Set3" palette
mycols <- colorRampPalette(brewer.pal(11, "Set3"))(11)

# Assign specific colors to each Order for consistency
order_colors <- setNames(mycols, c("Dermaptera", "Hymenoptera", "Psocodea", "Diptera", 
                                   "Coleoptera", "Hemiptera", "Orthoptera", "Lepidoptera", 
                                   "Unknown", "Pterygota", "Neuroptera"))

# Add a color for "Other" in case of lumping small categories
order_colors["Other"] <- "grey50"

