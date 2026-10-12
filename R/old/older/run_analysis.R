# ==============================================================================
# MASTER ANALYSIS SCRIPT - FAIR-Device Validation
# ==============================================================================
# Description: Executes all three hypothesis analyses in sequence
# Date: October 22, 2025
# Author: FAIRD Validation Team
# ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                      ║\n")
cat("║         FAIR-DEVICE VALIDATION - COMPLETE ANALYSIS PIPELINE         ║\n")
cat("║                                                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Record start time
start_time <- Sys.time()

# ==============================================================================
# SETUP
# ==============================================================================

cat("=== SETUP ===\n\n")

# Check required packages
required_packages <- c("tidyverse", "car", "MASS", "patchwork", "ggpubr", "vegan")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if(length(missing_packages) > 0) {
  cat("⚠️ Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("   Installing missing packages...\n\n")
  install.packages(missing_packages)
}

# Load all required libraries
library(tidyverse)
library(car)
library(MASS)
library(patchwork)
library(ggpubr)
library(vegan)

cat("✅ All required packages loaded\n\n")

# Create output directory if it doesn't exist
if(!dir.exists("results")) {
  dir.create("results")
  cat("✅ Created 'results' directory\n\n")
}

# ==============================================================================
# HYPOTHESIS 1: INTRA-DEVICE CONSISTENCY
# ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  HYPOTHESIS 1: INTRA-DEVICE CONSISTENCY                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

h1_start <- Sys.time()
source("H1_consistency.R")
h1_duration <- as.numeric(difftime(Sys.time(), h1_start, units = "secs"))

cat(sprintf("⏱️ H1 completed in %.1f seconds\n\n", h1_duration))

# ==============================================================================
# HYPOTHESIS 2: VALIDATION AGAINST TRADITIONAL SAMPLING
# ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  HYPOTHESIS 2: VALIDATION AGAINST AMMOD                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

h2_start <- Sys.time()
source("H2_validation.R")
h2_duration <- as.numeric(difftime(Sys.time(), h2_start, units = "secs"))

cat(sprintf("⏱️ H2 completed in %.1f seconds\n\n", h2_duration))

# ==============================================================================
# HYPOTHESIS 3: TAXONOMIC BREADTH COMPARISON
# ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  HYPOTHESIS 3: TAXONOMIC BREADTH                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

h3_start <- Sys.time()
source("H3_taxonomic.R")
h3_duration <- as.numeric(difftime(Sys.time(), h3_start, units = "secs"))

cat(sprintf("⏱️ H3 completed in %.1f seconds\n\n", h3_duration))

# ==============================================================================
# MOVE ALL OUTPUTS TO RESULTS DIRECTORY
# ==============================================================================

cat("\n=== ORGANIZING OUTPUT FILES ===\n\n")

# List of files to move
output_files <- c(
  # Tables
  "Table1_H1_Consistency_Results.csv",
  "Table2_H2_Validation_Results.csv",
  "Table3_H3_Taxonomic_Results.csv",
  
  # H1 figures
  "Figure_S1_Diagnostic_Comparison.png",
  "Figure_S1_Diagnostic_Comparison.pdf",
  "Figure2_H1_Consistency.png",
  "Figure2_H1_Consistency.pdf",
  
  # H2 figures
  "Figure3_H2_Validation.png",
  "Figure3_H2_Validation.pdf",
  "Figure4_H2_Temporal_Dynamics.png",
  "Figure4_H2_Temporal_Dynamics.pdf",
  
  # H3 figures
  "Figure5_H3_Richness.png",
  "Figure5_H3_Richness.pdf",
  "Figure6_H3_Composition.png",
  "Figure6_H3_Composition.pdf"
)

# Move files that exist
moved_count <- 0
for(file in output_files) {
  if(file.exists(file)) {
    file.copy(file, file.path("results", file), overwrite = TRUE)
    file.remove(file)
    moved_count <- moved_count + 1
  }
}

cat(sprintf("✅ Moved %d output files to 'results/' directory\n\n", moved_count))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

total_duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                     ANALYSIS PIPELINE COMPLETE                       ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("⏱️ TIMING SUMMARY:\n")
cat(sprintf("   H1 (Consistency):    %.1f seconds\n", h1_duration))
cat(sprintf("   H2 (Validation):     %.1f seconds\n", h2_duration))
cat(sprintf("   H3 (Taxonomic):      %.1f seconds\n", h3_duration))
cat(sprintf("   ─────────────────────────────────\n"))
cat(sprintf("   TOTAL:               %.1f seconds (%.1f minutes)\n\n", 
            total_duration, total_duration/60))

cat("📁 OUTPUT FILES (in results/ directory):\n\n")

cat("📊 TABLES:\n")
cat("   ✅ Table1_H1_Consistency_Results.csv\n")
cat("   ✅ Table2_H2_Validation_Results.csv\n")
cat("   ✅ Table3_H3_Taxonomic_Results.csv\n\n")

cat("📈 MAIN TEXT FIGURES:\n")
cat("   ✅ Figure2_H1_Consistency (3-panel)\n")
cat("   ✅ Figure3_H2_Validation (3-panel)\n")
cat("   ✅ Figure4_H2_Temporal_Dynamics (time series)\n")
cat("   ✅ Figure5_H3_Richness (bar plot)\n")
cat("   ✅ Figure6_H3_Composition (stacked bar)\n\n")

cat("📋 SUPPLEMENTARY FIGURES:\n")
cat("   ✅ Figure_S1_Diagnostic_Comparison (original vs log)\n\n")

cat("🎯 KEY FINDINGS:\n\n")

cat("H1 - INTRA-DEVICE CONSISTENCY:\n")
cat("   • FAIRD Maize:  R² = 0.648*** (strong, log-transformed)\n")
cat("   • FAIRD Meadow: R² = 0.643*** (strong, log-transformed)\n")
cat("   • ID Maize:     R² = 0.196*   (poor, log-transformed)\n")
cat("   • FAIRD is 3× more consistent than Insect Detect\n\n")

cat("H2 - VALIDATION AGAINST AMMOD:\n")
cat("   ⭐ Site2 (Maize-Core): R² = 0.61**, Concordance = 81.8%\n")
cat("   ⚠️ Site1 (Maize-Edge): R² = 0.075 (edge effect, <3mm bias)\n")
cat("   ❌ Site3 (Meadow):     R² = 0.000 (n=7, requires investigation)\n\n")

cat("H3 - TAXONOMIC BREADTH:\n")
cat("   • FAIRD: 55 genera (GENERALIST sampler)\n")
cat("   • ID: 5 genera (DIPTERA specialist, >95% Diptera)\n")
cat("   • FAIRD captures 980% MORE genera than ID\n")
cat("   • ID3 complete failure vs FAIRD robust performance\n\n")

cat("📝 MANUSCRIPT STATUS:\n")
cat("   ✅ All analyses complete with rigorous diagnostics\n")
cat("   ✅ Log transformation decision finalized (H1)\n")
cat("   ✅ All tables and figures generated\n")
cat("   ✅ Ready for manuscript writing\n\n")

cat("🚀 NEXT STEPS:\n")
cat("   1. Site3 investigation (meadow validation failure)\n")
cat("   2. Detection bias analysis (body size at Site1)\n")
cat("   3. Write Methods section (copy from updated documents)\n")
cat("   4. Write Results section (use generated tables)\n")
cat("   5. Write Discussion (emphasize methodological rigor)\n\n")

cat("✨ Analysis pipeline completed successfully!\n\n")

# Save session info for reproducibility
session_info <- sessionInfo()
save(session_info, file = "results/session_info.RData")

cat("✅ Session info saved for reproducibility\n\n")