# ==============================================================================
# H3: TAXONOMIC SIMILARITY ANALYSIS - FAIRD vs AMMOD (GENUS-LEVEL)
# ==============================================================================
# Description: Daily taxonomic similarity between FAIRD and AMMOD at GENUS level
#              FAIRD: ~43% identified to genus (real resolution)
#              AMMOD: Species collapsed to genus for fair comparison
# Date: November 26, 2025
# Author: Juan Chiavassa & Claude
# ==============================================================================
# Methodology: Based on Gemini Deep Research findings
#   - Taxonomic Similarity (ΔS) - Izsak & Price (2001)
#   - Beta diversity decomposition - Baselga framework
#   - Taxonomic distance matrices - Clarke & Warwick
#   - Analysis level: GENUS (reflects FAIRD's real taxonomic resolution)
# ==============================================================================

# --- Setup ---
library(tidyverse)
library(vegan)      # For taxa2dist(), taxondive()
library(betapart)   # For beta.pair()
library(patchwork)  # For combining plots
library(lubridate)

# ==============================================================================
# CONFIGURATION: SIZE THRESHOLD
# ==============================================================================
# Set minimum body length threshold for size-corrected analysis (in mm)
# Change this value to test different thresholds (e.g., 3, 4, 5 mm)
SIZE_THRESHOLD <- 5 # mm

cat("\n⚙️  SIZE THRESHOLD SET TO:", SIZE_THRESHOLD, "mm\n")
cat("   (All outputs will use this threshold)\n\n")

# Load data
cat("\n========================================\n")
cat("H3: TAXONOMIC SIMILARITY FAIRD vs AMMOD\n")
cat("========================================\n\n")

source("load_data.R")
load("data_cleaned.RData")

# Load custom colors
source("custom_colors.R")

cat("✅ Colors: FAIRD=", custom_colors["FAIRD"], 
    ", AMMOD=", custom_colors["AMMOD"], "\n\n")

# ==============================================================================
# SECTION 1: DATA PREPARATION
# ==============================================================================

cat("\n=== SECTION 1: DATA PREPARATION ===\n\n")

# Define temporal scope: ALL 22 days (both AMMOD bottle batches)
date_start <- as.Date("2023-08-23")
date_end <- as.Date("2023-09-13")  # CHANGED: was 2023-09-03

cat("Temporal scope: ", as.character(date_start), "to", as.character(date_end), "\n")
cat("Duration: 22 days (Batch 1: Aug 23-Sep 3 + Sep 4 excluded + Batch 2: Sep 5-13)\n\n")

# --- 1.1: Extract FAIRD data (individual-level with variable resolution) ---
cat("Loading FAIRD data...\n")

faird_data <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::filter(Date >= date_start & Date <= date_end) %>%
  dplyr::select(Date, Device, Site, Class, Order, Suborder, Infraorder, 
                Superfamily, Family, Subfamily, Genus, Species)

cat("\n✅ FAIRD data loaded successfully\n")
cat("   Total individuals:", nrow(faird_data), "\n")
cat("   Date range:", as.character(min(faird_data$Date)), "to", 
    as.character(max(faird_data$Date)), "\n")
cat("   Devices:", paste(sort(unique(faird_data$Device)), collapse=", "), "\n")
cat("   Sites:", paste(sort(unique(faird_data$Site)), collapse=", "), "\n")
cat("   Days with data:", length(unique(faird_data$Date)), "\n")

# Data quality checks
cat("\n📊 FAIRD data quality checks:\n")
cat("   Records with Class:", sum(!is.na(faird_data$Class) & faird_data$Class != ""), "\n")
cat("   Records with Order:", sum(!is.na(faird_data$Order) & faird_data$Order != ""), "\n")
cat("   Records with Family:", sum(!is.na(faird_data$Family) & faird_data$Family != ""), "\n")
cat("   Records with Genus:", sum(!is.na(faird_data$Genus) & faird_data$Genus != ""), "\n")
cat("   Records with Species:", sum(!is.na(faird_data$Species) & faird_data$Species != ""), "\n\n")

# --- 1.2: Extract AMMOD data (OTU-level, species resolution) ---
cat("Loading AMMOD data...\n")

ammod_data <- data_list$ammod_taxonomy %>%
  dplyr::filter(Date >= date_start & Date <= date_end) %>%
  dplyr::select(Date, Device, Site, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus, Species, Body_length, Body_width)
# Note: Count column not included - AMMOD taxonomy is presence/absence only
# Body_length & Body_width needed for size-corrected analysis

# --- CRITICAL: Exclude AMMOD4 mechanical failure period ---
# AMMOD4 bottle rotation system jammed Aug 28-Sep 3 (bottle 6 accumulated 7 days)
# These days must be excluded as they represent accumulated samples, not daily captures
ammod_data <- ammod_data %>%
  dplyr::filter(
    !(Device == "AMMOD4" & Date >= as.Date("2023-08-28") & Date <= as.Date("2023-09-03"))
  )

cat("⚠️  AMMOD4 MECHANICAL FAILURE EXCLUSION:\n")
cat("   Excluded AMMOD4 days: 2023-08-28 to 2023-09-03 (7 days)\n")
cat("   Reason: Bottle rotation jammed - accumulated sample not representative\n")
cat("   AMMOD4 available days: 2023-08-23 to 2023-08-27 (5 days only)\n\n")

cat("\n✅ AMMOD data loaded successfully\n")
cat("   Total OTUs:", nrow(ammod_data), "\n")
cat("   Date range:", as.character(min(ammod_data$Date)), "to",
    as.character(max(ammod_data$Date)), "\n")
cat("   Devices:", paste(sort(unique(ammod_data$Device)), collapse=", "), "\n")
cat("   Sites:", paste(sort(unique(ammod_data$Site)), collapse=", "), "\n")
cat("   Days with data:", length(unique(ammod_data$Date)), "\n")

# Data quality checks
cat("\n📊 AMMOD data quality checks:\n")
cat("   OTUs with Class:", sum(!is.na(ammod_data$Class) & ammod_data$Class != ""), "\n")
cat("   OTUs with Order:", sum(!is.na(ammod_data$Order) & ammod_data$Order != ""), "\n")
cat("   OTUs with Family:", sum(!is.na(ammod_data$Family) & ammod_data$Family != ""), "\n")
cat("   OTUs with Genus:", sum(!is.na(ammod_data$Genus) & ammod_data$Genus != ""), "\n")
cat("   OTUs with Species:", sum(!is.na(ammod_data$Species) & ammod_data$Species != ""), "\n")
cat("   Note: Data is presence/absence (Count column not used)\n\n")

# --- 1.3: Create taxonomic ID for each observation ---
# CRITICAL: Use GENUS level as comparison standard
# Rationale: FAIRD's real taxonomic resolution is predominantly at genus level
#            (only ~0.1% identified to species, ~43% to genus)
#            AMMOD species will be collapsed to genus level for fair comparison
#
# IMPORTANT: "#N/C" (No Clasificado) is treated as NA/missing, not valid taxonomy

create_taxon_id <- function(data) {
  data %>%
    dplyr::mutate(
      # Create unique taxon identifier at GENUS level (or higher if genus unavailable)
      # Check: not NA, not empty string, and not "#N/C"
      Taxon_ID = dplyr::case_when(
        !is.na(Genus) & Genus != "" & Genus != "#N/C" ~ Genus,  # Primary: Use Genus
        !is.na(Family) & Family != "" & Family != "#N/C" ~ Family,
        !is.na(Superfamily) & Superfamily != "" & Superfamily != "#N/C" ~ Superfamily,
        !is.na(Infraorder) & Infraorder != "" & Infraorder != "#N/C" ~ Infraorder,
        !is.na(Suborder) & Suborder != "" & Suborder != "#N/C" ~ Suborder,
        !is.na(Order) & Order != "" & Order != "#N/C" ~ Order,
        !is.na(Class) & Class != "" & Class != "#N/C" ~ Class,
        TRUE ~ "Unknown"
      ),
      
      # Store which taxonomic level was used
      Taxon_Level = dplyr::case_when(
        !is.na(Genus) & Genus != "" & Genus != "#N/C" ~ "Genus",
        !is.na(Family) & Family != "" & Family != "#N/C" ~ "Family",
        !is.na(Superfamily) & Superfamily != "" & Superfamily != "#N/C" ~ "Superfamily",
        !is.na(Infraorder) & Infraorder != "" & Infraorder != "#N/C" ~ "Infraorder",
        !is.na(Suborder) & Suborder != "" & Suborder != "#N/C" ~ "Suborder",
        !is.na(Order) & Order != "" & Order != "#N/C" ~ "Order",
        !is.na(Class) & Class != "" & Class != "#N/C" ~ "Class",
        TRUE ~ "Unknown"
      )
    )
}

faird_data <- create_taxon_id(faird_data)
ammod_data <- create_taxon_id(ammod_data)

cat("✅ Taxon IDs created successfully\n\n")

# --- 1.4: Create size-corrected AMMOD dataset (≥SIZE_THRESHOLD only) ---
cat("🔧 Creating size-corrected AMMOD dataset...\n")
cat("   Filtering AMMOD by Body_length >=", SIZE_THRESHOLD, "mm (FAIRD's detection threshold)\n\n")

# Check if Body_length column exists
if (!"Body_length" %in% names(ammod_data)) {
  cat("⚠️  WARNING: Body_length column not found in AMMOD data\n")
  cat("   Size-corrected analysis will be skipped\n\n")
  run_size_correction <- FALSE
  run_baseline_comparison <- FALSE
} else {
  # IMPORTANT: Create full dataset copy BEFORE any filtering
  ammod_data_full <- ammod_data
  
  # Filter AMMOD by size threshold
  ammod_data_size_corrected <- ammod_data %>%
    dplyr::filter(!is.na(Body_length) & Body_length >= SIZE_THRESHOLD)
  
  n_original <- nrow(ammod_data_full)
  n_filtered <- nrow(ammod_data_size_corrected)
  n_removed <- n_original - n_filtered
  pct_retained <- round(n_filtered / n_original * 100, 1)
  
  cat("📊 Size filtering results:\n")
  cat("   Original AMMOD OTUs:", n_original, "\n")
  cat("   OTUs ≥", SIZE_THRESHOLD, "mm:", n_filtered, "(", pct_retained, "%)\n")
  cat("   OTUs <", SIZE_THRESHOLD, "mm (excluded):", n_removed, "(", round(100-pct_retained, 1), "%)\n\n")
  
  # Check FAIRD size distribution for comparison
  if ("Body_length" %in% names(faird_data)) {
    faird_below_threshold <- sum(faird_data$Body_length < SIZE_THRESHOLD, na.rm = TRUE)
    faird_total <- sum(!is.na(faird_data$Body_length))
    cat("   FAIRD <", SIZE_THRESHOLD, "mm (for reference):", faird_below_threshold, "of", faird_total,
        "(", round(faird_below_threshold/faird_total*100, 2), "%)\n\n")
  }
  
  run_size_correction <- TRUE
  
  # If threshold is different from 3mm, also create 3mm baseline for comparison
  if (SIZE_THRESHOLD != 3) {
    cat("   📌 Additionally creating 3mm baseline for comparison...\n")
    ammod_data_size_corrected_3mm <- ammod_data_full %>%
      dplyr::filter(!is.na(Body_length) & Body_length >= 3.0)
    n_3mm <- nrow(ammod_data_size_corrected_3mm)
    cat("      OTUs ≥3mm (baseline):", n_3mm, "(", round(n_3mm/n_original*100, 1), "%)\n")
    cat("      ✅ ammod_data_size_corrected_3mm created with", n_3mm, "rows\n\n")
    run_baseline_comparison <- TRUE
    cat("      ✅ run_baseline_comparison set to TRUE\n\n")
  } else {
    run_baseline_comparison <- FALSE
    cat("      ℹ️  SIZE_THRESHOLD = 3, no baseline comparison needed\n\n")
  }
}

# NOW filter Unknown from all datasets
cat("🔍 Checking for '#N/C' (No Clasificado) entries:\n")

# Count "Unknown" in final Taxon_IDs
faird_unknown <- sum(faird_data$Taxon_ID == "Unknown", na.rm = TRUE)
ammod_unknown <- sum(ammod_data$Taxon_ID == "Unknown", na.rm = TRUE)

cat("   FAIRD: ", faird_unknown, "individuals classified as 'Unknown'\n")
cat("          (", round(faird_unknown/nrow(faird_data)*100, 1), "% of total)\n")
cat("   AMMOD: ", ammod_unknown, "OTUs classified as 'Unknown'\n")
cat("          (", round(ammod_unknown/nrow(ammod_data)*100, 1), "% of total)\n")

if (faird_unknown > 0 | ammod_unknown > 0) {
  cat("\n   ⚠️  'Unknown' taxa will be excluded from similarity calculations\n")
  cat("      (no valid taxonomic classification at any level)\n")
  
  # Filter out Unknown taxa from ALL datasets
  faird_data <- faird_data %>% dplyr::filter(Taxon_ID != "Unknown")
  ammod_data <- ammod_data %>% dplyr::filter(Taxon_ID != "Unknown")
  
  # Also filter from size-corrected datasets if they exist
  if (exists("ammod_data_full")) {
    ammod_data_full <- ammod_data_full %>% dplyr::filter(Taxon_ID != "Unknown")
  }
  if (exists("ammod_data_size_corrected")) {
    ammod_data_size_corrected <- ammod_data_size_corrected %>% dplyr::filter(Taxon_ID != "Unknown")
  }
  if (exists("ammod_data_size_corrected_3mm")) {
    ammod_data_size_corrected_3mm <- ammod_data_size_corrected_3mm %>% dplyr::filter(Taxon_ID != "Unknown")
  }
  
  cat("   ✅ Filtered data: FAIRD =", nrow(faird_data), "individuals | AMMOD =",
      nrow(ammod_data), "OTUs\n")
  
  # Verify size-corrected datasets still exist after filtering
  if (exists("ammod_data_full")) {
    cat("      ammod_data_full:", nrow(ammod_data_full), "rows (after Unknown filter)\n")
  }
  if (exists("ammod_data_size_corrected")) {
    cat("      ammod_data_size_corrected:", nrow(ammod_data_size_corrected), "rows (after Unknown filter)\n")
  }
  if (exists("ammod_data_size_corrected_3mm")) {
    cat("      ammod_data_size_corrected_3mm:", nrow(ammod_data_size_corrected_3mm), "rows (after Unknown filter)\n")
  }
  cat("\n")
} else {
  cat("\n   ✅ All taxa have valid taxonomic classification\n\n")
}

# Summary of taxonomic resolution
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📊 FAIRD TAXONOMIC RESOLUTION DISTRIBUTION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
faird_resolution <- table(faird_data$Taxon_Level)
print(faird_resolution)
cat("\nPercentages:\n")
print(round(prop.table(faird_resolution) * 100, 1))
cat("\nMost specific: ", names(faird_resolution)[which.max(faird_resolution)], 
    " (", max(faird_resolution), " individuals)\n", sep="")
cat("Unique taxa identified:", length(unique(faird_data$Taxon_ID)), "\n\n")

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📊 AMMOD TAXONOMIC RESOLUTION DISTRIBUTION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
ammod_resolution <- table(ammod_data$Taxon_Level)
print(ammod_resolution)
cat("\nPercentages:\n")
print(round(prop.table(ammod_resolution) * 100, 1))
cat("\nMost specific: ", names(ammod_resolution)[which.max(ammod_resolution)],
    " (", max(ammod_resolution), " OTUs)\n", sep="")
cat("Unique taxa identified:", length(unique(ammod_data$Taxon_ID)), "\n\n")

# ==============================================================================
# SECTION 2: BUILD TAXONOMIC DISTANCE MATRIX
# ==============================================================================

cat("\n=== SECTION 2: TAXONOMIC DISTANCE MATRIX ===\n\n")

# --- 2.1: Create combined taxonomy table ---
# We need all unique taxa from both FAIRD and AMMOD to build distance matrix

# Extract unique taxa from FAIRD
faird_taxa <- faird_data %>%
  dplyr::select(Taxon_ID, Class, Order, Suborder, Infraorder, Superfamily, 
                Family, Subfamily, Genus, Species) %>%
  dplyr::distinct()

# Extract unique taxa from AMMOD
ammod_taxa <- ammod_data %>%
  dplyr::select(Taxon_ID, Class, Order, Suborder, Infraorder, Superfamily,
                Family, Subfamily, Genus, Species) %>%
  dplyr::distinct()

# Combine and create master taxonomy table
all_taxa <- dplyr::bind_rows(
  faird_taxa %>% dplyr::mutate(Source = "FAIRD"),
  ammod_taxa %>% dplyr::mutate(Source = "AMMOD")
) %>%
  dplyr::distinct(Taxon_ID, .keep_all = TRUE)

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📊 COMBINED TAXONOMY TABLE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Total unique taxa across both methods:", nrow(all_taxa), "\n")
cat("  FAIRD unique taxa:", nrow(faird_taxa), "\n")
cat("  AMMOD unique taxa:", nrow(ammod_taxa), "\n")
cat("  Shared taxa (exact match):", nrow(faird_taxa) + nrow(ammod_taxa) - nrow(all_taxa), "\n")

shared_exact <- intersect(faird_taxa$Taxon_ID, ammod_taxa$Taxon_ID)
cat("  Shared % (of FAIRD):", round(length(shared_exact) / nrow(faird_taxa) * 100, 1), "%\n")
cat("  Shared % (of AMMOD):", round(length(shared_exact) / nrow(ammod_taxa) * 100, 1), "%\n")

cat("\nTaxonomic composition of combined taxa:\n")
cat("  Classes:", length(unique(all_taxa$Class[!is.na(all_taxa$Class) & all_taxa$Class != ""])), "\n")
cat("  Orders:", length(unique(all_taxa$Order[!is.na(all_taxa$Order) & all_taxa$Order != ""])), "\n")
cat("  Families:", length(unique(all_taxa$Family[!is.na(all_taxa$Family) & all_taxa$Family != ""])), "\n")
cat("  Genera:", length(unique(all_taxa$Genus[!is.na(all_taxa$Genus) & all_taxa$Genus != ""])), "\n\n")

# --- 2.2: Build taxonomic distance matrix using vegan::taxa2dist() ---

# Prepare taxonomy table in format required by taxa2dist()
# Must have rownames as taxon IDs and columns as hierarchical levels
taxonomy_table <- all_taxa %>%
  dplyr::select(Taxon_ID, Class, Order, Suborder, Infraorder, Superfamily,
                Family, Subfamily, Genus, Species) %>%
  as.data.frame()

rownames(taxonomy_table) <- taxonomy_table$Taxon_ID
taxonomy_table <- taxonomy_table %>% dplyr::select(-Taxon_ID)

# Create distance matrix with variable step lengths
# varstep = TRUE allows different weights for different hierarchical levels
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("🔧 CREATING TAXONOMIC DISTANCE MATRIX\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Using vegan::taxa2dist() with varstep=TRUE\n")
cat("Processing", nrow(taxonomy_table), "taxa...\n\n")

taxon_distances <- vegan::taxa2dist(taxonomy_table, varstep = TRUE, check = TRUE)

cat("✅ Distance matrix created successfully\n\n")
cat("📊 Distance matrix properties:\n")
cat("   Dimensions:", attr(taxon_distances, "Size"), "×", attr(taxon_distances, "Size"), "\n")
cat("   Total pairwise distances:", length(taxon_distances), "\n")
cat("   Min distance:", round(min(taxon_distances), 4), "\n")
cat("   Max distance:", round(max(taxon_distances), 4), "\n")
cat("   Mean distance:", round(mean(taxon_distances), 4), "\n")
cat("   Median distance:", round(median(taxon_distances), 4), "\n")
cat("   SD distance:", round(sd(taxon_distances), 4), "\n\n")

# Check for any issues
cat("🔍 Validation checks:\n")
if (any(is.na(taxon_distances))) {
  cat("   ⚠️  WARNING: NA values found in distance matrix!\n")
  cat("      Count:", sum(is.na(taxon_distances)), "\n")
} else {
  cat("   ✅ No NA values\n")
}

if (any(taxon_distances < 0)) {
  cat("   ⚠️  WARNING: Negative distances found! (Should not occur)\n")
} else {
  cat("   ✅ All distances non-negative\n")
}

if (max(taxon_distances) == 0) {
  cat("   ⚠️  WARNING: Max distance is zero! (All taxa identical?)\n")
} else {
  cat("   ✅ Distance range > 0\n")
}

cat("\n")

# Convert to matrix for easier indexing
dist_matrix <- as.matrix(taxon_distances)

# ==============================================================================
# SECTION 3: CALCULATE DAILY TAXONOMIC SIMILARITY (ΔS)
# ==============================================================================

cat("\n=== SECTION 3: DAILY TAXONOMIC SIMILARITY ===\n\n")

# --- 3.1: Custom implementation of Izsak-Price Taxonomic Similarity (ΔS) ---

calculate_taxonomic_similarity <- function(taxa_A, taxa_B, dist_matrix) {
  # Calculate taxonomic similarity between two taxon lists
  # Based on Izsak & Price (2001) - average minimum distance
  # 
  # Args:
  #   taxa_A: vector of taxon IDs from sample A
  #   taxa_B: vector of taxon IDs from sample B
  #   dist_matrix: matrix of pairwise taxonomic distances
  # 
  # Returns:
  #   ΔS: similarity value (0 = completely different, 1 = identical)
  
  # Remove any taxa not in distance matrix
  taxa_A <- taxa_A[taxa_A %in% rownames(dist_matrix)]
  taxa_B <- taxa_B[taxa_B %in% rownames(dist_matrix)]
  
  # Handle edge cases
  if (length(taxa_A) == 0 | length(taxa_B) == 0) {
    return(NA)
  }
  
  # For each taxon in A, find minimum distance to any taxon in B
  min_distances_A_to_B <- sapply(taxa_A, function(taxon_a) {
    if (taxon_a %in% taxa_B) {
      return(0)  # Exact match
    } else {
      distances_to_B <- dist_matrix[taxon_a, taxa_B]
      return(min(distances_to_B, na.rm = TRUE))
    }
  })
  
  # For each taxon in B, find minimum distance to any taxon in A
  min_distances_B_to_A <- sapply(taxa_B, function(taxon_b) {
    if (taxon_b %in% taxa_A) {
      return(0)  # Exact match
    } else {
      distances_to_A <- dist_matrix[taxon_b, taxa_A]
      return(min(distances_to_A, na.rm = TRUE))
    }
  })
  
  # Average minimum distance (bidirectional)
  avg_min_distance <- mean(c(min_distances_A_to_B, min_distances_B_to_A))
  
  # Maximum possible distance in the matrix
  max_distance <- max(dist_matrix, na.rm = TRUE)
  
  # Calculate similarity (inverse of normalized distance)
  delta_S <- 1 - (avg_min_distance / max_distance)
  
  return(delta_S)
}

# --- 3.2: Calculate similarity for each site and day ---

sites <- c("Site1", "Site2", "Site3", "Site4")
all_dates <- seq.Date(date_start, date_end, by = "day")

# Initialize results dataframe
similarity_results <- expand.grid(
  Date = all_dates,
  Site = sites,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    Delta_S_raw = NA_real_,
    Delta_S_size_corrected = NA_real_,
    Delta_S_size_corrected_3mm = NA_real_,
    N_taxa_FAIRD = NA_integer_,
    N_taxa_AMMOD_full = NA_integer_,
    N_taxa_AMMOD_filtered = NA_integer_,
    N_taxa_AMMOD_3mm = NA_integer_,
    N_shared_exact_raw = NA_integer_,
    N_shared_exact_size_corrected = NA_integer_,
    N_shared_exact_3mm = NA_integer_
  )

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("🔄 CALCULATING DAILY TAXONOMIC SIMILARITY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Total calculations:", nrow(similarity_results), "(", length(all_dates), "days ×", 
    length(sites), "sites )\n")

if (run_size_correction) {
  cat("\n   Computing TWO similarity metrics:\n")
  cat("   1. RAW: FAIRD vs all AMMOD taxa (original analysis)\n")
  cat("   2. SIZE-CORRECTED: FAIRD vs AMMOD ≥", SIZE_THRESHOLD, "mm (controlled comparison)\n")
}
cat("\n")

progress_counter <- 0
total_calcs <- nrow(similarity_results)

for (i in 1:nrow(similarity_results)) {
  current_date <- similarity_results$Date[i]
  current_site <- similarity_results$Site[i]
  
  # Progress indicator every 10 calculations
  progress_counter <- progress_counter + 1
  if (progress_counter %% 10 == 0 || progress_counter == total_calcs) {
    cat("   Progress:", progress_counter, "/", total_calcs, 
        "(", round(progress_counter/total_calcs*100, 0), "%)\n")
  }
  
  # Get FAIRD taxa for this site and date
  faird_taxa_day <- faird_data %>%
    dplyr::filter(Date == current_date, Site == current_site) %>%
    dplyr::pull(Taxon_ID) %>%
    unique()
  
  # Get AMMOD taxa (FULL dataset - original analysis)
  ammod_taxa_day_full <- ammod_data_full %>%
    dplyr::filter(Date == current_date, Site == current_site) %>%
    dplyr::pull(Taxon_ID) %>%
    unique()
  
  # Store counts
  similarity_results$N_taxa_FAIRD[i] <- length(faird_taxa_day)
  similarity_results$N_taxa_AMMOD_full[i] <- length(ammod_taxa_day_full)
  similarity_results$N_shared_exact_raw[i] <- length(intersect(faird_taxa_day, ammod_taxa_day_full))
  
  # Calculate RAW similarity (original analysis)
  if (length(faird_taxa_day) > 0 && length(ammod_taxa_day_full) > 0) {
    similarity_results$Delta_S_raw[i] <- calculate_taxonomic_similarity(
      faird_taxa_day, ammod_taxa_day_full, dist_matrix
    )
  }
  
  # Calculate SIZE-CORRECTED similarity if available
  if (run_size_correction) {
    # Get AMMOD taxa (FILTERED ≥SIZE_THRESHOLD)
    ammod_taxa_day_filtered <- ammod_data_size_corrected %>%
      dplyr::filter(Date == current_date, Site == current_site) %>%
      dplyr::pull(Taxon_ID) %>%
      unique()
    
    similarity_results$N_taxa_AMMOD_filtered[i] <- length(ammod_taxa_day_filtered)
    similarity_results$N_shared_exact_size_corrected[i] <- length(intersect(faird_taxa_day, ammod_taxa_day_filtered))
    
    if (length(faird_taxa_day) > 0 && length(ammod_taxa_day_filtered) > 0) {
      similarity_results$Delta_S_size_corrected[i] <- calculate_taxonomic_similarity(
        faird_taxa_day, ammod_taxa_day_filtered, dist_matrix
      )
    }
    
    # Also calculate 3mm baseline if different threshold is used
    if (run_baseline_comparison) {
      ammod_taxa_day_3mm <- ammod_data_size_corrected_3mm %>%
        dplyr::filter(Date == current_date, Site == current_site) %>%
        dplyr::pull(Taxon_ID) %>%
        unique()
      
      similarity_results$N_taxa_AMMOD_3mm[i] <- length(ammod_taxa_day_3mm)
      similarity_results$N_shared_exact_3mm[i] <- length(intersect(faird_taxa_day, ammod_taxa_day_3mm))
      
      if (length(faird_taxa_day) > 0 && length(ammod_taxa_day_3mm) > 0) {
        similarity_results$Delta_S_size_corrected_3mm[i] <- calculate_taxonomic_similarity(
          faird_taxa_day, ammod_taxa_day_3mm, dist_matrix
        )
      }
    }
  }
}

cat("\n✅ Similarity calculations complete\n\n")

# Quick summary of results
valid_results_raw <- sum(!is.na(similarity_results$Delta_S_raw))
cat("📊 Calculation summary:\n")
cat("   RAW similarity - Valid results:", valid_results_raw, "/", nrow(similarity_results), "\n")

if (run_size_correction) {
  valid_results_filtered <- sum(!is.na(similarity_results$Delta_S_size_corrected))
  cat("   SIZE-CORRECTED similarity - Valid results:", valid_results_filtered, "/", 
      nrow(similarity_results), "\n")
}

if (run_baseline_comparison) {
  valid_results_3mm <- sum(!is.na(similarity_results$Delta_S_size_corrected_3mm))
  cat("   3MM BASELINE similarity - Valid results:", valid_results_3mm, "/",
      nrow(similarity_results), "\n")
}

cat("   Missing (NA):", sum(is.na(similarity_results$Delta_S_raw)), "\n")
if (sum(is.na(similarity_results$Delta_S_raw)) > 0) {
  cat("   (NAs occur when one or both methods detected zero taxa)\n")
}
cat("\n")

# --- 3.3: Summary statistics by site ---

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📊 SUMMARY STATISTICS BY SITE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Summary for RAW similarity
summary_by_site_raw <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw)) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    N_days = dplyr::n(),
    Mean_Similarity_raw = mean(Delta_S_raw, na.rm = TRUE),
    SD_Similarity_raw = sd(Delta_S_raw, na.rm = TRUE),
    Min_Similarity_raw = min(Delta_S_raw, na.rm = TRUE),
    Max_Similarity_raw = max(Delta_S_raw, na.rm = TRUE),
    Mean_Taxa_FAIRD = mean(N_taxa_FAIRD, na.rm = TRUE),
    Mean_Taxa_AMMOD_full = mean(N_taxa_AMMOD_full, na.rm = TRUE),
    Mean_Shared_raw = mean(N_shared_exact_raw, na.rm = TRUE),
    .groups = "drop"
  )

# Summary for SIZE-CORRECTED similarity if available
if (run_size_correction) {
  summary_by_site_filtered <- similarity_results %>%
    dplyr::filter(!is.na(Delta_S_size_corrected)) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      Mean_Similarity_filtered = mean(Delta_S_size_corrected, na.rm = TRUE),
      SD_Similarity_filtered = sd(Delta_S_size_corrected, na.rm = TRUE),
      Min_Similarity_filtered = min(Delta_S_size_corrected, na.rm = TRUE),
      Max_Similarity_filtered = max(Delta_S_size_corrected, na.rm = TRUE),
      Mean_Taxa_AMMOD_filtered = mean(N_taxa_AMMOD_filtered, na.rm = TRUE),
      Mean_Shared_filtered = mean(N_shared_exact_size_corrected, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Combine both summaries
  summary_by_site <- summary_by_site_raw %>%
    dplyr::left_join(summary_by_site_filtered, by = "Site")
} else {
  summary_by_site <- summary_by_site_raw
}

print(summary_by_site)
cat("\n")

# Detailed breakdown by site
for (site in sites) {
  site_data <- summary_by_site %>% dplyr::filter(Site == site)
  
  if (nrow(site_data) > 0) {
    cat("▶", site, "\n")
    cat("   Days analyzed:", site_data$N_days, "/", length(all_dates), "\n\n")
    
    # RAW similarity
    cat("   RAW SIMILARITY (vs all AMMOD taxa):\n")
    cat("     Mean:", round(site_data$Mean_Similarity_raw * 100, 1), "% ±",
        round(site_data$SD_Similarity_raw * 100, 1), "%\n")
    cat("     Range:", round(site_data$Min_Similarity_raw * 100, 1), "% -",
        round(site_data$Max_Similarity_raw * 100, 1), "%\n")
    
    # Days below 50% similarity
    days_below_50 <- similarity_results %>%
      dplyr::filter(Site == site, !is.na(Delta_S_raw), Delta_S_raw < 0.5) %>%
      dplyr::arrange(Date) %>%
      dplyr::select(Date, Delta_S_raw)
    
    if (nrow(days_below_50) > 0) {
      cat("     Days below 50%:", nrow(days_below_50), "days\n")
      for (j in 1:nrow(days_below_50)) {
        cat("       ", as.character(days_below_50$Date[j]), ": ",
            round(days_below_50$Delta_S_raw[j] * 100, 1), "%\n", sep = "")
      }
    } else {
      cat("     Days below 50%: None (all days ≥50%)\n")
    }
    
    cat("     AMMOD taxa/day:", round(site_data$Mean_Taxa_AMMOD_full, 1), "\n")
    cat("     Exact matches/day:", round(site_data$Mean_Shared_raw, 1), 
        "(", round(site_data$Mean_Shared_raw / site_data$Mean_Taxa_FAIRD * 100, 1), "% of FAIRD)\n")
    
    # SIZE-CORRECTED similarity (if available)
    if (run_size_correction && !is.na(site_data$Mean_Similarity_filtered)) {
      cat("\n   SIZE-CORRECTED (vs AMMOD ≥", SIZE_THRESHOLD, "mm):\n")
      cat("     Mean:", round(site_data$Mean_Similarity_filtered * 100, 1), "% ±",
          round(site_data$SD_Similarity_filtered * 100, 1), "%\n")
      cat("     Range:", round(site_data$Min_Similarity_filtered * 100, 1), "% -",
          round(site_data$Max_Similarity_filtered * 100, 1), "%\n")
      
      # Days below 50% similarity (size-corrected)
      days_below_50_filtered <- similarity_results %>%
        dplyr::filter(Site == site, !is.na(Delta_S_size_corrected), Delta_S_size_corrected < 0.5) %>%
        dplyr::arrange(Date) %>%
        dplyr::select(Date, Delta_S_size_corrected)
      
      if (nrow(days_below_50_filtered) > 0) {
        cat("     Days below 50%:", nrow(days_below_50_filtered), "days\n")
        for (j in 1:nrow(days_below_50_filtered)) {
          cat("       ", as.character(days_below_50_filtered$Date[j]), ": ",
              round(days_below_50_filtered$Delta_S_size_corrected[j] * 100, 1), "%\n", sep = "")
        }
      } else {
        cat("     Days below 50%: None (all days ≥50%)\n")
      }
      
      cat("     AMMOD taxa/day:", round(site_data$Mean_Taxa_AMMOD_filtered, 1), "\n")
      cat("     Exact matches/day:", round(site_data$Mean_Shared_filtered, 1),
          "(", round(site_data$Mean_Shared_filtered / site_data$Mean_Taxa_FAIRD * 100, 1), "% of FAIRD)\n")
      
      # Calculate improvement
      improvement <- (site_data$Mean_Similarity_filtered - site_data$Mean_Similarity_raw) * 100
      cat("\n     ✨ IMPROVEMENT:", sprintf("%+.1f", improvement), "percentage points\n")
    }
    
    cat("\n   FAIRD taxa/day:", round(site_data$Mean_Taxa_FAIRD, 1), "\n\n")
  } else {
    cat("▶", site, ": No valid data\n\n")
  }
}

# Overall statistics
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📊 OVERALL STATISTICS (ALL SITES COMBINED)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# RAW similarity overall
overall_stats_raw <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw)) %>%
  dplyr::summarise(
    N_observations = dplyr::n(),
    Mean_Similarity = mean(Delta_S_raw, na.rm = TRUE),
    SD_Similarity = sd(Delta_S_raw, na.rm = TRUE),
    Median_Similarity = median(Delta_S_raw, na.rm = TRUE),
    Q25_Similarity = quantile(Delta_S_raw, 0.25, na.rm = TRUE),
    Q75_Similarity = quantile(Delta_S_raw, 0.75, na.rm = TRUE)
  )

cat("RAW SIMILARITY (vs all AMMOD taxa):\n")
cat("  Total observations:", overall_stats_raw$N_observations, "\n")
cat("  Mean similarity:", round(overall_stats_raw$Mean_Similarity * 100, 1), "%\n")
cat("  SD similarity:", round(overall_stats_raw$SD_Similarity * 100, 1), "%\n")
cat("  Median similarity:", round(overall_stats_raw$Median_Similarity * 100, 1), "%\n")
cat("  IQR:", round(overall_stats_raw$Q25_Similarity * 100, 1), "% -",
    round(overall_stats_raw$Q75_Similarity * 100, 1), "%\n")

# SIZE-CORRECTED similarity overall (if available)
if (run_size_correction) {
  overall_stats_filtered <- similarity_results %>%
    dplyr::filter(!is.na(Delta_S_size_corrected)) %>%
    dplyr::summarise(
      N_observations = dplyr::n(),
      Mean_Similarity = mean(Delta_S_size_corrected, na.rm = TRUE),
      SD_Similarity = sd(Delta_S_size_corrected, na.rm = TRUE),
      Median_Similarity = median(Delta_S_size_corrected, na.rm = TRUE),
      Q25_Similarity = quantile(Delta_S_size_corrected, 0.25, na.rm = TRUE),
      Q75_Similarity = quantile(Delta_S_size_corrected, 0.75, na.rm = TRUE)
    )
  
  cat("\nSIZE-CORRECTED SIMILARITY (vs AMMOD ≥", SIZE_THRESHOLD, "mm):\n")
  cat("  Total observations:", overall_stats_filtered$N_observations, "\n")
  cat("  Mean similarity:", round(overall_stats_filtered$Mean_Similarity * 100, 1), "%\n")
  cat("  SD similarity:", round(overall_stats_filtered$SD_Similarity * 100, 1), "%\n")
  cat("  Median similarity:", round(overall_stats_filtered$Median_Similarity * 100, 1), "%\n")
  cat("  IQR:", round(overall_stats_filtered$Q25_Similarity * 100, 1), "% -",
      round(overall_stats_filtered$Q75_Similarity * 100, 1), "%\n")
  
  improvement_overall <- (overall_stats_filtered$Mean_Similarity - overall_stats_raw$Mean_Similarity) * 100
  cat("\n✨ OVERALL IMPROVEMENT:", sprintf("%+.1f", improvement_overall), "percentage points\n")
}
cat("\n")

# ==============================================================================
# SECTION 4: BETA DIVERSITY DECOMPOSITION (OPTIONAL)
# ==============================================================================

cat("\n=== SECTION 4: BETA DIVERSITY DECOMPOSITION ===\n\n")

# This section uses betapart to decompose dissimilarity into turnover and nestedness
# Provides complementary information to taxonomic similarity

cat("Note: Beta diversity decomposition uses presence/absence at taxon level\n")
cat("      (not weighted by taxonomic distances like ΔS)\n\n")

cat("Processing beta diversity for each site...\n\n")

# For each site, create presence/absence matrices
beta_results <- list()

for (site in sites) {
  
  cat("━━━ Processing", site, "━━━\n")
  
  # Create date x taxon matrix for FAIRD
  faird_matrix <- faird_data %>%
    dplyr::filter(Site == site) %>%
    dplyr::select(Date, Taxon_ID) %>%
    dplyr::mutate(Present = 1) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = Taxon_ID, values_from = Present, values_fill = 0) %>%
    dplyr::arrange(Date)
  
  faird_dates <- faird_matrix$Date
  faird_matrix <- faird_matrix %>% dplyr::select(-Date) %>% as.matrix()
  rownames(faird_matrix) <- as.character(faird_dates)
  
  cat("  FAIRD matrix:", nrow(faird_matrix), "days ×", ncol(faird_matrix), "taxa\n")
  
  # Create date x taxon matrix for AMMOD
  ammod_matrix <- ammod_data %>%
    dplyr::filter(Site == site) %>%
    dplyr::select(Date, Taxon_ID) %>%
    dplyr::mutate(Present = 1) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = Taxon_ID, values_from = Present, values_fill = 0) %>%
    dplyr::arrange(Date)
  
  ammod_dates <- ammod_matrix$Date
  ammod_matrix <- ammod_matrix %>% dplyr::select(-Date) %>% as.matrix()
  rownames(ammod_matrix) <- as.character(ammod_dates)
  
  cat("  AMMOD matrix:", nrow(ammod_matrix), "days ×", ncol(ammod_matrix), "taxa\n")
  
  # Calculate beta diversity (if sufficient data)
  if (nrow(faird_matrix) >= 2) {
    beta_faird <- betapart::beta.pair(faird_matrix, index.family = "sorensen")
    beta_results[[paste0(site, "_FAIRD")]] <- beta_faird
    cat("  ✅ FAIRD beta diversity calculated\n")
  } else {
    cat("  ⚠️  FAIRD: Insufficient data (need ≥2 days)\n")
  }
  
  if (nrow(ammod_matrix) >= 2) {
    beta_ammod <- betapart::beta.pair(ammod_matrix, index.family = "sorensen")
    beta_results[[paste0(site, "_AMMOD")]] <- beta_ammod
    cat("  ✅ AMMOD beta diversity calculated\n")
  } else {
    cat("  ⚠️  AMMOD: Insufficient data (need ≥2 days)\n")
  }
  
  cat("\n")
}

cat("✅ Beta diversity decomposition complete\n")
cat("   Total results:", length(beta_results), "matrices\n\n")

# ==============================================================================
# SECTION 5: VISUALIZATION
# ==============================================================================

cat("\n=== SECTION 5: VISUALIZATION ===\n\n")

cat("Preparing data for plotting...\n")

# --- 5.1: Temporal concordance plot (style H1/H2) ---

# Prepare data for plotting - reshape to long format for both metrics
if (run_size_correction) {
  plot_data <- similarity_results %>%
    dplyr::filter(!is.na(Delta_S_raw) | !is.na(Delta_S_size_corrected)) %>%
    dplyr::mutate(
      Site_label = dplyr::case_when(
        Site == "Site1" ~ "Site 1 (Maize)",
        Site == "Site2" ~ "Site 2 (Maize)",
        Site == "Site3" ~ "Site 3 (Meadow)",
        Site == "Site4" ~ "Site 4 (Meadow)",
        TRUE ~ Site
      )
    ) %>%
    tidyr::pivot_longer(
      cols = c(Delta_S_raw, Delta_S_size_corrected),
      names_to = "Metric",
      values_to = "Similarity"
    ) %>%
    dplyr::filter(!is.na(Similarity)) %>%
    dplyr::mutate(
      Similarity_pct = Similarity * 100,
      Metric_label = dplyr::case_when(
        Metric == "Delta_S_raw" ~ "Raw (all AMMOD)",
        Metric == "Delta_S_size_corrected" ~ paste0("Size-corrected (≥", SIZE_THRESHOLD, "mm)"),
        TRUE ~ Metric
      )
    )
} else {
  plot_data <- similarity_results %>%
    dplyr::filter(!is.na(Delta_S_raw)) %>%
    dplyr::mutate(
      Similarity_pct = Delta_S_raw * 100,
      Site_label = dplyr::case_when(
        Site == "Site1" ~ "Site 1 (Maize)",
        Site == "Site2" ~ "Site 2 (Maize)",
        Site == "Site3" ~ "Site 3 (Meadow)",
        Site == "Site4" ~ "Site 4 (Meadow)",
        TRUE ~ Site
      ),
      Metric_label = "Raw (all AMMOD)"
    )
}

cat("  Plot data prepared:", length(unique(plot_data$Date)), "dates\n")
cat("  Sites in plot:", paste(unique(plot_data$Site_label), collapse=", "), "\n")
if (run_size_correction) {
  cat("  Metrics plotted: Raw similarity, Size-corrected similarity\n")
}
cat("\n")

cat("Creating temporal similarity plot...\n")

# Define colors for different metrics
metric_colors <- setNames(
  c(custom_colors["FAIRD"], "darkviolet"),
  c("Raw (all AMMOD)", paste0("Size-corrected (≥", SIZE_THRESHOLD, "mm)"))
)

# Create plot
if (run_size_correction) {
  p_similarity <- ggplot(plot_data, aes(x = Date, y = Similarity_pct, 
                                        color = Metric_label, group = Metric_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    facet_wrap(~Site_label, ncol = 2) +
    scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
    scale_color_manual(values = metric_colors) +
    labs(
      title = "Daily Taxonomic Similarity: FAIRD vs AMMOD (Genus-Level)",
      subtitle = paste0("Izsak-Price Index (ΔS) - Raw vs Size-Corrected (≥", SIZE_THRESHOLD, "mm) Comparison"),
      x = "Date",
      y = "Taxonomic Similarity (%)",
      color = "Metric",
      caption = "Dashed line = 50% similarity reference | Size correction tests if detection bias explains concordance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
} else {
  p_similarity <- ggplot(plot_data, aes(x = Date, y = Similarity_pct)) +
    geom_line(color = custom_colors["FAIRD"], linewidth = 0.8) +
    geom_point(color = custom_colors["FAIRD"], size = 2.5) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    facet_wrap(~Site_label, ncol = 2) +
    scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
    labs(
      title = "Daily Taxonomic Similarity: FAIRD vs AMMOD (Genus-Level)",
      subtitle = "Izsak-Price Index (ΔS) - Comparison at genus level (FAIRD's primary resolution)",
      x = "Date",
      y = "Taxonomic Similarity (%)",
      caption = "Dashed line = 50% similarity reference"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

print(p_similarity)

# Save plot
cat("\nSaving temporal similarity plot...\n")
ggsave("outputs/H3_FAIRD_AMMOD_similarity_temporal_GENUS.png", 
       plot = p_similarity, 
       width = 10, height = 8, dpi = 300)

cat("✅ Plot saved: outputs/H3_FAIRD_AMMOD_similarity_temporal_GENUS.png\n")
cat("   Dimensions: 10 × 8 inches (3000 × 2400 pixels @ 300 DPI)\n\n")

# --- 5.2: Threshold comparison plot (only if SIZE_THRESHOLD != 3) ---

cat("🔍 Checking if threshold comparison plot should be generated...\n")
cat("   SIZE_THRESHOLD =", SIZE_THRESHOLD, "mm\n")
cat("   run_baseline_comparison =", run_baseline_comparison, "\n")

if (run_baseline_comparison) {
  cat("   ✅ Generating threshold comparison plot (3mm vs", SIZE_THRESHOLD, "mm)\n\n")
  cat("Creating threshold comparison plot (3mm vs", SIZE_THRESHOLD, "mm)...\n")
  
  # Prepare data comparing 3mm baseline vs custom threshold
  plot_data_threshold <- similarity_results %>%
    dplyr::filter(!is.na(Delta_S_size_corrected_3mm) | !is.na(Delta_S_size_corrected)) %>%
    dplyr::mutate(
      Site_label = dplyr::case_when(
        Site == "Site1" ~ "Site 1 (Maize)",
        Site == "Site2" ~ "Site 2 (Maize)",
        Site == "Site3" ~ "Site 3 (Meadow)",
        Site == "Site4" ~ "Site 4 (Meadow)",
        TRUE ~ Site
      )
    ) %>%
    tidyr::pivot_longer(
      cols = c(Delta_S_size_corrected_3mm, Delta_S_size_corrected),
      names_to = "Threshold",
      values_to = "Similarity"
    ) %>%
    dplyr::filter(!is.na(Similarity)) %>%
    dplyr::mutate(
      Similarity_pct = Similarity * 100,
      Threshold_label = dplyr::case_when(
        Threshold == "Delta_S_size_corrected_3mm" ~ "Baseline (≥3mm)",
        Threshold == "Delta_S_size_corrected" ~ paste0("Custom (≥", SIZE_THRESHOLD, "mm)"),
        TRUE ~ Threshold
      )
    )
  
  # Define colors
  threshold_colors <- setNames(
    c(custom_colors["FAIRD"], "darkviolet"),
    c("Baseline (≥3mm)", paste0("Custom (≥", SIZE_THRESHOLD, "mm)"))
  )
  
  # Create plot
  p_threshold <- ggplot(plot_data_threshold, aes(x = Date, y = Similarity_pct, 
                                                 color = Threshold_label, 
                                                 group = Threshold_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    facet_wrap(~Site_label, ncol = 2) +
    scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
    scale_color_manual(values = threshold_colors) +
    labs(
      title = paste0("Threshold Comparison: 3mm (baseline) vs ", SIZE_THRESHOLD, "mm (custom)"),
      subtitle = "Izsak-Price Index (ΔS) - Effect of changing size threshold",
      x = "Date",
      y = "Taxonomic Similarity (%)",
      color = "Threshold",
      caption = paste0("Dashed line = 50% | Shows effect of ", SIZE_THRESHOLD, "mm vs standard 3mm threshold")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  
  print(p_threshold)
  
  # Save plot
  cat("\nSaving threshold comparison plot...\n")
  ggsave(paste0("outputs/H3_FAIRD_AMMOD_threshold_comparison_3mm_vs_", SIZE_THRESHOLD, "mm.png"),
         plot = p_threshold,
         width = 10, height = 8, dpi = 300)
  
  cat("✅ Plot saved: outputs/H3_FAIRD_AMMOD_threshold_comparison_3mm_vs_", SIZE_THRESHOLD, "mm.png\n")
  cat("   Dimensions: 10 × 8 inches (3000 × 2400 pixels @ 300 DPI)\n\n")
} else {
  cat("   ⏭️  Skipping threshold comparison plot (SIZE_THRESHOLD = 3mm, no comparison needed)\n\n")
}

# --- 5.3: Summary plot by site ---

cat("Creating summary plot by site...\n")

p_summary <- ggplot(summary_by_site, aes(x = Site, y = Mean_Similarity_raw * 100)) +
  geom_col(fill = custom_colors["FAIRD"], alpha = 0.7) +
  geom_errorbar(aes(ymin = (Mean_Similarity_raw - SD_Similarity_raw) * 100,
                    ymax = (Mean_Similarity_raw + SD_Similarity_raw) * 100),
                width = 0.2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  labs(
    title = "Mean Taxonomic Similarity by Site (Genus-Level)",
    subtitle = "Error bars = ±1 SD",
    x = "Site",
    y = "Mean Similarity (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold")
  )

print(p_summary)

cat("\nSaving summary plot...\n")
ggsave("outputs/H3_FAIRD_AMMOD_similarity_by_site_GENUS.png",
       plot = p_summary,
       width = 8, height = 6, dpi = 300)

cat("✅ Plot saved: outputs/H3_FAIRD_AMMOD_similarity_by_site_GENUS.png\n")
cat("   Dimensions: 8 × 6 inches (2400 × 1800 pixels @ 300 DPI)\n\n")

# ==============================================================================
# SECTION 6: STATISTICAL SUMMARY FOR MANUSCRIPT
# ==============================================================================

cat("\n=== SECTION 6: MANUSCRIPT SUMMARY ===\n\n")

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("H3: TAXONOMIC SIMILARITY ANALYSIS - FAIRD vs AMMOD (GENUS-LEVEL)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("TEMPORAL SCOPE:\n")
cat("  Period:", as.character(date_start), "to", as.character(date_end), "\n")
cat("  Duration: 12 days (first AMMOD bottle batch)\n\n")

cat("OVERALL TAXONOMIC SIMILARITY:\n")
cat("  Mean ΔS:", round(overall_stats_raw$Mean_Similarity * 100, 1), "%\n")
cat("  SD:", round(overall_stats_raw$SD_Similarity * 100, 1), "%\n")
cat("  Median:", round(overall_stats_raw$Median_Similarity * 100, 1), "%\n")
cat("  N observations:", overall_stats_raw$N_observations, "\n\n")

cat("SITE-SPECIFIC RESULTS:\n\n")
for (i in 1:nrow(summary_by_site)) {
  cat("  ", summary_by_site$Site[i], ":\n")
  cat("    Mean similarity:", round(summary_by_site$Mean_Similarity_raw[i] * 100, 1), "%\n")
  cat("    Range:", round(summary_by_site$Min_Similarity_raw[i] * 100, 1), "% -",
      round(summary_by_site$Max_Similarity_raw[i] * 100, 1), "%\n")
  cat("    Mean taxa detected - FAIRD:", round(summary_by_site$Mean_Taxa_FAIRD[i], 1), 
      "| AMMOD:", round(summary_by_site$Mean_Taxa_AMMOD_full[i], 1), "\n")
  cat("    Mean exact matches:", round(summary_by_site$Mean_Shared_raw[i], 1), "\n\n")
}

cat("KEY FINDINGS:\n\n")

# Identify best/worst performing sites
best_site <- summary_by_site %>% dplyr::slice_max(Mean_Similarity_raw, n = 1)
worst_site <- summary_by_site %>% dplyr::slice_min(Mean_Similarity_raw, n = 1)

cat("  🔹 Highest similarity:", best_site$Site, 
    "(", round(best_site$Mean_Similarity_raw * 100, 1), "%)\n")
cat("  🔹 Lowest similarity:", worst_site$Site,
    "(", round(worst_site$Mean_Similarity_raw * 100, 1), "%)\n")

# Check if there's a habitat effect (Maize vs Meadow)
maize_sites <- summary_by_site %>% dplyr::filter(Site %in% c("Site1", "Site2"))
meadow_sites <- summary_by_site %>% dplyr::filter(Site %in% c("Site3", "Site4"))

if (nrow(maize_sites) > 0 & nrow(meadow_sites) > 0) {
  cat("  🔹 Maize sites mean:", round(mean(maize_sites$Mean_Similarity_raw) * 100, 1), "%\n")
  cat("  🔹 Meadow sites mean:", round(mean(meadow_sites$Mean_Similarity_raw) * 100, 1), "%\n")
}

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("✅ Analysis complete - Results ready for manuscript integration\n\n")

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("💾 SAVING OUTPUT FILES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Save similarity results
cat("Saving detailed results...\n")
readr::write_csv(similarity_results, "outputs/H3_FAIRD_AMMOD_similarity_results_GENUS.csv")
cat("  ✅ outputs/H3_FAIRD_AMMOD_similarity_results_GENUS.csv\n")
cat("     Rows:", nrow(similarity_results), "| Columns:", ncol(similarity_results), "\n\n")

cat("Saving summary statistics...\n")
readr::write_csv(summary_by_site, "outputs/H3_FAIRD_AMMOD_similarity_summary_GENUS.csv")
cat("  ✅ outputs/H3_FAIRD_AMMOD_similarity_summary_GENUS.csv\n")
cat("     Rows:", nrow(summary_by_site), "| Columns:", ncol(summary_by_site), "\n\n")

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📁 COMPLETE OUTPUT SUMMARY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("FIGURES (PNG, 300 DPI):\n")
cat("  1. H3_FAIRD_AMMOD_similarity_temporal_GENUS.png (10×8 in)\n")
cat("  2. H3_FAIRD_AMMOD_similarity_by_site_GENUS.png (8×6 in)\n\n")

cat("DATA FILES (CSV):\n")
cat("  1. H3_FAIRD_AMMOD_similarity_results_GENUS.csv (daily values)\n")
cat("  2. H3_FAIRD_AMMOD_similarity_summary_GENUS.csv (site summaries)\n\n")

cat("✅ All files saved successfully\n")
cat("\n⚠️  NOTE: Analysis performed at GENUS level (FAIRD's primary resolution)\n")
cat("         AMMOD species collapsed to genera for fair comparison\n\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================