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
# CONFIGURATION: SIZE THRESHOLDS (fixed — 3mm and 5mm)
# ==============================================================================
# Both thresholds are always computed in a single run:
#   - 3mm: FAIRD physical detection limit (mesh size)
#   - 5mm: stricter sensitivity check
# Results feed Table 12 (main, Delta_S_raw) and Table S-X (sensitivity, 3mm vs 5mm).

SIZE_3MM <- 3  # mm — primary threshold (FAIRD detection limit)
SIZE_5MM <- 5  # mm — sensitivity threshold

cat("\n⚙️  FIXED SIZE THRESHOLDS: ", SIZE_3MM, "mm (primary) and", SIZE_5MM, "mm (sensitivity)\n\n")


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

# ==============================================================================
# TEMPORAL SCOPE - EXPANDED TO FULL COVERAGE (BATCH 1 + BATCH 2)
# ==============================================================================
# Date range: August 23 - September 13, 2023 (22 days)
# Exclusions:
#   1. September 4, 2023: Bottle change day (no sampling)
#   2. AMMOD4 August 30 - September 3: Mechanical failure period
#      - Rotation system stuck → non-representative sampling
#      - OTUs detected are contamination, not active captures
# ==============================================================================

date_start <- as.Date("2023-08-23")  # First day of sampling
date_end <- as.Date("2023-09-13")    # Last day of Batch 2
date_exclude <- as.Date("2023-09-04") # Bottle change day (no data)

# AMMOD4 mechanical failure period (rotation system stuck)
ammod4_failure_start <- as.Date("2023-08-28")
ammod4_failure_end <- as.Date("2023-09-03")

cat("\n=== TEMPORAL SCOPE (EXPANDED) ===\n\n")
cat("Date range:", as.character(date_start), "to", as.character(date_end), "\n")
cat("Duration: 22 days (Batch 1: Aug 23-Sep 3, Batch 2: Sep 5-13)\n")
cat("Excluded: Sep 4 (bottle change)\n")
cat("Excluded: AMMOD4", as.character(ammod4_failure_start), "to", 
    as.character(ammod4_failure_end), "(mechanical failure)\n")

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
  dplyr::filter(
    Date >= date_start & Date <= date_end,
    Date != date_exclude,  # Exclude bottle change day
    # Exclude AMMOD4 mechanical failure period
    !(Device == "AMMOD4" & Date >= ammod4_failure_start & Date <= ammod4_failure_end)
  ) %>%
  dplyr::select(Date, Device, Site, Class, Order, Suborder, Infraorder,
                Superfamily, Family, Subfamily, Genus, Species, Body_length, Body_width)
# Note: Count column not included - AMMOD taxonomy is presence/absence only
# Body_length & Body_width needed for size-corrected analysis

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

# --- 1.4: Create size-filtered AMMOD datasets (3mm and 5mm) ---
cat("🔧 Creating size-filtered AMMOD datasets...\n")

# Check if Body_length column exists
if (!"Body_length" %in% names(ammod_data)) {
  cat("⚠️  WARNING: Body_length column not found in AMMOD data\n")
  cat("   Size-filtered analyses will be skipped\n\n")
  stop("Body_length column required for H3 size-threshold analyses.")
}

# Full dataset (no size filter) — reference for Delta_S_raw
ammod_data_full <- ammod_data

# 3mm filter — primary threshold (FAIRD detection limit)
ammod_data_3mm <- ammod_data %>%
  dplyr::filter(!is.na(Body_length) & Body_length >= SIZE_3MM)

# 5mm filter — sensitivity threshold
ammod_data_5mm <- ammod_data %>%
  dplyr::filter(!is.na(Body_length) & Body_length >= SIZE_5MM)

n_original <- nrow(ammod_data_full)
n_3mm      <- nrow(ammod_data_3mm)
n_5mm      <- nrow(ammod_data_5mm)

cat("📊 Size filtering results:\n")
cat("   Original AMMOD OTUs:", n_original, "\n")
cat("   OTUs >= 3mm:", n_3mm, "(", round(n_3mm / n_original * 100, 1), "%)\n")
cat("   OTUs >= 5mm:", n_5mm, "(", round(n_5mm / n_original * 100, 1), "%)\n")
cat("   OTUs  < 3mm (excluded at primary threshold):",
    n_original - n_3mm, "(", round((n_original - n_3mm) / n_original * 100, 1), "%)\n\n")

# FAIRD size reference
if ("Body_length" %in% names(faird_data)) {
  faird_below_3mm <- sum(faird_data$Body_length < SIZE_3MM, na.rm = TRUE)
  faird_total     <- sum(!is.na(faird_data$Body_length))
  cat("   FAIRD < 3mm (for reference):", faird_below_3mm, "of", faird_total,
      "(", round(faird_below_3mm / faird_total * 100, 2), "%)\n\n")
}

# NOW filter Unknown from all datasets
cat("🔍 Checking for '#N/C' (No Clasificado) entries:\n")

faird_unknown <- sum(faird_data$Taxon_ID == "Unknown", na.rm = TRUE)
ammod_unknown <- sum(ammod_data$Taxon_ID == "Unknown", na.rm = TRUE)

cat("   FAIRD: ", faird_unknown, "individuals classified as 'Unknown'\n")
cat("          (", round(faird_unknown/nrow(faird_data)*100, 1), "% of total)\n")
cat("   AMMOD: ", ammod_unknown, "OTUs classified as 'Unknown'\n")
cat("          (", round(ammod_unknown/nrow(ammod_data)*100, 1), "% of total)\n")

if (faird_unknown > 0 | ammod_unknown > 0) {
  cat("\n   ⚠️  'Unknown' taxa will be excluded from similarity calculations\n")
  cat("      (no valid taxonomic classification at any level)\n")
  
  faird_data    <- faird_data    %>% dplyr::filter(Taxon_ID != "Unknown")
  ammod_data    <- ammod_data    %>% dplyr::filter(Taxon_ID != "Unknown")
  ammod_data_full <- ammod_data_full %>% dplyr::filter(Taxon_ID != "Unknown")
  ammod_data_3mm  <- ammod_data_3mm  %>% dplyr::filter(Taxon_ID != "Unknown")
  ammod_data_5mm  <- ammod_data_5mm  %>% dplyr::filter(Taxon_ID != "Unknown")
  
  cat("   ✅ Filtered data: FAIRD =", nrow(faird_data), "individuals | AMMOD =",
      nrow(ammod_data), "OTUs\n")
  cat("      ammod_data_full:", nrow(ammod_data_full), "rows (after Unknown filter)\n")
  cat("      ammod_data_3mm: ", nrow(ammod_data_3mm),  "rows (after Unknown filter)\n")
  cat("      ammod_data_5mm: ", nrow(ammod_data_5mm),  "rows (after Unknown filter)\n\n")
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
    # Primary metric: FAIRD vs all AMMOD (no size filter) — Table 12
    Delta_S_raw             = NA_real_,
    N_taxa_FAIRD            = NA_integer_,
    N_taxa_AMMOD_full       = NA_integer_,
    N_shared_exact_raw      = NA_integer_,
    # 3mm filter — primary size threshold (FAIRD detection limit) — Table S-X baseline
    Delta_S_3mm             = NA_real_,
    N_taxa_AMMOD_3mm        = NA_integer_,
    N_shared_exact_3mm      = NA_integer_,
    # 5mm filter — sensitivity threshold — Table S-X comparison
    Delta_S_5mm             = NA_real_,
    N_taxa_AMMOD_5mm        = NA_integer_,
    N_shared_exact_5mm      = NA_integer_
  )

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("🔄 CALCULATING DAILY TAXONOMIC SIMILARITY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Total calculations:", nrow(similarity_results), "(", length(all_dates), "days ×",
    length(sites), "sites )\n")
cat("\n   Computing THREE similarity metrics per day:\n")
cat("   1. RAW:  FAIRD vs all AMMOD (no size filter)  → Table 12\n")
cat("   2. 3mm:  FAIRD vs AMMOD >= 3mm               → Table S-X baseline\n")
cat("   3. 5mm:  FAIRD vs AMMOD >= 5mm               → Table S-X comparison\n\n")

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
  
  # Get AMMOD taxa — full dataset (no size filter)
  ammod_taxa_day_full <- ammod_data_full %>%
    dplyr::filter(Date == current_date, Site == current_site) %>%
    dplyr::pull(Taxon_ID) %>%
    unique()
  
  # Get AMMOD taxa — 3mm filter
  ammod_taxa_day_3mm <- ammod_data_3mm %>%
    dplyr::filter(Date == current_date, Site == current_site) %>%
    dplyr::pull(Taxon_ID) %>%
    unique()
  
  # Get AMMOD taxa — 5mm filter
  ammod_taxa_day_5mm <- ammod_data_5mm %>%
    dplyr::filter(Date == current_date, Site == current_site) %>%
    dplyr::pull(Taxon_ID) %>%
    unique()
  
  # Store taxon counts
  similarity_results$N_taxa_FAIRD[i]       <- length(faird_taxa_day)
  similarity_results$N_taxa_AMMOD_full[i]  <- length(ammod_taxa_day_full)
  similarity_results$N_taxa_AMMOD_3mm[i]   <- length(ammod_taxa_day_3mm)
  similarity_results$N_taxa_AMMOD_5mm[i]   <- length(ammod_taxa_day_5mm)
  similarity_results$N_shared_exact_raw[i] <- length(intersect(faird_taxa_day, ammod_taxa_day_full))
  similarity_results$N_shared_exact_3mm[i] <- length(intersect(faird_taxa_day, ammod_taxa_day_3mm))
  similarity_results$N_shared_exact_5mm[i] <- length(intersect(faird_taxa_day, ammod_taxa_day_5mm))
  
  # Calculate RAW similarity (FAIRD vs all AMMOD — Table 12)
  if (length(faird_taxa_day) > 0 && length(ammod_taxa_day_full) > 0) {
    similarity_results$Delta_S_raw[i] <- calculate_taxonomic_similarity(
      faird_taxa_day, ammod_taxa_day_full, dist_matrix
    )
  }
  
  # Calculate 3mm similarity (Table S-X baseline)
  if (length(faird_taxa_day) > 0 && length(ammod_taxa_day_3mm) > 0) {
    similarity_results$Delta_S_3mm[i] <- calculate_taxonomic_similarity(
      faird_taxa_day, ammod_taxa_day_3mm, dist_matrix
    )
  }
  
  # Calculate 5mm similarity (Table S-X comparison)
  if (length(faird_taxa_day) > 0 && length(ammod_taxa_day_5mm) > 0) {
    similarity_results$Delta_S_5mm[i] <- calculate_taxonomic_similarity(
      faird_taxa_day, ammod_taxa_day_5mm, dist_matrix
    )
  }
}

cat("\n✅ Similarity calculations complete\n\n")

# Quick summary of results
valid_results_raw <- sum(!is.na(similarity_results$Delta_S_raw))
valid_results_3mm <- sum(!is.na(similarity_results$Delta_S_3mm))
valid_results_5mm <- sum(!is.na(similarity_results$Delta_S_5mm))
cat("📊 Calculation summary:\n")
cat("   RAW (all AMMOD)    - Valid results:", valid_results_raw, "/", nrow(similarity_results), "\n")
cat("   3mm filter         - Valid results:", valid_results_3mm, "/", nrow(similarity_results), "\n")
cat("   5mm filter         - Valid results:", valid_results_5mm, "/", nrow(similarity_results), "\n")
cat("   Missing (NA) raw:  ", sum(is.na(similarity_results$Delta_S_raw)), "\n")
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
    N_days               = dplyr::n(),
    Mean_Similarity_raw  = mean(Delta_S_raw,       na.rm = TRUE),
    SD_Similarity_raw    = sd(Delta_S_raw,         na.rm = TRUE),
    Min_Similarity_raw   = min(Delta_S_raw,         na.rm = TRUE),
    Max_Similarity_raw   = max(Delta_S_raw,         na.rm = TRUE),
    Mean_Taxa_FAIRD      = mean(N_taxa_FAIRD,       na.rm = TRUE),
    Mean_Taxa_AMMOD_full = mean(N_taxa_AMMOD_full,  na.rm = TRUE),
    Mean_Shared_raw      = mean(N_shared_exact_raw, na.rm = TRUE),
    .groups = "drop"
  )

# Summary for 3mm filter
summary_by_site_3mm <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_3mm)) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    Mean_Similarity_3mm  = mean(Delta_S_3mm,       na.rm = TRUE),
    SD_Similarity_3mm    = sd(Delta_S_3mm,         na.rm = TRUE),
    Min_Similarity_3mm   = min(Delta_S_3mm,         na.rm = TRUE),
    Max_Similarity_3mm   = max(Delta_S_3mm,         na.rm = TRUE),
    Mean_Taxa_AMMOD_3mm  = mean(N_taxa_AMMOD_3mm,  na.rm = TRUE),
    Mean_Shared_3mm      = mean(N_shared_exact_3mm, na.rm = TRUE),
    .groups = "drop"
  )

# Summary for 5mm filter
summary_by_site_5mm <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_5mm)) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    Mean_Similarity_5mm  = mean(Delta_S_5mm,       na.rm = TRUE),
    SD_Similarity_5mm    = sd(Delta_S_5mm,         na.rm = TRUE),
    Min_Similarity_5mm   = min(Delta_S_5mm,         na.rm = TRUE),
    Max_Similarity_5mm   = max(Delta_S_5mm,         na.rm = TRUE),
    Mean_Taxa_AMMOD_5mm  = mean(N_taxa_AMMOD_5mm,  na.rm = TRUE),
    Mean_Shared_5mm      = mean(N_shared_exact_5mm, na.rm = TRUE),
    .groups = "drop"
  )

# Combine all summaries
summary_by_site <- summary_by_site_raw %>%
  dplyr::left_join(summary_by_site_3mm, by = "Site") %>%
  dplyr::left_join(summary_by_site_5mm, by = "Site")

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
    
    # 3mm filter summary
    if (!is.na(site_data$Mean_Similarity_3mm)) {
      cat("\n   SIZE-FILTERED (vs AMMOD >= 3mm):\n")
      cat("     Mean:", round(site_data$Mean_Similarity_3mm * 100, 1), "% ±",
          round(site_data$SD_Similarity_3mm * 100, 1), "%\n")
      cat("     Range:", round(site_data$Min_Similarity_3mm * 100, 1), "% -",
          round(site_data$Max_Similarity_3mm * 100, 1), "%\n")
      cat("     AMMOD taxa/day:", round(site_data$Mean_Taxa_AMMOD_3mm, 1), "\n")
    }
    
    # 5mm filter summary
    if (!is.na(site_data$Mean_Similarity_5mm)) {
      cat("\n   SIZE-FILTERED (vs AMMOD >= 5mm):\n")
      cat("     Mean:", round(site_data$Mean_Similarity_5mm * 100, 1), "% ±",
          round(site_data$SD_Similarity_5mm * 100, 1), "%\n")
      cat("     Range:", round(site_data$Min_Similarity_5mm * 100, 1), "% -",
          round(site_data$Max_Similarity_5mm * 100, 1), "%\n")
      cat("     AMMOD taxa/day:", round(site_data$Mean_Taxa_AMMOD_5mm, 1), "\n")
      improvement_5mm <- (site_data$Mean_Similarity_5mm - site_data$Mean_Similarity_3mm) * 100
      cat("\n     Improvement 3mm -> 5mm:", sprintf("%+.1f", improvement_5mm), "pp\n")
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

# SIZE-CORRECTED similarity overall — 3mm
overall_stats_3mm <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_3mm)) %>%
  dplyr::summarise(
    N_observations    = dplyr::n(),
    Mean_Similarity   = mean(Delta_S_3mm,   na.rm = TRUE),
    SD_Similarity     = sd(Delta_S_3mm,     na.rm = TRUE),
    Median_Similarity = median(Delta_S_3mm, na.rm = TRUE),
    Q25_Similarity    = quantile(Delta_S_3mm, 0.25, na.rm = TRUE),
    Q75_Similarity    = quantile(Delta_S_3mm, 0.75, na.rm = TRUE)
  )

cat("\nSIZE-CORRECTED SIMILARITY (vs AMMOD >= 3mm):\n")
cat("  Total observations:", overall_stats_3mm$N_observations, "\n")
cat("  Mean similarity:",   round(overall_stats_3mm$Mean_Similarity   * 100, 1), "%\n")
cat("  SD similarity:",     round(overall_stats_3mm$SD_Similarity     * 100, 1), "%\n")
cat("  Median similarity:", round(overall_stats_3mm$Median_Similarity * 100, 1), "%\n")
cat("  IQR:", round(overall_stats_3mm$Q25_Similarity * 100, 1), "% -",
    round(overall_stats_3mm$Q75_Similarity * 100, 1), "%\n")

# SIZE-CORRECTED similarity overall — 5mm
overall_stats_5mm <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_5mm)) %>%
  dplyr::summarise(
    N_observations    = dplyr::n(),
    Mean_Similarity   = mean(Delta_S_5mm,   na.rm = TRUE),
    SD_Similarity     = sd(Delta_S_5mm,     na.rm = TRUE),
    Median_Similarity = median(Delta_S_5mm, na.rm = TRUE),
    Q25_Similarity    = quantile(Delta_S_5mm, 0.25, na.rm = TRUE),
    Q75_Similarity    = quantile(Delta_S_5mm, 0.75, na.rm = TRUE)
  )

cat("\nSIZE-CORRECTED SIMILARITY (vs AMMOD >= 5mm):\n")
cat("  Total observations:", overall_stats_5mm$N_observations, "\n")
cat("  Mean similarity:",   round(overall_stats_5mm$Mean_Similarity   * 100, 1), "%\n")
cat("  SD similarity:",     round(overall_stats_5mm$SD_Similarity     * 100, 1), "%\n")
cat("  Median similarity:", round(overall_stats_5mm$Median_Similarity * 100, 1), "%\n")
cat("  IQR:", round(overall_stats_5mm$Q25_Similarity * 100, 1), "% -",
    round(overall_stats_5mm$Q75_Similarity * 100, 1), "%\n")

improvement_overall <- (overall_stats_5mm$Mean_Similarity - overall_stats_3mm$Mean_Similarity) * 100
cat("\nImprovement 3mm -> 5mm:", sprintf("%+.1f", improvement_overall), "pp\n")
cat("\n")

# ---- HABITAT-LEVEL STATISTICS (WEIGHTED AVERAGES) ----------------------------
# Calculates per-habitat statistics using individual daily observations (weighted)
# rather than simple means of site-level means (unweighted).
# This ensures each observation day contributes equally regardless of site sample size.

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📊 HABITAT-LEVEL STATISTICS (WEIGHTED BY OBSERVATIONS)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("ℹ️  Method: Weighted averages from individual daily observations\n")
cat("   (each observation day contributes equally, regardless of site sample size)\n\n")

# Assign habitat to each observation
habitat_data <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw)) %>%
  dplyr::mutate(
    Habitat = dplyr::case_when(
      Site %in% c("Site1", "Site2") ~ "Maize",
      Site %in% c("Site3", "Site4") ~ "Meadow",
      TRUE ~ NA_character_
    )
  )

# Calculate weighted stats per habitat
habitat_stats <- habitat_data %>%
  dplyr::group_by(Habitat) %>%
  dplyr::summarise(
    N_observations       = dplyr::n(),
    Mean_similarity      = mean(Delta_S_raw, na.rm = TRUE),
    SD_similarity        = sd(Delta_S_raw, na.rm = TRUE),
    Median_similarity    = median(Delta_S_raw, na.rm = TRUE),
    Days_below_50pct     = sum(Delta_S_raw < 0.5, na.rm = TRUE),
    Pct_below_50         = sum(Delta_S_raw < 0.5, na.rm = TRUE) / dplyr::n() * 100,
    Mean_FAIRD_genera    = mean(N_taxa_FAIRD, na.rm = TRUE),
    Mean_AMMOD_genera    = mean(N_taxa_AMMOD_full, na.rm = TRUE),
    Mean_shared_exact    = mean(N_shared_exact_raw, na.rm = TRUE),
    .groups = "drop"
  )

# Print habitat statistics
for (h in c("Maize", "Meadow")) {
  hs <- habitat_stats %>% dplyr::filter(Habitat == h)
  cat(h, " (", 
      paste(if (h == "Maize") "Site1 + Site2" else "Site3 + Site4"), 
      "):\n", sep = "")
  cat("  N observations:       ", hs$N_observations, "\n")
  cat("  Mean similarity:      ", round(hs$Mean_similarity * 100, 1), "% ±", 
      round(hs$SD_similarity * 100, 1), "%\n")
  cat("  Days below 50%:       ", hs$Days_below_50pct, "/", hs$N_observations,
      " (", round(hs$Pct_below_50, 1), "%)\n", sep = "")
  cat("  FAIRD genera/day:     ", round(hs$Mean_FAIRD_genera, 1), "\n")
  cat("  AMMOD genera/day:     ", round(hs$Mean_AMMOD_genera, 1), "\n")
  cat("  Exact matches/day:    ", round(hs$Mean_shared_exact, 1), "\n\n")
}

# Calculate and print differences
maize_row  <- habitat_stats %>% dplyr::filter(Habitat == "Maize")
meadow_row <- habitat_stats %>% dplyr::filter(Habitat == "Meadow")

diff_pp          <- (meadow_row$Mean_similarity - maize_row$Mean_similarity) * 100
diff_faird_gen   <- meadow_row$Mean_FAIRD_genera - maize_row$Mean_FAIRD_genera
diff_faird_pct   <- diff_faird_gen / maize_row$Mean_FAIRD_genera * 100
diff_ammod_gen   <- meadow_row$Mean_AMMOD_genera - maize_row$Mean_AMMOD_genera
diff_ammod_pct   <- diff_ammod_gen / maize_row$Mean_AMMOD_genera * 100
diff_days_below  <- meadow_row$Pct_below_50 - maize_row$Pct_below_50

cat("DIFFERENCE (Meadow - Maize):\n")
cat("  Similarity:           ", sprintf("%+.1f", diff_pp), "pp\n")
cat("  Days below 50%:       ", sprintf("%+.1f", diff_days_below), "pp\n")
cat("  FAIRD genera/day:     ", sprintf("%+.1f", diff_faird_gen), 
    " (", sprintf("%+.0f", diff_faird_pct), "%)\n", sep = "")
cat("  AMMOD genera/day:     ", sprintf("%+.1f", diff_ammod_gen),
    " (", sprintf("%+.0f", diff_ammod_pct), "%)\n", sep = "")

cat("\n")

# Also show comparison: weighted vs simple (site-level) means for transparency
cat("COMPARISON: Weighted vs Simple (site-level) means:\n")
maize_simple  <- mean(c(
  mean(habitat_data$Delta_S_raw[habitat_data$Site == "Site1"], na.rm = TRUE),
  mean(habitat_data$Delta_S_raw[habitat_data$Site == "Site2"], na.rm = TRUE)
))
meadow_simple <- mean(c(
  mean(habitat_data$Delta_S_raw[habitat_data$Site == "Site3"], na.rm = TRUE),
  mean(habitat_data$Delta_S_raw[habitat_data$Site == "Site4"], na.rm = TRUE)
))
cat("  Maize  - Weighted:", round(maize_row$Mean_similarity * 100, 1), 
    "% | Simple:", round(maize_simple * 100, 1), "%\n")
cat("  Meadow - Weighted:", round(meadow_row$Mean_similarity * 100, 1), 
    "% | Simple:", round(meadow_simple * 100, 1), "%\n")
cat("  Diff   - Weighted:", sprintf("%+.1f", diff_pp), 
    "pp | Simple:", sprintf("%+.1f", (meadow_simple - maize_simple) * 100), "pp\n")

cat("\n")

# ---- LINEAR TREND ANALYSIS ---------------------------------------------------
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📈 LINEAR TREND ANALYSIS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

trend_data <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw)) %>%
  dplyr::mutate(Day_numeric = as.numeric(Date - min(Date)))

trend_model <- lm(Delta_S_raw ~ Day_numeric, data = trend_data)
trend_summary <- summary(trend_model)

cat("Model: Delta_S_raw ~ Day_numeric (n =", nrow(trend_data), ")\n")
cat("  Slope:  ", sprintf("%+.4f", coef(trend_model)[2]),
    " (", sprintf("%+.2f", coef(trend_model)[2] * 100), "% per day)\n", sep = "")
cat("  t(", trend_summary$df[2], ") = ", 
    round(trend_summary$coefficients[2, 3], 2), 
    ", p = ", round(trend_summary$coefficients[2, 4], 3), "\n", sep = "")
cat("  p-value exact:", format.pval(trend_summary$coefficients[2, 4], digits = 3), "\n")
cat("  Significance: ", ifelse(trend_summary$coefficients[2, 4] < 0.05, 
                               "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
cat("  Interpretation: ",
    ifelse(coef(trend_model)[2] > 0, "Slight improvement", "Slight decline"),
    " trend but ",
    ifelse(trend_summary$coefficients[2, 4] < 0.05, "statistically robust", 
           "not statistically robust"), "\n\n")

# ---- DIAGNOSTIC: What drives the declining trend? ----------------------------
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("🔍 TREND DIAGNOSTIC: What drives the declining similarity?\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

diag_data <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw)) %>%
  dplyr::mutate(Day_numeric = as.numeric(Date - min(Date)))

# Trend in FAIRD genera/day
mod_faird <- lm(N_taxa_FAIRD ~ Day_numeric, data = diag_data)
s_faird <- summary(mod_faird)

# Trend in AMMOD genera/day
mod_ammod <- lm(N_taxa_AMMOD_full ~ Day_numeric, data = diag_data)
s_ammod <- summary(mod_ammod)

# Trend in exact matches
mod_shared <- lm(N_shared_exact_raw ~ Day_numeric, data = diag_data)
s_shared <- summary(mod_shared)

cat("1. FAIRD genera/day trend:\n")
cat("   Slope: ", sprintf("%+.3f", coef(mod_faird)[2]), " genera/day\n", sep = "")
cat("   t(", s_faird$df[2], ") = ", round(s_faird$coefficients[2, 3], 2),
    ", p = ", format.pval(s_faird$coefficients[2, 4], digits = 3), "\n", sep = "")
cat("   Direction: ", ifelse(coef(mod_faird)[2] > 0, "INCREASING", "DECLINING"), "\n\n")

cat("2. AMMOD genera/day trend:\n")
cat("   Slope: ", sprintf("%+.3f", coef(mod_ammod)[2]), " genera/day\n", sep = "")
cat("   t(", s_ammod$df[2], ") = ", round(s_ammod$coefficients[2, 3], 2),
    ", p = ", format.pval(s_ammod$coefficients[2, 4], digits = 3), "\n", sep = "")
cat("   Direction: ", ifelse(coef(mod_ammod)[2] > 0, "INCREASING", "DECLINING"), "\n\n")

cat("3. Exact matches/day trend:\n")
cat("   Slope: ", sprintf("%+.3f", coef(mod_shared)[2]), " matches/day\n", sep = "")
cat("   t(", s_shared$df[2], ") = ", round(s_shared$coefficients[2, 3], 2),
    ", p = ", format.pval(s_shared$coefficients[2, 4], digits = 3), "\n", sep = "")
cat("   Direction: ", ifelse(coef(mod_shared)[2] > 0, "INCREASING", "DECLINING"), "\n\n")

cat("INTERPRETATION:\n")
if (coef(mod_faird)[2] < 0 & s_faird$coefficients[2, 4] < 0.05) {
  cat("  → FAIRD detections decline significantly over time\n")
}
if (coef(mod_ammod)[2] < 0 & s_ammod$coefficients[2, 4] < 0.05) {
  cat("  → AMMOD richness also declines (seasonal effect)\n")
} else {
  cat("  → AMMOD richness does NOT decline significantly\n")
}
if (coef(mod_faird)[2] < 0 & coef(mod_ammod)[2] >= 0) {
  cat("  → Asymmetric decline: FAIRD loses taxa while AMMOD maintains richness\n")
  cat("    This explains declining similarity (fewer FAIRD taxa to match against stable AMMOD baseline)\n")
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

# --- 5.1: Temporal similarity plot — RAW vs 3mm vs 5mm ---

plot_data <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw) | !is.na(Delta_S_3mm) | !is.na(Delta_S_5mm)) %>%
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
    cols = c(Delta_S_raw, Delta_S_3mm, Delta_S_5mm),
    names_to  = "Metric",
    values_to = "Similarity"
  ) %>%
  dplyr::filter(!is.na(Similarity)) %>%
  dplyr::mutate(
    Similarity_pct = Similarity * 100,
    Metric_label = dplyr::case_when(
      Metric == "Delta_S_raw" ~ "Raw (all AMMOD)",
      Metric == "Delta_S_3mm" ~ "Size-filtered (>= 3mm)",
      Metric == "Delta_S_5mm" ~ "Size-filtered (>= 5mm)",
      TRUE ~ Metric
    )
  )

cat("  Plot data prepared:", length(unique(plot_data$Date)), "dates\n")
cat("  Sites in plot:", paste(unique(plot_data$Site_label), collapse = ", "), "\n\n")

cat("Creating temporal similarity plot...\n")

metric_colors <- c(
  "Raw (all AMMOD)"        = unname(custom_colors["FAIRD"]),
  "Size-filtered (>= 3mm)" = "darkorange",
  "Size-filtered (>= 5mm)" = "darkviolet"
)

p_similarity <- ggplot(plot_data,
                       aes(x = Date, y = Similarity_pct, color = Metric_label, group = Metric_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  facet_wrap(~Site_label, ncol = 2) +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  scale_color_manual(values = metric_colors) +
  labs(
    title    = "Daily Taxonomic Similarity: FAIRD vs AMMOD (Genus-Level)",
    subtitle = "Izsak-Price Index (Delta_S) — Raw vs size-filtered AMMOD (3mm and 5mm)",
    x        = "Date",
    y        = "Taxonomic Similarity (%)",
    color    = "Metric",
    caption  = "Dashed line = 50% similarity reference"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

print(p_similarity)

ggsave("outputs/H3_FAIRD_AMMOD_similarity_temporal_GENUS.png",
       plot = p_similarity, width = 10, height = 8, dpi = 300)
cat("✅ outputs/H3_FAIRD_AMMOD_similarity_temporal_GENUS.png\n\n")

# --- 5.1b: Clean publication figure — RAW similarity only ---
cat("Creating clean publication figure (Delta_S_raw, single line per site)...\n")

plot_data_clean <- similarity_results %>%
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
    Habitat = dplyr::if_else(Site %in% c("Site1", "Site2"), "Maize", "Meadow")
  )

site_colors_pub <- c(
  "Site 1 (Maize)"  = unname(custom_colors["Site1"]),
  "Site 2 (Maize)"  = unname(custom_colors["Site2"]),
  "Site 3 (Meadow)" = unname(custom_colors["Site3"]),
  "Site 4 (Meadow)" = unname(custom_colors["Site4"])
)

p_similarity_clean <- ggplot(plot_data_clean,
                             aes(x = Date, y = Similarity_pct, color = Site_label, group = Site_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  facet_wrap(~Site_label, ncol = 2) +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  scale_color_manual(values = site_colors_pub) +
  labs(
    title   = "Daily Taxonomic Similarity: FAIRD vs AMMOD (Genus-Level)",
    subtitle = "Izsak-Price Index (\u0394S) | FAIRD vs all AMMOD taxa | Genus-level comparison",
    x       = "Date",
    y       = "Taxonomic Similarity (%)",
    caption = "Dashed line = 50% similarity reference"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

print(p_similarity_clean)

ggsave("outputs/H3_FAIRD_AMMOD_similarity_temporal_GENUS_clean.png",
       plot = p_similarity_clean, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("outputs/H3_FAIRD_AMMOD_similarity_temporal_GENUS_clean.pdf",
       plot = p_similarity_clean, width = 10, height = 8, device = cairo_pdf)
cat("✅ outputs/H3_FAIRD_AMMOD_similarity_temporal_GENUS_clean.png/.pdf\n\n")

# --- 5.2: Threshold comparison plot (3mm vs 5mm) ---
cat("Creating threshold comparison plot (3mm vs 5mm)...\n")

plot_data_threshold <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_3mm) | !is.na(Delta_S_5mm)) %>%
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
    cols      = c(Delta_S_3mm, Delta_S_5mm),
    names_to  = "Threshold",
    values_to = "Similarity"
  ) %>%
  dplyr::filter(!is.na(Similarity)) %>%
  dplyr::mutate(
    Similarity_pct  = Similarity * 100,
    Threshold_label = dplyr::case_when(
      Threshold == "Delta_S_3mm" ~ "Baseline (>= 3mm)",
      Threshold == "Delta_S_5mm" ~ "Sensitivity (>= 5mm)",
      TRUE ~ Threshold
    )
  )

threshold_colors <- c(
  "Baseline (>= 3mm)"    = "darkorange",
  "Sensitivity (>= 5mm)" = "darkviolet"
)

p_threshold <- ggplot(plot_data_threshold,
                      aes(x = Date, y = Similarity_pct, color = Threshold_label, group = Threshold_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  facet_wrap(~Site_label, ncol = 2) +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  scale_color_manual(values = threshold_colors) +
  labs(
    title    = "Threshold Sensitivity: 3mm vs 5mm AMMOD filter",
    subtitle = "Izsak-Price Index (Delta_S) — Effect of size threshold on taxonomic similarity",
    x        = "Date",
    y        = "Taxonomic Similarity (%)",
    color    = "Threshold",
    caption  = "Dashed line = 50% similarity reference"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

print(p_threshold)

ggsave("outputs/H3_FAIRD_AMMOD_threshold_comparison_3mm_vs_5mm.png",
       plot = p_threshold, width = 10, height = 8, dpi = 300)
cat("✅ outputs/H3_FAIRD_AMMOD_threshold_comparison_3mm_vs_5mm.png\n\n")

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
       plot = p_summary, width = 8, height = 6, dpi = 300)

cat("✅ Plot saved: outputs/H3_FAIRD_AMMOD_similarity_by_site_GENUS.png\n")
cat("   Dimensions: 8 × 6 inches (2400 × 1800 pixels @ 300 DPI)\n\n")

# ==============================================================================
# SECTION 6: STATISTICAL SUMMARY
# ==============================================================================

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

cat("✅ Analysis complete \n\n")

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("💾 SAVING OUTPUT FILES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Save similarity results
cat("Saving detailed results...\n")

# Store the path in a variable to avoid repeating the paste0 logic
results_path <- "outputs/H3_FAIRD_AMMOD_similarity_results_GENUS.csv"

readr::write_csv(similarity_results, results_path)

cat("  ✅", results_path, "\n")
cat("      Rows:", nrow(similarity_results), "| Columns:", ncol(similarity_results), "\n\n")

# Save summary statistics
cat("Saving summary statistics...\n")

summary_path <- "outputs/H3_FAIRD_AMMOD_similarity_summary_GENUS.csv"

readr::write_csv(summary_by_site, summary_path)

cat("  ✅", summary_path, "\n")
cat("      Rows:", nrow(summary_by_site), "| Columns:", ncol(summary_by_site), "\n\n")

# ==============================================================================
# ==============================================================================
# SECTION 7a: TABLE 12 — MAIN H3 RESULTS
# ==============================================================================
# Publication-ready summary of daily FAIRD-AMMOD taxonomic similarity (H3).
# Always produced in every run.
# Uses Delta_S_raw (FAIRD vs ALL AMMOD taxa, no size filter), which is the
# primary H3 metric reported in the manuscript. Size-threshold sensitivity
# is handled separately in Section 7b (Table S-X).
# ==============================================================================

cat("\n=== SECTION 7a: TABLE 12 — MAIN H3 RESULTS (similarity + genera/day) ===\n\n")
cat("NOTE: Mean similarity = Delta_S_raw (FAIRD vs all AMMOD, no size filter).\n")
cat("      Size-threshold sensitivity (3mm vs 5mm) is in Table S-X (Section 7b).\n\n")

# --- Per-site rows ---
table12_by_site <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw)) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    Habitat       = dplyr::first(dplyr::case_when(
      Site %in% c("Site1", "Site2") ~ "Maize",
      TRUE                           ~ "Meadow")),
    n_days        = dplyr::n(),
    mean_DS       = round(mean(Delta_S_raw,        na.rm = TRUE) * 100, 1),
    sd_DS         = round(sd(Delta_S_raw,          na.rm = TRUE) * 100, 1),
    faird_genera  = round(mean(N_taxa_FAIRD,       na.rm = TRUE), 1),
    ammod_genera  = round(mean(N_taxa_AMMOD_full,  na.rm = TRUE), 1),
    exact_matches = round(mean(N_shared_exact_raw, na.rm = TRUE), 1),
    .groups = "drop"
  )

# --- Grand mean row (unweighted mean across all n=68 observations) ---
table12_grand <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_raw)) %>%
  dplyr::summarise(
    Site          = "All sites",
    Habitat       = "\u2014",
    n_days        = dplyr::n(),
    mean_DS       = round(mean(Delta_S_raw,        na.rm = TRUE) * 100, 1),
    sd_DS         = round(sd(Delta_S_raw,          na.rm = TRUE) * 100, 1),
    faird_genera  = round(mean(N_taxa_FAIRD,       na.rm = TRUE), 1),
    ammod_genera  = round(mean(N_taxa_AMMOD_full,  na.rm = TRUE), 1),
    exact_matches = round(mean(N_shared_exact_raw, na.rm = TRUE), 1)
  )

table12 <- dplyr::bind_rows(table12_by_site, table12_grand)

table12_display <- table12 %>%
  dplyr::rename(
    `Site`                = Site,
    `Habitat`             = Habitat,
    `n days`              = n_days,
    `Mean similarity (%)` = mean_DS,
    `SD`                  = sd_DS,
    `FAIRD genera/day`    = faird_genera,
    `AMMOD genera/day`    = ammod_genera,
    `Exact matches/day`   = exact_matches
  )

cat("--- Table 12: Daily taxonomic similarity (Izsak-Price, genus level) by site ---\n\n")
print(table12_display, n = Inf)

readr::write_csv(table12_display, "outputs/H3_Table12_similarity_by_site.csv")
cat("\n[OK] outputs/H3_Table12_similarity_by_site.csv\n\n")


# ==============================================================================
# SECTION 7b: TABLE S-X — SIZE THRESHOLD SENSITIVITY (3mm vs 5mm)
# ==============================================================================
# Supplementary table comparing FAIRD-AMMOD similarity and AMMOD taxon
# retention under the 3mm vs 5mm thresholds. Always produced.
# Note: n per site matches Table 12 because both 3mm and 5mm are always
#       computed on the same days.
# ==============================================================================

cat("\n=== SECTION 7b: TABLE S-X — THRESHOLD SENSITIVITY (3mm vs 5mm) ===\n\n")

# --- Per-site rows ---
tableS_by_site <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_3mm), !is.na(Delta_S_5mm)) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    Habitat        = dplyr::first(dplyr::case_when(
      Site %in% c("Site1", "Site2") ~ "Maize",
      TRUE                           ~ "Meadow")),
    n_days         = dplyr::n(),
    ds_3mm         = round(mean(Delta_S_3mm,       na.rm = TRUE) * 100, 1),
    ds_5mm         = round(mean(Delta_S_5mm,       na.rm = TRUE) * 100, 1),
    delta_pp       = round(mean(Delta_S_5mm,       na.rm = TRUE) * 100 -
                             mean(Delta_S_3mm,       na.rm = TRUE) * 100, 1),
    ammod_taxa_3mm = round(mean(N_taxa_AMMOD_3mm,  na.rm = TRUE), 1),
    ammod_taxa_5mm = round(mean(N_taxa_AMMOD_5mm,  na.rm = TRUE), 1),
    pct_retained   = round(mean(N_taxa_AMMOD_5mm,  na.rm = TRUE) /
                             mean(N_taxa_AMMOD_3mm,  na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

# --- Grand mean row ---
tableS_grand <- similarity_results %>%
  dplyr::filter(!is.na(Delta_S_3mm), !is.na(Delta_S_5mm)) %>%
  dplyr::summarise(
    Site           = "All sites",
    Habitat        = "\u2014",
    n_days         = dplyr::n(),
    ds_3mm         = round(mean(Delta_S_3mm,       na.rm = TRUE) * 100, 1),
    ds_5mm         = round(mean(Delta_S_5mm,       na.rm = TRUE) * 100, 1),
    delta_pp       = round(mean(Delta_S_5mm,       na.rm = TRUE) * 100 -
                             mean(Delta_S_3mm,       na.rm = TRUE) * 100, 1),
    ammod_taxa_3mm = round(mean(N_taxa_AMMOD_3mm,  na.rm = TRUE), 1),
    ammod_taxa_5mm = round(mean(N_taxa_AMMOD_5mm,  na.rm = TRUE), 1),
    pct_retained   = round(mean(N_taxa_AMMOD_5mm,  na.rm = TRUE) /
                             mean(N_taxa_AMMOD_3mm,  na.rm = TRUE) * 100, 1)
  )

tableS <- dplyr::bind_rows(tableS_by_site, tableS_grand)

tableS_display <- tableS %>%
  dplyr::rename(
    `Site`                    = Site,
    `Habitat`                 = Habitat,
    `n days`                  = n_days,
    `Mean DS >= 3mm (%)`      = ds_3mm,
    `Mean DS >= 5mm (%)`      = ds_5mm,
    `Delta (pp)`              = delta_pp,
    `AMMOD taxa >= 3mm`       = ammod_taxa_3mm,
    `AMMOD taxa >= 5mm`       = ammod_taxa_5mm,
    `AMMOD taxa retained (%)` = pct_retained
  )

cat("--- Table S-X: Threshold sensitivity (3mm vs 5mm) ---\n\n")
print(tableS_display, n = Inf)

cat("\nInterpretation:\n")
cat(sprintf("  Switching from 3mm to 5mm filter on AMMOD:\n"))
cat(sprintf("  Grand mean similarity change:  %+.1f pp\n",  tableS_grand$delta_pp))
cat(sprintf("  AMMOD taxa retained at 5mm:    %.1f%%\n",    tableS_grand$pct_retained))
cat(sprintf("  AMMOD taxa discarded at 5mm:   %.1f%%\n\n",  100 - tableS_grand$pct_retained))

readr::write_csv(tableS_display, "outputs/H3_TableSX_threshold_sensitivity_3mm_vs_5mm.csv")
cat("[OK] outputs/H3_TableSX_threshold_sensitivity_3mm_vs_5mm.csv\n\n")

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📁 COMPLETE OUTPUT SUMMARY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("FIGURES (PNG, 300 DPI):\n")
cat("  1. H3_FAIRD_AMMOD_similarity_temporal_GENUS.png (10×8 in)\n")
cat("  2. H3_FAIRD_AMMOD_similarity_by_site_GENUS.png (8×6 in)\n\n")

cat("DATA FILES (CSV):\n")
cat("  1. H3_FAIRD_AMMOD_similarity_results_GENUS.csv (daily values, all three metrics)\n")
cat("  2. H3_FAIRD_AMMOD_similarity_summary_GENUS.csv (site summaries)\n")
cat("  3. H3_Table12_similarity_by_site.csv (Table 12: main H3 results, Delta_S_raw)\n")
cat("  4. H3_TableSX_threshold_sensitivity_3mm_vs_5mm.csv (Table S-X: sensitivity)\n")
cat("\n")

cat("✅ All files saved successfully\n")
cat("\n⚠️  NOTE: Analysis performed at GENUS level (FAIRD's primary resolution)\n")
cat("         AMMOD species collapsed to genera for fair comparison\n\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================

# ============================================================
# CLOSE LOG
# ============================================================
sink()
cat(">> Output saved to:", output_file, "\n")