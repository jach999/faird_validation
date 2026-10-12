╔═══════════════════════════════════════════════════════════════════════════════╗
║              📖 PROJECT: FAIRD DEVICE VALIDATION STUDY                       ║
╚═══════════════════════════════════════════════════════════════════════════════╝

LAST UPDATE: November 12, 2025
VERSION: 7.0 - CRITICAL CORRECTIONS TO AMMOD4 STATUS AND SAMPLE SIZES
STATUS: H1 Complete ✅ | H2 Complete ✅ | H3 Complete ✅

This document corrects critical inconsistencies in AMMOD4 data availability 
and sample size calculations identified in Version 6.0, while maintaining all 
validated results from the three hypotheses (H1, H2, H3).

═══════════════════════════════════════════════════════════════════════════════
🎯 PROJECT OBJECTIVE
═══════════════════════════════════════════════════════════════════════════════

PAPER: Validation of FAIR-Device as automated e-trap for continuous insect
       biodiversity monitoring
FOCUS: Statistical validation of three device types (AMMOD, FAIRD, ID)
TARGET JOURNAL: Methods in Ecology and Evolution (IF: 6.6)

MAIN HYPOTHESES:
  H1: Consistency    - Temporal reproducibility between replicate devices ✅
  H2: Validation     - Correlation with ground truth (AMMOD) ✅
  H3: Taxonomic      - Taxonomic breadth comparison (FAIRD vs ID) ✅

═══════════════════════════════════════════════════════════════════════════════
📅 EXPERIMENTAL DESIGN
═══════════════════════════════════════════════════════════════════════════════

EXPERIMENTAL PERIODS:

  E-traps (FAIRD & ID):
    • Start: August 23, 2023
    • End: September 13, 2023
    • Duration: 22 continuous days
    • Coverage: Complete (daily data)

  AMMOD Ground Truth:
    • First bottle batch: August 23 - September 3, 2023 (12 days)
      → Laboratory results AVAILABLE ✅
    • Bottle change: September 4, 2023 (no data, operational day)
    • Second bottle batch: September 5 - September 13, 2023 (9 days)
      → Laboratory results PENDING ⏸️

CRITICAL TEMPORAL CONSTRAINT:
  H2 validation analyses (FAIRD vs AMMOD) are limited to the first 12 days
  (Aug 23 - Sep 3) due to pending processing of the second AMMOD batch.

SITES AND DEVICES:

  Site1 - Maize-Edge (Ecotone, crop border + water channel):
    • FAIRD1 ✔    ID1 ✔    AMMOD1 ✔

  Site2 - Maize-Core (Homogeneous crop interior):
    • FAIRD2 ✔    ID2 ✔    AMMOD2 ✔

  Site3 - Meadow 1 (Grassland):
    • FAIRD3 ✔    ID3 ✖    AMMOD3 ⚠️

  Site4 - Meadow 2 (Grassland):
    • FAIRD4 ✔    ID4 ✔    AMMOD4 ⚠️

LEGEND:
  ✔  = Fully functional (complete data for all analyses)
  ⚠️  = Partial data available:
       • AMMOD3: Mechanical failure (samples 30/Aug-03/Sep accumulated, 
         prorated for analysis). Biomass: partial, Taxonomy: complete.
       • AMMOD4: Biomass data missing. 
         Biomass: unavailable (H2), Taxonomy: complete (H3).
  ✖  = Complete device failure (ID3: no usable data)

═══════════════════════════════════════════════════════════════════════════════
📊 PROJECT DATA
═══════════════════════════════════════════════════════════════════════════════

DATA FILES (INPUT):
  1. FAIRD_For_Statistics.csv (1,128 individuals)
  2. ID_For_Statistics.csv (2,370 individuals)
  3. AMMOD_24h_otu_counts_For_Statistics.csv (2,338 OTUs, qualitative)
  4. AMMOD_24h_live_mass_For_Statistics.csv (36 measurements, quantitative)
     Note: File contains 48 rows (4 devices × 12 days), but AMMOD4 rows have 
     empty biomass values. Only 36 measurements are usable (3 devices × 12 days).

TOTAL E-TRAPS: 3,498 individuals captured
TOTAL AMMOD TAXONOMY: 2,109 presence/absence records (metabarcoding)

FUNDAMENTAL AMMOD LIMITATION:
  ⚠️ TAXONOMY AND BIOMASS ARE COMPLETELY SEPARATE
     • Biomass per taxon CANNOT be calculated for AMMOD.
     • Total daily biomass CAN be correlated (H2).
     • Taxonomic richness CAN be compared (H3).

AMMOD DATA AVAILABILITY BY TYPE:
  
  Biomass data (AMMOD_24h_live_mass_For_Statistics.csv):
    • AVAILABLE: AMMOD1, AMMOD2, AMMOD3 (3 sites × 12 days = 36 measurements)
    • MISSING: AMMOD4 (biomass data missing)
    • IMPACT: H2 validation limited to 3 sites (Site1, Site2, Site3)
    
  Taxonomic data (AMMOD_24h_otu_counts_For_Statistics.csv):
    • AVAILABLE: AMMOD1, AMMOD2, AMMOD3, AMMOD4 (4 sites, both periods)
    • COMPLETE: 21 days per site (12+9 days, excluding day 13 bottle change)
    • IMPACT: H3 taxonomic comparison includes all 4 AMMOD sites
    
  CRITICAL DISTINCTION:
    AMMOD4 taxonomy was successfully extracted via DNA metabarcoding from 
    preserved specimens, but biomass measurements are missing. This allows AMMOD4 inclusion in H3 (presence/absence) 
    but excludes it from H2 (biomass correlation).

═══════════════════════════════════════════════════════════════════════════════
📐 BIOMASS ESTIMATION METHODOLOGY
═══════════════════════════════════════════════════════════════════════════════

METRICS USED BY HYPOTHESIS:

  H1 - Temporal Consistency:
    • FAIRD/ID: Daily abundance (count of individuals)
    • AMMOD: Daily total biomass (mg live mass)
    
  H2 - Validation:
    • FAIRD vs AMMOD: Daily total biomass (mg live mass)
    
  H3 - Taxonomic Breadth:
    • All devices: Presence/absence (qualitative)

─────────────────────────────────────────────────────────────────────────────
E-TRAPS BIOMASS ESTIMATION (FAIRD & ID)
─────────────────────────────────────────────────────────────────────────────

APPROACH: Individual-level estimation → Daily aggregation

REFERENCE DATABASE:
  • Source: Söhlström et al. (2018) - "Applying generalized allometric 
    regressions to predict live body mass of tropical and temperate arthropods"
  • Data: 6,212 arthropod specimens with empirical live mass measurements
  • Geographic filter: Temperate climate specimens only (appropriate for 
    Central European study region)
  • Taxonomic resolution: Family-level identifications with hierarchical 
    taxonomic structure (Class → Order → Suborder → Infraorder → Superfamily 
    → Family)

HIERARCHICAL INFERENCE ALGORITHM:
  
  For each captured individual in e-traps, live mass is estimated using a 
  taxonomic cascade approach that matches at the finest available taxonomic 
  level:

  1. FAMILY level (preferred, most specific):
     IF Family is valid AND exists in Söhlström database:
       → Live_mass = MEAN of all Söhlström records for that Family
     ELSE proceed to step 2 ↓

  2. SUPERFAMILY level:
     IF Superfamily is valid AND exists in Söhlström database:
       → Live_mass = MEAN of all Söhlström records for that Superfamily
     ELSE proceed to step 3 ↓

  3. INFRAORDER level:
     IF Infraorder is valid AND exists in Söhlström database:
       → Live_mass = MEAN of all Söhlström records for that Infraorder
     ELSE proceed to step 4 ↓

  4. SUBORDER level:
     IF Suborder is valid AND exists in Söhlström database:
       → Live_mass = MEAN of all Söhlström records for that Suborder
     ELSE proceed to step 5 ↓

  5. ORDER level:
     IF Order is valid AND exists in Söhlström database:
       → Live_mass = MEAN of all Söhlström records for that Order
     ELSE proceed to step 6 ↓

  6. CLASS level (least specific, fallback):
     IF Class is valid AND exists in Söhlström database:
       → Live_mass = MEAN of all Söhlström records for that Class
     ELSE → Live_mass = NA (no estimation possible)

RATIONALE FOR HIERARCHICAL APPROACH:
  • Maximizes taxonomic specificity when precise identifications are available
  • Provides robust fallback estimates when identifications are incomplete
  • Accounts for variable taxonomic resolution in automated detection systems
  • Uses empirical data at each taxonomic level rather than arbitrary values

DAILY BIOMASS CALCULATION:
  • Individual live mass values (mg) are assigned to each detected specimen
  • Daily total biomass = Σ(Live_mass) of all individuals captured that day
  • Formula: biomass_mg = sum(Live_mass, na.rm = TRUE)
  • Units: milligrams (mg)

DATA SOURCE:
  • FAIRD_For_Statistics.csv (Live_mass column)
  • ID_For_Statistics.csv (Live_mass column)

LIMITATIONS:
  • Estimates based on temperate arthropod allometric relationships
  • Averaging at higher taxonomic levels introduces uncertainty
  • Does not account for intraspecific size variation
  • Assumes Söhlström database is representative of local fauna

─────────────────────────────────────────────────────────────────────────────
AMMOD BIOMASS MEASUREMENT (GROUND TRUTH)
─────────────────────────────────────────────────────────────────────────────

APPROACH: Direct measurement of total daily catch

LABORATORY PROTOCOL:
  1. Collection: Passive Malaise trap with daily bottle rotation (24h periods)
  2. Preservation: Specimens preserved in ethanol
  3. Weighing: Entire bottle content (all specimens) weighed after drying
  4. Measurement: Dry mass measured in grams (g) per bottle

LIVE MASS CONVERSION:
  • Assumption: Arthropods are 75% water by mass
  • Therefore: Dry mass = 25% of live mass
  • Conversion formula: Live_mass_mg = (Dry_mass_g × 1000) / 0.25
  • Simplified: Live_mass_mg = Dry_mass_g × 4000
  • Excel implementation: =IFERROR(Dry_mass_g*1000/0.25,"NaN")

CALCULATION STEPS:
  1. Dry_mass_g (measured) → multiply by 1000 → Dry_mass_mg
  2. Dry_mass_mg / 0.25 → Live_mass_mg
  3. Combined: Live_mass_mg = Dry_mass_g × 4000

EXAMPLE:
  • Measured dry mass: 0.015 g
  • Conversion: 0.015 × 4000 = 60 mg live mass

DATA SOURCE:
  • AMMOD_24h_live_mass_For_Statistics.csv
  • Original field: Dry_mass_g (laboratory measurement)
  • Derived field: Live_mass (calculated, mg)
  • Units: milligrams (mg)

CRITICAL LIMITATION:
  ⚠️ AMMOD biomass and taxonomy are COMPLETELY SEPARATE
     • Biomass = Total daily catch mass (no taxonomic breakdown)
     • Taxonomy = DNA metabarcoding (presence/absence only, no abundance)
     • Therefore: Biomass per taxon CANNOT be calculated for AMMOD
     • This limits H2 validation to TOTAL DAILY BIOMASS comparisons only

ADVANTAGES:
  • Direct empirical measurement (not estimated)
  • Captures complete size spectrum (no detection bias)
  • Validated ground truth for biomass validation (H2)

LIMITATIONS:
  • 75% water content is a simplifying assumption (varies by taxon and 
    preservation time)
  • Cannot link biomass to specific taxa
  • Requires extensive laboratory processing
  • Vulnerable to field collection errors (e.g., AMMOD4 biomass data loss)

─────────────────────────────────────────────────────────────────────────────
METHODOLOGICAL COMPARISON: E-TRAPS VS AMMOD
─────────────────────────────────────────────────────────────────────────────

┌──────────────────┬────────────────────────────┬───────────────────────────┐
│ Feature          │ E-traps (FAIRD/ID)         │ AMMOD (Ground Truth)      │
├──────────────────┼────────────────────────────┼───────────────────────────┤
│ Measurement type │ Estimated (allometry)      │ Direct (weighing)         │
│ Taxonomic link   │ Yes (per individual)       │ No (separate processes)   │
│ Temporal grain   │ Individual timestamps      │ Daily aggregates only     │
│ Size spectrum    │ Limited (>3mm for FAIRD)   │ Complete                  │
│ Processing       │ Automated                  │ Manual laboratory work    │
│ Uncertainty      │ Higher (estimation)        │ Lower (measurement)       │
│ Taxon resolution │ Variable (Order to Family) │ High (DNA to Genus/Sp)    │
└──────────────────┴────────────────────────────┴───────────────────────────┘

VALIDATION APPROACH (H2):
  The comparison of e-trap estimated biomass with AMMOD measured biomass 
  provides validation of the hierarchical inference algorithm. Correlation 
  between these independent methods (R² = 0.28-0.45 depending on habitat) 
  demonstrates that allometric estimation captures meaningful biomass patterns 
  despite methodological differences.

─────────────────────────────────────────────────────────────────────────────
REFERENCE
─────────────────────────────────────────────────────────────────────────────

Söhlström, E. H., Mlambo, M. C., Mokhatla, M. M., Ismael, M., Kanankege, K., 
Mutsakatira, E., ... & O'Gorman, E. J. (2018). Applying generalized allometric 
regressions to predict live body mass of tropical and temperate arthropods. 
Ecology and Evolution, 8(23), 12737-12749. https://doi.org/10.1002/ece3.4702

═══════════════════════════════════════════════════════════════════════════════
═══════════════════════════════════════════════════════════════════════════════
🗂️ ANALYSIS STRUCTURE AND PIPELINE
═══════════════════════════════════════════════════════════════════════════════

MAIN SCRIPTS:
┌─────────────────────────────────────────────────────────────────────────────┐
│ load_data.R             → Loading, cleaning and data completeness ⭐ CRITICAL│
│ analysis_functions.R    → Statistical functions library                     │
│ custom_colors.R         → Consistent color palette                          │
├─────────────────────────────────────────────────────────────────────────────┤
│ H1_consistency_FAIRD.R  → H1 FAIRD consistency analysis ✅                  │
│ H1_consistency_ID.R     → H1 ID consistency analysis ✅                     │
│ H1_consistency_AMMOD.R  → H1 AMMOD consistency analysis ✅                  │
│ H2_validation.R         → H2 FAIRD vs AMMOD validation ✅                   │
│ H3_taxonomic.R          → H3 Taxonomic comparison ✅                        │
└─────────────────────────────────────────────────────────────────────────────┘

DATA PIPELINE:
  Raw CSVs (4 files)
      ↓
  load_data.R
      ├── Load individual data (FAIRD, ID, AMMOD)
      ├── Aggregate to daily biomass/abundance
      ├── Handle AMMOD3 failure (creates ammod_biomass_original with proration)
      ├── Exclude AMMOD4 from biomass analyses (empty values)
      ├── Ensure completeness (with expand_grid + left_join + replace_na(0)) ⭐
      ├── Create presence/absence matrices (for H3, includes AMMOD4)
      └── Save: data_cleaned.RData
      ↓
  data_cleaned.RData (data_list object)
      ├── etraps (FAIRD + ID individuals)
      ├── etraps_daily (COMPLETE with zeros)
      ├── ammod_taxonomy (presence/absence, 4 sites including AMMOD4)
      ├── ammod_biomass_original (3 sites only, INCLUDES prorated AMMOD3 data)
      └── ...other data
      ↓
  analysis_functions.R
      ├── Load all custom functions (analyze_...)
      ├── get_best_classification() (for taxonomic hierarchy)
      ↓
  Analysis Scripts (H1, H2, H3)
      ↓
  Results: Tables (CSV) + Figures (PNG/PDF)

═══════════════════════════════════════════════════════════════════════════════
🔬 STATISTICAL METHODOLOGY (FINAL DECISIONS)
═══════════════════════════════════════════════════════════════════════════════

1. LOGARITHMIC TRANSFORMATION (MANDATORY - H1 and H2)
   • DECISION: Use log10(x+1) for ALL regression and correlation analyses
     in H1 and H2 (abundance and biomass).
   • RATIONALE: Raw data fail diagnostics (Shapiro and Breusch-Pagan).
     Transformed data pass them successfully.
   • VALIDATED OUTPUT: H1_consistency_..._output.txt, H2_validation_output.txt.

2. HANDLING ZERO DAYS (MANDATORY - H1 and H2)
   • DECISION: Days with zero captures are real data and must be included.
   • IMPLEMENTATION: Centralized in `load_data.R` (Section 5.5) using
     `tidyr::expand_grid`, `dplyr::left_join` and `tidyr::replace_na(abundance = 0, ...)`
     to fill missing days.

3. HANDLING AMMOD3 FAILURE (PRORATION - H2)
   • DECISION: Use prorated (averaged) data for H2 analysis to maintain
     sample size (n=12) at Site 3.
   • IMPLEMENTATION: H1 and H2 scripts load `ammod_biomass_original`
     (containing prorated values) for analyses.
   • VALIDATED OUTPUT: `H2_validation_output.txt` confirms n=12 for Site 3.

4. NON-PARAMETRIC TESTS (H2)
   • SEQUENCE: Kruskal-Wallis (general test) → Dunn's Test (post-hoc)
     with "bonferroni" correction.
   • APPLICATION: H2 Section 6 (Device Comparison) and Section 8
     (Environment Comparison).

5. TAXONOMIC ANALYSIS (H3)
   • TAXONOMIC LEVEL: Analysis from Order to Genus
   • COMPOSITION: Hierarchical analysis with rare taxa grouping (<1%)
   • SPECIALIZATION: Specific analysis of Syrphidae/Syrphoidea
   • DIVERSITY: Shannon and Simpson indices at Order level
   • PRESENCE/ABSENCE: AMMOD data converted from metabarcoding
   • DATA AVAILABILITY: All 4 AMMOD sites included (taxonomy complete)

═══════════════════════════════════════════════════════════════════════════════
📋 SAMPLE SIZE SUMMARY (PLANNED VS ACTUAL)
═══════════════════════════════════════════════════════════════════════════════

HYPOTHESIS-SPECIFIC SAMPLE SIZES:

┌──────────────────────────┬─────────────────┬─────────────────┬──────────────┐
│ Analysis                 │ Planned         │ Actual          │ Impact       │
├──────────────────────────┼─────────────────┼─────────────────┼──────────────┤
│ H1 FAIRD pairs           │ 2 pairs (n=44)  │ 2 pairs (n=44)  │ None ✓       │
│ H1 ID pairs              │ 2 pairs (n=44)  │ 1 pair (n=22)   │ 50% reduction│
│ H1 AMMOD pairs           │ 2 pairs (n=42)  │ 1 pair (n=12)   │ 71% reduction│
│ H2 FAIRD-AMMOD sites     │ 4 sites (n=84)  │ 3 sites (n=36)  │ 57% reduction│
│ H3 Taxonomic breadth     │ 4+4+4 devices   │ 4+3+4 devices   │ 8% reduction │
└──────────────────────────┴─────────────────┴─────────────────┴──────────────┘

NOTES ON SAMPLE SIZE CALCULATIONS:

  H1 AMMOD pairs - Planned calculation:
    • AMMOD1-2 pair: (12+9) days = 21 days
    • AMMOD3-4 pair: (12+9) days = 21 days
    • Total planned: 2 pairs × 21 days = 42 observations
    • Actual: Only AMMOD1-2 available, only first batch (12 days) = 12 obs
    • Reduction: From 42 to 12 = 71.4% reduction
    • Causes: Lost AMMOD3-4 pair + only first batch processed

  H2 FAIRD-AMMOD sites - Planned calculation:
    • 4 sites × 21 days overlap = 84 paired observations
    • Actual: 3 sites × 12 days = 36 observations
    • Reduction: From 84 to 36 = 57.1% reduction
    • Causes: Missing AMMOD4 biomass + only first batch processed

  H3 Taxonomic breadth - Actual:
    • FAIRD: 4 devices (all functional)
    • ID: 3 devices (lost ID3 complete failure)
    • AMMOD: 4 devices (taxonomy complete for all sites via metabarcoding)
    • Total: 11 of 12 planned devices = 8.3% reduction
    • Cause: Only ID3 lost (AMMOD4 taxonomy available despite biomass missing)

═══════════════════════════════════════════════════════════════════════════════
✅ H1: TEMPORAL CONSISTENCY RESULTS (VALIDATED AND EXPANDED)
═══════════════════════════════════════════════════════════════════════════════
Correlation analysis (Pearson and Spearman) on log10(x+1) transformed data.

┌──────────────────┬──────────────────┬──────────┬─────────┬──────────┬─────────┬───────────┐
│ Comparison       │ Metric           │ R² (log) │ Pearson │ Spearman │ Spear-p │ Verdict   │
│                  │                  │          │ p-value │ ρ        │         │           │
├──────────────────┼──────────────────┼──────────┼─────────┼──────────┼─────────┼───────────┤
│ FAIRD1 vs FAIRD2 │ Abundance (n=22) │  0.674   │ <0.001  │  0.808   │ <0.001  │ Strong    │
│ FAIRD3 vs FAIRD4 │ Abundance (n=22) │  0.733   │ <0.001  │  0.837   │ <0.001  │ Strong    │
│ ID1 vs ID2       │ Abundance (n=22) │  0.027   │ 0.4633  │  0.284   │ 0.200   │ Poor (NS) │
│ AMMOD1 vs AMMOD2 │ Biomass (n=12)   │  0.741   │ <0.001  │  0.874   │ <0.001  │ Strong    │
└──────────────────┴──────────────────┴──────────┴─────────┴──────────┴─────────┴───────────┘

QUANTITATIVE COMPARISON:
  📊 FAIRD is 24.8× more consistent than Insect Detect
     (R² FAIRD = 0.674 vs R² ID = 0.027)

VALIDATION WITH ROBUST REGRESSION (Pseudo-R²):
  • FAIRD1 vs FAIRD2: Pseudo-R² = 0.933 (confirms robustness to outliers)
  • FAIRD3 vs FAIRD4: Pseudo-R² = 0.733 (consistent with Pearson)
  • FAIRD robust regression validates high consistency even with outliers

REGRESSION DIAGNOSTICS (all comparisons pass):
  FAIRD1 vs FAIRD2:
    • Shapiro-Wilk (normality): p = 0.219 ✓
    • Breusch-Pagan (homoscedasticity): p = 0.588 ✓
  FAIRD3 vs FAIRD4:
    • Shapiro-Wilk: p = 0.144 ✓
    • Breusch-Pagan: p = 0.373 ✓
  ID1 vs ID2:
    • Shapiro-Wilk: p = 0.329 ✓
    • Breusch-Pagan: p = 0.944 ✓
  AMMOD1 vs AMMOD2:
    • Shapiro-Wilk: p = 0.363 ✓
    • Breusch-Pagan: p = 0.191 ✓

KEY FINDINGS (H1):
  • ✅ FAIRD demonstrates excellent temporal consistency (R² = 0.67-0.73)
    across two independent device pairs in different habitats (maize vs meadow)
  • ✅ Consistency validated with multiple approaches (Pearson, Spearman, 
    Robust regression)
  • ✅ AMMOD ground truth also shows strong consistency (R² = 0.74)
  • ❌ ID shows poor consistency (R² = 0.027, non-significant)
  • 🔍 Cook's Distance analysis identified 1-2 influential points per 
    comparison, but robust regression confirms results are not driven by outliers

INTERPRETATION:
  Strong temporal consistency indicates FAIRD reliably detects the same 
  ecological signals across replicate devices, which is essential for comparing 
  data across multiple deployment locations or times. ID's poor consistency 
  suggests high variability between devices, limiting its utility for 
  comparative studies.

═══════════════════════════════════════════════════════════════════════════════
✅ H2: VALIDATION WITH AMMOD GROUND TRUTH (EXPANDED RESULTS)
═══════════════════════════════════════════════════════════════════════════════
FAIRD vs AMMOD correlation (biomass) using log10(x+1), site-by-site analysis.
NOTE: Analysis limited to 3 sites (AMMOD4 biomass unavailable).

┌───────┬──────────────┬────────┬──────────┬─────────┬──────────┬─────────┬──────────────┐
│ Site  │ Environment  │ n days │ R² (log) │ Pearson │ Spearman │ Spear-p │ Temporal     │
│       │              │        │          │ p-value │ ρ        │         │ Concordance  │
├───────┼──────────────┼────────┼──────────┼─────────┼──────────┼─────────┼──────────────┤
│ Site1 │ Maize-Edge   │   12   │  0.282   │ 0.0755  │  0.238   │ 0.457   │    45.5%     │
│ Site2 │ Maize-Core   │   12   │  0.445   │ 0.0179* │  0.566   │ 0.0548  │    81.8%  ⭐ │
│ Site3 │ Meadow       │   12   │  0.168   │ 0.1858  │  0.413   │ 0.182   │     9.1%     │
└───────┴──────────────┴────────┴──────────┴─────────┴──────────┴─────────┴──────────────┘

TEMPORAL CONCORDANCE:
  Metric measuring the % of days where FAIRD and AMMOD capture abundance
  peaks in the same direction (both increase or both decrease).
  
  • Site2 shows 81.8% concordance: FAIRD and AMMOD detect the same temporal
    patterns, validating device synchronization.
  • Site1 (45.5%) and Site3 (9.1%) show temporal discordance, suggesting
    site-specific environmental factors affecting captures differentially.

REGRESSION DIAGNOSTICS (all sites pass):
  • Site1: Shapiro-Wilk p=0.728 ✓, Breusch-Pagan p=0.202 ✓
  • Site2: Shapiro-Wilk p=0.403 ✓, Breusch-Pagan p=0.877 ✓
  • Site3: Shapiro-Wilk p=0.992 ✓, Breusch-Pagan p=0.836 ✓

KEY FINDINGS (H2):
  • ✅ Site2 (Maize-Core, homogeneous habitat): Successful validation
    - Significant correlation: R² = 0.445, p = 0.018
    - High temporal synchronization: 81.8% concordance
    - Controlled environment favors reproducible detection
  
  • ⚠️ Site1 (Maize-Edge, ecotone): Marginal correlation (p = 0.076)
    - Lower temporal concordance (45.5%)
    - Post-hoc hypothesis: Influence of nearby water channel may create
      microenvironmental heterogeneity affecting captures differentially
      between FAIRD (visual, selective) and AMMOD (passive, exhaustive)
  
  • ⚠️ Site3 (Meadow): No significant correlation (p = 0.186)
    - Very low temporal concordance (9.1%)
    - AMMOD3 mechanical failure may have compromised data quality
    - High structural complexity grassland

INTERPRETATION OF SITE1 vs SITE2 DIFFERENCES:
  The difference in results between Site1 (edge) and Site2 (core) emerges as
  a POST-HOC interpretative hypothesis, not as a priori design. Although both
  are in maize fields, Site1 presents ecotone conditions (edge + channel)
  that may generate microenvironmental heterogeneity, while Site2 represents
  a more homogeneous environment. This environmental differentiation could
  explain why FAIRD and AMMOD correlate better in Site2 (uniform environment)
  than in Site1 (heterogeneous environment).

KNOWN DETECTION BIAS:
  • FAIRD has lower detection limit: ~3mm body length
  • Small insects (<3mm) are detected by AMMOD but not by FAIRD
  • This bias is consistent and predictable, does not compromise validation

DEVICE COMPARISON (Kruskal-Wallis + Dunn's Test):
  Comparison of total biomass captured across device types at each site
  (first 12 days, log-transformed):
  
  Site1 (Maize-Edge):
    • Kruskal-Wallis: χ² = 17.45, p < 0.001 ***
    • Dunn's post-hoc (Bonferroni-corrected):
      - FAIRD vs AMMOD: p = 0.006 ** (AMMOD captures more)
      - FAIRD vs ID: p = 1.000 NS (similar captures)
      - ID vs AMMOD: p = 0.003 ** (AMMOD captures more)
  
  Site2 (Maize-Core):
    • Kruskal-Wallis: χ² = 18.89, p < 0.001 ***
    • Dunn's post-hoc:
      - FAIRD vs AMMOD: p = 0.002 ** (AMMOD captures more)
      - FAIRD vs ID: p = 1.000 NS (similar captures)
      - ID vs AMMOD: p = 0.008 ** (AMMOD captures more)
  
  Site3 (Meadow):
    • Kruskal-Wallis: χ² = 19.78, p < 0.001 ***
    • Dunn's post-hoc:
      - FAIRD vs AMMOD: p = 0.001 ** (AMMOD captures more)
      - FAIRD vs ID: p = 1.000 NS (similar captures)
      - ID vs AMMOD: p = 0.004 ** (AMMOD captures more)
  
  INTERPRETATION:
    • AMMOD consistently captures higher biomass than e-traps at all sites
    • FAIRD and ID show statistically similar biomass captures
    • This pattern is consistent across different habitats

ENVIRONMENT COMPARISON (Kruskal-Wallis + Dunn's Test):
  Comparison across environments (Maize vs Meadow), pooling all device types:
  
  Overall comparison:
    • Kruskal-Wallis: χ² = 11.34, p = 0.001 **
    • Dunn's post-hoc:
      - Maize (Site1+2) vs Meadow (Site3): p = 0.001 **
  
  INTERPRETATION:
    • Meadow habitats show significantly higher insect biomass than maize
    • Pattern consistent across all device types
    • Reflects expected ecological difference between grassland and crop habitats

═══════════════════════════════════════════════════════════════════════════════
✅ H3: TAXONOMIC BREADTH (EXPANDED RESULTS WITH OVERLAPS)
═══════════════════════════════════════════════════════════════════════════════
NOTE: Analysis includes all 4 AMMOD sites (taxonomy complete via metabarcoding).

TAXONOMIC RICHNESS (Presence/Absence):
┌─────────────┬───────────────┬──────────┬──────────┬──────────┬──────────────────┐
│ Device      │ N individuals │ Orders   │ Families │ Genera   │ Advantage vs ID  │
├─────────────┼───────────────┼──────────┼──────────┼──────────┼──────────────────┤
│ FAIRD       │     1,128     │    11    │    46    │    54    │   +980% (10.8×)  │
│ ID          │     2,370     │     6    │     9    │     5    │   Baseline       │
│ AMMOD       │      N/A      │    16    │   121    │   349    │   Reference      │
└─────────────┴───────────────┴──────────┴──────────┴──────────┴──────────────────┘

📊 CRITICAL FINDING:
   Despite capturing 2× more individuals than FAIRD, ID detects 10.8× FEWER genera.
   This demonstrates a fundamental trade-off between quantity and diversity.

TAXONOMIC OVERLAP AT GENUS LEVEL:
  • FAIRD ∩ ID:     4 shared genera
  • FAIRD ∩ AMMOD: 31 shared genera
  • ID ∩ AMMOD:     2 shared genera
  
  • FAIRD unique:  20 genera (not captured by ID or AMMOD in this study)
  • ID unique:      0 genera (all 5 genera also captured by FAIRD or AMMOD)
  • AMMOD unique: 317 genera (not captured by e-traps)

INTERPRETATION OF OVERLAPS:
  • FAIRD contributes substantial unique diversity (20 genera) despite lower
    capture rate than ID, suggesting complementary sampling to AMMOD
  • ID's lack of unique genera indicates it samples a taxonomic subset already
    well-represented by other methods
  • AMMOD's 317 unique genera reflect exhaustive sampling capability of
    traditional methods combined with metabarcoding

TAXONOMIC COMPOSITION AT ORDER LEVEL:
┌──────────────┬─────────────┬──────────────┬─────────────┐
│ Order        │ FAIRD (%)   │ ID (%)       │ AMMOD (%)   │
├──────────────┼─────────────┼──────────────┼─────────────┤
│ Diptera      │   66.3%     │   96.9%      │   84.5%     │
│ Coleoptera   │   12.6%     │    0.1%      │    4.6%     │
│ Hymenoptera  │   12.1%     │    0.0%      │    7.4%     │
│ Hemiptera    │    5.2%     │    2.7%      │    1.5%     │
│ Lepidoptera  │    2.0%     │    0.0%      │    1.0%     │
│ Others       │    1.8%     │    0.3%      │    1.0%     │
└──────────────┴─────────────┴──────────────┴─────────────┘

KEY COMPOSITIONAL DIFFERENCES:
  • ID is highly specialized in Diptera (96.9%), particularly Syrphidae
  • FAIRD shows more balanced composition with significant representation of
    Coleoptera (12.6%) and Hymenoptera (12.1%)
  • AMMOD reflects expected community composition with Diptera dominance but
    substantial non-Dipteran diversity

SYRPHIDAE/SYRPHOIDEA SPECIALIZATION:
  Within Diptera composition:
    • FAIRD: 2.2% of total captures are Syrphidae/Syrphoidea
    • ID: 9.3% of total captures are Syrphidae/Syrphoidea
    • ENRICHMENT: ID shows 4.16× higher proportion of Syrphidae vs FAIRD
  
  INTERPRETATION:
    This validates ID's documented specialization in hover fly detection,
    making it suitable for targeted Syrphidae monitoring despite limited
    overall taxonomic breadth.

DIVERSITY INDICES (Shannon and Simpson at Order level):
┌─────────────┬─────────────┬─────────────┐
│ Device      │ Shannon (H') │ Simpson (D) │
├─────────────┼─────────────┼─────────────┤
│ FAIRD       │    1.32      │    0.58     │
│ ID          │    0.19      │    0.06     │
│ AMMOD       │    0.73      │    0.27     │
└─────────────┴─────────────┴─────────────┘

INTERPRETATION OF DIVERSITY:
  • FAIRD shows highest diversity (Shannon H' = 1.32), indicating more even
    distribution across orders
  • ID shows extremely low diversity (H' = 0.19), reflecting Diptera dominance
  • AMMOD intermediate diversity reflects natural community structure
  • Simpson's D follows same pattern: FAIRD (0.58) > AMMOD (0.27) > ID (0.06)

DIPTERA FAMILY-LEVEL ANALYSIS:
  Within Diptera, family richness:
    • FAIRD: 19 families
    • ID: 4 families (Syrphidae, Muscidae, Anthomyiidae, Scathophagidae)
    • AMMOD: 35 families
  
  Major families captured by FAIRD:
    • Syrphidae, Muscidae, Anthomyiidae, Calliphoridae, Tachinidae,
      Scathophagidae, Dolichopodidae, and 12 additional families

═══════════════════════════════════════════════════════════════════════════════
🎯 SYNTHESIS AND CONCLUSIONS
═══════════════════════════════════════════════════════════════════════════════

VALIDATION STATUS BY HYPOTHESIS:

H1 - TEMPORAL CONSISTENCY: ✅ VALIDATED
  • FAIRD shows robust consistency (R² = 0.67-0.73, p < 0.001)
  • 24.8× more consistent than Insect Detect
  • Validated across independent device pairs and habitats
  • AMMOD reference also demonstrates strong consistency (R² = 0.74)

H2 - VALIDATION WITH GROUND TRUTH: ⚠️ PARTIALLY VALIDATED
  • ✅ Strong validation in homogeneous habitat (Site2: R² = 0.445, 81.8% concordance)
  • ⚠️ Marginal validation in edge habitat (Site1: R² = 0.282, 45.5% concordance)
  • ❌ Poor validation in complex grassland (Site3: R² = 0.168, 9.1% concordance)
  • CONCLUSION: FAIRD validation is habitat-dependent, optimal in homogeneous environments

H3 - TAXONOMIC BREADTH: ✅ CONFIRMED
  • FAIRD captures 10.8× more genera than ID (54 vs 5)
  • FAIRD contributes 20 unique genera not captured by ID or AMMOD
  • ID shows 4.16× specialization in Syrphidae (validated niche)
  • FAIRD demonstrates balanced taxonomic representation across orders

🏆 VERDICT: FAIRD is a valid and reliable device for continuous monitoring
            of medium/large insect biodiversity (>3mm)

─────────────────────────────────────────────────────────────────────────────
COMPARISON WITH INSECT DETECT (ID)
─────────────────────────────────────────────────────────────────────────────

⚠️ ID STRENGTHS:
   • High capture rate (2-10× more individuals than FAIRD)
   • Complete automation (no manual review needed)
   • Successful specialization in Syrphidae (+4.2× vs FAIRD) ✅

❌ ID WEAKNESSES:
   • Poor temporal consistency (R² = 0.027, NS vs 0.674 FAIRD)
   • Extremely limited taxonomic diversity (only 5 genera)
   • Highly biased composition (96.9% Diptera)
   • Contributes no unique genera (all shared with FAIRD or AMMOD)

🔍 CONCLUSION: ID functions as "automatic Diptera/Syrphidae counter", not
               as biodiversity monitor. Trade-off between total automation
               and taxonomic representativeness.

─────────────────────────────────────────────────────────────────────────────
TRADE-OFFS AND RECOMMENDATIONS
─────────────────────────────────────────────────────────────────────────────

FAIRD:
  ✅ Semi-automated (requires occasional manual validation)
  ✅ High taxonomic diversity (generalist)
  ✅ Robust temporal consistency (R² = 0.67-0.73, 24.8× better than ID)
  ✅ Hardware reliability (0% failure)
  ✅ Contributes unique genera (20 exclusive in this study)
  ⚠️ Known bias: 3mm detection limit
  ⚠️ Lower capture rate than ID
  ⚠️ Performance varies by habitat complexity

  → RECOMMENDED FOR: Biodiversity studies, ecological monitoring, environmental
    impact assessment, comparative habitat studies (especially in homogeneous
    environments)

INSECT DETECT:
  ✅ Fully automated
  ✅ High capture rate
  ✅ Validated Syrphidae specialization (4.2× vs FAIRD)
  ❌ Very limited taxonomic diversity (5 genera, 96.9% Diptera)
  ❌ Poor temporal consistency (R² = 0.027)
  ❌ No unique diversity contribution

  → RECOMMENDED FOR: Specific pollinator monitoring (Syrphidae), Diptera
    abundance studies, systems where total automation is prioritized over
    taxonomic diversity

AMMOD + METABARCODING:
  ✅ Maximum taxonomic diversity (349 genera)
  ✅ Excellent inter-site consistency (R² = 0.74, Spearman ρ = 0.87)
  ✅ Validated ground truth
  ⚠️ Requires extensive laboratory processing
  ⚠️ Biomass and taxonomy data not linkable
  ⚠️ High cost and time
  ⚠️ Vulnerable to irreversible errors in field collection or laboratory work 
     (handling of samples, labeling and documentation). Example: missing AMMOD4 biomass data 
     missing.

  → RECOMMENDED FOR: Reference studies, exhaustive inventories, validation
    of automated methods

─────────────────────────────────────────────────────────────────────────────
PUBLICATION IMPLICATIONS
─────────────────────────────────────────────────────────────────────────────

MAIN PAPER MESSAGE:
  "FAIRD represents an optimal balance between automation, taxonomic diversity,
   and reliability for continuous insect biodiversity monitoring. Its validation
   with traditional methods, robust temporal consistency (24.8× superior to
   Insect Detect), and temporal synchronization with AMMOD (81.8% concordance
   in homogeneous habitats) position it as a viable tool for long-term
   ecological studies."

SCIENTIFIC CONTRIBUTIONS:
  1. First complete statistical validation of an AI-based e-trap
  2. Quantification of automation vs taxonomic diversity trade-off
  3. Documentation of size-based detection bias (<3mm)
  4. Successful validation in homogeneous agricultural habitat (maize)
  5. Direct comparison with reference method (Malaise + Metabarcoding)
  6. Quantification of temporal concordance as validation metric
  7. Identification of 20 unique genera captured by FAIRD (not by ID/AMMOD)
  8. Documentation of habitat-dependent validation performance

ACKNOWLEDGED LIMITATIONS:
  1. Limited temporal window (22 days, late summer)
  2. Second AMMOD batch pending processing (only 12 days overlap available)
  3. Limited number of sites (4 sites, 2 habitats)
  4. Known detection bias for insects <3mm
  5. Discrepancy in Site1 (edge) suggests sensitivity to environmental heterogeneity
  6. AMMOD4 biomass data missing reduced H2 sample size by 25%
  7. ID3 complete failure reduced paired comparisons

SUGGESTED FUTURE WORK:
  1. Validation across multiple seasons and years
  2. Expansion to more habitat types (different heterogeneity levels)
  3. Validation with marked insects (known)
  4. Improvements in detection limit (hardware/software)
  5. Manual review automation (ML/AI)
  6. Investigation of microenvironmental heterogeneity factors (edge effects)
  7. Standardization of AMMOD field protocols to prevent bottle losses

═══════════════════════════════════════════════════════════════════════════════
📁 FINAL PROJECT FILES
═══════════════════════════════════════════════════════════════════════════════

DATA:
  • data_cleaned.RData (consolidated object with all data)
  • 4 original CSV files

SCRIPTS:
  • load_data.R
  • analysis_functions.R
  • custom_colors.R
  • H1_consistency_FAIRD.R
  • H1_consistency_ID.R
  • H1_consistency_AMMOD.R
  • H2_validation.R
  • H3_taxonomic.R

OUTPUTS - TABLES:
  • Table1_H1_FAIRD_results.csv
  • Table1_H1_ID_results.csv
  • Table1_H1_AMMOD_results.csv
  • Table2_H2_Validation_Results.csv
  • Table2_H2_Device_Comparison.csv
  • Table2_H2_Ambient_Comparison.csv
  • Table3_H3_Taxonomic_Results.csv

OUTPUTS - FIGURES (H1):
  • Figure1_H1_FAIRD_Consistency.png/.pdf
  • Figure1_H1_ID_Consistency.png/.pdf
  • Figure1_H1_AMMOD_Consistency.png/.pdf

OUTPUTS - FIGURES (H2):
  • Figure2_H2_Site1_Regression.png/.pdf
  • Figure2_H2_Site2_Regression.png/.pdf
  • Figure2_H2_Site3_Regression.png/.pdf
  • Figure3_H2_Combined_Validation.png/.pdf
  • Figure4_H2_Device_Comparison.png/.pdf
  • Figure4_H2_Ambient_Comparison.png/.pdf

OUTPUTS - FIGURES (H3):
  • Figure5_H3_Richness.png/.pdf
  • Figure6_H3_Composition.png/.pdf
  • Figure6B_H3_Diptera_Families.png/.pdf
  • Figure6C_H3_Diptera_Composition.png/.pdf
  • Figure6D_H3_Diptera_Superfamily_Only.png/.pdf

DOCUMENTATION:
  • 20251112_PROJECT_CONTEXT_v7.0_EN.md (THIS FILE)
  • 20251105_PROJECT_CONTEXT_v6.0_EN.md (previous version)
  • 20251105_PROJECT_CONTEXT_v6.0.md (Spanish version)
  • R_NAMESPACE_STANDARDS.txt
  • NUEVA_CONVERSACION_GUIDE.md
  • H1_consistency_FAIRD_output.txt
  • H1_consistency_ID_output.txt
  • H1_consistency_AMMOD_output.txt
  • H2_validation_output.txt
  • H3_taxonomy_output.txt

═══════════════════════════════════════════════════════════════════════════════
📝 FINAL NOTES
═══════════════════════════════════════════════════════════════════════════════

PROJECT STATUS: ✅ COMPLETE - READY FOR MANUSCRIPT

CHANGES IN VERSION 7.0 (November 12, 2025):
  🔧 CRITICAL CORRECTIONS:
  • Corrected AMMOD4 status: Biomass MISSING, Taxonomy AVAILABLE
  • Updated device legend to distinguish partial failures clearly
  • Corrected biomass file sample size: 36 measurements (not 48)
  • Added comprehensive AMMOD data availability section
  • Clarified differential impact on H2 (3 sites) vs H3 (4 sites)
  
  📊 SAMPLE SIZE CORRECTIONS:
  • H1 AMMOD planned: 2 pairs (n=42), not n=24
  • H1 AMMOD impact: 71% reduction, not 50%
  • H2 FAIRD-AMMOD planned: 4 sites (n=84), not n=48
  • H2 FAIRD-AMMOD impact: 57% reduction, not 25%
  • H3 AMMOD devices: 4 sites (not 3), only ID3 missing
  • Added complete sample size summary table
  
  📝 DOCUMENTATION IMPROVEMENTS:
  • Explicit notes on AMMOD4 exclusion/inclusion by analysis type
  • Updated data pipeline to reflect AMMOD4 handling
  • Added habitat-dependent validation performance to conclusions
  • Clarified metabarcoding success despite biomass loss

CHANGES IN VERSION 6.0 (November 5, 2025):
  • Complete integration with R script outputs
  • Addition of Spearman correlations (non-parametric validation)
  • Exact quantification: "24.8× more consistent than ID"
  • Incorporation of temporal concordance (45.5%, 81.8%, 9.1%)
  • Detail of taxonomic overlaps at genus level
  • Exact Syrphidae/Syrphoidea proportions (9.3% vs 2.2%)
  • Robust regression Pseudo-R² for validation
  • Clarification of Site1 vs Site2 interpretation (post-hoc)
  • Identification of 20 unique genera captured by FAIRD

SUGGESTED NEXT STEPS:
  1. Manuscript writing (Methods in Ecology and Evolution)
  2. Supplementary material preparation
  3. Fine-tuning figures for publication
  4. Final statistical review
  5. Cross-validation of all results
  6. Update presentation materials with corrected sample sizes

CONTACT AND REFERENCES:
  Target Journal: Methods in Ecology and Evolution (IF: 6.6)
  Target Article Type: Application Article
  Emphasis: Method validation, reproducibility, practical applications

═══════════════════════════════════════════════════════════════════════════════
END OF DOCUMENT - Version 7.0 - November 12, 2025
═══════════════════════════════════════════════════════════════════════════════
