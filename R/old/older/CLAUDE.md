# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based statistical analysis project for validating the FAIR-Device (FAIRD) automated e-trap for continuous insect population monitoring. The project compares three device types (FAIRD, Insect Detect, AMMOD) across 4 field sites over a 22-day period (Aug 23 - Sep 13, 2023).

**Key Context Document**: Always read `PROJECT_CONTEXT_2025_10_28_INTEGRATED.txt` before starting work - it contains critical methodological decisions, established standards, and project state.

## Running Analyses

### Complete Analysis Pipeline

```r
# Run all three hypotheses in sequence
source("run_analysis.R")
```

This executes H1 (consistency) → H2 (validation) → H3 (taxonomic) and organizes outputs into `results/` directory.

### Individual Hypothesis Scripts

```r
# H1: Intra-device consistency (FAIRD vs FAIRD, ID vs ID)
source("H1_consistency.R")

# H2: FAIRD-AMMOD validation (correlations, temporal synchronization)
source("H2_validation.R")

# H3: Taxonomic breadth comparison (FAIRD vs ID)
source("H3_taxonomic.R")
```

Each script requires `load_data.R` to be run first (or loads it automatically).

### Data Loading

```r
# Load and clean all data (required before running analyses)
source("load_data.R")
load("data_cleaned.RData")

# Access datasets via list
data_list$etraps_daily           # E-traps daily biomass (COMPLETE with zeros)
data_list$ammod_biomass          # AMMOD biomass (for stats, excludes AMMOD3 malfunction period)
data_list$ammod_biomass_original # AMMOD biomass (includes prorated AMMOD3 values)
data_list$faird                  # FAIRD individual-level data
data_list$id                     # Insect Detect individual-level data
```

## Architecture

### Data Pipeline Flow

```
Raw CSVs (4 files)
    ↓
load_data.R
    ├── Individual data (faird_data, id_data, etraps_data)
    ├── Daily aggregation (etraps_daily_biomass)
    ├── Data completeness enforcement (tidyr::complete) ⭐ CRITICAL
    ├── AMMOD taxonomy (presence/absence only)
    ├── AMMOD biomass (two versions: stats vs original)
    └── Presence/absence matrices
    ↓
data_cleaned.RData (data_list object)
    ↓
H1/H2/H3 analysis scripts
    ↓
Tables (CSV) + Figures (PNG/PDF)
```

### Critical Data Processing (load_data.R:244-295)

The **data completeness section** (added Oct 28, 2025) is essential:

- Uses `tidyr::complete()` to create a master grid of all device×date combinations
- Fills missing days with `abundance = 0`, `biomass_mg = 0`
- This ensures all analyses include biologically meaningful zero-capture days
- **NEVER remove or modify this section** - it's foundational for H1, H2, H3

### AMMOD3 Mechanical Failure Handling

AMMOD3 experienced mechanical failure on Aug 30, 2023:

- **Two versions created**:
  - `ammod_biomass_for_stats`: Excludes malfunction period (Aug 30 - Sep 3) - **use for statistical tests**
  - `ammod_biomass_original`: Includes prorated values - **use for H2 visualization** (see `H2_validation.R:40`)
- Default `data_list$ammod_biomass` points to stats version
- **Justification**: Prorated values (total/5 days) assume uniform capture rate - documented assumption in manuscript

## Established Methodological Standards

### 1. Log Transformation (FINALIZED Oct 28, 2025)

**Decision**: Use `log₁₀(x + 1)` transformation for **ALL biomass/abundance analyses**

```r
# Standard transformation pattern
data_log <- data %>%
  mutate(
    biomass_log = log10(biomass_mg + 1),
    abundance_log = log10(abundance + 1)
  )

model <- lm(biomass_log ~ predictor_log, data = data_log)
```

**Justification**:
- Raw biomass is right-skewed and heteroscedastic (fails Shapiro-Wilk, Breusch-Pagan)
- `log(x+1)` preserves zeros (log(0+1) = 0)
- Stabilizes variance and meets regression assumptions
- Applied to **H1** (consistency) and **H2** (validation)

**Validation**: All transformed analyses pass Shapiro-Wilk (p > 0.05) and Breusch-Pagan tests

### 2. Zero-Days Handling (FINALIZED Oct 28, 2025)

**Decision**: **ALWAYS include** zero-capture days as real observations

**Implementation**: Centralized in `load_data.R` via `tidyr::complete()` - benefits all scripts

**Biological Rationale**:
- Zeros represent real conditions (adverse weather, low activity periods)
- NOT missing data or technical failures
- FAIRD captures ~10× more than AMMOD → FAIRD zeros align with AMMOD low-biomass days

**Impact**: Including 1 zero-day at Site1 improved R² by 61.7% (0.075 → 0.121)

### 3. Namespace Standards (MANDATORY)

**ALWAYS use explicit namespace for tidyverse functions** to avoid conflicts with MASS package:

```r
# CORRECT
dplyr::filter()
dplyr::select()
dplyr::mutate()
tidyr::complete()
tidyr::pivot_wider()
readr::read_csv()

# WRONG
filter()   # Conflicts with MASS::filter() if MASS loaded after tidyverse
select()   # Conflicts with MASS::select()
```

**Package loading order**:
```r
library(tidyverse)   # Load FIRST
library(car)
library(MASS)        # Load LAST (highest conflict risk)
```

### 4. Non-Parametric Tests

Use **Kruskal-Wallis + Dunn's post-hoc** for Device×Site comparisons:

```r
kruskal.test(abundance ~ DevicexSite, data = data)
dunn.test::dunn.test(data$abundance, data$DevicexSite, method = "bonferroni")
```

**Effect size**: Report epsilon-squared (ε²) via `rstatix::epsilon_squared()`

### 5. Results Must Be Dynamic (NO Hardcoding)

```r
# CORRECT - extract from model
r2 <- summary(model)$r.squared
p_val <- summary(model)$coefficients[2, 4]

# WRONG - hardcoded values
annotate(label = "R² = 0.61, p = 0.003")  # ❌ Will become stale
```

## Key Results (Oct 28, 2025 Status)

**H1 - Consistency** ✅ COMPLETE:
- FAIRD Maize: R² = 0.67***, FAIRD Meadow: R² = 0.73*** (log-transformed)
- ID Maize: R² = 0.03 (not significant)
- **FAIRD is 3× more consistent than ID**

**H2 - Validation** ✅ COMPLETE:
- Site2 (Maize-Core): R² = 0.45, p = 0.018* (log-transformed, best validation)
- Site1 (Maize-Edge): R² = 0.28, p = 0.076 (weak, edge effects)
- Site3 (Meadow): R² = 0.17, p = 0.186 (very weak, only n=7 days)

**H3 - Taxonomic** ⏸️ PENDING:
- FAIRD: 55 genera (generalist)
- ID: 5 genera (Diptera specialist)
- **980% more genera captured by FAIRD**

## File Structure

### Core Scripts (Root Directory)
- `load_data.R` - Data loading + cleaning + completeness enforcement ⭐
- `run_analysis.R` - Master pipeline
- `analysis_functions.R` - Shared statistical functions
- `custom_colors.R` - Consistent color palette
- `H1_consistency.R` - Hypothesis 1 analysis
- `H2_validation.R` - Hypothesis 2 analysis
- `H3_taxonomic.R` - Hypothesis 3 analysis

### Input Data (Root Directory)
- `FAIRD_For_Statistics.csv` - FAIRD individual-level data
- `ID_For_Statistics.csv` - Insect Detect individual-level data
- `AMMOD_24h_otu_counts_For_Statistics.csv` - AMMOD taxonomy (presence/absence)
- `AMMOD_24h_live_mass_For_Statistics.csv` - AMMOD biomass (gravimetric)

### Key Outputs
- `data_cleaned.RData` - Master data object (load with `load()`)
- `Table1_H1_Consistency_Results.csv` - H1 correlation results
- `Table2_H2_Validation_Results.csv` - H2 validation metrics
- `Figure_S1_Diagnostic_Comparison.png` - Original vs log diagnostics
- `results/` - Organized outputs from `run_analysis.R`

### Documentation (Root Directory)
- `PROJECT_CONTEXT_2025_10_28_INTEGRATED.txt` ⭐ **READ THIS FIRST**
- `executive_summary_updated.md` - Detailed experimental design

## Critical Data Constraints

### AMMOD Data Is NOT Linked
- **Taxonomy**: Metabarcoding (16 orders, 349 genera) - presence/absence ONLY
- **Biomass**: Gravimetric (Live_mass column) - quantitative, NO taxonomy
- **Cannot connect taxonomy to biomass** for AMMOD
- Use AMMOD taxonomy for H3 (presence/absence comparisons only)
- Use AMMOD biomass for H2 (validation correlations only)

### Device Failures
- **ID3**: Complete failure (no data) - excluded from analyses
- **AMMOD4**: Complete failure (all NaN) - Site4 cannot use AMMOD
- **AMMOD3**: Partial failure (Aug 30 onwards) - see handling strategy above

### Valid Sites for FAIRD-AMMOD Correlation
Only **Site1, Site2, Site3** have overlapping FAIRD + AMMOD biomass data (Site4 excluded due to AMMOD4 failure)

## Working with Previous Session Work

This project has been developed collaboratively between **Gemini** and **Claude**:

- **Gemini** (Oct 28): Systematic implementation (log transform, tidyr::complete, H1/H2 scripts)
- **Claude** (Oct 28): Validation analysis (zero-days impact, documentation)

**When continuing work**:
1. Read `PROJECT_CONTEXT_2025_10_28_INTEGRATED.txt` for established decisions
2. Do NOT re-debate log transformation or zero-inclusion (finalized)
3. Maintain consistency with existing H1/H2 approaches for H3
4. Follow namespace standards strictly
5. Document all new decisions in project context files

## Package Dependencies

**Core Analysis**:
```r
library(tidyverse)   # Data manipulation, visualization
library(car)         # Diagnostic tests (ncvTest for Breusch-Pagan)
library(MASS)        # Robust regression (load LAST)
library(broom)       # Tidy model outputs
library(patchwork)   # Combine plots
library(ggpubr)      # stat_cor, stat_regline_equation
library(vegan)       # Diversity indices (H3)
```

**Post-Hoc Tests**:
```r
library(dunn.test)   # Dunn's test for Kruskal-Wallis
library(rstatix)     # epsilon_squared effect size
```

**Managed with renv**: Dependencies are tracked in `renv.lock`

## Next Priority Tasks

1. ⏳ **H3 Taxonomic Analysis**: Apply established standards (log if needed, data completeness automatic)
2. ⏳ **Integrate Section 6, 8**: Merge Device×Site and Environment analyses into main H2 script
3. ⏳ **Site3 Investigation**: Understand why Meadow validation is very weak (only n=7 days)
4. ⏳ **Manuscript Preparation**: Update Methods/Results with new log-transformed values

## Common Pitfalls to Avoid

1. **Never modify `tidyr::complete()` section in load_data.R** - foundational for all analyses
2. **Never use raw scale for H1/H2** - log transformation is finalized and validated
3. **Never use MASS functions before declaring intent** - load MASS last, use `::` notation
4. **Never hardcode R²/p-values in figures** - extract dynamically from models
5. **Never exclude zero-days without consultation** - established as real data
6. **Never use `ammod_biomass_for_stats` in H2 validation figures** - use `ammod_biomass_original` (includes AMMOD3 prorated) for visualization consistency

## Diagnostic Standards

All regressions must pass:
1. **Shapiro-Wilk test** (residuals): p > 0.05 (normality)
2. **Breusch-Pagan test** (ncvTest): p > 0.05 (homoscedasticity)
3. **Cook's Distance**: < 4/n threshold (influential points check)

If diagnostics fail on raw scale → apply log₁₀(x + 1) transformation
