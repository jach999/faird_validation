# CLAUDE.md — FAIRD Validation Project

## Project

R-based statistical analysis validating the FAIR-Device (FAIRD) against AMMOD (DNA metabarcoding + gravimetric biomass) for automated insect monitoring. Target journal: *Methods in Ecology and Evolution*. Manuscript: Draft V5.

Study: 22 days, 2 habitats (Maize/Meadow), 4 sites, 12 devices total (Triesdorf, Bavaria, Aug–Sep 2023).

## Repository structure

### Target (post-refactoring)
```
FAIRD_validation/
├── CLAUDE.md
├── FAIRD_ID_AMMOD.Rproj          # Root-level .Rproj (CWD = FAIRD_validation/) — NOT used for curated scripts
├── docs/
│   ├── section_mapping.md        # Script → paper section mapping (READ FIRST for refactoring)
│   └── namespace_standards.md    # R namespace rules (MANDATORY)
├── R/                            # CURATED scripts — open R/FAIRD_ID_AMMOD.Rproj to run
│   ├── FAIRD_ID_AMMOD.Rproj      # ← USE THIS .Rproj (CWD = R/)
│   ├── 00_load_data.R            # Data loading + Truth Calendar + completeness enforcement
│   ├── 00_custom_colors.R        # Shared color palette
│   ├── 00_analysis_functions.R   # Shared functions (cleaned)
│   ├── H1_quanti_consistency.R   # §3.2.1 — Quantitative intra-system consistency
│   ├── H1_quali_consistency.R    # §3.2.2 — Qualitative intra-system consistency
│   ├── H2_quanti_validation.R    # §3.3.1 — Quantitative inter-system agreement
│   ├── H2_quali_validation.R     # §3.3.2 — Qualitative inter-system agreement
│   └── outputs/                  # All outputs from curated scripts
│       ├── tables/               # CSV summary tables (Tables 8–19)
│       ├── figures/              # PNG/PDF figures (Figs 3–27)
│       └── console/              # Console capture *_output_*.txt
├── source_tables/                # 4 CSV input files (read-only)
├── SAS/                          # LMM analyses (PROC MIXED) — §3.1
├── diagnostics/                  # H2_divergence_diagnostic.R
├── internal/                     # ID-related scripts (excluded from paper)
└── old/                          # Archived/superseded scripts
```

### During transition
Original scripts remain in the project root untouched. New curated scripts are generated
in `R/` with paths relative to the project root. Once a new script is validated, the
original scripts it replaces are moved to `old/`. Do NOT modify original scripts.

## Running analyses

**Open `R/FAIRD_ID_AMMOD.Rproj`** (the one inside `R/`, not the root one). This sets the working directory to `R/`, so all plain relative paths resolve correctly within that directory.

Each analysis script must begin with:
```r
library(here)
here::i_am("R/<script_name>.R")   # anchors here root to FAIRD_validation/
```

`here::i_am()` walks upward from the script until it finds the directory that contains `R/<script_name>.R` as a relative path. It sets **that** directory (`FAIRD_validation/`) as the `here` root, without changing the working directory. This is needed so `00_load_data.R` can resolve data file paths via `here::here()` regardless of how the script is invoked.

```r
# Inside R/H1_quanti_consistency.R:
library(here)
here::i_am("R/H1_quanti_consistency.R")

source("00_load_data.R")      # plain name — works because CWD = R/
source("00_custom_colors.R")  # plain name — works because CWD = R/

# Output paths: relative to CWD (R/) — land in R/outputs/
sink("outputs/console/H1_quanti_consistency_output.txt", split = TRUE)
ggsave("outputs/figures/Fig04_H1_FAIRD_Maize.png", width = 10, height = 6)
readr::write_csv(table_9, "outputs/tables/Table09_H1_intra_habitat.csv")
```

Inside `00_load_data.R`, source data files use `here::here()` so they resolve to the project root regardless of CWD:
```r
here::i_am("R/00_load_data.R")
readr::read_csv(here::here("source_tables", "FAIRD_For_Statistics.csv"))
```

Access datasets via `data_list$faird`, `data_list$etraps_daily`, `data_list$ammod_biomass_for_stats`, etc.

**Never `setwd()` inside scripts. Never use the root-level `.Rproj`for curated scripts.**

## Refactoring corrections (apply when creating curated scripts)

1. **`R/00_load_data.R`:** Remove the three `write_csv()` calls (lines 539–541 in original) that export intermediate CSVs (`etraps_daily_biomass_clean.csv`, `ammod_biomass_for_stats.csv`, `ammod_taxonomy_presence_clean.csv`) to the project root. All data is already available via `data_list` + `data_cleaned.RData` — intermediate CSVs are unnecessary and clutter the workspace.
2. **All analysis scripts:** Never read intermediate CSVs. Always use `data_list$` objects loaded via `source("00_load_data.R")`. The original `H2_validation_counts_vs_drymass.R` reads CSVs directly — fix this in `R/H2_quanti_validation.R`.
3. **Output files only in `outputs/`:** Scripts must never write files to the project root. All outputs go to `outputs/tables/`, `outputs/figures/`, or `outputs/console/`.

## Critical data constraints

- **AMMOD taxonomy is presence/absence ONLY.** The Count column = DNA reads, NOT individual counts. Never use it for abundance.
- **AMMOD biomass and taxonomy are NOT linked.** Biomass = gravimetric; taxonomy = metabarcoding. Cannot connect them per-taxon.
- **AMMOD4:** No biomass data (mechanical failure). Site4 excluded from H2.
- **AMMOD3:** 8 of 12 biomass days valid (processing method changed in cursu). OTU data available for all 21 days.
- **ID3:** Complete device failure. ID excluded from manuscript entirely.
- **Zero-capture days are real observations**, not missing data. Enforced via `tidyr::complete()` in `00_load_data.R`. Never remove this section.

## Analytical hierarchy

- **H1 (intra-system):** Intra-Habitat Level → Cross-Habitat Level (daily averaging only; pooling invalid, site-pair covariate p=0.013)
- **H2 (inter-system):** Site Level → Intra-Habitat Level → Overall Level. Cross-Habitat Level excluded (analytically invalid for inter-system comparison).
- **Overall Level** is the primary level for generalizability (Scherber request). R²=0.321***.

## Established methodological decisions (DO NOT re-debate)

1. **Log transformation:** `log10(x + 1)` for ALL biomass/abundance regressions. Validated via Shapiro-Wilk + Breusch-Pagan.
2. **Zero-days inclusion:** Always. Centralized in `00_load_data.R` via `tidyr::complete()`.
3. **Biomass inference > raw counts:** Allometric approach (Söhlström et al. 2018) produces stronger correlations (ΔR²=+0.01 to +0.11).
4. **3-day windows** for taxonomic concordance indices (Morisita-Horn for FAIRD, Izsak-Price for FAIRD↔AMMOD).
5. **Superfamily** as operative taxonomic level for H2 qualitative (§3.3.2).
6. **CCF aggregation:** Sum daily biomass across Maize sites before computing single CCF (treats Maize as habitat unit).

## R coding standards (MANDATORY)

**Always use explicit namespaces** for tidyverse functions. See `docs/namespace_standards.md` for complete list.

```r
# CORRECT
dplyr::filter()  dplyr::select()  dplyr::mutate()  dplyr::summarise()
dplyr::group_by()  dplyr::arrange()  dplyr::rename()  dplyr::pull()
dplyr::left_join()  dplyr::inner_join()  dplyr::bind_rows()
dplyr::case_when()  dplyr::if_else()  dplyr::n()  dplyr::lag()
tidyr::pivot_wider()  tidyr::pivot_longer()  tidyr::drop_na()
tidyr::complete()  readr::read_csv()  readr::write_csv()

# WRONG — will fail when MASS is loaded
filter()  select()  mutate()  summarise()
```

**Additional conventions:**
- Comments, variable names, and print outputs in **English**
- `source("R/00_load_data.R")` and `source("R/00_custom_colors.R")` at script start
- `while (sink.number() > 0) sink()` before every `sink()` call
- `custom_colors` passed as full named vector to `scale_color_manual()` (subsetting breaks ggplot2)
- No hardcoded statistical values in figures — extract dynamically from models
- Console output capture via `sink(output_file, split = TRUE)`

## Diagnostic standards

All regressions must report:
1. **Shapiro-Wilk** (residual normality): p > 0.05
2. **Breusch-Pagan** via `car::ncvTest()` (homoscedasticity): p > 0.05
3. **Durbin-Watson** via `lmtest::dwtest()` (autocorrelation)
4. **Non-parametric:** Spearman ρ + Kendall τ (on daily changes, not absolute values)

## Package dependencies

```r
# Core
library(tidyverse)    # dplyr, tidyr, ggplot2, readr, stringr
library(patchwork)    # Multi-panel figures
library(scales)       # Axis formatting
# Diagnostics
library(car)          # ncvTest (Breusch-Pagan)
library(lmtest)       # dwtest (Durbin-Watson)
library(MASS)         # rlm (robust regression) — LOAD LAST
# Taxonomic
library(vegan)        # Diversity indices (Morisita-Horn)
# Visualization
library(ggrepel)      # Non-overlapping labels
```

## Terminology

- "System" for FAIRD/AMMOD/ID (not "method" in device contexts)
- Maize/Meadow capitalized as experimental categories
- "Inter-" reserved for CCF statistical term; "Inter-System/Inter-Habitat" for device comparisons
- Superfamily as operative taxonomic level for H2 qualitative
- "Intra-Habitat Level", "Overall Level" — not "Habitat Level" or "Global Level"

## Key collaborators

- **Christoph Scherber:** Thesis supervisor. Requested generalization beyond individual devices.
- **Max Sittinger:** Insect Detect developer. Language about ID must remain diplomatic/non-critical.