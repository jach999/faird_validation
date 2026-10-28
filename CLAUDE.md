# CLAUDE.md, FAIRD validation project

Replaces the V5-era CLAUDE.md (archived as `docs/old/CLAUDE_V5_refactoring.md` on
2026-10-27). Its refactoring instructions, the exclusion of Insect Detect from the
manuscript and the fixed V5 analysis rules no longer apply.

## Project

Field validation of the FAIR-Device (FAIRD), a non-lethal camera trap, against the AMMOD
multisampler (Malaise trap, gravimetric biomass and DNA metabarcoding) and Insect Detect
(ID). Triesdorf, Bavaria, 23 August to 13 September 2023; two habitats (Maize, Meadow),
four sites, one unit per system and site. Target journal: Methods in Ecology and
Evolution. Current draft: V5 (in the claude.ai project); V6 is being planned from
`docs/V6_OUTLINE.md`.

## Where things are

- Start-of-session routine, role of each document and closing routine:
  `claude/NEW_CONVERSATION_GUIDE.md` in the claude.ai project (not in this repository).
- Results with status and sources: `docs/RESULTS_REGISTRY.md`. Results come from here,
  not from V5 prose or from memory.
- V6 structure and claims: `docs/V6_OUTLINE.md`.
- Scherber's guidelines: `claude/20261013_SCHERBER_GUIDELINES_EN.md` (project).
- Design, incidents and methods: `EXPERIMENT_CONTEXT_STATIC_V1_5_EN.md` (project).
- Branches: work on `revision-scherber`; `main` holds Scherber's commits; tag
  `v5-supplementary` is the frozen V5 state.

## Repository layout

```
R/            one script per V6 section (since 2026-10-28):
              00_load_data.R, 00_custom_colors.R, 00_trajectory_functions.R (shared)
              2-4_data_availability.R        Methods 2.4, Table 1
              3-1a_quantity_trajectories.R   Results 3.1, Figure 3, trajectory tests
              3-1b_quantity_concordance.R    Results 3.1, concordance, simple correlations
              3-2a_composition_orders.R      Results 3.2, Figure 4, Table 2
              3-2b_composition_overlap.R     Results 3.2, Izsak-Price, Table 3
R/supplementary/  S1 to S4, one script per supplementary section
R/old/        archived scripts, not part of the published repository
              (stage2_v5/: V5 framework; stage3_versions/: superseded Stage 3 versions)
outputs/      console/ and tables/, file names prefixed with the script name;
              figures not versioned; old/ subfolders (and outputs/old/) hold outputs
              of archived scripts, not part of the published repository
source_tables/  the four input CSV files (read-only)
SAS/          V5 diagnostic mixed models (program and results)
diagnostics/  throwaway diagnostic scripts
docs/         registry, outline, namespace standards, old/ for retired documents
```

## Running scripts

- Juan runs the scripts and shares the outputs; Claude does not run them.
- Working directory = repository root (root `FAIRD_validation.Rproj`).
- Scripts source `R/00_load_data.R` and `R/00_custom_colors.R`; `00_load_data.R` anchors
  paths with `here::i_am("R/00_load_data.R")`. Paths via `here::here()`, never
  `setwd()` or hardcoded paths.
- Console output to `outputs/console/`, tables to `outputs/tables/`, each prefixed with
  the script name (`<script>_output.txt`, `<script>_<table>.csv`).
- File names follow the V6 sections and contain no person names; the same holds for
  outputs and printed text. Only `R/old/` and the `old/` output folders are exempt.

## R conventions

- Explicit namespaces for tidyverse functions (`dplyr::`, `tidyr::`, ...); MASS masks
  dplyr. Full list: `docs/namespace_standards.md`.
- Comments, variable names and printed output in English; no apostrophes in comments
  (Windows encoding).
- `while (sink.number() > 0) sink()` before the main `sink()`; main sink with
  `split = TRUE`, nested diagnostic sinks with `split = FALSE`.
- No hardcoded statistical values; scale factors and labels derived from the data.
- Scripts meant for publication contain no editorial conclusions.

## Data facts

- AMMOD taxonomy is presence/absence only: the `Count` column holds DNA reads, not
  individuals. AMMOD biomass (gravimetric) and taxonomy (metabarcoding) are not linked
  per taxon.
- AMMOD biomass: `Live_mass`, the unit of the camera systems' allometric estimates.
  AMMOD1 and AMMOD2 12 days, AMMOD3 8 days, AMMOD4 none.
- ID3 recorded nothing (SD card). ID power outages from 28 August to 1 September give
  partial, unknown-effort days; analyses use the full series and the window excluded,
  the exclusion symmetric for both systems of a comparison.
- Zero-capture days are real observations, completed in `00_load_data.R`.

## Analysis and figures

- Current models: glmmTMB Tweedie on the raw response (log link), Date spline by
  Device_type, Device as random effect. V5 methods (OLS on log10(x+1), pairwise
  correlations, SAS LMM) are legacy.
- Main-text comparisons are between systems; comparisons between individual units go to
  supplementary material.
- Figures on linear axes, never logarithmic; AMMOD on a secondary axis divided by a
  data-derived factor.

## Terminology and text

- "System" for FAIRD, AMMOD and ID; "unit" for one device; "method" only in statistical
  contexts. Maize and Meadow capitalised. Superfamily is the operative taxonomic level.
- No em dashes.
- Language about Insect Detect stays descriptive and non-competitive (its developer, Max
  Sittinger, is a close colleague).
