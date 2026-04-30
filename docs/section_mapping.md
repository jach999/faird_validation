# Section Mapping: R Scripts → Draft V5

This document maps each section of the manuscript to its source scripts and outputs.
Used during refactoring to consolidate ~16,000 lines across 23 scripts into 4 core scripts.

**Source scripts** are the originals in the project root (e.g., `H1_consistency_FAIRD.R`).
**New curated scripts** are generated in `R/` (e.g., `R/H1_quanti_consistency.R`).
**Outputs** go to `outputs/tables/`, `outputs/figures/`, `outputs/console/`.

### Known issues to fix during refactoring
1. **`load_data.R` lines 539–541:** Remove `write_csv()` calls that export intermediate CSVs to root. All data available via `data_list`.
2. **`H2_validation_counts_vs_drymass.R` lines 54, 58:** Reads intermediate CSVs directly instead of using `data_list`. Fix in curated script.
3. **`analysis_functions.R`:** `analyze_device_consistency_robust()` is duplicated (lines 163–210 and 579–626). Remove duplicate. Add namespaces to `drop_na()`, `count()`, `column_to_rownames()`.

## Paper structure overview

```
§3    Results
├── §3.1    Diagnostic Mixed Model (SAS — not in R scripts)
├── §3.2    H1 — Intra-System Temporal Consistency
│   ├── §3.2.1   Quantitative Consistency  →  R/H1_quanti_consistency.R
│   └── §3.2.2   Qualitative Consistency   →  R/H1_quali_consistency.R
└── §3.3    H2 — Inter-System Temporal Agreement
    ├── §3.3.1   Quantitative Agreement    →  R/H2_quanti_validation.R
    └── §3.3.2   Qualitative Agreement     →  R/H2_quali_validation.R
```

---

## R/H1_quanti_consistency.R — §3.2.1

### Section 1: Diagnostic heteroscedasticity (§2.6, Figure 3)
- **Source:** `H1_consistency_FAIRD.R` lines ~1–120
- **Output:** Figure 3 (raw vs log diagnostic panels)
- **Note:** This figure is referenced in Methods (§2.6) but generated here

### Section 2: Intra-Habitat Level — FAIRD (§3.2.1.1)
- **Source:** `H1_consistency_FAIRD.R`
- **Analysis:** OLS regression FAIRD1↔FAIRD2 (Maize, n=22), FAIRD3↔FAIRD4 (Meadow, n=22)
- **Output:** Figures 4–5, Table 9 (FAIRD rows)
- **Metrics:** R², Spearman ρ, Kendall τ (on daily changes), concordance %, regression diagnostics

### Section 3: Intra-Habitat Level — AMMOD (§3.2.1.1)
- **Source:** `H1_consistency_AMMOD.R`
- **Analysis:** OLS regression AMMOD1↔AMMOD2 (Maize only, n=12). No Meadow pair (AMMOD4 failed).
- **Output:** Figure 6, Table 9 (AMMOD rows)

### Section 4: Supplementary — Robust regression (§3.2.1.2)
- **Source:** `H1_consistency_FAIRD.R` (robust regression section)
- **Analysis:** MASS::rlm() validation of OLS results
- **Output:** Table 9 supplement (OLS vs robust comparison)

### Section 5: Inter-Habitat Level — FAIRD (§3.2.1.3)
- **Source:** `H1_consistency_FAIRD_Maize_Meadow.R`
- **Analysis:** mean(FAIRD1+FAIRD2)/2 vs mean(FAIRD3+FAIRD4)/2, n=22 days
- **Output:** Figure 9, Table 10 (FAIRD row)
- **Note:** Uses AVERAGING approach (pooling invalid — site-pair covariate p=0.013)

### Section 6: Inter-Habitat Level — AMMOD (§3.2.1.3)
- **Source:** `H1_consistency_AMMOD_Maize_Meadow.R`
- **Analysis:** AMMOD1↔AMMOD3, AMMOD2↔AMMOD3 (no pairing possible; AMMOD4 failed)
- **Output:** Figures 10–11, Table 10 (AMMOD rows)

### Section 7: Summary tables (§3.2.1)
- **Source:** `H1_consistency_summary.R` (table consolidation logic only)
- **Output:** Tables 8, 9, 10 (final consolidated versions)

---

## R/H1_quali_consistency.R — §3.2.2

### Section 1: Intra-Habitat Level — FAIRD (§3.2.2.1)
- **Source:** `H1_consistency_FAIRD_taxonomic.R`
- **Analysis:** Morisita-Horn similarity in 3-day windows at Superfamily/Order/Genus levels
- **Output:** Figures 12–13, Table 11 (FAIRD rows)
- **Index:** Morisita-Horn (abundance-based, appropriate for FAIRD count data)

### Section 2: Intra-Habitat Level — AMMOD (§3.2.2.1)
- **Source:** `H1_consistency_AMMOD_taxonomic.R`
- **Analysis:** Izsak-Price similarity daily + 3-day windowed at Species/Genus/Superfamily
- **Output:** Figures 14–15, Tables 11 (AMMOD rows), 12
- **Index:** Izsak-Price (presence/absence-based, required for AMMOD OTU data)

### Section 3: Inter-Habitat Level — FAIRD (§3.2.2.2)
- **Source:** `H1_H2_H3_taxonomic_size_analysis.R` (H1 FAIRD portions)
- **Analysis:** Cross-habitat Morisita-Horn: Pool(FAIRD1+2) vs Pool(FAIRD3+4)
- **Output:** Table 11 (inter-habitat rows)

### Section 4: Inter-Habitat Level — AMMOD (§3.2.2.2)
- **Source:** `H1_H2_H3_taxonomic_size_analysis.R` (H1 AMMOD portions)
- **Analysis:** Cross-habitat Izsak-Price
- **Output:** Table 11 (inter-habitat rows)

---

## R/H2_quanti_validation.R — §3.3.1

### Section 1: Site Level (§3.3.1.1)
- **Source:** `H2_validation_FAIRD_AMMOD.R`
- **Analysis:** Per-site OLS: Site1 (n=12), Site2 (n=12), Site3 (n=8). Site4 excluded.
- **Output:** Figures 16–18, Table 13 (site-level rows)
- **Metrics:** R², slope [95% CI], Spearman ρ, Kendall τ (daily changes), concordance %
- **Key results:** Site1 R²=0.282 ns, Site2 R²=0.445*, Site3 R²=0.058 ns

### Section 2: Intra-Habitat Level (§3.3.1.2)
- **Source:** `H2_validation_aggregated_levels.R` (Maize section)
- **Analysis:** Maize pooled (Sites 1+2, n=24). Site covariate p=0.861 → pooling valid.
- **Output:** Figure 19, Table 13 (Maize Habitat row)
- **Key result:** R²=0.311**

### Section 3: Overall Level (§3.3.1.3)
- **Source:** `H2_validation_aggregated_levels.R` (Overall section)
- **Analysis:** All sites pooled (n=32). Site covariate all p>0.64 → pooling valid. Bootstrap CIs.
- **Output:** Figure 20, Table 13 (Overall row)
- **Key result:** R²=0.321*** (primary result per Scherber request)

### Section 4: Cross-correlation analysis (§3.3.1.4)
- **Source:** `H2_validation_cross_corr.R` (1581 lines — longest single script)
- **Analysis:** CCF at lag-0 per site + Intra-Habitat summed series
- **Output:** Figures 21–24, Table 13 (CCF column)
- **Note:** Overall Level CCF excluded from manuscript (n=8 at Site3 insufficient). Reported in table with methodological note.

### Section 5: Detection metric sensitivity (§3.3.1.5)
- **Source:** `H2_validation_counts_vs_drymass.R`
- **Analysis:** R² comparison: raw counts vs dry mass vs allometric biomass inference
- **Output:** Table 14
- **Key finding:** Allometric biomass > dry mass > raw counts (ΔR²=+0.01 to +0.11)
- **⚠️ FIX:** Original reads `etraps_daily_biomass_clean.csv` and `ammod_biomass_for_stats.csv` directly. Replace with `data_list$etraps_daily` and `data_list$ammod_biomass_for_stats`.

### Section 6: Body size distributions (§3.3.1.6)
- **Source:** `H2_H3_extra_figures.R`
- **Analysis:** Body length distributions FAIRD vs AMMOD across sites
- **Output:** Table 15, supplementary figures

---

## R/H2_quali_validation.R — §3.3.2

### Section 1: Site Level (§3.3.2.1)
- **Source:** `H3_taxonomic_FAIRD_AMMOD_superfamily.R`
- **Analysis:** Izsak-Price ΔS at Superfamily level, 3-day windows, per site
- **Output:** Figures 25–27, Table 16 (site rows)

### Section 2: Intra-Habitat Level (§3.3.2.2)
- **Source:** `H3_taxonomic_FAIRD_AMMOD_Maize_Meadow.R`
- **Analysis:** t-test, Cohen's d, Levene's for Maize vs Meadow I-P distributions
- **Output:** Figure 26 (habitat comparison), Table 16 (habitat rows)

### Section 3: Overall Level (§3.3.2.3)
- **Source:** `H3_taxonomic_FAIRD_AMMOD_superfamily.R`
- **Analysis:** Pooled I-P across all sites
- **Output:** Table 16 (Overall row)
- **Key result:** Overall ΔS = 67.8%

### Section 4: Size-threshold sensitivity (§3.3.2.4)
- **Source:** `H2_validation_size_threshold.R` (→ renamed H2_quali_validation_size_threshold.R)
- **Analysis:** Three-level sensitivity: raw → ≥3mm → ≥5mm. Sub-3mm diagnostic.
- **Output:** Table 17
- **Key findings:** Near-zero raw→3mm delta (+0.2 pp); 10 exclusively sub-3mm Superfamilies present in 26/26 windows (mean 1.9/window)

### Section 5: Taxonomic breadth (§3.3.2.5)
- **Source:** `H3_taxonomic.R` (relevant portions only — script contains obsolete analyses)
- **Analysis:** Diptera dominance (91.8% Maize, 84.5% Meadow). Taxonomic resolution by Order. Classification success rates.
- **Output:** Tables 18–19

---

## SAS analyses — §3.1

- **Preparation:** `SAS/FAIRD_ID_AMMOD_prepare_sas_data.R` (generates filtered CSVs from data_cleaned.RData)
- **Execution:** `SAS/FAIRD_ID_AMMOD.sas` (PROC MIXED, REML, Kenward-Roger df, AR(1))
- **Analysis A:** Biomass validation (Table 6) — device_type × site, days 1–12, Sites 1–3
- **Analysis B:** Abundance comparison FAIRD vs ID (Table 7) — device_type × site, days 1–22, Sites 1,2,4
- **Results transcribed manually to manuscript**

---

## Scripts NOT in paper (retained for documentation)

| Script | Location | Purpose |
|--------|----------|---------|
| `H1_consistency_ID.R` | `internal/` | ID intra-system consistency (excluded from paper) |
| `H2_quanti_validation_ID_AMMOD.R` | `internal/` | ID vs AMMOD validation (excluded from paper) |
| `H2_divergence_diagnostic.R` | `diagnostics/` | AMMOD-exclusive taxa on high-residual days (feeds §4.2 Discussion) |
