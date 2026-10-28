# FAIR-DEVICE VALIDATION STUDY, RESULTS REGISTRY

**Version:** 0.7.4 (October 28, 2026). R-TB-INT-01: Maize-pair LRT values corrected
(transcription error; reading unchanged), checked against the rerun of
`3-1a_quantity_trajectories.R`. The reorganised scripts 2-4, 3-1a, 3-1b, 3-2a and 3-2b
reproduce every value in this registry.
v0.7.3 (October 28, 2026). Conventions: which script files produced the
current results, and superseded versions moved to `R/old/stage3_versions/`; adds
R-CC-DESC-01 (Section 2.2, simple correlations of daily means, biomass and abundance);
R-DC-BM-02: factor k printed (11); R-DC-TX-03: denominators; both Section 6 rows
removed. No earlier result changed.
v0.7.2 (October 27, 2026). File moved from ADP (`claude/20261012_RESULTS_REGISTRY.md`)
to the repository as `docs/RESULTS_REGISTRY.md`; update rule, location of Scherber's
original script and reference to his criticisms (now `claude/20261013_SCHERBER_GUIDELINES_EN.md`)
updated in Conventions. No result changed.
v0.7.1 (October 13, 2026): diagnostic mixed models A and B (Section 2.5)
verified against the SAS results file in the repository; their row removed from Section
6; path of PC V16 updated (`docs/old/`). Otherwise as v0.7.
v0.7 (October 12, 2026, second update): populates Section 2.3 (temporal
abundance, FAIRD vs ID: R-TA-FI-01, R-TA-FI-02 planned) and adds R-DX-FI-ABUND-01, from
the run of 2026-09-02 and its session summary (summary retired from ADP); adds two rows
to Section 6; updates file locations after the move of scripts and outputs to the GitHub
repository (Conventions, common notes of Section 2.2); everything else as in v0.6.
v0.6 (October 12, 2026): adds Section 2.6 (Scherber's device-comparison
script of 2026-10-08: multinomial order composition and Tweedie biomass model, original
and input-corrected runs), branch `DC`, observations R-OBS-05 to -07, a manuscript-level
flag on device-pair references (Scherber criticism C2), the 2026-10-12 rerun note of the
concordance script, and new rows in Section 6; everything else as in v0.5.
v0.5 (October 8, 2026): adds the results of the registry diagnostics
(Section 8 of `H2_quanti_validation.R`, run of 2026-10-08) to R-V5-H2Q-01 to -04,
R-OBS-02 and Sections 5 and 6; everything else as in v0.4.
v0.4 (October 7, 2026): review of every LEGACY entry against Draft V5,
PROJECT_CONTEXT V16 and the 2026-07-22 script outputs completed: Sections 1.1, 1.2, 2.2,
2.4 and 2.5 (new), OUTPUT ONLY entries from the H1 and H2 outputs, PC V16 tests dropped
in V5 (entered as SUPERSEDED). New sections: 4 (observations), 5 (manuscript
corrections for V6), 6 (values pending resolution). Stage 3 entries (R-TB, R-CC, R-DX)
unchanged from v0.3.
**Scope:** current and superseded results of every analysis, one entry per analysis.
Methods live in EXPERIMENT_CONTEXT_STATIC (V1.5); the chronology of decisions lives in
the session logs. This file answers one question: what is the result of analysis X, and
is it still valid.

## Conventions

- **ID scheme:** `R-<branch>-<pair or system>-<nn>`. Branches: `TB` temporal biomass
  model, `TA` temporal abundance model, `CC` system-level concordance, `TX` taxonomic,
  `DX` diagnostic or data fact, `OBS` observation (Section 4), `DC` Scherber's
  device-comparison script (Section 2.6; sub-tags `TX` taxon composition, `BM`
  biomass), `V5` legacy result from Draft V5, `V16` legacy result
  from PROJECT_CONTEXT V16 (Stage 1) not carried into V5, `OUT` result present only in
  a script output. Pairs: `AF` AMMOD vs
  FAIRD, `AI` AMMOD vs ID, `FI` FAIRD vs ID, `ALL` three systems. Sub-tags after `V5`,
  `V16` and `OUT`: `H1Q`, `H1T`, `H2Q`, `H2T`, `H3T` (Stage 1 or 2 hypothesis and
  dimension), `MM` diagnostic mixed model, `LV` level comparison between systems.
- **Status:** `CURRENT`, `SUPERSEDED by <ID>`, `LEGACY V5` (Stage 2 device-to-device
  framework, kept for reference, not for new claims), `LEGACY PC V16` (Stage 1 result
  absent from V5, same restrictions), `OUTPUT ONLY` (result printed by a Stage 2 script
  output of 2026-07-22 but reported neither in V5 nor in PC V16; absence checked by
  searching both documents for the values; same restrictions as LEGACY), `PLANNED`.
  `PROVISIONAL` (v0.7, added to CURRENT): result of a run whose design is known to
  deviate from a decision taken later; kept until the rerun, not for the manuscript.
- **Flags** (added to the status):
  - `DISCREPANCY V5/V16`: V5 and PC V16 report different values for the same test. Both
    values are recorded, V5 is kept, and the entry states `Cause: verified` or
    `Cause: hypothesis` with its evidence.
  - `V5 NOT REPRODUCED`: the V5 value differs from the output of the script that
    produced it (outputs of 2026-07-22). The V5 value is kept, the output value is
    recorded beside it, and the item is listed under "Manuscript corrections for V6".
  - `NOT FOR MAIN TEXT (C2)` (v0.6): the result rests on a comparison between
    individual devices (device pair as unit). Scherber's criticism C2
    (`claude/SCHERBER_FEEDBACK_V5_EN.md`): device-device comparisons are not
    scientifically admissible; the main argument uses system-level comparisons over
    both blocks. Kept as diagnostic or appendix material only (R-OBS-05).
  - Criticism IDs C1 to C9 (v0.7.2 note): `SCHERBER_FEEDBACK_V5_EN.md` was replaced on
    2026-10-13 by `claude/20261013_SCHERBER_GUIDELINES_EN.md`, which numbers his
    guidelines G1 to G18 and maps the old IDs in its Section 4 (C1 = G6, C2 = G8 and
    G11, C4 = G14, C9 = G3). The flag keeps its name.
- **Sources:** V5 = `FAIRD_Paper_V5.docx`, cited by Results section number and table;
  PC V16 = `docs/old/20260317_PROJECT_CONTEXT_V16_EN_UPDATE_NEEDED.md` (repository), cited by
  block heading (its table of contents line numbers are not reliable); outputs =
  `<script>_output.txt` with line numbers, now in `outputs/console/` of the repository
  (except `H1_consistency_summary_output_*.txt`, in `outputs/`). Line numbers refer to
  the files as of 2026-10-06, except where an entry or common note states otherwise.
- **File locations (v0.7):** scripts and outputs live in the GitHub repository
  `jach999/faird_validation`, branch `revision-scherber`, synced to ADP. Scripts cited by
  bare name are in `R/` (current) or `R/old/` (archived). Earlier entries that say a
  file is "in ADP" or "under `claude/`" refer to the situation at the time of the entry;
  the file is now in the repository unless stated otherwise.
- **Superseded entries are never deleted.** They are marked and linked to their
  replacement, so that discarded analyses are not rerun.
- **Reading line:** one line, statistical reading only. Ecological interpretation belongs
  to the manuscript.
- **Update rule:** at the end of each session, from its session log. ONE cumulative
  document: entries are added or marked in place, never moved to a separate file. Since
  v0.7.2 the file lives in the repository as `docs/RESULTS_REGISTRY.md`, without a date
  in its name: the version and date are in the header, and git keeps the history.
  Until v0.7.1 it lived in ADP and was renamed with the date of each update.
- **Scripts:** Section 2e is unchanged from `Scherber_temporal_biomass_3.R` to
  `Scherber_temporal_biomass_5.R`; version 5 adds Sections 2f (identity link) and 2g
  (internal consistency). Concordance: `Scherber_temporal_concordance_biomass.R`.
  Section 2.6: Scherber's `Christoph_device_comparisons_2026-10-08_taxon_and_biomass_analyses.R`
  (with `Christoph_Device_Comparisons.RData`), run through two versions prepared on
  2026-10-12: `Christoph_device_comparisons_2026-10-08_runnable.R` (runs with
  `source()`, analyses unchanged) and `Christoph_device_comparisons_2026-10-08_corrected.R`
  (same models, corrected inputs). Original script in `R/` (the file is named
  `Christoph_device_comparisons_2026-10-08_taxon and biomass analyses.R`, with spaces;
  byte-identical to the copy on `main`, brought to `revision-scherber` on 2026-10-13,
  the ADP copy deleted); runnable and corrected versions in `R/`; run outputs in
  `outputs/console/`.
- **Script files (v0.7.3).** Checked on 2026-10-28 against the outputs in the
  repository. Entries citing `Scherber_temporal_biomass_5.R` for Sections 2e to 2g refer
  to the file `R/Scherber_temporal_biomass_6.R`, whose header also says "Version 5" and
  which is the only version with Section 2g (it writes
  `Scherber_biomass_internal_consistency.csv`); the file named `_5.R` has Sections 2e and
  2f only. Entries citing `Scherber_temporal_concordance_biomass.R` refer to the current
  `R/Scherber_temporal_concordance_biomass_4.R`, which wrote
  `outputs/tables/Scherber_concordance_biomass.csv`; the console output of the
  2026-10-12 rerun was written by `Scherber_temporal_concordance_biomass (1).R`, the
  same analysis plus the internal parallelism LRT, with identical values (checked by
  R-CC-DESC-01 against the table). The file without suffix is the earlier version with
  the bootstrap p-value that never finished. Superseded versions (`Scherber_temporal_biomass.R`
  and `_2` to `_5`, concordance without suffix, `(1)` and `_2`, `H2_bridge_exploratory.R`,
  `Scherber_temporal_biomass_block_comparison.R`, `diagnostic_per_site_biomass.R`) moved to
  `R/old/stage3_versions/`; `diagnostic_per_site_biomass_2.R` contains Table 3c cited by
  R-DX-INTRA-01. File names are kept until the repository reorganisation.
- **Repository reorganisation (2026-10-28, after v0.7.3).** Scripts renamed after the V6
  sections; entries keep their old citations until the new scripts have been run and
  checked against the old outputs, then the citations are replaced in one pass. Mapping:

  | Old script (now in `R/old/`) | New script |
  | --- | --- |
  | `Scherber_temporal_biomass_6.R`, Sections 1, 2a to 2e (primary and complementary pairs), 2g | `3-1a_quantity_trajectories.R` (functions in `00_trajectory_functions.R`) |
  | `Scherber_temporal_biomass_6.R`, Block random, nlme cross-check, pairing sensitivities, post hoc, 2f | `supplementary/S1_pair_model_sensitivity.R` |
  | `Christoph_device_comparisons_2026-10-08_corrected.R`, biomass part (model, k, dual-axis figure) | `3-1a_quantity_trajectories.R` |
  | `Christoph_device_comparisons_2026-10-08_corrected.R`, taxon part (multinomial, Figure 1 of the original) | `3-2a_composition_orders.R` |
  | `Scherber_temporal_concordance_biomass_4.R` and `descriptive_correlations_daily_means.R` | `3-1b_quantity_concordance.R` |
  | `H2_quali_validation.R` | `3-2b_composition_overlap.R` (renamed, analysis unchanged) |
  | `diagnostic_per_site_biomass_2.R`, `diagnostic_per_device_operational.R` | `supplementary/S3_site_diagnostics.R` (abundance part now keeps all days with the outage window flagged for both systems) |
  | `Scherber_temporal_abundance_FAIRD_ID.R` | `supplementary/S4_abundance_cameras.R` (descriptive only; the provisional tests of R-TA-FI-01 stay in `R/old/`) |
  | `H1_quanti_consistency.R`, `H2_quanti_validation.R` (pair correlations only) | `supplementary/S2_unit_pair_correlations.R` |
  | (new) | `2-4_data_availability.R` (Table 1, valid unit-days, checks of the validity calendar) |
- **Extraction marks:** the `[cite: ...]` marks of the v0.2 LEGACY rows were removed in
  v0.4; every value carries its exact source.
- **Observations:** Section 4 holds cross-entry syntheses (`R-OBS-nn`) with the fields
  Rests on, Statement, Basis (verified or hypothesis) and Implication. An observation is
  reviewed whenever an entry it rests on changes status. The Reading line of each entry
  stays statistical and entry-specific.
- **Corrections and pending values:** every V5 NOT REPRODUCED item is listed in Section
  5. Values whose cause is not determined, that no output in ADP prints, or whose
  hypothesis can only be checked with a source not in ADP are listed in Section 6 with
  what would resolve them. An item can appear in both.
- **Stage 1 scripts** named in PC V16: `H1_consistency_ID.R`,
  `H2_validation_aggregated_levels.R` and `H2_validation_cross_corr.R` are in `R/old/`
  of the repository (since 2026-10-10 in ADP, since 2026-10-12 in the repository);
  `BodySize_Distribution.R` and `TaxConsistency_IntraCrossHabitat.R` are not available.
- p-values are LRT unless stated; "Holm" = Holm-adjusted over the primary family.

---

## 1. INTERNAL CONSISTENCY (within each system)

### 1.1 Quantitative (biomass, abundance)

| ID | Status | Note |
| --- | --- | --- |
| R-V5-H1Q-FAIRD-01 | LEGACY V5 | FAIRD replicate pairs within habitat, daily abundance (entry below) |
| R-V5-H1Q-FAIRD-02 | LEGACY V5, V5 NOT REPRODUCED | FAIRD Maize vs Meadow, daily abundance (entry below; was part of FAIRD-01 in v0.3) |
| R-V5-H1Q-AMMOD-01 | LEGACY V5, DISCREPANCY V5/V16 | AMMOD1 vs AMMOD2, daily biomass (entry below) |
| R-V5-H1Q-AMMOD-02 | LEGACY V5, V5 NOT REPRODUCED | AMMOD Maize devices vs AMMOD3, daily biomass (entry below) |
| R-V5-H1Q-ID-01 | LEGACY V5, DISCREPANCY V5/V16, V5 NOT REPRODUCED | ID1 vs ID2, conservative and operational, daily abundance (entry below) |
| R-V5-H1Q-ROB-01 | LEGACY V5 | Robust regression check, log scale (entry below) |
| R-V16-H1Q-ROB-01 | LEGACY PC V16 | Robust regression check, original scale (entry below) |
| R-V16-H1Q-DIAG-01 | LEGACY PC V16 | Raw vs log10(x+1) diagnostic behind the transformation decision (entry below) |
| R-OUT-H1Q-01 | OUTPUT ONLY | Spearman rho on daily log changes, all H1 pairs (entry below) |
| R-OUT-H1Q-02 | OUTPUT ONLY | OLS residual diagnostics of the H1 pairs (entry below) |
| R-CC-INT-* | DONE, see R-TB-INT-01 and R-CC-INT-01 below | System-level concordance among devices of the same system |

Common to the R-V5-H1Q entries (Stage 2 framework, `H1_quanti_consistency.R`):
- Method: pairwise OLS on log-transformed daily values; Pearson R² (primary), Spearman
  rho, temporal concordance (% of consecutive day-pairs moving in the same direction),
  Kendall tau on daily changes. Source: V5, Methods, "H1, Quantitative Consistency".
- Metric: daily abundance (individual counts) for FAIRD and ID; daily live-mass biomass
  (mg) for AMMOD (AMMOD has no abundance). Source: same V5 section. This confirms the
  v0.3 metric note.
- Transformation: V5 states only "log-transformed"; log10(x+1) per PC V16, H1
  Methodology, and EXPERIMENT_CONTEXT_STATIC V1.5, Section 14, Stage 1.
- Zero-capture days retained.
- Comparison with PC V16 made against its H1 summary table and H1 result blocks
  (PC V16, "H1: TEMPORAL CONSISTENCY").
- Verification against script output: `H1_quanti_consistency_output.txt` (run of
  2026-07-22). Spearman rho is on absolute log values and Kendall tau on raw daily
  changes (output, l.331). The script has no ID section, so ID entries cannot be
  verified against it.

**R-V5-H1Q-FAIRD-01, FAIRD replicate pairs within habitat** | LEGACY V5
- Level: Intra-Habitat (device pair within one habitat). Window: full 22-day series.
- Source: V5, Results 3.2.1.1 (Intra-Habitat Level) and Table 8.
- Result: FAIRD1 vs FAIRD2 (Maize, n = 22): R² = 0.674, p < 0.001; slope 0.84
  [0.56, 1.11]; concordance 57.1%; Spearman rho 0.808, p < 0.001; Kendall tau 0.445,
  p = 0.006. FAIRD3 vs FAIRD4 (Meadow, n = 22): R² = 0.733, p < 0.001; slope 1.15
  [0.83, 1.48]; concordance 71.4%; rho 0.837, p < 0.001; tau 0.601, p < 0.001.
- PC V16: same values (H1 summary table).
- Output: same values (l.177 to 195).
- Reading: strong replicate agreement in both habitats.

**R-V5-H1Q-FAIRD-02, FAIRD Maize vs Meadow** | LEGACY V5
- Level: habitat (Inter-Habitat; daily habitat series = mean of the two replicates).
  Window: full 22-day series, n = 22.
- Source: V5, Results 3.2.1.3 (Inter-Habitat Level), Tables 8 and 10; method in V5,
  Methods, "H1, Quantitative Consistency", second paragraph.
- Result: R² = 0.719, p < 0.001; slope 0.60 [0.42, 0.77]; concordance 81.0%; rho 0.795,
  p < 0.001; tau 0.640, p < 0.001; pseudo-R² (log) 0.717, low-weight days 3/22.
- PC V16: same values (H1, "Supplementary: Cross-Habitat Consistency"), which adds
  Pearson r = 0.848 and 3/22 days with divergence > 0.5 log units.
- Output: all values identical except low-weight days (l.264 to 272).
- V5 NOT REPRODUCED: low-weight days 3/22 (V5, Table 10) vs 1/22 in the output
  (threshold weight < 0.5). Cause: hypothesis, V5 used a different weight threshold.
- Note: Table 9 and the text of 3.2.1.2 report zero low-weight days only for the
  within-habitat pairs; the between-habitat count is a separate figure.
- Reading: the two habitats share the daily FAIRD abundance signal.

**R-V5-H1Q-AMMOD-01, AMMOD1 vs AMMOD2** | LEGACY V5, DISCREPANCY V5/V16
- Level: Intra-Habitat (Maize). Window: 12 days (first AMMOD batch), n = 12.
- Source: V5, Results 3.2.1.1 and Table 8.
- Result (V5): R² = 0.741, p < 0.001; slope 1.13 [0.66, 1.61]; concordance 63.6%;
  rho 0.874, p < 0.001; tau 0.345, ns.
- PC V16: slope 1.14 [0.66, 1.61] (H1 summary table); all other values identical.
- DISCREPANCY: slope 1.13 (V5) vs 1.14 (PC V16), same interval. Cause: verified,
  rounding of 1.135 (output l.219: 1.135 +/- 0.212 SE, [0.663, 1.607]). Both values
  are correct roundings; V5 kept.
- Output: all other values identical (l.217 to 225).
- Reading: reference level of replicate agreement for the H1 comparison.

**R-V5-H1Q-AMMOD-02, AMMOD Maize devices vs AMMOD3** | LEGACY V5
- Level: Inter-Habitat, device against device (no averaging possible on the Meadow
  side). Window: 8 days (AMMOD3 limited by ATL processing), n = 8.
- Source: V5, Results 3.2.1.3, Tables 8 and 10.
- Result: AMMOD1 vs AMMOD3: R² = 0.472, p = 0.060; slope 0.53 [-0.03, 1.10];
  concordance 57.1%; rho 0.667, ns; tau -0.048, ns; pseudo-R² 0.470, low-weight days 2.
  AMMOD2 vs AMMOD3: R² = 0.351, p = 0.122; slope 0.35 [-0.12, 0.81]; concordance 42.9%;
  rho 0.571, ns; tau -0.333, ns; pseudo-R² 0.350, low-weight days 1.
- PC V16: same values (H1, "Supplementary: Cross-Habitat Consistency"); adds p for rho
  (0.083, 0.151) and tau (0.881, 0.293).
- Output: all values identical except low-weight days (l.288 to 307).
- V5 NOT REPRODUCED: low-weight days 2 and 1 (V5, Table 10) vs 0/8 and 0/8 in the output
  (threshold weight < 0.5). Cause: hypothesis, as in R-V5-H1Q-FAIRD-02. Note that
  PC V16 reports the same 2/8 and 1/8 as V5.
- Reading: no detectable between-habitat agreement for AMMOD at n = 8; not interpretable.

**R-V5-H1Q-ID-01, ID1 vs ID2** | LEGACY V5, DISCREPANCY V5/V16, V5 NOT REPRODUCED
- Level: Intra-Habitat (Maize; ID3 lost, so no Meadow pair).
- Source: V5, Results 3.2.1.1 (text) and Table 8; scenarios defined in V5, Methods,
  "H1, Quantitative Consistency".
- Result, conservative (all 22 days, battery-failure days included as zeros, n = 22):
  R² = 0.027, p = 0.461; slope 0.19 [-0.34, 0.73]; concordance 61.9%; rho 0.284, ns;
  tau 0.202, ns.
- Result, operational (Aug 28 to Sep 1 excluded, n = 17): R² = 0.434, p = 0.004; slope
  0.49 [0.18, 0.80]; concordance 68.8%; rho 0.625, p = 0.007; tau 0.438, p = 0.020.
- PC V16: conservative p = 0.463 (H1, "Results, Insect Detect"); all other values in the
  H1 summary table identical.
- DISCREPANCY: conservative p 0.461 (V5) vs 0.463 (PC V16); no effect on the conclusion.
  Cause: hypothesis, rounding of R² before reporting. With n = 22, any R² that rounds
  to 0.027 gives p between 0.461 and 0.469, so both values are compatible. Not
  verifiable: the ID rows come from the Stage 1 script `H1_consistency_ID.R` (PC V16,
  Final project files and changelog), which is not in ADP; no ADP script computes ID
  OLS pairs (`H1_quanti_consistency.R` assembles Table 8 from FAIRD and AMMOD only).
- Verification: `H1_consistency_summary_output_20261006_114514.txt` (run of 2026-10-06,
  current data). Conservative: R² = 0.0272, slope 0.19 [-0.34, 0.73], p = 0.4633,
  rho 0.284 (p = 0.200), concordance 61.9%, tau 0.202 (p = 0.204) (l.74, 129 to 132,
  216). Operational: R² = 0.434, slope 0.489 [0.181, 0.796], p = 0.0040, rho 0.625,
  concordance 68.8%, tau 0.438 (p = 0.020) (l.217, 220, 230).
- V5 NOT REPRODUCED: conservative p 0.461 (V5) vs 0.4633 (output, = PC V16). Cause:
  verified that the script gives 0.463; origin of 0.461 unknown. No effect on the
  conclusion. All other ID values reproduced.
- Scope note: V5 treats the battery window as zero days; the window has since been
  reclassified as partial, unknown-effort days (R-DX-ID-WIN-01). The conservative
  scenario therefore mixes failure and genuine low activity.
- Reading: ID replicates agree only when the battery window is removed.

**R-V5-H1Q-ROB-01, robust regression check, log scale** | LEGACY V5
- Method: Huber M-estimation (`MASS::rlm`); pseudo-R² = 1 - SS_res,robust / SS_tot,OLS,
  compared with OLS R². Same pairs, metrics and windows as the entries above.
- Source: V5, Results 3.2.1.2 and Table 9 (within-habitat pairs); Table 10
  (between-habitat comparisons, reported in R-V5-H1Q-FAIRD-02 and -AMMOD-02).
- Result (pseudo-R² vs R², delta, low-weight days): FAIRD1-FAIRD2 0.674 vs 0.674, 0.000,
  0/22; FAIRD3-FAIRD4 0.731 vs 0.733, 0.002, 0/22; AMMOD1-AMMOD2 0.740 vs 0.741, 0.001,
  0/12; ID operational 0.433 vs 0.434, 0.001, 0/17; ID conservative -0.013 vs 0.027,
  0.037, 3/22 (the battery-failure days).
- PC V16: same log-scale values for the four pairs it lists (H1, "Robust Regression
  Validation"); it does not give the log-scale conservative ID value.
- Output: FAIRD and AMMOD rows identical (Table09 block, l.240 to 250); ID rows not
  verifiable (script `H1_consistency_ID.R`, not in ADP).
- V5 INTERNAL INCONSISTENCY, ID conservative row: 0.027 - (-0.013) = 0.040, not the
  reported delta 0.037. PC V16 changelog lists only three ID pseudo-R² streams
  (conservative original 0.013, operational original 0.589, operational log 0.433), so
  the log-scale conservative value in V5 has no documented source. Cause: hypothesis,
  the original-scale 0.013 entered Table 9 with a sign change. Supporting evidence:
  `H1_consistency_summary.R` (provided 2026-10-06), which computes the ID pairs, runs
  the within-habitat robust regression only on the original scale (l.184, `rlm(y ~ x)`
  on raw values); log-scale robust regression appears only in the cross-habitat
  function (l.838). Its output of 2026-10-06 gives ID conservative pseudo-R² = 0.013
  (original scale, l.201). No code found that yields a log-scale conservative ID
  pseudo-R². V5 NOT REPRODUCED (add to V6 corrections).
- Reading: OLS R² not driven by influential days, except ID conservative.

**R-V16-H1Q-ROB-01, robust regression check, original scale** | LEGACY PC V16
- Not in V5. Source: PC V16, H1, "Robust Regression Validation".
- Metric: same pairs as above, untransformed daily values.
- Result (pseudo-R², low-weight days): FAIRD1-FAIRD2 0.933; ID conservative 0.013, 3/22;
  ID operational 0.589, 3/17; AMMOD1-AMMOD2 0.588, 3/12. FAIRD3-FAIRD4 not reported.
- Reading: on the raw scale low-weight days appear for ID and AMMOD; the log scale
  removes them. Not to be confused with R-V5-H1Q-ROB-01 (ID conservative -0.013 there is
  log scale).

**R-TB-INT-01, parallelism within each system (constant ratio)** | CURRENT
- Source: `Scherber_temporal_biomass_5.R`, Section 2g. Same GLMM, cascade, PAIRED_BLOCK
  pairing and LRT as Section 2e; units of one system replace Device_type. Maize-only
  comparisons may use last-resort levels without Site/Device intercepts. Unadjusted.
- Result (LRT p): FAIRD halves {1,3} vs {2,4} 0.94 / 0.78 / not computed and {1,4} vs
  {2,3} 0.47 / 0.70 / 0.50 (A_12d / B_22d full / B_22d excl); FAIRD1 vs FAIRD2 0.59 /
  0.51 / 0.14; AMMOD1 vs AMMOD2 0.31 (A_12d, level fixed, Wald 0.27); ID1 vs ID2 0.71 /
  0.99 / 0.99 (B_22d excl Wald 0.56).
- Correction (2026-10-28): until this date the entry gave FAIRD1 vs FAIRD2 0.59 / 0.99 /
  0.99 and ID1 vs ID2 0.71 / 0.14 / 0.99, which match neither
  `Scherber_biomass_internal_consistency.csv` nor the 2026-10-28 run of
  `3-1a_quantity_trajectories.R` (both give the values above). Transcription error;
  the reading is unchanged (no internal comparison rejects parallelism).
- Flag (v0.6): the Maize device-pair rows (FAIRD1 vs FAIRD2, AMMOD1 vs AMMOD2, ID1 vs
  ID2) are NOT FOR MAIN TEXT (C2). The FAIRD halves combine devices across both blocks.
- Reading: no internal comparison rejects parallelism; the non-proportionality between
  systems (R-TB-AF-01, R-TB-FI-22F, R-TB-FI-22X) does not appear within any system.

**R-CC-INT-01, concordance within each system** | CURRENT
- Source: `Scherber_temporal_concordance_biomass.R`. LOW = r between fitted unit curves
  (raw scale, coefficient-simulation interval, lower limit pessimistic); HIGH = r between
  unit series of Pearson residuals (no-AR1 model), lag 0; TOTAL = r of daily unit means.
- Result: FAIRD halves LOW 0.85 to 0.99, HIGH 0.70 to 0.88, TOTAL 0.90 to 0.95;
  FAIRD1 vs FAIRD2 LOW 0.63 to 0.99, HIGH 0.70 to 0.86, TOTAL 0.88 to 0.94;
  AMMOD1 vs AMMOD2 (A_12d) LOW 0.85 [0.12, 0.99] (level fixed), HIGH 0.67, TOTAL 0.80;
  ID1 vs ID2 LOW 0.09 / 0.99 / 0.95, HIGH 0.20 / 0.16 / 0.42, TOTAL 0.08 / 0.25 / 0.69
  (A_12d / B_22d full / B_22d excl).
- Rerun of 2026-10-12: all values reproduced (see the common note of Section 2.2).
- Flag (v0.6): the Maize device-pair rows are NOT FOR MAIN TEXT (C2); FAIRD halves
  over both blocks are the system-level internal reference.
- Reading: FAIRD consistent in all components; AMMOD reasonably consistent (Maize,
  12 days); ID poorly consistent day to day, mainly on partial-operation days
  (consistent with R-V5-H1Q-ID-01).

**R-V16-H1Q-DIAG-01, transformation diagnostic, FAIRD1 vs FAIRD2** | LEGACY PC V16
- Not reported with values in V5. Source: PC V16, "Statistical Methodology", point 1
  (Logarithmic transformation); reproduced by `H1_quanti_consistency_output.txt`,
  Section 1 (l.159 to 172).
- Metric: daily abundance, Maize pair, n = 22.
- Result: raw scale R² = 0.934, Shapiro-Wilk p = 0.134, Breusch-Pagan p = 0.028, DW
  1.521, Cook's D 3/22 influential; log10(x+1) scale R² = 0.674, Shapiro-Wilk p = 0.216,
  Breusch-Pagan p = 0.392, DW 1.645. PC V16 adds 1/22 influential points on the log
  scale (not printed in the output).
- Reading: the raw scale fails homoscedasticity; basis of the log10(x+1) decision for
  all Stage 2 regressions.

**R-OUT-H1Q-01, Spearman rho on daily log changes** | OUTPUT ONLY
- Source: `H1_quanti_consistency_output.txt`, l.182, 192, 222, 269, 293, 304.
- Metric: Spearman rho between day-to-day changes of the log10(x+1) series (V5 reports
  rho on absolute log values instead).
- Result: FAIRD1-FAIRD2 0.573 (p < 0.01); FAIRD3-FAIRD4 0.727 (p < 0.001); AMMOD1-AMMOD2
  0.409, ns; FAIRD Maize vs Meadow 0.765 (p < 0.001); AMMOD1-AMMOD3 -0.071, ns;
  AMMOD2-AMMOD3 -0.500, ns.
- Caution: the "Intra-Habitat Summary" lines of the output (l.204 to 205, 230, 277)
  print this daily-change rho under the label "Spearman rho", while the manuscript
  tables use rho on absolute values. Quote with the metric stated.

**R-OUT-H1Q-02, OLS residual diagnostics of the H1 pairs** | OUTPUT ONLY
- Source: `H1_quanti_consistency_output.txt`, robust regression table (l.240 to 250).
- Metric: log10(x+1) OLS models of R-V5-H1Q-FAIRD-01 and -AMMOD-01.
- Result (Shapiro-Wilk p / Breusch-Pagan p / DW, DW p): FAIRD1-FAIRD2 0.216 / 0.392 /
  1.65, 0.172; FAIRD3-FAIRD4 0.411 / 0.174 / 1.42, 0.056; AMMOD1-AMMOD2 0.702 / 0.225 /
  2.03, 0.442.
- Reading: no assumption violation detected; FAIRD3-FAIRD4 residual autocorrelation
  borderline.

### 1.2 Taxonomic

| ID | Status | Note |
| --- | --- | --- |
| R-V5-H1T-01 | LEGACY V5 | Intra-habitat similarity, 3-day windows, Superfamily (entry below) |
| R-V5-H1T-02 | LEGACY V5 | Inter-habitat similarity, 3-day windows, Superfamily (entry below) |
| R-V5-H1T-03 | LEGACY V5 | AMMOD daily intra-habitat baseline, Species to Superfamily (entry below) |
| R-V5-H1T-04 | LEGACY V5 | AMMOD daily inter-habitat similarity, pooled habitats (entry below) |
| R-OUT-H1T-01 | OUTPUT ONLY | FAIRD M-H at Order and Genus levels, daily and windowed (entry below) |
| R-OUT-H1T-02 | OUTPUT ONLY | Temporal trend of FAIRD windowed M-H (entry below) |
| R-OUT-H1T-03 | OUTPUT ONLY | AMMOD daily inter-habitat I-P at Superfamily level (entry below) |
| R-OUT-H1T-04 | OUTPUT ONLY | Intra vs inter-habitat and habitat tests, Wilcoxon (entry below) |

Common to the R-V5-H1T entries (Stage 2 framework, `H1_quali_consistency.R`):
- Indices: Morisita-Horn (M-H) for FAIRD, on abundance counts; Izsak-Price (I-P) for
  AMMOD, on presence/absence with taxonomic distances (`vegan::taxa2dist`). Source: V5,
  Methods, "H1, Qualitative Consistency". The two indices are not on the same scale;
  FAIRD vs AMMOD values are compared side by side, not tested against each other.
- Windowed entries: non-overlapping 3-day windows, Superfamily level; windows with fewer
  than 3 FAIRD individuals or 3 AMMOD OTUs per device excluded (2 of 14 device-windows
  for FAIRD, 2 of 14 for AMMOD). Source: same V5 section.
- Daily entries: AMMOD only, no windowing.
- PC V16 has no counterpart: in Stage 1, H1 was quantitative only and AMMOD taxonomic
  consistency was not evaluated (PC V16, "H1: TEMPORAL CONSISTENCY", methodology note).
- Verification: `H1_quali_consistency_output.txt` (run of 2026-07-22). All values in
  this section are reproduced by the output.

**R-V5-H1T-01, intra-habitat similarity, 3-day windows** | LEGACY V5
- Level: Intra-Habitat (replicate pair within one habitat). Metric: community
  composition, Superfamily. Window: full season in 3-day windows.
- Source: V5, Results 3.2.2.1 and Table 11.
- Result: FAIRD (M-H) Maize 88.6% +/- 17.3 (n = 5 windows), Meadow 82.9% +/- 13.6
  (n = 7); Student t(10) = 0.641, p = 0.536 (Levene p = 0.797). AMMOD (I-P) Maize
  87.9% +/- 5.3 (n = 7), Meadow 87.8% +/- 4.2 (n = 5).
- Output: same values (FAIRD l.229 to 230 and l.251 to 254; AMMOD windowed summary).
- Reading: high replicate similarity for both systems; no habitat difference for FAIRD.

**R-V5-H1T-02, inter-habitat similarity, 3-day windows** | LEGACY V5
- Level: habitat (Maize pool vs Meadow pool). Metric and window as R-V5-H1T-01.
- Source: V5, Results 3.2.2.2, last paragraph (values in text only, no table).
- Result: FAIRD (M-H) 86.1% +/- 12.1 (n = 7 windows); AMMOD (I-P) 86.6% +/- 2.3
  (n = 7).
- Output: same values (l.392 and l.399).
- Note: this is the value cited in the draft email to Scherber (v0.3 note).
- Reading: between-habitat similarity of the same order as within-habitat similarity,
  for both systems.

**R-V5-H1T-03, AMMOD daily intra-habitat baseline** | LEGACY V5
- Level: Intra-Habitat (AMMOD1 vs AMMOD2 Maize; AMMOD3 vs AMMOD4 Meadow). Metric: I-P,
  presence/absence. Window: daily; Maize n = 21 days, Meadow n = 13 (AMMOD4 limited by
  the mechanical failure, Aug 28 to Sep 3).
- Source: V5, Results 3.2.2.1, last sentence, and Table 12.
- Result (Maize / Meadow): Species 74.1% +/- 6.8 / 78.8% +/- 3.5; Genus 76.5% +/- 6.8 /
  80.9% +/- 3.8; Superfamily 84.3% +/- 5.8 / 87.1% +/- 3.6.
- Output: same values (daily summary block).
- Reading: baseline of natural spatial turnover for the reference system; similarity
  rises with coarser taxonomic level.

**R-V5-H1T-04, AMMOD daily inter-habitat similarity** | LEGACY V5
- Level: habitat (daily detections pooled: AMMOD1+2 vs AMMOD3+4). Metric: I-P. Window:
  daily, n = 21 days.
- Source: V5, Results 3.2.2.2, first paragraph.
- Result: Species 75.0% +/- 5.0; Genus 77.1% +/- 4.8.
- Output: same values (Inter-Habitat Daily Summary); the Superfamily value of the same
  run is recorded in R-OUT-H1T-03.
- Reading: inter-habitat similarity close to the intra-Maize baseline (74.1% and 76.5%
  in R-V5-H1T-03).

Common to the R-OUT-H1T entries: source `H1_quali_consistency_output.txt`; indices,
windows and filters as in the R-V5-H1T common note.

**R-OUT-H1T-01, FAIRD M-H at Order and Genus levels** | OUTPUT ONLY
- Source: output, Sections 1.5 and windowed summary (l.210 to 232).
- Result, 3-day windows filtered (mean +/- SD, n windows; Maize / Meadow): Order 91.1
  +/- 17.3, n = 6 / 98.5 +/- 1.9, n = 7; Genus 83.7 +/- 23.5, n = 5 / 63.8 +/- 23.4,
  n = 7. Superfamily as R-V5-H1T-01.
- Result, daily unfiltered (mean +/- SD, n days; Maize / Meadow): Order 87.7 +/- 27.0,
  n = 15 / 93.8 +/- 7.0, n = 22; Superfamily 78.7 +/- 37.2, n = 12 / 61.5 +/- 33.7,
  n = 22; Genus 71.5 +/- 37.6, n = 11 / 51.0 +/- 35.5, n = 18.
- Reading: FAIRD replicate similarity falls with finer taxonomic level and is very
  variable at daily resolution, especially in Meadow.

**R-OUT-H1T-02, temporal trend of FAIRD windowed M-H, Superfamily** | OUTPUT ONLY
- Source: output, Section 1.6 (l.242 to 247). Linear regression of M-H on window index.
- Result: pooled slope -6.19 pp per window, R² = 0.653, p = 0.0015, n = 12; FAIRD1-FAIRD2
  -8.43 pp, R² = 0.595, p = 0.127, n = 5; FAIRD3-FAIRD4 -5.44 pp, R² = 0.746, p = 0.012,
  n = 7.
- Reading: replicate similarity declines over the season (consistent with falling
  FAIRD capture counts, V5 Methods, H1 Qualitative Consistency).

**R-OUT-H1T-03, AMMOD daily inter-habitat I-P, Superfamily** | OUTPUT ONLY
- Source: output, Section 4.1, Inter-Habitat Daily Summary. Pooled AMMOD1+2 vs AMMOD3+4.
- Result: 83.6% +/- 5.2, n = 21 days (Species and Genus in R-V5-H1T-04).

**R-OUT-H1T-04, intra vs inter-habitat and habitat tests (Wilcoxon)** | OUTPUT ONLY
- Source: output, Section 5 (l.389 to 410). 3-day windows, Superfamily.
- Result: FAIRD intra (n = 12, mean 85.2%) vs inter (n = 7, 86.1%): W = 45, p = 0.833,
  d = -0.064. AMMOD intra (n = 12, 87.9%) vs inter (n = 7, 86.6%): W = 52, p = 0.422,
  d = 0.322. FAIRD Maize vs Meadow: W = 22, p = 0.516, d = 0.376. AMMOD Maize vs Meadow:
  W = 18, p = 1.000, d = 0.016.
- Note: V5 reports the FAIRD habitat comparison with Student's t (R-V5-H1T-01); the
  Wilcoxon version is not in V5. Same conclusion.
- Reading: no intra vs inter-habitat difference for either system; supports the
  side-by-side reading of R-V5-H1T-01 and -02.

---

## 2. AGREEMENT BETWEEN SYSTEMS

### 2.1 Temporal biomass trajectories, pair models (parallelism, log link)

Common to all entries: `Scherber_temporal_biomass_3.R`, Section 2e, run of 2026-10-02.
Model bs(Date.num, 3) * Device_type + Block fixed, (1|Site/Device), AR1, Tweedie, log
link; interaction 3 df. Parallel on the log link = constant ratio between systems. Full
table with all variants: `outputs/tables/Scherber_biomass_pairs.csv`. The same results
are reproduced by `Scherber_temporal_biomass_5.R` (Section 2e unchanged).

---

**R-TB-ALL-01, omnibus three systems** | CURRENT
- Question: are the three biomass trajectories parallel?
- Data: 12-day AMMOD window (Aug 23 to Sep 3), FULL, ID full series.
- Result: Date x Device_type p = 0.0063 (Wald, 6 df).
- Reading: not all three trajectories are parallel.

**R-TB-AF-01, AMMOD vs FAIRD, primary** | CURRENT
- Data: A_12d, PAIRED_BLOCK.
- Result: LRT p = 0.013, Holm 0.039.
- Sensitivities: every variant p < 0.05 (range 0.005 to 0.049); post hoc removal of the
  first 3 days p = 0.049.
- Reading: not parallel; robust to the pairing criterion; the only pair surviving Holm.

**R-TB-AF-00, AMMOD vs FAIRD, co-located site-days** | SUPERSEDED by R-TB-AF-01 as
primary (kept as sensitivity)
- Data: 12 days, PAIRED by site (now PAIRED_SITE).
- Result: Wald p = 0.015, LRT p = 0.028 (LRT at level `simple`). Reproduced exactly on
  2026-10-02.

**R-TB-AI-01, AMMOD vs ID, primary** | CURRENT
- Data: A_12d, PAIRED_BLOCK, 64 obs (Maize: AMMOD1-2 vs ID1-2, 12 days; Meadow: AMMOD3 vs
  ID4, 8 days).
- Result: LRT p = 0.083, Holm 0.085; Wald 0.034 at the same level.
- Sensitivities: range 0.011 to 0.049 (PAIRED_SITE is Maize only, Block not estimable).
- Reading: not significant after correction; Wald anti-conservative at this n; not
  evidence of parallelism.

**R-TB-FI-01, FAIRD vs ID, 12 days, primary** | CURRENT
- Data: A_12d, PAIRED_BLOCK (identical to FULL), ID full series.
- Result: LRT p = 0.043, Holm 0.085; PAIRED_SITE p = 0.091.
- Reading: not significant after correction; fragile to the pairing criterion.

**R-TB-FI-22F, FAIRD vs ID, 22 days, full series** | CURRENT (complementary)
- Data: B_22d, PAIRED_BLOCK, all days.
- Result: LRT p = 0.0008.
- Reading: not parallel over the season.

**R-TB-FI-22X, FAIRD vs ID, 22 days, window excluded** | CURRENT (complementary)
- Data: B_22d, PAIRED_BLOCK, Aug 28 to Sep 1 removed for both systems (17 days); AR1
  valid.
- Result: LRT p < 0.0001.
- Reading: not parallel; the adverse window does not explain the non-parallelism.
- Scale note: on a linear axis the fitted curves look nearly identical (difference about
  100 to 140 mg) while their ratio goes from about 1.4 to 10. Proportionality is
  rejected; shared temporal shape is not tested here (see 2.2).
  [SUPERSEDED in part (v0.3): the difference is not constant either, see
  R-TB-ALL-ID-01; shared shape tested in R-CC-FI-22-01.]

**R-TB-ALL-ID-01, constant-difference test (identity link), all pairs** | CURRENT
- Source: `Scherber_temporal_biomass_5.R`, Section 2f. Same runs and data as above; only
  the link changes. AMMOD divided by the dual-axis figure factor k. Cascade: Tweedie
  identity, then Gaussian identity with dispersion by system.
- Result (LRT, Holm over the A_12d primary pairs): AMMOD-FAIRD 0.018 (Holm 0.054;
  Tweedie, no Block, no AR1); AMMOD-ID 0.058 (0.116; Tweedie, no AR1); FAIRD-ID A_12d
  0.106 (0.116; Tweedie, AR1); FAIRD-ID B_22d full 0.80 (Gaussian, no AR1); B_22d excl
  0.45 (Gaussian, AR1).
- Caveats: FAIRD-ID B_22d on the Gaussian fallback, which predicts negative FAIRD biomass
  in the second half; fitted difference FAIRD minus ID from -71 to -259 mg. Only FAIRD-ID
  A_12d isolates the link change cleanly. Identity-link figures not for the manuscript.
- Reading: neither constant ratio nor constant difference describes the systems;
  non-rejection under the identity link is lack of precision, not parallelism.

### 2.2 Temporal concordance, system level (shared signal)

| ID | Status | Note |
| --- | --- | --- |
| R-CC-AF/AI/FI-* | DONE, see entries below | Two components per pair: shared low-frequency shape, shared day-to-day fluctuation. Scale to be declared per claim |
| R-CC-DESC-01 | CURRENT, DESCRIPTIVE | Simple correlations of daily system means, biomass and abundance, Pearson and Spearman (entry below) |
| R-V5-H2Q-01 | LEGACY V5, DISCREPANCY V5/V16 | FAIRD vs AMMOD, site-level OLS, daily biomass (entry below) |
| R-V5-H2Q-02 | LEGACY V5, DISCREPANCY V5/V16, V5 NOT REPRODUCED | FAIRD vs AMMOD, Maize pooled, n = 24 (entry below) |
| R-V5-H2Q-03 | LEGACY V5, DISCREPANCY V5/V16, V5 NOT REPRODUCED | FAIRD vs AMMOD, Overall pooled, n = 32 (entry below) |
| R-V5-H2Q-04 | LEGACY V5, DISCREPANCY V5/V16 | FAIRD vs AMMOD cross-correlation (entry below) |
| R-V5-H2Q-05 | LEGACY V5 | Detection metric sensitivity, counts and dry mass (entry below) |
| R-V5-H2Q-06 | LEGACY V5 | AMMOD body length and taxonomic DeltaS by site, V5 Table 15 (entry below) |
| R-V16-H2Q-01 | LEGACY PC V16 | Pooling-validity covariate models (entry below) |
| R-V16-H2Q-02 | LEGACY PC V16 | OLS residual diagnostics and adjusted R², all levels (entry below) |
| R-V16-H2Q-03 | LEGACY PC V16 | Taxonomic cross-habitat tests behind V5 Table 15 (entry below) |
| R-V16-H2Q-AI-01 | SUPERSEDED by R-TB-AI-01 and R-CC-AI-01 | ID vs AMMOD exploratory correlation, n = 7 (entry below) |
| R-OUT-H2Q-01 | OUTPUT ONLY | Spearman rho on daily log changes, all levels (entry below) |
| R-OUT-H2Q-02 | OUTPUT ONLY | Robust regression check, all levels (entry below) |
| R-OUT-H2Q-03 | OUTPUT ONLY | Overall CCF and Fisher z average (entry below) |
| R-OUT-H2Q-04 | OUTPUT ONLY | FAIRD body-length classes by site (entry below) |
| (moved) | see Section 2.5 | Diagnostic mixed models A and B (R-V5-H2Q-MM-01, -02) |

Common to the entries below: `Scherber_temporal_concordance_biomass.R`, spline df 3,
PAIRED_BLOCK unless stated; definitions as in R-CC-INT-01; HIGH intervals by Fisher z
(indicative). A bootstrap p-value for HIGH was implemented and removed (run never
finished); its partial values are not for reporting.
Rerun of 2026-10-12 (`outputs/console/Scherber_temporal_concordance_biomass_output.txt`,
3.8 min; until v0.6 cited as `claude/Scherber_temporal_concordance_biomass_output_20261012.txt`,
same run): every
point estimate of R-CC-INT-01, R-CC-AF-01 to R-CC-SENS-01 and every parallelism LRT of
R-TB-INT-01 reproduced. One interval differs from the logged value: ID1 vs ID2 B_22d
full LOW 0.991 [-0.558, 0.982] here vs [-0.60, 0.988] in the session log of 2026-10-05.
The script fixes the seed (`set.seed(20261002)`, 1000 coefficient draws), so the cause
is not determined (possibly a different script version on 2026-10-05). The estimate
above the upper limit persists (Section 6).
Level of comparison for the manuscript (v0.6): the two-block rows are the primary,
system-level results. Maize-only rows compare systems within one block (secondary).
Disattenuation rests on Maize device-pair references (C2, R-OBS-05).

**R-CC-AF-01, AMMOD vs FAIRD, A_12d** | CURRENT
- Result: LOW 0.25 [-0.39, 0.78]; HIGH 0.61 [0.05, 0.88]; TOTAL 0.38. Maize only: LOW
  0.05; HIGH 0.87 [0.59, 0.96]; TOTAL 0.66.
- Reading: shape over 12 days not estimable; shared day-to-day fluctuation, part of it
  the common response to the cold spell. Primary manuscript value: HIGH two-block 0.61
  (R-OBS-05).

**R-CC-AI-01, AMMOD vs ID, A_12d** | CURRENT
- Result: LOW -0.20; HIGH 0.65; TOTAL 0.20. Maize only: LOW -0.72; HIGH 0.53; TOTAL 0.05.
- Reading: daily agreement similar to AMMOD-FAIRD; shape not estimable.

**R-CC-FI-01, FAIRD vs ID, A_12d** | CURRENT
- Result: LOW -0.46; HIGH 0.12; TOTAL -0.07. Maize only: LOW -0.44; HIGH 0.33; TOTAL 0.13.
- Reading: no shared signal over 12 days.

**R-CC-FI-22-01, FAIRD vs ID, B_22d** | CURRENT, LOW NOT ROBUST (see R-CC-SENS-01)
- Result: full LOW 0.92 [-0.45, 0.97], HIGH -0.02, TOTAL 0.15; window excluded LOW 0.93
  [0.51, 0.98], HIGH -0.34, TOTAL 0.11; excluded, Maize only LOW 0.94 [0.65, 0.97],
  HIGH -0.13, TOTAL 0.35.
- Reading: only the timing of the rise and peak is shared (clearest without the
  partial-operation days); curves diverge in the second half; no shared day-to-day
  fluctuation.

**R-CC-DIS-01, disattenuated HIGH, Maize only (Spearman-Brown)** | CURRENT, NOT FOR MAIN
TEXT (C2)
- Result: AMMOD-FAIRD 1.01; AMMOD-ID 1.04 (unstable, ID reliability 0.33); FAIRD-ID 0.61 /
  0.23 / -0.18 (A_12d / B_22d full / B_22d excl).
- Reading rules: > 1 can reflect site-level signal shared by co-located systems; < 1
  mixes detection with recording-coverage differences.
- Flag (v0.6): reliabilities come from device pairs (AMMOD1 vs AMMOD2, FAIRD1 vs
  FAIRD2, ID1 vs ID2); diagnostic or appendix only (R-OBS-05).
- Reading: AMMOD-FAIRD share the daily signal as far as reliability allows; FAIRD-ID
  disagreement not explained by ID noise alone.

**R-CC-SENS-01, sensitivity spline df 4** | CURRENT
- Result: LOW FAIRD-ID B_22d 0.92 to 0.09 (full), 0.93 to 0.40 (excl); AMMOD-FAIRD 0.25
  to 0.53; FAIRD-ID A_12d -0.46 to -0.02. HIGH AMMOD-FAIRD 0.61 to 0.53; AMMOD-ID 0.65 to
  0.59.
- Reading: LOW depends on spline flexibility; HIGH is stable.

**R-CC-DESC-01, simple correlations of daily system means** | CURRENT, DESCRIPTIVE
- Source: `descriptive_correlations_daily_means.R`, run of 2026-10-28; output
  `outputs/console/descriptive_correlations_daily_means_output.txt`, table
  `outputs/tables/descriptive_correlations_daily_means.csv`. No models, no tests.
- Metric: daily unit means averaged within block, then across blocks (definition of
  TOTAL); Pearson (raw scale) and Spearman. Same comparisons, windows, block-day
  pairing and symmetric window exclusion as R-CC-*. Biomass for all systems, abundance
  for FAIRD and ID.
- Check: biomass Pearson reproduces r_total of `Scherber_concordance_biomass.csv` for
  all 23 comparisons (largest difference 1.1e-16).
- Result, Spearman (Pearson), biomass / abundance:
  - FAIRD halves, mean of two splits, A_12d / B_22d full / B_22d excl: biomass 0.74 /
    0.85 / 0.82 (0.91 / 0.93 / 0.95); abundance 0.83 / 0.88 / 0.84 (0.96 / 0.97 / 0.98).
  - FAIRD1 vs FAIRD2: biomass 0.71 / 0.79 / 0.79 (0.88 / 0.90 / 0.93); abundance 0.82 /
    0.81 / 0.78 (0.96 / 0.97 / 0.98).
  - AMMOD1 vs AMMOD2, A_12d: 0.87 (0.80).
  - ID1 vs ID2, A_12d / B_22d full / B_22d excl: biomass 0.16 / 0.27 / 0.60 (0.08 /
    0.25 / 0.69); abundance 0.21 / 0.28 / 0.63 (0.09 / 0.24 / 0.78).
  - AMMOD vs FAIRD, A_12d: two blocks 0.33 (0.38); Maize 0.33 (0.66).
  - AMMOD vs ID, A_12d: two blocks 0.42 (0.20); Maize 0.21 (0.05).
  - FAIRD vs ID, A_12d: two blocks biomass 0.03 (-0.07), abundance 0.06 (-0.08).
  - FAIRD vs ID, B_22d full / excl, two blocks: biomass 0.24 / 0.28 (0.15 / 0.11);
    abundance 0.22 / 0.19 (0.13 / 0.11). Maize: biomass 0.58 / 0.43 (0.31 / 0.35);
    abundance 0.58 / 0.44 (0.26 / 0.36).
- Reading: raw-scale Pearson is dominated by the shared seasonal peak (FAIRD internal
  Pearson close to LOW); Spearman weights all days equally. Within FAIRD, abundance is
  more consistent than biomass (the mass estimates add noise). Between FAIRD and ID,
  abundance agrees no better than biomass.
- Manuscript use (Juan's decision, 2026-10-28): Spearman of daily means in parentheses
  next to the model-based results, for every comparison; Pearson and both measures in
  Supplementary Material S2.

**H2 quantitative agreement, FAIRD vs AMMOD (Stage 2 framework)**

Common to R-V5-H2Q-01 to -05 and the R-OUT-H2Q entries (`H2_quanti_validation.R`):
- Pair: FAIRD vs the co-located AMMOD; Site 4 excluded (no AMMOD4 biomass). Window:
  Sites 1 and 2, 12 days (Aug 23 to Sep 3); Site 3, 8 days (Aug 23 to 30).
- Metric: daily biomass (mg), FAIRD inferred live mass (allometric), AMMOD live mass.
- Transformation: V5 states "log-transformed"; log10(x+1), AMMOD as predictor and FAIRD
  as response (script l.28 and l.67 to 68).
- Statistics: OLS R² and slope; Spearman rho on absolute log values; Kendall tau and
  concordance on raw daily changes (script l.106); for pooled levels, daily changes
  computed within site before pooling; bootstrap R² 1000 resamples, seed 42, quantile
  CI (script l.302).
- Comparison with PC V16 made against its H2 summary table, "Section 1-3", "Section 9",
  "Section 10" and "Section 11" (Stage 1 scripts `H2_validation_*.R`, not in ADP).
- Verification: `H2_quanti_validation_output.txt` (run of 2026-07-22). The rerun of
  2026-10-08 reproduces every value of its Sections 1 to 7; line references are kept to
  the 2026-07-22 file.
- File in use (v0.7): the 2026-07-22 file is no longer in ADP; the repository holds the
  2026-10-08 rerun (`outputs/console/H2_quanti_validation_output.txt`). Its Sections 1
  to 7 are the 2026-07-22 lines shifted by +3 (three Truth Calendar lines added by
  `00_load_data.R`; checked at l.170, 205, 209, 228 and 231, cited here as l.167, 202,
  206, 225 and 228); add 3 to every line reference below when reading that file. Its
  Section 8 (registry diagnostics) also appears there from l.359, because the run
  predates the `split = FALSE` fix; the separate diagnostics file is
  `outputs/console/H2_quanti_registry_diagnostics_output.txt`.
- Diagnostics: `H2_quanti_registry_diagnostics_output.txt` ("diag" below; Section 8 of
  the script, run of 2026-10-08), tests of the hypotheses recorded in v0.4.

**R-V5-H2Q-01, FAIRD vs AMMOD, site level** | LEGACY V5, DISCREPANCY V5/V16
- Level: Site (one FAIRD vs one AMMOD).
- Source: V5, Results 3.3.1.1 and Table 13.
- Result: Site 1 (FAIRD1 vs AMMOD1, Maize, n = 12): R² = 0.282, p = 0.076; slope 1.44
  [-0.18, 3.05]; concordance 45.5%; rho 0.238, ns; tau 0.200, ns. Site 2 (Maize,
  n = 12): R² = 0.445, p = 0.018; slope 0.86 [0.18, 1.54]; concordance 81.8%; rho
  0.566, ns; tau 0.673, p = 0.004. Site 3 (Meadow, n = 8): R² = 0.058, p = 0.566; slope
  0.33 [-1.01, 1.68]; concordance 28.6%; rho 0.214, ns; tau -0.238, ns.
- PC V16: same values in the H2 summary table and "Section 1-3" (adds rho p 0.457,
  0.0548, 0.610). The "Summary table, all aggregation levels" (end of "Section 11")
  gives other values: tau 0.152 ns / 0.455 * / 0.071 ns and Site 3 slope 0.34.
- DISCREPANCY: (a) tau 0.200 / 0.673 / -0.238 (V5) vs 0.152 / 0.455 / 0.071 (PC V16,
  Section 11 table). Cause: verified, the PC V16 site rows are Kendall tau on absolute
  log values, not on daily changes (diag 8B rows 3, 6, 9: 0.152, 0.455, 0.071); its
  pooled rows (0.359, 0.187) are on raw daily changes, as in V5. The output reproduces
  V5. (b) Site 3 slope 0.33 (V5) vs 0.34 (PC V16).
  Cause: verified, rounding of 0.335 (output l.182).
- Output: all V5 values reproduced (l.167 to 186; Table 13 assembly l.307 to 309).
- Reading: significant covariation of daily log biomass only at Site 2; Site 1
  inconclusive; Site 3 underpowered (n = 8).

**R-V5-H2Q-02, FAIRD vs AMMOD, Maize pooled** | LEGACY V5, DISCREPANCY V5/V16, V5 NOT
REPRODUCED
- Level: Intra-Habitat (Sites 1 and 2 stacked, n = 24 site-days).
- Source: V5, Results 3.3.1.2 and Table 13.
- Result (V5): R² = 0.311, p = 0.005; bootstrap 95% CI [0.033, 0.554]; slope 1.08
  [0.37, 1.79]; rho 0.482, p = 0.017; tau 0.359, p = 0.011; concordance 63.6%.
- PC V16 ("Section 10"; "H2, New: Habitat and Overall Aggregation Levels"): same R²,
  slope, rho, concordance and CI (bootstrap point estimate 0.312); p = 0.006; tau
  p = 0.019; DW 1.385, p = 0.041.
- DISCREPANCY: p 0.005 (V5) vs 0.006 (PC V16). Cause: verified, PC V16 reports the
  log_ammod p of the covariate model (output l.201: 0.0060), not the simple OLS p
  (output l.317: 0.0047, = V5).
- V5 NOT REPRODUCED: (a) tau p 0.011 (V5) vs 0.0193 (output l.324, = PC V16); tau
  identical. Cause: not determined. Variants tested (diag 8B rows 10 to 15): normal
  approximation 0.0193, with continuity correction 0.0208, exact 0.0194, one-sided
  0.0096, log changes 0.0259, absolute log values 0.0131 (tau 0.362); none gives
  0.011. (b) bootstrap CI [0.033, 0.554] (V5, = PC V16) vs [0.041, 0.556]
  (output l.206). Cause: hypothesis, Stage 1 bootstrap with another seed or resampling
  scheme.
- PC V16 vs output (no V5 value, not a V6 correction): DW 1.385, p = 0.041 (PC V16) vs
  1.845, p = 0.299 (output l.211). Cause: verified, row order: sorted by Date then Site
  the same model gives DW 1.38, p = 0.0409 (diag 8A row 2); the script sorts by Site
  then Date (script l.310). See R-OBS-02.
- Reading: pooled Maize daily log biomass covaries positively; R² interval wide.

**R-V5-H2Q-03, FAIRD vs AMMOD, Overall pooled** | LEGACY V5, DISCREPANCY V5/V16, V5 NOT
REPRODUCED
- Level: Overall (Sites 1 to 3 stacked, n = 32 site-days).
- Source: V5, Results 3.3.1.3 and Table 13.
- Result (V5): R² = 0.321, p < 0.001; bootstrap 95% CI [0.063, 0.543]; slope 1.03
  [0.47, 1.58]; rho 0.491, p = 0.004; tau 0.187, p = 0.154; concordance 55.2%; DW
  0.917, p < 0.001 (reported in V5 as a limitation).
- PC V16 ("Section 11"): same values, except p = 0.0036.
- DISCREPANCY: p < 0.001 (V5) vs 0.0036 (PC V16). Cause: verified, as in R-V5-H2Q-02
  (covariate model output l.224: 0.00364; simple OLS l.318: 0.0007).
- V5 NOT REPRODUCED: (a) bootstrap CI [0.063, 0.543] (V5, = PC V16) vs [0.064, 0.541]
  (output l.228). Cause: hypothesis, as in R-V5-H2Q-02. (b) DW 0.917, p < 0.001 (V5,
  = PC V16) vs 1.864, p = 0.306 (output l.233). Cause: verified, row order: sorted by
  Date then Site the same model gives DW 0.917, p = 0.0002 (diag 8A row 4); with Site
  then Date (script l.458), the order valid for temporal autocorrelation, 1.864. The V5
  limitation sentence rests on the 0.917 value (R-OBS-02).
- Reading: positive covariation across the three sites; directional (tau) agreement
  diluted by Site 3.

**R-V5-H2Q-04, FAIRD vs AMMOD cross-correlation** | LEGACY V5, DISCREPANCY V5/V16
- Method: z-scored raw (untransformed) daily biomass, `stats::acf` on the paired series,
  CI +/- 1.96/sqrt(n) (script l.585 to 597, l.612 to 625). Lag 0 primary.
- Source: V5, Results 3.3.1.4 and Table 13.
- Result (V5): lag 0 Site 1 0.347, ns (CI +/- 0.566); Site 2 0.781, p < 0.05 (+/- 0.566);
  Site 3 -0.002, ns (+/- 0.693); Maize summed series (n = 12) 0.657, p < 0.05. Lag +1
  0.548 / 0.622 / 0.575, none significant; Maize lag +1 0.589, exceeds the CI. Overall
  not computed.
- PC V16 ("Section 9"; H2 summary table): lag 0 0.382 / 0.696 / -0.010; lag +1 (FAIRD
  leading) 0.621 / 0.658 / 0.487; Site 3 lag +2 0.792; AMMOD lag-1 autocorrelation
  0.381 to 0.530. No Maize summed series.
- DISCREPANCY: lag 0 and lag +1 site values. Cause: not determined. The v0.4 hypothesis
  (CCF on log10 series) is rejected: log10(x+1) series give lag 0 0.531 / 0.667 / 0.241
  and lag +1 0.389 / 0.658 / 0.587 (diag 8C); only Site 2 lag +1 coincides.
- Output: lag 0 values reproduced (l.241 to 252); lag +1 reproduced by diag 8C (raw
  series): 0.548 / 0.622 / 0.575, Maize summed 0.589.
- Direction: the script's `acf[, 1, 2]` equals `stats::ccf(faird, ammod)` at every lag
  >= 0 (diag 8C), so lag k = cor(FAIRD[t+k], AMMOD[t]): a positive lag means FAIRD
  later than AMMOD. The PC V16 label "FAIRD leading" cannot be checked, its values not
  being reproduced.
- Not in V5: Site 3 raw lag +2 0.769 exceeds its CI (+/- 0.693) (diag 8C; PC V16
  0.792).
- Reading: same-day synchrony at Site 2 and in the Maize summed series; lag +1 values
  compatible with autocorrelation, no lead-lag claim.

**R-V5-H2Q-05, detection metric sensitivity** | LEGACY V5
- Metric: same OLS, with FAIRD raw counts and AMMOD dry mass ("Alternative") against
  inferred biomass and live mass ("Standard"); DeltaR² = R²(alternative) minus
  R²(standard).
- Source: V5, Results 3.3.1.5 and Table 14.
- Result: R² standard vs alternative: Site 1 0.282 vs 0.171 (DeltaR² -0.111); Site 2
  0.445 vs 0.365 (-0.079); Site 3 0.058 vs 0.029 (-0.029). Alternative slope +/- SE,
  p, rho: Site 1 2.643 +/- 1.838, 0.181, 0.267; Site 2 2.015 +/- 0.840, 0.037, 0.542;
  Site 3 0.398 +/- 0.941, 0.687, 0.119.
- PC V16 (block of `H2_validation_counts_vs_drymass.R`, after "Section 11"): same values.
- Output: R² reproduced (l.265 to 275); slopes, p and rho of the alternative not
  printed. DeltaR² printed with opposite sign convention, and Site 2 as 0.080 instead
  of 0.079. Cause: verified, the script rounds R² to 3 decimals (script l.854) before
  subtracting (l.915), so 0.445 minus 0.365; the V5 value is correct.
- Reading: the inferred-biomass metrics give higher R² than raw counts and dry mass at
  every site.

**R-V5-H2Q-06, AMMOD body length and taxonomic DeltaS by site** | LEGACY V5
- Metric: mean daily family body length (mm) of the AMMOD-detected taxon pool; pairwise
  Mann-Whitney, Bonferroni over 6 pairs. Table 15 also lists daily Izsak-Price DeltaS
  (taxonomic level not stated in V5; genus in PC V16).
- Source: V5, Results 3.3.1.6 (Supplementary) and Table 15.
- Result: body length Site 1 5.670 +/- 0.834; Site 2 5.719 +/- 0.874; Site 3 5.457
  +/- 0.457; Site 4 5.243 +/- 0.355; all p_bonf > 0.16 (closest Site 1 vs Site 4,
  p_raw = 0.027, p_bonf = 0.165). DeltaS: intra-Maize 0.734 +/- 0.078 (n = 21),
  intra-Meadow 0.793 +/- 0.043 (n = 13), Site 1 vs Meadow 0.716 +/- 0.050 (n = 34),
  Site 2 vs Meadow 0.709 +/- 0.080 (n = 34).
- PC V16 ("Supplementary: Body Size Distribution by Site", `BodySize_Distribution.R`;
  "Section 7C", which adds medians and W): same values.
- Output: not verifiable. Section 6 of the output, also labelled "Table 15", is a
  different analysis (R-OUT-H2Q-04); the Stage 1 scripts are not in ADP.
- Reading: no body-length difference between sites in the AMMOD taxon pool.

**R-V16-H2Q-01, pooling-validity covariate models** | LEGACY PC V16
- Source: PC V16, "Section 10" and "Section 11"; output l.198 to 204 and l.221 to 226.
- Model: log_faird ~ log_ammod + Site, log10(x+1).
- Result: Maize SiteSite2 p = 0.861; Overall SiteSite2 p = 0.830, SiteSite3 p = 0.645;
  log_ammod 1.072 (p = 0.0060) and 0.956 (p = 0.0036).
- Output: same values.
- Reading: no site offset detected at fixed AMMOD biomass; used to justify pooling.

**R-V16-H2Q-02, OLS residual diagnostics and adjusted R²** | LEGACY PC V16
- Source: PC V16, "Section 1-3" (sites), "Section 10" and "Section 11" (pooled); output
  l.172, 179, 186, 211, 233; Table 13 assembly l.307 to 311.
- Result (Shapiro-Wilk p / Breusch-Pagan p / DW, DW p): Site 1 0.728 / 0.202 / 1.988,
  0.410; Site 2 0.403 / 0.877 / 1.429, 0.096; Site 3 0.974 / 0.559 / 1.800, 0.264;
  Maize 0.278 / 0.193; Overall 0.089 / 0.082 (pooled DW in R-V5-H2Q-02 and -03).
  Adjusted R²: 0.210 / 0.389 / -0.099; Maize 0.279; Overall 0.298.
- Output: same values.
- Note: per-site DW had been listed as an OUTPUT ONLY candidate; it is in PC V16.
- Reading: no assumption violation detected at site level.

**R-V16-H2Q-03, taxonomic cross-habitat tests** | LEGACY PC V16
- Source: PC V16, "Supplementary: Taxonomic Cross-Habitat Consistency"
  (`TaxConsistency_IntraCrossHabitat.R`), quoted in the Site 1 and Site 3 blocks of
  "Section 1-3". No output in ADP.
- Metric: daily Izsak-Price DeltaS, genus level, `vegan::taxa2dist`; Wilcoxon rank-sum
  and Cohen's d.
- Result: AMMOD Site 1 vs Meadow (0.716) vs Site 2 vs Meadow (0.709): W = 572.0,
  p = 0.946, d = 0.107. AMMOD intra-Maize (0.734) vs intra-Meadow (0.793): W = 60.0,
  p = 0.007, d = 0.886. AMMOD cross-habitat vs intra-habitat (0.756 vs 0.756):
  d = 0.008. FAIRD counterparts: p = 0.749, d = 0.070; p = 0.438, d = 0.195; p = 0.306,
  d = 0.262.
- Relation to V5: Table 15 shows the AMMOD DeltaS means without these tests.
- Reading: higher compositional heterogeneity within Maize than within Meadow (AMMOD);
  no difference between the two Maize sites relative to Meadow.

**R-V16-H2Q-AI-01, ID vs AMMOD, exploratory correlation** | SUPERSEDED by R-TB-AI-01 and
R-CC-AI-01
- Source: PC V16, "H2, Supplementary: ID vs AMMOD Validation (Exploratory)"
  (`H2_validation_ID.R`, not in ADP). Dropped in V5.
- Data: Sites 1 and 2 (ID1 vs AMMOD1, ID2 vs AMMOD2), n = 7 paired days per site
  (Aug 23 to 27, Sep 2 to 3; battery window excluded); Meadow pair AMMOD3 vs ID4 not
  used. Metric: daily biomass, log10.
- Result: Site 1 R² = 0.072, ns; slope 0.62 [-1.93, 3.16]; concordance 83.3%; tau 0.333,
  ns. Site 2 R² = 0.052, ns; slope 0.23 [-0.91, 1.38]; concordance 50.0%; tau 0.333,
  ns. CCF (CI +/- 0.741) maximum at lag 0 for Site 1 (0.269) and at lag -1 for Site 2
  (0.245; lag 0 0.228).
- Replaced by: R-TB-AI-01 and R-CC-AI-01 (block pairing, AMMOD3 vs ID4 included,
  64 obs).

**R-OUT-H2Q-01, Spearman rho on daily log changes** | OUTPUT ONLY
- Source: output l.169, 176, 183, 208, 230.
- Result: Site 1 0.509 (p = 0.110); Site 2 0.636 (p = 0.035); Site 3 0.071 (p = 0.879);
  Maize 0.504 (p < 0.05); Overall 0.370 (p < 0.05).
- Caution: V5 reports rho on absolute log values; quote with the metric stated.

**R-OUT-H2Q-02, robust regression check** | OUTPUT ONLY
- Source: output l.171, 178, 185, 210, 232. Method: `MASS::rlm`, Huber psi, log10(x+1);
  pseudo-R² = 1 minus SS(residual) / SS(total); low-weight days = weight < 0.5 (script
  l.118 to 122).
- Result: pseudo-R² Site 0.280 / 0.445 / 0.048; Maize 0.302; Overall 0.314; low-weight
  days 0 at every level.
- Reading: OLS results not driven by influential days.

**R-OUT-H2Q-03, Overall CCF and Fisher z average** | OUTPUT ONLY
- Source: output l.254 to 257.
- Method: daily mean of FAIRD and of AMMOD biomass across the sites with data that day
  (Sites 1 to 3; Site 3 only up to Aug 30), 12 days, then the CCF of R-V5-H2Q-04
  (script l.779 to 788).
- Result: Overall lag 0 0.520, CI +/- 0.566, ns (labelled descriptive). Fisher z
  average of Site 1 and Site 2 lag 0: 0.608 (descriptive).
- Caution: V5 states the Overall CCF was not computed (Methods, H2 Quantitative; Table
  13 footnote).

**R-OUT-H2Q-04, FAIRD body-length classes by site** | OUTPUT ONLY
- Source: output Section 6, "Table 15", l.283 to 291 (labelled Table 15 but unrelated to
  V5 Table 15).
- Result (% of individuals < 3 mm / 3 to 5 mm / > 5 mm): Site 1 (n = 265) 0 / 48.3 /
  51.7; Site 2 (186) 0 / 34.4 / 65.6; Site 3 (279) 0.7 / 25.8 / 73.5; Site 4 (398) 0.5 /
  39.4 / 60.1; Overall (1128) 0.4 / 37.3 / 62.3.
- Reading: FAIRD records almost no individuals below 3 mm at any site.

### 2.3 Temporal abundance, FAIRD vs ID

| ID | Status | Note |
| --- | --- | --- |
| R-TA-FI-01 | CURRENT, PROVISIONAL | Parallelism of daily abundance trajectories, Gaussian lme vs Tweedie glmmTMB, conservative and operational runs (entry below) |
| R-DX-FI-ABUND-01 | see Section 3 | Per-site check of the central ID abundance peak |
| R-TA-FI-02 | PLANNED | Rerun with symmetric window exclusion and full `times` levels; then the concordance measures of 2.2 (LOW, HIGH, TOTAL) and internal consistency in abundance |

Common to the R-TA-FI entries:
- Source: `Scherber_temporal_abundance_FAIRD_ID.R`, run of 2026-09-02; output
  `outputs/console/Scherber_temporal_abundance_FAIRD_ID_output.txt` (line numbers
  below); tables `outputs/tables/Scherber_abundance_FAIRD_ID_*.csv`. Session summary of
  the same date (`resumen_sesion_abundancia_FAIRD_ID.txt`, retired from ADP on
  2026-10-12): its results are entered here and in R-DX-FI-ABUND-01, its methodological
  findings are in EXPERIMENT_CONTEXT V1.5, Sections 14 and 15.
- Data: daily abundance (individuals per device-day) from `data_list$etraps_daily`,
  camera systems only (FAIRD1 to FAIRD4, ID1, ID2, ID4), zero days included.
  Conservative: 154 device-days, 22 days. Operational: 139 device-days; the ID battery
  window (Aug 28 to Sep 1, flag `Field_Valid_Operational_ID`) is removed for ID only,
  FAIRD keeps all 22 days (output l.22 to 27).
- Models: BX `nlme::lme(abundance ~ bs(Date.num, 3) * Device_type, random = ~ Date.num |
  Block/Site/Device, correlation = corAR1())`, Type II ANOVA (`car::Anova`); glmmTMB
  Tweedie B0 `(1 | Block/Site/Device)`, B1 = B0 + `ar1(times + 0 | Device)`, B_simple
  `(1 | Site/Device)`, Type III Wald ANOVA on the plotting model (cascade B1, B0,
  B_simple). A glmmTMB fit enters only if `tmb_is_valid()` (positive-definite Hessian,
  finite AIC).
- PROVISIONAL because the operational exclusion is applied to ID only, against the rule
  of October 2026 that the window is excluded for both systems of a comparison
  (EXPERIMENT_CONTEXT V1.5, Section 8.5); the planned rerun (R-TA-FI-02) also completes
  the `times` levels. Not for the manuscript until the rerun.

**R-TA-FI-01, parallelism of FAIRD and ID abundance trajectories** | CURRENT, PROVISIONAL
- Result, conservative (n = 154): BX interaction chi² = 0.89, 3 df, p = 0.828 (Date
  p = 0.029, Device_type p = 0.003) (l.60 to 71). glmmTMB: B0, B1 and B_simple valid;
  B1 best by AIC (1198.4 vs 1215.3 and 1217.3); interaction chi² = 20.13, 3 df,
  p = 0.0002 (Device_type p = 0.24) (l.75 to 112).
- Result, operational (n = 139): BX interaction chi² = 1.65, 3 df, p = 0.647 (Date and
  Device_type p < 0.001) (l.155 to 166). glmmTMB: B1 invalid (non-positive-definite
  Hessian, NA AIC; the 5-day ID gap breaks the AR1 structure); B0 and B_simple valid;
  interaction chi² = 24.25, 3 df, p < 0.0001 (Device_type p = 0.21), cell labelled B0
  (B_simple best by AIC, 1039.7 vs 1041.1) (l.170 to 212, 233 to 239).
- Mean abundance per device-day: FAIRD 12.8 (both runs); ID 35.9 (conservative), 32.4
  (operational) (l.36 to 41, 131 to 136).
- Correction of the session summary: it states that the conservative glmmTMB p came
  from an invalid B1. The output shows B1 valid in the conservative run (l.78, 95) and
  invalid only in the operational run (l.173 to 177, 195).
- Figures (`outputs/figures/`, not versioned): linear y axis capped at 100 (98th
  percentile 90.9 rounded up, `coord_cartesian`, the fit uses all data); 3
  (conservative) and 1 (operational) points above the cap (l.252 to 261). On that axis
  the operational run separates the curves more in the central period than the
  conservative run (session summary; EXPERIMENT_CONTEXT V1.5, Section 14).
- Reading: parallel under the additive Gaussian model in both runs, not parallel under
  the Tweedie log-link model in both runs, with a valid glmmTMB in each; the
  disagreement is not an AR1 artefact.
- Open: reading of the disagreement as a difference of distributional framework is
  pending validation by Scherber (Section 6).

### 2.4 Taxonomic agreement

| ID | Status | Note |
| --- | --- | --- |
| R-V5-H2T-01 | LEGACY V5 | FAIRD vs AMMOD Izsak-Price, 3-day windows, Superfamily: site, intra-habitat, overall (entry below) |
| R-V5-H2T-02 | LEGACY V5, V5 NOT REPRODUCED | Size-threshold sensitivity and sub-3 mm diagnostic (entry below) |
| R-V5-H2T-03 | LEGACY V5, DISCREPANCY V5/V16, V5 NOT REPRODUCED | Taxonomic breadth, FAIRD, ID, AMMOD (entry below) |
| R-V5-H2T-04 | LEGACY V5, DISCREPANCY V5/V16, V5 NOT REPRODUCED | FAIRD identification success by Order (entry below) |
| R-V5-H2T-05 | LEGACY V5, DISCREPANCY V5/V16 | Genus-level daily values quoted in V5 as comparison (entry below) |
| R-V16-H3T-01 | SUPERSEDED by R-V5-H2T-01 and -02 | Daily genus-level Izsak-Price, FAIRD vs AMMOD (entry below) |
| R-V16-H3T-02 | SUPERSEDED by R-V5-H2T-03 | FAIRD vs ID taxonomic comparison (entry below) |
| R-OUT-H2T-01 | OUTPUT ONLY | FAIRD Order composition by site (entry below) |
| R-OUT-H2T-02 | OUTPUT ONLY | Identification success of the minor Orders (entry below) |
| R-OUT-H2T-03 | OUTPUT ONLY | Sub-3 mm diagnostic by site (entry below) |
| R-DC-TX-* | see Section 2.6 | Order composition, three systems, Scherber's multinomial model |

Common to the R-V5-H2T entries (Stage 2 framework, `H2_quali_validation.R`):
- Index: Izsak-Price DeltaS at Superfamily level, taxonomic distances from
  `vegan::taxa2dist`; FAIRD counts reduced to presence/absence; bidirectional
  similarity averaged per site-window. Source: V5, Methods, "H2, Qualitative Agreement".
- Windows: non-overlapping 3-day windows, up to 7 per site (28 site-windows computed);
  windows with fewer than 3 FAIRD individuals resolved to Superfamily excluded, leaving
  24 (output l.172 to 173). FAIRD individuals with Superfamily: 768/1128 = 68.1% (output
  l.162; V5 Methods "~68%").
- PC V16 has no counterpart for the windowed Superfamily analysis: Stage 1 used daily
  genus-level Izsak-Price (H3), dropped in V5 (R-V16-H3T-01).
- Verification: `H2_quali_validation_output.txt` (run of 2026-07-22).

**R-V5-H2T-01, FAIRD vs AMMOD similarity, 3-day windows** | LEGACY V5
- Level: Site (FAIRD vs co-located AMMOD), Intra-Habitat (site-windows pooled by
  habitat), Overall (all 24 site-windows). Metric: community composition,
  presence/absence, Superfamily.
- Source: V5, Results 3.3.2 (Site, Intra-Habitat and Overall Level) and Table 16.
- Result: Site 1 65.2% +/- 2.5 (n = 7 windows); Site 2 61.3 +/- 3.0 (5); Site 3 72.2
  +/- 4.9 (7); Site 4 71.8 +/- 3.2 (5). Maize 63.6 +/- 3.3 (12), Meadow 72.0 +/- 4.1
  (12); Student t(22) = -5.619, p < 0.001; Cohen's d = 2.294; Levene p > 0.05. Overall
  67.8 +/- 5.6 (n = 24), 95% CI [65.5, 70.0]. Superfamilies per window, FAIRD 3.9 /
  2.6 / 7.4 / 8.2, AMMOD 21.1 / 19.0 / 30.6 / 25.8; exact matches 2.9 / 1.8 / 6.0 / 5.0.
- Output: reproduced (l.177 to 217); Levene F = 0.045, p = 0.834 (l.192); d printed
  with negative sign, Maize minus Meadow (l.194). The intra-habitat and overall
  per-window counts of Table 16 are not printed; they agree with the window-weighted
  means of the site rows within rounding. FAIRD individuals per window not printed, not
  verified.
- Reading: FAIRD-AMMOD compositional similarity higher in Meadow than in Maize; all
  sites between 61 and 72%.

**R-V5-H2T-02, size-threshold sensitivity and sub-3 mm diagnostic** | LEGACY V5, V5 NOT
REPRODUCED
- Metric: DeltaS of R-V5-H2T-01 with three AMMOD pools by body length (raw, >= 3 mm,
  >= 5 mm), same FAIRD pool; differences in percentage points (pp).
- Source: V5, Results 3.3.2 (Size-threshold sensitivity analysis) and Table 17.
- Result: AMMOD pool 2023 OTUs, 1783 at >= 3 mm (-11.9%), 771 at >= 5 mm. Overall
  DeltaS 67.8 / 68.0 / 72.3 (+0.2 and +4.3 pp); site differences raw to 3 mm -1.0 to
  +1.1 pp, 3 to 5 mm +3.1 to +5.4 pp (Table 17). Superfamilies present only below 3 mm
  in 26/26 windows, mean 1.9 per window (about 8% of the raw pool); ten Superfamilies
  (Aleyrodoidea, Aphidoidea, Araneoidea, Chalcidoidea, Cucujoidea, Diaprioidea,
  Phoroidea, Platygastroidea, Platypezoidea, Psylloidea).
- PC V16 (H3 block, "Results - Size threshold sensitivity analysis"): same pool sizes;
  56.8% of the >= 3 mm pool lost between 3 and 5 mm. Its DeltaS values are genus-level
  daily, not comparable (R-V16-H3T-01).
- Output: reproduced (l.224 to 277). Table 17 AMMOD Superfamilies per window at >= 3 and
  >= 5 mm only partly printed (l.240).
- V5 NOT REPRODUCED: (a) "a further 56.7% of the >= 3 mm pool" vs 1 - 771/1783 = 56.8%
  (output l.225 to 226; = PC V16). Cause: verified, mis-rounding of 56.76%. (b) "all 26
  analysis windows": the diagnostic uses 26 windows without the FAIRD minimum filter
  (Site 2 has 7 there, 5 in the analysis); the analysed set is 24 (output l.228, 256,
  265). Value reproduced, description wrong.
- Reading: removing the sub-3 mm fraction barely changes DeltaS; restricting AMMOD to
  >= 5 mm raises it by about 4 pp.

**R-V5-H2T-03, taxonomic breadth** | LEGACY V5, DISCREPANCY V5/V16, V5 NOT REPRODUCED
- Metric: distinct Orders, Families and Genera per system over the full study period,
  #N/C excluded; coverage = genera of the e-trap / AMMOD genera.
- Source: V5, Results 3.3.2 (Taxonomic breadth) and Table 18.
- Result: FAIRD 1128 individuals (1102 classified), 10 Orders, 45 Families, 53 Genera;
  ID 2370 individuals, 5 / 8 / 4; AMMOD 15 / 118 / 341. FAIRD coverage: Orders 66.7%
  (10/15), Families 38.1% (45/118), Genera 15.5% (53/341); about 13 times the genera of
  ID. ID coverage "1.5% (5/341)".
- PC V16 (H3 block, "Results - Richness (Table 3)"): FAIRD 11 / 46 / 54, ID 6 / 9 / 5,
  AMMOD identical.
- DISCREPANCY: Cause: verified for V5 (#N/C excluded, output l.136 to 140); hypothesis
  for PC V16 (#N/C counted as one taxon, +1 at every level for both e-traps).
- Output: V5 counts reproduced (l.138 to 140).
- V5 NOT REPRODUCED: (a) ID coverage "1.5% (5/341)" (Table 18 and text) vs 4/341 = 1.2%
  (output l.140); 5 is the PC V16 count. (b) "10 of 15 orders, spanning Diptera,
  Hymenoptera, Coleoptera, Lepidoptera, Hemiptera, and six additional orders" names 11;
  the output lists 10 (l.359 to 368), so five additional. Araneae is one of the 10.
- Reading: FAIRD covers a broader taxonomic range than ID; its coverage of AMMOD
  decreases from Order to Genus.

**R-V5-H2T-04, FAIRD identification success by Order** | LEGACY V5, DISCREPANCY V5/V16,
V5 NOT REPRODUCED
- Metric: % of FAIRD individuals resolved to Genus, to Family only, and to Superfamily
  or lower; denominator 1102 (#N/C excluded).
- Source: V5, Results 3.3.2 (Taxonomic breadth) and Table 19.
- Result: Hymenoptera (n = 44) 81.8 / 18.2 / 0.0; Coleoptera (53) 66.0 / 34.0 / 0.0;
  Diptera (986) 41.5 / 4.1 / 54.5; Lepidoptera (3) 33.3 / 33.3 / 33.3; Overall (1102)
  44.5 / 6.3 / 50.4. Diptera 93.7% of individuals in Maize, 86.7% in Meadow, 89.5%
  overall (986/1102).
- PC V16 (H3 block, "Results - Taxonomic resolution by Order"): same Order rows; Overall
  row n = 1128, 43.4 / 6.3 / 50.4.
- DISCREPANCY: Overall Genus 44.5% (V5, 490/1102) vs 43.4% (PC V16, 490/1128). Cause:
  verified, denominator (output l.392). The same pair explains 43% in V5 Methods vs
  44.5% in V5 Results.
- V5 NOT REPRODUCED: Overall Family 6.3 vs 6.4 (71/1102) and Superfamily or lower 50.4
  vs 49.1 (541/1102) (output l.392); the V5 row sums to 101.2%. Cause: verified in part,
  V5 recomputed Genus over 1102 but carried the other two values from the PC V16 row
  (6.3 = 71/1128); 50.4 matches no count over 1128 either (567/1128 = 50.3), hypothesis
  of a transcription error in PC V16.
- Output: Order rows reproduced (l.376 to 381); Diptera shares l.307 to 308 and l.349.
- Not verifiable: genus share "37 to 51%" across devices (V5 Methods, H2 Qualitative);
  not printed.
- Reading: Diptera dominate FAIRD records and have the lowest genus resolution among the
  main Orders.

**R-V5-H2T-05, genus-level daily values quoted in V5** | LEGACY V5, DISCREPANCY V5/V16
- Source: V5, Results 3.3.2: Site Level (SD "4.8 to 12.0% in the genus-level daily
  analysis"), Intra-Habitat Level (d = 0.507), Overall Level (50.4% +/- 6.0%, 95% CI
  49.0 to 51.8%). Origin: R-V16-H3T-01, dropped in V5.
- PC V16 (H3 block, "Results - Overall concordance"): 50.4 +/- 6.0, CI 49.0 to 51.8
  (n = 68 site-days); per-site SD 4.8 to 6.8%; habitat test t(66) = 2.091, p = 0.040,
  d = 0.507.
- DISCREPANCY: SD range 4.8 to 12.0% (V5) vs 4.8 to 6.8% (PC V16). Cause: not determined.
- Output: no counterpart in the 2026-07-22 outputs; not verifiable.
- Reading: V5 compares the windowed Superfamily results with a superseded analysis;
  these are not current results.

**R-V16-H3T-01, daily genus-level Izsak-Price, FAIRD vs AMMOD** | SUPERSEDED by
R-V5-H2T-01 and -02
- Source: PC V16, "H3: Daily Taxonomic Similarity" (H3 block). Dropped in V5 as a
  hypothesis; some values still quoted (R-V5-H2T-05).
- Data: genus level, daily, 21 days (Sep 4 excluded, AMMOD bottle change), n = 68
  site-days; AMMOD >= 3 mm pool as baseline.
- Result: overall 50.4% +/- 6.0, 95% CI 49.0 to 51.8; Site 1 50.9 +/- 4.8, Site 2 46.8
  +/- 6.8, Site 3 51.9 +/- 5.9, Site 4 52.0 +/- 5.4. Maize 48.9 +/- 6.1 vs Meadow 51.9
  +/- 5.6: Student t(66) = 2.091, p = 0.040, d = 0.507 [0.024, 0.990] (Wilcoxon
  p = 0.102). Linear trend -0.41% per day, t(66) = -4.3, p < 0.001. >= 5 mm threshold:
  53.5% overall (+3.5 pp).
- Replaced by: R-V5-H2T-01 and -02 (Superfamily, 3-day windows, which address the low
  daily FAIRD counts behind the declining trend).

**R-V16-H3T-02, FAIRD vs ID taxonomic comparison** | SUPERSEDED by R-V5-H2T-03
- Source: PC V16, "Additional Analysis: FAIRD vs ID Taxonomic Comparison
  (Supplementary)" (H3 block). Dropped in V5 except for the breadth counts.
- Metric: genus level, 22 days; #N/C counted as a taxon (see R-V5-H2T-03).
- Result: genera FAIRD 54, ID 5, AMMOD 341; exclusive genera FAIRD 20, ID 0, AMMOD 309;
  overlaps FAIRD and AMMOD 31, FAIRD and ID 4, ID and AMMOD 2; genera per g of biomass
  FAIRD 4.49, ID 0.22, AMMOD 4.54; Diptera share 87.4% (FAIRD) vs 96.9% (ID);
  Syrphoidea share of Diptera 2.23% vs 9.28% (4.16-fold).
- Note: the FAIRD Diptera share 87.4% differs from 89.5% in V5 (R-V5-H2T-04), which
  excludes #N/C.
- Replaced by: R-V5-H2T-03.

**R-OUT-H2T-01, FAIRD Order composition by site** | OUTPUT ONLY
- Source: output l.314 to 343.
- Result: Diptera Site 1 91.5% (236/258), Site 2 96.7% (178/184), Site 3 84.8%
  (234/276), Site 4 88.0% (338/384); remaining Orders by site in the same lines.
- Note: site totals here exclude #N/C and differ from those of R-OUT-H2Q-04 (265, 186,
  279, 398).

**R-OUT-H2T-02, identification success of the minor Orders** | OUTPUT ONLY
- Source: output l.379 to 385. Metric as in R-V5-H2T-04.
- Result: Araneae (n = 7) 28.6 / 28.6 / 42.9; Hemiptera (5) 80 / 20 / 0; Dermaptera,
  Neuroptera, Orthoptera (1 each) to Genus; Trichoptera (1) to Family.
- Note: these 16 individuals complete the 1102 of the Overall row of Table 19.

**R-OUT-H2T-03, sub-3 mm diagnostic by site** | OUTPUT ONLY
- Source: output l.255 to 265.
- Result: mean Superfamilies present only below 3 mm per window, Site 1 1.86, Site 2
  1.43, Site 3 2.14, Site 4 2.40 (maximum 3 / 3 / 3 / 4); present in every window at
  every site.

### 2.5 Level differences between systems (diagnostic mixed models)

| ID | Status | Note |
| --- | --- | --- |
| R-V5-H2Q-MM-01 | LEGACY V5 | Diagnostic mixed model A, biomass level by device type, Sites 1 to 3 (entry below; moved from 2.2) |
| R-V5-H2Q-MM-02 | LEGACY V5 | Diagnostic mixed model B, abundance level FAIRD vs ID, Sites 1, 2, 4 (entry below; moved from 2.2) |
| R-V16-LV-01 | SUPERSEDED by R-V5-H2Q-MM-01 | Three-system Kruskal-Wallis and Dunn (entry below) |
| R-V16-LV-02 | SUPERSEDED by R-V5-H2Q-MM-02 | FAIRD vs ID by habitat, Mann-Whitney (entry below) |
| R-V16-LV-03 | SUPERSEDED by R-V5-H2Q-MM-01 and -MM-02 | Habitat x device Kruskal-Wallis (entry below) |

Common to R-V5-H2Q-MM-01 and -02:
- Scope: compare levels between systems (device type main effect and its consistency
  across sites). They do not test temporal concordance; for the temporal question see
  2.1 to 2.3.
- Model: LMM, fixed device type, site and their interaction; repeated measures with
  AR(1), subject = device x site; REML, Kenward-Roger df; SAS PROC MIXED. Source: PC
  V16, "Mixed Model Validation", Methodology. V5 says only "on the log scale".
- Verification (v0.7.1): `SAS/FAIRD_ID_AMMOD-results.html` (repository; program
  `SAS/FAIRD_ID_AMMOD.sas`, inputs `SAS/FAIRD_ID_AMMOD_sas_biomass_12day.csv` and
  `SAS/FAIRD_ID_AMMOD_sas_counts_22day.csv`). REML, Kenward-Roger df. Every V5 value of
  Tables 6 and 7 listed below is printed there (F, df, LS means, contrasts, slices; null
  LRT chi² = 5.88, p = 0.0153 for Analysis B). The back-transformed LS means and the
  17.7-fold ratio are not printed (computed outside SAS, PC V16). An older results file
  (`SAS/old/`, not tracked) is superseded by this one.

**R-V5-H2Q-MM-01, Analysis A, biomass level by device type** | LEGACY V5
- Level: device type x site. Metric: daily biomass (mg), log10 (PC V16). Window: 12-day
  AMMOD overlap (Aug 23 to Sep 3); Sites 1 to 3; n = 92.
- Source: V5, Results 3.1 (Analysis A) and Table 6.
- Result: device type F(2, 27.6) = 32.67, p < 0.001; site F(2, 27.7) = 1.16, p = 0.327;
  device x site F(3, 27.6) = 0.29, p = 0.834. LS means (log10) AMMOD 3.278 +/- 0.120,
  FAIRD 2.030 +/- 0.111, ID non-estimable (no ID at Site 3); contrast AMMOD minus FAIRD
  1.248 +/- 0.163, p < 0.001 (17.7-fold ratio). Slices: Site 1 F(2, 27.4) = 13.30,
  p < 0.001; Site 2 F(2, 27.4) = 11.67, p = 0.0002; Site 3 F(1, 28.3) = 16.93,
  p = 0.0003.
- PC V16: same values; adds back-transformed LS means AMMOD 1898 mg/day (95% CI 1080 to
  3337), FAIRD 107 mg/day (64 to 181). Its breakdown of n = 92 sums to 88; n = 92 is
  the SAS "N Used".
- Reading: AMMOD biomass level above FAIRD at every site, with no detectable
  site-dependence of the gap (see R-OBS-03 on the 17.7 ratio).

**R-V5-H2Q-MM-02, Analysis B, abundance level FAIRD vs ID** | LEGACY V5
- Level: device type x site. Metric: daily abundance, log10 (PC V16). Window: 22 days;
  Sites 1, 2, 4; n = 132.
- Source: V5, Results 3.1 (Analysis B) and Table 7.
- Result: device type F(1, 29.5) = 17.91, p = 0.0002; site F(2, 29.5) = 2.32,
  p = 0.116; device x site F(2, 29.5) = 2.31, p = 0.117. LS means (log10) FAIRD 0.847
  +/- 0.077, ID 1.309 +/- 0.077; contrast FAIRD minus ID -0.462 +/- 0.109, p = 0.0002.
- PC V16: same values; adds back-transformed LS means FAIRD 7.0 individuals/day (4.9 to
  10.1), ID 20.4 (14.1 to 29.4); AR(1) 0.220, null LRT chi² = 5.88, p = 0.015.
- Reading: ID records more individuals per day than FAIRD, with no detectable
  site-dependence of the gap.

**R-V16-LV-01, three-system level comparison, Kruskal-Wallis and Dunn** | SUPERSEDED by
R-V5-H2Q-MM-01
- Source: PC V16, H2, "Section 6B: Three-way device comparison". Dropped in V5.
- Data: Sites 1 to 3 for AMMOD (n = 32) and FAIRD (n = 66), Sites 1 and 2 for ID
  (n = 44); n = 142. Metric: daily captures in each system's own unit (AMMOD mg/day,
  FAIRD and ID individuals/day); the converted-unit table is descriptive only.
- Result: chi² = 77.914, df = 2, p < 0.0001, epsilon² = 0.546. Dunn (Bonferroni):
  AMMOD vs FAIRD Z = 8.826, p < 0.0001; AMMOD vs ID 5.427, p < 0.0001; FAIRD vs ID
  -3.291, p = 0.0015.
- Limitation recorded in PC V16 itself: ranks pooled across groups measured in different
  units; repeated daily observations treated as independent.
- Replaced by: R-V5-H2Q-MM-01 (common unit, AR(1), device x site).

**R-V16-LV-02, FAIRD vs ID by habitat, Mann-Whitney** | SUPERSEDED by R-V5-H2Q-MM-02
- Source: PC V16, "Section 8: E-trap performance by habitat type". Dropped in V5.
- Metric: daily abundance (individuals/day), 22 days.
- Result: Maize (Sites 1 and 2, n = 44 each) FAIRD 10.2 +/- 18.0, ID 36.3 +/- 51.2,
  p < 0.0001; Meadow (Site 4, n = 22 each) FAIRD 18.1 +/- 19.6, ID 35.1 +/- 32.7,
  p = 0.0219.
- Replaced by: R-V5-H2Q-MM-02; abundance trajectories in 2.3 (R-TA-FI).

**R-V16-LV-03, habitat x device, Kruskal-Wallis** | SUPERSEDED by R-V5-H2Q-MM-01 and
-MM-02
- Source: PC V16, "Section 8B: Habitat x device comparison". Dropped in V5.
- Data: five groups (Maize AMMOD, ID, FAIRD; Meadow ID, FAIRD), n = 156; units mixed as
  in R-V16-LV-01.
- Result: chi² = 73.659, df = 4, p < 0.0001, epsilon² = 0.461. Dunn (adjusted): FAIRD
  Maize vs Meadow Z = -2.131, p = 0.165; ID Maize vs Meadow -0.949, p = 1.000; Maize
  FAIRD vs Maize ID -3.768, p = 0.0008; Meadow FAIRD vs Meadow ID -1.641, p = 0.504.
- Note: the Meadow FAIRD vs ID contrast is ns here but p = 0.0219 in R-V16-LV-02
  (different test and adjustment).
- Replaced by: R-V5-H2Q-MM-01 and -MM-02 (site and device x site terms).

### 2.6 Scherber's device-comparison script (main figures proposed by Scherber, 2026-10-08)

| ID | Status | Note |
| --- | --- | --- |
| R-DC-TX-01 | SUPERSEDED by R-DC-TX-02 | Multinomial order composition, original input (AMMOD genus-level records) |
| R-DC-TX-02 | CURRENT, DESCRIPTIVE | Multinomial order composition, order-level presence for all systems (Figure 1 proposed) |
| R-DC-TX-03 | CURRENT | Order detection frequency by system, device-days with the order present |
| R-DC-BM-01 | SUPERSEDED by R-DC-BM-02 | Tweedie biomass model, AMMOD dry mass |
| R-DC-BM-02 | CURRENT | Tweedie biomass model, AMMOD live mass (Figure 2 proposed) |

Common to the R-DC entries:
- Origin: Scherber's email of 2026-10-08 proposes the two figures of this script as the
  main figures of the paper (lean MEE manuscript: no detailed pairwise comparisons; a
  Venn diagram of composition and a device comparison table suggested as additions).
  The model of R-DC-BM is the one he named earlier (email E5,
  `claude/SCHERBER_FEEDBACK_V5_EN.md`).
- Scripts: original `Christoph_device_comparisons_2026-10-08_taxon_and_biomass_analyses.R`;
  runnable version (pasted console output commented, unused packages disabled,
  exploratory `detach()`/`update.packages()` disabled, plots wrapped in `print()`, sink
  log); corrected version (same models, `[INPUT FIX]` tags). Input objects from
  `00_load_data.R` (`etraps_presence_daily`, `ammod_presence_daily`,
  `etraps_daily_overlap`, `ammod_biomass_overlap`).
- Outputs: `Christoph_device_comparisons_output_20261012_143600.txt` (runnable, original
  inputs); `Christoph_device_comparisons_corrected_output_20261012_153222.txt`
  (corrected inputs); figures `Device effects 2026-10-08-2.pdf` (Scherber) and
  `Device effects 2026-10-08-2_corrected.pdf` (2 pages; run before the dual-axis page
  was added).
- Input corrections (all on inputs, models unchanged): (1) AMMOD presence deduplicated
  to Date x Device x Order, the unit of `etraps_presence_daily`; (2) AMMOD biomass =
  `Live_mass` (mg) instead of `Dry_mass_g * 1000` (unit homogeneity with the allometric
  live mass of the e-traps; in the source table Live_mass = 4 x dry mass); (3) top
  orders selected by name (six most frequent, "Indet." excluded) instead of the rank
  index `c(1:5,7)`; (4) optional sensitivity without the ID battery window, off by
  default, not run; (5) additional dual-axis biomass figure (see R-DC-BM-02).
- No random effects in the multinomial model: Scherber's attempts with random effects
  (`mclogit`, `mixcat`, `multpois`) did not run and `brms` could not be installed (his
  email of 2026-10-08). Its p-values ignore Site/Device and the non-independence of
  records and are anti-conservative.

**R-DC-TX-01, multinomial order composition, original input** | SUPERSEDED by R-DC-TX-02
- Model: `nnet::multinom(Order ~ bs(Date.num, 3) * Device_type)`, six orders (Diptera,
  Hymenoptera, Coleoptera, Lepidoptera, Hemiptera, Araneae); `car::Anova` Type II LR.
- Data: AMMOD 2023 records (one per Date x Device x Order x Family x Genus, i.e. genera);
  e-traps 271 records (one per Date x Device x Order).
- Result (runnable output l.32 to 39): Date chi² = 41.04, 15 df, p = 0.0003;
  Device_type 70.63, 10 df, p = 3.4e-11; interaction 72.05, 30 df, p = 2.6e-5.
- Reason superseded: the two record units differ, so the AMMOD panel shows genus
  richness by order and the e-trap panels show order presence; the proportions are not
  comparable. The rank index `c(1:5,7)` excluded "Indet." (rank 6), so the order set is
  unchanged in R-DC-TX-02.

**R-DC-TX-02, multinomial order composition, order-level presence** | CURRENT, DESCRIPTIVE
- Model: as R-DC-TX-01. Data: AMMOD 319 records, e-traps 271, all one per
  Date x Device x Order. Abundance enters for no system (one record per order present
  per device-day).
- Result (corrected output l.67 to 73): Date chi² = 17.16, 15 df, p = 0.31;
  Device_type 86.09, 10 df, p = 3.2e-14; interaction 59.00, 30 df, p = 0.0012.
- Figure (corrected PDF, page 1): AMMOD panel broad (Diptera about 0.2 to 0.3,
  Hymenoptera about 0.25, Coleoptera about 0.2, Lepidoptera up to about 0.25, Hemiptera
  about 0.08); FAIRD and ID panels dominated by Diptera (unchanged from the original,
  same data). FAIRD Coleoptera and Araneae bulge near Sep 7 rests on few records (5
  Araneae), spline artefact.
- Reading: the systems differ in which orders they detect; no common temporal trend in
  composition once the record unit is harmonised. Descriptive only (p-values
  anti-conservative).

**R-DC-TX-03, order detection frequency by system** | CURRENT
- Source: corrected output l.18 to 36 (order-level presence records, i.e. device-days
  with the order present).
- Result (AMMOD / FAIRD / ID): Diptera 76 / 79 / 61; Hymenoptera 67 / 24 / 14;
  Coleoptera 51 / 30 / 11; Lepidoptera 43 / 3 / 2; Hemiptera 34 / 5 / 3; Araneae 8 / 5 /
  0; Opiliones 10 / 0 / 0; Psocoptera 10 / 0 / 0; Trichoptera 6 / 1 / 0; Orthoptera 5 /
  1 / 0; Dermaptera, Isopoda, Neuroptera, Odonata 2 / 0 to 1 / 0; Trombidiformes 1 / 0 /
  0; Indet. 0 / 13 / 17.
- Denominators (v0.7.3): valid taxonomic unit-days from the `effort_base` rules
  (`00_load_data.R`, `Taxo_Valid`): AMMOD 76 (AMMOD1 to AMMOD3, 21 each; AMMOD4 13),
  FAIRD 88, ID 66 (ID3 excluded; outage days are valid unit-days). Derived from the
  rules, not printed by a run; equal to EXPERIMENT_CONTEXT V1.5, Section 9.
- Reading: Diptera detected on almost every device-day by all systems; Lepidoptera and
  Hemiptera detected by AMMOD on many device-days and by the e-traps on very few.
  Causes: hypotheses only (R-OBS-06).

**R-DC-BM-01, Tweedie biomass model, AMMOD dry mass** | SUPERSEDED by R-DC-BM-02
- Model: `glmmTMB(biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site/Device),
  family = tweedie)`, log link; no Block, no AR1. Data: 116 rows, 12 dates (Aug 23 to
  Sep 3), 4 sites, 10 devices. AMMOD = `Dry_mass_g * 1000`.
- Result (runnable output l.48 to 80): intercept 6.85; Device_type FAIRD -2.04
  (p = 0.0002), ID -0.76 (p = 0.148); spline x FAIRD 2.33 / 2.01 / 0.43, spline x ID
  -1.21 / 2.81 (p = 0.015) / 0.81; Device:Site variance 0.046, Site 1.7e-9; dispersion
  12.6; AIC 1577.0.
- Reason superseded: dry mass for AMMOD against allometric live mass for the e-traps.
  The non-significant ID level contrast was an artefact of the units (R-OBS-07).

**R-DC-BM-02, Tweedie biomass model, AMMOD live mass** | CURRENT
- Model and data: as R-DC-BM-01, AMMOD = `Live_mass` (mg).
- Result (corrected output l.87 to 119): intercept 8.24; Device_type FAIRD -3.42
  (p = 7.6e-11), ID -2.14 (p = 2.5e-5); spline x FAIRD 2.30 / 2.04 / 0.40, spline x ID
  -1.22 / 2.83 (p = 0.007) / 0.79; Device:Site variance 0.038, Site 2.6e-9 (Site
  intercept collapses); dispersion 12.1; AIC 1657.4 (not comparable with R-DC-BM-01,
  response rescaled).
- Check: intercept and Device_type shifts = 1.386 = ln 4; interaction coefficients
  unchanged within 0.03. Level contrasts refer to the reference date (Aug 23), because
  of the interaction.
- No test of the Date x Device_type interaction printed by this script (Wald
  coefficients only); the three-system test with Block and AR1 is R-TB-ALL-01.
- Figures: page 2, Scherber's `effects` plot on the response scale (e-trap curves
  compressed near zero by the live-mass AMMOD level); page 3 (added after the run),
  dual-axis figure: linear axes, FAIRD/ID left, AMMOD right divided by a data-derived
  factor k (sum AMMOD / sum FAIRD over co-located site-days, 2 significant digits).
  Juan's decision: figure kept as designed, k data-derived (not rounded to 10).
  k = 11 over 32 co-located site-days (corrected output of 2026-10-12 16:02:55, l.122;
  the 15:32 output predates the dual-axis page and is otherwise identical).
- Exploratory checks in the script: `fitdistrplus::fitdist` gamma fails (zero-biomass
  days); `tweedie_profile` result printed as NULL (`t1$p.max` not in the object;
  Scherber's note: power about 1.6).
- Reading: same temporal shape as R-DC-BM-01; with live mass both e-traps record
  significantly less biomass than AMMOD at the reference date.

---

## 3. DIAGNOSTICS AND DATA FACTS

**R-DX-AF-RATIO-01, AMMOD/FAIRD ratio drift** | CURRENT
- Source: Section 2e, ratio of device means per block-day (PAIRED_BLOCK).
- Result: Spearman rho(ratio, date) = -0.08; AMMOD-scale factor for figures about 10.
- Reading: no monotonic drift; divergence is of shape, concentrated at the ends of the
  12-day window.

**R-DX-AI-RATIO-01 and R-DX-FI-RATIO-01, ratio drift** | CURRENT
- Result: AMMOD/ID rho = -0.77; FAIRD/ID rho = -0.57.
- Reading: FAIRD loses relative signal in the second half of the season (consistent with
  context Section 13).

**R-DX-INTRA-01, intra-FAIRD check of the PAIRED_BLOCK assumption** | CURRENT
- Source: `diagnostic_per_site_biomass.R`, Table 3c. Diagnostic only, not for the
  manuscript.
- Result: FAIRD3/FAIRD4 on AMMOD3 days 0.34 to 1.16 (median 0.87); FAIRD1/FAIRD2 0 to 6.1.
- Reading: Meadow sites share their dynamics (supports AMMOD3 vs ID4 pairing); Maize
  instability is irrelevant because Maize pairs are co-located.

**R-DX-ID-WIN-01, structure of the ID battery window** | CURRENT
- Result: 15 device-days, 5 zeros and 10 with biomass; by session calendar 3 fully off,
  8 partial, 4 full effort; Sep 2 partial for ID2 and ID4. ID1 Aug 29 is a genuine zero
  under full effort. Details: context Section 8.5.

**R-DX-ID-BOUT-01, temporal clustering of ID detections** | CURRENT (not used in models)
- Result: bouts/detections (same class, gap < 30 s) 0.37 and 0.29 for ID2 Aug 30 and 31,
  0.23 for ID4 Aug 24, against device medians about 0.6; a session-level Poisson model
  predicts much of it.
- Reading: compatible with revisits, not separable from dense genuine activity.

**R-DX-WX-01, weather during the campaign** | CURRENT
- Result: cold front Aug 27 to Sep 1 (daytime mean 14 C on Aug 28 and 29, radiation
  4 to 5 MJ m-2 on Aug 27 to 29, rain most days); battery failures lag the radiation
  minimum by about one day. Details: context Section 8.7.

**R-DX-AR1-01, behaviour of the AR1 term in the pair models** | CURRENT
- Result: rho 0.02 to 0.21, SD 0.38 to 0.80 (log scale) in all valid pair models; Site
  and Device:Site variances about 0.
- Reading: AR1 acts as a device-day effect; its BLUPs are the day-to-day deviations
  (excluded when computing the HIGH component).

**R-DX-COVER-01, daily recording coverage** | CURRENT (data fact)
- FAIRD 07-19 h continuous; ID 9 sessions of about 40 min within 07-09, 12-14 and
  17-19 h; AMMOD 00-24 h. Shared calendar-day boundary, lag 0 by design.

**R-DX-FAIL-01, attribution of ID failures** | CURRENT (data fact, for the reliability
analysis)
- ID3 total loss: corrupt or badly flashed SD card, noticed at the end of the experiment
  (operator). ID power failures: two devices on one solar panel and one battery in the
  maize field (deployment configuration, corrected mid-experiment). Not device-intrinsic.
- AMMOD losses (AMMOD4 no biomass; AMMOD3 from Aug 30, ATL processing and jam): to be
  classified with the same criterion.

**R-DX-PRES-01, record units of the presence objects in `00_load_data.R`** | CURRENT
(data fact, v0.6)
- `etraps_presence_daily`: one row per Date x Device x Order (`present = 1`), abundance
  collapsed; 271 rows. `ammod_presence_daily`: from `ammod_taxonomy_presence`,
  deduplicated by Date x Device x Order x Family x Genus; 2023 rows (319 at order
  level).
- Consequence: any cross-system comparison of these two objects must first bring them
  to a common unit (done inside the corrected script, not in `00_load_data.R`, which
  other scripts use at genus level).

**R-DX-FI-ABUND-01, per-site check of the central ID abundance peak** | CURRENT
(diagnostic, v0.7), NOT FOR MAIN TEXT (C2)
- Source: `diagnostic_per_device_operational.R` (throwaway diagnostic, no output file in
  the repository); values from the session summary of 2026-09-02. Operational period
  (ID battery window removed), sites with both camera systems (Site1, Site2, Site4),
  central window Aug 27 to Sep 1 as defined in the script.
- Result: Site2, marked central ID peak (central mean ID 20 vs FAIRD 12.5; ID2 maximum
  91); Site4, ID starts high and declines (ID4 maximum 141, the highest daily value of
  the dataset; central mean FAIRD 30.2 vs ID 19); Site1, nearly parallel.
- Reading: the aggregated central ID peak of R-TA-FI-01 is the superposition of three
  different site dynamics, not a consistent system pattern; with one device per system
  per site it cannot be separated from device idiosyncrasy.

---


## 4. OBSERVATIONS (cross-entry syntheses)

Format and update rule: see Conventions. Observations condense reasoning across entries;
they are not ecological interpretation.

**R-OBS-01, pooled p-values in PC V16 come from another model**
- Rests on: R-V5-H2Q-02, R-V5-H2Q-03, R-V16-H2Q-01.
- Statement: for the Maize (n = 24) and Overall (n = 32) levels, PC V16 reports as the
  regression p the p of log_ammod in the model with a Site covariate (0.006, 0.0036),
  not the p of the simple OLS whose R² it reports (0.0047, 0.0007).
- Basis: verified (output l.201, 224, 317, 318).
- Implication: no pooled-level p is taken from PC V16; V5 values are the correct ones.

**R-OBS-02, the V5 Durbin-Watson limitation is probably an ordering artefact**
- Rests on: R-V5-H2Q-02, R-V5-H2Q-03, R-V16-H2Q-02.
- Statement: on a stacked multi-site series DW depends on row order. Sorted by site then
  date, consecutive residuals are consecutive days of one site, which is the temporal
  autocorrelation DW is meant to test: 1.845 (Maize) and 1.864 (Overall), with
  per-site values 1.43 to 1.99. Sorted by date then site, consecutive residuals are
  different sites on the same day, so the low DW (0.917, 1.385) measures the shared
  daily signal between sites, not a violated assumption.
- Basis: verified; the same models sorted by date then site give exactly the V5 and PC
  V16 values (diag 8A: 0.917, p = 0.0002; 1.38, p = 0.0409).
- Implication: the V5 limitation sentence in 3.3.1.3 rests on a doubtful value; V6
  drops or recomputes it (Section 5).

**R-OBS-03, the 17.7-fold AMMOD/FAIRD ratio is not a conversion factor**
- Rests on: R-V5-H2Q-MM-01, R-TB-AF-01, R-DX-AF-RATIO-01.
- Statement: 17.7 is the ratio of geometric means (LS-mean difference 1.248 on log10)
  averaged over the 12-day window. R-TB-AF-01 rejects a constant ratio over time, so it
  does not convert FAIRD biomass to AMMOD biomass for a given day. The figure factor of
  about 10 (R-DX-AF-RATIO-01) is a different estimand (ratio of device means per
  block-day) and is not interchangeable with 17.7.
- Basis: verified for the definitions; the size of the gap between 17.7 and 10 is not
  explained here.
- Implication: if the manuscript cites both numbers, it states what each one is.

**R-OBS-04, Stage 2 and Stage 3 agree on what AMMOD and FAIRD share**
- Rests on: R-V5-H2Q-01, R-V5-H2Q-04, R-CC-AF-01, R-TB-AF-01.
- Statement: both frameworks show shared day-to-day fluctuation, mainly in Maize (CCF
  lag 0 0.781 at Site 2 and 0.657 for the Maize summed series; HIGH 0.87 Maize only).
  Neither shows a shared seasonal shape (LOW not estimable over 12 days; non-parallel on
  the log link). Part of the daily synchrony is the common response to the cold spell.
- Basis: verified (values of the cited entries).
- Implication: "FAIRD tracks AMMOD" is qualified in the manuscript as day-to-day
  synchrony, not proportional tracking. (v0.6: the main-text value is the two-block HIGH
  0.61, R-OBS-05; Site 2 CCF is LEGACY and device-level.)

**R-OBS-05, level of comparison for the manuscript (Scherber C2)**
- Rests on: R-CC-AF-01, R-CC-AI-01, R-CC-FI-01, R-CC-INT-01, R-CC-DIS-01, R-TB-INT-01,
  R-TB-AF-01, `claude/SCHERBER_FEEDBACK_V5_EN.md` (C1, C2).
- Statement: main-text comparisons are between systems over both blocks. Admissible:
  pair models and concordance at two-block level (e.g. AMMOD vs FAIRD HIGH 0.61
  [0.05, 0.88], parallelism LRT 0.013, Holm 0.039); Maize-only rows as single-block
  secondary values; FAIRD halves over both blocks as the internal reference. Not
  admissible in the main text: Maize device pairs (AMMOD1 vs AMMOD2, FAIRD1 vs FAIRD2,
  ID1 vs ID2) and everything derived from them (disattenuation), and all LEGACY
  site-level device pairs.
- Basis: verified against Scherber's written criticism (email E4) and Juan's decision
  of 2026-10-12.
- Implication: the device-pair references stay as diagnostics or appendix; the
  manuscript does not use them to calibrate between-system agreement.

**R-OBS-06, Figure 1 shows detection frequency, not abundance composition**
- Rests on: R-DC-TX-02, R-DC-TX-03, R-DX-PRES-01, R-V5-H2T-04, R-DX-COVER-01,
  R-OUT-H2Q-04.
- Statement: in R-DC-TX-02 each system contributes one record per order present per
  device-day, so the Diptera dominance of the e-trap panels reflects that Diptera are
  present on almost every device-day while other orders are rare, not the 89.5% Diptera
  share of FAIRD individuals. The low e-trap detection of Lepidoptera and Hemiptera has
  candidate causes that these data cannot separate: recording coverage (AMMOD 24 h,
  FAIRD 07-19 h, ID three 2 h windows; nocturnal and crepuscular Lepidoptera, Juan's
  hypothesis); FAIRD size threshold of about 3 mm (many Hemiptera are smaller);
  classification to "Indet." (30 e-trap records, none in AMMOD); metabarcoding
  sensitivity to single small individuals and fragments; ID attractant platform
  selecting flower visitors.
- Basis: verified for the record unit; hypothesis for the causes. Partial check
  possible: hour of `Checkin` of the FAIRD Lepidoptera and Hemiptera records (Section 6).
- Implication: the figure and caption speak of detection frequency; causes go to the
  Discussion as candidates.

**R-OBS-07, the biomass unit correction changes levels, not shapes**
- Rests on: R-DC-BM-01, R-DC-BM-02.
- Statement: replacing AMMOD dry mass by live mass multiplies AMMOD by 4; on the log
  link this shifts the intercept and the Device_type contrasts by ln 4 and leaves the
  interaction coefficients unchanged. The level of ID relative to AMMOD goes from not
  significant (p = 0.148) to significant (p = 2.5e-5).
- Basis: verified (coefficients of both runs).
- Implication: Figure 2 and any level statement use live mass; the original figure
  suggested comparable AMMOD and ID levels by an artefact of units.

---


## 5. MANUSCRIPT CORRECTIONS FOR V6

All V5 NOT REPRODUCED items. Correct value = output of the script that produced the V5
value (or the count it implies). Details and causes in the cited entries.

| Entry | V5 location | V5 value | Correct value | Cause |
| --- | --- | --- | --- | --- |
| R-V5-H1Q-FAIRD-02 | Table 10, FAIRD Maize vs Meadow | low-weight days 3/22 | 1/22 (weight < 0.5) | hypothesis, other weight threshold |
| R-V5-H1Q-AMMOD-02 | Table 10, AMMOD1 and AMMOD2 vs AMMOD3 | low-weight days 2 and 1 | 0/8 and 0/8 | hypothesis, as above |
| R-V5-H1Q-ID-01 | 3.2.1.1 and Table 8, ID conservative | p = 0.461 | p = 0.463 | origin of 0.461 unknown |
| R-V5-H1Q-ROB-01 | Table 9, ID conservative row | pseudo-R² -0.013, delta 0.037 | no documented log-scale value; recompute or remove the row (original scale 0.013) | hypothesis, original-scale value entered with a sign change |
| R-V5-H2Q-02 | 3.3.1.2, Maize pooled | Kendall p = 0.011 | p = 0.019 | not determined (six variants excluded) |
| R-V5-H2Q-02 | 3.3.1.2, Maize pooled | bootstrap CI [0.033, 0.554] | [0.041, 0.556] | hypothesis, Stage 1 seed or scheme |
| R-V5-H2Q-03 | 3.3.1.3, Overall pooled | bootstrap CI [0.063, 0.543] | [0.064, 0.541] | hypothesis, as above |
| R-V5-H2Q-03 | 3.3.1.3, Overall pooled, limitation sentence | DW = 0.917, p < 0.001 | DW = 1.864, p = 0.306; drop the limitation (R-OBS-02) | verified, row order (date then site) |
| R-V5-H2T-02 | 3.3.2, size-threshold text | "a further 56.7%" | 56.8% | verified, mis-rounding |
| R-V5-H2T-02 | 3.3.2, sub-3 mm diagnostic text | "all 26 analysis windows" | 26 diagnostic windows (no FAIRD filter); analysed set 24 | verified |
| R-V5-H2T-03 | Table 18 and text | ID coverage 1.5% (5/341) | 1.2% (4/341) | verified, PC V16 count carried over |
| R-V5-H2T-03 | 3.3.2, Taxonomic breadth text | "and six additional orders" | five additional orders | verified, 10 Orders listed in output |
| R-V5-H2T-04 | Table 19, Overall row | Family 6.3%; Superfamily or lower 50.4% | 6.4% (71/1102); 49.1% (541/1102) | verified in part, PC V16 row carried over |

Notes:
- Table 13 significance stars are unchanged by the R-V5-H2Q-02 corrections (p < 0.05 in
  both values).
- Not V6 corrections (no V5 value differs from the output): DW Maize in PC V16
  (R-V5-H2Q-02); the H2 quantitative DeltaR² rounding (R-V5-H2Q-05, V5 correct).

---

## 6. VALUES PENDING RESOLUTION

Values whose cause is not determined, that no output in ADP prints, or whose hypothesis
needs a source not in ADP. "Testable in ADP" means a rerun or an extra print in an ADP
script would settle it.

| Entry | Value | Problem | Needed to resolve |
| --- | --- | --- | --- |
| R-V5-H1Q-FAIRD-02, R-V5-H1Q-AMMOD-02 | low-weight days 3/22; 2 and 1 (V5) | cause hypothesis (weight threshold) | threshold used for V5 Table 10; testable in ADP by rerunning with other thresholds |
| R-V5-H1Q-ID-01 | conservative p 0.461 (V5) | origin unknown | Stage 1 `H1_consistency_ID.R` output |
| R-V5-H1Q-ROB-01 | ID operational log pseudo-R² 0.433; ID conservative -0.013 | not printed by any ADP output; -0.013 has no documented source | log-scale `rlm` on the ID pairs (extra print in `H1_consistency_summary.R`) |
| R-V5-H2Q-02 | Kendall p 0.011 (V5) | cause not determined; six variants excluded (diag 8B) | Stage 1 `H2_validation_aggregated_levels.R`; or not worth resolving (V6 uses 0.019) |
| R-V5-H2Q-02, -03 | bootstrap CI (V5 = PC V16) | cause hypothesis (seed) | Stage 1 script, or not worth resolving (V6 uses the output values) |
| R-V5-H2Q-04 | CCF lag 0 and lag +1 in PC V16 | cause not determined; log10 hypothesis rejected (diag 8C) | Stage 1 `H2_validation_*.R`; not needed for V6 (V5 values reproduced) |
| R-V5-H2Q-06, R-V16-H2Q-03 | body length and DeltaS tests | not verifiable | Stage 1 `BodySize_Distribution.R`, `TaxConsistency_IntraCrossHabitat.R` |
| R-V5-H2T-01 | FAIRD individuals per window (Table 16) | not printed | testable in ADP: extra print in `H2_quali_validation.R` |
| R-V5-H2T-04 | Superfamily or lower 50.4% in PC V16 | transcription hypothesis | not needed for V6 (output value used) |
| R-V5-H2T-04 | genus share 37 to 51% across devices (V5 Methods) | not printed | testable in ADP: genus share by device |
| R-V5-H2T-05 | SD range 4.8 to 12.0% (V5) | cause not determined (PC V16 4.8 to 6.8%) | Stage 1 H3 output; or remove the comparison in V6 |
| R-CC-INT-01 | ID1 vs ID2 B_22d full LOW 0.991, interval [-0.558, 0.982] | estimate above the upper limit; interval differs from the 2026-10-05 log ([-0.60, 0.988]) although the seed is fixed | change the LOW interval method for the paper, or report LOW without interval |
| R-OBS-06 | hour of `Checkin` of FAIRD Lepidoptera and Hemiptera records | not computed | `table()` of the hour by Order in `faird_data` |
| R-DC-BM-02 | Tweedie power from `tweedie_profile` | printed as NULL | print the correct element of the profile object |
| R-TA-FI-01 | interaction p BX 0.828 / 0.647 vs Tweedie 0.0002 / < 0.0001 | operational exclusion asymmetric (ID only); reading of the BX vs Tweedie disagreement (distributional framework) not validated | rerun R-TA-FI-02 with symmetric exclusion; Scherber's reading of the disagreement |
| R-DX-FI-ABUND-01 | per-site central means and maxima | values only from the session summary; no output file | rerun `diagnostic_per_device_operational.R` with console capture, if the values are needed |
