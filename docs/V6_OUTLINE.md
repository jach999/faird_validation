# FAIR-DEVICE VALIDATION STUDY, OUTLINE OF MANUSCRIPT V6

**Version:** 0.2 (2026-10-29, close of session 2026-10-27/29). Decisions D1, D2, D4 to D7
taken; all analyses of Section 7 of v0.1 done or dropped; figure numbers of the V6 draft;
supplementary sections and script map (Sections 5 and 8). v0.1: agreed in session
2026-10-13/14, written down 2026-10-27.
**Purpose:** structure of the V6 rewrite and one sentence per claim, to be used as
paragraph headings when drafting. Each claim cites the registry entry that supports it
(`docs/RESULTS_REGISTRY.md`); the registry, not this file, holds the values and their
status. Scherber's guidelines (G1 to G18, emails E1 to E7) are in
`claude/20261013_SCHERBER_GUIDELINES_EN.md`.
**Draft:** Claude Doc "FAIRD Paper V6" (claude.ai artifact
`Rh8dWfkNim5wZrbBsbaQ82`), written 2026-10-27/28 from this outline. Juan merges the V5
Introduction and Methods into it by hand (2026-10-29).

Status marks used below: **[confirmed]** agreed by Juan; **[proposed]** proposed by
Claude, not yet confirmed; **[pending analysis]** needs a run before the text can be
written; **[open]** decision not taken; **[dropped]** decided against.

---

## 1. Principles of the rewrite

- No hypothesis-driven structure (no H1/H2/H3) and no V5 level hierarchy (Site,
  Intra-Habitat, Overall) as section titles. The V5 terms are too rigid for the
  analyses now planned, which take a general view rather than an exhaustive test of
  consistency. [confirmed]
- Results follow Scherber's two main figures (G15): one section on quantity, one on
  composition. [confirmed]
- System-level comparisons only in the main text (G6, G8); unit-pair comparisons go to
  the supplementary material, except the one-sentence statement of AMMOD and ID
  internal consistency (D1). [confirmed]
- Figures on linear axes, AMMOD on a secondary axis with a data-derived factor (G14).
  [confirmed]
- Only analyses already done; no new tests (Juan, 2026-10-27: a compact paper that tells
  one story and follows Scherber). [confirmed]

### Vocabulary

Defined once in Methods and used consistently. [confirmed]

| Term | Meaning | Registry / script name |
| --- | --- | --- |
| seasonal trajectory | shape of a system's curve over the period | LOW |
| day-to-day fluctuation | daily deviations around that curve | HIGH |
| proportional | two systems whose curves keep a constant ratio (parallel on the log link) | parallelism test |
| across both habitats / in Maize | Maize and Meadow together / Maize only | two-block / Maize-only rows |
| system / unit | FAIRD, AMMOD or ID / one device of a system | Device_type / Device |

---

## 2. Materials and Methods

Figures 1 (the three systems in the field) and 2 (study site and layout of the four
sites): photographs and map, to be supplied by Juan.

- **2.1 Devices.** FAIRD, AMMOD, ID; power supply. [placeholder in the draft: power-supply
  contradiction to resolve]
- **2.2 Study site and field design.** Four sites, two habitats, 12 units, 22 days; one
  unit per system and site, stated as the replication limit (G10).
- **2.3 Data processing.** Per system; FAIRD and ID biomass estimated from counts and
  allometry (Sohlström et al. 2018), AMMOD live mass; daily recording coverage (AMMOD 00-24 h, FAIRD
  07-19 h, ID three 2 h windows) [R-DX-COVER-01]; abundance exists only for the camera
  systems (AMMOD taxonomy is presence/absence). [placeholders: metabarcoding protocol,
  source of AMMOD body length]
- **2.4 Data availability and exclusions.** The whole failure story in one place
  [confirmed]: Table 1, data losses by system with cause and treatment, unit-days of a
  planned 88 per system (`2-4_data_availability.R`); the ID power outages coincided with
  the cold spell of 28 August to 1 September [R-DX-WX-01, R-DX-ID-WIN-01, R-DX-FAIL-01];
  exclusion rule for the ID outage window (full series and window excluded, symmetric for
  both systems of a comparison); effective sample sizes. The Discussion refers to the
  table and does not repeat the facts. [done]
- **2.5 Statistical analysis.** Tweedie GLMM: Figure 3 from Scherber's model, tests from
  the same model with Block and AR1 (D2); three-system Wald test and pair LRTs with Holm
  over the three 12-day pairs; FAIRD vs ID also over 22 days. The two concordance
  components. Spearman rho of daily system means as a descriptive measure (D7). Multinomial
  model of order composition (no unit random effect, G13, stated with the replication
  limit). Izsak-Price similarity in 3-day windows at Superfamily level. Sensitivities in
  the supplementary material (S1).

---

## 3. Results

Spearman rho of daily system means in parentheses next to every model-based comparison
[R-CC-DESC-01; D7]. [confirmed]

### 3.1 Quantity

Figure 3: biomass of the three systems over time, 12 days with AMMOD biomass (Scherber's
model, corrected inputs, dual axis, k = 11; `3-1a_quantity_trajectories.R`). [confirmed]

1. *The two halves of the FAIRD network followed proportional trajectories and shared
   their day-to-day fluctuation.* [R-TB-INT-01, R-CC-INT-01, R-CC-DESC-01; FAIRD halves
   across both habitats] Opening paragraph, as the internal reference. One sentence on
   AMMOD and ID, whose internal consistency exists only as the Maize unit pairs, details
   in S2 (D1). [confirmed]
2. *AMMOD recorded substantially higher daily biomass than both camera systems.*
   [R-DC-BM-02, R-OBS-07; contrast at the reference date, see R-OBS-03 before quoting a
   single ratio]
3. *FAIRD biomass was not proportional to AMMOD biomass over time.* [R-TB-AF-01,
   R-TB-ALL-01, R-TB-ALL-02]
4. *FAIRD and AMMOD shared the day-to-day fluctuation of biomass, partly through a
   common response to the cold spell.* [R-CC-AF-01, R-CC-SENS-01, R-OBS-04]
5. *ID showed a similar day-to-day agreement with AMMOD, but whether its biomass was
   proportional to AMMOD could not be determined.* [R-CC-AI-01, R-TB-AI-01]
6. *FAIRD and ID shared the timing of the seasonal peak, but diverged in the second half
   of the period.* [R-TB-FI-22F, R-TB-FI-22X, R-CC-FI-22-01, R-CC-SENS-01]
7. Closing paragraph on abundance, FAIRD vs ID only, no figure in the main text (S4).
   [confirmed]
   - *ID recorded more individuals per unit-day than FAIRD.* [R-TA-FI-01, S4 means:
     35.9 vs 12.8 over 22 days, 32.4 vs 12.6 with the window excluded for both]
   - *The two camera systems agreed as little in daily abundance as in biomass, so their
     divergence does not come from the mass estimates.* [R-CC-DESC-01] The trajectory
     tests of abundance (R-TA-FI-01, provisional) are not used. [rerun R-TA-FI-02
     dropped, no new tests]

### 3.2 Composition

Figure 4: order composition over time, one panel per system (Scherber's multinomial model,
corrected inputs; `3-2a_composition_orders.R`). [confirmed]

1. *The three systems differed in which orders they detected, with no common temporal
   trend in composition.* [R-DC-TX-02; descriptive, p-values anti-conservative]
2. *Diptera were detected on almost every unit-day by all systems; Lepidoptera and
   Hemiptera mainly by AMMOD.* [R-DC-TX-03, R-OBS-06] Table 2, unit-days per order out of
   the valid taxonomic unit-days (AMMOD 76, FAIRD 88, ID 66). [done]
3. *FAIRD and AMMOD communities overlapped by about two thirds at Superfamily level, more
   in Meadow than in Maize.* [R-V5-H2T-01, CURRENT since D4; `3-2b_composition_overlap.R`]
4. *Excluding AMMOD taxa smaller than 3 mm barely changed this overlap; FAIRD recorded
   almost no individuals below that size.* [R-V5-H2T-02, R-OUT-H2Q-04; the FAIRD
   body-length shares are printed by no current script, a print is to be added to 3-2b]
5. *FAIRD covered a broader taxonomic range than ID, consistent with their different
   classification approaches, but a small fraction of the genera detected by AMMOD.*
   Table 3. [R-V5-H2T-03; wording on ID descriptive of its design, not a deficiency]
6. *Diptera made up about 90% of FAIRD records and had the lowest identification success
   to genus.* [R-V5-H2T-04]
7. Venn (Euler) diagram of composition (G9). [dropped, D6]

---

## 4. Discussion

1. **What a camera system captures of the Malaise signal.** FAIRD follows the
   day-to-day fluctuation, not the amount; a complement to Malaise sampling, not a
   substitute. Qualify E6 point 6 ("high concordance"), which the data do not support as
   stated.
2. **Limits of taxonomic coverage.** Poorly resolved Diptera, the 3 mm threshold,
   recording coverage; candidate causes for the low detection of Lepidoptera and
   Hemiptera [R-OBS-06].
3. **Device reliability in context.** Refers to Table 1; causes, silent failures (ID3),
   non-random data loss; trade-off between field processing and data transmission.
   Juan's position, not a ranking of systems.
4. **Choosing a system.** Table 4, comparison of the three systems as deployed (G4;
   placeholder: time to data), and Figure 5, decision tree (G3; placeholder).
5. **Limitations.** One unit per system and site; AMMOD biomass over 12 days only;
   asymptotic tests.
6. **Outlook.** Light-weight devices combining optical methods and eDNA (E6 point 7).

---

## 5. Supplementary material

| Section | Content | Scripts |
| --- | --- | --- |
| S1 | Pair models in detail: pairing criteria, Wald and LRT, outage window excluded; Block random and Gaussian model as checks of the three-system test; identity-link tests; concordance with a spline of 4 df | `supplementary/S1_pair_model_sensitivity.R`; df 4 concordance in `3-1b_quantity_concordance.R` |
| S2 | Internal consistency and simple correlations: Maize unit pairs of AMMOD and ID (proportionality, both concordance components); Pearson and Spearman of daily means for every comparison; conventional log-scale correlations between paired units | `3-1a` (internal parallelism), `3-1b` (concordance, Section 9 correlations), `supplementary/S2_unit_pair_correlations.R` |
| S3 | Per-site diagnostics, biomass and abundance; ID detection classes on outage days | `supplementary/S3_site_diagnostics.R` |
| S4 | Abundance of the camera systems, descriptive figure (caption: the highest ID2 values fall on outage days) | `supplementary/S4_abundance_cameras.R` |
| S5 | Taxonomic details: Izsak-Price by site and habitat, size threshold by site, FAIRD identification success by order | `3-2b_composition_overlap.R` |
| S6 | Code and data: repository `jach999/faird_validation` (Zenodo DOI at submission) | all |

---

## 6. Decisions

| ID | Decision | Outcome | Status |
| --- | --- | --- | --- |
| D1 | Internal consistency of AMMOD and ID exists only as Maize unit pairs (G8) | Revised by Juan (2026-10-27): one sentence in the main text, details in S2. Reason: whether FAIRD resembles AMMOD or ID loses weight if the reader does not see whether those systems are consistent within themselves. AMMOD3 and ID4 not added as extra pairs (no new tests) | [confirmed] |
| D2 | Figure and test from the same model | Figure 3 from Scherber's model; tests from the same model with Block and AR1 (3-1a) | [confirmed] |
| D3 | Place of abundance | Closing paragraph of 3.1 | [confirmed] |
| D4 | Taxonomic entries used in V6 were LEGACY or OUTPUT ONLY | Reclassified as CURRENT in registry v0.7.8 | [confirmed] |
| D5 | Izsak-Price at system level | Not recomputed (no new tests); the unit-level comparison goes with the replication statement | [dropped] |
| D6 | Venn diagram (G9) | Left out (Juan: unclear what Scherber meant) | [dropped] |
| D7 | Simple correlations | Spearman of daily system means in parentheses in the main text, Pearson and both in S2. Reasons (Juan): a reviewer expects a simple correlation and would suspect cherry-picking without it; Scherber objects to figure slots and sections on simple correlations, not to a mention; Spearman is what most precedent studies used and is not driven by the peak days | [confirmed] |

## 7. Analyses needed before drafting

All items of v0.1 done (k, interaction test, Table 1 counts, Table 2 denominators) or
dropped (abundance rerun R-TA-FI-02). Remaining, not needed for the draft:
1. Print of FAIRD body-length classes in 3-2b (support of 3.2 claim 4).
2. Run of S3 Part C (ID detection classes on outage days).

## 8. Script map

| V6 section | Script | Main outputs (`outputs/tables/`) |
| --- | --- | --- |
| Shared | `00_load_data.R`, `00_custom_colors.R`, `00_trajectory_functions.R` | |
| 2.4, Table 1 | `2-4_data_availability.R` | `2-4_data_availability_table1_incidents.csv`, `_valid_by_system.csv`, `_valid_by_unit.csv` |
| 3.1, Figure 3, trajectory tests | `3-1a_quantity_trajectories.R` | `3-1a_quantity_trajectories_pair_models.csv`, `_internal_parallelism.csv`, `_test_model_ANOVA_III.csv`, `_figure_model_coefficients.csv` |
| 3.1, concordance, Spearman | `3-1b_quantity_concordance.R` | `3-1b_quantity_concordance_components.csv`, `_disattenuation.csv`, `_daily_mean_correlations.csv` |
| 3.2, Figure 4, Table 2 | `3-2a_composition_orders.R` | `3-2a_composition_orders_table2_order_unit_days.csv`, `_multinomial_LR_tests.csv` |
| 3.2, Izsak-Price, Table 3 | `3-2b_composition_overlap.R` | `3-2b_overlap_summary.csv`, `3-2b_size_threshold.csv`, `3-2b_diptera_share.csv`, `3-2b_identification_*.csv` |
| S1 to S4 | `R/supplementary/S1_...` to `S4_...` | `S1_pair_model_sensitivity_*.csv`, `S2_unit_pair_correlations.csv`, `S3_site_diagnostics_id_class_*.csv`, `S4_abundance_cameras_*.csv` |
