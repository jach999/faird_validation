# FAIR-DEVICE VALIDATION STUDY, OUTLINE OF MANUSCRIPT V6

**Version:** 0.1 (agreed in session 2026-10-13/14; written down 2026-10-27).
**Purpose:** structure of the V6 rewrite and one sentence per claim, to be used as
paragraph headings when drafting. Each claim cites the registry entry that supports it
(`docs/RESULTS_REGISTRY.md`); the registry, not this file, holds the values and their
status. Scherber's guidelines (G1 to G18, emails E1 to E7) are in
`claude/20261013_SCHERBER_GUIDELINES_EN.md`.

Status marks used below: **[confirmed]** agreed by Juan; **[proposed]** proposed by
Claude, not yet confirmed; **[pending analysis]** needs a run before the text can be
written; **[open]** decision not taken.

---

## 1. Principles of the rewrite

- No hypothesis-driven structure (no H1/H2/H3) and no V5 level hierarchy (Site,
  Intra-Habitat, Overall) as section titles. The V5 terms are too rigid for the
  analyses now planned, which take a general view rather than an exhaustive test of
  consistency. [confirmed]
- Results follow Scherber's two main figures (G15): one section on quantity, one on
  composition. [confirmed]
- System-level comparisons only in the main text (G6, G8); device-pair comparisons go to
  the supplementary material. [confirmed]
- Figures on linear axes, AMMOD on a secondary axis with a data-derived factor (G14).
  [confirmed]

### Vocabulary

Defined once in Methods and used consistently. [confirmed]

| Term | Meaning | Registry / script name |
| --- | --- | --- |
| seasonal trajectory | shape of a system's curve over the period | LOW |
| day-to-day fluctuation | daily deviations around that curve | HIGH |
| proportional | two systems whose curves keep a constant ratio (parallel on the log link) | parallelism test, Section 2e |
| across both habitats / in Maize | Maize and Meadow together / Maize only | two-block / Maize-only rows |
| system / unit | FAIRD, AMMOD or ID / one device of a system | Device_type / Device |

---

## 2. Materials and Methods

- **2.1 Devices.** FAIRD, AMMOD, ID; power supply.
- **2.2 Study site and field design.** Four sites, two habitats, 12 units, 22 days; one
  unit per system and site, stated as the replication limit (G10).
- **2.3 Data processing.** Per system; FAIRD and ID biomass estimated from counts and
  allometry (Söhlström), AMMOD live mass; daily recording coverage (AMMOD 00-24 h, FAIRD
  07-19 h, ID three 2 h windows) [R-DX-COVER-01]; abundance exists only for the camera
  systems (AMMOD taxonomy is presence/absence).
- **2.4 Data availability and exclusions.** The whole failure story in one place
  [confirmed; not split between Methods and Results]: one table per system with
  unit-days planned and lost, type of failure and known cause (answers G5); the ID power
  outages coincided with the cold spell of 28 August to 1 September [R-DX-WX-01,
  R-DX-ID-WIN-01, R-DX-FAIL-01]; exclusion rule for the ID battery window (full series
  vs window excluded, symmetric for both systems of a comparison); effective sample
  sizes. The Discussion refers to the table and does not repeat the facts.
  [pending analysis: counts for the table; AMMOD losses classified with the same
  criterion as ID]
- **2.5 Statistical analysis.** Tweedie GLMM of trajectories (G11); the two concordance
  components; multinomial model of order composition; Izsak-Price similarity in 3-day
  windows at Superfamily level. Details in the supplementary material. The multinomial
  model has no unit random effect (G13); the text says so as part of the replication
  statement.

---

## 3. Results

### 3.1 Quantity

Figure 2: biomass of the three systems over time (Scherber's figure, corrected inputs,
dual-axis version). [confirmed]

1. *The two halves of the FAIRD network followed proportional trajectories and shared
   their day-to-day fluctuation.* [R-TB-INT-01, R-CC-INT-01; FAIRD halves across both
   habitats] Opening paragraph, as the internal reference.
2. *AMMOD recorded substantially higher daily biomass than both camera systems.*
   [R-DC-BM-02, R-OBS-07; contrast at the reference date, see R-OBS-03 before quoting a
   single ratio]
3. *FAIRD biomass was not proportional to AMMOD biomass over time.* [R-TB-AF-01]
4. *FAIRD and AMMOD shared the day-to-day fluctuation of biomass, partly through a
   common response to the cold spell.* [R-CC-AF-01, R-CC-SENS-01, R-OBS-04]
5. *ID showed a similar day-to-day agreement with AMMOD, but whether its biomass was
   proportional to AMMOD could not be determined.* [R-CC-AI-01, R-TB-AI-01]
6. *FAIRD and ID shared the timing of the seasonal peak, but diverged in the second half
   of the period.* [R-TB-FI-22F, R-TB-FI-22X, R-CC-FI-22-01, R-CC-SENS-01]
7. Closing paragraph on abundance, FAIRD vs ID only, as a check on the camera
   comparison: abundance is the native camera measure and separates detection from
   taxonomic resolution, which camera biomass mixes. No figure in the main text.
   [confirmed as a paragraph, not a section; pending analysis: rerun R-TA-FI-02]
   - *ID recorded more individuals per day than FAIRD.* [R-TA-FI-01, PROVISIONAL]
   - Claim on trajectories: to be written after the rerun. If the rerun adds nothing to
     claim 6, the paragraph shrinks to one sentence with a reference to the
     supplementary material.

### 3.2 Composition

Figure 1: order composition over time, one panel per system (Scherber's figure,
corrected inputs). [confirmed]

1. *The three systems differed in which orders they detected, with no common temporal
   trend in composition.* [R-DC-TX-02; descriptive, p-values anti-conservative]
2. *Diptera were detected on almost every unit-day by all systems; Lepidoptera and
   Hemiptera mainly by AMMOD.* [R-DC-TX-03, R-OBS-06; pending analysis: valid unit-days
   per system as denominators]
3. *FAIRD and AMMOD communities overlapped by about two thirds at Superfamily level, more
   in Meadow than in Maize.* [R-V5-H2T-01; kept by Juan's decision, confirmed; to be
   reclassified in the registry, see D4]
4. *Excluding AMMOD taxa smaller than 3 mm barely changed this overlap; FAIRD recorded
   almost no individuals below that size.* [R-V5-H2T-02, R-OUT-H2Q-04]
5. *FAIRD covered a broader taxonomic range than ID, consistent with their different
   classification approaches, but a small fraction of the genera detected by AMMOD.*
   [R-V5-H2T-03; wording on ID descriptive of its design, not a deficiency]
6. *Diptera made up about 90% of FAIRD records and had the lowest identification success
   to genus.* [R-V5-H2T-04]
7. Venn (Euler) diagram of composition (G9). [open: common taxonomic level for the three
   systems]

---

## 4. Discussion

1. **What camera systems capture of the reference signal.** FAIRD follows the
   day-to-day fluctuation, not the amount; a complement to Malaise sampling, not a
   substitute. Qualify E6 point 6 ("high concordance"), which the data do not support as
   stated.
2. **Limits of taxonomic coverage.** Poorly resolved Diptera, the 3 mm threshold,
   recording coverage; candidate causes for the low detection of Lepidoptera and
   Hemiptera [R-OBS-06].
3. **Device reliability in context.** Refers to the table of 2.4; causes, silent failures
   (ID3), non-random data loss; trade-off between field processing and data
   transmission. Juan's position, not a ranking of systems.
4. **Choosing a system.** Device comparison table (G4) and decision tree (G3).
5. **Limitations.** One unit per system and site; AMMOD biomass over 12 days only;
   asymptotic tests.
6. **Outlook.** Light-weight devices combining optical methods and eDNA (E6 point 7).

---

## 5. Supplementary material

AMMOD and ID unit pairs (internal consistency), pair models in detail, identity-link
tests, sensitivities (spline with 4 df, ID battery window), per-site diagnostics,
abundance figure.

---

## 6. Open decisions

| ID | Decision | Proposal | Status |
| --- | --- | --- | --- |
| D1 | Internal consistency of AMMOD and ID exists only as Maize unit pairs (G8) | Supplementary material; main text uses FAIRD halves only | [proposed] |
| D2 | Figure 2 and its test should come from the same model | Scherber's model (G11, G13), adding a test of the Date x Device_type interaction; the Block/AR1 pair models as sensitivity | [proposed] |
| D3 | Place of abundance | Closing paragraph of 3.1 | [confirmed] |
| D4 | Taxonomic entries used in V6 are LEGACY in the registry (R-V5-H2T-01 to -04, R-OUT-H2Q-04) | Reclassify as CURRENT; all are reproduced by the 2026-07-22 outputs | [proposed] |
| D5 | Izsak-Price compares each FAIRD with the co-located AMMOD and treats windows as independent, which Scherber could read as a device-device comparison (G8) | Optionally recompute at system level (union of the units of each system per habitat and window) in `H2_quali_validation.R` | [open] |
| D6 | Taxonomic level of the Venn diagram | Family or Superfamily, AMMOD OTUs aggregated to it, unassigned share stated | [open] |

## 7. Analyses needed before drafting

1. Corrected device-comparison script: prints the factor k of Figure 2 (and the
   interaction test, if D2 is accepted).
2. Abundance rerun R-TA-FI-02 (symmetric exclusion, full `times` levels).
3. Counts for the table of 2.4 (unit-days planned and lost per system and cause).
4. Denominators for R-DC-TX-03 (valid unit-days per system, from `effort_base`).
