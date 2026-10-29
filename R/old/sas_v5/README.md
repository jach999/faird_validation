# V5 diagnostic mixed models (SAS)

Archived on 2026-10-29 from the former `SAS/` folder of the repository root. Not part of
the published repository.

- `FAIRD_ID_AMMOD.sas`: PROC MIXED models A and B of Draft V5 (REML, Kenward-Roger df,
  AR(1)); registry entries R-V5-H2Q-MM-01 and -02.
- `FAIRD_ID_AMMOD-results.html`: SAS results; prints every V5 value of Tables 6 and 7.
- `FAIRD_ID_AMMOD_prepare_sas_data.R`, `load_data.R`: preparation of the SAS inputs
  `FAIRD_ID_AMMOD_sas_biomass_12day.csv` and `FAIRD_ID_AMMOD_sas_counts_22day.csv`.
- `source_tables/` (removed): byte-identical copies of the four CSV files in the root
  `source_tables/`. To rerun `load_data.R` here, copy them back into
  `R/old/sas_v5/source_tables/` (paths anchored by `.here` in this folder).
