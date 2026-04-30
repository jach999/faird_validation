/* -------------------------------------------------------------------------
   FINAL ANALYSIS - PRE-FILTERED DATA STRATEGY
   
   This script relies on CSV files that were already filtered in R.
   No WHERE clauses are applied to date ranges to avoid logic errors.
   ------------------------------------------------------------------------- */

/* 0. CLEANUP MEMORY */
PROC DATASETS LIBRARY=WORK KILL NOLIST; RUN; QUIT;


/* =========================================================================
   ANALYSIS A: BIOMASS VALIDATION
   Input: sas_final_biomass_12day.csv (Days 1-12 -> AMMODs only have data the first 12 days, Sites 1-3 -> AMMOD4 missing)
   Expected Result: F ~ 29.53
   ========================================================================= */

/* IMPORTANT: Check that the path matches your uploaded file location */
FILENAME REF_BIO '/home/u64392136/FAIRD_ID_AMMOD_sas_biomass_12day.csv';

PROC IMPORT DATAFILE=REF_BIO
	DBMS=CSV
	OUT=WORK.biomass_data;
	GETNAMES=YES;
RUN;

TITLE "Analysis A: Biomass Validation (Strict - Sites 1-3)";

PROC MIXED data=WORK.biomass_data plots=studentpanel method=reml;
    CLASS day_num site device_type;
    
    /* MODEL:
       Testing differences in Log Biomass between devices.
       Site is a fixed effect.
    */
    MODEL log_biomass = device_type site device_type*site / ddfm=kr;
    
    /* REPEATED MEASURES:
       AR(1) structure for temporal autocorrelation.
    */
    REPEATED day_num / subject=device_type*site type=ar(1);
    
    /* POST-HOC:
       Pairwise comparisons (AMMOD vs FAIRD, etc.)
    */
    LSMEANS device_type / adjust=bon diff cl;
    LSMEANS device_type*site / slice=site adjust=bon diff cl;
RUN;


/* =========================================================================
   ANALYSIS B: ABUNDANCE COMPARISON (FAIRD vs ID)
   Input: sas_final_counts_22day.csv (Days 1-22, All Sites, No AMMOD)
   Expected Result: F ~ 21
   ========================================================================= */

FILENAME REF_CNT '/home/u64392136/FAIRD_ID_AMMOD_sas_counts_22day.csv';

PROC IMPORT DATAFILE=REF_CNT
	DBMS=CSV
	OUT=WORK.counts_data;
	GETNAMES=YES;
RUN;

TITLE "Analysis B: Abundance Comparison (FAIRD vs ID - All Days)";

PROC MIXED data=WORK.counts_data plots=studentpanel method=reml;
    CLASS day_num site device_type;
    
    /* MODEL:
       Testing differences in Log Abundance (Counts).
    */
    MODEL log_abundance = device_type site device_type*site / ddfm=kr;
    
    REPEATED day_num / subject=device_type*site type=ar(1);
    
    LSMEANS device_type / adjust=bon diff cl;
RUN;

TITLE; /* Clear titles */