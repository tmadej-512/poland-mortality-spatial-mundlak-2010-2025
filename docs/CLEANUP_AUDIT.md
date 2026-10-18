# Cleanup audit — curated SSM R repository

This distribution is assembled from the latest provided R zip and the earlier GitHub preparation package. The original files were not modified.

## Scope

**17 R files retained**: core configuration/helpers; ridge, GAM, SVC uncertainty, annual validation, perturbation, core figure generation, model-specific mortality profiling, outcome-defined mortality profiles, SVC reporting, zeroing sensitivity, collection, two basis checks, audit and one canonical `RUN_ALL.R`.

**Not included** (retained in the author's original local archive):

- `R/00_preflight.R` — auxiliary preflight script
- `R/05c_three_era_mortality_profiles.R` — alternate profile implementation
- `R/05f_patch_fig1_fig4.R` — obsolete version-specific patch runner
- `R/RUN_AFTER_GAM.R` — alternate runner / redundant execution pathway
- `R/RUN_FIGURES_ONLY.R` — alternate runner / redundant execution pathway
- `R/RUN_GAM_ONLY.R` — alternate runner / redundant execution pathway
- `R/RUN_PATCH_ONLY_v40_13.R` — obsolete version-specific patch runner
- `R/RUN_PATCH_ONLY_v40_15.R` — obsolete version-specific patch runner
- `R/RUN_RIDGE_ONLY.R` — alternate runner / redundant execution pathway
- `R/RUN_SVC_REPORTING_ONLY.R` — alternate runner / redundant execution pathway
- `R/09_export_ssm_deduplicated_figures.R` — alternate runner / redundant execution pathway

## Cautions

- `RUN_ALL.R` was historically the active path, and its dependency on `05c_four_period_mortality_profiles.R` is preserved. The latter explicitly includes a fourth all-years reporting window alongside three study eras. Resolve any mismatch with the **published** three-era tables before changing analysis code.
- No GAMs were re-estimated. No data and saved model objects were provided; full execution cannot be validated here.
- This is a curated *code release*, not a fully verified one-click reproduction package. The reader must supply the final county panel.
- Before a public release, compare the `R/` folder against the latest original analysis source and inspect any absolute paths, environment configuration, and source data permissions.
- Never delete the source archive until the published calculations have been independently reproduced.
