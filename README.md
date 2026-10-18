# Place-based structural disadvantage and mortality inequality in Poland: spatial Mundlak GAM, 2010–2025

R analysis pipeline for the manuscript *Place-based structural disadvantage and mortality inequality across 379 Polish counties, 2010–2025: a spatial Mundlak analysis* (submitted/prepared for **Social Science & Medicine**, Health Geography).

> **Publication note:** This is a prepared source-code update, not a certified end-to-end reproducibility release. Model estimation has **not** been rerun for this GitHub update. Make repository changes public only after checking double-anonymized review requirements.

---

## Overview

The study analyses age–sex-standardised all-cause mortality (ASMR per 100,000) in **379 harmonised Polish counties**. Three **overlapping** analysis windows are used: **M1 2010–2019**, **M2 2010–2022**, and **M3 2010–2025** (estimation observations begin in 2011 after lagging). The spatial Mundlak GAM separates persistent **between-county** contrasts from lagged **within-county** ecological co-movement, with non-linear national functions and spatially varying coefficient (SVC) modifiers. Neither component is interpreted as establishing individual-level causation.

Main features:

1. Nine structural between-county components and eight lagged within-county components.
2. `mgcv::bam`, fREML, `gamma = 1.4`.
3. Retained basis dimensions: national 5; SVC 25; baseline spatial 175; calendar time 9.
4. A linear Mundlak ridge benchmark, residual diagnostics, strict leave-one-year-out temporal validation, local coefficient-covariance uncertainty, same-fit SVC-zeroing comparisons, and 50 structured-noise perturbations per model.
5. Descriptive mortality profiles for 2010–2019, 2020–2022 and 2023–2025.

---

## Repository structure

```text
R/                       Analysis, validation, reporting and figure scripts
  RUN_ALL.R               Full pipeline; computationally intensive—do NOT run for documentation-only updates
  RUN_FIGURES_ONLY.R      Figure-specific workflow; requires precomputed local outputs
  RUN_SVC_REPORTING_ONLY.R  SVC reporting workflow; requires precomputed local outputs
  00_config.R             Settings and period definitions
  02_spatial_mundlak_gam.R  Spatial Mundlak GAM
  03_validation.R         Validation routines
  04_perturbation.R       Structured-noise sensitivity
  05_figures.R            Main figure generation
  06_collect_results.R    Reporting and table assembly
data/                     Dataset access notes; actual county-year data not committed
docs/                     Reproducibility and data-deposit notes
output*/                  Generated locally; not tracked by Git
```

The original R scripts are retained from the supplied latest `R (3)(1).zip`; they may rely on a matching local project directory, source data and saved model objects.

## Reproducibility

**Do not rerun the full pipeline to publish this code.** To reproduce the estimates from scratch, obtain the study data and necessary spatial files, install R and required packages, confirm paths in `R/00_config.R`, and then run from the project root:

```r
source("R/RUN_ALL.R")
```

This command can be computationally expensive. It is not required to upload the source code to GitHub. Recorded model fits and reported findings should be independently validated against the complete released replication materials once available.

## Data

The local master entry point is `data/finalny_panel_przestrzenny_379_powiatow.rds`, with an optional alternative `data/finalny_panel_przestrzenny_379_powiatow.csv.gz`. These data files were **not supplied in the GitHub-preparation package**. The public data-deposit record and DOI must be inserted **only after a deposit exists**. See `data/README_data.md`.

Source datasets include Statistics Poland BDL and national censuses, CAMS EAC4 and CAMS European air-quality reanalyses, and county GIS boundaries. Redistribution of all derived data and source geometry must be checked against their original conditions.

## Interpretation and limitations

The lagged within-county components capture conditional ecological temporal coupling, not exogenous causal shocks. Education trajectories interpolate/extrapolate the 2011/2021 census anchors. The overlapping models are sensitivity windows, not independent pre- and post-COVID comparisons. Strict temporal leave-one-year-out validation showed severe instability in some pandemic-inclusive folds; the GAM should not be represented as a forecasting model.

## Citation

Manuscript citation and data/code DOI: **to be added once registered**. Do not insert placeholder DOI as if it were valid.

## License

The **existing GitHub repository has an MIT license**. Preserve its LICENSE and history when updating; no LICENSE replacement is included in this patch. Dataset license/rights are separate from the code license.


## Curated research code (2026 SSM submission)

`R/RUN_ALL.R` is the canonical end-to-end entry point. It is computationally expensive; **do not run it merely to upload this repository**. Scripts `R/07_basis_dimension_grid_48.R` and `R/08_baseline_boundary_check_175_200.R` document the reported basis checks and are intentionally not invoked by the standard entry point. `R/05e_svc_zeroing_sensitivity.R` documents the separately executed sensitivity analysis.

One potential provenance mismatch requires checking before claiming full one-click replication: the currently active `RUN_ALL.R` invokes `05c_four_period_mortality_profiles.R`, whose reporting additionally includes the all-years window, while the manuscript focuses on three descriptive eras. The manuscript and the exact final CSVs are the authority for the published results.

See `docs/CLEANUP_AUDIT.md` for retained and excluded scripts. The final analysis panel is **not included** in this code repository.
