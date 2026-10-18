# Place-based mortality inequality in Poland (2010–2025)

Reproducible R workflow accompanying the study *Place-based structural disadvantage and mortality inequality across 379 Polish counties, 2010–2025: a spatial Mundlak analysis*.

## Overview

This study examines age- and sex-standardised all-cause mortality rates (ASMR, per 100,000) for 379 harmonised Polish counties. Three overlapping modelling windows are considered: 2010–2019 (M1), 2010–2022 (M2), and 2010–2025 (M3); estimation starts in 2011 due to lagged predictors.

Spatial Mundlak generalised additive models distinguish persistent between-county associations from lagged within-county deviations. National nonlinear relationships and spatially varying coefficients are estimated using `mgcv::bam` (fREML, `gamma = 1.4`). The selected basis dimensions are 5 (national determinant smooths), 25 (spatial modifiers), 175 (baseline spatial surface), and 9 (calendar time).

Analyses include a linear ridge benchmark, residual and temporal validation, local uncertainty, same-fit spatial-modifier sensitivity, and structured-noise perturbation (50 replicates per model). These are ecological associations, not causal individual-level effects.

## Repository structure

```text
R/
  00_config.R                         Model settings and study periods
  00_helpers.R                        Shared functions
  00_data_audit.R                     Input-data checks
  01_mundlak_ridge.R                  Linear Mundlak benchmark
  02_spatial_mundlak_gam.R            Main spatial Mundlak GAM
  02b_svc_uncertainty.R               Local coefficient-covariance uncertainty
  03_validation.R                     Residual diagnostics and temporal validation
  04_perturbation.R                   Structured-noise sensitivity
  05_figures.R                        Model effect and map figures
  05b_mortality_profiles.R            Mortality-profile analysis
  05c_four_period_mortality_profiles.R  Descriptive mortality-profile reporting
  05d_svc_reporting.R                 Spatially varying coefficient reporting
  05e_svc_zeroing_sensitivity.R       Same-fit spatial-modifier sensitivity
  06_collect_results.R                Tables and aggregate results
  07_basis_dimension_grid_48.R        Basis-dimension sensitivity grid
  08_baseline_boundary_check_175_200.R  Baseline spatial basis check
  RUN_ALL.R                            Full analysis workflow
data/
  README_data.md                      Data inputs, provenance and access requirements
docs/
  QUICK_START.md                      How to prepare and run the analysis
  REPRODUCIBILITY_NOTES.md            Scope, dependencies and validation notes
zenodo_draft/                         Draft dataset metadata, not a published deposit
```

Generated results and the analytical panel are not tracked in this code repository.

## Reproducing the analysis

1. Obtain the final county-level analysis panel and county boundary geometries; see [`data/README_data.md`](data/README_data.md).
2. Install the required R packages and confirm the paths and settings in `R/00_config.R` and `R/RUN_ALL.R`.
3. Set the R working directory to the repository root.
4. To intentionally rebuild the study analyses, run:

```r
source("R/RUN_ALL.R")
```

The full workflow is computationally intensive. The basis-dimension scripts and SVC-zeroing sensitivity are retained as separately executed methodological checks; they are not all invoked by the main runner. See [`docs/REPRODUCIBILITY_NOTES.md`](docs/REPRODUCIBILITY_NOTES.md). The complete workflow has not been independently rerun from this curated release.

## Data availability

The primary input is `data/finalny_panel_przestrzenny_379_powiatow.rds`; the workflow also recognises `data/finalny_panel_przestrzenny_379_powiatow.csv.gz` when available. The data deposit and its persistent identifier will be added after publication and rights verification. Source material includes Statistics Poland, Copernicus Atmosphere Monitoring Service and administrative spatial data.

## Interpretation

Between-county associations describe persistent geographical structure. Within-county terms describe lagged ecological co-movement, including census-derived educational trajectories that are interpolated or extrapolated from 2011 and 2021 observations. Overlapping windows are sensitivity analyses, not independent pre-/post-pandemic cohorts. Temporal validation shows instability across some held-out years; the models should not be presented as forecasting tools.

## Citation

Bibliographic details and permanent identifiers will be added after registration/publication.

## License

Code: MIT, as specified by the existing repository's `LICENSE`. Data rights and licensing are separate.
