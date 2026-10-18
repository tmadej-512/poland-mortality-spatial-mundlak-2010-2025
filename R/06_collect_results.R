# ============================================================
# 06_collect_results.R — zbiorcze tabele do manuskryptu
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)

suppressPackageStartupMessages({library(dplyr); library(purrr); library(readr); library(tibble)})

OUT <- file.path(OUTPUT_ROOT, "publication_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
pps <- names(PERIOD_SPECS)

read_if <- function(path) if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else tibble()

aps <- map_dfr(pps, ~ read_if(file.path(OUTPUT_ROOT, .x, "spatial_mundlak_gam", "tables", "aps_all.csv")))
ridge <- map_dfr(pps, ~ read_if(file.path(OUTPUT_ROOT, .x, "mundlak_ridge", "tables", "coef_final.csv")))
gam_fit <- map_dfr(pps, ~ read_if(file.path(OUTPUT_ROOT, .x, "spatial_mundlak_gam", "tables", "gam_fit_summary.csv")))
cv <- map_dfr(pps, ~ read_if(file.path(OUTPUT_ROOT, .x, "validation", "tables", "strict_loo_summary.csv")))
pert <- map_dfr(pps, ~ read_if(file.path(OUTPUT_ROOT, .x, "perturbation", "tables", "snr_table.csv")))
moran <- map_dfr(pps, ~ read_if(file.path(OUTPUT_ROOT, .x, "validation", "tables", "moran_summary.csv")))
acf <- map_dfr(pps, ~ read_if(file.path(OUTPUT_ROOT, .x, "validation", "tables", "temporal_acf1_summary.csv")))
profiles <- map_dfr(pps, ~ read_if(file.path(
  OUTPUT_ROOT, .x, "mortality_profiles", "tables", "profile_table_publication.csv"
)) %>% dplyr::mutate(period = .x))
profile_outcomes <- map_dfr(pps, ~ read_if(file.path(
  OUTPUT_ROOT, .x, "mortality_profiles", "tables", "profile_outcome_summary.csv"
)))
profile_validation <- map_dfr(pps, ~ read_if(file.path(
  OUTPUT_ROOT, .x, "mortality_profiles", "tables", "profile_validation_against_observed_ASMR.csv"
)))

if (nrow(aps)) readr::write_csv(aps, file.path(OUT, "APS_all_periods.csv"))
if (nrow(ridge)) readr::write_csv(ridge, file.path(OUT, "ridge_all_periods.csv"))
if (nrow(gam_fit)) readr::write_csv(gam_fit, file.path(OUT, "GAM_fit_all_periods.csv"))
if (nrow(cv)) readr::write_csv(cv, file.path(OUT, "strict_CV_all_periods.csv"))
if (nrow(pert)) readr::write_csv(pert, file.path(OUT, "perturbation_all_periods.csv"))
if (nrow(moran)) readr::write_csv(moran, file.path(OUT, "Moran_residuals_all_periods.csv"))
if (nrow(acf)) readr::write_csv(acf, file.path(OUT, "ACF1_residuals_all_periods.csv"))
if (nrow(profiles)) readr::write_csv(profiles, file.path(OUT, "mortality_profiles_all_periods.csv"))
if (nrow(profile_outcomes)) readr::write_csv(profile_outcomes, file.path(OUT, "mortality_profile_outcomes_all_periods.csv"))
if (nrow(profile_validation)) readr::write_csv(profile_validation, file.path(OUT, "mortality_profile_validation_all_periods.csv"))

model_manifest <- map_dfr(pps, function(pp) {
  sp <- PERIOD_SPECS[[pp]]
  tibble(
    period = pp,
    model = sp$model_label,
    years = sp$year_label,
    year_min = sp$year_min,
    year_max = sp$year_max,
    controls = ifelse(length(sp$controls), paste(sp$controls, collapse = ";"), "none"),
    covid_dummy_years = paste(range(COVID_YEARS), collapse = "–"),
    post_covid_dummy_years = ifelse("post_covid_dummy" %in% sp$controls,
                                    paste(range(POST_COVID_YEARS), collapse = "–"), "not included"),
    within_lag = WITHIN_LAG_N,
    scaling_mode = SCALING_MODE,
    scaling_reference_period = SCALING_REFERENCE_PERIOD,
    gam_architecture = GAM_ARCHITECTURE,
    k_global = K_GLOBAL,
    k_svc = K_SVC,
    k_baseline = K_BASELINE,
    k_year = K_YEAR,
    ridge_bootstrap_B = BOOT_B,
    perturbation_B = B_PERTURB,
    strict_loo = RUN_STRICT_LOO,
    between_within_variables = paste(MODEL_VARS, collapse = ";"),
    between_only_variables = paste(ONLY_BETWEEN_VARS, collapse = ";")
  )
})
readr::write_csv(model_manifest, file.path(OUT, "analysis_manifest.csv"))

# Jedna tabela interpretacyjna łącząca GAM + ridge + perturbację.
if (nrow(aps)) {
  combined <- aps %>%
    dplyr::left_join(ridge %>% dplyr::select(period, variable, component, beta_debiased, conf.low, conf.high),
              by = c("period", "variable", "component")) %>%
    dplyr::left_join(
      pert %>% dplyr::select(
        period, variable, component, snr, r_shape, mean_r_shape,
        pct_aps_same_sign, p10_surface_cor, mean_county_sign_agreement,
        spatial_stable, final_class
      ),
      by = c("period", "variable", "component")
    )
  readr::write_csv(combined, file.path(OUT, "combined_effect_results.csv"))
}

cat("\nCollected publication tables:", OUT, "\n")


# PRIMARY four-period observed-ASMR/model-contribution profile outputs.
profile4_dir <- file.path(
  OUTPUT_ROOT,
  "four_period_observed_ASMR_model_profiles",
  "tables"
)

copy_if <- function(src, dst) {
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
  }
}

copy_if(
  file.path(
    profile4_dir,
    "model_contribution_HIGH_ASMR_minus_LOW_ASMR.csv"
  ),
  file.path(
    OUT,
    "model_contribution_HIGH_ASMR_minus_LOW_ASMR_four_periods.csv"
  )
)

copy_if(
  file.path(
    profile4_dir,
    "raw_determinant_profiles_LOW_HIGH_observed_ASMR.csv"
  ),
  file.path(
    OUT,
    "raw_determinant_profiles_LOW_HIGH_observed_ASMR_four_periods.csv"
  )
)

copy_if(
  file.path(
    profile4_dir,
    "observed_ASMR_profile_outcomes_and_model_burden.csv"
  ),
  file.path(
    OUT,
    "observed_ASMR_profiles_model_burden_summary.csv"
  )
)

copy_if(
  file.path(
    profile4_dir,
    "profile_methods_record.csv"
  ),
  file.path(
    OUT,
    "observed_ASMR_model_profile_methods_record.csv"
  )
)


# SVC uncertainty and reporting tables.
svc_uncertainty <- map_dfr(
  pps,
  ~ read_if(file.path(
    OUTPUT_ROOT, .x,
    "spatial_mundlak_gam", "tables",
    "svc_uncertainty_summary.csv"
  ))
)
aps_uncertainty <- map_dfr(
  pps,
  ~ read_if(file.path(
    OUTPUT_ROOT, .x,
    "spatial_mundlak_gam", "tables",
    "average_derivative_uncertainty.csv"
  ))
)

if (nrow(svc_uncertainty)) {
  readr::write_csv(
    svc_uncertainty,
    file.path(OUT, "SVC_uncertainty_all_periods.csv")
  )
}
if (nrow(aps_uncertainty)) {
  readr::write_csv(
    aps_uncertainty,
    file.path(OUT, "GAM_average_derivative_CI_all_periods.csv")
  )
}

for (nm in c(
  "SVC_cross_period_surface_stability.csv",
  "SVC_uncertainty_coverage_all_periods.csv"
)) {
  src <- file.path(OUTPUT_ROOT, "publication_tables", nm)
  # File may already be in OUT; no copy needed if identical destination.
  if (file.exists(src)) {
    message("SVC reporting table available: ", src)
  }
}

# v40.13: focal-SVC zeroing sensitivity (post-estimation; no refit).
svc_zeroing <- map_dfr(
  pps,
  ~ read_if(file.path(
    OUTPUT_ROOT, .x,
    "spatial_mundlak_gam", "tables",
    "svc_zeroing_sensitivity.csv"
  ))
)
svc_ribbon <- map_dfr(
  pps,
  ~ read_if(file.path(
    OUTPUT_ROOT, .x,
    "spatial_mundlak_gam", "tables",
    "svc_zeroing_curve_diagnostics.csv"
  ))
)
if (nrow(svc_zeroing)) {
  readr::write_csv(svc_zeroing, file.path(OUT, "SVC_zeroing_sensitivity_all_models.csv"))
}
if (nrow(svc_ribbon)) {
  readr::write_csv(svc_ribbon, file.path(OUT, "SVC_curve_ribbon_diagnostics_all_models.csv"))
}
