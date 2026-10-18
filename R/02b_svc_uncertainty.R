# SSM v6 namespace patch: dplyr verbs explicitly qualified to prevent select()/filter() masking.
# ============================================================
# 02b_svc_uncertainty.R
#
# Post-estimation uncertainty for the Spatial Mundlak SVC-GAM.
#
# No refit is performed.
#
# For every county-specific marginal effect and accumulated contribution,
# this script constructs the corresponding linear contrast from the mgcv
# lpmatrix and propagates the full coefficient covariance matrix:
#
#   L = X(newdata_1) - X(newdata_0)
#   estimate = L beta
#   SE = sqrt(diag(L V L'))
#
# For finite-difference marginal effects, L is divided by SLOPE_EPS.
#
# Pointwise intervals are approximate Bayesian intervals in the usual mgcv
# sense. If vcov(..., unconditional=TRUE) is available, smoothing-parameter
# uncertainty is included; otherwise the script falls back to the fitted
# coefficient covariance matrix. These are not multiplicity-adjusted
# local significance tests.
#
# Outputs:
#   effect_maps_uncertainty.csv
#   svc_uncertainty_summary.csv
#   effect_curve_mean_uncertainty.csv
#   average_derivative_uncertainty.csv
# and an updated gam_results.rds containing these objects.
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) {
  source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
}
if (!exists("make_reference_counties")) {
  source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(mgcv)
})

OUT_DIR <- file.path(OUT_BASE, "spatial_mundlak_gam")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDS_PATH <- file.path(OUT_DIR, "gam_results.rds")

if (!file.exists(RDS_PATH)) {
  stop("Brak gotowego GAM-u: ", RDS_PATH)
}

g <- readRDS(RDS_PATH)
m <- g$m
dat_model <- g$dat_model
between_terms <- g$between_terms
within_terms <- g$within_terms
all_terms <- g$all_terms
controls <- g$settings$CONTROL_VARS

ci_level <- if (exists("SVC_CI_LEVEL")) SVC_CI_LEVEL else 0.95
zcrit <- stats::qnorm(1 - (1 - ci_level) / 2)

cat("=== SVC UNCERTAINTY | ", PERIOD, " ===\n", sep = "")
cat("CI level:", ci_level, "\n")

# ------------------------------------------------------------
# Covariance matrix.
# ------------------------------------------------------------
get_coef_vcov <- function(mod) {
  requested_unconditional <- isTRUE(SVC_USE_UNCONDITIONAL_VCOV)

  if (requested_unconditional) {
    V <- suppressWarnings(
      tryCatch(
        stats::vcov(mod, unconditional = TRUE),
        error = function(e) NULL
      )
    )
    if (!is.null(V) && all(is.finite(V))) {
      return(list(
        V = as.matrix(V),
        mode = "mgcv vcov, unconditional=TRUE"
      ))
    }
  }

  V <- tryCatch(
    stats::vcov(mod, unconditional = FALSE),
    error = function(e) NULL
  )

  if (!is.null(V) && all(is.finite(V))) {
    return(list(
      V = as.matrix(V),
      mode = "mgcv vcov, conditional on smoothing parameters"
    ))
  }

  if (!is.null(mod$Vp) && all(is.finite(mod$Vp))) {
    return(list(
      V = as.matrix(mod$Vp),
      mode = "model$Vp fallback"
    ))
  }

  stop("Nie udało się uzyskać macierzy kowariancji współczynników GAM.")
}

vc <- get_coef_vcov(m)
V <- vc$V
beta <- stats::coef(m)

if (nrow(V) != length(beta) || ncol(V) != length(beta)) {
  stop("Niezgodny wymiar V i coef(m).")
}

cat("Covariance:", vc$mode, "\n")

# ------------------------------------------------------------
# Helpers.
# ------------------------------------------------------------
classify_sign <- function(est, lo, hi) {
  dplyr::case_when(
    !is.finite(est) | !is.finite(lo) | !is.finite(hi) ~ "not estimable",
    lo > 0 ~ "robust positive",
    hi < 0 ~ "robust negative",
    TRUE ~ "uncertain"
  )
}

contrast_from_L <- function(L) {
  if (is.null(dim(L))) L <- matrix(L, nrow = 1L)

  est <- as.numeric(L %*% beta)
  LV <- L %*% V
  v <- rowSums(LV * L)
  v[is.finite(v) & v < 0 & abs(v) < 1e-10] <- 0
  se <- sqrt(pmax(v, 0))

  lo <- est - zcrit * se
  hi <- est + zcrit * se

  tibble(
    value = est,
    se = se,
    ci_low = lo,
    ci_high = hi,
    ci_width = hi - lo,
    ci_excludes_zero = is.finite(lo) & is.finite(hi) &
      ((lo > 0) | (hi < 0)),
    robust_sign = classify_sign(est, lo, hi)
  )
}

make_L <- function(nd0, nd1, scale = 1) {
  X0 <- predict(m, newdata = nd0, type = "lpmatrix")
  X1 <- predict(m, newdata = nd1, type = "lpmatrix")
  (X1 - X0) / scale
}

aggregate_L_by_group <- function(L, group) {
  group <- as.character(group)
  sums <- rowsum(L, group = group, reorder = FALSE)
  counts <- as.numeric(table(factor(group, levels = rownames(sums))))
  sums / counts
}

ref_counties <- make_reference_counties(
  dat_model,
  all_terms,
  controls
)

# ------------------------------------------------------------
# County maps with pointwise uncertainty.
# ------------------------------------------------------------
compute_map_inference <- function(var) {
  comp <- comp_from_term(var)

  if (comp == "between") {
    vals <- dat_model %>%
      dplyr::select(code, code_pow, name, dplyr::all_of(var)) %>%
      dplyr::distinct(code, .keep_all = TRUE) %>%
      dplyr::rename(v_actual = dplyr::all_of(var))

    base <- ref_counties %>%
      dplyr::inner_join(vals %>% dplyr::select(code, v_actual), by = "code")

    # Local finite-difference marginal effect.
    nd0 <- base
    nd1 <- base
    nd0[[var]] <- base$v_actual
    nd1[[var]] <- base$v_actual + SLOPE_EPS

    inf_deriv <- contrast_from_L(
      make_L(nd0, nd1, scale = SLOPE_EPS)
    )

    # Accumulated contribution: observed structural level vs 0 SD.
    nd_ref <- base
    nd_act <- base
    nd_ref[[var]] <- 0
    nd_act[[var]] <- base$v_actual

    inf_delta <- contrast_from_L(
      make_L(nd_ref, nd_act, scale = 1)
    )

    return(dplyr::bind_rows(
      dplyr::bind_cols(
        base %>%
          dplyr::transmute(code, code_pow, name, x, y, v_actual),
        inf_deriv
      ) %>%
        dplyr::mutate(
          map_type = "deriv",
          contribution_definition = NA_character_
        ),
      dplyr::bind_cols(
        base %>%
          dplyr::transmute(code, code_pow, name, x, y, v_actual),
        inf_delta
      ) %>%
        dplyr::mutate(
          map_type = "delta",
          contribution_definition =
            "observed county structural level vs 0 SD"
        )
    ) %>%
      dplyr::mutate(
        var = var,
        variable = clean_var(var),
        var_label = label_var(variable),
        component = comp,
        period = PERIOD,
        year_label = YEAR_LABEL,
        ci_level = ci_level,
        covariance_mode = vc$mode,
        interval_type = "approximate pointwise coefficient-covariance interval"
      ))
  }

  # WITHIN:
  # marginal effect = county mean of annual local derivatives evaluated at
  # observed lagged within deviations. Uncertainty is propagated for the
  # county-mean linear contrast, not averaged from annual SEs.
  meta <- dat_model %>%
    dplyr::select(code, code_pow, name, year, x, y, dplyr::all_of(var)) %>%
    dplyr::rename(v_actual = dplyr::all_of(var))

  nd_base <- meta
  nd_base$year_c <- 0
  for (cc in controls) nd_base[[cc]] <- 0
  for (tt in all_terms) nd_base[[tt]] <- 0

  nd0 <- nd_base
  nd1 <- nd_base
  nd0[[var]] <- meta$v_actual
  nd1[[var]] <- meta$v_actual + SLOPE_EPS

  L_it <- make_L(nd0, nd1, scale = SLOPE_EPS)
  L_county <- aggregate_L_by_group(L_it, meta$code)

  group_order <- rownames(L_county)
  meta_county <- meta %>%
    dplyr::group_by(code) %>%
    dplyr::summarise(
      code_pow = dplyr::first(code_pow),
      name = dplyr::first(name),
      x = dplyr::first(x),
      y = dplyr::first(y),
      v_actual = mean(v_actual, na.rm = TRUE),
      n_years = sum(is.finite(v_actual)),
      .groups = "drop"
    ) %>%
    dplyr::slice(match(group_order, code))

  inf_deriv <- contrast_from_L(L_county)

  rm(L_it, L_county)
  invisible(gc())

  # Standardized accumulated within contrast: +1 SD vs 0 SD.
  nd_ref <- ref_counties
  nd_plus1 <- ref_counties
  nd_ref[[var]] <- 0
  nd_plus1[[var]] <- 1

  inf_delta <- contrast_from_L(
    make_L(nd_ref, nd_plus1, scale = 1)
  )

  nyrs <- meta_county %>%
    dplyr::select(code, n_years)

  delta_meta <- ref_counties %>%
    dplyr::transmute(
      code, code_pow, name, x, y,
      v_actual = 1
    ) %>%
    dplyr::left_join(nyrs, by = "code")

  dplyr::bind_rows(
    dplyr::bind_cols(
      meta_county %>%
        dplyr::select(code, code_pow, name, x, y, v_actual, n_years),
      inf_deriv
    ) %>%
      dplyr::mutate(
        map_type = "deriv",
        contribution_definition = NA_character_
      ),
    dplyr::bind_cols(
      delta_meta,
      inf_delta
    ) %>%
      dplyr::mutate(
        map_type = "delta",
        contribution_definition =
          "+1 SD lag-1 within deviation vs 0 SD"
      )
  ) %>%
    dplyr::mutate(
      var = var,
      variable = clean_var(var),
      var_label = label_var(variable),
      component = comp,
      period = PERIOD,
      year_label = YEAR_LABEL,
      ci_level = ci_level,
      covariance_mode = vc$mode,
      interval_type = "approximate pointwise coefficient-covariance interval"
    )
}

cat("County map uncertainty...\n")
effect_maps_uncertainty <- purrr::map_dfr(
  all_terms,
  ~ tryCatch(
    compute_map_inference(.x),
    error = function(e) {
      warning("Uncertainty map ", .x, ": ", e$message)
      NULL
    }
  )
)

if (!nrow(effect_maps_uncertainty)) {
  stop("Nie udało się policzyć uncertainty dla effect maps.")
}

readr::write_csv(
  effect_maps_uncertainty,
  file.path(TAB_DIR, "effect_maps_uncertainty.csv")
)

svc_uncertainty_summary <- effect_maps_uncertainty %>%
  dplyr::group_by(
    period, year_label,
    variable, var_label, component, map_type
  ) %>%
  dplyr::summarise(
    n_counties = sum(is.finite(value)),
    mean_estimate = mean(value, na.rm = TRUE),
    median_estimate = median(value, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    median_se = median(se, na.rm = TRUE),
    pct_robust_positive =
      100 * mean(robust_sign == "robust positive", na.rm = TRUE),
    pct_robust_negative =
      100 * mean(robust_sign == "robust negative", na.rm = TRUE),
    pct_uncertain =
      100 * mean(robust_sign == "uncertain", na.rm = TRUE),
    pct_ci_excludes_zero =
      100 * mean(ci_excludes_zero, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  svc_uncertainty_summary,
  file.path(TAB_DIR, "svc_uncertainty_summary.csv")
)

# ------------------------------------------------------------
# Mean effect-curve contrast with pointwise uncertainty.
#
# Existing effect_curves.csv remains the descriptive spatial distribution:
# median/IQR/P10-P90 across county-specific contrasts.
# Here we additionally report the mean county contrast and its model-based CI.
# ------------------------------------------------------------
compute_curve_inference <- function(var) {
  dc <- g$curves_all %>%
    dplyr::filter(var == !!var) %>%
    dplyr::arrange(var_value)

  if (!nrow(dc)) return(NULL)

  purrr::map_dfr(dc$var_value, function(vv) {
    nd0 <- ref_counties
    nd1 <- ref_counties
    nd0[[var]] <- 0
    nd1[[var]] <- vv

    L <- make_L(nd0, nd1, scale = 1)
    L_mean <- matrix(colMeans(L), nrow = 1L)
    inf <- contrast_from_L(L_mean)

    dplyr::bind_cols(
      tibble(var_value = vv),
      inf
    )
  }) %>%
    dplyr::rename(
      effect_mean = value,
      effect_mean_se = se,
      effect_mean_ci_low = ci_low,
      effect_mean_ci_high = ci_high,
      effect_mean_ci_width = ci_width,
      effect_mean_ci_excludes_zero = ci_excludes_zero,
      effect_mean_robust_sign = robust_sign
    ) %>%
    dplyr::mutate(
      var = var,
      variable = clean_var(var),
      var_label = label_var(variable),
      component = comp_from_term(var),
      period = PERIOD,
      year_label = YEAR_LABEL,
      ci_level = ci_level,
      covariance_mode = vc$mode
    )
}

cat("Mean curve uncertainty...\n")
curve_mean_uncertainty <- purrr::map_dfr(
  all_terms,
  ~ tryCatch(
    compute_curve_inference(.x),
    error = function(e) {
      warning("Curve uncertainty ", .x, ": ", e$message)
      NULL
    }
  )
)

readr::write_csv(
  curve_mean_uncertainty,
  file.path(TAB_DIR, "effect_curve_mean_uncertainty.csv")
)

# ------------------------------------------------------------
# Average derivative over grid AND counties, with model-based uncertainty.
# This is a linear mean derivative and therefore differs deliberately from
# the manuscript's descriptive median county APS.
# ------------------------------------------------------------
compute_average_derivative_inference <- function(var) {
  grid <- g$curves_all %>%
    dplyr::filter(var == !!var) %>%
    dplyr::distinct(var_value) %>%
    dplyr::arrange(var_value) %>%
    dplyr::pull(var_value)

  if (!length(grid)) return(NULL)

  L_sum <- rep(0, length(beta))

  for (vv in grid) {
    nd0 <- ref_counties
    nd1 <- ref_counties
    nd0[[var]] <- vv
    nd1[[var]] <- vv + SLOPE_EPS
    L <- make_L(nd0, nd1, scale = SLOPE_EPS)
    L_sum <- L_sum + colMeans(L) / length(grid)
  }

  inf <- contrast_from_L(matrix(L_sum, nrow = 1L))

  inf %>%
    dplyr::rename(
      average_derivative_mean = value,
      average_derivative_se = se
    ) %>%
    dplyr::mutate(
      var = var,
      variable = clean_var(var),
      var_label = label_var(variable),
      component = comp_from_term(var),
      period = PERIOD,
      year_label = YEAR_LABEL,
      n_grid = length(grid),
      n_counties = nrow(ref_counties),
      ci_level = ci_level,
      covariance_mode = vc$mode,
      estimand =
        "mean finite-difference derivative averaged over P2-P98 grid and counties"
    )
}

cat("Average derivative uncertainty...\n")
aps_inference <- purrr::map_dfr(
  all_terms,
  ~ tryCatch(
    compute_average_derivative_inference(.x),
    error = function(e) {
      warning("Average derivative uncertainty ", .x, ": ", e$message)
      NULL
    }
  )
)

readr::write_csv(
  aps_inference,
  file.path(TAB_DIR, "average_derivative_uncertainty.csv")
)

# ------------------------------------------------------------
# Save exact reporting metadata.
# ------------------------------------------------------------
uncertainty_manifest <- tibble(
  period = PERIOD,
  year_label = YEAR_LABEL,
  ci_level = ci_level,
  covariance_mode = vc$mode,
  unconditional_requested = isTRUE(SVC_USE_UNCONDITIONAL_VCOV),
  interval_scope = "pointwise",
  multiplicity_adjusted = FALSE,
  robust_sign_definition =
    "approximate pointwise CI excludes zero",
  marginal_estimand =
    "finite-difference local derivative; within maps average the lpmatrix contrast over observed lagged within deviations within county",
  between_contribution_estimand =
    "observed county structural level vs 0 SD",
  within_contribution_estimand =
    "+1 SD lag-1 within deviation vs 0 SD",
  important_caveat =
    "robust sign is an uncertainty display, not a family-wise local significance test"
)

readr::write_csv(
  uncertainty_manifest,
  file.path(TAB_DIR, "svc_uncertainty_manifest.csv")
)

# ------------------------------------------------------------
# Update saved GAM bundle. Point estimates are replaced by the lpmatrix
# estimates so the estimate and its SE come from the same contrast.
# ------------------------------------------------------------
g$effect_maps <- effect_maps_uncertainty
g$effect_maps_uncertainty <- effect_maps_uncertainty
g$svc_uncertainty_summary <- svc_uncertainty_summary
g$curve_mean_uncertainty <- curve_mean_uncertainty
g$aps_inference <- aps_inference
g$uncertainty_manifest <- uncertainty_manifest

saveRDS(g, RDS_PATH)

cat("\nSVC uncertainty DONE:", OUT_DIR, "\n")
print(
  svc_uncertainty_summary %>%
    dplyr::select(
      var_label, component, map_type,
      pct_robust_positive, pct_robust_negative, pct_uncertain
    ),
  n = 100
)
