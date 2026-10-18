# ============================================================
# 02_spatial_mundlak_gam.R — Spatial Mundlak GAM / SVC
# v40.5: selected architecture/K + M1-reference scaling
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("build_mundlak_panel")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(stringr); library(mgcv); library(readr)
})

OUT_DIR <- file.path(OUT_BASE, "spatial_mundlak_gam")
MOD_DIR <- file.path(OUT_DIR, "models")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(MOD_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== SPATIAL MUNDLAK GAM v40.5 ===\n")
cat(MODEL_LABEL, "|", YEAR_LABEL, "| controls:",
    ifelse(length(CONTROL_VARS), paste(CONTROL_VARS, collapse = ", "), "none"), "\n")
cat("Architecture:", GAM_ARCHITECTURE,
    "| K global/SVC/baseline/year =", K_GLOBAL, "/", K_SVC, "/", K_BASELINE, "/", K_YEAR, "\n")
cat("Bases global/SVC/baseline/year =", BS_GLOBAL, "/", BS_SVC, "/", BS_BASELINE, "/", BS_YEAR, "\n")
cat("Scaling:", SCALING_MODE, "| reference:", SCALING_REFERENCE_PERIOD, "\n")

prep <- build_mundlak_panel(dat, period_id = PERIOD, scaling_mode = SCALING_MODE)
dat_all <- prep$data
between_terms <- prep$between_terms
within_terms <- prep$within_terms
all_terms <- prep$all_terms
readr::write_csv(prep$scaler_table, file.path(TAB_DIR, "scaling_parameters.csv"))

dat_model <- dat_all %>%
  model_complete_cases(OUTCOME, all_terms, CONTROL_VARS, require_xy = TRUE) %>%
  arrange(code, year)

if (nrow(dat_model) == 0L) stop("Brak kompletnych obserwacji do GAM.")
cat("n model =", nrow(dat_model), "| counties =", n_distinct(dat_model$code),
    "| years =", min(dat_model$year), "–", max(dat_model$year), "\n")

fml <- build_gam_formula(between_terms, within_terms, CONTROL_VARS)
cat("Formula:\n")
print(fml)

m <- mgcv::bam(
  fml,
  data = dat_model,
  family = FAMILY,
  method = BAM_METHOD,
  discrete = BAM_DISCRETE,
  gamma = GAMMA,
  nthreads = BAM_NTHREADS,
  gc.level = BAM_GC_LEVEL,
  na.action = na.exclude,
  control = gam.control(trace = FALSE)
)

sm <- summary(m)
cat("R² =", round(sm$r.sq, 3), "| dev.expl =", round(sm$dev.expl, 3), "\n")
saveRDS(m, file.path(MOD_DIR, "gam_model.rds"))

# ------------------------------------------------------------
# EDF / smooth diagnostics table
# ------------------------------------------------------------
edf_tbl <- as_tibble(sm$s.table, rownames = "term") %>%
  mutate(
    term_type = case_when(
      str_detect(term, "year_c") ~ "year",
      str_detect(term, "s\\(x,y\\)") & str_detect(term, ":") ~ "SVC",
      str_detect(term, "s\\(x,y\\)") ~ "baseline",
      TRUE ~ "global"
    ),
    raw_var = case_when(
      term_type == "SVC" ~ str_remove(term, "^.*:"),
      TRUE ~ str_extract(term, "(?<=s\\()[^,)]+")
    ),
    component = comp_from_term(raw_var),
    variable = clean_var(raw_var),
    var_label = label_var(variable)
  )
readr::write_csv(edf_tbl, file.path(TAB_DIR, "edf_table.csv"))

# ------------------------------------------------------------
# Prediction helpers: one shared reference regime
# other determinant terms = 0; year_c = 0; period dummies = 0.
# ------------------------------------------------------------
ref_counties <- make_reference_counties(dat_model, all_terms, CONTROL_VARS)

predict_delta <- function(mod, nd0, nd1) {
  p0 <- tryCatch(as.numeric(predict(mod, newdata = nd0, type = "response")),
                 error = function(e) rep(NA_real_, nrow(nd0)))
  p1 <- tryCatch(as.numeric(predict(mod, newdata = nd1, type = "response")),
                 error = function(e) rep(NA_real_, nrow(nd1)))
  p1 - p0
}

term_grid <- function(var, n_grid = N_GRID) {
  q <- quantile(dat_model[[var]], c(.02, .98), na.rm = TRUE)
  if (!all(is.finite(q)) || q[1] == q[2]) return(numeric(0))
  seq(q[1], q[2], length.out = n_grid)
}

# ------------------------------------------------------------
# Average partial slopes (county-specific -> median/IQR)
# ------------------------------------------------------------
compute_aps <- function(var) {
  vg <- term_grid(var)
  if (length(vg) == 0L) return(NULL)

  slope_grid <- map_dfr(vg, function(vv) {
    nd0 <- ref_counties
    nd1 <- ref_counties
    nd0[[var]] <- vv
    nd1[[var]] <- vv + SLOPE_EPS
    tibble(
      code = ref_counties$code,
      code_pow = ref_counties$code_pow,
      name = ref_counties$name,
      grid_value = vv,
      partial_slope = predict_delta(m, nd0, nd1) / SLOPE_EPS
    )
  })

  county <- slope_grid %>%
    group_by(code, code_pow, name) %>%
    summarise(aps_county = mean(partial_slope, na.rm = TRUE), .groups = "drop")

  summary_tbl <- county %>%
    summarise(
      aps = median(aps_county, na.rm = TRUE),
      aps_mean = mean(aps_county, na.rm = TRUE),
      aps_q25 = quantile(aps_county, .25, na.rm = TRUE),
      aps_q75 = quantile(aps_county, .75, na.rm = TRUE),
      aps_p10 = quantile(aps_county, .10, na.rm = TRUE),
      aps_p90 = quantile(aps_county, .90, na.rm = TRUE),
      n_counties = sum(is.finite(aps_county))
    ) %>%
    mutate(
      var = var,
      variable = clean_var(var),
      var_label = label_var(variable),
      component = comp_from_term(var),
      period = PERIOD,
      year_label = YEAR_LABEL
    )

  list(summary = summary_tbl,
       county = county %>% mutate(var = var, variable = clean_var(var),
                                  var_label = label_var(clean_var(var)),
                                  component = comp_from_term(var), period = PERIOD),
       grid = slope_grid %>% mutate(var = var, variable = clean_var(var),
                                    component = comp_from_term(var), period = PERIOD))
}

cat("APS...\n")
aps_res <- map(all_terms, ~ tryCatch(compute_aps(.x), error = function(e) {
  warning("APS ", .x, ": ", e$message); NULL
}))
aps_all <- map_dfr(aps_res, "summary") %>%
  mutate(direction = case_when(
    aps_q75 < 0 ~ "lower outcome",
    aps_q25 > 0 ~ "higher outcome",
    TRUE ~ "heterogeneous / crosses zero"
  ))
aps_county_all <- map_dfr(aps_res, "county")
aps_grid_county_all <- map_dfr(aps_res, "grid")
readr::write_csv(aps_all, file.path(TAB_DIR, "aps_all.csv"))
readr::write_csv(aps_county_all, file.path(TAB_DIR, "aps_county_all.csv"))

# ------------------------------------------------------------
# Effect curves
# ------------------------------------------------------------
compute_curve <- function(var) {
  vg <- term_grid(var)
  if (length(vg) == 0L) return(NULL)
  nd_ref <- ref_counties
  nd_ref[[var]] <- 0
  p_ref <- as.numeric(predict(m, newdata = nd_ref, type = "response"))

  map_dfr(vg, function(vv) {
    nd <- ref_counties
    nd[[var]] <- vv
    eff <- as.numeric(predict(m, newdata = nd, type = "response")) - p_ref
    tibble(
      var_value = vv,
      effect_median = median(eff, na.rm = TRUE),
      effect_q25 = quantile(eff, .25, na.rm = TRUE),
      effect_q75 = quantile(eff, .75, na.rm = TRUE),
      effect_p10 = quantile(eff, .10, na.rm = TRUE),
      effect_p90 = quantile(eff, .90, na.rm = TRUE)
    )
  }) %>%
    mutate(var = var,
           variable = clean_var(var),
           var_label = label_var(variable),
           component = comp_from_term(var),
           period = PERIOD,
           year_label = YEAR_LABEL)
}

cat("Effect curves...\n")
curves_all <- map_dfr(all_terms, ~ tryCatch(compute_curve(.x), error = function(e) {
  warning("Curve ", .x, ": ", e$message); NULL
}))
readr::write_csv(curves_all, file.path(TAB_DIR, "effect_curves.csv"))

# ------------------------------------------------------------
# Marginal and accumulated-effect maps
# Between: county structural level vs 0.
# Within: mean over observed lagged within deviations vs 0.
# ------------------------------------------------------------
compute_effect_maps <- function(var) {
  comp <- comp_from_term(var)

  if (comp == "between") {
    vals <- dat_model %>%
      select(code, code_pow, name, all_of(var)) %>%
      distinct(code, .keep_all = TRUE) %>%
      rename(v_actual = all_of(var))
    base <- ref_counties %>%
      inner_join(vals %>% select(code, v_actual), by = "code")

    nd0 <- base; nd1 <- base
    nd0[[var]] <- base$v_actual
    nd1[[var]] <- base$v_actual + SLOPE_EPS
    deriv <- predict_delta(m, nd0, nd1) / SLOPE_EPS

    nd_ref <- base; nd_act <- base
    nd_ref[[var]] <- 0
    nd_act[[var]] <- base$v_actual
    delta <- predict_delta(m, nd_ref, nd_act)

    return(bind_rows(
      base %>% transmute(code, code_pow, name, x, y, v_actual,
                         value = deriv, map_type = "deriv",
                         contribution_definition = NA_character_),
      base %>% transmute(code, code_pow, name, x, y, v_actual,
                         value = delta, map_type = "delta",
                         contribution_definition = "observed county structural level vs 0 SD")
    ) %>% mutate(var = var, variable = clean_var(var), var_label = label_var(variable),
                 component = comp, period = PERIOD, year_label = YEAR_LABEL))
  }

  # within
  # Marginal-effect map: average local derivative over the county's observed
  # lagged within deviations.
  meta <- dat_model %>%
    select(code, code_pow, name, year, x, y, all_of(var)) %>%
    rename(v_actual = all_of(var))
  nd_base <- meta
  nd_base$year_c <- 0
  for (cc in CONTROL_VARS) nd_base[[cc]] <- 0
  for (tt in all_terms) nd_base[[tt]] <- 0

  nd0 <- nd_base; nd1 <- nd_base
  nd0[[var]] <- meta$v_actual
  nd1[[var]] <- meta$v_actual + SLOPE_EPS
  deriv_it <- predict_delta(m, nd0, nd1) / SLOPE_EPS

  pow_deriv <- meta %>%
    mutate(deriv_it = deriv_it) %>%
    group_by(code, code_pow, name) %>%
    summarise(
      x = first(x), y = first(y),
      v_actual = mean(v_actual, na.rm = TRUE),
      value = mean(deriv_it, na.rm = TRUE),
      n_years = sum(is.finite(v_actual)),
      .groups = "drop"
    ) %>%
    mutate(
      map_type = "deriv",
      contribution_definition = NA_character_
    )

  # Accumulated contribution for a within term cannot be the mean of
  # f(within_it)-f(0) over time: within deviations are centred around zero,
  # so that signed average mechanically cancels for near-linear effects.
  # Use a common standardized scenario instead: +1 SD lagged within deviation
  # versus 0 SD, evaluated separately for each county.
  ref <- ref_counties
  nd_ref <- ref
  nd_plus1 <- ref
  nd_ref[[var]] <- 0
  nd_plus1[[var]] <- 1
  delta_plus1 <- predict_delta(m, nd_ref, nd_plus1)

  nyrs <- meta %>%
    group_by(code, code_pow, name) %>%
    summarise(n_years = sum(is.finite(v_actual)), .groups = "drop")

  pow_delta <- ref %>%
    transmute(
      code, code_pow, name, x, y,
      v_actual = 1,
      value = delta_plus1,
      map_type = "delta",
      contribution_definition = "+1 SD lag-1 within deviation vs 0 SD"
    ) %>%
    left_join(nyrs, by = c("code", "code_pow", "name"))

  bind_rows(pow_deriv, pow_delta) %>%
    mutate(
      var = var,
      variable = clean_var(var),
      var_label = label_var(variable),
      component = comp,
      period = PERIOD,
      year_label = YEAR_LABEL
    )
}

cat("Effect maps...\n")
effect_maps <- map_dfr(all_terms, ~ tryCatch(compute_effect_maps(.x), error = function(e) {
  warning("Map ", .x, ": ", e$message); NULL
}))
readr::write_csv(effect_maps %>% select(-any_of(c("geometry"))),
                 file.path(TAB_DIR, "effect_maps.csv"))

# ------------------------------------------------------------
# In-sample fitted summary by county
# ------------------------------------------------------------
dat_fit <- dat_model %>%
  mutate(
    fitted = as.numeric(predict(m, newdata = dat_model, type = "response")),
    residual = .data[[OUTCOME]] - fitted
  )
pred_powiat <- dat_fit %>%
  group_by(code, code_pow, name) %>%
  summarise(
    observed_mean = mean(.data[[OUTCOME]], na.rm = TRUE),
    predicted_mean = mean(fitted, na.rm = TRUE),
    residual_mean = mean(residual, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(period = PERIOD, year_label = YEAR_LABEL)
readr::write_csv(pred_powiat, file.path(TAB_DIR, "powiat_predictions_summary.csv"))

fit_summary <- tibble(
  period = PERIOD,
  year_label = YEAR_LABEL,
  n_model = nrow(dat_model),
  n_counties = n_distinct(dat_model$code),
  year_min_model = min(dat_model$year),
  year_max_model = max(dat_model$year),
  r_sq = sm$r.sq,
  dev_expl = sm$dev.expl,
  rmse_insample = sqrt(mean((dat_model[[OUTCOME]] - fitted(m))^2, na.rm = TRUE)),
  aic = AIC(m),
  bic = BIC(m),
  total_smooth_edf = sum(sm$s.table[, "edf"], na.rm = TRUE),
  freml_score = if (!is.null(m$gcv.ubre) && length(m$gcv.ubre)) unname(m$gcv.ubre[1]) else NA_real_,
  controls = paste(CONTROL_VARS, collapse = ";"),
  scaling_mode = prep$scaling_mode,
  scaling_source = prep$scaler_source,
  architecture = GAM_ARCHITECTURE,
  bs_global = BS_GLOBAL,
  bs_svc = BS_SVC,
  bs_baseline = BS_BASELINE,
  bs_year = BS_YEAR,
  k_global = K_GLOBAL,
  k_svc = K_SVC,
  k_baseline = K_BASELINE,
  k_year = K_YEAR,
  gamma = GAMMA
)
readr::write_csv(fit_summary, file.path(TAB_DIR, "gam_fit_summary.csv"))

saveRDS(list(
  m = m,
  formula = fml,
  dat_all = dat_all,
  dat_model = dat_model,
  prep = prep,
  between_terms = between_terms,
  within_terms = within_terms,
  all_terms = all_terms,
  edf_tbl = edf_tbl,
  aps_all = aps_all,
  aps_county_all = aps_county_all,
  aps_grid_county_all = aps_grid_county_all,
  curves_all = curves_all,
  effect_maps = effect_maps,
  pred_powiat = pred_powiat,
  fit_summary = fit_summary,
  VAR_LABELS = VAR_LABELS,
  VAR_UNITS = VAR_UNITS,
  period = PERIOD,
  year_label = YEAR_LABEL,
  settings = list(
    OUTCOME = OUTCOME,
    CONTROL_VARS = CONTROL_VARS,
    WITHIN_LAG_N = WITHIN_LAG_N,
    YEAR_MIN = YEAR_MIN, YEAR_MAX = YEAR_MAX,
    YEAR_LABEL = YEAR_LABEL,
    GAMMA = GAMMA,
    K_GLOBAL = K_GLOBAL, K_SVC = K_SVC, K_BASELINE = K_BASELINE, K_YEAR = K_YEAR,
    GAM_ARCHITECTURE = GAM_ARCHITECTURE,
    BS_GLOBAL = BS_GLOBAL, BS_SVC = BS_SVC, BS_BASELINE = BS_BASELINE, BS_YEAR = BS_YEAR,
    SCALING_MODE = SCALING_MODE, SCALING_REFERENCE_PERIOD = SCALING_REFERENCE_PERIOD
  )
), file.path(OUT_DIR, "gam_results.rds"))

cat("\nSpatial Mundlak GAM DONE:", OUT_DIR, "\n")
print(aps_all %>% select(var_label, component, aps, aps_q25, aps_q75), n = 50)
