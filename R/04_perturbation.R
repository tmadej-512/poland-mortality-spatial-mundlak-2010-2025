# ============================================================
# 04_perturbation.R — FAST structured residual-noise stability test
# v40.5 FAST:
#   * exact same structured-noise idea as before;
#   * exact curve and APS definitions retained;
#   * prediction calls vectorized over the complete evaluation grid;
#   * spatial stability evaluated as the county-specific local slope
#     at the common reference z = 0, rather than integrating a full
#     county × grid surface in every perturbation;
#   * progress is printed AFTER each completed replicate, with ETA.
#
# This is a stability stress test, not bootstrap inference.
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("make_reference_counties")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(mgcv); library(readr); library(ggplot2); library(forcats)
})

set.seed(20260908)

GAM_RDS <- file.path(OUT_BASE, "spatial_mundlak_gam", "gam_results.rds")
if (!file.exists(GAM_RDS)) stop("Brak: ", GAM_RDS)

res <- readRDS(GAM_RDS)
m <- res$m
dat_model <- res$dat_model
all_terms <- res$all_terms
aps_all <- res$aps_all

OUT_DIR <- file.path(OUT_BASE, "perturbation")
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== PERTURBATION STABILITY | FAST FINAL ===\n")
cat(MODEL_LABEL, "|", YEAR_LABEL,
    "| B =", B_PERTURB,
    "| noise scale =", NOISE_SCALE, "× residual SD\n")
cat("bam threads =", BAM_NTHREADS,
    "| gc.level =", BAM_GC_LEVEL,
    "| curve/APS grid =", N_GRID,
    "| spatial reference z =", PERT_SPATIAL_REFERENCE_Z, "\n")
cat("Spatial stability = county local-slope pattern at z=0; APS remains integrated over the full analysis grid.\n")

if (!exists("powiaty_spdf") || is.null(powiaty_spdf)) {
  stop("Finalna perturbacja wymaga powiaty_spdf.")
}
if (!requireNamespace("spdep", quietly = TRUE)) {
  stop("Finalna perturbacja wymaga pakietu 'spdep'.")
}

# ------------------------------------------------------------
# Residual structure
# ------------------------------------------------------------
fitted0 <- as.numeric(predict(m, newdata = dat_model, type = "response"))
dat_r <- dat_model %>%
  mutate(
    fitted0 = fitted0,
    resid = .data[[OUTCOME]] - fitted0
  ) %>%
  arrange(code, year)

sigma_resid <- sd(dat_r$resid, na.rm = TRUE)

rho_t <- dat_r %>%
  group_by(code) %>%
  summarise(
    rho = if (n() >= 3L)
      suppressWarnings(cor(resid[-n()], resid[-1], use = "complete.obs"))
    else NA_real_,
    .groups = "drop"
  ) %>%
  summarise(rho = median(rho, na.rm = TRUE)) %>%
  pull(rho)

if (!is.finite(rho_t)) rho_t <- 0
rho_t <- max(-.95, min(.95, rho_t))

lambda_s <- 0
W_mat <- NULL

if (!is.null(powiaty_spdf) && requireNamespace("spdep", quietly = TRUE)) {
  pow <- sf::st_as_sf(powiaty_spdf)
  code_col <- intersect(c(POW_CODE_COL, "JPT_KOD_JE", "code_pow"), names(pow))[1]

  pow[[code_col]] <- as.character(pow[[code_col]])
  pow[[code_col]] <- ifelse(
    nchar(pow[[code_col]]) == 3L,
    paste0("0", pow[[code_col]]),
    pow[[code_col]]
  )
  pow <- pow %>% arrange(.data[[code_col]])

  nb <- tryCatch(spdep::poly2nb(pow, queen = TRUE), error = function(e) NULL)

  if (!is.null(nb)) {
    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    W_mat <- tryCatch(spdep::listw2mat(lw), error = function(e) NULL)

    if (!is.null(W_mat)) {
      rownames(W_mat) <- colnames(W_mat) <- pow[[code_col]]
    }

    lambda_by_year <- map_dbl(sort(unique(dat_r$year)), function(yr) {
      dd <- dat_r %>% filter(year == yr)
      idx <- match(pow[[code_col]], dd$code_pow)
      rr <- dd$resid[idx]

      if (sum(is.finite(rr)) < 20L) return(NA_real_)

      wr <- tryCatch(
        spdep::lag.listw(lw, rr, zero.policy = TRUE, NAOK = TRUE),
        error = function(e) rep(NA_real_, length(rr))
      )

      ok <- is.finite(rr) & is.finite(wr)
      if (sum(ok) < 20L) return(NA_real_)

      unname(coef(lm(rr[ok] ~ wr[ok]))[2])
    })

    lambda_s <- median(lambda_by_year, na.rm = TRUE)
    if (!is.finite(lambda_s)) lambda_s <- 0
    lambda_s <- max(-.8, min(.8, lambda_s))
  }
}

cat(sprintf(
  "Residual SD=%.3f | median rho_t=%.3f | spatial lambda proxy=%.3f\n",
  sigma_resid, rho_t, lambda_s
))

# ------------------------------------------------------------
# Structured noise generator
# ------------------------------------------------------------
gen_noise <- function(dat_in, scale_val = NOISE_SCALE) {
  counties <- sort(unique(dat_in$code_pow))
  years <- sort(unique(dat_in$year))
  sig <- sigma_resid * scale_val

  innov <- matrix(
    rnorm(length(counties) * length(years), 0, sig),
    nrow = length(counties),
    ncol = length(years),
    dimnames = list(counties, as.character(years))
  )

  eps <- innov

  if (length(years) >= 2L) {
    for (jj in 2:length(years)) {
      eps[, jj] <-
        rho_t * eps[, jj - 1] +
        sqrt(max(0, 1 - rho_t^2)) * innov[, jj]
    }
  }

  if (!is.null(W_mat) && abs(lambda_s) > .01) {
    w_codes <- rownames(W_mat)

    for (jj in seq_along(years)) {
      idx <- match(counties, w_codes)
      ok <- !is.na(idx)

      e_full <- rep(0, length(w_codes))
      e_full[idx[ok]] <- eps[ok, jj]

      w_e <- as.numeric(W_mat %*% e_full)
      eps[ok, jj] <- eps[ok, jj] + lambda_s * w_e[idx[ok]]
    }
  }

  as.data.frame(eps) %>%
    tibble::rownames_to_column("code_pow") %>%
    pivot_longer(-code_pow, names_to = "year", values_to = "noise") %>%
    mutate(year = as.integer(year))
}

# ------------------------------------------------------------
# Common reference regime and baseline objects
# ------------------------------------------------------------
ref_counties <- make_reference_counties(dat_model, all_terms, CONTROL_VARS) %>%
  mutate(.county_index = row_number())

# Exactly the same grid as the unperturbed GAM curves.
base_grids <- res$curves_all %>%
  select(
    var, variable, var_label, component,
    var_value, effect_base = effect_median
  )

base_aps <- aps_all %>%
  select(var, variable, var_label, component, aps_base = aps)

safe_predict <- function(mod, nd) {
  tryCatch(
    as.numeric(predict(mod, newdata = nd, type = "response")),
    error = function(e) {
      warning("Prediction failed: ", e$message)
      rep(NA_real_, nrow(nd))
    }
  )
}

# Baseline spatial local-effect pattern at the common reference z=0.
# This includes the global derivative + the SVC modifier and is therefore
# directly interpretable as the county-specific local marginal slope at z=0.
compute_spatial_z0 <- function(mod, var) {
  nd0 <- ref_counties
  nd1 <- ref_counties
  nd0[[var]] <- PERT_SPATIAL_REFERENCE_Z
  nd1[[var]] <- PERT_SPATIAL_REFERENCE_Z + SLOPE_EPS

  p0 <- safe_predict(mod, nd0)
  p1 <- safe_predict(mod, nd1)

  tibble(
    var = var,
    variable = clean_var(var),
    var_label = label_var(clean_var(var)),
    component = comp_from_term(var),
    code = ref_counties$code,
    code_pow = ref_counties$code_pow,
    name = ref_counties$name,
    local_slope_base = (p1 - p0) / SLOPE_EPS
  )
}

cat("Preparing unperturbed spatial reference patterns at z=0...\n")
base_spatial_z0 <- map_dfr(
  all_terms,
  ~ compute_spatial_z0(m, .x)
)

readr::write_csv(
  base_spatial_z0,
  file.path(TAB_DIR, "baseline_spatial_z0_local_slopes.csv")
)

# ------------------------------------------------------------
# FAST term evaluation
# ------------------------------------------------------------
# One term is evaluated with TWO vectorized prediction calls:
#   pred(z) and pred(z + eps)
# for all 379 counties and all grid values, including z=0.
#
# From the same prediction block we recover:
#   1) the exact median/IQR effect curve;
#   2) county-specific APS averaged across the full N_GRID grid,
#      then the national median APS;
#   3) the county local-slope pattern at z=0.
#
# This preserves the substantive definitions while avoiding hundreds of
# predict() calls per term.
evaluate_term_fast <- function(mod, var, grid, run_id) {
  grid <- sort(unique(as.numeric(grid)))
  grid <- grid[is.finite(grid)]
  if (length(grid) == 0L) return(NULL)

  eval_values <- unique(c(PERT_SPATIAL_REFERENCE_Z, grid))

  nd0 <- bind_rows(lapply(eval_values, function(vv) {
    nd <- ref_counties
    nd[[var]] <- vv
    nd$.eval_value <- vv
    nd
  }))

  nd1 <- nd0
  nd1[[var]] <- nd1[[var]] + SLOPE_EPS

  pred0 <- safe_predict(mod, nd0)
  pred1 <- safe_predict(mod, nd1)

  ev <- nd0 %>%
    transmute(
      .county_index,
      code, code_pow, name,
      eval_value = .eval_value,
      pred0 = pred0,
      slope = (pred1 - pred0) / SLOPE_EPS
    )

  ref_pred <- ev %>%
    filter(eval_value == PERT_SPATIAL_REFERENCE_Z) %>%
    select(.county_index, pred_ref = pred0)

  curve <- ev %>%
    filter(eval_value %in% grid) %>%
    left_join(ref_pred, by = ".county_index") %>%
    mutate(effect = pred0 - pred_ref) %>%
    group_by(eval_value) %>%
    summarise(
      effect = median(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      var_value = eval_value,
      effect = effect,
      var = var,
      variable = clean_var(var),
      var_label = label_var(clean_var(var)),
      component = comp_from_term(var),
      run = run_id
    )

  county_aps <- ev %>%
    filter(eval_value %in% grid) %>%
    group_by(.county_index, code, code_pow, name) %>%
    summarise(
      aps_county = mean(slope, na.rm = TRUE),
      .groups = "drop"
    )

  aps_run <- county_aps %>%
    summarise(aps_pert = median(aps_county, na.rm = TRUE)) %>%
    mutate(
      var = var,
      variable = clean_var(var),
      var_label = label_var(clean_var(var)),
      component = comp_from_term(var),
      run = run_id
    ) %>%
    select(var, variable, var_label, component, run, aps_pert)

  spatial_z0 <- ev %>%
    filter(eval_value == PERT_SPATIAL_REFERENCE_Z) %>%
    transmute(
      var = var,
      variable = clean_var(var),
      var_label = label_var(clean_var(var)),
      component = comp_from_term(var),
      run = run_id,
      code, code_pow, name,
      local_slope = slope
    )

  list(
    curve = curve,
    aps = aps_run,
    spatial_z0 = spatial_z0
  )
}

summarise_spatial_run <- function(spatial_z0_b) {
  spatial_z0_b %>%
    inner_join(
      base_spatial_z0,
      by = c(
        "var", "variable", "var_label", "component",
        "code", "code_pow", "name"
      )
    ) %>%
    group_by(var, variable, var_label, component, run) %>%
    summarise(
      surface_sd_base = sd(local_slope_base, na.rm = TRUE),

      surface_cor =
        if (
          sum(is.finite(local_slope_base) & is.finite(local_slope)) >= 10L &&
          sd(local_slope_base, na.rm = TRUE) > 0 &&
          sd(local_slope, na.rm = TRUE) > 0
        ) {
          suppressWarnings(
            cor(local_slope_base, local_slope, use = "complete.obs")
          )
        } else {
          NA_real_
        },

      surface_rank_cor =
        if (sum(is.finite(local_slope_base) & is.finite(local_slope)) >= 10L) {
          suppressWarnings(
            cor(
              local_slope_base, local_slope,
              method = "spearman",
              use = "complete.obs"
            )
          )
        } else {
          NA_real_
        },

      surface_nrmse = {
        denom <- sd(local_slope_base, na.rm = TRUE)
        if (!is.finite(denom) || denom < SURFACE_SD_MIN) {
          NA_real_
        } else {
          sqrt(mean((local_slope - local_slope_base)^2, na.rm = TRUE)) / denom
        }
      },

      county_sign_agreement = {
        ok <- is.finite(local_slope_base) & is.finite(local_slope)
        if (sum(ok) == 0L) {
          NA_real_
        } else {
          mean(sign(local_slope_base[ok]) == sign(local_slope[ok]))
        }
      },

      spatial_reference_z = PERT_SPATIAL_REFERENCE_Z,
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# One perturbation replicate
# ------------------------------------------------------------
run_single <- function(b) {
  noise <- gen_noise(dat_model)

  dp <- dat_model %>%
    left_join(noise, by = c("code_pow", "year")) %>%
    mutate(noise = replace_na(noise, 0))

  dp[[OUTCOME]] <- dp[[OUTCOME]] + dp$noise

  mp <- tryCatch(
    mgcv::bam(
      formula(m),
      data = dp,
      family = m$family,
      method = BAM_METHOD,
      discrete = BAM_DISCRETE,
      gamma = GAMMA,
      nthreads = BAM_NTHREADS,
      gc.level = BAM_GC_LEVEL,
      na.action = na.exclude,
      control = gam.control(trace = FALSE)
    ),
    error = function(e) {
      warning("Perturbation ", b, ": ", e$message)
      NULL
    }
  )

  if (is.null(mp)) return(NULL)

  term_res <- map(all_terms, function(v) {
    grid <- base_grids %>%
      filter(var == v) %>%
      arrange(var_value) %>%
      pull(var_value)

    if (length(grid) == 0L) return(NULL)

    tryCatch(
      evaluate_term_fast(mp, v, grid, b),
      error = function(e) {
        warning("Perturbation ", b, ", term ", v, ": ", e$message)
        NULL
      }
    )
  })

  term_res <- compact(term_res)

  if (length(term_res) == 0L) {
    rm(mp, dp, noise)
    invisible(gc(full = FALSE))
    return(NULL)
  }

  curves_b <- bind_rows(map(term_res, "curve"))
  aps_b <- bind_rows(map(term_res, "aps"))
  spatial_z0_b <- bind_rows(map(term_res, "spatial_z0"))
  spatial_metrics_b <- summarise_spatial_run(spatial_z0_b)

  rm(mp, dp, noise, term_res, spatial_z0_b)
  invisible(gc(full = FALSE))

  list(
    curves = curves_b,
    aps = aps_b,
    spatial = spatial_metrics_b
  )
}

# ------------------------------------------------------------
# Run perturbations with timing / ETA
# ------------------------------------------------------------
cat("Running perturbations...\n")
run_objects <- vector("list", B_PERTURB)
run_minutes <- rep(NA_real_, B_PERTURB)

for (b in seq_len(B_PERTURB)) {
  t_b <- proc.time()[["elapsed"]]

  run_objects[[b]] <- run_single(b)

  run_minutes[b] <- (proc.time()[["elapsed"]] - t_b) / 60
  mean_min <- mean(run_minutes[seq_len(b)], na.rm = TRUE)
  eta_min <- mean_min * (B_PERTURB - b)

  status <- if (is.null(run_objects[[b]])) "FAILED" else "OK"

  cat(sprintf(
    "  %d / %d | %s | %.2f min | mean %.2f min | ETA %.1f min\n",
    b, B_PERTURB, status, run_minutes[b], mean_min, eta_min
  ))
}

ok_run <- vapply(run_objects, function(z) !is.null(z), logical(1))
n_success <- sum(ok_run)

min_success <- max(20L, floor(0.90 * B_PERTURB))
if (n_success < min_success) {
  stop(
    "Za mało poprawnych replikacji perturbacji: ",
    n_success, "/", B_PERTURB,
    "; wymagane co najmniej ", min_success, "."
  )
}

curves_pert <- bind_rows(lapply(run_objects[ok_run], `[[`, "curves"))
aps_by_run <- bind_rows(lapply(run_objects[ok_run], `[[`, "aps"))
spatial_by_run <- bind_rows(lapply(run_objects[ok_run], `[[`, "spatial"))

rm(run_objects)
invisible(gc(full = TRUE))

# ------------------------------------------------------------
# Curve summaries
# ------------------------------------------------------------
readr::write_csv(
  curves_pert,
  file.path(TAB_DIR, "curves_perturbed.csv")
)

pert_summary <- curves_pert %>%
  group_by(var, variable, var_label, component, var_value) %>%
  summarise(
    effect_median = median(effect, na.rm = TRUE),
    effect_q25 = quantile(effect, .25, na.rm = TRUE),
    effect_q75 = quantile(effect, .75, na.rm = TRUE),
    effect_p10 = quantile(effect, .10, na.rm = TRUE),
    effect_p90 = quantile(effect, .90, na.rm = TRUE),
    effect_sd = sd(effect, na.rm = TRUE),
    n_runs = n_distinct(run),
    .groups = "drop"
  )

shape_by_run <- base_grids %>%
  inner_join(
    curves_pert,
    by = c(
      "var", "variable", "var_label", "component", "var_value"
    )
  ) %>%
  group_by(var, variable, var_label, component, run) %>%
  summarise(
    r_run =
      if (n() >= 3L) {
        suppressWarnings(cor(effect_base, effect, use = "complete.obs"))
      } else {
        NA_real_
      },
    .groups = "drop"
  )

# ------------------------------------------------------------
# APS stability
# APS itself is still calculated over the complete N_GRID grid and then
# summarized across counties exactly as in the baseline GAM.
# ------------------------------------------------------------
aps_by_run <- aps_by_run %>%
  left_join(
    base_aps,
    by = c("var", "variable", "var_label", "component")
  ) %>%
  mutate(
    aps_sign_same = case_when(
      !is.finite(aps_base) | !is.finite(aps_pert) ~ NA,
      aps_base == 0 & aps_pert == 0 ~ TRUE,
      TRUE ~ sign(aps_base) == sign(aps_pert)
    )
  )

aps_stability <- aps_by_run %>%
  group_by(var, variable, var_label, component) %>%
  summarise(
    aps_pert_median = median(aps_pert, na.rm = TRUE),
    aps_pert_q025 = quantile(aps_pert, .025, na.rm = TRUE),
    aps_pert_q975 = quantile(aps_pert, .975, na.rm = TRUE),
    pct_aps_same_sign = mean(aps_sign_same, na.rm = TRUE),
    n_runs_aps = n_distinct(run),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Spatial local-effect stability at z=0
# Keep legacy output column names (surface_cor etc.) for downstream
# compatibility, but record explicitly that the compared spatial object
# is the local marginal-slope pattern at z=0.
# ------------------------------------------------------------
spatial_stability <- spatial_by_run %>%
  group_by(var, variable, var_label, component) %>%
  summarise(
    surface_sd_base = first(surface_sd_base),
    surface_meaningful =
      is.finite(first(surface_sd_base)) &&
      first(surface_sd_base) >= SURFACE_SD_MIN,
    mean_surface_cor = mean(surface_cor, na.rm = TRUE),
    p10_surface_cor = quantile(surface_cor, .10, na.rm = TRUE),
    mean_surface_rank_cor = mean(surface_rank_cor, na.rm = TRUE),
    median_surface_nrmse = median(surface_nrmse, na.rm = TRUE),
    p90_surface_nrmse = quantile(surface_nrmse, .90, na.rm = TRUE),
    mean_county_sign_agreement = mean(county_sign_agreement, na.rm = TRUE),
    p10_county_sign_agreement = quantile(county_sign_agreement, .10, na.rm = TRUE),
    n_runs_spatial = n_distinct(run),
    spatial_reference_z = first(spatial_reference_z),
    .groups = "drop"
  ) %>%
  mutate(
    spatial_stable =
      !surface_meaningful |
      (
        p10_surface_cor >= SURFACE_R_THRESHOLD &
        mean_county_sign_agreement >= COUNTY_SIGN_THRESHOLD
      )
  )

# ------------------------------------------------------------
# SNR + final stability class
# ------------------------------------------------------------
snr_table <- base_grids %>%
  group_by(var, variable, var_label, component) %>%
  summarise(
    signal_range = diff(range(effect_base, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(
    pert_summary %>%
      group_by(var, variable, var_label, component) %>%
      summarise(
        noise_sd = mean(effect_sd, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("var", "variable", "var_label", "component")
  ) %>%
  left_join(
    base_grids %>%
      inner_join(
        pert_summary %>%
          select(var, variable, component, var_value, effect_median),
        by = c("var", "variable", "component", "var_value")
      ) %>%
      group_by(var, variable, component) %>%
      summarise(
        r_shape = suppressWarnings(
          cor(effect_base, effect_median, use = "complete.obs")
        ),
        .groups = "drop"
      ),
    by = c("var", "variable", "component")
  ) %>%
  left_join(
    shape_by_run %>%
      group_by(var, variable, component) %>%
      summarise(
        mean_r_shape = mean(r_run, na.rm = TRUE),
        sd_r_shape = sd(r_run, na.rm = TRUE),
        pct_positive_r = mean(r_run > 0, na.rm = TRUE),
        p10_r_shape = quantile(r_run, .10, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("var", "variable", "component")
  ) %>%
  left_join(
    aps_all %>%
      select(variable, component, aps, aps_q25, aps_q75),
    by = c("variable", "component")
  ) %>%
  left_join(
    aps_stability,
    by = c("var", "variable", "var_label", "component")
  ) %>%
  left_join(
    spatial_stability,
    by = c("var", "variable", "var_label", "component")
  ) %>%
  mutate(
    snr = signal_range / noise_sd,
    meaningful = abs(aps) >= MEANINGFUL_APS_THRESHOLD,

    shape_stable =
      r_shape >= SHAPE_R_THRESHOLD &
      mean_r_shape >= SHAPE_R_THRESHOLD &
      pct_positive_r >= .90,

    aps_sign_stable =
      pct_aps_same_sign >= APS_SIGN_THRESHOLD,

    final_class = case_when(
      !meaningful ~ "negligible_effect",

      snr > SNR_THRESHOLD &
        shape_stable &
        aps_sign_stable &
        spatial_stable ~ "perturbation_stable",

      TRUE ~ "unstable"
    ),

    period = PERIOD,
    year_label = YEAR_LABEL
  ) %>%
  arrange(component, desc(snr))

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
readr::write_csv(
  pert_summary,
  file.path(TAB_DIR, "perturbation_curve_summary.csv")
)
readr::write_csv(
  shape_by_run,
  file.path(TAB_DIR, "shape_by_run.csv")
)
readr::write_csv(
  aps_by_run,
  file.path(TAB_DIR, "aps_by_run.csv")
)
readr::write_csv(
  aps_stability,
  file.path(TAB_DIR, "aps_stability.csv")
)

# New explicit filenames:
readr::write_csv(
  spatial_by_run,
  file.path(TAB_DIR, "spatial_z0_stability_by_run.csv")
)
readr::write_csv(
  spatial_stability,
  file.path(TAB_DIR, "spatial_z0_stability.csv")
)

# Backward-compatible filenames for 05_figures.R / existing workflows:
readr::write_csv(
  spatial_by_run,
  file.path(TAB_DIR, "spatial_stability_by_run.csv")
)
readr::write_csv(
  spatial_stability,
  file.path(TAB_DIR, "spatial_stability.csv")
)

readr::write_csv(
  snr_table,
  file.path(TAB_DIR, "snr_table.csv")
)

noise_params <- tibble(
  period = PERIOD,
  year_label = YEAR_LABEL,
  sigma_resid = sigma_resid,
  rho_t = rho_t,
  lambda_s = lambda_s,
  noise_scale = NOISE_SCALE,
  B_requested = B_PERTURB,
  B_successful = n_success,
  mean_minutes_per_successful_run = mean(run_minutes[ok_run], na.rm = TRUE),
  total_perturbation_minutes = sum(run_minutes, na.rm = TRUE),
  bam_nthreads = BAM_NTHREADS,
  bam_gc_level = BAM_GC_LEVEL,
  analysis_grid_n = N_GRID,
  slope_eps = SLOPE_EPS,
  spatial_reference_z = PERT_SPATIAL_REFERENCE_Z,
  spatial_stability_definition = "county local marginal slope at common z=0",
  aps_definition = "county APS integrated across complete analysis grid; national median across counties",
  implementation = "fast_vectorized_predictions_v2",
  snr_threshold = SNR_THRESHOLD,
  shape_r_threshold = SHAPE_R_THRESHOLD,
  surface_r_threshold = SURFACE_R_THRESHOLD,
  county_sign_threshold = COUNTY_SIGN_THRESHOLD,
  aps_sign_threshold = APS_SIGN_THRESHOLD,
  surface_sd_min = SURFACE_SD_MIN,
  meaningful_aps_threshold = MEANINGFUL_APS_THRESHOLD
)

readr::write_csv(
  noise_params,
  file.path(TAB_DIR, "perturbation_settings.csv")
)

# Compact diagnostic plot.
p <- snr_table %>%
  mutate(
    label = forcats::fct_reorder(
      paste0(var_label, " (", component, ")"),
      snr
    )
  ) %>%
  ggplot(aes(snr, label, shape = component)) +
  geom_vline(
    xintercept = SNR_THRESHOLD,
    linetype = "dashed",
    colour = LANCET$grey_mid
  ) +
  geom_point(size = 2.2, colour = LANCET$navy) +
  labs(
    title = paste0("Perturbation stability | ", YEAR_LABEL),
    subtitle = paste0(
      "Structured residual-noise stress test; curve + APS + ",
      "county local-slope stability at z=0"
    ),
    x = "Signal-to-noise ratio",
    y = NULL
  ) +
  theme_lancet(base_size = 9)

save_hp_plot(
  file.path(FIG_DIR, "snr_by_term.png"),
  p,
  FIG_W_FULL,
  max(4.8, nrow(snr_table) * .22)
)

saveRDS(
  list(
    curves_pert = curves_pert,
    pert_summary = pert_summary,
    shape_by_run = shape_by_run,
    aps_by_run = aps_by_run,
    aps_stability = aps_stability,
    base_spatial_z0 = base_spatial_z0,
    spatial_by_run = spatial_by_run,
    spatial_stability = spatial_stability,
    snr_table = snr_table,
    noise_params = noise_params,
    run_minutes = run_minutes,
    spatial_reference_z = PERT_SPATIAL_REFERENCE_Z,
    implementation = "fast_vectorized_predictions_v2",
    period = PERIOD,
    year_label = YEAR_LABEL
  ),
  file.path(OUT_DIR, "perturbation_results.rds")
)

cat("\nPerturbation DONE:", OUT_DIR, "\n")
cat(sprintf(
  "Successful runs: %d/%d | mean %.2f min/run | total %.1f min\n",
  n_success,
  B_PERTURB,
  mean(run_minutes[ok_run], na.rm = TRUE),
  sum(run_minutes, na.rm = TRUE)
))

print(
  snr_table %>%
    select(
      var_label, component, aps, snr,
      r_shape, mean_r_shape,
      pct_aps_same_sign,
      p10_surface_cor,
      mean_county_sign_agreement,
      final_class
    ),
  n = 50
)
