# ============================================================
# 03_validation.R — diagnostyka + strict leave-one-year-out CV
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("strict_fold_data")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(mgcv); library(ggplot2); library(readr)
})

GAM_RDS <- file.path(OUT_BASE, "spatial_mundlak_gam", "gam_results.rds")
if (!file.exists(GAM_RDS)) stop("Brak: ", GAM_RDS)
res <- readRDS(GAM_RDS)
m <- res$m
dat_model <- res$dat_model
between_terms <- res$between_terms
within_terms <- res$within_terms
all_terms <- res$all_terms

OUT_DIR <- file.path(OUT_BASE, "validation")
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== VALIDATION ===\n")
cat(MODEL_LABEL, "|", YEAR_LABEL, "\n")

if (!exists("powiaty_spdf") || is.null(powiaty_spdf)) {
  stop("Finalna walidacja wymaga powiaty_spdf do Moran's I.")
}
if (!requireNamespace("spdep", quietly = TRUE)) {
  stop("Finalna walidacja wymaga pakietu 'spdep'.")
}

rmse_fn <- function(e) sqrt(mean(e^2, na.rm = TRUE))
mae_fn  <- function(e) mean(abs(e), na.rm = TRUE)

fitted_v <- as.numeric(predict(m, newdata = dat_model, type = "response"))
obs_v <- as.numeric(dat_model[[OUTCOME]])
resid_v <- obs_v - fitted_v
resid_sd <- sd(resid_v, na.rm = TRUE)
if (!is.finite(resid_sd) || resid_sd == 0) resid_sd <- 1

dat_val <- dat_model %>%
  mutate(observed = obs_v,
         fitted = fitted_v,
         resid = resid_v,
         resid_std = resid_v / resid_sd)

val_global <- tibble(
  period = PERIOD,
  year_label = YEAR_LABEL,
  n = nrow(dat_val),
  n_counties = n_distinct(dat_val$code),
  r2_insample = 1 - sum((obs_v - fitted_v)^2) / sum((obs_v - mean(obs_v))^2),
  rmse_insample = rmse_fn(resid_v),
  mae_insample = mae_fn(resid_v),
  bias_insample = mean(fitted_v - obs_v)
)
readr::write_csv(val_global, file.path(TAB_DIR, "validation_global.csv"))

val_by_year <- dat_val %>%
  group_by(year) %>%
  summarise(
    n = n(),
    rmse = rmse_fn(resid),
    mae = mae_fn(resid),
    bias = mean(fitted - observed, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(period = PERIOD)
readr::write_csv(val_by_year, file.path(TAB_DIR, "validation_by_year.csv"))

# ------------------------------------------------------------
# Temporal residual ACF1 within counties
# ------------------------------------------------------------
temporal_acf <- dat_val %>%
  arrange(code, year) %>%
  group_by(code) %>%
  summarise(
    n = n(),
    acf1 = if (n() >= 3L) suppressWarnings(cor(resid[-n()], resid[-1], use = "complete.obs")) else NA_real_,
    .groups = "drop"
  )
acf_summary <- temporal_acf %>%
  summarise(
    median_acf1 = median(acf1, na.rm = TRUE),
    mean_acf1 = mean(acf1, na.rm = TRUE),
    pct_abs_gt_02 = 100 * mean(abs(acf1) > .2, na.rm = TRUE)
  )
readr::write_csv(temporal_acf, file.path(TAB_DIR, "temporal_acf1_by_county.csv"))
readr::write_csv(acf_summary, file.path(TAB_DIR, "temporal_acf1_summary.csv"))

# ------------------------------------------------------------
# Moran's I by year
# ------------------------------------------------------------
moran_by_year <- tibble()
moran_summary <- tibble()
if (!is.null(powiaty_spdf) && requireNamespace("spdep", quietly = TRUE)) {
  pow <- sf::st_as_sf(powiaty_spdf)
  code_col <- intersect(c(POW_CODE_COL, "JPT_KOD_JE", "code_pow"), names(pow))[1]
  pow[[code_col]] <- as.character(pow[[code_col]])
  pow[[code_col]] <- ifelse(nchar(pow[[code_col]]) == 3L, paste0("0", pow[[code_col]]), pow[[code_col]])
  pow <- pow %>% arrange(.data[[code_col]])
  nb <- spdep::poly2nb(pow, queen = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  moran_by_year <- map_dfr(sort(unique(dat_val$year)), function(yr) {
    dd <- dat_val %>% filter(year == yr)
    idx <- match(pow[[code_col]], dd$code_pow)
    rr <- dd$resid[idx]
    ok <- is.finite(rr)
    if (sum(ok) < 10L) {
      return(tibble(year = yr, moran_i = NA_real_,
                    p_greater = NA_real_, p_two_sided = NA_real_))
    }
    mt_g <- tryCatch(
      spdep::moran.test(rr, lw, zero.policy = TRUE, na.action = na.exclude,
                        alternative = "greater"),
      error = function(e) NULL
    )
    mt_2 <- tryCatch(
      spdep::moran.test(rr, lw, zero.policy = TRUE, na.action = na.exclude,
                        alternative = "two.sided"),
      error = function(e) NULL
    )
    if (is.null(mt_g)) {
      return(tibble(year = yr, moran_i = NA_real_,
                    p_greater = NA_real_, p_two_sided = NA_real_))
    }
    tibble(
      year = yr,
      moran_i = unname(mt_g$estimate[["Moran I statistic"]]),
      p_greater = mt_g$p.value,
      p_two_sided = if (is.null(mt_2)) NA_real_ else mt_2$p.value
    )
  }) %>%
    mutate(
      sig_positive = is.finite(moran_i) & moran_i > 0 & p_greater < .05,
      sig_two_sided = p_two_sided < .05,
      period = PERIOD
    )
  readr::write_csv(moran_by_year, file.path(TAB_DIR, "moran_by_year.csv"))

  moran_summary <- moran_by_year %>%
    summarise(
      period = PERIOD,
      year_label = YEAR_LABEL,
      n_years = sum(is.finite(moran_i)),
      median_moran_i = median(moran_i, na.rm = TRUE),
      max_abs_moran_i = max(abs(moran_i), na.rm = TRUE),
      n_positive_p05 = sum(sig_positive, na.rm = TRUE),
      n_two_sided_p05 = sum(sig_two_sided, na.rm = TRUE)
    )
  readr::write_csv(moran_summary, file.path(TAB_DIR, "moran_summary.csv"))
} else {
  warning("spdep lub geometria niedostępna — Moran's I pominięty.")
}

# ------------------------------------------------------------
# k-check / concurvity
# ------------------------------------------------------------
k_check <- tryCatch(as.data.frame(mgcv::k.check(m), stringsAsFactors = FALSE),
                    error = function(e) data.frame(error = e$message))
if (!"error" %in% names(k_check)) k_check <- tibble::rownames_to_column(k_check, "term")
readr::write_csv(k_check, file.path(TAB_DIR, "k_check.csv"))

concurv <- tryCatch(mgcv::concurvity(m, full = TRUE), error = function(e) NULL)
if (!is.null(concurv)) {
  concurv_df <- as.data.frame(concurv) %>% tibble::rownames_to_column("metric")
  readr::write_csv(concurv_df, file.path(TAB_DIR, "concurvity.csv"))
}

# ------------------------------------------------------------
# STRICT LOO-year CV
# Preprocessing jest fitowane WYŁĄCZNIE na training years:
# county means, scaling oraz allowed lag source years.
# Dodatkowo obserwacja t+1 w treningu wypada, jeśli jej lag pochodzi
# z held-out year t. To zapobiega transductive leakage.
# ------------------------------------------------------------
loo_results <- tibble()
if (isTRUE(RUN_STRICT_LOO)) {
  cv_years <- seq.int(YEAR_MIN + WITHIN_LAG_N, YEAR_MAX)
  cat("Strict LOO-year CV:", paste(cv_years, collapse = ", "), "\n")

  loo_results <- map_dfr(cv_years, function(yr) {
    cat("  held-out:", yr, "\n")
    fd <- tryCatch(strict_fold_data(dat, yr, PERIOD), error = function(e) NULL)
    if (is.null(fd) || nrow(fd$train) == 0L || nrow(fd$test) == 0L) {
      return(tibble(year = yr, n_train = if (is.null(fd)) NA_integer_ else nrow(fd$train),
                    n_test = if (is.null(fd)) NA_integer_ else nrow(fd$test),
                    k_year_cv = NA_integer_,
                    rmse = NA_real_, mae = NA_real_, r2 = NA_real_, bias = NA_real_))
    }

    # Held-out folds, especially M1, contain fewer unique outcome years than
    # the final full-sample model. Thin-plate k cannot exceed the number of
    # unique time points, so cap only the CV temporal basis at the largest
    # feasible value. All other basis dimensions remain identical.
    spec_cv <- GAM_SPEC
    n_time_cv <- dplyr::n_distinct(fd$train$year_c[is.finite(fd$train$year_c)])
    spec_cv$k_year <- min(as.integer(GAM_SPEC$k_year), as.integer(n_time_cv))
    if (!is.finite(spec_cv$k_year) || spec_cv$k_year < 3L) {
      return(tibble(year = yr, n_train = nrow(fd$train), n_test = nrow(fd$test),
                    k_year_cv = as.integer(spec_cv$k_year),
                    rmse = NA_real_, mae = NA_real_, r2 = NA_real_, bias = NA_real_))
    }

    fml_cv <- build_gam_formula(
      fd$prep$between_terms, fd$prep$within_terms, CONTROL_VARS,
      spec = spec_cv
    )
    m_cv <- tryCatch(
      mgcv::bam(fml_cv, data = fd$train, family = FAMILY, method = BAM_METHOD,
                discrete = BAM_DISCRETE, gamma = GAMMA, nthreads = BAM_NTHREADS,
                gc.level = BAM_GC_LEVEL,
                na.action = na.exclude, control = gam.control(trace = FALSE)),
      error = function(e) { warning("CV ", yr, ": ", e$message); NULL }
    )
    if (is.null(m_cv)) {
      return(tibble(year = yr, n_train = nrow(fd$train), n_test = nrow(fd$test),
                    k_year_cv = as.integer(spec_cv$k_year),
                    rmse = NA_real_, mae = NA_real_, r2 = NA_real_, bias = NA_real_))
    }
    pred <- tryCatch(as.numeric(predict(m_cv, newdata = fd$test, type = "response")),
                     error = function(e) rep(NA_real_, nrow(fd$test)))
    obs <- as.numeric(fd$test[[OUTCOME]])
    ok <- is.finite(obs) & is.finite(pred)
    obs <- obs[ok]; pred <- pred[ok]
    denom <- sum((obs - mean(obs))^2)
    tibble(
      year = yr,
      n_train = nrow(fd$train),
      n_test = length(obs),
      k_year_cv = as.integer(spec_cv$k_year),
      rmse = rmse_fn(obs - pred),
      mae = mae_fn(obs - pred),
      r2 = if (is.finite(denom) && denom > 0) 1 - sum((obs - pred)^2) / denom else NA_real_,
      bias = mean(pred - obs)
    )
  }) %>% mutate(period = PERIOD)

  readr::write_csv(loo_results, file.path(TAB_DIR, "strict_loo_year_cv.csv"))
}

loo_summary <- if (nrow(loo_results)) {
  loo_results %>% summarise(
    period = PERIOD,
    year_label = YEAR_LABEL,
    n_years = sum(is.finite(rmse)),
    k_year_cv_min = if (any(is.finite(k_year_cv))) min(k_year_cv[is.finite(k_year_cv)]) else NA_integer_,
    k_year_cv_max = if (any(is.finite(k_year_cv))) max(k_year_cv[is.finite(k_year_cv)]) else NA_integer_,
    rmse_mean = mean(rmse, na.rm = TRUE),
    rmse_sd = sd(rmse, na.rm = TRUE),
    mae_mean = mean(mae, na.rm = TRUE),
    r2_mean = mean(r2, na.rm = TRUE),
    r2_median = median(r2, na.rm = TRUE),
    bias_mean = mean(bias, na.rm = TRUE)
  )
} else tibble(period = PERIOD, year_label = YEAR_LABEL)
readr::write_csv(loo_summary, file.path(TAB_DIR, "strict_loo_summary.csv"))

# ------------------------------------------------------------
# Diagnostyczne figury
# ------------------------------------------------------------
p1 <- ggplot(dat_val, aes(observed, fitted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = LANCET$grey_mid) +
  geom_point(alpha = .25, size = .55, colour = LANCET$navy) +
  coord_equal() + labs(title = "Observed vs fitted", x = "Observed", y = "Fitted") +
  theme_lancet(base_size = 9)
save_hp_plot(file.path(FIG_DIR, "observed_vs_fitted.png"), p1, FIG_W_MID, FIG_W_MID * .72)

if (nrow(loo_results)) {
  p2 <- loo_results %>%
    select(year, RMSE = rmse, `R²` = r2) %>%
    pivot_longer(c(RMSE, `R²`), names_to = "metric", values_to = "value") %>%
    ggplot(aes(year, value)) + geom_line() + geom_point(size = 1.6) +
    facet_wrap(~metric, scales = "free_y") +
    labs(title = "Strict leave-one-year-out validation", x = "Held-out year", y = NULL) +
    theme_lancet(base_size = 9)
  save_hp_plot(file.path(FIG_DIR, "strict_loo_year_cv.png"), p2, FIG_W_FULL, FIG_W_FULL * .45)
}

saveRDS(list(
  val_global = val_global,
  val_by_year = val_by_year,
  temporal_acf = temporal_acf,
  acf_summary = acf_summary,
  moran_by_year = moran_by_year,
  moran_summary = moran_summary,
  k_check = k_check,
  concurvity = concurv,
  loo_results = loo_results,
  loo_summary = loo_summary,
  period = PERIOD,
  year_label = YEAR_LABEL
), file.path(OUT_DIR, "validation_results.rds"))

cat("\nValidation DONE:", OUT_DIR, "\n")
print(val_global)
print(loo_summary)
