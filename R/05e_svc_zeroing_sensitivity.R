# SSM v6 namespace patch: dplyr verbs explicitly qualified to prevent select()/filter() masking.
# ============================================================
# 05e_svc_zeroing_sensitivity.R
# Post-estimation SVC zeroing sensitivity check — NO MODEL REFIT.
#
# Compares the reported focal-determinant effect (global smooth + focal SVC)
# with the SAME fitted GAM evaluated after suppressing ONLY the focal SVC
# coefficient block in the lpmatrix contrast.
#
# This is not a refit without SVC. It answers a narrower reporting question:
# how much do fitted conclusions change when the local spatial modifier of the
# focal determinant is removed while the fitted global smooth is retained?
#
# Outputs per model:
#   spatial_mundlak_gam/tables/svc_zeroing_sensitivity.csv
#   spatial_mundlak_gam/tables/svc_zeroing_curve_diagnostics.csv
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("make_reference_counties")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

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
if (!file.exists(RDS_PATH)) stop("Brak gotowego GAM-u: ", RDS_PATH)

g <- readRDS(RDS_PATH)
m <- g$m
dat_model <- g$dat_model
between_terms <- g$between_terms
within_terms <- g$within_terms
all_terms <- g$all_terms
controls <- g$settings$CONTROL_VARS
beta <- stats::coef(m)

ref_counties <- make_reference_counties(dat_model, all_terms, controls)

# Identify coefficient columns belonging specifically to s(x,y, by=focal_var).
# mgcv stores this robustly in each smooth object's $by and parameter range.
svc_coef_index <- function(mod, var) {
  idx <- integer(0)
  matched_labels <- character(0)
  for (sm in mod$smooth) {
    by_nm <- if (!is.null(sm$by)) as.character(sm$by) else "NA"
    lab <- if (!is.null(sm$label)) as.character(sm$label) else ""
    is_xy <- grepl("s\\(x,y", gsub(" ", "", lab), fixed = FALSE)
    if (identical(by_nm, var) && is_xy) {
      idx <- c(idx, seq.int(sm$first.para, sm$last.para))
      matched_labels <- c(matched_labels, lab)
    }
  }
  list(index = sort(unique(idx)), labels = unique(matched_labels))
}

svc_edf_for <- function(mod, labels) {
  if (!length(labels)) return(0)
  st <- tryCatch(summary(mod)$s.table, error = function(e) NULL)
  if (is.null(st) || !nrow(st)) return(NA_real_)
  rr <- rownames(st)
  ii <- unique(unlist(lapply(labels, function(lbl) {
    which(rr == lbl | grepl(lbl, rr, fixed = TRUE))
  })))
  if (!length(ii)) return(NA_real_)
  sum(as.numeric(st[ii, "edf"]), na.rm = TRUE)
}

predict_contrast <- function(nd0, nd1, zero_idx = integer(0), scale = 1) {
  X0 <- predict(m, newdata = nd0, type = "lpmatrix")
  X1 <- predict(m, newdata = nd1, type = "lpmatrix")
  L <- (X1 - X0) / scale
  total <- as.numeric(L %*% beta)
  if (!length(zero_idx)) return(list(total = total, zeroed = total))
  Lz <- L
  Lz[, zero_idx] <- 0
  zeroed <- as.numeric(Lz %*% beta)
  list(total = total, zeroed = zeroed)
}

safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L || stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok], method = "pearson"))
}

sign_agreement <- function(a, b, tol = 1e-8) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  sa <- ifelse(abs(a[ok]) <= tol, 0, sign(a[ok]))
  sb <- ifelse(abs(b[ok]) <= tol, 0, sign(b[ok]))
  mean(sa == sb)
}

summary_map_pair <- function(total, zeroed, prefix) {
  tibble(
    !!paste0(prefix, "_pearson_r") := safe_cor(total, zeroed),
    !!paste0(prefix, "_sign_agreement_pct") := 100 * sign_agreement(total, zeroed),
    !!paste0(prefix, "_median_abs_difference") := median(abs(total - zeroed), na.rm = TRUE),
    !!paste0(prefix, "_p90_abs_difference") := as.numeric(quantile(abs(total - zeroed), .90, na.rm = TRUE, names = FALSE))
  )
}

compute_one <- function(var) {
  comp <- comp_from_term(var)
  si <- svc_coef_index(m, var)
  zero_idx <- si$index
  svc_edf <- svc_edf_for(m, si$labels)

  vals <- dat_model[[var]]
  q <- as.numeric(stats::quantile(vals, c(.02, .98), na.rm = TRUE, names = FALSE))
  grid <- seq(q[1], q[2], length.out = 25L)

  # ----- Effect curve and county APS -----
  curve_total <- matrix(NA_real_, nrow(ref_counties), length(grid))
  curve_zero <- matrix(NA_real_, nrow(ref_counties), length(grid))
  slope_total <- matrix(NA_real_, nrow(ref_counties), length(grid))
  slope_zero <- matrix(NA_real_, nrow(ref_counties), length(grid))

  nd_ref <- ref_counties
  nd_ref[[var]] <- 0

  for (kk in seq_along(grid)) {
    vv <- grid[kk]
    nd <- ref_counties
    nd[[var]] <- vv
    cc <- predict_contrast(nd_ref, nd, zero_idx = zero_idx)
    curve_total[, kk] <- cc$total
    curve_zero[, kk] <- cc$zeroed

    nd0 <- ref_counties
    nd1 <- ref_counties
    nd0[[var]] <- vv
    nd1[[var]] <- vv + SLOPE_EPS
    ss <- predict_contrast(nd0, nd1, zero_idx = zero_idx, scale = SLOPE_EPS)
    slope_total[, kk] <- ss$total
    slope_zero[, kk] <- ss$zeroed
  }

  aps_county_total <- rowMeans(slope_total, na.rm = TRUE)
  aps_county_zero <- rowMeans(slope_zero, na.rm = TRUE)

  med_curve_total <- apply(curve_total, 2, median, na.rm = TRUE)
  med_curve_zero <- apply(curve_zero, 2, median, na.rm = TRUE)
  iqr_width_total <- apply(curve_total, 2, function(z) diff(quantile(z, c(.25, .75), na.rm = TRUE)))
  p1090_width_total <- apply(curve_total, 2, function(z) diff(quantile(z, c(.10, .90), na.rm = TRUE)))
  iqr_width_zero <- apply(curve_zero, 2, function(z) diff(quantile(z, c(.25, .75), na.rm = TRUE)))

  # ----- Marginal and accumulated maps -----
  if (comp == "between") {
    v_actual <- dat_model %>%
      dplyr::select(code, dplyr::all_of(var)) %>%
      dplyr::distinct(code, .keep_all = TRUE) %>%
      dplyr::rename(v_actual = dplyr::all_of(var))
    base <- ref_counties %>% dplyr::inner_join(v_actual, by = "code")

    nd0 <- base; nd1 <- base
    nd0[[var]] <- base$v_actual
    nd1[[var]] <- base$v_actual + SLOPE_EPS
    mm <- predict_contrast(nd0, nd1, zero_idx = zero_idx, scale = SLOPE_EPS)

    nr <- base; na <- base
    nr[[var]] <- 0
    na[[var]] <- base$v_actual
    aa <- predict_contrast(nr, na, zero_idx = zero_idx)

    marginal_total <- mm$total
    marginal_zero <- mm$zeroed
    accum_total <- aa$total
    accum_zero <- aa$zeroed
  } else {
    meta <- dat_model %>%
      dplyr::select(code, code_pow, name, year, x, y, dplyr::all_of(var)) %>%
      dplyr::rename(v_actual = dplyr::all_of(var))
    nd0 <- meta
    nd1 <- meta
    nd0$year_c <- 0; nd1$year_c <- 0
    for (cc in controls) { nd0[[cc]] <- 0; nd1[[cc]] <- 0 }
    for (tt in all_terms) { nd0[[tt]] <- 0; nd1[[tt]] <- 0 }
    nd0[[var]] <- meta$v_actual
    nd1[[var]] <- meta$v_actual + SLOPE_EPS
    mm_it <- predict_contrast(nd0, nd1, zero_idx = zero_idx, scale = SLOPE_EPS)

    mm_df <- tibble(code = meta$code, total = mm_it$total, zeroed = mm_it$zeroed) %>%
      dplyr::group_by(code) %>%
      dplyr::summarise(total = mean(total, na.rm = TRUE), zeroed = mean(zeroed, na.rm = TRUE), .groups = "drop")
    marginal_total <- mm_df$total
    marginal_zero <- mm_df$zeroed

    nr <- ref_counties; na <- ref_counties
    nr[[var]] <- 0
    na[[var]] <- 1
    aa <- predict_contrast(nr, na, zero_idx = zero_idx)
    accum_total <- aa$total
    accum_zero <- aa$zeroed
  }

  aps_total <- median(aps_county_total, na.rm = TRUE)
  aps_zero <- median(aps_county_zero, na.rm = TRUE)

  dplyr::bind_cols(
    tibble(
      period = PERIOD,
      year_label = YEAR_LABEL,
      variable = clean_var(var),
      var = var,
      var_label = label_var(clean_var(var)),
      component = comp,
      svc_term_found = length(zero_idx) > 0,
      svc_smooth_label = paste(si$labels, collapse = "; "),
      svc_edf = svc_edf,
      aps_current = aps_total,
      aps_svc_zeroed = aps_zero,
      delta_aps_current_minus_zeroed = aps_total - aps_zero,
      abs_delta_aps = abs(aps_total - aps_zero),
      aps_sign_changed = ifelse(abs(aps_total) <= 1e-8 || abs(aps_zero) <= 1e-8,
                                sign(aps_total) != sign(aps_zero),
                                sign(aps_total) != sign(aps_zero)),
      county_aps_sign_agreement_pct = 100 * sign_agreement(aps_county_total, aps_county_zero),
      median_curve_pearson_r = safe_cor(med_curve_total, med_curve_zero),
      current_curve_median_IQR_width = median(iqr_width_total, na.rm = TRUE),
      current_curve_median_P10_P90_width = median(p1090_width_total, na.rm = TRUE),
      svc_zeroed_curve_median_IQR_width = median(iqr_width_zero, na.rm = TRUE),
      ribbon_interpretation = ifelse(
        length(zero_idx) == 0 || (!is.na(svc_edf) && svc_edf < 0.05),
        "SVC absent/effectively penalised to zero; little spatial effect heterogeneity expected",
        "Non-zero ribbon width reflects cross-county heterogeneity induced by the focal SVC term"
      )
    ),
    summary_map_pair(marginal_total, marginal_zero, "marginal_map"),
    summary_map_pair(accum_total, accum_zero, "accumulated_map")
  )
}

cat("=== SVC ZEROING SENSITIVITY | ", PERIOD, " ===\n", sep = "")
res <- purrr::map_dfr(all_terms, function(v) {
  tryCatch(compute_one(v), error = function(e) {
    warning("SVC zeroing failed for ", v, ": ", e$message)
    tibble(
      period = PERIOD, year_label = YEAR_LABEL,
      variable = clean_var(v), var = v,
      var_label = label_var(clean_var(v)), component = comp_from_term(v),
      error = e$message
    )
  })
})

readr::write_csv(res, file.path(TAB_DIR, "svc_zeroing_sensitivity.csv"))

curve_diag <- res %>%
  dplyr::select(dplyr::any_of(c(
    "period", "year_label", "variable", "var", "var_label", "component",
    "svc_term_found", "svc_smooth_label", "svc_edf",
    "current_curve_median_IQR_width", "current_curve_median_P10_P90_width",
    "svc_zeroed_curve_median_IQR_width", "median_curve_pearson_r",
    "ribbon_interpretation"
  )))
readr::write_csv(curve_diag, file.path(TAB_DIR, "svc_zeroing_curve_diagnostics.csv"))

cat("SVC zeroing sensitivity written to: ", TAB_DIR, "\n", sep = "")
print(res %>% dplyr::select(dplyr::any_of(c(
  "var_label", "component", "svc_edf", "aps_current", "aps_svc_zeroed",
  "delta_aps_current_minus_zeroed", "aps_sign_changed",
  "marginal_map_sign_agreement_pct", "accumulated_map_sign_agreement_pct"
))), n = 100)
