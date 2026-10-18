# ============================================================
# 08_baseline_boundary_check_175_200.R
# ============================================================
# Ostatni boundary check dla K_BASELINE na modelu głównym M1.
#
# Stałe:
#   K_GLOBAL   = 5
#   K_SVC      = 25
#   K_YEAR     = 9
#
# Porównanie:
#   K_BASELINE = 175 vs 200
#
# Cel:
# - sprawdzić, czy zwiększenie baseline spatial basis 175 -> 200
#   przestaje substantywnie zmieniać APS, krzywe i powierzchnie SVC;
# - sprawdzić EDF/k', k.check, AIC/BIC/RMSE, Moran's I i ACF(1);
# - BEZ perturbacji.
#
# Dwa modele są celowo estymowane w tym samym skrypcie, nawet jeśli
# model BASELINE=175 był już liczony w gridzie 48. Dzięki temu oba
# fit'y i signature effects są wyliczone identycznymi ustawieniami,
# co daje czyste, samodzielne porównanie do raportowania.
# ============================================================

# -------------------------
# USTAWIENIA
# -------------------------
KGRID_PERIOD       <- "m1_2010_2019"
KGRID_ARCHITECTURE <- "C_all_svc"

KGRID_GLOBAL_VALUES   <- 5L
KGRID_SVC_VALUES      <- 25L
KGRID_BASELINE_VALUES <- c(175L, 200L)
KGRID_YEAR_VALUES     <- 9L

# Dla finalnego boundary check zwiększamy rozdzielczość diagnostyczną.
KGRID_KCHECK_NREP   <- 200L
KGRID_EFFECT_GRID_N <- 15L
KGRID_SLOPE_EPS     <- 0.01

# Robocze progi STABILNOŚCI — heurystyki opisowe, nie formalne testy.
KGRID_STAB_APS_COR_MIN       <- 0.99
KGRID_STAB_APS_SIGN_MIN      <- 1.00
KGRID_STAB_CURVE_COR_MIN     <- 0.97
KGRID_STAB_SURFACE_COR_MIN   <- 0.95
KGRID_MEANINGFUL_APS         <- 1.0
KGRID_SURFACE_SD_MIN         <- 0.25
KGRID_CURVE_RANGE_MIN        <- 0.50

KGRID_RESUME      <- TRUE
KGRID_SAVE_MODELS <- FALSE

# -------------------------
# ŚCIEŻKI / PAKIETY
# -------------------------
PROJECT_DIR <- get0(
  "PROJECT_DIR",
  ifnotfound = normalizePath(getwd(), winslash = "/", mustWork = TRUE)
)
SCRIPTS_DIR <- get0("SCRIPTS_DIR", ifnotfound = file.path(PROJECT_DIR, "R"))
DATA_DIR    <- get0("DATA_DIR",    ifnotfound = file.path(PROJECT_DIR, "data"))
PANEL_RDS   <- get0("PANEL_RDS",   ifnotfound = file.path(DATA_DIR, "finalny_panel_przestrzenny_379_powiatow.rds"))
PANEL_CSV   <- get0("PANEL_CSV",   ifnotfound = file.path(DATA_DIR, "finalny_panel_przestrzenny_379_powiatow.csv.gz"))
POWIAT_SHP_PATH <- get0("POWIAT_SHP_PATH", ifnotfound = NULL)

if (!dir.exists(PROJECT_DIR)) stop("Nie istnieje PROJECT_DIR: ", PROJECT_DIR)
if (!file.exists(file.path(SCRIPTS_DIR, "00_config.R"))) stop("Brak R/00_config.R")
if (!file.exists(file.path(SCRIPTS_DIR, "00_helpers.R"))) stop("Brak R/00_helpers.R")
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(mgcv)
  library(sf)
  library(ggplot2)
})

# Zapobiega załadowaniu starego selected_gam_spec.rds.
RUNNING_GAM_SELECTION <- TRUE
PERIOD <- KGRID_PERIOD
source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

# ------------------------------------------------------------
# COMPATIBILITY / STANDALONE DEFAULTS
# ------------------------------------------------------------
# Ten skrypt ma działać także wtedy, gdy lokalny 00_config.R pochodzi
# z v40.2/v40.3 i nie definiuje jeszcze wszystkich nazw z v40.4/v40.5.
# Nie nadpisujemy wartości istniejących w configu.
if (!exists("OUTPUT_ROOT", inherits = TRUE)) OUTPUT_ROOT <- "output_v40_5_final"
if (!exists("SCALING_MODE", inherits = TRUE)) SCALING_MODE <- "m1_reference"
if (!exists("GAMMA", inherits = TRUE)) GAMMA <- 1.4
if (!exists("BAM_METHOD", inherits = TRUE)) BAM_METHOD <- "fREML"
if (!exists("BAM_DISCRETE", inherits = TRUE)) BAM_DISCRETE <- TRUE
if (!exists("BAM_NTHREADS", inherits = TRUE)) {
  BAM_NTHREADS <- min(2L, max(1L, parallel::detectCores() - 1L))
}
if (!exists("BAM_GC_LEVEL", inherits = TRUE)) BAM_GC_LEVEL <- 2L
if (!exists("FAMILY", inherits = TRUE)) FAMILY <- gaussian()

.required_helpers <- c(
  "normalize_panel_columns", "resolve_outcome", "to_numeric_safe",
  "prepare_spatial_objects", "validate_raw_panel", "build_mundlak_panel",
  "model_complete_cases", "make_reference_counties",
  "clean_var", "label_var", "comp_from_term"
)
.missing_helpers <- .required_helpers[!vapply(.required_helpers, exists, logical(1),
                                               mode = "function", inherits = TRUE)]
if (length(.missing_helpers)) {
  stop(
    "Twój R/00_helpers.R jest starszy niż skrypt 08_baseline_boundary_check_175_200.R. ",
    "Brakuje funkcji: ", paste(.missing_helpers, collapse = ", "), ". ",
    "Podmień 00_helpers.R i 00_config.R na wersję z v40.5 albo użyj pełnej paczki v40.5."
  )
}

if (!KGRID_ARCHITECTURE %in% c("A_global_only", "B_between_svc", "C_all_svc")) {
  stop("Nieznana KGRID_ARCHITECTURE: ", KGRID_ARCHITECTURE)
}
if (length(KGRID_GLOBAL_VALUES) != 1L ||
    length(KGRID_SVC_VALUES) != 1L ||
    length(KGRID_BASELINE_VALUES) != 2L ||
    length(KGRID_YEAR_VALUES) != 1L) {
  stop("Ten skrypt oczekuje dokładnie: GLOBAL=5, SVC=25, BASELINE=175/200, YEAR=9.")
}

OUT_DIR <- file.path(OUTPUT_ROOT, "k_baseline_boundary_175_200", KGRID_PERIOD)
CHK_DIR <- file.path(OUT_DIR, "checkpoints")
MOD_DIR <- file.path(OUT_DIR, "models_optional")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHK_DIR, recursive = TRUE, showWarnings = FALSE)
if (KGRID_SAVE_MODELS) dir.create(MOD_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n=== K_BASELINE BOUNDARY CHECK 175 -> 200 | NO PERTURBATION ===\n")
cat("Period:", KGRID_PERIOD, "|", YEAR_LABEL, "\n")
cat("Architecture:", KGRID_ARCHITECTURE, "\n")
cat("K_GLOBAL:", paste(KGRID_GLOBAL_VALUES, collapse = ", "), "\n")
cat("K_SVC:", paste(KGRID_SVC_VALUES, collapse = ", "), "\n")
cat("K_BASELINE:", paste(KGRID_BASELINE_VALUES, collapse = ", "), "\n")
cat("K_YEAR:", paste(KGRID_YEAR_VALUES, collapse = ", "), "\n")
cat("Total fits:", length(KGRID_BASELINE_VALUES), "\n")
cat("Perturbation: OFF (dopiero po wyborze wystarczających baz)\n\n")

# -------------------------
# WCZYTANIE PANELU
# -------------------------
if (file.exists(PANEL_RDS)) {
  dat0 <- readRDS(PANEL_RDS)
} else if (file.exists(PANEL_CSV)) {
  dat0 <- readr::read_csv(PANEL_CSV, show_col_types = FALSE)
} else {
  stop("Nie znaleziono panelu RDS/CSV.")
}

dat0 <- normalize_panel_columns(dat0)
OUTCOME <- resolve_outcome(dat0, OUTCOME)
assign("OUTCOME", OUTCOME, envir = .GlobalEnv)
for (nm in intersect(c(OUTCOME, ALL_RAW_VARS), names(dat0))) dat0[[nm]] <- to_numeric_safe(dat0[[nm]])

powiaty0 <- NULL
if (!inherits(dat0, "sf")) {
  shp <- POWIAT_SHP_PATH
  if (is.null(shp) || !file.exists(shp)) {
    hits <- list.files(PROJECT_DIR, pattern = "powiaty\\.shp$", recursive = TRUE,
                       full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) shp <- hits[1]
  }
  if (!is.null(shp) && file.exists(shp)) powiaty0 <- sf::st_read(shp, quiet = TRUE)
}
sp_obj <- prepare_spatial_objects(dat0, powiaty0)
dat <- sp_obj$dat
powiaty_spdf <- sp_obj$powiaty_spdf
rm(dat0, powiaty0, sp_obj)
validate_raw_panel(dat, OUTCOME)
invisible(gc(full = TRUE))

# Compatibility: v40.4/v40.5 helper accepts `scaling_mode`, while v40.2/v40.3 does not.
# For the selection grid we work on M1, so the old period-specific scaling is
# numerically the same reference period as `m1_reference`.
.bp_args <- names(formals(build_mundlak_panel))
if ("scaling_mode" %in% .bp_args) {
  prep <- build_mundlak_panel(dat, period_id = KGRID_PERIOD, scaling_mode = SCALING_MODE)
} else {
  if (!identical(KGRID_PERIOD, "m1_2010_2019") && identical(SCALING_MODE, "m1_reference")) {
    stop("Lokalny 00_helpers.R jest starszy i nie obsługuje scaling_mode='m1_reference' poza M1. ",
         "Dla gridu użyj M1 albo podmień 00_helpers.R na v40.5.")
  }
  message("[compat] build_mundlak_panel() bez argumentu scaling_mode — używam interfejsu v40.2/v40.3.")
  prep <- build_mundlak_panel(dat, period_id = KGRID_PERIOD)
}
dat_model <- prep$data %>%
  model_complete_cases(OUTCOME, prep$all_terms, CONTROL_VARS, require_xy = TRUE) %>%
  arrange(code, year)

between_terms <- prep$between_terms
within_terms  <- prep$within_terms
all_terms     <- prep$all_terms
ref_counties  <- make_reference_counties(dat_model, all_terms, CONTROL_VARS)

cat("Model sample: n=", nrow(dat_model),
    " | counties=", n_distinct(dat_model$code),
    " | years=", min(dat_model$year), "-", max(dat_model$year), "\n\n", sep = "")

# -------------------------
# DWA PORÓWNYWALNE FIT'Y: różnią się wyłącznie K_BASELINE
# -------------------------
grid <- tidyr::crossing(
  k_global = sort(unique(KGRID_GLOBAL_VALUES)),
  k_svc = sort(unique(KGRID_SVC_VALUES)),
  k_baseline = sort(unique(KGRID_BASELINE_VALUES)),
  k_year = sort(unique(KGRID_YEAR_VALUES))
) %>%
  mutate(
    candidate_id = sprintf("glob%02d_svc%02d_base%03d_year%02d",
                           k_global, k_svc, k_baseline, k_year),
    # tylko ranking kosztu/złożoności, nie kryterium statystyczne
    complexity = k_global * 2 + k_svc * 17 + k_baseline + k_year
  ) %>%
  arrange(k_global, k_year, k_baseline, k_svc)

readr::write_csv(grid, file.path(OUT_DIR, "grid_definition.csv"))

signature <- paste(
  "script", "kbaseline_boundary_175_200_v1",
  "period", KGRID_PERIOD,
  "arch", KGRID_ARCHITECTURE,
  "kg", paste(KGRID_GLOBAL_VALUES, collapse = ","),
  "ksvc", paste(KGRID_SVC_VALUES, collapse = ","),
  "kbase", paste(KGRID_BASELINE_VALUES, collapse = ","),
  "kyear", paste(KGRID_YEAR_VALUES, collapse = ","),
  "scaling", SCALING_MODE,
  "gamma", GAMMA,
  "kcheck_nrep", KGRID_KCHECK_NREP,
  "effect_grid_n", KGRID_EFFECT_GRID_N,
  "slope_eps", KGRID_SLOPE_EPS,
  "pred", paste(c(MODEL_VARS, ONLY_BETWEEN_VARS), collapse = "|"),
  sep = "::"
)

make_spec <- function(k_global, k_svc, k_baseline, k_year) {
  list(
    architecture = KGRID_ARCHITECTURE,
    k_global = as.integer(k_global),
    k_svc = as.integer(k_svc),
    k_baseline = as.integer(k_baseline),
    k_year = as.integer(k_year),
    bs_global = "ts",
    bs_svc = "ts",
    bs_baseline = "tp",
    bs_year = "tp"
  )
}

# Local formula builder. Do not use the helper's build_gam_formula here because
# v40.2/v40.3 hard-coded K_* and m=1 and has no `spec=` argument.
# This guarantees that every candidate really uses the requested grid values.
build_kgrid_formula <- function(between_terms, within_terms, controls, spec) {
  arch <- spec$architecture
  allowed <- c("A_global_only", "B_between_svc", "C_all_svc")
  if (!arch %in% allowed) stop("Nieznana architektura GAM: ", arch)

  year_term <- sprintf("s(year_c,bs='%s',k=%d)",
                       spec$bs_year, as.integer(spec$k_year))
  baseline_term <- sprintf("s(x,y,bs='%s',k=%d)",
                           spec$bs_baseline, as.integer(spec$k_baseline))

  global_terms <- c(
    vapply(between_terms, function(v)
      sprintf("s(%s,bs='%s',k=%d)", v, spec$bs_global, as.integer(spec$k_global)),
      character(1)),
    vapply(within_terms, function(v)
      sprintf("s(%s,bs='%s',k=%d)", v, spec$bs_global, as.integer(spec$k_global)),
      character(1))
  )

  svc_between <- if (arch %in% c("B_between_svc", "C_all_svc")) {
    vapply(between_terms, function(v)
      sprintf("s(x,y,by=%s,bs='%s',k=%d)", v, spec$bs_svc, as.integer(spec$k_svc)),
      character(1))
  } else character(0)

  svc_within <- if (identical(arch, "C_all_svc")) {
    vapply(within_terms, function(v)
      sprintf("s(x,y,by=%s,bs='%s',k=%d)", v, spec$bs_svc, as.integer(spec$k_svc)),
      character(1))
  } else character(0)

  rhs <- c(year_term, controls, baseline_term, global_terms, svc_between, svc_within)
  stats::as.formula(paste(OUTCOME, "~", paste(rhs, collapse = " + ")))
}

# -------------------------
# HELPERY DIAGNOSTYCZNE
# -------------------------
classify_smooth <- function(term) {
  if (grepl("year_c", term, fixed = TRUE)) return("year")
  if (identical(gsub(" ", "", term), "s(x,y)")) return("baseline")
  if (grepl("s(x,y):", gsub(" ", "", term), fixed = TRUE)) return("svc")
  "global"
}

safe_skew <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^3)
}

safe_kurtosis_excess <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^4) - 3
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, names = FALSE, type = 7))
}

safe_cor <- function(x, y, min_n = 3L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_n) return(NA_real_)
  sx <- stats::sd(x[ok]); sy <- stats::sd(y[ok])
  if (!is.finite(sx) || !is.finite(sy) || sx < 1e-10 || sy < 1e-10) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok]))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

median_resid_acf1 <- function(d, resid_vec) {
  z <- d %>% mutate(.resid = resid_vec) %>% arrange(code, year) %>%
    group_by(code) %>%
    summarise(
      acf1 = if (n() >= 3L) suppressWarnings(cor(.resid[-n()], .resid[-1], use = "complete.obs")) else NA_real_,
      .groups = "drop"
    )
  stats::median(z$acf1, na.rm = TRUE)
}

# Moran I na średniej reszcie powiatowej. kNN=8 jest tani i niezależny od shapefile.
resid_moran_knn <- function(d, resid_vec, k = 8L) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    return(tibble(moran_i = NA_real_, moran_p = NA_real_, moran_method = "spdep_not_installed"))
  }
  dd <- d %>% mutate(.resid = resid_vec) %>%
    group_by(code, code_pow) %>%
    summarise(x = first(x), y = first(y), resid_mean = mean(.resid, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(x), is.finite(y), is.finite(resid_mean))
  if (nrow(dd) < (k + 2L)) return(tibble(moran_i = NA_real_, moran_p = NA_real_, moran_method = "too_few_counties"))
  coords <- as.matrix(dd[, c("x", "y")])
  knn <- tryCatch(spdep::knearneigh(coords, k = min(k, nrow(dd) - 1L)), error = function(e) NULL)
  if (is.null(knn)) return(tibble(moran_i = NA_real_, moran_p = NA_real_, moran_method = "knn_error"))
  nb <- spdep::knn2nb(knn)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  mt <- tryCatch(spdep::moran.test(dd$resid_mean, lw, zero.policy = TRUE), error = function(e) NULL)
  if (is.null(mt)) return(tibble(moran_i = NA_real_, moran_p = NA_real_, moran_method = "moran_error"))
  tibble(
    moran_i = unname(mt$estimate[["Moran I statistic"]]),
    moran_p = mt$p.value,
    moran_method = "mean county residual; 8-NN"
  )
}


# Moran's I rok po roku z użyciem geometrii powiatów (queen contiguity),
# zgodnie z głównym skryptem walidacyjnym. Jeśli geometria/spdep są
# niedostępne, zwracamy pustą tabelę i grid nadal działa.
resid_moran_by_year <- function(d, resid_vec, powiaty_spdf) {
  if (is.null(powiaty_spdf) || !requireNamespace("spdep", quietly = TRUE)) {
    return(tibble())
  }

  pow <- tryCatch(sf::st_as_sf(powiaty_spdf), error = function(e) NULL)
  if (is.null(pow)) return(tibble())

  code_candidates <- unique(c(
    if (exists("POW_CODE_COL", inherits = TRUE)) get("POW_CODE_COL", inherits = TRUE) else NULL,
    "JPT_KOD_JE", "code_pow", "code"
  ))
  code_col <- intersect(code_candidates, names(pow))[1]
  if (is.na(code_col) || length(code_col) == 0L) return(tibble())

  pow[[code_col]] <- as.character(pow[[code_col]])
  pow[[code_col]] <- ifelse(nchar(pow[[code_col]]) == 3L,
                            paste0("0", pow[[code_col]]),
                            pow[[code_col]])
  pow <- pow %>% arrange(.data[[code_col]])

  nb <- tryCatch(spdep::poly2nb(pow, queen = TRUE), error = function(e) NULL)
  if (is.null(nb)) return(tibble())
  lw <- tryCatch(spdep::nb2listw(nb, style = "W", zero.policy = TRUE),
                 error = function(e) NULL)
  if (is.null(lw)) return(tibble())

  dd_all <- d %>% mutate(.resid = resid_vec)

  purrr::map_dfr(sort(unique(dd_all$year)), function(yr) {
    dd <- dd_all %>% filter(year == yr)
    idx <- match(pow[[code_col]], dd$code_pow)
    rr <- dd$.resid[idx]
    ok <- is.finite(rr)

    if (sum(ok) < 10L) {
      return(tibble(
        year = yr, moran_i = NA_real_,
        p_greater = NA_real_, p_two_sided = NA_real_
      ))
    }

    mt_g <- tryCatch(
      spdep::moran.test(rr, lw, zero.policy = TRUE,
                        na.action = na.exclude, alternative = "greater"),
      error = function(e) NULL
    )
    mt_2 <- tryCatch(
      spdep::moran.test(rr, lw, zero.policy = TRUE,
                        na.action = na.exclude, alternative = "two.sided"),
      error = function(e) NULL
    )

    tibble(
      year = yr,
      moran_i = if (is.null(mt_g)) NA_real_
                else unname(mt_g$estimate[["Moran I statistic"]]),
      p_greater = if (is.null(mt_g)) NA_real_ else mt_g$p.value,
      p_two_sided = if (is.null(mt_2)) NA_real_ else mt_2$p.value
    )
  })
}

extract_k_diag <- function(mod, candidate_id, spec) {
  kd <- tryCatch(
    as.data.frame(mgcv::k.check(mod, n.rep = KGRID_KCHECK_NREP)),
    error = function(e) NULL
  )
  if (is.null(kd) || nrow(kd) == 0L || ncol(kd) < 4L) return(tibble())
  kd <- tibble::rownames_to_column(kd, "term")
  names(kd)[2:5] <- c("k_prime", "edf", "k_index", "p_value")
  kd %>%
    mutate(
      term_group = vapply(term, classify_smooth, character(1)),
      edf_ratio = ifelse(is.finite(k_prime) & k_prime > 0, edf / k_prime, NA_real_),
      # soft = sygnał do obejrzenia; hard = mocniejszy sygnał niedowymiarowania.
      soft_k_flag = is.finite(k_index) & is.finite(p_value) & k_index < 0.85 & p_value < 0.05,
      hard_k_flag = is.finite(k_index) & is.finite(p_value) & is.finite(edf_ratio) &
        k_index < 0.80 & p_value < 0.01 & edf_ratio > 0.70,
      saturation_flag = is.finite(edf_ratio) & edf_ratio >= 0.90,
      candidate_id = candidate_id,
      k_global = spec$k_global,
      k_svc = spec$k_svc,
      k_baseline = spec$k_baseline,
      k_year = spec$k_year
    )
}

extract_edf_table <- function(mod, candidate_id, spec) {
  sm <- summary(mod)
  st <- as.data.frame(sm$s.table)
  if (!nrow(st)) return(tibble())
  tibble::rownames_to_column(st, "term") %>%
    as_tibble() %>%
    mutate(
      term_group = vapply(term, classify_smooth, character(1)),
      candidate_id = candidate_id,
      k_global = spec$k_global,
      k_svc = spec$k_svc,
      k_baseline = spec$k_baseline,
      k_year = spec$k_year
    )
}

# Stałe gridy dla wszystkich kandydatów, żeby porównania były punkt-w-punkt.
term_grids <- setNames(lapply(all_terms, function(v) {
  q <- quantile(dat_model[[v]], c(.05, .95), na.rm = TRUE)
  if (!all(is.finite(q)) || q[1] == q[2]) return(numeric(0))
  seq(q[1], q[2], length.out = KGRID_EFFECT_GRID_N)
}), all_terms)

# Mini-APS: średnia po 7 punktach gridu z lokalnych pochodnych, osobno dla powiatu.
# To jest screening stability, nie zamiennik finalnego APS z 02_spatial_mundlak_gam.R.
effect_signature <- function(mod, candidate_id) {
  county_all <- list()
  curve_all <- list()
  cc <- 0L
  for (v in all_terms) {
    vg <- term_grids[[v]]
    if (!length(vg)) next
    slope_mat <- matrix(NA_real_, nrow = nrow(ref_counties), ncol = length(vg))
    curve_mat <- matrix(NA_real_, nrow = nrow(ref_counties), ncol = length(vg))

    nd_ref <- ref_counties
    nd_ref[[v]] <- 0
    pref <- tryCatch(as.numeric(predict(mod, newdata = nd_ref, type = "response")),
                     error = function(e) rep(NA_real_, nrow(nd_ref)))

    for (j in seq_along(vg)) {
      vv <- vg[j]
      nd0 <- ref_counties; nd1 <- ref_counties; ndc <- ref_counties
      nd0[[v]] <- vv
      nd1[[v]] <- vv + KGRID_SLOPE_EPS
      ndc[[v]] <- vv
      p0 <- tryCatch(as.numeric(predict(mod, newdata = nd0, type = "response")),
                     error = function(e) rep(NA_real_, nrow(nd0)))
      p1 <- tryCatch(as.numeric(predict(mod, newdata = nd1, type = "response")),
                     error = function(e) rep(NA_real_, nrow(nd1)))
      pc <- tryCatch(as.numeric(predict(mod, newdata = ndc, type = "response")),
                     error = function(e) rep(NA_real_, nrow(ndc)))
      slope_mat[, j] <- (p1 - p0) / KGRID_SLOPE_EPS
      curve_mat[, j] <- pc - pref
    }

    cc <- cc + 1L
    county_all[[cc]] <- tibble(
      candidate_id = candidate_id,
      var = v,
      variable = clean_var(v),
      var_label = label_var(clean_var(v)),
      component = comp_from_term(v),
      code = ref_counties$code,
      code_pow = ref_counties$code_pow,
      mini_aps_county = rowMeans(slope_mat, na.rm = TRUE)
    )
    curve_all[[cc]] <- tibble(
      candidate_id = candidate_id,
      var = v,
      variable = clean_var(v),
      var_label = label_var(clean_var(v)),
      component = comp_from_term(v),
      var_value = vg,
      effect_median = apply(curve_mat, 2, median, na.rm = TRUE)
    )
  }
  list(county = bind_rows(county_all), curve = bind_rows(curve_all))
}

fit_one_candidate <- function(row) {
  id <- row$candidate_id[[1]]
  spec <- make_spec(row$k_global[[1]], row$k_svc[[1]], row$k_baseline[[1]], row$k_year[[1]])
  chk_file <- file.path(CHK_DIR, paste0(id, ".rds"))

  if (KGRID_RESUME && file.exists(chk_file)) {
    old <- tryCatch(readRDS(chk_file), error = function(e) NULL)
    if (!is.null(old) && identical(old$signature, signature)) {
      cat("[resume]", id, "\n")
      return(old)
    }
  }

  cat("\n------------------------------------------------------------\n")
  cat("FIT", id, "| K_GLOBAL=", spec$k_global,
      " K_SVC=", spec$k_svc,
      " K_BASELINE=", spec$k_baseline,
      " K_YEAR=", spec$k_year, "\n", sep = "")

  fml <- build_kgrid_formula(between_terms, within_terms, CONTROL_VARS, spec = spec)
  t0 <- proc.time()[["elapsed"]]
  err <- NULL
  mod <- tryCatch(
    mgcv::bam(
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
    ),
    error = function(e) { err <<- conditionMessage(e); NULL }
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  if (is.null(mod)) {
    obj <- list(
      signature = signature,
      summary = tibble(
        candidate_id = id, k_global = spec$k_global, k_svc = spec$k_svc,
        k_baseline = spec$k_baseline, k_year = spec$k_year,
        fit_ok = FALSE, error = err, elapsed_sec = elapsed
      ),
      kdiag = tibble(), edf = tibble(), effect_county = tibble(), effect_curve = tibble(),
      moran_by_year = tibble()
    )
    saveRDS(obj, chk_file)
    return(obj)
  }

  fitted0 <- as.numeric(predict(mod, newdata = dat_model, type = "response"))
  resid0 <- as.numeric(dat_model[[OUTCOME]]) - fitted0
  sm <- summary(mod)
  kd <- extract_k_diag(mod, id, spec)
  edf <- extract_edf_table(mod, id, spec)
  mor <- resid_moran_knn(dat_model, resid0)
  mor_year <- resid_moran_by_year(dat_model, resid0, powiaty_spdf) %>%
    mutate(candidate_id = id,
           k_global = spec$k_global,
           k_svc = spec$k_svc,
           k_baseline = spec$k_baseline,
           k_year = spec$k_year)
  eff <- effect_signature(mod, id)

  kd_sum <- if (nrow(kd)) {
    kd %>% summarise(
      min_k_index = suppressWarnings(min(k_index, na.rm = TRUE)),
      max_edf_ratio = suppressWarnings(max(edf_ratio, na.rm = TRUE)),
      n_soft_k_flags = sum(soft_k_flag, na.rm = TRUE),
      n_hard_k_flags = sum(hard_k_flag, na.rm = TRUE),
      n_saturation_flags = sum(saturation_flag, na.rm = TRUE),
      year_edf_ratio = suppressWarnings(max(edf_ratio[term_group == "year"], na.rm = TRUE)),
      baseline_edf_ratio = suppressWarnings(max(edf_ratio[term_group == "baseline"], na.rm = TRUE)),
      svc_max_edf_ratio = suppressWarnings(max(edf_ratio[term_group == "svc"], na.rm = TRUE)),
      svc_min_k_index = suppressWarnings(min(k_index[term_group == "svc"], na.rm = TRUE))
    )
  } else tibble(
    min_k_index = NA_real_, max_edf_ratio = NA_real_, n_soft_k_flags = NA_integer_,
    n_hard_k_flags = NA_integer_, n_saturation_flags = NA_integer_,
    year_edf_ratio = NA_real_, baseline_edf_ratio = NA_real_,
    svc_max_edf_ratio = NA_real_, svc_min_k_index = NA_real_
  )

  # zamień Inf/-Inf po pustych grupach na NA
  kd_sum <- kd_sum %>% mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))

  summary_row <- tibble(
    candidate_id = id,
    k_global = spec$k_global,
    k_svc = spec$k_svc,
    k_baseline = spec$k_baseline,
    k_year = spec$k_year,
    fit_ok = TRUE,
    error = NA_character_,
    n = nrow(dat_model),
    n_coef = length(coef(mod)),
    n_sp = length(mod$sp),
    model_size_mb = as.numeric(object.size(mod)) / 1024^2,
    aic = tryCatch(as.numeric(AIC(mod)), error = function(e) NA_real_),
    bic = tryCatch(as.numeric(BIC(mod)), error = function(e) NA_real_),
    freml_score = tryCatch(as.numeric(mod$gcv.ubre), error = function(e) NA_real_),
    r_sq = sm$r.sq,
    dev_expl = sm$dev.expl,
    total_smooth_edf = sum(sm$s.table[, "edf"], na.rm = TRUE),
    rmse_in = sqrt(mean(resid0^2, na.rm = TRUE)),
    mae_in = mean(abs(resid0), na.rm = TRUE),
    resid_sd = sd(resid0, na.rm = TRUE),
    resid_skew = safe_skew(resid0),
    resid_kurtosis_excess = safe_kurtosis_excess(resid0),
    resid_acf1_median = median_resid_acf1(dat_model, resid0),
    moran_i = mor$moran_i[[1]],
    moran_p = mor$moran_p[[1]],
    moran_year_median = if (nrow(mor_year)) median(mor_year$moran_i, na.rm = TRUE) else NA_real_,
    moran_year_max_abs = if (nrow(mor_year)) max(abs(mor_year$moran_i), na.rm = TRUE) else NA_real_,
    n_moran_year_pos_p05 = if (nrow(mor_year)) sum(mor_year$moran_i > 0 & mor_year$p_greater < .05, na.rm = TRUE) else NA_integer_,
    n_moran_year_twosided_p05 = if (nrow(mor_year)) sum(mor_year$p_two_sided < .05, na.rm = TRUE) else NA_integer_,
    elapsed_sec = elapsed
  ) %>% bind_cols(kd_sum)

  obj <- list(
    signature = signature,
    summary = summary_row,
    kdiag = kd,
    edf = edf,
    effect_county = eff$county,
    effect_curve = eff$curve,
    moran_by_year = mor_year
  )
  saveRDS(obj, chk_file)
  if (KGRID_SAVE_MODELS) saveRDS(mod, file.path(MOD_DIR, paste0(id, ".rds")))

  cat(sprintf("  AIC=%.1f | R2=%.4f | dev=%.4f | max EDF/k'=%.3f | min k-index=%.3f | %.1f min\n",
              summary_row$aic, summary_row$r_sq, summary_row$dev_expl,
              summary_row$max_edf_ratio, summary_row$min_k_index, elapsed / 60))

  rm(mod, fitted0, resid0, eff, mor_year)
  invisible(gc(full = TRUE))
  obj
}

# -------------------------
# STAGE 1: DWA PEŁNE FIT'Y
# -------------------------
results <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  results[[i]] <- fit_one_candidate(grid[i, ])

  # progres po każdym kandydacie
  tmp_sum <- bind_rows(lapply(results[seq_len(i)], `[[`, "summary"))
  readr::write_csv(tmp_sum, file.path(OUT_DIR, "grid_summary_progress.csv"))
}

summary_tbl <- bind_rows(lapply(results, `[[`, "summary")) %>%
  left_join(grid %>% select(candidate_id, complexity), by = "candidate_id") %>%
  mutate(
    delta_aic = ifelse(fit_ok & is.finite(aic), aic - min(aic[fit_ok & is.finite(aic)], na.rm = TRUE), NA_real_),
    delta_bic = ifelse(fit_ok & is.finite(bic), bic - min(bic[fit_ok & is.finite(bic)], na.rm = TRUE), NA_real_),
    # Tylko flaga diagnostyczna do przeglądu; NIE jest automatycznym kryterium odrzucenia.
    diagnostic_attention = fit_ok & (
      (!is.na(max_edf_ratio) & max_edf_ratio >= 0.95) |
      (!is.na(n_hard_k_flags) & n_hard_k_flags > 0)
    )
  ) %>%
  arrange(k_global, k_year, k_baseline, k_svc)

kdiag_tbl <- bind_rows(lapply(results, `[[`, "kdiag"))
edf_tbl <- bind_rows(lapply(results, `[[`, "edf"))
effect_county_tbl <- bind_rows(lapply(results, `[[`, "effect_county"))
effect_curve_tbl <- bind_rows(lapply(results, `[[`, "effect_curve"))
moran_by_year_tbl <- bind_rows(lapply(results, function(x) {
  if (!is.null(x$moran_by_year)) x$moran_by_year else tibble()
}))

readr::write_csv(summary_tbl, file.path(OUT_DIR, "grid_summary.csv"))
readr::write_csv(kdiag_tbl, file.path(OUT_DIR, "kcheck_by_term.csv"))
readr::write_csv(edf_tbl, file.path(OUT_DIR, "edf_by_term.csv"))
readr::write_csv(moran_by_year_tbl, file.path(OUT_DIR, "moran_by_year.csv"))

# Term-by-term / group-level view: pokazuje CO generuje max EDF/k' i min k-index.
kcheck_group_summary <- kdiag_tbl %>%
  group_by(candidate_id, k_global, k_svc, k_baseline, k_year, term_group) %>%
  summarise(
    n_terms = n(),
    max_edf_ratio = suppressWarnings(max(edf_ratio, na.rm = TRUE)),
    min_k_index = suppressWarnings(min(k_index, na.rm = TRUE)),
    n_edf_ratio_ge_090 = sum(edf_ratio >= 0.90, na.rm = TRUE),
    n_edf_ratio_ge_095 = sum(edf_ratio >= 0.95, na.rm = TRUE),
    n_kindex_lt_085 = sum(k_index < 0.85, na.rm = TRUE),
    n_kindex_lt_080 = sum(k_index < 0.80, na.rm = TRUE),
    n_kcheck_p_lt_005 = sum(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    max_edf_ratio = ifelse(is.infinite(max_edf_ratio), NA_real_, max_edf_ratio),
    min_k_index = ifelse(is.infinite(min_k_index), NA_real_, min_k_index)
  )

critical_terms <- kdiag_tbl %>%
  group_by(candidate_id) %>%
  summarise(
    worst_edf_term = if (any(is.finite(edf_ratio))) term[which.max(replace(edf_ratio, !is.finite(edf_ratio), -Inf))] else NA_character_,
    worst_edf_group = if (any(is.finite(edf_ratio))) term_group[which.max(replace(edf_ratio, !is.finite(edf_ratio), -Inf))] else NA_character_,
    worst_edf_ratio = if (any(is.finite(edf_ratio))) max(edf_ratio, na.rm = TRUE) else NA_real_,
    worst_k_term = if (any(is.finite(k_index))) term[which.min(replace(k_index, !is.finite(k_index), Inf))] else NA_character_,
    worst_k_group = if (any(is.finite(k_index))) term_group[which.min(replace(k_index, !is.finite(k_index), Inf))] else NA_character_,
    worst_k_index = if (any(is.finite(k_index))) min(k_index, na.rm = TRUE) else NA_real_,
    worst_k_p = if (any(is.finite(k_index))) p_value[which.min(replace(k_index, !is.finite(k_index), Inf))] else NA_real_,
    .groups = "drop"
  )

readr::write_csv(kcheck_group_summary, file.path(OUT_DIR, "kcheck_group_summary.csv"))
readr::write_csv(critical_terms, file.path(OUT_DIR, "critical_terms.csv"))
readr::write_csv(effect_county_tbl, file.path(OUT_DIR, "mini_aps_county_signatures.csv"))
readr::write_csv(effect_curve_tbl, file.path(OUT_DIR, "curve_signatures.csv"))

# -------------------------
# STABILNOŚĆ WZGLĘDEM ZWIĘKSZANIA k
# porównujemy tylko sąsiednie wartości jednej osi przy pozostałych stałych.
# -------------------------
make_adjacent_pairs <- function(grid, dim_name, values) {
  values <- sort(unique(values))
  pair_list <- list()
  z <- 0L
  for (j in seq_len(length(values) - 1L)) {
    lo <- values[j]; hi <- values[j + 1L]
    other <- setdiff(c("k_global", "k_svc", "k_baseline", "k_year"), dim_name)
    low_tbl <- grid %>% filter(.data[[dim_name]] == lo)
    high_tbl <- grid %>% filter(.data[[dim_name]] == hi)
    pp <- inner_join(
      low_tbl %>% select(all_of(other), candidate_low = candidate_id),
      high_tbl %>% select(all_of(other), candidate_high = candidate_id),
      by = other
    ) %>% mutate(dimension = dim_name, k_low = lo, k_high = hi)
    z <- z + 1L; pair_list[[z]] <- pp
  }
  bind_rows(pair_list)
}

pairs <- make_adjacent_pairs(grid, "k_baseline", KGRID_BASELINE_VALUES)

compare_pair_term <- function(low_id, high_id, dim, klo, khi) {
  surf <- inner_join(
    effect_county_tbl %>% filter(candidate_id == low_id) %>%
      select(var, variable, var_label, component, code, low = mini_aps_county),
    effect_county_tbl %>% filter(candidate_id == high_id) %>%
      select(var, code, high = mini_aps_county),
    by = c("var", "code")
  )
  curve <- inner_join(
    effect_curve_tbl %>% filter(candidate_id == low_id) %>%
      select(var, var_value, curve_low = effect_median),
    effect_curve_tbl %>% filter(candidate_id == high_id) %>%
      select(var, var_value, curve_high = effect_median),
    by = c("var", "var_value")
  )

  surf_sum <- surf %>% group_by(var, variable, var_label, component) %>%
    summarise(
      aps_low = median(low, na.rm = TRUE),
      aps_high = median(high, na.rm = TRUE),
      abs_aps_diff = abs(aps_high - aps_low),
      rel_aps_diff = abs_aps_diff / max(1, abs(aps_low), abs(aps_high), na.rm = TRUE),
      aps_meaningful = max(abs(aps_low), abs(aps_high), na.rm = TRUE) >= KGRID_MEANINGFUL_APS,
      aps_sign_same = ifelse(aps_meaningful, sign(aps_low) == sign(aps_high), NA),
      surface_sd_low = sd(low, na.rm = TRUE),
      surface_sd_high = sd(high, na.rm = TRUE),
      surface_active = max(surface_sd_low, surface_sd_high, na.rm = TRUE) >= KGRID_SURFACE_SD_MIN,
      surface_cor = safe_cor(low, high, min_n = 10L),
      county_sign_agreement = mean(sign(low) == sign(high), na.rm = TRUE),
      surface_rmse = sqrt(mean((high - low)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  curve_sum <- curve %>% group_by(var) %>%
    summarise(
      curve_range_low = diff(range(curve_low, na.rm = TRUE)),
      curve_range_high = diff(range(curve_high, na.rm = TRUE)),
      curve_active = max(curve_range_low, curve_range_high, na.rm = TRUE) >= KGRID_CURVE_RANGE_MIN,
      curve_cor = safe_cor(curve_low, curve_high, min_n = 3L),
      curve_rmse = sqrt(mean((curve_high - curve_low)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  surf_sum %>% left_join(curve_sum, by = "var") %>%
    mutate(
      dimension = dim, k_low = klo, k_high = khi,
      candidate_low = low_id, candidate_high = high_id
    )
}

adjacent_term <- pmap_dfr(
  pairs %>% select(candidate_low, candidate_high, dimension, k_low, k_high),
  ~ compare_pair_term(..1, ..2, ..3, ..4, ..5)
)

adjacent_candidate <- adjacent_term %>%
  group_by(dimension, k_low, k_high, candidate_low, candidate_high) %>%
  summarise(
    n_terms = n(),
    aps_vector_cor = safe_cor(aps_low, aps_high, min_n = 5L),
    n_meaningful_aps = sum(aps_meaningful, na.rm = TRUE),
    aps_sign_agreement = ifelse(n_meaningful_aps > 0,
                                mean(aps_sign_same[aps_meaningful], na.rm = TRUE), NA_real_),
    median_rel_aps_diff = median(rel_aps_diff, na.rm = TRUE),
    max_rel_aps_diff = max(rel_aps_diff, na.rm = TRUE),
    n_active_surfaces = sum(surface_active, na.rm = TRUE),
    median_surface_cor = median(surface_cor[surface_active], na.rm = TRUE),
    p10_surface_cor = safe_quantile(surface_cor[surface_active], .10),
    n_active_curves = sum(curve_active, na.rm = TRUE),
    median_curve_cor = median(curve_cor[curve_active], na.rm = TRUE),
    p10_curve_cor = safe_quantile(curve_cor[curve_active], .10),
    mean_county_sign_agreement = mean(county_sign_agreement, na.rm = TRUE),
    median_abs_aps_diff = median(abs_aps_diff, na.rm = TRUE),
    max_abs_aps_diff = max(abs_aps_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Brak aktywnych powierzchni/krzywych nie jest porażką: shrinkage może je wyzerować.
    surface_ok = n_active_surfaces == 0L | is.na(p10_surface_cor) |
      p10_surface_cor >= KGRID_STAB_SURFACE_COR_MIN,
    curve_ok = n_active_curves == 0L | is.na(p10_curve_cor) |
      p10_curve_cor >= KGRID_STAB_CURVE_COR_MIN,
    aps_ok = is.na(aps_vector_cor) | aps_vector_cor >= KGRID_STAB_APS_COR_MIN,
    aps_sign_ok = n_meaningful_aps == 0L | is.na(aps_sign_agreement) |
      aps_sign_agreement >= KGRID_STAB_APS_SIGN_MIN,
    stability_pass = aps_ok & aps_sign_ok & curve_ok & surface_ok
  )

# Interesuje nas szczególnie, czy mniejsza baza jest już stabilna
# względem następnego większego k przy pozostałych wymiarach stałych.
lower_stability <- adjacent_candidate %>%
  group_by(candidate_id = candidate_low) %>%
  summarise(
    n_upward_pairs = n(),
    n_stable_upward = sum(stability_pass, na.rm = TRUE),
    n_unstable_upward = sum(!stability_pass, na.rm = TRUE),
    min_aps_vector_cor = suppressWarnings(min(aps_vector_cor, na.rm = TRUE)),
    min_aps_sign_agreement = suppressWarnings(min(aps_sign_agreement, na.rm = TRUE)),
    min_p10_surface_cor = suppressWarnings(min(p10_surface_cor, na.rm = TRUE)),
    min_p10_curve_cor = suppressWarnings(min(p10_curve_cor, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))

readr::write_csv(adjacent_term, file.path(OUT_DIR, "adjacent_k_stability_by_term.csv"))
readr::write_csv(adjacent_candidate, file.path(OUT_DIR, "adjacent_k_stability_summary.csv"))
readr::write_csv(lower_stability, file.path(OUT_DIR, "candidate_upward_stability.csv"))

# -------------------------
# PROSTE PLOTY BOUNDARY CHECK
# -------------------------
plot_metrics <- summary_tbl %>%
  filter(fit_ok) %>%
  select(k_baseline, aic, bic, rmse_in, r_sq, dev_expl,
         baseline_edf_ratio, min_k_index,
         moran_year_median, moran_year_max_abs) %>%
  pivot_longer(-k_baseline, names_to = "metric", values_to = "value")

p_metrics <- ggplot(plot_metrics,
                    aes(x = factor(k_baseline), y = value, group = 1)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  labs(
    title = "K_BASELINE boundary check: 175 vs 200",
    x = "K_BASELINE", y = NULL
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(OUT_DIR, "baseline_175_200_metrics.png"),
       p_metrics, width = 10, height = 8, dpi = 180)

# -------------------------
# RAPORT KOŃCOWY / META — BEZ PERTURBACJI
# -------------------------
# Zestawienie zmian fitu przy zwiększaniu pojedynczego k.
fit_delta_pairs <- pairs %>%
  left_join(summary_tbl %>% select(candidate_low = candidate_id,
                                   aic_low = aic, bic_low = bic, dev_low = dev_expl,
                                   r2_low = r_sq, rmse_low = rmse_in,
                                   moran_low = moran_i, acf_low = resid_acf1_median,
                                   max_edf_low = max_edf_ratio, min_k_low = min_k_index),
            by = "candidate_low") %>%
  left_join(summary_tbl %>% select(candidate_high = candidate_id,
                                   aic_high = aic, bic_high = bic, dev_high = dev_expl,
                                   r2_high = r_sq, rmse_high = rmse_in,
                                   moran_high = moran_i, acf_high = resid_acf1_median,
                                   max_edf_high = max_edf_ratio, min_k_high = min_k_index),
            by = "candidate_high") %>%
  mutate(
    delta_aic_high_minus_low = aic_high - aic_low,
    delta_bic_high_minus_low = bic_high - bic_low,
    delta_dev_high_minus_low = dev_high - dev_low,
    delta_r2_high_minus_low = r2_high - r2_low,
    delta_rmse_high_minus_low = rmse_high - rmse_low,
    delta_moran_high_minus_low = moran_high - moran_low,
    delta_acf_high_minus_low = acf_high - acf_low
  ) %>%
  left_join(adjacent_candidate,
            by = c("dimension", "k_low", "k_high", "candidate_low", "candidate_high"))

readr::write_csv(fit_delta_pairs, file.path(OUT_DIR, "adjacent_k_fit_and_effect_changes.csv"))

meta <- c(
  paste0("created: ", Sys.time()),
  paste0("period: ", KGRID_PERIOD),
  paste0("architecture: ", KGRID_ARCHITECTURE),
  paste0("K_GLOBAL: ", paste(KGRID_GLOBAL_VALUES, collapse = ",")),
  paste0("K_SVC: ", paste(KGRID_SVC_VALUES, collapse = ",")),
  paste0("K_BASELINE: ", paste(KGRID_BASELINE_VALUES, collapse = ",")),
  paste0("K_YEAR: ", paste(KGRID_YEAR_VALUES, collapse = ",")),
  paste0("family: gaussian; method: ", BAM_METHOD, "; discrete: ", BAM_DISCRETE, "; gamma: ", GAMMA),
  paste0("scaling: ", SCALING_MODE),
  paste0("k.check n.rep: ", KGRID_KCHECK_NREP),
  paste0("effect signature grid points: ", KGRID_EFFECT_GRID_N),
  "PERTURBATION: OFF by design. Run only after basis dimensions stabilize.",
  "Boundary-check rationale: K_GLOBAL=5, K_SVC=25 and K_YEAR=9 are held fixed; only K_BASELINE is increased from 175 to 200.",
  "Concurvity is intentionally not used as a selection/fail criterion.",
  "Interpret k.check term-by-term. max EDF/k' and min k-index are screening summaries only.",
  "Selection principle: retain K_BASELINE=175 only if 175->200 is substantively stable; otherwise perform one further targeted check at K_BASELINE=225.",
  "Mini-APS/curve signatures are screening diagnostics; final APS must be recomputed by the main GAM script."
)
writeLines(meta, file.path(OUT_DIR, "RUN_INFO.txt"))


# -------------------------
# AUTOMATYCZNE PODSUMOWANIE BOUNDARY CHECK
# -------------------------
boundary_pair <- fit_delta_pairs %>%
  filter(dimension == "k_baseline", k_low == 175, k_high == 200)

boundary_terms <- adjacent_term %>%
  filter(dimension == "k_baseline", k_low == 175, k_high == 200) %>%
  arrange(desc(aps_meaningful), desc(rel_aps_diff))

readr::write_csv(boundary_pair,
                 file.path(OUT_DIR, "baseline_175_200_decision_summary.csv"))
readr::write_csv(boundary_terms,
                 file.path(OUT_DIR, "baseline_175_200_by_term.csv"))

if (nrow(boundary_pair) == 1L) {
  bp <- boundary_pair[1, ]

  recommendation <- if (isTRUE(bp$stability_pass)) {
    "RETAIN_K_BASELINE_175"
  } else {
    "TEST_K_BASELINE_225"
  }

  decision_text <- c(
    "K_BASELINE boundary check: 175 -> 200",
    paste0("screening_recommendation: ", recommendation),
    paste0("stability_pass: ", bp$stability_pass),
    paste0("delta_AIC_200_minus_175: ", round(bp$delta_aic_high_minus_low, 3)),
    paste0("delta_BIC_200_minus_175: ", round(bp$delta_bic_high_minus_low, 3)),
    paste0("delta_RMSE_200_minus_175: ", round(bp$delta_rmse_high_minus_low, 4)),
    paste0("APS_vector_correlation: ", round(bp$aps_vector_cor, 5)),
    paste0("APS_sign_agreement: ", round(bp$aps_sign_agreement, 5)),
    paste0("median_relative_APS_difference: ", round(bp$median_rel_aps_diff, 5)),
    paste0("maximum_relative_APS_difference: ", round(bp$max_rel_aps_diff, 5)),
    paste0("p10_SVC_surface_correlation: ", round(bp$p10_surface_cor, 5)),
    paste0("p10_curve_correlation: ", round(bp$p10_curve_cor, 5)),
    paste0("mean_county_sign_agreement: ", round(bp$mean_county_sign_agreement, 5)),
    "",
    "Interpretation:",
    "- stability_pass is a screening heuristic, not a formal statistical test;",
    "- the final choice should combine effect stability, term-by-term EDF/k', residual spatial structure and parsimony;",
    "- concurvity is intentionally NOT a selection/fail criterion;",
    "- if 175->200 is stable, retain the smaller basis K_BASELINE=175;",
    "- if it is not stable, run a single additional K_BASELINE=225 fit and compare 200->225."
  )
  writeLines(decision_text, file.path(OUT_DIR, "BOUNDARY_DECISION.txt"))
}

cat("\n============================================================\n")
cat("K_BASELINE BOUNDARY CHECK 175 -> 200 DONE\n")
cat("Output:", OUT_DIR, "\n")
cat("Najważniejsze pliki:\n")
cat("  grid_summary.csv\n")
cat("  kcheck_by_term.csv\n")
cat("  critical_terms.csv\n")
cat("  edf_by_term.csv\n")
cat("  moran_by_year.csv\n")
cat("  baseline_175_200_decision_summary.csv\n")
cat("  baseline_175_200_by_term.csv\n")
cat("  adjacent_k_fit_and_effect_changes.csv\n")
cat("  BOUNDARY_DECISION.txt\n")
cat("============================================================\n\n")

print(summary_tbl %>%
        select(candidate_id, k_global, k_svc, k_baseline, k_year,
               aic, delta_aic, bic, delta_bic, r_sq, dev_expl, rmse_in,
               max_edf_ratio, year_edf_ratio, baseline_edf_ratio,
               svc_max_edf_ratio, min_k_index, n_hard_k_flags,
               resid_acf1_median, moran_i, moran_year_median,
               moran_year_max_abs, n_moran_year_pos_p05,
               elapsed_sec), n = 10)

cat("\nPorównanie 175 -> 200:\n")
print(boundary_pair, n = 10)

cat("\nNajwiększe zmiany APS term-by-term:\n")
print(boundary_terms %>%
        select(var_label, component, aps_low, aps_high,
               abs_aps_diff, rel_aps_diff, aps_sign_same,
               surface_cor, curve_cor) %>%
        arrange(desc(rel_aps_diff)), n = 40)

if (file.exists(file.path(OUT_DIR, "BOUNDARY_DECISION.txt"))) {
  cat("\n")
  cat(paste(readLines(file.path(OUT_DIR, "BOUNDARY_DECISION.txt")),
            collapse = "\n"))
  cat("\n")
}

rm(RUNNING_GAM_SELECTION)
