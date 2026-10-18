# ============================================================
# 00_config.R — konfiguracja Spatial Mundlak GAM, 3 okresy
# SSM final pipeline v40.14
# ============================================================

# PERIOD można ustawić przed source("R/00_config.R").
if (!exists("PERIOD")) PERIOD <- "m1_2010_2019"

PERIOD_SPECS <- list(
  m1_2010_2019 = list(
    year_min = 2010L,
    year_max = 2019L,
    year_label = "2010–2019",
    model_label = "Model 1",
    controls = character(0)
  ),
  m2_2010_2022 = list(
    year_min = 2010L,
    year_max = 2022L,
    year_label = "2010–2022",
    model_label = "Model 2",
    controls = c("covid_dummy")
  ),
  m3_2010_2025 = list(
    year_min = 2010L,
    year_max = 2025L,
    year_label = "2010–2025",
    model_label = "Model 3",
    controls = c("covid_dummy", "post_covid_dummy")
  )
)

if (!PERIOD %in% names(PERIOD_SPECS)) {
  stop("Nieznany PERIOD: ", PERIOD,
       ". Dozwolone: ", paste(names(PERIOD_SPECS), collapse = ", "))
}

PERIOD_SPEC  <- PERIOD_SPECS[[PERIOD]]
YEAR_MIN     <- PERIOD_SPEC$year_min
YEAR_MAX     <- PERIOD_SPEC$year_max
YEAR_LABEL   <- PERIOD_SPEC$year_label
MODEL_LABEL  <- PERIOD_SPEC$model_label
CONTROL_VARS <- PERIOD_SPEC$controls

# ------------------------------------------------------------
# Definicja okresów wskaźnikowych.
# Założenie v40:
#   COVID dummy      = 2020–2022
#   post-COVID dummy = 2023–2025
# Dzięki temu modele są zagnieżdżone i nakładające się:
# M1 2010–2019; M2 2010–2022; M3 2010–2025.
# Jeśli chcesz inną definicję, zmień tylko te dwa wektory.
# ------------------------------------------------------------
COVID_YEARS      <- 2020:2022
POST_COVID_YEARS <- 2023:2025

# ---- Kolumny panelu ----
ID_COL       <- "code"
NAME_COL     <- "name"
YEAR_COL     <- "year"
POW_CODE_COL <- "code_pow"

# Outcome. Jeżeli tej nazwy nie ma, RUN_ALL spróbuje kilku typowych aliasów.
OUTCOME <- "asmr_na_100tys"
OUTCOME_CANDIDATES <- c(
  "asmr_na_100tys", "asmr_2010_2019std", "asmr", "ASMR", "asmr_rate", "asmr_per_100k"
)
OUTCOME_LABEL <- "ASMR / 100,000"

# ---- Predyktory ----
# Between + within (lag 1)
MODEL_VARS <- c(
  "wyksztalcenie_wyzsze_pct",
  "wyksztalcenie_gimnazjalne_i_podstawowe_pct",
  "dlugotrwale_bezrobotni_pct_aktywnych",
  "wynagrodzenie_brutto",
  "lekarze_na_10_tys",
  "pm25_mean_year",
  "saldo_migracji_na_1000",
  "density_pop_km2"
)

# Between-only
ONLY_BETWEEN_VARS <- c(
  "distance_to_city_km"
)

ALL_RAW_VARS <- unique(c(MODEL_VARS, ONLY_BETWEEN_VARS))
WITHIN_LAG_N <- 1L

VAR_LABELS <- c(
  wyksztalcenie_wyzsze_pct = "Higher education",
  wyksztalcenie_gimnazjalne_i_podstawowe_pct = "Primary / lower-secondary education",
  dlugotrwale_bezrobotni_pct_aktywnych = "Long-term unemployment",
  wynagrodzenie_brutto = "Gross wages",
  lekarze_na_10_tys = "Physicians",
  pm25_mean_year = "PM2.5",
  saldo_migracji_na_1000 = "Net migration",
  density_pop_km2 = "Population density",
  distance_to_city_km = "Distance to city with county rights"
)

VAR_UNITS <- c(
  wyksztalcenie_wyzsze_pct = "%",
  wyksztalcenie_gimnazjalne_i_podstawowe_pct = "%",
  dlugotrwale_bezrobotni_pct_aktywnych = "% of economically active",
  wynagrodzenie_brutto = "PLN",
  lekarze_na_10_tys = "per 10,000",
  pm25_mean_year = "µg/m³",
  saldo_migracji_na_1000 = "per 1,000",
  density_pop_km2 = "persons/km²",
  distance_to_city_km = "km"
)

# ---- Polityka braków ----
# Braki nie są imputowane medianą. Model używa kompletnych obserwacji
# po zbudowaniu decomposition/lag. Audyt braków jest zapisywany osobno.
MISSING_POLICY <- "complete_case"

# ---- Standaryzacja porównawcza M1/M2/M3 ----
# Finalne modele używają jednej referencyjnej skali 2010–2019.
# County means Mundlaka nadal są liczone osobno dla okresu danego modelu;
# zamrażamy wyłącznie centrum/SD używane do standaryzacji efektów.
SCALING_MODE <- "m1_reference"
SCALING_REFERENCE_PERIOD <- "m1_2010_2019"
# Strict CV musi pozostać leakage-free, dlatego w foldach skala jest fitowana
# wyłącznie na latach treningowych.
CV_SCALING_MODE <- "training_only"

# ---- Finalna specyfikacja GAM ----
# Basis dimensions zostały zamrożone po pełnej sensitivity analysis na M1:
# 48-model grid + targeted K_BASELINE 175 -> 200 boundary check.
# Ta sama specyfikacja jest używana w M1/M2/M3 dla bezpośredniej porównywalności.
FINAL_GAM_SPEC <- list(
  architecture = "C_all_svc",
  k_global = 5L,
  k_svc = 25L,
  k_baseline = 175L,
  k_year = 9L,
  bs_global = "ts",
  bs_svc = "ts",
  bs_baseline = "tp",
  bs_year = "tp"
)

GAM_SPEC <- FINAL_GAM_SPEC
GAM_ARCHITECTURE <- GAM_SPEC$architecture
K_GLOBAL    <- as.integer(GAM_SPEC$k_global)
K_SVC       <- as.integer(GAM_SPEC$k_svc)
K_BASELINE  <- as.integer(GAM_SPEC$k_baseline)
K_YEAR      <- as.integer(GAM_SPEC$k_year)
BS_GLOBAL   <- GAM_SPEC$bs_global
BS_SVC      <- GAM_SPEC$bs_svc
BS_BASELINE <- GAM_SPEC$bs_baseline
BS_YEAR     <- GAM_SPEC$bs_year

FAMILY       <- gaussian()
GAMMA        <- 1.4
BAM_METHOD   <- "fREML"
BAM_DISCRETE <- TRUE
BAM_NTHREADS <- max(1L, parallel::detectCores() - 1L)
BAM_GC_LEVEL <- 0L


SLOPE_EPS <- 0.01
N_GRID    <- 25L


# ---- SVC uncertainty / reporting ----
# Approximate Bayesian covariance from mgcv is used for linear contrasts.
# When available, unconditional=TRUE adds smoothing-parameter uncertainty.
SVC_CI_LEVEL <- 0.95
SVC_USE_UNCONDITIONAL_VCOV <- TRUE

# "Robust sign" is descriptive shorthand only:
# the approximate pointwise CI for the local contrast does not include zero.
# It is NOT a multiplicity-adjusted local significance test.
SVC_FADE_UNCERTAIN_MAPS <- TRUE
SVC_UNCERTAIN_ALPHA <- 0.24

# Cross-period comparison is a nested-period stability diagnostic,
# not an independent replication study.
SVC_STABILITY_MIN_COMMON_COUNTIES <- 30L

# ---- Mundlak ridge ----
CV_FOLDS    <- 10L
LAMBDA_RULE <- "lambda.min"
BOOT_B      <- 1000L
BOOT_RHO    <- 0.5
BOOT_TYPE   <- "dependent_multiplier_fixed_operator"

# ---- Walidacja ----
RUN_STRICT_LOO <- TRUE
MIN_TRAIN_YEARS_STRICT_CV <- 5L

# ---- Perturbacja ----
B_PERTURB      <- 50L
NOISE_SCALE    <- 0.50
SNR_THRESHOLD  <- 2.0
SHAPE_R_THRESHOLD <- 0.80
SURFACE_R_THRESHOLD <- 0.95
COUNTY_SIGN_THRESHOLD <- 0.90
APS_SIGN_THRESHOLD <- 0.90
# Spatial perturbation stability is evaluated at the common reference z = 0.
PERT_SPATIAL_REFERENCE_Z <- 0
SURFACE_SD_MIN <- 0.25
MEANINGFUL_APS_THRESHOLD <- 1.0
MAIN_MAX_BETWEEN <- 5L
MAIN_MAX_WITHIN  <- 5L

# ---- Output ----
OUTPUT_ROOT <- "output_v40_6_periods_2010_2019_2020_2022_2023_2025"
OUT_BASE <- file.path(OUTPUT_ROOT, PERIOD)

# ---- Figury / Elsevier ----
# User-facing switch for figure export in RUN_ALL / RUN_AFTER_GAM / reporting-only runs.
# Allowed values: "png", "pdf", "both". You may also set FIG_EXPORT_FORMATS directly.
if (!exists("FIG_EXPORT_MODE")) FIG_EXPORT_MODE <- "both"
if (!exists("FIG_EXPORT_FORMATS")) {
  mode_norm <- tolower(trimws(as.character(FIG_EXPORT_MODE)[1]))
  FIG_EXPORT_FORMATS <- switch(
    mode_norm,
    "png" = "png",
    "pdf" = "pdf",
    "both" = c("png", "pdf"),
    stop('FIG_EXPORT_MODE must be one of: "png", "pdf", "both".')
  )
}
if (!exists("FIG_DPI")) FIG_DPI <- 500
if (!exists("FIG_FONT_FAMILY")) FIG_FONT_FAMILY <- "sans"

# ---- Publication cartogram policy (SSM) ----
# Effect cartograms use one symmetric colour range for the marginal and
# accumulated-contribution columns. The range is locked across M1/M2/M3
# within each Mundlak layer (between / within), so the same colour always
# represents the same ASMR contrast. Extreme values are squished only for
# display; model estimates and exported numeric tables are unchanged.
if (!exists("FIG_EFFECT_MAP_CLIP_PROB")) FIG_EFFECT_MAP_CLIP_PROB <- 0.95
if (!exists("FIG_EFFECT_MAP_LOCK_ACROSS_PERIODS")) FIG_EFFECT_MAP_LOCK_ACROSS_PERIODS <- TRUE
if (!exists("FIG_EFFECT_MAP_PAIR_SHARED_SCALE")) FIG_EFFECT_MAP_PAIR_SHARED_SCALE <- TRUE

# Larger maps: curve | marginal map | accumulated-contribution map.
if (!exists("FIG_EFFECT_GRID_WIDTHS")) FIG_EFFECT_GRID_WIDTHS <- c(1.12, 1.00, 1.00)
if (!exists("FIG_UNCERTAINTY_GRID_WIDTHS")) FIG_UNCERTAINTY_GRID_WIDTHS <- c(0.64, 1.00, 1.00)

# Cartographic line weights.
if (!exists("FIG_MAP_COUNTY_BORDER_WIDTH")) FIG_MAP_COUNTY_BORDER_WIDTH <- 0.045
if (!exists("FIG_MAP_COUNTRY_BORDER_WIDTH")) FIG_MAP_COUNTRY_BORDER_WIDTH <- 0.52

FIG_W_SINGLE <- 90 / 25.4
FIG_W_MID    <- 140 / 25.4
FIG_W_FULL   <- 190 / 25.4

HAP_COLS <- list(
  blue       = "#0072B2",
  sky        = "#56B4E9",
  orange     = "#E69F00",
  vermillion = "#D55E00",
  green      = "#009E73",
  purple     = "#CC79A7",
  yellow     = "#F0E442",
  black      = "#1A1A1A",
  grey_dark  = "#4D4D4D",
  grey_mid   = "#7F7F7F",
  grey_light = "#E6E6E6",
  grey_panel = "#F7F7F7",
  white      = "#FFFFFF"
)

LANCET <- list(
  red        = HAP_COLS$vermillion,
  navy       = HAP_COLS$blue,
  grey_dark  = HAP_COLS$grey_dark,
  grey_mid   = HAP_COLS$grey_mid,
  grey_light = HAP_COLS$grey_light,
  white      = HAP_COLS$white,
  black      = HAP_COLS$black
)

label_var <- function(x) {
  out <- VAR_LABELS[x]
  out[is.na(out)] <- x[is.na(out)]
  as.character(out)
}

unit_var <- function(x) {
  out <- VAR_UNITS[x]
  out[is.na(out)] <- ""
  as.character(out)
}

clean_var <- function(t) {
  t <- as.character(t)
  t <- stringr::str_remove(t, "_between_z$")
  t <- stringr::str_remove(t, paste0("_within_lag", WITHIN_LAG_N, "_z$"))
  t
}

comp_from_term <- function(t) {
  dplyr::case_when(
    stringr::str_detect(t, "_between_z$") ~ "between",
    stringr::str_detect(t, paste0("_within_lag", WITHIN_LAG_N, "_z$")) ~ "within",
    TRUE ~ "other"
  )
}

make_period_controls <- function(d) {
  d %>%
    dplyr::mutate(
      covid_dummy = as.integer(.data[[YEAR_COL]] %in% COVID_YEARS),
      post_covid_dummy = as.integer(.data[[YEAR_COL]] %in% POST_COVID_YEARS)
    )
}

theme_health_place <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = FIG_FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = HAP_COLS$black, family = FIG_FONT_FAMILY),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0,
                                         margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = HAP_COLS$grey_mid,
                                            hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.caption = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.line = ggplot2::element_line(colour = HAP_COLS$grey_dark, linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = HAP_COLS$grey_dark, linewidth = 0.30),
      axis.ticks.length = ggplot2::unit(2.5, "pt"),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = HAP_COLS$black),
      axis.title = ggplot2::element_text(size = base_size, colour = HAP_COLS$black),
      panel.grid.major.x = ggplot2::element_line(colour = HAP_COLS$grey_light, linewidth = 0.22),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(4, 6, 4, 4)
    )
}

theme_health_place_map <- function(base_size = 8) {
  ggplot2::theme_void(base_size = base_size, base_family = FIG_FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = HAP_COLS$black, family = FIG_FONT_FAMILY),
      plot.title = ggplot2::element_text(
        size = base_size + 0.1, face = "bold", hjust = 0.5,
        lineheight = 0.92, margin = ggplot2::margin(b = 2.5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1, colour = HAP_COLS$grey_mid,
        hjust = 0.5, lineheight = 1.0
      ),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.key.width = ggplot2::unit(1.2, "cm"),
      legend.key.height = ggplot2::unit(0.22, "cm"),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

theme_lancet <- theme_health_place
theme_lancet_map <- theme_health_place_map

save_hp_plot <- function(filename, plot, width, height,
                         dpi = FIG_DPI,
                         formats = FIG_EXPORT_FORMATS,
                         bg = "white") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  base <- tools::file_path_sans_ext(filename)
  formats <- unique(tolower(formats))
  if ("png" %in% formats) {
    ggplot2::ggsave(paste0(base, ".png"), plot,
                    width = width, height = height, units = "in",
                    dpi = dpi, bg = bg)
  }
  if ("pdf" %in% formats) {
    ggplot2::ggsave(paste0(base, ".pdf"), plot,
                    width = width, height = height, units = "in",
                    device = grDevices::cairo_pdf, bg = bg)
  }
  if ("tiff" %in% formats || "tif" %in% formats) {
    ggplot2::ggsave(paste0(base, ".tiff"), plot,
                    width = width, height = height, units = "in",
                    dpi = dpi, bg = bg, device = "tiff", compression = "lzw")
  }
  invisible(filename)
}
