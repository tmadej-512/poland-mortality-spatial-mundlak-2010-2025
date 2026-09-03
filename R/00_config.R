# ============================================================
# 00_config.R — wspólna konfiguracja dla wszystkich skryptów
#
# SOURCE ten plik na początku każdego skryptu:
#   source("00_config.R")
#
# Zmień PERIOD żeby przełączyć między okresami:
#   "pre_covid"  → 2010–2019  (main analysis)
#   "full"       → 2010–2023  (includes COVID)
# ============================================================

# ---- Okres analizy ----
# PERIOD mozna ustawic przed sourcowaniem tego pliku
# Jesli nie ustawione z zewnatrz, domyslnie pre_covid
if (!exists("PERIOD")) PERIOD <- "pre_covid"

YEAR_MIN <- 2010

YEAR_MAX <- switch(PERIOD,
  pre_covid = 2019,
  full      = 2023
)

YEAR_LABEL <- switch(PERIOD,
  pre_covid = "2010\u20132019",
  full      = "2010\u20132023"
)

# ---- Outcome ----
OUTCOME <- "asmr_2010_2019std"

# ---- Zmienne ----
MODEL_VARS <- c(
  "inc", "unmp", "docs", "beds",
  "med_exp", "dens", "mig", "pm25"
)
ONLY_BETWEEN_VARS <- c(
  "sec", "uni", "distance_to_city"
)

WITHIN_LAG_N <- 1

# ---- Parametry GAM ----
K_GLOBAL    <- 3
K_SVC       <- 50
K_BASELINE  <- 150
K_YEAR      <- 3
BS_GLOBAL   <- "ts"
BS_SVC      <- "ts"
BS_BASELINE <- "ts"
BS_YEAR     <- "tp"
M_GLOBAL    <- 1
M_SVC       <- 1
M_YEAR      <- 1

FAMILY       <- gaussian()
GAMMA        <- 1.4
BAM_METHOD   <- "fREML"
BAM_DISCRETE <- TRUE
BAM_NTHREADS <- max(1L, parallel::detectCores() - 1L)

SLOPE_EPS <- 0.01
# Backward-compatible alias used by older helper names.
AME_EPS <- SLOPE_EPS
N_GRID  <- 10

# ---- Parametry Mundlak ----
CV_FOLDS   <- 10
LAMBDA_RULE <- "lambda.min"
BOOT_B     <- 1000
BOOT_RHO   <- 0.5

# ---- Output ----
OUT_BASE   <- file.path("output", PERIOD)

# ---- Etykiety ----
VAR_LABELS <- c(
  inc              = "Income / wages",
  unmp             = "Unemployment",
  docs             = "Physicians per 100k",
  beds             = "Hospital beds per 100k",
  med_exp          = "Health expenditure",
  dens             = "Population density",
  mig              = "Net migration",
  pm25             = "PM2.5",
  sec              = "Low education",
  uni              = "Higher education",
  distance_to_city = "Distance to city"
)


# ---- Ustawienia figur pod Health & Place / Elsevier ----
# Domyslnie zapisujemy wersje rastrowa do szybkiego podgladu
# oraz wektorowa PDF do submission/artwork. TIFF mozna wlaczyc
# przez FIG_EXPORT_FORMATS <- c("png", "pdf", "tiff") przed sourcowaniem.
# Submission (peer review): PNG 300 DPI — mniejsze pliki, akceptowane przez Elsevier
# Po akceptacji zmien na: FIG_EXPORT_FORMATS <- c("tiff"); FIG_DPI <- 500
if (!exists("FIG_EXPORT_FORMATS")) FIG_EXPORT_FORMATS <- c("png")
if (!exists("FIG_DPI"))            FIG_DPI <- 300
if (!exists("FIG_FONT_FAMILY")) FIG_FONT_FAMILY <- "sans"

# Paleta Okabe-Ito / color-blind safe. Stale nazwy LANCET zostaja
# jako aliasy, zeby nie przepisywac reszty pipeline.
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

# ---- Paleta figur ----
LANCET <- list(
  # aliasy historyczne: red/navy sa teraz barwami color-blind safe
  red        = HAP_COLS$vermillion,
  navy       = HAP_COLS$blue,
  grey_dark  = HAP_COLS$grey_dark,
  grey_mid   = HAP_COLS$grey_mid,
  grey_light = HAP_COLS$grey_light,
  white      = HAP_COLS$white,
  black      = HAP_COLS$black
)

# ---- Helpery ----
label_var <- function(x) {
  out <- VAR_LABELS[x]
  out[is.na(out)] <- x[is.na(out)]
  as.character(out)
}

clean_var <- function(t) {
  t %>%
    stringr::str_remove("_between_z$") %>%
    stringr::str_remove(
      paste0("_within_lag", WITHIN_LAG_N, "_z$")
    )
}

comp_from_term <- function(t) {
  dplyr::case_when(
    stringr::str_detect(t, "_between_z$") ~ "between",
    stringr::str_detect(
      t, paste0("_within_lag", WITHIN_LAG_N, "_z$")
    ) ~ "within",
    TRUE ~ "other"
  )
}

make_code_pow <- function(code) {
  s <- as.character(code)
  s <- ifelse(nchar(s) == 6L, paste0("0", s), s)
  substring(s, 1L, 4L)
}

impute_median <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) med <- 0
  x[is.na(x)] <- med
  x
}

scale_safe <- function(x) {
  x <- impute_median(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

covid_period_factor <- function(year) {
  factor(
    dplyr::case_when(
      year <= 2019            ~ "pre_covid",
      year %in% c(2020, 2021) ~ "covid_2020_2021",
      year >= 2022            ~ "post_covid",
      TRUE                    ~ NA_character_
    ),
    levels = c("pre_covid", "covid_2020_2021",
               "post_covid")
  )
}

# Dla predykcji — użyj mediany okresu jako referencji
covid_factor_ref <- function(period = PERIOD) {
  ref <- if (period == "pre_covid") "pre_covid"
         else "pre_covid"   # zawsze pre jako referencja
  factor(ref,
         levels = c("pre_covid", "covid_2020_2021",
                    "post_covid"))
}

# Motyw Health & Place / Elsevier-style
# - biale tlo
# - bez pionowych linii i ciezkich ramek
# - font generyczny sans; przy eksporcie PDF system wybierze Arial/Helvetica-like
# - minimalne gridy tylko tam, gdzie pomagaja czytac wartosci
# - domyslna legenda na dole

theme_health_place <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = FIG_FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(
        colour = HAP_COLS$black, family = FIG_FONT_FAMILY),
      plot.title = ggplot2::element_text(
        size = base_size + 1, face = "bold", hjust = 0,
        margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1, colour = HAP_COLS$grey_mid,
        hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.caption = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(
        size = base_size + 2, face = "bold"),
      axis.line = ggplot2::element_line(
        colour = HAP_COLS$grey_dark, linewidth = 0.35),
      axis.ticks = ggplot2::element_line(
        colour = HAP_COLS$grey_dark, linewidth = 0.30),
      axis.ticks.length = ggplot2::unit(2.5, "pt"),
      axis.text = ggplot2::element_text(
        size = base_size - 0.5, colour = HAP_COLS$black),
      axis.title = ggplot2::element_text(
        size = base_size, colour = HAP_COLS$black),
      panel.background = ggplot2::element_rect(
        fill = HAP_COLS$white, colour = NA),
      panel.grid.major.x = ggplot2::element_line(
        colour = HAP_COLS$grey_light, linewidth = 0.22),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(size = base_size - 1),
      legend.text      = ggplot2::element_text(
        size = base_size - 1, colour = HAP_COLS$black),
      legend.key       = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", colour = HAP_COLS$black),
      plot.margin = ggplot2::margin(4, 6, 4, 4)
    )
}

theme_health_place_map <- function(base_size = 8) {
  ggplot2::theme_void(base_size = base_size, base_family = FIG_FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(
        colour = HAP_COLS$black, family = FIG_FONT_FAMILY),
      plot.title = ggplot2::element_text(
        size = base_size + 0.5, face = "bold",
        hjust = 0.5, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1, colour = HAP_COLS$grey_mid,
        hjust = 0.5),
      legend.position   = "bottom",
      legend.text       = ggplot2::element_text(size = base_size - 1),
      legend.key.width  = ggplot2::unit(1.25, "cm"),
      legend.key.height = ggplot2::unit(0.22, "cm"),
      legend.title      = ggplot2::element_text(size = base_size - 1),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

# Zostawiamy stare nazwy jako kompatybilne aliasy.
theme_lancet <- theme_health_place
theme_lancet_map <- theme_health_place_map

# Jednolity eksport figur: PNG do podgladu + PDF wektorowy do submission.
# Opcjonalnie TIFF LZW: ustaw FIG_EXPORT_FORMATS <- c("png", "pdf", "tiff").
save_hp_plot <- function(filename, plot, width, height,
                         dpi = FIG_DPI,
                         formats = FIG_EXPORT_FORMATS,
                         bg = "white") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  base <- tools::file_path_sans_ext(filename)
  formats <- unique(tolower(formats))

  if ("png" %in% formats) {
    ggplot2::ggsave(
      paste0(base, ".png"), plot,
      width = width, height = height, units = "in",
      dpi = dpi, bg = bg
    )
  }

  if ("pdf" %in% formats) {
    ggplot2::ggsave(
      paste0(base, ".pdf"), plot,
      width = width, height = height, units = "in",
      device = grDevices::cairo_pdf, bg = bg
    )
  }

  if ("tiff" %in% formats || "tif" %in% formats) {
    ggplot2::ggsave(
      paste0(base, ".tiff"), plot,
      width = width, height = height, units = "in",
      dpi = dpi, bg = bg,
      device = "tiff", compression = "lzw"
    )
  }

  invisible(filename)
}

# ---- Walidacja zakresu danych ----
validate_period_data <- function(d, period = PERIOD, year_min = YEAR_MIN,
                                 year_max = YEAR_MAX, strict_full = TRUE) {
  if (!"year" %in% names(d))
    stop("Brak kolumny 'year' w obiekcie dat.")

  yrs <- sort(unique(as.integer(d$year)))
  yrs <- yrs[is.finite(yrs)]
  if (length(yrs) == 0)
    stop("Brak poprawnych lat w obiekcie dat.")

  data_min <- min(yrs, na.rm = TRUE)
  data_max <- max(yrs, na.rm = TRUE)

  if (data_min > year_min) {
    warning(sprintf(
      "Dane zaczynają się od %s, a konfiguracja od %s.",
      data_min, year_min
    ))
  }

  if (period == "full" && strict_full && data_max < year_max) {
    stop(sprintf(
      paste0(
        "PERIOD='full' wymaga danych do %s, ale obiekt dat kończy się na %s.\n",
        "Najpierw odśwież/wygeneruj data/zbior_asmr_5yr_age_sex_2010_2019std.rds ",
        "tak, aby zawierał lata 2020–2023. Nazwa pliku oznacza standard ASMR 2010–2019, ",
        "a nie zakres lat analizy.\n",
        "Kontrola: range(dat$year) oraz table(covid_period_factor(dat$year))."
      ),
      year_max, data_max
    ))
  }

  invisible(tibble::tibble(
    period = period, year_min_config = year_min,
    year_max_config = year_max, data_min = data_min, data_max = data_max,
    n_years = dplyr::n_distinct(yrs)
  ))
}

use_covid_control <- function(d, period = PERIOD) {
  period == "full" &&
    "covid_period" %in% names(d) &&
    nlevels(droplevels(d$covid_period)) >= 2
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

cat(sprintf(
  "Config loaded: PERIOD=%s | %s\n",
  PERIOD, YEAR_LABEL
))
