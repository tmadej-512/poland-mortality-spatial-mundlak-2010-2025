# ============================================================
# RUN_ALL.R — finalny pipeline v40.5, panel 2010–2025, 3 modele
# ============================================================

# ------------------------------------------------------------
# 1. ŚCIEŻKI — zgodne z podaną lokalizacją nowego panelu
# ------------------------------------------------------------
PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# Robust script directory discovery.
script_candidates <- unique(c(
  file.path(PROJECT_DIR, "R"),
  file.path(PROJECT_DIR, "SSM_v40_6_SVC_REPORTING", "R"),
  file.path(getwd(), "R")
))
required_core <- c("00_config.R", "00_helpers.R")
ok_candidate <- vapply(
  script_candidates,
  function(dd) all(file.exists(file.path(dd, required_core))),
  logical(1)
)
if (!any(ok_candidate)) {
  stop(
    "Nie znaleziono kompletnego katalogu skryptów. Oczekiwano m.in. 00_config.R i 00_helpers.R w:\n",
    paste(script_candidates, collapse = "\n")
  )
}
SCRIPTS_DIR <- script_candidates[which(ok_candidate)[1]]
cat("Scripts: ", normalizePath(SCRIPTS_DIR, winslash = "/", mustWork = TRUE), "\n", sep = "")
DATA_DIR    <- file.path(PROJECT_DIR, "data")

PANEL_RDS <- file.path(DATA_DIR, "finalny_panel_przestrzenny_379_powiatow.rds")
PANEL_CSV <- file.path(DATA_DIR, "finalny_panel_przestrzenny_379_powiatow.csv.gz")

# ------------------------------------------------------------
# 1b. FORMAT EKSPORTU FIGUR — JEDYNY SWITCH DO ZMIANY
# ------------------------------------------------------------
# Dozwolone wartości:
#   "png"  = eksportuj tylko PNG
#   "pdf"  = eksportuj tylko PDF
#   "both" = eksportuj PNG i PDF
# Dotyczy wszystkich figur zapisywanych przez save_hp_plot().
FIG_EXPORT_MODE <- "both"

FIG_EXPORT_MODE <- tolower(trimws(FIG_EXPORT_MODE))
if (!FIG_EXPORT_MODE %in% c("png", "pdf", "both")) {
  stop('FIG_EXPORT_MODE musi być jednym z: "png", "pdf", "both".')
}
cat("Figure export mode: ", FIG_EXPORT_MODE, "\n", sep = "")

# Opcjonalnie ustaw ręcznie pełną ścieżkę do .shp.
# Jeśli NULL, skrypt spróbuje znaleźć *Powiaty.shp rekurencyjnie w PROJECT_DIR.
POWIAT_SHP_PATH <- NULL

setwd(PROJECT_DIR)
cat("Project:", getwd(), "\n")

# ------------------------------------------------------------
# 2. PAKIETY
# ------------------------------------------------------------
required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "tibble", "stringr", "ggplot2",
  "patchwork", "forcats", "scales", "mgcv", "glmnet", "sf", "spdep", "readr"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Brakuje pakietów R: ", paste(missing_pkgs, collapse = ", "),
       "\nZainstaluj je przed uruchomieniem pipeline.")
}

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(readr)
})

# ------------------------------------------------------------
# 3. SKRYPTY / CONFIG
# ------------------------------------------------------------
required_scripts <- c(
  "00_config.R", "00_helpers.R", "00_data_audit.R",
  "01_mundlak_ridge.R", "02_spatial_mundlak_gam.R",
  "02b_svc_uncertainty.R",
  "03_validation.R", "04_perturbation.R",
  "05_figures.R", "05b_mortality_profiles.R",
  "05c_four_period_mortality_profiles.R",
  "05d_svc_reporting.R", "06_collect_results.R"
)
missing_scripts <- file.path(SCRIPTS_DIR, required_scripts)
missing_scripts <- missing_scripts[!file.exists(missing_scripts)]
if (length(missing_scripts)) {
  stop("Brak skryptów w R/:\n", paste(missing_scripts, collapse = "\n"))
}

PERIOD <- "m1_2010_2019"
source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

cat("Final fixed GAM: ", GAM_ARCHITECTURE,
    " | k global/SVC/baseline/year = ", K_GLOBAL, "/", K_SVC, "/", K_BASELINE, "/", K_YEAR, "\n", sep = "")
stopifnot(
  identical(GAM_ARCHITECTURE, "C_all_svc"),
  K_GLOBAL == 5L, K_SVC == 25L, K_BASELINE == 175L, K_YEAR == 9L
)

# ------------------------------------------------------------
# 4. PANEL
# ------------------------------------------------------------
if (file.exists(PANEL_RDS)) {
  cat("Loading RDS:", PANEL_RDS, "\n")
  dat0 <- readRDS(PANEL_RDS)
} else if (file.exists(PANEL_CSV)) {
  cat("RDS missing; loading CSV.GZ:", PANEL_CSV, "\n")
  dat0 <- readr::read_csv(PANEL_CSV, show_col_types = FALSE)
} else {
  stop("Nie znaleziono panelu. Oczekiwano:\n", PANEL_RDS, "\nlub\n", PANEL_CSV)
}

dat0 <- normalize_panel_columns(dat0)
map0 <- attr(dat0, "column_mapping")
if (!is.null(map0)) {
  cat("Column mapping:\n")
  for (nm in names(map0)) cat("  ", nm, " <- ", unname(map0[[nm]]), "\n", sep = "")
}
OUTCOME <- resolve_outcome(dat0, OUTCOME)
OUTCOME_RESOLVED <- OUTCOME
assign("OUTCOME", OUTCOME_RESOLVED, envir = .GlobalEnv)
for (nm in intersect(c(OUTCOME_RESOLVED, ALL_RAW_VARS), names(dat0))) dat0[[nm]] <- to_numeric_safe(dat0[[nm]])
cat("Outcome:", OUTCOME_RESOLVED, "\n")

# ------------------------------------------------------------
# 5. GEOMETRIA / CENTROIDY
# ------------------------------------------------------------
powiaty0 <- NULL
if (!inherits(dat0, "sf")) {
  shp <- POWIAT_SHP_PATH
  if (is.null(shp) || !file.exists(shp)) {
    hits <- list.files(PROJECT_DIR, pattern = "powiaty\\.shp$",
                       recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) shp <- hits[1]
  }
  if (!is.null(shp) && file.exists(shp)) {
    cat("Shapefile:", shp, "\n")
    powiaty0 <- sf::st_read(shp, quiet = TRUE)
  } else {
    stop(
      "Finalny pipeline wymaga geometrii powiatów do Moran's I, perturbacji przestrzennej i figur. ",
      "Nie znaleziono Powiaty.shp. Ustaw POWIAT_SHP_PATH w RUN_ALL.R."
    )
  }
}

sp <- prepare_spatial_objects(dat0, powiaty0)
dat <- sp$dat
powiaty_spdf <- sp$powiaty_spdf
rm(dat0, powiaty0, sp)

validate_raw_panel(dat, OUTCOME_RESOLVED)
dat_input <- dat
cat("Panel:", nrow(dat_input), "rows |", n_distinct(dat_input$code), "counties | years",
    min(dat_input$year, na.rm = TRUE), "–", max(dat_input$year, na.rm = TRUE), "\n")

# ------------------------------------------------------------
# 6. AUDYT DANYCH
# ------------------------------------------------------------
PERIOD <- "m3_2010_2025"
source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)
OUTCOME <- OUTCOME_RESOLVED
dat <- dat_input
source(file.path(SCRIPTS_DIR, "00_data_audit.R"), local = FALSE)
invisible(gc(full = TRUE))

# ------------------------------------------------------------
# 7. FUNKCJA URUCHAMIAJĄCA KROK DLA OKRESU
# ------------------------------------------------------------
run_period_step <- function(script, period_id, label) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(label, "|", period_id, "\n")
  cat(strrep("=", 72), "\n")
  t0 <- proc.time()

  # Każdy kosztowny krok działa w osobnym środowisku. Dzięki temu duże
  # obiekty bam/glmnet/bootstrap nie zostają w .GlobalEnv między modelami.
  e <- new.env(parent = .GlobalEnv)
  e$PERIOD <- period_id
  e$dat <- dat_input
  e$powiaty_spdf <- powiaty_spdf
  e$SCRIPTS_DIR <- SCRIPTS_DIR
  sys.source(file.path(SCRIPTS_DIR, "00_config.R"), envir = e)
  e$OUTCOME <- OUTCOME_RESOLVED
  sys.source(file.path(SCRIPTS_DIR, "00_helpers.R"), envir = e)
  sys.source(file.path(SCRIPTS_DIR, script), envir = e)
  rm(e)
  invisible(gc(full = TRUE))

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("DONE: %s (%.1f min)\n", label, elapsed / 60))
}

periods <- c("m1_2010_2019", "m2_2010_2022", "m3_2010_2025")
t_start <- proc.time()

# Najpierw oba modele dla każdego okresu.
for (pp in periods) run_period_step("01_mundlak_ridge.R", pp, "Mundlak ridge")
for (pp in periods) run_period_step("02_spatial_mundlak_gam.R", pp, "Spatial Mundlak GAM")

# Post-estimation SVC uncertainty; no refit.
for (pp in periods) run_period_step(
  "02b_svc_uncertainty.R",
  pp,
  "SVC uncertainty"
)

# Model-associated structural mortality profiles by M1/M2/M3
# (retained as supplementary post-estimation output).
for (pp in periods) run_period_step(
  "05b_mortality_profiles.R",
  pp,
  "Mortality-linked structural profiles"
)

# PRIMARY four-period profile analysis:
# LOW/HIGH groups are independently defined by observed ASMR,
# then direct M3 term contributions are compared between those groups.
PERIOD <- "m1_2010_2019"
dat <- dat_input
source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
OUTCOME <- OUTCOME_RESOLVED
source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)
source(
  file.path(SCRIPTS_DIR, "05c_four_period_mortality_profiles.R"),
  local = FALSE
)

# Walidacja i perturbacje po dopasowaniu modeli.
for (pp in periods) run_period_step("03_validation.R", pp, "Validation")
for (pp in periods) run_period_step("04_perturbation.R", pp, "Perturbation stability")

# ------------------------------------------------------------
# 8. FIGURY I TABELE ZBIORCZE — po wszystkich 3 modelach
# ------------------------------------------------------------
PERIOD <- "m1_2010_2019"
dat <- dat_input
source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
OUTCOME <- OUTCOME_RESOLVED
source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)
source(file.path(SCRIPTS_DIR, "05_figures.R"), local = FALSE)
source(file.path(SCRIPTS_DIR, "05d_svc_reporting.R"), local = FALSE)
source(file.path(SCRIPTS_DIR, "06_collect_results.R"), local = FALSE)

elapsed <- (proc.time() - t_start)["elapsed"]
cat("\n", strrep("=", 72), "\n", sep = "")
cat("FINAL PIPELINE v40.6 SVC REPORTING DONE\n")
cat(sprintf("Total elapsed: %.1f min\n", elapsed / 60))
cat("Outputs:\n")
cat("  ", file.path(OUTPUT_ROOT, "m1_2010_2019"), "/\n", sep = "")
cat("  ", file.path(OUTPUT_ROOT, "m2_2010_2022"), "/\n", sep = "")
cat("  ", file.path(OUTPUT_ROOT, "m3_2010_2025"), "/\n", sep = "")
cat("  ", file.path(OUTPUT_ROOT, "data_audit"), "/\n", sep = "")
cat("  ", file.path(OUTPUT_ROOT, "publication_figures"), "/\n", sep = "")
cat("  ", file.path(OUTPUT_ROOT, "publication_tables"), "/\n", sep = "")
cat(strrep("=", 72), "\n")
