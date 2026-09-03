# ============================================================
# RUN_ALL.R — główny skrypt uruchamiający cały pipeline
#
# Ustaw ścieżkę projektu poniżej, potem uruchom:
#   source("RUN_ALL.R")
# lub z linii poleceń:
#   Rscript RUN_ALL.R
# ============================================================

# ---- Ścieżka projektu ----
PROJECT_DIR <- "C:/Users/panzc/Documents/Instytut/analiza_przestrzenna"

# ---- Załaduj podstawowe pakiety ----
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

# ---- Przejdź do katalogu projektu ----
setwd(PROJECT_DIR)
cat("Katalog roboczy:", getwd(), "\n")

# ---- Ścieżka do skryptów ----
SCRIPTS_DIR <- file.path(PROJECT_DIR, "R")

# ---- Sprawdź czy skrypty istnieją ----
required_scripts <- file.path(SCRIPTS_DIR, c(
  "00_config.R",
  "01_mundlak_ridge.R",
  "02_stvc_gam.R",
  "03_figures.R",
  "04_validation.R",
  "05_perturbation.R",
  "07_descriptive_asmr_context.R"
))

missing <- required_scripts[
  !file.exists(required_scripts)
]
if (length(missing) > 0) {
  cat("Brakujące skrypty:\n")
  cat(paste(" -", missing, collapse="\n"), "\n")
  stop("Skopiuj skrypty do katalogu projektu.")
}

# ---- Wczytaj dane z data/ ----
DATA_DIR <- file.path(PROJECT_DIR, "data")

# Panel powiatowy
df_path <- file.path(DATA_DIR,
                      "zbior_asmr_5yr_age_sex_2010_2019std.rds")
if (!exists("dat")) {
  if (!file.exists(df_path))
    stop("Brak pliku: ", df_path)
  cat("Wczytywanie dat...\n")
  dat <- readRDS(df_path)
  cat("dat wczytane: n =", nrow(dat), "\n")
} else {
  cat("dat już w pamięci: n =", nrow(dat), "\n")
}

# Shapefile powiatów
shp_path <- file.path(DATA_DIR, "Powiaty.shp")
if (!exists("powiaty_spdf")) {
  if (!file.exists(shp_path))
    stop("Brak pliku: ", shp_path)
  cat("Wczytywanie powiaty_spdf...\n")
  suppressPackageStartupMessages(library(sf))
  powiaty_spdf <- sf::st_read(shp_path, quiet = TRUE)
  cat("powiaty_spdf wczytane:",
      nrow(powiaty_spdf), "powiatów\n")
} else {
  cat("powiaty_spdf już w pamięci:",
      nrow(powiaty_spdf), "powiatów\n")
}

# Zachowaj niezmienioną kopię panelu wejściowego. Część skryptów wczytuje własny
# obiekt `dat` z RDS wyników; przed każdym kolejnym krokiem przywracamy panel źródłowy.
dat_input <- dat

cat("dat: n =", nrow(dat_input),
    "| lata:", min(dat_input$year, na.rm=TRUE),
    "-", max(dat_input$year, na.rm=TRUE),
    "| powiaty:", n_distinct(dat_input$code), "\n")

# ============================================================
# PIPELINE
# ============================================================

# Czas całkowity
t_start_all <- proc.time()

run_step <- function(script, period_val,
                      step_name) {
  cat("\n", strrep("=", 60), "\n", sep="")
  cat("KROK:", step_name, "| period =", period_val, "\n")
  cat(strrep("=", 60), "\n")

  script_path <- file.path(SCRIPTS_DIR, script)

  t0 <- proc.time()

  tryCatch({
    # Ustaw PERIOD i re-sourcuj config
    # zeby YEAR_MAX, YEAR_LABEL itp. byly aktualne
    assign("PERIOD", period_val, envir = .GlobalEnv)
    assign("dat", dat_input, envir = .GlobalEnv)
    source(file.path(SCRIPTS_DIR, "00_config.R"),
           local = FALSE)
    validate_period_data(dat)

    source(script_path, local = FALSE)

    elapsed <- (proc.time() - t0)["elapsed"]
    cat(sprintf("\n%s DONE (%.1f min)\n",
                step_name, elapsed / 60))

  }, error = function(e) {
    cat("\nBŁĄD w kroku:", step_name, "\n")
    cat("Komunikat:", e$message, "\n")
    cat("Skrypt:", script, "\n")
    stop(e)
  })
}


# ============================================================
# KROK 0: Opisowy kontekst ASMR — 2010–2023
# ============================================================
run_step("07_descriptive_asmr_context.R", "full",
          "Descriptive ASMR context (levels and temporal change)")

# ============================================================
# KROK 1: Mundlak ridge — pre-COVID
# ============================================================
run_step("01_mundlak_ridge.R", "pre_covid",
          "Mundlak ridge (pre-COVID 2010-2019)")

# ============================================================
# KROK 2: Mundlak ridge — pełny okres
# ============================================================
run_step("01_mundlak_ridge.R", "full",
          "Mundlak ridge (full 2010-2023)")

# ============================================================
# KROK 3: Spatial Mundlak GAM — pre-COVID
# ============================================================
run_step("02_stvc_gam.R", "pre_covid",
          "Spatial Mundlak GAM (pre-COVID 2010-2019)")

# ============================================================
# KROK 4: Spatial Mundlak GAM — pełny okres
# ============================================================
run_step("02_stvc_gam.R", "full",
          "Spatial Mundlak GAM (full 2010-2023)")

# ============================================================
# KROK 5: Wykresy publikacyjne
# ============================================================
run_step("03_figures.R", "pre_covid",
          "Publication figures")

# ============================================================
# KROK 6: Walidacja — pre-COVID
# ============================================================
run_step("04_validation.R", "pre_covid",
          "Validation (pre-COVID)")

# ============================================================
# KROK 7: Walidacja — pełny okres
# ============================================================
run_step("04_validation.R", "full",
          "Validation (full)")

# ============================================================
# KROK 8: Perturbacja — pre-COVID
# ============================================================
run_step("05_perturbation.R", "pre_covid",
          "Perturbation analysis (pre-COVID)")

# ============================================================
# KROK 9: Perturbacja — pełny okres
# ============================================================
run_step("05_perturbation.R", "full",
          "Perturbation analysis (full)")

# ============================================================
# PODSUMOWANIE
# ============================================================

total_elapsed <- (proc.time() - t_start_all)["elapsed"]

cat("\n", strrep("=", 60), "\n", sep="")
cat("PIPELINE DONE\n")
cat(strrep("=", 60), "\n")
cat(sprintf("Łączny czas: %.1f min\n",
            total_elapsed / 60))
cat("\nOutput:\n")
cat("  output/pre_covid/   — pre-COVID (2010-2019)\n")
cat("  output/full/        — pełny okres (2010-2023)\n")
cat("  output/figures/     — wykresy publikacyjne\n")
cat("  output/descriptive_asmr/ — opisowe mapy i tabele ASMR\n")
cat(strrep("=", 60), "\n")
