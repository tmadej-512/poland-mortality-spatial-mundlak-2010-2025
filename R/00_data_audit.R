# ============================================================
# 00_data_audit.R — audyt nowego panelu 2010–2025
# ============================================================

if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("build_mundlak_panel")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

OUT_DIR <- file.path(OUTPUT_ROOT, "data_audit")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

validate_raw_panel(dat, OUTCOME)

missing_tbl <- make_missingness_audit(dat, OUTCOME)
readr::write_csv(missing_tbl, file.path(OUT_DIR, "missingness.csv"))

panel_counts <- dat %>%
  drop_geom() %>%
  group_by(year) %>%
  summarise(
    n_rows = n(),
    n_counties = n_distinct(code),
    n_outcome = sum(is.finite(as.numeric(.data[[OUTCOME]]))),
    .groups = "drop"
  )
readr::write_csv(panel_counts, file.path(OUT_DIR, "panel_counts_by_year.csv"))

range_tbl <- dat %>%
  drop_geom() %>%
  select(all_of(c(OUTCOME, ALL_RAW_VARS))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    n = sum(is.finite(as.numeric(value))),
    min = min(as.numeric(value), na.rm = TRUE),
    p02 = quantile(as.numeric(value), .02, na.rm = TRUE),
    p25 = quantile(as.numeric(value), .25, na.rm = TRUE),
    median = median(as.numeric(value), na.rm = TRUE),
    mean = mean(as.numeric(value), na.rm = TRUE),
    p75 = quantile(as.numeric(value), .75, na.rm = TRUE),
    p98 = quantile(as.numeric(value), .98, na.rm = TRUE),
    max = max(as.numeric(value), na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(range_tbl, file.path(OUT_DIR, "variable_ranges.csv"))

period_audit <- purrr::map_dfr(names(PERIOD_SPECS), function(pp) {
  prep <- build_mundlak_panel(dat, period_id = pp)
  spec <- PERIOD_SPECS[[pp]]
  dm <- prep$data %>%
    model_complete_cases(OUTCOME, prep$all_terms, spec$controls, require_xy = TRUE)
  tibble(
    period = pp,
    year_label = spec$year_label,
    raw_rows = sum(dat$year >= spec$year_min & dat$year <= spec$year_max),
    raw_counties = n_distinct(dat$code[dat$year >= spec$year_min & dat$year <= spec$year_max]),
    model_rows = nrow(dm),
    model_counties = n_distinct(dm$code),
    model_year_min = if (nrow(dm)) min(dm$year) else NA_integer_,
    model_year_max = if (nrow(dm)) max(dm$year) else NA_integer_,
    controls = paste(spec$controls, collapse = ";")
  )
})
readr::write_csv(period_audit, file.path(OUT_DIR, "period_model_sample_audit.csv"))

cat("\n=== DATA AUDIT ===\n")
print(period_audit)
cat("\nMissingness:\n")
print(missing_tbl)
cat("\nZapisano:", OUT_DIR, "\n")
