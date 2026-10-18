# ============================================================
# 05c_four_period_mortality_profiles.R
#
# FINAL FOUR-PERIOD PROFILE ANALYSIS
#
# PRIMARY QUESTION
#   Which fitted determinant contributions distinguish counties with
#   genuinely high vs low observed mortality?
#
# GROUP DEFINITION (independent of the contribution score)
#   Within each reporting window, counties are ranked by their mean
#   OBSERVED ASMR:
#     LOW  = bottom observed-ASMR decile
#     HIGH = top observed-ASMR decile
#
# REPORTING WINDOWS
#   1) 2010–2019
#   2) COVID 2020–2022
#   3) Post-COVID 2023–2025
#   4) All 2010–2025
#
# MODEL CONTRIBUTIONS
#   Contributions are extracted from ONE common M3 Spatial Mundlak GAM
#   (2010–2025; effective model years 2011–2025).
#
# For every determinant term j and county-year it:
#
#   C_jit =
#     prediction(term j at observed z)
#       - prediction(term j at 0 SD),
#
# with all other determinant terms set to 0 SD, year_c=0, and period
# controls=0. The common spatial baseline/intercept cancel in the contrast.
#
# PRIMARY HEATMAP
#   mean(C_j | HIGH observed ASMR) - mean(C_j | LOW observed ASMR)
#   in ASMR / 100,000.
#
# This avoids circularity: HIGH/LOW groups are defined from the outcome,
# not from the same model-contribution score being decomposed.
#
# Raw determinant values of HIGH/LOW observed-ASMR counties are exported
# as median [IQR] in original units.
#
# NOTE ON 2010
#   Observed-ASMR group membership and raw profiles use the full nominal
#   windows (e.g. 2010–2019). Model contributions for windows containing
#   2010 necessarily begin in 2011 because within terms are lagged by one
#   year and M3 is estimated on 2011–2025.
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) {
  source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
}
if (!exists("drop_geom")) {
  source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(patchwork)
})

if (!exists("dat")) {
  stop("Brak pełnego panelu `dat`.")
}

M3_RDS <- file.path(
  OUTPUT_ROOT,
  "m3_2010_2025",
  "spatial_mundlak_gam",
  "gam_results.rds"
)
if (!file.exists(M3_RDS)) {
  stop(
    "Brak finalnego M3 GAM: ", M3_RDS,
    "\nNajpierw dopasuj Model 3."
  )
}

g <- readRDS(M3_RDS)
m <- g$m
dm <- g$dat_model %>% drop_geom()

m3_controls <- g$settings$CONTROL_VARS
m3_terms <- g$all_terms
m3_between <- g$between_terms
m3_within <- g$within_terms

if (length(m3_between) != 9L || length(m3_within) != 8L) {
  stop(
    "Oczekiwano 9 between i 8 within terms w M3; znaleziono ",
    length(m3_between), " i ", length(m3_within), "."
  )
}

OUT_DIR <- file.path(
  OUTPUT_ROOT,
  "four_period_observed_ASMR_model_profiles"
)
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
PUB_TAB_DIR <- file.path(OUTPUT_ROOT, "publication_tables")
PUB_FIG_DIR <- file.path(OUTPUT_ROOT, "publication_figures")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== FOUR-PERIOD OBSERVED-ASMR GROUPS + MODEL CONTRIBUTIONS ===\n")
cat("Grouping: bottom/top decile of observed county mean ASMR.\n")
cat("Contribution source: one common M3 Spatial Mundlak GAM.\n")

PROFILE_WINDOWS <- tibble(
  period_id = c(
    "pre_2010_2019",
    "covid_2020_2022",
    "post_2023_2025",
    "all_2010_2025"
  ),
  period_label = c(
    "2010–2019",
    "COVID 2020–2022",
    "Post-COVID 2023–2025",
    "All 2010–2025"
  ),
  nominal_year_min = c(2010L, 2020L, 2023L, 2010L),
  nominal_year_max = c(2019L, 2022L, 2025L, 2025L),
  effective_year_min = c(2011L, 2020L, 2023L, 2011L),
  effective_year_max = c(2019L, 2022L, 2025L, 2025L)
)

safe_pred <- function(mod, newdata) {
  tryCatch(
    as.numeric(predict(mod, newdata = newdata, type = "response")),
    error = function(e) {
      stop("Prediction failed in profile module: ", e$message)
    }
  )
}

# ------------------------------------------------------------
# 1. Define LOW/HIGH groups ONLY from observed ASMR in nominal windows.
# ------------------------------------------------------------
raw_panel <- dat %>%
  drop_geom() %>%
  mutate(
    year = as.integer(year),
    across(all_of(c(OUTCOME, ALL_RAW_VARS)), to_numeric_safe)
  )

observed_groups <- pmap_dfr(
  PROFILE_WINDOWS,
  function(
    period_id, period_label,
    nominal_year_min, nominal_year_max,
    effective_year_min, effective_year_max
  ) {
    raw_panel %>%
      filter(
        year >= nominal_year_min,
        year <= nominal_year_max
      ) %>%
      group_by(code, code_pow, name) %>%
      summarise(
        observed_asmr = mean(.data[[OUTCOME]], na.rm = TRUE),
        n_observed_years = sum(is.finite(.data[[OUTCOME]])),
        .groups = "drop"
      ) %>%
      mutate(
        period_id = period_id,
        period_label = period_label,
        nominal_year_min = nominal_year_min,
        nominal_year_max = nominal_year_max,
        effective_year_min = effective_year_min,
        effective_year_max = effective_year_max
      )
  }
) %>%
  group_by(period_id) %>%
  arrange(observed_asmr, .by_group = TRUE) %>%
  mutate(
    observed_asmr_decile = ntile(observed_asmr, 10L),
    profile = case_when(
      observed_asmr_decile == 1L  ~ "Low observed ASMR",
      observed_asmr_decile == 10L ~ "High observed ASMR",
      TRUE ~ "Middle 80%"
    )
  ) %>%
  ungroup()

readr::write_csv(
  observed_groups,
  file.path(TAB_DIR, "county_observed_ASMR_groups_four_periods.csv")
)

# ------------------------------------------------------------
# 2. Extract direct term-specific M3 contributions for every county-year.
# ------------------------------------------------------------
nd_ref <- dm
nd_ref$year_c <- 0
for (cc in m3_controls) nd_ref[[cc]] <- 0
for (tt in m3_terms) nd_ref[[tt]] <- 0

cat("Reference prediction...\n")
eta_ref <- safe_pred(m, nd_ref)

term_contributions <- map_dfr(m3_terms, function(tt) {
  cat("  contribution:", tt, "\n")
  nd <- nd_ref
  nd[[tt]] <- dm[[tt]]
  eta_tt <- safe_pred(m, nd)

  tibble(
    code = dm$code,
    code_pow = dm$code_pow,
    name = dm$name,
    year = dm$year,
    term = tt,
    variable = clean_var(tt),
    var_label = label_var(clean_var(tt)),
    component = comp_from_term(tt),
    z_observed = dm[[tt]],
    contribution = eta_tt - eta_ref
  )
})

readr::write_csv(
  term_contributions,
  file.path(TAB_DIR, "M3_term_contributions_by_county_year.csv")
)

# ------------------------------------------------------------
# 3. Average each model term within county x reporting window,
#    using model-estimable years.
# ------------------------------------------------------------
term_window_county <- pmap_dfr(
  PROFILE_WINDOWS,
  function(
    period_id, period_label,
    nominal_year_min, nominal_year_max,
    effective_year_min, effective_year_max
  ) {
    term_contributions %>%
      filter(
        year >= effective_year_min,
        year <= effective_year_max
      ) %>%
      group_by(
        code, code_pow, name,
        term, variable, var_label, component
      ) %>%
      summarise(
        mean_term_contribution = mean(contribution, na.rm = TRUE),
        median_term_contribution = median(contribution, na.rm = TRUE),
        n_model_years = n(),
        .groups = "drop"
      ) %>%
      mutate(
        period_id = period_id,
        period_label = period_label,
        nominal_year_min = nominal_year_min,
        nominal_year_max = nominal_year_max,
        effective_year_min = effective_year_min,
        effective_year_max = effective_year_max
      )
  }
) %>%
  left_join(
    observed_groups %>%
      select(
        code, code_pow, period_id,
        observed_asmr, observed_asmr_decile, profile
      ),
    by = c("code", "code_pow", "period_id")
  )

readr::write_csv(
  term_window_county,
  file.path(TAB_DIR, "county_term_contributions_by_observed_ASMR_group.csv")
)

# ------------------------------------------------------------
# 4. Main contribution contrast: HIGH observed ASMR - LOW observed ASMR.
#
# Use group means for the main decomposition because the GAM is additive
# and the means have a direct decomposition interpretation.
# Medians/IQRs are saved too.
# ------------------------------------------------------------
contribution_profile_summary <- term_window_county %>%
  filter(profile %in% c(
    "Low observed ASMR",
    "High observed ASMR"
  )) %>%
  group_by(
    period_id, period_label,
    component, term, variable, var_label, profile
  ) %>%
  summarise(
    n_counties = n(),
    mean_contribution = mean(mean_term_contribution, na.rm = TRUE),
    median_contribution = median(mean_term_contribution, na.rm = TRUE),
    q25_contribution = quantile(mean_term_contribution, .25, na.rm = TRUE),
    q75_contribution = quantile(mean_term_contribution, .75, na.rm = TRUE),
    .groups = "drop"
  )

contribution_contrast <- contribution_profile_summary %>%
  pivot_wider(
    names_from = profile,
    values_from = c(
      n_counties,
      mean_contribution,
      median_contribution,
      q25_contribution,
      q75_contribution
    ),
    names_sep = "__"
  )

names(contribution_contrast) <- names(contribution_contrast) %>%
  str_replace_all("Low observed ASMR", "LOW_ASMR") %>%
  str_replace_all("High observed ASMR", "HIGH_ASMR")

contribution_contrast <- contribution_contrast %>%
  mutate(
    high_minus_low_mean_contribution =
      mean_contribution__HIGH_ASMR - mean_contribution__LOW_ASMR,
    high_minus_low_median_contribution =
      median_contribution__HIGH_ASMR - median_contribution__LOW_ASMR,
    period_order = match(period_id, PROFILE_WINDOWS$period_id),
    component_order = if_else(component == "between", 1L, 2L),
    variable_order = match(variable, ALL_RAW_VARS)
  ) %>%
  arrange(period_order, component_order, variable_order) %>%
  select(-period_order, -component_order, -variable_order)

readr::write_csv(
  contribution_profile_summary,
  file.path(TAB_DIR, "model_contributions_LOW_HIGH_observed_ASMR_long.csv")
)

readr::write_csv(
  contribution_contrast,
  file.path(TAB_DIR, "model_contribution_HIGH_ASMR_minus_LOW_ASMR.csv")
)

readr::write_csv(
  contribution_contrast,
  file.path(
    PUB_TAB_DIR,
    "Table_model_contribution_HIGH_ASMR_minus_LOW_ASMR_four_periods.csv"
  )
)

# ------------------------------------------------------------
# 5. Determinant-only total burden is retained as a descriptive MODEL score,
#    but it no longer defines the groups.
# ------------------------------------------------------------
burden_it <- term_contributions %>%
  group_by(code, code_pow, name, year) %>%
  summarise(
    between_contribution =
      sum(contribution[component == "between"], na.rm = FALSE),
    within_contribution =
      sum(contribution[component == "within"], na.rm = FALSE),
    total_determinant_contribution =
      between_contribution + within_contribution,
    .groups = "drop"
  )

burden_window <- pmap_dfr(
  PROFILE_WINDOWS,
  function(
    period_id, period_label,
    nominal_year_min, nominal_year_max,
    effective_year_min, effective_year_max
  ) {
    burden_it %>%
      filter(
        year >= effective_year_min,
        year <= effective_year_max
      ) %>%
      group_by(code, code_pow, name) %>%
      summarise(
        determinant_burden =
          mean(total_determinant_contribution, na.rm = TRUE),
        between_burden =
          mean(between_contribution, na.rm = TRUE),
        within_burden =
          mean(within_contribution, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        period_id = period_id,
        period_label = period_label
      )
  }
) %>%
  left_join(
    observed_groups %>%
      select(
        code, code_pow, period_id,
        observed_asmr, profile
      ),
    by = c("code", "code_pow", "period_id")
  )

# Complete-model prediction for context/validation.
pred_full <- safe_pred(m, g$dat_model)
pred_dm <- dm %>%
  transmute(
    code, code_pow, name, year,
    predicted_asmr = pred_full
  )

pred_window <- pmap_dfr(
  PROFILE_WINDOWS,
  function(
    period_id, period_label,
    nominal_year_min, nominal_year_max,
    effective_year_min, effective_year_max
  ) {
    pred_dm %>%
      filter(
        year >= effective_year_min,
        year <= effective_year_max
      ) %>%
      group_by(code, code_pow, name) %>%
      summarise(
        predicted_asmr = mean(predicted_asmr, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(period_id = period_id)
  }
)

burden_window <- burden_window %>%
  left_join(
    pred_window %>% select(code, code_pow, period_id, predicted_asmr),
    by = c("code", "code_pow", "period_id")
  )

readr::write_csv(
  burden_window,
  file.path(TAB_DIR, "observed_ASMR_groups_with_model_burden.csv")
)

profile_outcomes <- burden_window %>%
  filter(profile %in% c("Low observed ASMR", "High observed ASMR")) %>%
  group_by(period_id, period_label, profile) %>%
  summarise(
    n_counties = n(),
    observed_asmr_median = median(observed_asmr, na.rm = TRUE),
    observed_asmr_q25 = quantile(observed_asmr, .25, na.rm = TRUE),
    observed_asmr_q75 = quantile(observed_asmr, .75, na.rm = TRUE),
    predicted_asmr_median = median(predicted_asmr, na.rm = TRUE),
    determinant_burden_mean = mean(determinant_burden, na.rm = TRUE),
    determinant_burden_median = median(determinant_burden, na.rm = TRUE),
    between_burden_mean = mean(between_burden, na.rm = TRUE),
    within_burden_mean = mean(within_burden, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  profile_outcomes,
  file.path(TAB_DIR, "observed_ASMR_profile_outcomes_and_model_burden.csv")
)

# ------------------------------------------------------------
# 6. Raw determinant table for LOW/HIGH OBSERVED-ASMR groups.
#    Uses the exact nominal windows, including 2010 where applicable.
# ------------------------------------------------------------
county_raw_window <- pmap_dfr(
  PROFILE_WINDOWS,
  function(
    period_id, period_label,
    nominal_year_min, nominal_year_max,
    effective_year_min, effective_year_max
  ) {
    raw_panel %>%
      filter(
        year >= nominal_year_min,
        year <= nominal_year_max
      ) %>%
      group_by(code, code_pow, name) %>%
      summarise(
        across(
          all_of(ALL_RAW_VARS),
          ~ mean(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      mutate(
        period_id = period_id,
        period_label = period_label
      )
  }
)

raw_long <- county_raw_window %>%
  pivot_longer(
    cols = all_of(ALL_RAW_VARS),
    names_to = "variable",
    values_to = "raw_value"
  ) %>%
  mutate(
    var_label = label_var(variable),
    unit = unit_var(variable)
  ) %>%
  left_join(
    observed_groups %>%
      select(code, code_pow, period_id, profile),
    by = c("code", "code_pow", "period_id")
  ) %>%
  filter(profile %in% c(
    "Low observed ASMR",
    "High observed ASMR"
  ))

raw_profile_long <- raw_long %>%
  group_by(
    period_id, period_label, profile,
    variable, var_label, unit
  ) %>%
  summarise(
    n_counties = sum(is.finite(raw_value)),
    raw_median = median(raw_value, na.rm = TRUE),
    raw_q25 = quantile(raw_value, .25, na.rm = TRUE),
    raw_q75 = quantile(raw_value, .75, na.rm = TRUE),
    .groups = "drop"
  )

format_raw <- function(x, unit) {
  if (!is.finite(x)) return(NA_character_)

  if (unit %in% c("PLN", "persons/km²")) {
    return(scales::label_number(
      accuracy = 1,
      big.mark = " ",
      decimal.mark = ".",
      trim = TRUE
    )(x))
  }

  scales::label_number(
    accuracy = .1,
    big.mark = " ",
    decimal.mark = ".",
    trim = TRUE
  )(x)
}

fmt_iqr <- function(med, q25, q75, unit) {
  paste0(
    format_raw(med, unit),
    " [",
    format_raw(q25, unit),
    "–",
    format_raw(q75, unit),
    "]"
  )
}

raw_profile_table <- raw_profile_long %>%
  rowwise() %>%
  mutate(
    raw_median_IQR =
      fmt_iqr(raw_median, raw_q25, raw_q75, unit)
  ) %>%
  ungroup() %>%
  select(
    period_id, period_label,
    variable, var_label, unit,
    profile, raw_median_IQR
  ) %>%
  pivot_wider(
    names_from = profile,
    values_from = raw_median_IQR,
    names_sep = "__"
  )

names(raw_profile_table) <- names(raw_profile_table) %>%
  str_replace_all(
    "Low observed ASMR",
    "LOW_observed_ASMR"
  ) %>%
  str_replace_all(
    "High observed ASMR",
    "HIGH_observed_ASMR"
  )

raw_profile_table <- raw_profile_table %>%
  mutate(
    period_order = match(period_id, PROFILE_WINDOWS$period_id),
    variable_order = match(variable, ALL_RAW_VARS)
  ) %>%
  arrange(period_order, variable_order) %>%
  select(-period_order, -variable_order)

readr::write_csv(
  raw_profile_table,
  file.path(TAB_DIR, "raw_determinant_profiles_LOW_HIGH_observed_ASMR.csv")
)

readr::write_csv(
  raw_profile_table,
  file.path(
    PUB_TAB_DIR,
    "Table_raw_determinant_profiles_LOW_HIGH_observed_ASMR_four_periods.csv"
  )
)

# ------------------------------------------------------------
# 7. PRIMARY PUBLICATION HEATMAP A:
#    LOW and HIGH observed-ASMR groups shown SIDE BY SIDE.
#
# Each cell is the mean direct M3 term contribution in that group,
# expressed in ASMR / 100,000.
#
# Interpretation:
#   positive / orange = this term raises modeled ASMR relative to its 0-SD reference;
#   negative / blue   = this term lowers modeled ASMR relative to its 0-SD reference.
#
# The LOW/HIGH grouping is still independent of the model contribution score:
# it comes from observed ASMR deciles.
# ------------------------------------------------------------
period_levels <- PROFILE_WINDOWS$period_label

paired_plot_data <- contribution_profile_summary %>%
  mutate(
    period_label = factor(period_label, levels = period_levels),
    profile = factor(
      profile,
      levels = c("Low observed ASMR", "High observed ASMR")
    ),
    column_label = paste(
      period_label,
      ifelse(
        profile == "Low observed ASMR",
        "LOW ASMR",
        "HIGH ASMR"
      ),
      sep = "\n"
    ),
    column_label = factor(
      column_label,
      levels = as.vector(rbind(
        paste(period_levels, "LOW ASMR", sep = "\n"),
        paste(period_levels, "HIGH ASMR", sep = "\n")
      ))
    ),
    var_label = factor(
      var_label,
      levels = rev(vapply(ALL_RAW_VARS, label_var, character(1)))
    )
  )

pretty_lim <- function(x) {
  if (!is.finite(x) || x <= 0) return(1)
  if (x <= 2) return(ceiling(x * 2) / 2)
  if (x <= 10) return(ceiling(x))
  if (x <= 50) return(ceiling(x / 5) * 5)
  ceiling(x / 10) * 10
}

between_level_lim <- paired_plot_data %>%
  filter(component == "between") %>%
  summarise(v = max(abs(mean_contribution), na.rm = TRUE)) %>%
  pull(v) %>%
  pretty_lim()

within_level_lim <- paired_plot_data %>%
  filter(component == "within") %>%
  summarise(v = max(abs(mean_contribution), na.rm = TRUE)) %>%
  pull(v) %>%
  pretty_lim()

make_paired_heat <- function(comp, lim, title_text) {
  dd <- paired_plot_data %>% filter(component == comp)

  ggplot(
    dd,
    aes(
      x = column_label,
      y = var_label,
      fill = mean_contribution
    )
  ) +
    geom_tile(colour = "white", linewidth = .75) +
    geom_text(
      aes(label = sprintf("%+.1f", mean_contribution)),
      size = 2.45,
      colour = "black",
      family = FIG_FONT_FAMILY
    ) +
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "#F7F7F7",
      high = "#D55E00",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      name = "Mean model\ncontribution\n(ASMR / 100,000)"
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_size = 8.6,
      base_family = FIG_FONT_FAMILY
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        size = 6.7,
        face = "bold",
        colour = "black",
        lineheight = .95
      ),
      axis.text.y = element_text(
        size = 7.4,
        colour = "black"
      ),
      plot.title = element_text(
        size = 10,
        face = "bold",
        hjust = 0
      ),
      legend.position = "right",
      legend.title = element_text(size = 7.2),
      legend.text = element_text(size = 6.7),
      plot.margin = margin(5, 8, 5, 5)
    )
}

p_between_levels <- make_paired_heat(
  "between",
  between_level_lim,
  "A. Between-county determinant contributions"
)

p_within_levels <- make_paired_heat(
  "within",
  within_level_lim,
  "B. Within-county lag-1 determinant contributions"
)

p_paired_main <- (
  p_between_levels /
  p_within_levels
) +
  plot_layout(heights = c(9, 8)) +
  plot_annotation(
    title = "Model contributions in counties with low and high observed mortality",
    subtitle = paste0(
      "LOW/HIGH groups are bottom/top deciles of observed county mean ASMR within each period. ",
      "Cells show the mean direct M3 term contribution in each group (ASMR per 100,000). ",
      "Orange = positive contribution to modeled ASMR; blue = negative contribution. ",
      "For periods containing 2010, contribution averages begin in 2011 because within terms are lagged."
    ),
    theme = theme(
      plot.title = element_text(
        size = 12.2,
        face = "bold",
        hjust = 0,
        family = FIG_FONT_FAMILY
      ),
      plot.subtitle = element_text(
        size = 7.8,
        colour = HAP_COLS$grey_mid,
        hjust = 0,
        family = FIG_FONT_FAMILY
      )
    )
  )

save_hp_plot(
  file.path(
    FIG_DIR,
    "MAIN_model_contributions_LOW_HIGH_observed_ASMR_side_by_side.png"
  ),
  p_paired_main,
  width = 13.2,
  height = 10.2
)

save_hp_plot(
  file.path(
    PUB_FIG_DIR,
    "Figure_model_contributions_LOW_HIGH_observed_ASMR_four_periods.png"
  ),
  p_paired_main,
  width = 13.2,
  height = 10.2
)

# ------------------------------------------------------------
# 7b. SECONDARY PUBLICATION HEATMAP:
#     HIGH minus LOW observed-ASMR contribution contrast.
#     Kept because it gives the cleanest direct group contrast.
# ------------------------------------------------------------

plot_data <- contribution_contrast %>%
  mutate(
    period_label = factor(period_label, levels = period_levels),
    var_label = factor(
      var_label,
      levels = rev(vapply(ALL_RAW_VARS, label_var, character(1)))
    )
  )

between_lim <- plot_data %>%
  filter(component == "between") %>%
  summarise(v = max(abs(high_minus_low_mean_contribution), na.rm = TRUE)) %>%
  pull(v) %>%
  pretty_lim()

within_lim <- plot_data %>%
  filter(component == "within") %>%
  summarise(v = max(abs(high_minus_low_mean_contribution), na.rm = TRUE)) %>%
  pull(v) %>%
  pretty_lim()

make_heat <- function(comp, lim, title_text) {
  dd <- plot_data %>% filter(component == comp)

  ggplot(
    dd,
    aes(
      x = period_label,
      y = var_label,
      fill = high_minus_low_mean_contribution
    )
  ) +
    geom_tile(colour = "white", linewidth = .8) +
    geom_text(
      aes(label = sprintf("%+.1f", high_minus_low_mean_contribution)),
      size = 2.9,
      colour = "black",
      family = FIG_FONT_FAMILY
    ) +
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "#F7F7F7",
      high = "#D55E00",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      name = "HIGH ASMR − LOW ASMR\nmodel contribution\n(ASMR / 100,000)"
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_size = 8.8,
      base_family = FIG_FONT_FAMILY
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        size = 7.8,
        face = "bold",
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 7.5,
        colour = "black"
      ),
      plot.title = element_text(
        size = 10,
        face = "bold",
        hjust = 0
      ),
      legend.position = "right",
      legend.title = element_text(size = 7.3),
      legend.text = element_text(size = 6.8),
      plot.margin = margin(5, 8, 5, 5)
    )
}

p_between <- make_heat(
  "between",
  between_lim,
  "A. Between-county determinant contributions"
)

p_within <- make_heat(
  "within",
  within_lim,
  "B. Within-county lag-1 determinant contributions"
)

p_main <- (
  p_between /
  p_within
) +
  plot_layout(heights = c(9, 8)) +
  plot_annotation(
    title = "Model contributions distinguishing counties with high vs low observed mortality",
    subtitle = paste0(
      "Counties are defined independently by observed ASMR: top vs bottom decile within each period. ",
      "Cells show the difference in mean M3 term contribution (HIGH observed ASMR − LOW observed ASMR), ",
      "in ASMR per 100,000. Orange = more positive / less protective contribution in high-mortality counties; ",
      "blue = more negative / more protective contribution. ",
      "For windows containing 2010, contribution averages begin in 2011 because within terms are lagged."
    ),
    theme = theme(
      plot.title = element_text(
        size = 12.2,
        face = "bold",
        hjust = 0,
        family = FIG_FONT_FAMILY
      ),
      plot.subtitle = element_text(
        size = 7.8,
        colour = HAP_COLS$grey_mid,
        hjust = 0,
        family = FIG_FONT_FAMILY
      )
    )
  )

save_hp_plot(
  file.path(
    FIG_DIR,
    "SECONDARY_model_contributions_HIGH_minus_LOW_observed_ASMR.png"
  ),
  p_main,
  width = 9.4,
  height = 10.2
)

save_hp_plot(
  file.path(
    PUB_FIG_DIR,
    "Figure_model_contributions_HIGH_vs_LOW_observed_ASMR_four_periods.png"
  ),
  p_main,
  width = 9.4,
  height = 10.2
)

# ------------------------------------------------------------
# 8. Supplement: model-burden-defined decomposition, retained only to show
#    how the determinant-only score itself is constructed.
# ------------------------------------------------------------
supp_groups <- burden_window %>%
  group_by(period_id) %>%
  mutate(
    burden_decile = ntile(determinant_burden, 10L),
    burden_profile = case_when(
      burden_decile == 1L  ~ "Low model burden",
      burden_decile == 10L ~ "High model burden",
      TRUE ~ "Middle 80%"
    )
  ) %>%
  ungroup() %>%
  select(code, code_pow, period_id, burden_profile)

supp_contrast <- term_window_county %>%
  select(
    -profile, -observed_asmr, -observed_asmr_decile
  ) %>%
  left_join(
    supp_groups,
    by = c("code", "code_pow", "period_id")
  ) %>%
  filter(burden_profile %in% c(
    "Low model burden",
    "High model burden"
  )) %>%
  group_by(
    period_id, period_label,
    component, term, variable, var_label,
    burden_profile
  ) %>%
  summarise(
    mean_contribution = mean(mean_term_contribution, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = burden_profile,
    values_from = mean_contribution
  ) %>%
  mutate(
    high_model_burden_minus_low_model_burden =
      .data[["High model burden"]] -
      .data[["Low model burden"]]
  )

readr::write_csv(
  supp_contrast,
  file.path(
    TAB_DIR,
    "SUPP_model_burden_defined_contribution_decomposition.csv"
  )
)

# LOW/HIGH MODEL-BURDEN groups shown side by side, as in the earlier profile
# presentation. This is retained for decomposition/interpretation, not as an
# independent validation, because group membership is defined by the same
# determinant burden being summarized.
supp_levels_long <- term_window_county %>%
  select(
    -profile, -observed_asmr, -observed_asmr_decile
  ) %>%
  left_join(
    supp_groups,
    by = c("code", "code_pow", "period_id")
  ) %>%
  filter(burden_profile %in% c(
    "Low model burden",
    "High model burden"
  )) %>%
  group_by(
    period_id, period_label,
    component, term, variable, var_label,
    burden_profile
  ) %>%
  summarise(
    mean_contribution = mean(mean_term_contribution, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    period_label = factor(period_label, levels = period_levels),
    burden_profile = factor(
      burden_profile,
      levels = c("Low model burden", "High model burden")
    ),
    column_label = paste(
      period_label,
      ifelse(
        burden_profile == "Low model burden",
        "LOW burden",
        "HIGH burden"
      ),
      sep = "\n"
    ),
    column_label = factor(
      column_label,
      levels = as.vector(rbind(
        paste(period_levels, "LOW burden", sep = "\n"),
        paste(period_levels, "HIGH burden", sep = "\n")
      ))
    ),
    var_label = factor(
      var_label,
      levels = rev(vapply(ALL_RAW_VARS, label_var, character(1)))
    )
  )

supp_between_lim <- supp_levels_long %>%
  filter(component == "between") %>%
  summarise(v = max(abs(mean_contribution), na.rm = TRUE)) %>%
  pull(v) %>%
  pretty_lim()

supp_within_lim <- supp_levels_long %>%
  filter(component == "within") %>%
  summarise(v = max(abs(mean_contribution), na.rm = TRUE)) %>%
  pull(v) %>%
  pretty_lim()

make_supp_burden_heat <- function(comp, lim, title_text) {
  dd <- supp_levels_long %>% filter(component == comp)

  ggplot(
    dd,
    aes(
      x = column_label,
      y = var_label,
      fill = mean_contribution
    )
  ) +
    geom_tile(colour = "white", linewidth = .75) +
    geom_text(
      aes(label = sprintf("%+.1f", mean_contribution)),
      size = 2.45,
      colour = "black",
      family = FIG_FONT_FAMILY
    ) +
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "#F7F7F7",
      high = "#D55E00",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      name = "Mean model\ncontribution\n(ASMR / 100,000)"
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_size = 8.6,
      base_family = FIG_FONT_FAMILY
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        size = 6.7,
        face = "bold",
        colour = "black",
        lineheight = .95
      ),
      axis.text.y = element_text(
        size = 7.4,
        colour = "black"
      ),
      plot.title = element_text(
        size = 10,
        face = "bold",
        hjust = 0
      ),
      legend.position = "right",
      legend.title = element_text(size = 7.2),
      legend.text = element_text(size = 6.7)
    )
}

p_supp_burden_between <- make_supp_burden_heat(
  "between",
  supp_between_lim,
  "A. Between-county determinant contributions"
)

p_supp_burden_within <- make_supp_burden_heat(
  "within",
  supp_within_lim,
  "B. Within-county lag-1 determinant contributions"
)

p_supp_burden <- (
  p_supp_burden_between /
  p_supp_burden_within
) +
  plot_layout(heights = c(9, 8)) +
  plot_annotation(
    title = "Model contributions in LOW vs HIGH model-burden counties",
    subtitle = paste0(
      "Groups are bottom/top deciles of the determinant-only M3 burden score. ",
      "Cells show mean model term contributions in each group. ",
      "This figure decomposes the model-defined burden itself and is therefore supplementary."
    ),
    theme = theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0,
        family = FIG_FONT_FAMILY
      ),
      plot.subtitle = element_text(
        size = 7.8,
        colour = HAP_COLS$grey_mid,
        hjust = 0,
        family = FIG_FONT_FAMILY
      )
    )
  )

save_hp_plot(
  file.path(
    FIG_DIR,
    "SUPP_model_contributions_LOW_HIGH_model_burden_side_by_side.png"
  ),
  p_supp_burden,
  width = 13.2,
  height = 10.2
)

# ------------------------------------------------------------
# 9. Methods record.
# ------------------------------------------------------------
methods_record <- tibble(
  primary_group_definition =
    "bottom/top decile of observed county mean ASMR within each nominal reporting window",
  contribution_model =
    "single common M3 2010-2025 Spatial Mundlak GAM",
  contribution_statistic =
    "mean term contribution in HIGH observed-ASMR counties minus mean term contribution in LOW observed-ASMR counties",
  contribution_units =
    "ASMR per 100,000",
  n_between_terms = length(m3_between),
  n_within_terms = length(m3_within),
  contribution_reference =
    "observed standardized term value versus 0 SD; all other determinant terms at 0 SD",
  excluded_from_term_contrasts =
    "common intercept/baseline spatial surface cancels; year_c and period dummies held at 0",
  nominal_group_windows =
    "2010-2019; 2020-2022; 2023-2025; 2010-2025",
  effective_contribution_windows =
    "2011-2019; 2020-2022; 2023-2025; 2011-2025",
  reason_for_effective_window_shift =
    "M3 within terms use lag 1, so no 2010 model row exists",
  raw_profile_table =
    "median [IQR] raw determinant values in observed-ASMR-defined LOW/HIGH groups",
  supplementary_decomposition =
    "model-burden-defined HIGH-minus-LOW contribution decomposition retained only as supplement"
)

readr::write_csv(
  methods_record,
  file.path(TAB_DIR, "profile_methods_record.csv")
)

cat("\nFOUR-PERIOD observed-ASMR/model-contribution profiles DONE:\n")
cat("  ", file.path(
  FIG_DIR,
  "MAIN_model_contributions_LOW_HIGH_observed_ASMR_side_by_side.png"
), "\n", sep = "")
cat("  ", file.path(
  FIG_DIR,
  "SECONDARY_model_contributions_HIGH_minus_LOW_observed_ASMR.png"
), "\n", sep = "")
cat("  ", file.path(
  TAB_DIR,
  "model_contribution_HIGH_ASMR_minus_LOW_ASMR.csv"
), "\n", sep = "")
cat("  ", file.path(
  TAB_DIR,
  "raw_determinant_profiles_LOW_HIGH_observed_ASMR.csv"
), "\n", sep = "")
cat("  ", file.path(
  TAB_DIR,
  "observed_ASMR_profile_outcomes_and_model_burden.csv"
), "\n", sep = "")

print(profile_outcomes, n = 20)
print(contribution_contrast, n = 50)
print(raw_profile_table, n = 50)
