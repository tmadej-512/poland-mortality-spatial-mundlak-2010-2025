# ============================================================
# 05b_mortality_profiles.R
# Structural determinant profiles at low vs high model-associated
# mortality burden.
#
# Definition:
#   total between-determinant contribution for county i =
#     sum_j [ predicted effect of observed county structural level z_ij
#             relative to the common 0-SD reference ],
#   using ALL 9 between-county determinants from the fitted Spatial
#   Mundlak GAM.
#
# Profiles:
#   LOW  = bottom decile of total between contribution
#   HIGH = top decile of total between contribution
#
# This is a model-associated descriptive synthesis, not a causal
# classification and not a clustering algorithm.
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("drop_geom")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
})

GAM_RDS <- file.path(OUT_BASE, "spatial_mundlak_gam", "gam_results.rds")
if (!file.exists(GAM_RDS)) stop("Brak: ", GAM_RDS)

g <- readRDS(GAM_RDS)

OUT_DIR <- file.path(OUT_BASE, "mortality_profiles")
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
PUB_TAB_DIR <- file.path(OUTPUT_ROOT, "publication_tables")
PUB_FIG_DIR <- file.path(OUTPUT_ROOT, "publication_figures")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== MORTALITY-LINKED STRUCTURAL PROFILES ===\n")
cat(MODEL_LABEL, "|", YEAR_LABEL, "\n")

# ------------------------------------------------------------
# 1. Total structural between-determinant contribution by county
# ------------------------------------------------------------
between_delta <- g$effect_maps %>%
  filter(component == "between", map_type == "delta") %>%
  mutate(
    variable = as.character(variable),
    var_label = label_var(variable)
  )

expected_between <- unique(vapply(g$between_terms, clean_var, character(1)))
missing_between <- setdiff(expected_between, unique(between_delta$variable))
if (length(missing_between)) {
  stop(
    "Brakuje between contribution dla: ",
    paste(missing_between, collapse = ", ")
  )
}

score_by_county <- between_delta %>%
  filter(variable %in% expected_between) %>%
  group_by(code, code_pow, name) %>%
  summarise(
    n_between_terms = n_distinct(variable),
    total_between_contribution =
      if (all(is.finite(value))) sum(value) else NA_real_,
    .groups = "drop"
  )

if (any(score_by_county$n_between_terms != length(expected_between))) {
  stop("Nie każdy powiat ma komplet 9 between contributions.")
}

# Observed/predicted mortality summaries from the same fitted sample.
mortality_by_county <- g$pred_powiat %>%
  select(code, code_pow, name, observed_mean, predicted_mean, residual_mean)

county_scores <- score_by_county %>%
  left_join(mortality_by_county, by = c("code", "code_pow", "name")) %>%
  filter(is.finite(total_between_contribution)) %>%
  mutate(
    burden_decile = dplyr::ntile(total_between_contribution, 10L),
    profile = case_when(
      burden_decile == 1L  ~ "Low model-associated burden",
      burden_decile == 10L ~ "High model-associated burden",
      TRUE ~ "Middle 80%"
    ),
    observed_asmr_decile = dplyr::ntile(observed_mean, 10L),
    observed_profile = case_when(
      observed_asmr_decile == 1L  ~ "Low observed ASMR",
      observed_asmr_decile == 10L ~ "High observed ASMR",
      TRUE ~ "Middle 80%"
    ),
    period = PERIOD,
    year_label = YEAR_LABEL
  )

readr::write_csv(
  county_scores,
  file.path(TAB_DIR, "county_profile_scores.csv")
)

# ------------------------------------------------------------
# 2. Profile-level outcome summary
# ------------------------------------------------------------
profile_outcome_summary <- county_scores %>%
  filter(profile != "Middle 80%") %>%
  group_by(profile) %>%
  summarise(
    n_counties = n(),
    contribution_median = median(total_between_contribution, na.rm = TRUE),
    contribution_q25 = quantile(total_between_contribution, .25, na.rm = TRUE),
    contribution_q75 = quantile(total_between_contribution, .75, na.rm = TRUE),
    observed_asmr_median = median(observed_mean, na.rm = TRUE),
    observed_asmr_q25 = quantile(observed_mean, .25, na.rm = TRUE),
    observed_asmr_q75 = quantile(observed_mean, .75, na.rm = TRUE),
    predicted_asmr_median = median(predicted_mean, na.rm = TRUE),
    predicted_asmr_q25 = quantile(predicted_mean, .25, na.rm = TRUE),
    predicted_asmr_q75 = quantile(predicted_mean, .75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(period = PERIOD, year_label = YEAR_LABEL)

readr::write_csv(
  profile_outcome_summary,
  file.path(TAB_DIR, "profile_outcome_summary.csv")
)

# ------------------------------------------------------------
# 3. County structural raw/z values for all 9 between determinants
# ------------------------------------------------------------
between_values <- purrr::map_dfr(g$between_terms, function(var_z) {
  raw_term <- sub("_z$", "", var_z)
  variable <- clean_var(var_z)

  if (!raw_term %in% names(g$dat_model)) {
    stop("Brak raw between term w dat_model: ", raw_term)
  }
  if (!var_z %in% names(g$dat_model)) {
    stop("Brak standardized between term w dat_model: ", var_z)
  }

  g$dat_model %>%
    select(code, code_pow, name, all_of(raw_term), all_of(var_z)) %>%
    distinct(code, .keep_all = TRUE) %>%
    transmute(
      code, code_pow, name,
      variable = variable,
      var_label = label_var(variable),
      unit = unit_var(variable),
      raw_value = as.numeric(.data[[raw_term]]),
      z_value = as.numeric(.data[[var_z]])
    )
})

profile_values <- between_values %>%
  left_join(
    county_scores %>%
      select(code, code_pow, profile, total_between_contribution,
             observed_mean, predicted_mean),
    by = c("code", "code_pow")
  ) %>%
  filter(profile %in% c(
    "Low model-associated burden",
    "High model-associated burden"
  ))

profile_determinants_long <- profile_values %>%
  group_by(profile, variable, var_label, unit) %>%
  summarise(
    n_counties = sum(is.finite(raw_value) & is.finite(z_value)),
    raw_median = median(raw_value, na.rm = TRUE),
    raw_q25 = quantile(raw_value, .25, na.rm = TRUE),
    raw_q75 = quantile(raw_value, .75, na.rm = TRUE),
    z_median = median(z_value, na.rm = TRUE),
    z_q25 = quantile(z_value, .25, na.rm = TRUE),
    z_q75 = quantile(z_value, .75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(period = PERIOD, year_label = YEAR_LABEL)

readr::write_csv(
  profile_determinants_long,
  file.path(TAB_DIR, "profile_determinants_long.csv")
)

# ------------------------------------------------------------
# 4. Publication-ready table: low vs high burden structural profile
# ------------------------------------------------------------
format_raw <- function(x, unit) {
  if (!is.finite(x)) return(NA_character_)
  if (unit %in% c("PLN", "persons/km²")) {
    return(scales::label_number(
      accuracy = 1, big.mark = " ", decimal.mark = ".", trim = TRUE
    )(x))
  }
  if (unit %in% c("%", "% of economically active", "per 10,000",
                  "µg/m³", "per 1,000", "km")) {
    return(scales::label_number(
      accuracy = .1, big.mark = " ", decimal.mark = ".", trim = TRUE
    )(x))
  }
  scales::label_number(accuracy = .1, big.mark = " ", trim = TRUE)(x)
}

fmt_iqr <- function(med, q25, q75, unit) {
  paste0(
    format_raw(med, unit),
    " [", format_raw(q25, unit), "–", format_raw(q75, unit), "]"
  )
}

profile_table_long <- profile_determinants_long %>%
  rowwise() %>%
  mutate(
    raw_median_iqr = fmt_iqr(raw_median, raw_q25, raw_q75, unit)
  ) %>%
  ungroup()

profile_table_wide <- profile_table_long %>%
  select(variable, var_label, unit, profile, raw_median_iqr, z_median) %>%
  pivot_wider(
    names_from = profile,
    values_from = c(raw_median_iqr, z_median),
    names_sep = "__"
  )

# Stable readable column names.
names(profile_table_wide) <- names(profile_table_wide) %>%
  stringr::str_replace_all("raw_median_iqr__Low model-associated burden",
                           "low_burden_raw_median_IQR") %>%
  stringr::str_replace_all("raw_median_iqr__High model-associated burden",
                           "high_burden_raw_median_IQR") %>%
  stringr::str_replace_all("z_median__Low model-associated burden",
                           "low_burden_median_z") %>%
  stringr::str_replace_all("z_median__High model-associated burden",
                           "high_burden_median_z")

if (all(c("low_burden_median_z", "high_burden_median_z") %in%
        names(profile_table_wide))) {
  profile_table_wide <- profile_table_wide %>%
    mutate(
      high_minus_low_z =
        high_burden_median_z - low_burden_median_z
    )
}

# Preserve manuscript variable order.
profile_table_wide <- tibble(variable = ALL_RAW_VARS) %>%
  left_join(profile_table_wide, by = "variable")

readr::write_csv(
  profile_table_wide,
  file.path(TAB_DIR, "profile_table_publication.csv")
)

# Also copy a period-labelled version to publication_tables.
readr::write_csv(
  profile_table_wide %>% mutate(period = PERIOD, year_label = YEAR_LABEL),
  file.path(PUB_TAB_DIR, paste0("mortality_profile_", PERIOD, ".csv"))
)

# ------------------------------------------------------------
# 5. Diagnostic link with actual mortality
# ------------------------------------------------------------
high_model <- county_scores$burden_decile == 10L
low_model  <- county_scores$burden_decile == 1L
high_obs   <- county_scores$observed_asmr_decile == 10L
low_obs    <- county_scores$observed_asmr_decile == 1L

profile_validation <- tibble(
  period = PERIOD,
  year_label = YEAR_LABEL,
  n_counties = nrow(county_scores),
  spearman_contribution_vs_observed_asmr =
    suppressWarnings(cor(
      county_scores$total_between_contribution,
      county_scores$observed_mean,
      method = "spearman",
      use = "complete.obs"
    )),
  pearson_contribution_vs_observed_asmr =
    suppressWarnings(cor(
      county_scores$total_between_contribution,
      county_scores$observed_mean,
      method = "pearson",
      use = "complete.obs"
    )),
  spearman_contribution_vs_predicted_asmr =
    suppressWarnings(cor(
      county_scores$total_between_contribution,
      county_scores$predicted_mean,
      method = "spearman",
      use = "complete.obs"
    )),
  high_model_n = sum(high_model, na.rm = TRUE),
  high_observed_n = sum(high_obs, na.rm = TRUE),
  high_overlap_n = sum(high_model & high_obs, na.rm = TRUE),
  high_overlap_share_of_model_high =
    sum(high_model & high_obs, na.rm = TRUE) / sum(high_model, na.rm = TRUE),
  low_model_n = sum(low_model, na.rm = TRUE),
  low_observed_n = sum(low_obs, na.rm = TRUE),
  low_overlap_n = sum(low_model & low_obs, na.rm = TRUE),
  low_overlap_share_of_model_low =
    sum(low_model & low_obs, na.rm = TRUE) / sum(low_model, na.rm = TRUE)
)

readr::write_csv(
  profile_validation,
  file.path(TAB_DIR, "profile_validation_against_observed_ASMR.csv")
)

# ------------------------------------------------------------
# 6. Heatmap
# Fill = median standardized structural level (M1-reference SD units).
# Cell label = median raw value.
# Blue = below the M1-reference mean; orange = above.
# ------------------------------------------------------------
profile_levels <- c(
  "Low model-associated burden",
  "High model-associated burden"
)

heat <- profile_determinants_long %>%
  mutate(
    profile = factor(profile, levels = profile_levels),
    variable = factor(variable, levels = rev(ALL_RAW_VARS)),
    row_label = paste0(var_label, "\n", unit)
  )

# Keep row order in manuscript order, top-to-bottom.
row_map <- tibble(
  variable = factor(rev(ALL_RAW_VARS), levels = rev(ALL_RAW_VARS)),
  row_order = seq_along(ALL_RAW_VARS)
)

heat <- heat %>%
  mutate(
    row_label = factor(
      row_label,
      levels = rev(
        paste0(
          vapply(ALL_RAW_VARS, label_var, character(1)),
          "\n",
          vapply(ALL_RAW_VARS, unit_var, character(1))
        )
      )
    ),
    raw_label = purrr::map2_chr(raw_median, unit, format_raw)
  )

z_lim <- max(abs(heat$z_median), na.rm = TRUE)
if (!is.finite(z_lim) || z_lim <= 0) z_lim <- 1
z_lim <- max(1, ceiling(z_lim * 2) / 2)

p_heat <- ggplot(
  heat,
  aes(x = profile, y = row_label, fill = z_median)
) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_text(
    aes(label = raw_label),
    size = 3.0,
    colour = "black",
    family = FIG_FONT_FAMILY
  ) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "#F7F7F7",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-z_lim, z_lim),
    oob = scales::squish,
    name = "Median structural\nlevel (SD)"
  ) +
  scale_x_discrete(
    labels = c(
      "Low model-associated burden" = "Lowest 10%\nmodel-associated burden",
      "High model-associated burden" = "Highest 10%\nmodel-associated burden"
    ),
    position = "top"
  ) +
  labs(
    title = paste0(
      "Structural determinant profiles at low and high model-associated mortality burden | ",
      YEAR_LABEL
    ),
    subtitle = paste0(
      "Groups are the bottom and top deciles of the sum of all 9 between-determinant ",
      "accumulated contributions. Cell text = median in original units; colour = median in M1-reference SD units."
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 9, base_family = FIG_FONT_FAMILY) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 8.2, colour = "black", lineheight = .95),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8.3, colour = HAP_COLS$grey_mid, hjust = 0),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    plot.margin = margin(8, 8, 8, 8)
  )

heat_name <- paste0("mortality_profile_heatmap_", PERIOD, ".png")
save_hp_plot(
  file.path(FIG_DIR, heat_name),
  p_heat,
  width = 8.0,
  height = 6.8
)

# Copy to publication figures with a stable descriptive filename.
save_hp_plot(
  file.path(PUB_FIG_DIR, heat_name),
  p_heat,
  width = 8.0,
  height = 6.8
)

cat("\nMortality-profile outputs:", OUT_DIR, "\n")
cat("  profile_table_publication.csv\n")
cat("  profile_outcome_summary.csv\n")
cat("  profile_validation_against_observed_ASMR.csv\n")
cat("  mortality_profile_heatmap_", PERIOD, ".png/.pdf\n", sep = "")

print(profile_outcome_summary, n = 10)
print(profile_table_wide, n = 20)
print(profile_validation, n = 10)
