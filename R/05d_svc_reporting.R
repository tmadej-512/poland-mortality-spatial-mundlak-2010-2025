# ============================================================
# 05d_svc_reporting.R
#
# Reporting/stability layer for SVC results.
#
# This script operationalizes the main reporting ideas discussed in the
# modern SVC visualization/specification literature:
#   - distinguish point estimates from uncertainty;
#   - use a zero-centred diverging palette;
#   - expose the reference contrast explicitly;
#   - summarize local uncertainty rather than mapping estimates alone;
#   - assess stability across alternative/nested analysis windows.
#
# The M1/M2/M3 windows overlap, so "cross-period stability" below is NOT
# independent replication. It is a specification/time-window sensitivity
# diagnostic.
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
})

wrap_fig_text <- function(x, width = 90) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(x)
  stringr::str_wrap(x, width = width)
}

PERIOD_ORDER <- c(
  "m1_2010_2019",
  "m2_2010_2022",
  "m3_2010_2025"
)

OUT <- file.path(OUTPUT_ROOT, "svc_reporting")
TAB_DIR <- file.path(OUT, "tables")
FIG_DIR <- file.path(OUT, "figures")
PUB_TAB_DIR <- file.path(OUTPUT_ROOT, "publication_tables")
PUB_FIG_DIR <- file.path(OUTPUT_ROOT, "publication_figures")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

read_period_maps <- function(pp) {
  f <- file.path(
    OUTPUT_ROOT, pp,
    "spatial_mundlak_gam", "tables",
    "effect_maps_uncertainty.csv"
  )
  if (!file.exists(f)) {
    warning("Brak uncertainty table: ", f)
    return(tibble())
  }
  readr::read_csv(f, show_col_types = FALSE) %>%
    dplyr::mutate(
      period_id = pp,
      period_label = PERIOD_SPECS[[pp]]$year_label
    )
}

maps_all <- map_dfr(PERIOD_ORDER, read_period_maps)

if (!nrow(maps_all)) {
  stop("Brak effect_maps_uncertainty.csv. Uruchom najpierw 02b_svc_uncertainty.R.")
}

# ------------------------------------------------------------
# 1. Cross-period surface stability.
# ------------------------------------------------------------
pair_specs <- tribble(
  ~period_a, ~period_b,
  "m1_2010_2019", "m2_2010_2022",
  "m1_2010_2019", "m3_2010_2025",
  "m2_2010_2022", "m3_2010_2025"
)

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

surface_stability <- pmap_dfr(
  pair_specs,
  function(period_a, period_b) {
    # Resolve labels before entering dplyr verbs. Inside dplyr::mutate(), newly created
    # columns named period_a/period_b would mask these pmap() arguments.
    label_a_value <- PERIOD_SPECS[[period_a]]$year_label
    label_b_value <- PERIOD_SPECS[[period_b]]$year_label

    aa <- maps_all %>%
      dplyr::filter(period_id == period_a) %>%
      dplyr::select(
        code, variable, var_label, component, map_type,
        value_a = value,
        robust_a = robust_sign,
        ci0_a = ci_excludes_zero
      )

    bb <- maps_all %>%
      dplyr::filter(period_id == period_b) %>%
      dplyr::select(
        code, variable, var_label, component, map_type,
        value_b = value,
        robust_b = robust_sign,
        ci0_b = ci_excludes_zero
      )

    aa %>%
      dplyr::inner_join(
        bb,
        by = c("code", "variable", "var_label", "component", "map_type")
      ) %>%
      dplyr::group_by(variable, var_label, component, map_type) %>%
      dplyr::summarise(
        n_common = n(),
        pearson_r = safe_cor(value_a, value_b, "pearson"),
        spearman_r = safe_cor(value_a, value_b, "spearman"),
        raw_sign_agreement =
          mean(sign(value_a) == sign(value_b), na.rm = TRUE),
        robust_share_a = mean(ci0_a, na.rm = TRUE),
        robust_share_b = mean(ci0_b, na.rm = TRUE),
        both_robust_share = mean(ci0_a & ci0_b, na.rm = TRUE),
        robust_sign_agreement_both =
          ifelse(
            sum(ci0_a & ci0_b, na.rm = TRUE) >= 3L,
            mean(
              robust_a[ci0_a & ci0_b] ==
                robust_b[ci0_a & ci0_b],
              na.rm = TRUE
            ),
            NA_real_
          ),
        mean_abs_difference =
          mean(abs(value_b - value_a), na.rm = TRUE),
        rmse_difference =
          sqrt(mean((value_b - value_a)^2, na.rm = TRUE)),
        nrmse_by_pooled_sd = {
          ss <- stats::sd(c(value_a, value_b), na.rm = TRUE)
          ifelse(
            is.finite(ss) && ss > 0,
            sqrt(mean((value_b - value_a)^2, na.rm = TRUE)) / ss,
            NA_real_
          )
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        period_a = .env$period_a,
        period_b = .env$period_b,
        label_a = .env$label_a_value,
        label_b = .env$label_b_value
      )
  }
)

readr::write_csv(
  surface_stability,
  file.path(TAB_DIR, "cross_period_surface_stability.csv")
)
readr::write_csv(
  surface_stability,
  file.path(PUB_TAB_DIR, "SVC_cross_period_surface_stability.csv")
)

# ------------------------------------------------------------
# 2. Uncertainty coverage summary across periods.
# ------------------------------------------------------------
uncertainty_coverage <- maps_all %>%
  dplyr::group_by(
    period_id, period_label,
    variable, var_label, component, map_type
  ) %>%
  dplyr::summarise(
    pct_robust_positive =
      100 * mean(robust_sign == "robust positive", na.rm = TRUE),
    pct_robust_negative =
      100 * mean(robust_sign == "robust negative", na.rm = TRUE),
    pct_uncertain =
      100 * mean(robust_sign == "uncertain", na.rm = TRUE),
    pct_ci_excludes_zero =
      100 * mean(ci_excludes_zero, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  uncertainty_coverage,
  file.path(TAB_DIR, "svc_uncertainty_coverage_all_periods.csv")
)
readr::write_csv(
  uncertainty_coverage,
  file.path(PUB_TAB_DIR, "SVC_uncertainty_coverage_all_periods.csv")
)

# ------------------------------------------------------------
# 3. Compact uncertainty-coverage figure.
#    Separate panels keep marginal and accumulated contrasts distinct.
# ------------------------------------------------------------
pd <- uncertainty_coverage %>%
  dplyr::mutate(
    period_label = factor(
      period_label,
      levels = vapply(
        PERIOD_ORDER,
        function(pp) PERIOD_SPECS[[pp]]$year_label,
        character(1)
      )
    ),
    component = factor(
      component,
      levels = c("between", "within"),
      labels = c("Between-county", "Within-county lag-1")
    ),
    map_type = factor(
      map_type,
      levels = c("deriv", "delta"),
      labels = c("Marginal effect", "Accumulated contribution")
    ),
    var_label = factor(
      var_label,
      levels = rev(label_var(c(MODEL_VARS, ONLY_BETWEEN_VARS)))
    )
  )

p_unc <- ggplot(
  pd,
  aes(
    x = period_label,
    y = var_label,
    fill = pct_ci_excludes_zero
  )
) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(
    aes(label = sprintf("%.0f%%", pct_ci_excludes_zero)),
    size = 2.35,
    family = FIG_FONT_FAMILY
  ) +
  facet_grid(component ~ map_type, scales = "free_y", space = "free_y") +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = HAP_COLS$blue,
    limits = c(0, 100),
    name = "Counties with\n95% CI excluding 0"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(
    base_size = 8.2,
    base_family = FIG_FONT_FAMILY
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 6.6, face = "bold", lineheight = .92, margin = margin(t = 3)),
    axis.text.y = element_text(size = 6.8),
    strip.text = element_text(size = 7.4, face = "bold"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(t = 6, r = 8, b = 6, l = 6),
    legend.position = "right"
  )

save_hp_plot(
  file.path(
    PUB_FIG_DIR,
    "Figure_S_SVC_uncertainty_coverage.png"
  ),
  p_unc,
  width = FIG_W_FULL,
  height = 8.4
)

# ------------------------------------------------------------
# 4. Cross-period stability heatmap.
# ------------------------------------------------------------
ps <- surface_stability %>%
  dplyr::mutate(
    comparison = paste0(label_a, "\nvs\n", label_b),
    component = factor(
      component,
      levels = c("between", "within"),
      labels = c("Between-county", "Within-county lag-1")
    ),
    map_type = factor(
      map_type,
      levels = c("deriv", "delta"),
      labels = c("Marginal effect", "Accumulated contribution")
    ),
    var_label = factor(
      var_label,
      levels = rev(label_var(c(MODEL_VARS, ONLY_BETWEEN_VARS)))
    )
  )

p_stab <- ggplot(
  ps,
  aes(
    x = comparison,
    y = var_label,
    fill = pearson_r
  )
) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(
    aes(label = ifelse(is.finite(pearson_r), sprintf("%.2f", pearson_r), "NA")),
    size = 2.3,
    family = FIG_FONT_FAMILY
  ) +
  facet_grid(component ~ map_type, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = HAP_COLS$vermillion,
    mid = HAP_COLS$white,
    high = HAP_COLS$blue,
    midpoint = 0,
    limits = c(-1, 1),
    name = "County-surface\nPearson r"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(
    base_size = 8.2,
    base_family = FIG_FONT_FAMILY
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 6.1, face = "bold", lineheight = .90, margin = margin(t = 3)),
    axis.text.y = element_text(size = 6.8),
    strip.text = element_text(size = 7.4, face = "bold"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(t = 6, r = 8, b = 8, l = 6),
    legend.position = "right"
  )

save_hp_plot(
  file.path(
    PUB_FIG_DIR,
    "Figure_S_SVC_cross_period_surface_stability.png"
  ),
  p_stab,
  width = FIG_W_FULL,
  height = 8.8
)

# ------------------------------------------------------------
# 5. Companion uncertainty cartograms (robust-sign classification).
#
# These maps answer a different question from the point-estimate maps:
# where does the approximate pointwise 95% coefficient-covariance interval
# support a negative sign, a positive sign, or neither?  They are NOT
# multiplicity-adjusted significance maps and must not be described as such.
#
# One figure is exported for each period x Mundlak component.  Rows are
# determinants and columns are marginal effect vs accumulated contribution.
# ------------------------------------------------------------

ROBUST_SIGN_LEVELS <- c(
  "robust negative",
  "uncertain",
  "robust positive",
  "not estimable"
)
ROBUST_SIGN_COLS <- c(
  "robust negative" = HAP_COLS$blue,
  "uncertain" = "#D9D9D9",
  "robust positive" = HAP_COLS$vermillion,
  "not estimable" = "#FFFFFF"
)

make_robust_sign_map <- function(df, title = NULL, show_legend = FALSE) {
  map_vals <- df %>%
    dplyr::select(code_pow, robust_sign) %>%
    dplyr::distinct(code_pow, .keep_all = TRUE) %>%
    dplyr::mutate(
      code_pow = as.character(code_pow),
      robust_sign = factor(robust_sign, levels = ROBUST_SIGN_LEVELS)
    )

  sf_d <- pow_geom %>%
    dplyr::inner_join(map_vals, by = "code_pow")

  p <- ggplot(sf_d) +
    geom_sf(
      aes(fill = robust_sign),
      colour = "#D9D9D9",
      linewidth = FIG_MAP_COUNTY_BORDER_WIDTH
    ) +
    geom_sf(
      data = country_border,
      fill = NA,
      colour = HAP_COLS$black,
      linewidth = FIG_MAP_COUNTRY_BORDER_WIDTH,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = ROBUST_SIGN_COLS,
      limits = ROBUST_SIGN_LEVELS,
      drop = FALSE,
      name = "Approx. 95% interval",
      labels = c(
        "Robust negative",
        "Includes zero",
        "Robust positive",
        "Not estimable"
      ),
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(colour = "#8C8C8C", linewidth = .25)
      )
    ) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title) +
    theme_void(base_family = FIG_FONT_FAMILY) +
    theme(
      plot.title = element_text(
        size = 7.3, face = "bold", hjust = .5,
        margin = margin(b = 2)
      ),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 7.0, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key.width = grid::unit(.42, "cm"),
      legend.key.height = grid::unit(.30, "cm"),
      legend.key = element_rect(fill = "white", colour = "#A6A6A6", linewidth = .25)
    )

  p
}

make_robust_sign_legend <- function() {
  leg_df <- tibble::tibble(
    x = seq_along(ROBUST_SIGN_LEVELS),
    y = 1,
    robust_sign = factor(ROBUST_SIGN_LEVELS, levels = ROBUST_SIGN_LEVELS)
  )

  p_leg <- ggplot(leg_df, aes(x = x, y = y, fill = robust_sign)) +
    geom_tile(colour = "#8C8C8C", linewidth = .25, show.legend = TRUE) +
    scale_fill_manual(
      values = ROBUST_SIGN_COLS,
      limits = ROBUST_SIGN_LEVELS,
      drop = FALSE,
      name = "Approx. 95% interval",
      labels = c(
        "Robust negative",
        "Includes zero",
        "Robust positive",
        "Not estimable"
      ),
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          colour = "#8C8C8C",
          linewidth = .25,
          alpha = 1
        )
      )
    ) +
    theme_void(base_family = FIG_FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 7.0, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key.width = grid::unit(.42, "cm"),
      legend.key.height = grid::unit(.30, "cm"),
      legend.key = element_rect(fill = "white", colour = "#A6A6A6", linewidth = .25),
      plot.margin = margin(0, 0, 0, 0)
    )

  leg_gt <- ggplotGrob(p_leg)
  guide_idx <- which(vapply(
    leg_gt$grobs,
    function(g) identical(g$name, "guide-box"),
    logical(1)
  ))
  if (!length(guide_idx)) {
    stop("Nie udało się wydobyć legendy robust-sign.")
  }
  patchwork::wrap_elements(full = leg_gt$grobs[[guide_idx[1]]])
}

make_uncertainty_cartogram_grid <- function(pp, component) {
  dd <- maps_all %>%
    dplyr::filter(period_id == pp, component == !!component)

  if (!nrow(dd)) {
    warning("Brak map uncertainty dla ", pp, " / ", component)
    return(invisible(NULL))
  }

  desired_vars <- if (component == "between") {
    c(MODEL_VARS, ONLY_BETWEEN_VARS)
  } else {
    MODEL_VARS
  }
  desired_labels <- label_var(desired_vars)

  rows <- lapply(desired_labels, function(vlab) {
    dm <- dd %>% dplyr::filter(var_label == vlab, map_type == "deriv")
    da <- dd %>% dplyr::filter(var_label == vlab, map_type == "delta")

    if (!nrow(dm) || !nrow(da)) return(NULL)

    pm <- make_robust_sign_map(dm, title = NULL, show_legend = FALSE)
    pa <- make_robust_sign_map(da, title = NULL, show_legend = FALSE)

    label_panel <- ggplot() +
      annotate(
        "text", x = 0, y = .5, label = vlab,
        hjust = 0, size = 2.65, fontface = "bold",
        family = FIG_FONT_FAMILY
      ) +
      xlim(0, 1) + ylim(0, 1) + theme_void()

    label_panel | pm | pa
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(invisible(NULL))

  # Apply consistent widths to every determinant row.
  rows <- lapply(rows, function(x) x + plot_layout(widths = FIG_UNCERTAINTY_GRID_WIDTHS))

  header <- wrap_plots(
    ggplot() + theme_void(),
    ggplot() +
      annotate("text", .5, .5, label = "Marginal effect", fontface = "bold", size = 3) +
      theme_void(),
    ggplot() +
      annotate("text", .5, .5, label = "Accumulated contribution", fontface = "bold", size = 3) +
      theme_void(),
    ncol = 3,
    widths = FIG_UNCERTAINTY_GRID_WIDTHS
  )

  # Dedicated legend is constructed from a dummy dataset containing all
  # robust-sign classes, so colours never disappear when a source map happens
  # to contain only a subset of classes.
  legend_plot <- make_robust_sign_legend()

  body <- wrap_plots(rows, ncol = 1)
  out <- header / body / legend_plot +
    plot_layout(heights = c(.060, 1, .075)) &
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(t = 3, r = 5, b = 4, l = 4)
    )

  fn <- paste0(
    "Figure_S_SVC_uncertainty_cartograms_",
    toupper(sub("_.*$", "", pp)), "_", component, ".png"
  )
  save_hp_plot(
    file.path(PUB_FIG_DIR, fn),
    out,
    width = FIG_W_FULL,
    height = 1.6 + 1.18 * length(rows)
  )
  invisible(out)
}

for (pp in PERIOD_ORDER) {
  make_uncertainty_cartogram_grid(pp, "between")
  make_uncertainty_cartogram_grid(pp, "within")
}

# Machine-readable index for Supplement assembly / provenance.
uncertainty_cartogram_index <- tidyr::expand_grid(
  period_id = PERIOD_ORDER,
  component = c("between", "within")
) %>%
  dplyr::mutate(
    period_label = vapply(
      period_id,
      function(pp) PERIOD_SPECS[[pp]]$year_label,
      character(1)
    ),
    filename = paste0(
      "Figure_S_SVC_uncertainty_cartograms_",
      toupper(sub("_.*$", "", period_id)), "_", component, ".png"
    ),
    definition = paste0(
      "robust sign if approximate pointwise ",
      round(100 * SVC_CI_LEVEL),
      "% coefficient-covariance interval excludes zero; not multiplicity-adjusted"
    )
  )

readr::write_csv(
  uncertainty_cartogram_index,
  file.path(PUB_TAB_DIR, "SVC_uncertainty_cartogram_index.csv")
)

# ------------------------------------------------------------
# 6. Reporting checklist / provenance.
# ------------------------------------------------------------
reporting_checklist <- tribble(
  ~domain, ~pipeline_implementation, ~status,
  "Global comparator",
  "Homogeneous linear Mundlak ridge benchmark fitted for every period",
  "implemented",
  "SVC model specification",
  "Global smooth + baseline spatial smooth + determinant-specific spatially varying coefficient smooth",
  "implemented",
  "Basis adequacy",
  "48-model basis grid plus targeted baseline boundary check; effect stability prioritized over AIC alone",
  "implemented",
  "Residual spatial dependence",
  "Year-specific queen-contiguity Moran's I",
  "implemented",
  "Residual temporal dependence",
  "County residual ACF1",
  "implemented",
  "Out-of-sample behavior",
  "Leakage-controlled leave-one-year-out validation",
  "implemented",
  "Structured robustness",
  "Spatial-temporal residual perturbation refits",
  "implemented",
  "Local uncertainty",
  "lpmatrix contrasts with full coefficient covariance; approximate pointwise intervals",
  "implemented",
  "Uncertainty-aware mapping",
  "Effect maps fade counties whose approximate pointwise 95% interval includes zero",
  "implemented",
  "Reference contrasts",
  "Between contribution = observed structural level vs 0 SD; within contribution = +1 SD lagged deviation vs 0 SD",
  "implemented",
  "Map colour",
  "Zero-centred colourblind-friendly diverging palette",
  "implemented",
  "Cross-period stability",
  "County-surface correlation/sign agreement across nested M1/M2/M3 windows",
  "implemented",
  "Multiplicity caveat",
  "Robust sign is explicitly not described as a multiplicity-adjusted local significance test",
  "implemented",
  "Independent replication",
  "Not available: M1/M2/M3 overlap and are therefore sensitivity windows rather than independent datasets",
  "not applicable / limitation"
)

readr::write_csv(
  reporting_checklist,
  file.path(TAB_DIR, "svc_reporting_checklist.csv")
)

cat("\nSVC reporting outputs DONE:", OUT, "\n")
