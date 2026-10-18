# ============================================================
# 05_figures.R
# PLOTFIX: national outline uses full shapefile with polygon holes stripped;
# effect-grid legends are extracted guide grobs placed directly below maps. — publication figures, ALL estimated effects
# SSM final pipeline v40.7 — comprehensive figure + uncertainty reporting patch
#
# Main rules:
#   * plot ALL 9 between and ALL 8 within-lag1 effects — no perturbation/top-N filtering;
#   * curve x ticks are one model SD apart but LABELLED in original units;
#   * y scales are locked across periods within component (between / within);
#   * map legends use sparse, non-overlapping ticks and appear once per effect-grid column;
#   * every county map has county boundaries + black Poland outline;
#   * within accumulated contribution is defined as a standardized +1 SD vs 0 SD
#     scenario (signed observed within deviations average to ~0 by construction);
#   * dominant-contribution maps use ALL determinants, not stability filtering;
#   * when SVC uncertainty is available, effect maps fade counties whose approximate
#     95% pointwise interval includes zero (estimate hue is unchanged).
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("drop_geom")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(patchwork); library(forcats); library(scales)
  library(sf); library(readr)
})

wrap_fig_text <- function(x, width = 90) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(x)
  stringr::str_wrap(x, width = width)
}

raw_panel_title <- function(variable, unit = NULL, suffix = NULL, width = 20) {
  lbl <- stringr::str_wrap(label_var(variable), width = width)
  parts <- c(lbl, suffix, unit)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  paste(parts, collapse = "\n")
}

PUB_DIR <- file.path(OUTPUT_ROOT, "publication_figures")
dir.create(PUB_DIR, recursive = TRUE, showWarnings = FALSE)

PERIOD_ORDER <- c("m1_2010_2019", "m2_2010_2022", "m3_2010_2025")

load_period <- function(pp) {
  gam_path <- file.path(OUTPUT_ROOT, pp, "spatial_mundlak_gam", "gam_results.rds")
  ridge_path <- file.path(OUTPUT_ROOT, pp, "mundlak_ridge", "mundlak_results.rds")
  pert_path <- file.path(OUTPUT_ROOT, pp, "perturbation", "perturbation_results.rds")
  val_path <- file.path(OUTPUT_ROOT, pp, "validation", "validation_results.rds")
  if (!file.exists(gam_path)) stop("Brak: ", gam_path)
  list(
    gam = readRDS(gam_path),
    ridge = if (file.exists(ridge_path)) readRDS(ridge_path) else NULL,
    perturb = if (file.exists(pert_path)) readRDS(pert_path) else NULL,
    validation = if (file.exists(val_path)) readRDS(val_path) else NULL
  )
}

R <- setNames(lapply(PERIOD_ORDER, load_period), PERIOD_ORDER)

# Return only determinant variables that are actually present in the saved GAM
# object for the requested Mundlak layer.  Keeping this helper local to the
# figure script makes RUN_FIGURES_ONLY.R independent of model-fitting scripts.
active_vars <- function(bundle, component = c("between", "within")) {
  component <- match.arg(component)
  g <- bundle$gam
  if (is.null(g$curves_all) || !all(c("component", "variable") %in% names(g$curves_all))) {
    stop("Saved gam_results.rds does not contain curves_all with component/variable columns.")
  }

  wanted <- if (component == "between") {
    c(MODEL_VARS, ONLY_BETWEEN_VARS)
  } else {
    MODEL_VARS
  }

  available <- g$curves_all %>%
    dplyr::filter(.data$component == component) %>%
    dplyr::pull(.data$variable) %>%
    as.character() %>%
    unique()

  vars <- wanted[wanted %in% available]
  missing <- setdiff(wanted, vars)
  if (length(missing)) {
    warning(
      "Omitting variables not present in saved ", component,
      " curves: ", paste(missing, collapse = ", ")
    )
  }
  vars
}

# ------------------------------------------------------------
# Geometry: canonical county layer + single black country outline
# ------------------------------------------------------------
if (is.null(powiaty_spdf)) stop("05_figures.R wymaga powiaty_spdf.")
pow_geom <- sf::st_as_sf(powiaty_spdf)
if (!POW_CODE_COL %in% names(pow_geom)) {
  cc <- intersect(c("JPT_KOD_JE", "code_pow"), names(pow_geom))[1]
  if (is.na(cc)) stop("Brak code_pow/JPT_KOD_JE w geometrii.")
  names(pow_geom)[names(pow_geom) == cc] <- POW_CODE_COL
}
pow_geom[[POW_CODE_COL]] <- as.character(pow_geom[[POW_CODE_COL]])
pow_geom <- pow_geom %>% dplyr::select(dplyr::all_of(POW_CODE_COL), geometry)

# All downstream joins use code_pow explicitly.
if (POW_CODE_COL != "code_pow") {
  names(pow_geom)[names(pow_geom) == POW_CODE_COL] <- "code_pow"
}
pow_geom$code_pow <- as.character(pow_geom$code_pow)

# Keep a full-geography copy for the OUTER national outline.
# The fill layer is still restricted to the 379 harmonised analytical counties.
pow_geom_full <- pow_geom

analysis_codes <- dat %>%
  drop_geom() %>%
  dplyr::distinct(code_pow) %>%
  dplyr::mutate(code_pow = as.character(code_pow))

pow_geom <- pow_geom %>%
  dplyr::semi_join(analysis_codes, by = "code_pow")

# Remove interior polygon rings from a dissolved national geometry.
# This prevents city-county / harmonisation holes from being redrawn as
# small black "dots" or loops inside Poland by the country-outline layer.
strip_polygon_holes <- function(g) {
  if (inherits(g, "POLYGON")) {
    return(sf::st_polygon(list(g[[1]])))
  }

  if (inherits(g, "MULTIPOLYGON")) {
    return(sf::st_multipolygon(
      lapply(g, function(poly) list(poly[[1]]))
    ))
  }

  if (inherits(g, "GEOMETRYCOLLECTION")) {
    polys <- list()
    for (part in g) {
      if (inherits(part, "POLYGON")) {
        polys[[length(polys) + 1L]] <- list(part[[1]])
      } else if (inherits(part, "MULTIPOLYGON")) {
        for (poly in part) {
          polys[[length(polys) + 1L]] <- list(poly[[1]])
        }
      }
    }
    if (length(polys)) return(sf::st_multipolygon(polys))
  }

  g
}

country_union_full <- sf::st_union(
  sf::st_geometry(sf::st_make_valid(pow_geom_full))
)

country_outer_sfg <- strip_polygon_holes(country_union_full[[1]])

country_border <- sf::st_sf(
  geometry = sf::st_sfc(
    country_outer_sfg,
    crs = sf::st_crs(pow_geom_full)
  )
)

compact_lab <- function(x) {
  x <- as.numeric(x)
  finite <- abs(x[is.finite(x)])
  mx <- if (length(finite)) max(finite) else 1
  acc <- dplyr::case_when(
    mx >= 100 ~ 1,
    mx >= 10  ~ 0.1,
    mx >= 1   ~ 0.1,
    mx >= .1  ~ 0.01,
    TRUE      ~ 0.001
  )
  scales::label_number(accuracy = acc, big.mark = " ", trim = TRUE)(x)
}

raw_axis_lab <- function(x, unit = "") {
  x <- as.numeric(x)
  if (unit == "PLN") {
    return(scales::label_number(accuracy = 1, big.mark = " ", trim = TRUE)(x))
  }
  if (unit == "persons/km²") {
    return(scales::label_number(accuracy = 1, big.mark = " ", trim = TRUE)(x))
  }
  if (unit %in% c("%", "% of economically active", "per 10,000", "µg/m³", "per 1,000", "km")) {
    return(scales::label_number(accuracy = 0.1, big.mark = " ", trim = TRUE)(x))
  }
  scales::label_number(accuracy = 0.1, big.mark = " ", trim = TRUE)(x)
}

three_scale_breaks <- function(limits, symmetric = FALSE) {
  limits <- as.numeric(limits)
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[1] == limits[2]) {
    return(c(-1, 0, 1))
  }
  if (symmetric) {
    lim <- max(abs(limits))
    return(c(-lim, 0, lim))
  }
  c(limits[1], mean(limits), limits[2])
}

# High-contrast sequential palette for descriptive/raw maps.
SEQ_COLS <- c("#F7FBFF", "#DCEAF6", "#A9C9E5", "#6BAED6", "#3182BD", "#08519C", "#08306B")
# High-contrast, colour-blind-friendly diverging palette for effect cartograms.
# The scale remains linear in the data; stronger intermediate colours improve
# discrimination of moderate effects without changing the underlying estimates.
DIV_COLS <- c("#005A9C", "#2C7FB8", "#9ECAE1", HAP_COLS$white,
              "#F6C28B", "#E6783A", "#B63C1A")

robust_abs_limit <- function(x, prob = .98) {
  z <- abs(x[is.finite(x)])
  if (!length(z)) return(1)
  lim <- as.numeric(stats::quantile(z, prob, na.rm = TRUE, names = FALSE))
  if (!is.finite(lim) || lim <= 0) lim <- max(z, na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  lim
}

robust_seq_limits <- function(x, probs = c(.02, .98)) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  q <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  if (length(unique(q)) < 2L || any(!is.finite(q))) q <- range(x, na.rm = TRUE)
  if (q[1] == q[2]) q <- q + c(-.5, .5)
  q
}

nice_symmetric_scale <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(list(lim = 1, breaks = c(-1, 0, 1)))
  m <- max(abs(x), na.rm = TRUE) * 1.03
  if (!is.finite(m) || m <= 0) m <- 1
  br <- pretty(c(-m, m), n = 5)
  lim <- max(abs(br), na.rm = TRUE)
  br <- br[br >= -lim & br <= lim]
  list(lim = lim, breaks = br)
}

# One locked curve-y scale for all periods within a component.
# Between and within remain separate because their effect magnitudes differ strongly.
shared_curve_y_info <- function(component) {
  yy <- purrr::map_dfr(R, function(bundle) {
    bundle$gam$curves_all %>%
      dplyr::filter(component == !!component) %>%
      dplyr::select(effect_p10, effect_p90, effect_q25, effect_q75, effect_median)
  })
  nice_symmetric_scale(unlist(yy, use.names = FALSE))
}

# Convert model-z positions to original units using the stored M1-reference scaler.
# Tick positions stay exactly 1 SD apart in z-space.
term_axis_spec <- function(bundle, var) {
  g <- bundle$gam
  raw_term <- sub("_z$", "", var)
  st <- g$prep$scaler_table
  rr <- st %>% dplyr::filter(raw_term == !!raw_term) %>% dplyr::slice_head(n = 1)

  dc <- g$curves_all %>% dplyr::filter(var == !!var)
  zr <- range(dc$var_value, na.rm = TRUE)
  zb <- seq(ceiling(zr[1]), floor(zr[2]), by = 1)
  if (!length(zb)) zb <- pretty(zr, n = 3)

  base_var <- clean_var(var)
  unit <- unit_var(base_var)

  if (nrow(rr) == 1L && is.finite(rr$center[1]) && is.finite(rr$scale[1])) {
    raw_breaks <- rr$center[1] + zb * rr$scale[1]
    labs <- raw_axis_lab(raw_breaks, unit)
  } else {
    labs <- compact_lab(zb)
  }

  axis_title <- if (comp_from_term(var) == "between") {
    paste0("County mean (", unit, "; tick spacing = 1 SD)")
  } else {
    paste0("Lag-1 within-county deviation (", unit, "; tick spacing = 1 SD)")
  }

  list(breaks = zb, labels = labs, title = axis_title)
}

# Rebuild only the problematic WITHIN contribution map at figure time.
# Averaging signed observed within deviations over time is near zero by construction.
# For visualization we therefore use the standardized accumulated contrast:
# predicted outcome at +1 SD lagged within deviation minus predicted outcome at 0 SD.
# No GAM refit is needed; this uses the already fitted model.
plot_effect_maps <- function(bundle, component) {
  g <- bundle$gam
  stored <- g$effect_maps %>% dplyr::filter(component == !!component)
  if (component != "within") return(stored)

  # New pipeline stores the standardized within contribution directly.
  if ("contribution_definition" %in% names(stored)) {
    def <- stored$contribution_definition[stored$map_type == "delta"]
    if (length(def) && any(grepl("\\+1 SD", def))) return(stored)
  }

  # Backward compatibility for GAM objects fitted before this figure patch.
  deriv <- stored %>% dplyr::filter(map_type == "deriv")
  controls <- if (!is.null(g$settings$CONTROL_VARS)) g$settings$CONTROL_VARS else CONTROL_VARS
  ref <- make_reference_counties(g$dat_model, g$all_terms, controls)

  delta_plus1 <- purrr::map_dfr(g$within_terms, function(var) {
    nd0 <- ref
    nd1 <- ref
    nd0[[var]] <- 0
    nd1[[var]] <- 1
    p0 <- as.numeric(predict(g$m, newdata = nd0, type = "response"))
    p1 <- as.numeric(predict(g$m, newdata = nd1, type = "response"))

    ref %>% dplyr::transmute(
      code, code_pow, name, x, y,
      v_actual = 1,
      value = p1 - p0,
      map_type = "delta",
      var = var,
      variable = clean_var(var),
      var_label = label_var(variable),
      component = "within",
      period = g$period,
      year_label = g$year_label,
      contribution_definition = "+1 SD lag-1 within deviation vs 0 SD"
    )
  })

  dplyr::bind_rows(deriv, delta_plus1)
}

# ------------------------------------------------------------
# Continuous county map used by effect maps and raw predictor maps.
# All maps receive county boundaries + black outer border.
# ------------------------------------------------------------
make_cont_map <- function(df, value_col = "value", title = NULL,
                          legend_title = NULL, limits = NULL,
                          symmetric = TRUE, show_legend = TRUE,
                          sequential = FALSE,
                          uncertainty_aware = FALSE) {
  keep_cols <- unique(c(
    "code_pow",
    value_col,
    intersect(c("ci_excludes_zero", "robust_sign", "se", "ci_low", "ci_high"), names(df))
  ))

  map_vals <- df %>%
    drop_geom() %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    dplyr::distinct(code_pow, .keep_all = TRUE)

  missing_map_codes <- setdiff(pow_geom$code_pow, map_vals$code_pow)
  if (length(missing_map_codes)) {
    warning(
      "Map data missing for ", length(missing_map_codes),
      " analytical counties. They are omitted from the fill layer: ",
      paste(head(missing_map_codes, 12), collapse = ", ")
    )
  }

  sf_d <- pow_geom %>%
    dplyr::inner_join(map_vals, by = "code_pow")

  vals <- sf_d[[value_col]]
  if (is.null(limits)) {
    if (symmetric && !sequential) {
      lim <- robust_abs_limit(vals, FIG_EFFECT_MAP_CLIP_PROB)
      limits <- c(-lim, lim)
    } else {
      limits <- robust_seq_limits(vals, c(.02, .98))
    }
  }
  breaks <- if (symmetric && !sequential) {
    lim_b <- max(abs(limits))
    c(-lim_b, -lim_b / 2, 0, lim_b / 2, lim_b)
  } else {
    three_scale_breaks(limits, symmetric = FALSE)
  }

  use_uncertainty <- isTRUE(uncertainty_aware) &&
    "ci_excludes_zero" %in% names(sf_d)

  if (use_uncertainty) {
    sf_d$map_alpha <- ifelse(
      sf_d$ci_excludes_zero,
      1,
      if (exists("SVC_UNCERTAIN_ALPHA")) SVC_UNCERTAIN_ALPHA else 0.24
    )

    p <- ggplot(sf_d) +
      geom_sf(
        aes(
          fill = .data[[value_col]],
          alpha = map_alpha
        ),
        colour = "#D9D9D9",
        linewidth = FIG_MAP_COUNTY_BORDER_WIDTH
      ) +
      scale_alpha_identity()
  } else {
    p <- ggplot(sf_d) +
      geom_sf(
        aes(fill = .data[[value_col]]),
        colour = "#D9D9D9",
        linewidth = FIG_MAP_COUNTY_BORDER_WIDTH
      )
  }

  if (symmetric && !sequential) {
    p <- p + scale_fill_gradientn(
      colours = DIV_COLS,
      values = scales::rescale(c(-1, -.55, -.20, 0, .20, .55, 1), from = c(-1, 1)),
      limits = limits, breaks = breaks, oob = scales::squish,
      name = legend_title,
      labels = compact_lab,
      guide = guide_colorbar(
        barwidth = grid::unit(2.55, "cm"),
        barheight = grid::unit(.20, "cm"),
        title.position = "top",
        ticks = TRUE
      )
    )
  } else {
    p <- p + scale_fill_gradientn(
      colours = SEQ_COLS,
      limits = limits, breaks = breaks, oob = scales::squish,
      name = legend_title,
      labels = compact_lab,
      guide = guide_colorbar(
        barwidth = grid::unit(2.55, "cm"),
        barheight = grid::unit(.20, "cm"),
        title.position = "top",
        ticks = TRUE
      )
    )
  }

  p +
    geom_sf(data = country_border, fill = NA, colour = "black",
            linewidth = FIG_MAP_COUNTRY_BORDER_WIDTH, inherit.aes = FALSE) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = wrap_fig_text(title, 24)) +
    theme_health_place_map(base_size = 7.2) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 6.5, face = "bold", hjust = .5),
      legend.text = element_text(size = 5.3),
      legend.margin = margin(t = 1, r = 2, b = 1, l = 2),
      plot.title = element_text(size = 6.9, face = "bold", hjust = .5, lineheight = .92,
                                margin = margin(b = 2.5)),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# ------------------------------------------------------------
# Percentile-informed categorical map, shared breaks across panels.
# Used for Figure 1 to improve visual discrimination of high ASMR.
# Numeric class labels preserve the ASMR scale.
# ------------------------------------------------------------
pooled_class_spec <- function(breaks) {
  br <- unique(as.numeric(breaks))
  if (length(br) < 3L) stop("Za mało unikalnych progów dla mapy kwantylowej.")
  class_labels <- vapply(seq_len(length(br) - 1L), function(i) {
    paste0(compact_lab(br[i]), "–", compact_lab(br[i + 1L]))
  }, character(1))
  pal <- grDevices::colorRampPalette(SEQ_COLS)(length(class_labels))
  names(pal) <- class_labels
  list(breaks = br, labels = class_labels, palette = pal)
}

make_pooled_quantile_map <- function(df, breaks, title = NULL,
                                     legend_title = OUTCOME_LABEL,
                                     show_legend = FALSE) {
  stopifnot(length(breaks) >= 3L)
  sp <- pooled_class_spec(breaks)
  br <- sp$breaks

  dd <- df %>% drop_geom() %>% dplyr::select(code_pow, value) %>%
    dplyr::mutate(
      value_clip = pmin(pmax(value, br[1]), br[length(br)]),
      q_class = cut(
        value_clip,
        breaks = br,
        include.lowest = TRUE,
        labels = sp$labels,
        ordered_result = TRUE
      ),
      q_class = factor(q_class, levels = sp$labels, ordered = TRUE)
    )

  # LEFT JOIN is deliberate: counties without a descriptive value remain on the
  # map and are shown in grey rather than disappearing into the white background.
  sf_d <- pow_geom %>% dplyr::left_join(dd, by = "code_pow")

  ggplot(sf_d) +
    geom_sf(aes(fill = q_class), colour = "#D9D9D9", linewidth = FIG_MAP_COUNTY_BORDER_WIDTH) +
    geom_sf(data = country_border, fill = NA, colour = "black",
            linewidth = FIG_MAP_COUNTRY_BORDER_WIDTH, inherit.aes = FALSE) +
    coord_sf(datum = NA, expand = FALSE) +
    scale_fill_manual(
      values = sp$palette,
      limits = sp$labels,
      drop = FALSE,
      na.value = "#BDBDBD",
      na.translate = TRUE,
      name = legend_title,
      guide = guide_legend(
        nrow = 2, byrow = TRUE,
        keywidth = grid::unit(.55, "cm"),
        keyheight = grid::unit(.24, "cm"),
        title.position = "top"
      )
    ) +
    labs(title = title) +
    theme_health_place_map(base_size = 7.7) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 6.7, face = "bold", hjust = .5),
      legend.text = element_text(size = 5.2),
      legend.spacing.x = grid::unit(.05, "cm"),
      plot.title = element_text(size = 8.7, face = "bold", hjust = .5),
      plot.margin = margin(1, 1, 1.5, 1)
    )
}

make_pooled_quantile_legend <- function(breaks, legend_title = OUTCOME_LABEL) {
  sp <- pooled_class_spec(breaks)
  dd <- tibble(
    q_class = factor(sp$labels, levels = sp$labels, ordered = TRUE),
    x = seq_along(sp$labels),
    y = 1
  )
  p_leg <- ggplot(dd, aes(x, y, fill = q_class)) +
    geom_point(shape = 22, size = 0, alpha = 0, show.legend = TRUE) +
    scale_fill_manual(
      values = sp$palette,
      limits = sp$labels,
      drop = FALSE,
      name = legend_title,
      guide = guide_legend(
        nrow = 2, byrow = TRUE,
        override.aes = list(alpha = 1, size = 4, shape = 22),
        keywidth = grid::unit(.55, "cm"),
        keyheight = grid::unit(.24, "cm"),
        title.position = "top"
      )
    ) +
    theme_void(base_family = FIG_FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 6.7, face = "bold", hjust = .5),
      legend.text = element_text(size = 5.2),
      legend.spacing.x = grid::unit(.05, "cm"),
      legend.margin = margin(0, 2, 0, 2),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  gt <- ggplotGrob(p_leg)
  guide_idx <- which(vapply(gt$grobs, function(g) identical(g$name, "guide-box"), logical(1)))
  if (!length(guide_idx)) return(plot_spacer())
  patchwork::wrap_elements(full = gt$grobs[[guide_idx[1]]])
}

# ------------------------------------------------------------
# Shared SSM cartogram scale.
# One colour range is used for BOTH marginal and accumulated-contribution
# maps, and (by default) is locked across M1/M2/M3 within the same Mundlak
# layer. This makes hue directly comparable across paired maps and periods.
# ------------------------------------------------------------
shared_effect_map_limit <- function(component, vars = NULL) {
  if (!component %in% c("between", "within")) stop("Unknown component: ", component)

  if (isTRUE(FIG_EFFECT_MAP_LOCK_ACROSS_PERIODS)) {
    dd <- purrr::map_dfr(R, function(bundle) {
      plot_effect_maps(bundle, component) %>%
        dplyr::select(variable, map_type, value)
    })
  } else {
    stop("shared_effect_map_limit() requires explicit current map data when period locking is disabled.")
  }

  if (!is.null(vars)) dd <- dd %>% dplyr::filter(variable %in% vars)
  robust_abs_limit(dd$value, FIG_EFFECT_MAP_CLIP_PROB)
}

make_diverging_legend_strip <- function(lim, legend_title) {
  lim <- as.numeric(lim)[1]
  if (!is.finite(lim) || lim <= 0) lim <- 1
  dd <- tibble::tibble(x = c(-lim, 0, lim), y = 1, value = c(-lim, 0, lim))
  p_leg <- ggplot(dd, aes(x, y, colour = value)) +
    geom_point(size = 0, alpha = 0, show.legend = TRUE) +
    scale_colour_gradientn(
      colours = DIV_COLS,
      limits = c(-lim, lim),
      breaks = c(-lim, -lim / 2, 0, lim / 2, lim),
      labels = compact_lab,
      oob = scales::squish,
      name = legend_title,
      guide = guide_colorbar(
        barwidth = grid::unit(3.15, "cm"),
        barheight = grid::unit(.22, "cm"),
        title.position = "top",
        ticks = TRUE
      )
    ) +
    theme_void(base_family = FIG_FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 6.4, face = "bold", hjust = .5),
      legend.text = element_text(size = 5.3),
      legend.margin = margin(0, 2, 0, 2),
      legend.box.margin = margin(0, 0, 0, 0)
    )

  gt <- ggplotGrob(p_leg)
  ii <- which(vapply(gt$grobs, function(g) identical(g$name, "guide-box"), logical(1)))
  if (!length(ii)) return(plot_spacer())
  patchwork::wrap_elements(full = gt$grobs[[ii[1]]])
}

# ------------------------------------------------------------
# Grid: curve + marginal map + accumulated contribution map.
# Y-axis is LOCKED across all rows of a given figure.
# Map scales are shared per map column; legends appear once at the bottom.
# ------------------------------------------------------------
make_effect_grid <- function(bundle, component, vars = NULL, filename,
                             stable_only = FALSE, max_n = Inf,
                             title = NULL) {
  g <- bundle$gam
  if (is.null(vars)) vars <- active_vars(bundle, component)
  if (length(vars) == 0L) {
    warning("Brak zmiennych do figury: ", filename)
    return(invisible(NULL))
  }

  curve_df <- g$curves_all %>% dplyr::filter(component == !!component, variable %in% vars)
  map_df <- plot_effect_maps(bundle, component) %>% dplyr::filter(variable %in% vars)

  # SSM cartogram policy: paired marginal/contribution maps use the SAME
  # symmetric colour range. By default the range is also locked across M1/M2/M3
  # within the same Mundlak layer. Robust clipping affects display only.
  if (isTRUE(FIG_EFFECT_MAP_PAIR_SHARED_SCALE)) {
    if (isTRUE(FIG_EFFECT_MAP_LOCK_ACROSS_PERIODS)) {
      lim_pair <- shared_effect_map_limit(component, vars = vars)
    } else {
      lim_pair <- robust_abs_limit(map_df$value, FIG_EFFECT_MAP_CLIP_PROB)
    }
    lim_deriv <- lim_pair
    lim_delta <- lim_pair
  } else {
    lim_deriv <- robust_abs_limit(
      map_df$value[map_df$map_type == "deriv"], FIG_EFFECT_MAP_CLIP_PROB
    )
    lim_delta <- robust_abs_limit(
      map_df$value[map_df$map_type == "delta"], FIG_EFFECT_MAP_CLIP_PROB
    )
  }

  # LOCKED across M1/M2/M3 for the same component.
  y_info <- shared_curve_y_info(component)

  # Keep a fixed, manuscript-friendly determinant order rather than filtering.
  desired <- label_var(c(MODEL_VARS, ONLY_BETWEEN_VARS))
  ord_tbl <- tibble(variable = c(MODEL_VARS, ONLY_BETWEEN_VARS), ord = seq_along(desired))
  if (component == "within") ord_tbl <- ord_tbl %>% dplyr::filter(variable %in% MODEL_VARS)
  vars <- ord_tbl %>% dplyr::filter(variable %in% vars) %>% dplyr::arrange(ord) %>% dplyr::pull(variable)

  rows <- lapply(seq_along(vars), function(i) {
    vv <- vars[i]
    lab <- label_var(vv)
    dc <- curve_df %>% dplyr::filter(variable == vv)
    dm <- map_df %>% dplyr::filter(variable == vv, map_type == "deriv")
    dd <- map_df %>% dplyr::filter(variable == vv, map_type == "delta")

    line_col <- if (component == "between") LANCET$navy else LANCET$red
    term_name <- if (component == "between") paste0(vv, "_between_z") else paste0(vv, "_within_lag", WITHIN_LAG_N, "_z")
    ax <- term_axis_spec(bundle, term_name)

    pc <- ggplot(dc, aes(var_value, effect_median)) +
      geom_hline(yintercept = 0, linewidth = .30, colour = LANCET$grey_mid) +
      geom_ribbon(aes(ymin = effect_p10, ymax = effect_p90), alpha = .07, fill = line_col) +
      geom_ribbon(aes(ymin = effect_q25, ymax = effect_q75), alpha = .16, fill = line_col) +
      geom_line(linewidth = .75, colour = line_col) +
      scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
      scale_y_continuous(breaks = y_info$breaks, labels = compact_lab) +
      coord_cartesian(ylim = c(-y_info$lim, y_info$lim), clip = "on") +
      labs(title = lab, x = ax$title, y = "Δ ASMR / 100,000") +
      theme_lancet(base_size = 7.2) +
      theme(
        plot.title = element_text(size = 7.6, face = "bold"),
        axis.text.x = element_text(size = 5.8),
        axis.text.y = element_text(size = 6.0),
        axis.title.x = element_text(size = 5.9),
        axis.title.y = element_text(size = 6.4),
        panel.grid.major.y = element_line(colour = HAP_COLS$grey_light, linewidth = .18)
      )

    pm <- make_cont_map(
      dm, title = NULL, legend_title = NULL,
      limits = c(-lim_deriv, lim_deriv), symmetric = TRUE,
      show_legend = FALSE,
      uncertainty_aware = isTRUE(SVC_FADE_UNCERTAIN_MAPS)
    )
    pd <- make_cont_map(
      dd, title = NULL, legend_title = NULL,
      limits = c(-lim_delta, lim_delta), symmetric = TRUE,
      show_legend = FALSE,
      uncertainty_aware = isTRUE(SVC_FADE_UNCERTAIN_MAPS)
    )

    (pc | pm | pd) + plot_layout(widths = FIG_EFFECT_GRID_WIDTHS)
  })

  contribution_header <- if (component == "between") {
    "Accumulated contribution\n(observed structural level vs 0 SD)"
  } else {
    "Accumulated contribution\n(+1 SD lag-1 deviation vs 0 SD)"
  }

  header <- wrap_plots(
    ggplot() + annotate("text", .5, .5, label = "Effect curve", fontface = "bold", size = 3) + theme_void(),
    ggplot() + annotate("text", .5, .5, label = "Marginal effect", fontface = "bold", size = 3) + theme_void(),
    ggplot() + annotate("text", .5, .5, label = contribution_header, fontface = "bold", size = 2.7) + theme_void(),
    ncol = 3, widths = FIG_EFFECT_GRID_WIDTHS
  )

  body <- wrap_plots(rows, ncol = 1)

  delta_legend_title <- if (component == "between") {
    "Contribution: observed level vs 0 SD"
  } else {
    "Contribution: +1 SD vs 0 SD"
  }

  legend_row <- wrap_plots(
    plot_spacer(),
    make_diverging_legend_strip(lim_deriv, "Marginal effect per +1 SD"),
    make_diverging_legend_strip(lim_delta, delta_legend_title),
    ncol = 3, widths = FIG_EFFECT_GRID_WIDTHS
  )

  uncertainty_subtitle <- if (
    "ci_excludes_zero" %in% names(map_df) &&
      isTRUE(SVC_FADE_UNCERTAIN_MAPS)
  ) {
    paste0(
      "Maps are uncertainty-aware: pale counties have approximate 95% pointwise intervals that include zero. ",
      "Curve ribbons show cross-county heterogeneity generated by the focal determinant-specific SVC term ",
      "(P10–P90 and IQR), not confidence intervals. Paired cartograms share one colour range, locked across periods ",
      "within the Mundlak layer; display limits use the pooled P", round(100 * FIG_EFFECT_MAP_CLIP_PROB),
      " of |effect| and squish more extreme values to the endpoints."
    )
  } else {
    NULL
  }

  # SSM artwork: no figure-level title/subtitle inside the image.
  # Methodological explanation belongs in the manuscript caption.
  # SSM artwork: no figure-level title/subtitle is added here.
  # IMPORTANT: do not blank plot.title recursively with patchwork `&`,
  # because the determinant names are row-level plot titles on the effect curves.
  p_out <- header / body / legend_row +
    plot_layout(heights = c(.060, 1, .075)) &
    theme(
      plot.margin = margin(t = 3, r = 5, b = 4, l = 4)
    )

  # Full all-term export. Tall by design: no estimated determinant is omitted.
  h <- 0.95 + length(vars) * 1.48
  save_hp_plot(file.path(PUB_DIR, filename), p_out, FIG_W_FULL, h)
  invisible(p_out)
}

# ------------------------------------------------------------
# Dominant accumulated contribution maps.
# For every county and sign, select the determinant with the largest
# absolute accumulated contribution. All estimated terms are used.
# ------------------------------------------------------------
DOM_PALETTE <- c(
  "Higher education" = "#0072B2",
  "Primary / lower-secondary education" = "#E69F00",
  "Long-term unemployment" = "#CC79A7",
  "Gross wages" = "#56B4E9",
  "Physicians" = "#009E73",
  "PM2.5" = "#F0E442",
  "Net migration" = "#D55E00",
  "Population density" = "#6F4C9B",
  "Distance to city with county rights" = "#8C564B",
  "No dominant contribution" = "#EEEEEE"
)

get_dominant_vars <- function(bundle, component, stable_only = FALSE) {
  # ALWAYS use all estimated determinants for dominant-contribution maps.
  plot_effect_maps(bundle, component) %>%
    dplyr::filter(map_type == "delta") %>%
    dplyr::distinct(variable) %>% dplyr::pull(variable)
}

make_dominant_legend <- function(levels_use, palette_use) {
  dd <- tibble(
    category = factor(levels_use, levels = levels_use),
    x = seq_along(levels_use), y = 1
  )
  ggplot(dd, aes(x, y, fill = category)) +
    geom_point(shape = 22, size = 0, alpha = 0, show.legend = TRUE) +
    scale_fill_manual(
      values = palette_use, drop = FALSE, name = "Dominant determinant",
      guide = guide_legend(
        nrow = 3, byrow = TRUE,
        override.aes = list(alpha = 1, size = 4, shape = 22),
        keywidth = grid::unit(.34, "cm"),
        keyheight = grid::unit(.28, "cm")
      )
    ) +
    theme_void(base_family = FIG_FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6.6),
      legend.spacing.x = grid::unit(.10, "cm"),
      legend.box.margin = margin(0, 2, 0, 2)
    )
}

make_dominant_contributor_figure <- function(bundle, filename,
                                             title = NULL,
                                             stable_only = FALSE,
                                             zero_thresh = 0) {
  em <- dplyr::bind_rows(
    plot_effect_maps(bundle, "between"),
    plot_effect_maps(bundle, "within")
  )

  vars_by_comp <- list(
    between = get_dominant_vars(bundle, "between", FALSE),
    within = get_dominant_vars(bundle, "within", FALSE)
  )

  labels_present <- em %>%
    dplyr::filter(map_type == "delta") %>%
    dplyr::transmute(var_label = label_var(variable)) %>%
    dplyr::distinct() %>% dplyr::pull(var_label)

  legend_levels <- c(labels_present, "No dominant contribution")
  legend_levels <- unique(legend_levels)
  unknown <- setdiff(legend_levels, names(DOM_PALETTE))
  if (length(unknown)) {
    extra <- setNames(grDevices::hcl.colors(length(unknown), palette = "Dark 3"), unknown)
    pal_use <- c(DOM_PALETTE, extra)[legend_levels]
  } else {
    pal_use <- DOM_PALETTE[legend_levels]
  }

  one_map <- function(comp, sign, panel_title) {
    vars_use <- vars_by_comp[[comp]]
    d <- em %>%
      dplyr::filter(component == !!comp, map_type == "delta", variable %in% vars_use) %>%
      dplyr::mutate(var_label = label_var(variable))

    all_codes <- pow_geom %>% st_drop_geometry() %>% dplyr::distinct(code_pow)
    if (sign == "positive") {
      signed <- d %>% dplyr::filter(is.finite(value), value > zero_thresh) %>% dplyr::mutate(magnitude = value)
    } else {
      signed <- d %>% dplyr::filter(is.finite(value), value < -zero_thresh) %>% dplyr::mutate(magnitude = abs(value))
    }

    best <- signed %>%
      dplyr::group_by(code_pow) %>%
      dplyr::slice_max(order_by = magnitude, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(code_pow, var_label, magnitude)

    best_full <- all_codes %>%
      dplyr::left_join(best, by = "code_pow") %>%
      dplyr::mutate(
        var_label = ifelse(is.na(var_label), "No dominant contribution", var_label),
        var_label = factor(var_label, levels = legend_levels)
      )

    sf_d <- pow_geom %>% dplyr::left_join(best_full, by = "code_pow")

    ggplot(sf_d) +
      geom_sf(aes(fill = var_label), colour = "#D9D9D9", linewidth = FIG_MAP_COUNTY_BORDER_WIDTH) +
      geom_sf(data = country_border, fill = NA, colour = "black",
              linewidth = FIG_MAP_COUNTRY_BORDER_WIDTH, inherit.aes = FALSE) +
      coord_sf(datum = NA, expand = FALSE) +
      scale_fill_manual(values = pal_use, drop = FALSE, na.translate = FALSE, guide = "none") +
      labs(title = panel_title) +
      theme_health_place_map(base_size = 7.8) +
      theme(
        plot.title = element_text(size = 8.1, face = "bold", hjust = .5),
        plot.margin = margin(2, 2, 2, 2)
      )
  }

  p1 <- one_map("between", "positive", "Between — largest positive")
  p2 <- one_map("within", "positive", "Within — largest positive")
  p3 <- one_map("between", "negative", "Between — largest negative")
  p4 <- one_map("within", "negative", "Within — largest negative")
  leg <- make_dominant_legend(legend_levels, pal_use)

  # SSM artwork: retain only short panel headings; full explanation goes in caption.
  p_out <- ((p1 | p2) / (p3 | p4)) / leg +
    plot_layout(heights = c(1, 1, .28)) &
    theme(
      plot.title = element_text(size = 8.1, face = "bold", hjust = .5),
      plot.subtitle = element_blank(),
      plot.margin = margin(t = 3, r = 5, b = 4, l = 4)
    )

  save_hp_plot(file.path(PUB_DIR, filename), p_out, FIG_W_FULL, FIG_W_FULL * .92)
  invisible(p_out)
}

# ------------------------------------------------------------
# Figure 1: descriptive outcome geography, 3 non-overlapping eras.
# Shared pooled quantile breakpoints make high-ASMR counties visible and
# keep the classification identical across periods.
# ------------------------------------------------------------
period_eras <- list(
  `2010–2019` = 2010:2019,
  `2020–2022` = 2020:2022,
  `2023–2025` = 2023:2025
)

desc_outcome <- map_dfr(names(period_eras), function(lbl) {
  yy <- period_eras[[lbl]]
  dat %>% drop_geom() %>% dplyr::filter(year %in% yy) %>%
    dplyr::group_by(code_pow) %>%
    dplyr::summarise(value = mean(.data[[OUTCOME]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(era = lbl)
})

# Tail-sensitive pooled percentiles: more visual resolution in the tails,
# especially the high-mortality end. Same numeric cut points in all panels.
fig1_probs <- c(0, .10, .25, .50, .75, .90, 1)
fig1_breaks <- as.numeric(quantile(desc_outcome$value, probs = fig1_probs,
                                   na.rm = TRUE, names = FALSE))
fig1_breaks <- unique(fig1_breaks)
if (length(fig1_breaks) < 4L) {
  fig1_breaks <- pretty(desc_outcome$value, n = 6)
}

# Explicit missingness audit: missing counties are grey, never white/blank.
fig1_missingness <- tidyr::expand_grid(
  era = names(period_eras),
  code_pow = pow_geom$code_pow
) %>%
  dplyr::left_join(desc_outcome, by = c("era", "code_pow")) %>%
  dplyr::group_by(era) %>%
  dplyr::summarise(
    n_counties = n(),
    n_missing = sum(!is.finite(value)),
    pct_missing = 100 * n_missing / n_counties,
    .groups = "drop"
  )
readr::write_csv(
  fig1_missingness,
  file.path(PUB_DIR, "Figure_1_descriptive_missingness.csv")
)

plots_desc <- lapply(names(period_eras), function(lbl) {
  make_pooled_quantile_map(
    desc_outcome %>% dplyr::filter(era == lbl),
    breaks = fig1_breaks,
    title = lbl,
    legend_title = paste0(OUTCOME_LABEL, "\npooled percentile-informed classes"),
    show_legend = FALSE
  )
})

# One deliberately constructed legend for all three maps. This avoids patchwork
# duplicating guides when one panel does not contain every observed class.
fig1_leg <- make_pooled_quantile_legend(
  fig1_breaks,
  paste0(OUTCOME_LABEL, "\npooled percentile-informed classes")
)

fig1_body <- wrap_plots(plots_desc, nrow = 1)
# SSM artwork: panel-year labels stay; the full figure title belongs in the caption.
# Extra legend height/bottom margin prevents clipping in PNG/PDF export.
fig1 <- (fig1_body / fig1_leg) +
  plot_layout(heights = c(1, .24)) &
  theme(
    plot.title = element_text(size = 8.7, face = "bold", hjust = .5),
    plot.margin = margin(t = 3, r = 5, b = 7, l = 4)
  )

save_hp_plot(file.path(PUB_DIR, "Figure_1_descriptive_outcome.png"),
             fig1, FIG_W_FULL, FIG_W_FULL * .52)

# Export the exact colour limits used by SVC effect cartograms for QA/captions.
scale_meta <- dplyr::bind_rows(lapply(c("between", "within"), function(cc) {
  lim <- shared_effect_map_limit(cc)
  tibble::tibble(
    component = cc,
    lower = -lim,
    upper = lim,
    pooled_abs_quantile = FIG_EFFECT_MAP_CLIP_PROB,
    locked_across_periods = FIG_EFFECT_MAP_LOCK_ACROSS_PERIODS,
    paired_map_columns_share_scale = FIG_EFFECT_MAP_PAIR_SHARED_SCALE
  )
}))
readr::write_csv(scale_meta, file.path(PUB_DIR, "Figure_effect_cartogram_scale_metadata.csv"))

# ------------------------------------------------------------
# Main Model 1 figures
# ------------------------------------------------------------
M1 <- R[["m1_2010_2019"]]
make_effect_grid(
  M1, "between",
  filename = "Figure_2_M1_between.png",
  stable_only = FALSE, max_n = Inf,
  title = "Between-county effects | Model 1 (2010–2019)"
)
make_effect_grid(
  M1, "within",
  filename = "Figure_3_M1_within.png",
  stable_only = FALSE, max_n = Inf,
  title = "Within-county lag-1 effects | Model 1 (2010–2019)"
)

# ------------------------------------------------------------
# Figure 4: linear Mundlak vs GAM, Model 1
# ------------------------------------------------------------
if (!is.null(M1$ridge)) {
  rdf <- M1$ridge$coef_final %>%
    dplyr::transmute(variable, var_label, component,
              estimator = "Linear Mundlak β",
              effect = beta_debiased, lo = conf.low, hi = conf.high,
              interval = "95% bootstrap CI")
  if (!is.null(M1$gam$aps_inference) && nrow(M1$gam$aps_inference)) {
    gdf <- M1$gam$aps_inference %>%
      dplyr::transmute(
        variable, var_label, component,
        estimator = "Spatial Mundlak GAM",
        effect = average_derivative_mean,
        lo = ci_low,
        hi = ci_high,
        interval = "approx. 95% pointwise CI"
      )
  } else {
    stop(
      "Figure 4 requires M1 GAM coefficient-covariance intervals (aps_inference). ",
      "Run 02b_svc_uncertainty.R for M1; no GAM refit is needed."
    )
  }
  fd <- dplyr::bind_rows(rdf, gdf)
  order_v <- M1$gam$aps_all %>% dplyr::filter(component == "between") %>%
    dplyr::arrange(aps) %>% dplyr::pull(var_label)
  order_all <- c(order_v, setdiff(unique(fd$var_label), order_v))
  fd$var_label <- factor(fd$var_label, levels = order_all)
  est_cols <- c("Linear Mundlak β" = LANCET$navy,
                "Spatial Mundlak GAM" = LANCET$red)

  p4 <- ggplot(fd, aes(effect, var_label, colour = estimator, shape = estimator)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = LANCET$grey_mid, linewidth = .3) +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(.55), linewidth = .75) +
    geom_point(position = position_dodge(.55), size = 2) +
    facet_wrap(~component, scales = "free_x") +
    scale_colour_manual(values = est_cols, name = NULL) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    scale_y_discrete(drop = FALSE) +
    labs(
      x = "Δ ASMR / 100,000 per +1 SD", y = NULL
    ) +
    theme_lancet(base_size = 8.2) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(t = 4, r = 6, b = 4, l = 4)
    )

  save_hp_plot(file.path(PUB_DIR, "Figure_4_M1_linear_vs_GAM.png"),
               p4, FIG_W_FULL, max(5.2, length(order_all) * .36 + 2.2))
}

# ------------------------------------------------------------
# Figure 5: dominant accumulated contribution, Model 1
# ------------------------------------------------------------
make_dominant_contributor_figure(
  M1,
  filename = "Figure_5_M1_dominant_contribution.png",
  title = "Dominant accumulated determinant contribution | Model 1 (2010–2019)",
  stable_only = FALSE,
  zero_thresh = 0
)

# ------------------------------------------------------------
# Supplement: all three periods APS comparison
# ------------------------------------------------------------
aps_periods <- map_dfr(PERIOD_ORDER, function(pp) {
  gg <- R[[pp]]$gam

  # IMPORTANT: build the branch result first, then mutate it.
  # Piping directly from `if (...) { ... } else { ... } %>% dplyr::mutate(...)`
  # attaches the pipe only to the `else` branch in R. Once aps_inference
  # exists, that silently drops period_label/period_id and later factor()
  # assignment fails with a 51-row vs 0-row size mismatch.
  out <- if (!is.null(gg$aps_inference) && nrow(gg$aps_inference)) {
    gg$aps_inference %>%
      dplyr::transmute(
        variable, var_label, component,
        aps = average_derivative_mean,
        aps_q25 = ci_low,
        aps_q75 = ci_high
      )
  } else {
    gg$aps_all %>%
      dplyr::select(variable, var_label, component, aps, aps_q25, aps_q75)
  }

  out %>%
    dplyr::mutate(
      period_label = PERIOD_SPECS[[pp]]$year_label,
      period_id = pp
    )
})

required_aps_cols <- c(
  "variable", "var_label", "component", "aps", "aps_q25", "aps_q75",
  "period_label", "period_id"
)
missing_aps_cols <- setdiff(required_aps_cols, names(aps_periods))
if (!nrow(aps_periods) || length(missing_aps_cols)) {
  stop(
    "Nie udało się zbudować tabeli APS across periods. Brakujące kolumny: ",
    paste(missing_aps_cols, collapse = ", "),
    ". Sprawdź gam_results.rds / aps_inference dla M1-M3."
  )
}

aps_periods$period_label <- factor(
  aps_periods$period_label,
  levels = vapply(PERIOD_ORDER, function(pp) PERIOD_SPECS[[pp]]$year_label, character(1))
)
var_order <- aps_periods %>%
  dplyr::filter(component == "between", period_id == "m1_2010_2019") %>%
  dplyr::arrange(aps) %>% dplyr::pull(var_label)
aps_periods$var_label <- factor(
  aps_periods$var_label,
  levels = c(var_order, setdiff(unique(aps_periods$var_label), var_order))
)
period_cols <- c(
  "2010–2019" = HAP_COLS$blue,
  "2010–2022" = HAP_COLS$orange,
  "2010–2025" = HAP_COLS$green
)

pS1 <- ggplot(aps_periods,
              aes(aps, var_label, colour = period_label, shape = period_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = LANCET$grey_mid, linewidth = .3) +
  geom_linerange(aes(xmin = aps_q25, xmax = aps_q75),
                 position = position_dodge(.65), linewidth = .7) +
  geom_point(position = position_dodge(.65), size = 1.9) +
  facet_wrap(~component, scales = "free_x") +
  scale_colour_manual(values = period_cols, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(
    x = "Δ ASMR / 100,000 per +1 SD", y = NULL
  ) +
  theme_lancet(base_size = 8.2) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(t = 4, r = 6, b = 4, l = 4)
  )

save_hp_plot(file.path(PUB_DIR, "Figure_S1_APS_period_comparison.png"),
             pS1, FIG_W_FULL, max(5.2, length(levels(aps_periods$var_label)) * .36 + 2.4))

# Full effect grids for all periods/components. The same locked-y / single-legend
# rules apply automatically because they use make_effect_grid().
supp_idx <- 2L
for (pp in PERIOD_ORDER) {
  for (cc in c("between", "within")) {
    nm <- sprintf("Figure_S%d_%s_%s_all.png", supp_idx, pp, cc)
    make_effect_grid(
      R[[pp]], cc, filename = nm,
      stable_only = FALSE, max_n = Inf,
      title = paste0(
        ifelse(cc == "between", "Between-county", "Within-county lag-1"),
        " effects | ", PERIOD_SPECS[[pp]]$year_label
      )
    )
    supp_idx <- supp_idx + 1L
  }
}

# ------------------------------------------------------------
# Raw predictor maps: each variable gets a robust P2-P98 numeric scale.
# Every map retains black country outline and county boundaries.
# ------------------------------------------------------------
raw_means <- dat %>% drop_geom() %>%
  dplyr::filter(year >= 2010, year <= 2025) %>%
  dplyr::group_by(code_pow) %>%
  dplyr::summarise(across(dplyr::all_of(ALL_RAW_VARS), ~ mean(to_numeric_safe(.x), na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(dplyr::all_of(ALL_RAW_VARS), names_to = "variable", values_to = "value")

plots_raw <- lapply(ALL_RAW_VARS, function(vv) {
  dd <- raw_means %>% dplyr::filter(variable == vv)
  make_cont_map(
    dd,
    title = raw_panel_title(vv, unit = unit_var(vv), width = 20),
    legend_title = unit_var(vv),
    symmetric = FALSE, sequential = TRUE,
    show_legend = TRUE
  )
})
raw_grid <- wrap_plots(plots_raw, ncol = 3) &
  theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 4))

save_hp_plot(file.path(PUB_DIR, sprintf("Figure_S%d_raw_predictor_means.png", supp_idx)),
             raw_grid, FIG_W_FULL, ceiling(length(plots_raw) / 3) * 2.58)
supp_idx <- supp_idx + 1L

# Annual raw changes only for time-varying predictors.
raw_changes <- dat %>% drop_geom() %>%
  dplyr::filter(year >= 2010, year <= 2025) %>%
  dplyr::arrange(code, year) %>%
  dplyr::group_by(code, code_pow) %>%
  dplyr::summarise(across(dplyr::all_of(MODEL_VARS), ~ {
    xx <- to_numeric_safe(.x)
    yy <- year
    ok <- is.finite(xx) & is.finite(yy)
    if (sum(ok) < 2L) NA_real_ else unname(coef(lm(xx[ok] ~ yy[ok]))[2])
  }), .groups = "drop") %>%
  dplyr::group_by(code_pow) %>%
  dplyr::summarise(across(dplyr::all_of(MODEL_VARS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(dplyr::all_of(MODEL_VARS), names_to = "variable", values_to = "value")

plots_change <- lapply(MODEL_VARS, function(vv) {
  dd <- raw_changes %>% dplyr::filter(variable == vv)
  make_cont_map(
    dd,
    title = raw_panel_title(vv, suffix = "annual change", width = 20),
    legend_title = paste0(unit_var(vv), "/year"),
    symmetric = TRUE, sequential = FALSE,
    show_legend = TRUE
  )
})
change_grid <- wrap_plots(plots_change, ncol = 3) &
  theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 4))

save_hp_plot(file.path(PUB_DIR, sprintf("Figure_S%d_raw_predictor_changes.png", supp_idx)),
             change_grid, FIG_W_FULL, ceiling(length(plots_change) / 3) * 2.45)
supp_idx <- supp_idx + 1L

# Dominant-contribution sensitivity maps for Models 2 and 3.
for (pp in c("m2_2010_2022", "m3_2010_2025")) {
  nm <- sprintf("Figure_S%d_%s_dominant_contribution.png", supp_idx, pp)
  make_dominant_contributor_figure(
    R[[pp]],
    filename = nm,
    title = paste0("Dominant accumulated determinant contribution | ", PERIOD_SPECS[[pp]]$year_label),
    stable_only = FALSE,
    zero_thresh = 0
  )
  supp_idx <- supp_idx + 1L
}

# ------------------------------------------------------------
# Manifest — single source of figure numbering.
# ------------------------------------------------------------
manifest <- tibble(
  figure = c(
    "Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5",
    paste0("Figure S", 1:11)
  ),
  file = c(
    "Figure_1_descriptive_outcome.png",
    "Figure_2_M1_between.png",
    "Figure_3_M1_within.png",
    "Figure_4_M1_linear_vs_GAM.png",
    "Figure_5_M1_dominant_contribution.png",
    "Figure_S1_APS_period_comparison.png",
    "Figure_S2_m1_2010_2019_between_all.png",
    "Figure_S3_m1_2010_2019_within_all.png",
    "Figure_S4_m2_2010_2022_between_all.png",
    "Figure_S5_m2_2010_2022_within_all.png",
    "Figure_S6_m3_2010_2025_between_all.png",
    "Figure_S7_m3_2010_2025_within_all.png",
    "Figure_S8_raw_predictor_means.png",
    "Figure_S9_raw_predictor_changes.png",
    "Figure_S10_m2_2010_2022_dominant_contribution.png",
    "Figure_S11_m3_2010_2025_dominant_contribution.png"
  ),
  role = c(
    "Descriptive mortality geography: 2010–2019, 2020–2022, 2023–2025; pooled percentile-informed ASMR classes",
    "Model 1 between-county effects; all 9 estimated between determinants",
    "Model 1 within-county lag-1 effects; all 8 estimated within determinants",
    "Model 1 linear Mundlak benchmark vs spatial GAM",
    "Model 1 dominant contribution by county/sign using all determinants; within uses +1 SD vs 0 SD",
    "GAM partial-slope sensitivity across all three overlapping periods",
    "Full Model 1 between-county effect grid",
    "Full Model 1 within-county effect grid",
    "Full Model 2 between-county effect grid",
    "Full Model 2 within-county effect grid",
    "Full Model 3 between-county effect grid",
    "Full Model 3 within-county effect grid",
    "Mean raw determinant values, 2010–2025",
    "County-specific linear annual changes in time-varying raw determinants, 2010–2025",
    "Model 2 dominant contribution by county/sign using all determinants; within uses +1 SD vs 0 SD",
    "Model 3 dominant contribution by county/sign using all determinants; within uses +1 SD vs 0 SD"
  )
)
readr::write_csv(manifest, file.path(PUB_DIR, "figure_manifest.csv"))

cat("\nPublication figures DONE:", PUB_DIR, "\n")
cat("Figure patch: ALL effects plotted; raw-unit x labels at 1-SD spacing; locked y-scales;\n")
cat("sparse legends; black country outline; within contribution = +1 SD vs 0; all-term dominant maps.\n")
