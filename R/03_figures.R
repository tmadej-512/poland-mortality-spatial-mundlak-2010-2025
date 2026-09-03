# ============================================================
# 03_figures.R
# Wykresy publikacyjne — Main + Supplement
#
# Wymaga: powiaty_spdf w pamięci
# Wczytuje: output/pre_covid/ + output/full/
#
# Produkuje:
#   output/figures/
#     Figure_1.png / Figure_2.png / Figure_3.png / Figure_4.png
#     fig4_covid_comparison.png / fig5_predicted_total_powiat_pre_full.png
#     suppl_fig1_mean_asmr_periods.png
#     suppl_fig2_desc_means_pre.png / suppl_fig3_desc_changes_pre.png
#     suppl_fig4_desc_means_full.png / suppl_fig5_desc_changes_full.png
#     suppl_fig6_observed_vs_predicted_total_powiat.png
#     suppl_fig7_prediction_error_total_powiat.png
#     suppl_fig8_predicted_mean_powiat_pre_full.png
#     suppl_fig9_between_full.png / suppl_fig10_within_full.png
#     suppl_fig11_gi_between.png / suppl_fig12_gi_within.png
#     suppl_fig13_between_all.png / suppl_fig14_within_all.png
#     suppl_fig15_gi_between_full.png / suppl_fig16_gi_within_full.png
# ============================================================

if (!exists("LANCET")) source(file.path(getwd(), "R", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(ggplot2); library(patchwork)
  library(forcats); library(tibble); library(sf)
  library(scales)
})

FIG_DIR <- "output/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. WCZYTAJ WYNIKI OBU OKRESÓW
# ============================================================

load_results <- function(period) {
  gam_path     <- file.path("output", period,
                             "stvc_gam", "gam_results.rds")
  mundlak_path <- file.path("output", period,
                             "mundlak_ridge",
                             "mundlak_results.rds")

  if (!file.exists(gam_path))
    stop("Brak: ", gam_path,
         "\nUruchom najpierw 02_stvc_gam.R z PERIOD='",
         period, "'")
  if (!file.exists(mundlak_path))
    stop("Brak: ", mundlak_path,
         "\nUruchom najpierw 01_mundlak_ridge.R z PERIOD='",
         period, "'")

  gam     <- readRDS(gam_path)
  mundlak <- readRDS(mundlak_path)

  pred_path <- file.path("output", period,
                         "stvc_gam", "tables",
                         "powiat_predictions_summary.csv")
  pred_powiat <- if (!is.null(gam$pred_powiat)) {
    gam$pred_powiat
  } else if (file.exists(pred_path)) {
    read.csv(pred_path, stringsAsFactors = FALSE)
  } else {
    tibble()
  }

  snr_path <- file.path("output", period,
                         "perturbation", "tables",
                         "snr_table.csv")
  snr_table <- if (file.exists(snr_path)) {
    read.csv(snr_path, stringsAsFactors = FALSE)
  } else {
    warning("Brak snr_table dla ", period,
            " — active_only filtr nie będzie działał")
    NULL
  }

  list(
    period        = period,
    year_label    = gam$settings$YEAR_LABEL,
    m             = gam$m,
    dat_all       = gam$dat_all,
    between_terms = gam$between_terms,
    within_terms  = gam$within_terms,
    all_terms     = gam$all_terms,
    edf_tbl       = gam$edf_tbl,
    aps_all       = gam$aps_all %||% {
      warning("gam_results.rds nie zawiera aps_all; używam legacy ame_all jako fallback.")
      gam$ame_all %>%
        dplyr::rename(
          aps = ame,
          aps_mean = ame,
          aps_q25 = ame_q25,
          aps_q75 = ame_q75,
          aps_p10 = ame_p10,
          aps_p90 = ame_p90
        )
    },
    aps_county_all = gam$aps_county_all %||% tibble(),
    curves_all    = gam$curves_all,
    effect_maps   = gam$effect_maps,
    gi_all        = gam$gi_all,
    pred_powiat   = pred_powiat,
    coef_final    = mundlak$coef_final,
    snr_table     = snr_table
  )
}

cat("Wczytywanie wyników...\n")
res_pre  <- load_results("pre_covid")
res_full <- load_results("full")

# ============================================================
# Figure 4 simulation removed by design.
# ============================================================


VAR_LABELS_USE <- res_pre$aps_all %>%
  dplyr::distinct(variable, var_label) %>%
  deframe()

# ============================================================
# 2. GEOMETRIA
# ============================================================

powiaty_geom <- powiaty_spdf %>%
  dplyr::mutate(
    JPT_KOD_JE = as.character(JPT_KOD_JE),
    JPT_KOD_JE = ifelse(nchar(JPT_KOD_JE) == 3,
                        paste0("0", JPT_KOD_JE),
                        JPT_KOD_JE)
  ) %>%
  dplyr::rename(code_pow = JPT_KOD_JE) %>%
  dplyr::select(code_pow, geometry)

# Skrócone etykiety osi i legend: 1k, 1.2M itd.
format_compact_num <- function(x, accuracy_small = 0.1) {
  fmt_one <- function(z) {
    if (is.na(z)) return(NA_character_)
    if (!is.finite(z)) return(as.character(z))
    az <- abs(z)
    out <- if (az >= 1e9) {
      sprintf("%.1fB", z / 1e9)
    } else if (az >= 1e6) {
      sprintf("%.1fM", z / 1e6)
    } else if (az >= 1e3) {
      sprintf("%.1fk", z / 1e3)
    } else {
      digits <- if (accuracy_small < 1) nchar(strsplit(as.character(accuracy_small), ".", fixed = TRUE)[[1]][2]) else 0
      formatC(round(z / accuracy_small) * accuracy_small, format = "f", digits = digits)
    }
    out <- sub("\\.0([kMB])$", "\\1", out)
    out <- sub("\\.0$", "", out)
    out
  }
  vapply(x, fmt_one, character(1))
}

label_compact_num <- function(accuracy_small = 0.1) {
  force(accuracy_small)
  function(x) format_compact_num(x, accuracy_small = accuracy_small)
}

make_map <- function(df_map, fill_col,
                      title = NULL, subtitle = NULL,
                      legend_name = "value",
                      palette = "RdBu",
                      symmetric = TRUE,
                      limits = NULL,
                      show_border = TRUE,
                      show_legend = FALSE) {
  sf_data <- df_map %>%
    dplyr::left_join(powiaty_geom, by = "code_pow") %>%
    st_as_sf()

  country_border <- if (show_border) {
    tryCatch(
      sf::st_union(sf_data) %>%
        sf::st_simplify(dTolerance = 500),
      error = function(e) NULL
    )
  } else NULL
  vals <- sf_data[[fill_col]]
  lims <- if (!is.null(limits)) limits else {
    if (symmetric) {
      lim <- max(abs(vals), na.rm = TRUE)
      if (!is.finite(lim) || lim == 0) lim <- 1
      c(-lim, lim)
    } else {
      r <- range(vals, na.rm = TRUE)
      if (r[1] == r[2]) r + c(-1, 1) else r
    }
  }
  p <- ggplot(sf_data) +
    geom_sf(aes(fill = .data[[fill_col]]), color = NA)

  if (symmetric || palette == "RdBu") {
    p <- p + scale_fill_gradient2(
      low = HAP_COLS$blue,
      mid = HAP_COLS$white,
      high = HAP_COLS$vermillion,
      midpoint = 0,
      limits = lims,
      name = legend_name,
      labels = label_compact_num(accuracy_small = 0.1),
      oob = squish,
      guide = guide_colorbar(
        barwidth = unit(4.8, "cm"),
        barheight = unit(0.28, "cm")
      )
    )
  } else {
    p <- p + scale_fill_gradient(
      low = HAP_COLS$grey_panel,
      high = HAP_COLS$blue,
      limits = lims,
      name = legend_name,
      labels = label_compact_num(accuracy_small = 0.1),
      oob = squish,
      guide = guide_colorbar(
        barwidth = unit(4.8, "cm"),
        barheight = unit(0.28, "cm")
      )
    )
  }
  if (!is.null(country_border)) {
    p <- p + geom_sf(data = country_border,
                     fill = NA, colour = "grey25",
                     linewidth = 0.35,
                     inherit.aes = FALSE)
  }
  p <- p + labs(title = title, subtitle = subtitle) +
    theme_health_place_map() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(hjust = 0, size = 6, face = "bold"),
      legend.text = element_text(size = 5.5),
      legend.box.margin = margin(t = 2),
      plot.caption = element_blank()
    )
  if (!show_legend) p <- p + guides(fill = "none")
  p
}

# ============================================================
# 3. FUNKCJA: GRID KRZYWA + MAPA POCHODNA + MAPA DELTA
# ============================================================

# ============================================================
# 3. FUNKCJA: GRID KRZYWA + ACCUMULATED MAP
#    Between: burden contrast vs P2
#    Within:  average temporal-deviation contribution vs 0
#    Układ:   2 kolumny × n_vars wierszy
#    Submission: plot_annotation usuwany (caption w manuskrypcie)
# ============================================================


# ============================================================
# HELPER: Mapa kwantylowa (globalna — używana przez wiele funkcji)
# ============================================================
# make_sym_breaks: symetryczne biny dla map z wartosciami ujemnymi i pozytywnymi
# n_q musi byc parzyste (np. 10 = 5 neg + 5 pos)
make_sym_breaks <- function(vals, n_q = 10L) {
  lim <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  n_half <- n_q %/% 2L
  neg_breaks <- seq(-lim, 0, length.out = n_half + 1)
  pos_breaks <- seq(0,  lim, length.out = n_half + 1)
  unique(c(neg_breaks, pos_breaks))
}

make_map_quantile <- function(df_m, fill_col,
                               title        = NULL,
                               subtitle     = NULL,
                               legend_name  = "\u0394 ASMR",
                               n_q          = 10L,
                               breaks       = NULL,
                               pal_cols     = NULL,
                               pal          = NULL,
                               show_legend  = TRUE) {
  if (nrow(df_m) == 0) return(ggplot() + theme_void())

  sf_data <- df_m %>%
    dplyr::left_join(powiaty_geom, by = "code_pow") %>%
    sf::st_as_sf()

  vals <- sf_data[[fill_col]]

  # Breaks: jesli nie podane z zewnatrz, generuj symetrycznie
  if (is.null(breaks)) {
    has_neg <- any(vals < 0, na.rm = TRUE)
    if (has_neg) {
      breaks <- make_sym_breaks(vals, n_q)
    } else {
      breaks <- unique(quantile(vals,
                                 probs = seq(0, 1, length.out = n_q + 1),
                                 na.rm = TRUE))
      if (length(breaks) < 3) breaks <- pretty(vals, n = n_q)
    }
  }

  n_cls <- length(breaks) - 1
  has_neg <- any(breaks < 0)

  # Paleta: jesli nie podana z zewnatrz. Argument `pal` jest
  # zachowany dla zgodnosci z wywolaniami map opisowych.
  if (is.null(pal_cols)) {
    pal_name <- if (!is.null(pal)) pal else if (has_neg) "RdBu" else "Reds"
    base_n <- if (pal_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      min(RColorBrewer::brewer.pal.info[pal_name, "maxcolors"], 11L)
    } else {
      9L
    }
    base_cols <- RColorBrewer::brewer.pal(base_n, pal_name)
    if (identical(pal_name, "RdBu")) base_cols <- rev(base_cols)
    pal_cols <- colorRampPalette(base_cols)(n_cls)
  }

  lbl <- sapply(seq_len(n_cls), function(k)
    sprintf("%.0f\u2013%.0f", breaks[k], breaks[k+1]))

  sf_data$q_class  <- cut(vals, breaks = breaks,
                           include.lowest = TRUE, labels = FALSE)
  sf_data$q_factor <- factor(sf_data$q_class,
                              levels = seq_len(n_cls),
                              labels = lbl)

  country_border <- tryCatch(
    sf::st_union(sf_data) %>% sf::st_simplify(dTolerance = 500),
    error = function(e) NULL
  )

  p <- ggplot(sf_data) +
    geom_sf(aes(fill = q_factor), colour = NA) +
    scale_fill_manual(values   = pal_cols,
                      name     = legend_name,
                      na.value = "grey85",
                      drop     = FALSE) +
    { if (!is.null(country_border))
        geom_sf(data = country_border, fill = NA,
                colour = "grey25", linewidth = 0.35,
                inherit.aes = FALSE)
      else NULL } +
    labs(title = title, subtitle = subtitle) +
    guides(fill = if (show_legend) {
      guide_legend(
        nrow      = 1L,
        reverse   = FALSE,
        keywidth  = unit(0.5, "cm"),
        keyheight = unit(0.25, "cm"),
        label.position = "bottom"
      )
    } else {
      "none"
    }) +
    theme_health_place_map() +
    theme(
      legend.position   = "bottom",
      legend.text       = element_text(size = 5.5),
      legend.title      = element_text(size = 6, face = "bold", hjust = 0),
      legend.box.margin = margin(t = 2),
      plot.caption      = element_blank()
    )
  p
}

make_grid_figure <- function(res, comp_name,
                              filename,
                              zero_thresh    = 1.0,
                              active_only    = TRUE,
                              add_annotation = FALSE,
                              n_quantile     = 5L) {

  terms_comp <- if (comp_name == "between")
    res$between_terms else res$within_terms
  col <- if (comp_name == "between")
    LANCET$navy else LANCET$red

  # Filtr obowiązkowy: zmienne o wyzerowanym efekcie nie są pokazywane
  # w żadnej figurze (main ani supplement).
  if (!is.null(res$curves_all)) {
    active_vars <- res$curves_all %>%
      dplyr::filter(component == comp_name) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        max_abs = max(abs(effect_median), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(max_abs > zero_thresh) %>%
      dplyr::pull(variable)
    terms_comp <- terms_comp[clean_var(terms_comp) %in% active_vars]
  }

  n_vars <- length(terms_comp)
  if (n_vars == 0) {
    message("Brak termów: ", comp_name); return(NULL)
  }

  # Wspolna skala Y krzywych
  ylim_df <- res$curves_all %>%
    dplyr::filter(component == comp_name,
                  variable %in% clean_var(terms_comp))
  ylim_range <- range(c(ylim_df$effect_p10, ylim_df$effect_p90),
                       na.rm = TRUE)
  pad      <- diff(ylim_range) * 0.05
  ylim_use <- ylim_range + c(-pad, pad)

  # Stara ciągła paleta: jedna wspólna, symetryczna skala dla wszystkich
  # kartogramów (marginal + contribution) w obrębie całej figury.
  all_maps <- res$effect_maps %>%
    dplyr::filter(component == comp_name,
                  map_type %in% c("deriv", "delta"),
                  variable %in% clean_var(terms_comp))

  # Wspólny, odporny na wartości odstające zakres kolorów.
  # Funkcja make_grid_figure() jest wywoływana osobno dla "between"
  # i "within", dlatego 98. percentyl jest liczony niezależnie dla obu figur.
  # Wartości poza zakresem nie znikają: scale_fill_gradient2() używa
  # oob = scales::squish i nadaje im skrajny kolor palety.
  finite_map_values <- abs(all_maps$value[is.finite(all_maps$value)])

  map_lim <- if (length(finite_map_values) > 0L) {
    as.numeric(stats::quantile(
      finite_map_values,
      probs = 0.98,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    ))
  } else {
    1
  }

  if (!is.finite(map_lim) || map_lim <= 0) map_lim <- 1
  limits_maps <- c(-map_lim, map_lim)

  # Buduj panele per zmienna
  panels <- list()

  for (i in seq_along(terms_comp)) {
    var    <- terms_comp[i]
    vclean <- clean_var(var)
    vlabel <- label_var(vclean)

    df_c <- res$curves_all %>%
      dplyr::filter(variable == vclean, component == comp_name)

    is_zero <- nrow(df_c) == 0 ||
      max(abs(df_c$effect_median), na.rm = TRUE) < zero_thresh

    bg <- if (is_zero) "#F5F5F5" else LANCET$white
    lc <- if (is_zero) LANCET$grey_mid else col

    # 1. Effect curve
    p_c <- if (nrow(df_c) > 0) {
      ggplot(df_c, aes(x = var_value)) +
        annotate("rect", xmin=-Inf, xmax=Inf,
                 ymin=-Inf, ymax=Inf, fill=bg) +
        geom_hline(yintercept = 0,
                   colour = LANCET$grey_mid, linewidth = 0.3) +
        geom_ribbon(aes(ymin=effect_p10, ymax=effect_p90),
                    fill=lc, alpha=0.10, colour=NA) +
        geom_ribbon(aes(ymin=effect_q25, ymax=effect_q75),
                    fill=lc, alpha=0.22, colour=NA) +
        geom_line(aes(y=effect_median),
                  colour=lc, linewidth=0.9) +
        coord_cartesian(ylim = ylim_use) +
        labs(title = vlabel,
             x = "Determinant value (SD)",
             y = "Δ ASMR / 100k") +
        theme_lancet(base_size = 8) +
        theme(panel.background = element_rect(fill=bg, colour=NA))
    } else ggplot() + theme_void()

    # 2. Marginal map
    df_mg <- res$effect_maps %>%
      dplyr::filter(variable == vclean,
                    component == comp_name,
                    map_type  == "deriv")

    p_mg <- make_map(
      df_mg, "value",
      title        = NULL,
      legend_name  = "Effect value",
      palette      = "RdBu",
      symmetric    = TRUE,
      limits       = limits_maps,
      show_legend  = (i == 1L)
    )

    # 3. Contribution map (accumulated)
    df_ac <- res$effect_maps %>%
      dplyr::filter(variable == vclean,
                    component == comp_name,
                    map_type  == "delta")

    p_ac <- make_map(
      df_ac, "value",
      title        = NULL,
      legend_name  = "Effect value",
      palette      = "RdBu",
      symmetric    = TRUE,
      limits       = limits_maps,
      show_legend  = FALSE
    )

    panels[[i]] <- list(curve = p_c, marginal = p_mg, contrib = p_ac,
                        label = vlabel)
  }

  # ============================================================
  # Uklad 2xN: between → pary zmiennych obok siebie (2 kolumny
  # zmiennych, kazda zmienna = 3 panele)
  # Within → 1 kolumna zmiennych, 3 panele per zmienna
  # ============================================================

  hdr <- function(txt, size=3.0)
    ggplot() +
    annotate("text", x=.5, y=.5, label=txt,
             fontface="bold", size=size,
             colour=LANCET$black) +
    theme_void() +
    theme(plot.margin=margin(0,0,1,0))

  # Uklad: between = 2xN (6 kolumn), within = 1xN (3 kolumny)

  hdr <- function(txt, size=3.0)
    ggplot() +
    annotate("text", x=.5, y=.5, label=txt,
             fontface="bold", size=size,
             colour=LANCET$black) +
    theme_void() +
    theme(plot.margin=margin(0,0,1,0))

  if (comp_name == "between") {
    hdrs <- list(
      hdr("Effect curve"), hdr("Marginal"), hdr("Contribution"),
      hdr("Effect curve"), hdr("Marginal"), hdr("Contribution")
    )
    all_plots <- hdrs
    n_pairs   <- ceiling(n_vars / 2)

    for (p in seq_len(n_pairs)) {
      i1 <- (p-1)*2 + 1
      i2 <- i1 + 1
      pan1 <- panels[[i1]]
      if (i2 <= n_vars) {
        pan2 <- panels[[i2]]
        all_plots <- c(all_plots,
          list(pan1$curve, pan1$marginal, pan1$contrib,
               pan2$curve, pan2$marginal, pan2$contrib))
      } else {
        all_plots <- c(all_plots,
          list(pan1$curve, pan1$marginal, pan1$contrib,
               ggplot()+theme_void(),
               ggplot()+theme_void(),
               ggplot()+theme_void()))
      }
    }

    p_grid <- wrap_plots(
      all_plots,
      ncol    = 6,
      byrow   = TRUE,
      widths  = c(1.2, 0.9, 0.9, 1.2, 0.9, 0.9),
      heights = c(0.12, rep(1, n_pairs))
    )
    fig_width  <- 14.0
    fig_height <- 0.5 + n_pairs * 3.5

  } else {
    # Within: 1 kolumna zmiennych, 3 panele per zmienna
    hdrs <- list(
      hdr("Effect curve"), hdr("Marginal effect"), hdr("Avg. contribution")
    )
    all_plots <- hdrs
    for (pan in panels) {
      all_plots <- c(all_plots,
        list(pan$curve, pan$marginal, pan$contrib))
    }

    p_grid <- wrap_plots(
      all_plots,
      ncol    = 3,
      byrow   = TRUE,
      widths  = c(1.2, 0.9, 0.9),
      heights = c(0.12, rep(1, n_vars))
    )
    fig_width  <- 7.09
    fig_height <- 0.5 + n_vars * 3.5
  }

  # Jedna wspólna, ciągła legenda pod całą figurą. Tylko pierwszy
  # kartogram ma aktywny guide, a patchwork przenosi go pod kompozycję.
  p_grid <- p_grid & theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    plot.margin = margin(4, 6, 12, 4),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

  save_hp_plot(file.path(FIG_DIR, filename), p_grid,
               width  = fig_width,
               height = fig_height,
               dpi    = 300)
  cat("Zapisano:", filename,
      sprintf("(%.1f x %.1f in)\n", fig_width, fig_height))
  invisible(p_grid)
}


# ============================================================
# 4. FOREST PLOT — Linear Mundlak beta vs GAM APS
# ============================================================

make_forest_fig <- function(res_main, res_comp = NULL,
                            filename,
                            title_label    = NULL,
                            add_annotation = FALSE) {

  # ----------------------------------------------------------
  # Linear Mundlak: debiased beta
  # ----------------------------------------------------------

  mdf <- res_main$coef_final %>%
    dplyr::filter(!is.na(beta_debiased)) %>%
    dplyr::mutate(
      estimator = "Linear Mundlak \u03b2",
      effect    = beta_debiased,
      ci_low    = conf.low,
      ci_high   = conf.high,
      layer     = component,
      interval_type = "95% CI"
    ) %>%
    dplyr::select(
      variable, var_label, layer,
      estimator, effect, ci_low, ci_high, interval_type
    )

  # ----------------------------------------------------------
  # Spatial Mundlak GAM: average partial slope
  # ----------------------------------------------------------

  adf <- res_main$aps_all %>%
    dplyr::mutate(
      estimator = "Spatial Mundlak GAM partial slope",
      effect    = aps,
      ci_low    = aps_q25,
      ci_high   = aps_q75,
      layer     = component,
      interval_type = "IQR across counties"
    ) %>%
    dplyr::select(
      variable, var_label, layer,
      estimator, effect, ci_low, ci_high, interval_type
    )

  # ----------------------------------------------------------
  # Optional comparison period
  # ----------------------------------------------------------

  if (!is.null(res_comp)) {
    adf_comp <- res_comp$aps_all %>%
      dplyr::mutate(
        estimator = paste0(
          "Spatial Mundlak GAM partial slope (",
          res_comp$year_label,
          ")"
        ),
        effect    = ame,
        ci_low    = ame_q25,
        ci_high   = ame_q75,
        layer     = component,
        interval_type = "IQR across counties"
      ) %>%
      dplyr::select(
        variable, var_label, layer,
        estimator, effect, ci_low, ci_high, interval_type
      )

    adf <- bind_rows(adf, adf_comp)
  }

  # ----------------------------------------------------------
  # Keep only terms present in both Mundlak and GAM
  # ----------------------------------------------------------

  keys_m <- mdf %>%
    dplyr::distinct(variable, layer) %>%
    dplyr::mutate(k = paste0(variable, "_", layer))

  keys_g <- adf %>%
    dplyr::distinct(variable, layer) %>%
    dplyr::mutate(k = paste0(variable, "_", layer))

  keys_common <- intersect(keys_m$k, keys_g$k)

  forest_df <- bind_rows(
    mdf %>%
      dplyr::filter(paste0(variable, "_", layer) %in% keys_common),
    adf %>%
      dplyr::filter(paste0(variable, "_", layer) %in% keys_common)
  ) %>%
    dplyr::mutate(
      estimator = factor(estimator, levels = unique(estimator)),
      interval_type = factor(
        interval_type,
        levels = c("95% CI", "IQR across counties")
      )
    )

  # ----------------------------------------------------------
  # Variable order: by Linear Mundlak beta
  # ----------------------------------------------------------

  # Sortuj po GAM APS between
  var_order_b <- forest_df %>%
    dplyr::filter(
      layer == "between",
      estimator == "Spatial Mundlak GAM partial slope"
    ) %>%
    dplyr::arrange(effect) %>%
    dplyr::pull(var_label) %>%
    unique()

  var_order_b <- c(
    var_order_b,
    setdiff(
      unique(forest_df$var_label[forest_df$layer == "between"]),
      var_order_b
    )
  )

  # Oba panele muszą mieć dokładnie ten sam zestaw i kolejność zmiennych.
  # Jeżeli dana zmienna nie występuje w warstwie within, pozostaje na osi Y,
  # ale bez punktu i przedziału oszacowania.
  var_order_all <- c(
    var_order_b,
    setdiff(unique(forest_df$var_label), var_order_b)
  )
  var_order_w <- var_order_all
  var_order_b <- var_order_all

  # ----------------------------------------------------------
  # Colours and shapes
  # ----------------------------------------------------------

  est_levels <- levels(forest_df$estimator)

  est_cols <- c(
    "Linear Mundlak \u03b2" = LANCET$navy,
    "Spatial Mundlak GAM partial slope" = LANCET$red,
    setNames(
      c(HAP_COLS$green, HAP_COLS$orange, HAP_COLS$purple)[
        seq_len(max(0, length(est_levels) - 2))
      ],
      est_levels[
        seq_len(max(0, length(est_levels) - 2)) + 2
      ]
    )
  )

  est_cols <- est_cols[est_levels]

  est_shps <- c(16, 17, 15, 18)[seq_along(est_levels)]
  names(est_shps) <- est_levels

  # ----------------------------------------------------------
  # X-axis limits
  # ----------------------------------------------------------

  xlim_b <- forest_df %>%
    dplyr::filter(layer == "between") %>%
    dplyr::summarise(
      lo = min(ci_low, na.rm = TRUE) * 1.08,
      hi = max(ci_high, na.rm = TRUE) * 1.08
    )

  xlim_w <- forest_df %>%
    dplyr::filter(layer == "within") %>%
    dplyr::summarise(
      lo = min(ci_low, na.rm = TRUE) * 1.08,
      hi = max(ci_high, na.rm = TRUE) * 1.08
    )

  # ----------------------------------------------------------
  # Panel builder
  # ----------------------------------------------------------

  make_panel <- function(layer_name, var_order, xlim) {
    dat <- forest_df %>%
      dplyr::filter(layer == layer_name) %>%
      dplyr::mutate(
        var_label = factor(var_label, levels = var_order)
      )

    panel_title <- if (layer_name == "between") {
      "Between-county means"
    } else {
      "Within-county change (lag 1)"
    }

    ggplot(
      dat,
      aes(
        x = effect,
        y = var_label,
        colour = estimator,
        shape = estimator
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = LANCET$grey_mid,
        linewidth = 0.35
      ) +
      geom_linerange(
        aes(
          xmin = ci_low,
          xmax = ci_high
        ),
        linewidth = 0.85,
        position = position_dodge(width = 0.6)
      ) +
      geom_point(
        size = 2.4,
        position = position_dodge(width = 0.6)
      ) +
      scale_colour_manual(
        values = est_cols,
        name = NULL
      ) +
      scale_shape_manual(
        values = est_shps,
        name = NULL
      ) +
      scale_y_discrete(
        limits = var_order,
        drop = FALSE
      ) +
      coord_cartesian(
        xlim = c(xlim$lo, xlim$hi)
      ) +
      scale_x_continuous(
        labels = number_format(accuracy = 1)
      ) +
      labs(
        title = NULL,
        subtitle = NULL,
        caption = NULL,
        x = "\u0394 ASMR per 100k at +1 SD",
        y = NULL
      ) +
      theme_lancet(base_size = 9) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 7.5, colour = LANCET$grey_mid),
        plot.caption = element_blank()
      )
  }

  # ----------------------------------------------------------
  # Compose equal-width panels
  # ----------------------------------------------------------

  p_b <- make_panel(
    "between",
    var_order_b,
    xlim_b
  )

  p_w <- make_panel(
    "within",
    var_order_w,
    xlim_w
  )

  p_out <- (p_b | p_w) +
    plot_layout(
      widths = c(1, 1),
      guides = "collect"
    )

  if (add_annotation) {
    p_out <- p_out +
    plot_annotation(
      title = title_label %||% paste0(
        "Linear Mundlak \u03b2 and GAM average partial slopes | ",
        res_main$year_label
      ),
      subtitle = paste0(
        "For the linear Mundlak reference model, \u03b2 is equivalent ",
        "to a constant slope for standardised predictors; GAM slopes summarise ",
        "the reference partial surface used for curves and maps."
      ),
      theme = theme(
        plot.title = element_text(
          size = 11,
          face = "bold",
          hjust = 0,
          colour = LANCET$black,
          margin = margin(b = 4)
        ),
        plot.subtitle = element_text(
          size = 8,
          colour = LANCET$grey_mid,
          hjust = 0
        ),
        legend.position = "bottom"
      )
    )
  }

  p_out <- p_out &
    theme(
      legend.position = "bottom",
      plot.caption = element_blank()
    )



  n_vars <- max(
    length(var_order_b),
    length(var_order_w)
  )

  h <- max(
    6,
    n_vars * 0.55 + 2.7
  )

  save_hp_plot(
    file.path(FIG_DIR, filename),
    p_out,
    width = 14.5,
    height = h
  )

  cat("Zapisano:", filename, "\n")

  invisible(p_out)
}

make_dominant_contributor_figure <- function(res,
                                             filename    = "Figure_4.png",
                                             active_only = TRUE,
                                             zero_thresh = 1.0,
                                             add_annotation = FALSE) {
  components <- c("between", "within")

  get_active_vars <- function(comp_name) {
    vars <- res$effect_maps %>%
      dplyr::filter(component == comp_name, map_type == "delta") %>%
      dplyr::pull(variable) %>%
      unique()
    if (!isTRUE(active_only)) return(vars)
    active_vars <- res$curves_all %>%
      dplyr::filter(component == comp_name) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(max_abs = max(abs(effect_median), na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::filter(is.finite(max_abs), max_abs >= zero_thresh) %>%
      dplyr::pull(variable)
    vars[vars %in% active_vars]
  }

  # Stale kolory per zmienna
  fixed_pal <- c(
    "Higher education"       = HAP_COLS$blue,
    "Low education"          = HAP_COLS$vermillion,
    "PM2.5"                  = HAP_COLS$green,
    "Unemployment"           = HAP_COLS$purple,
    "Population density"     = HAP_COLS$orange,
    "Income / wages"         = HAP_COLS$sky,
    "Physicians per 100k"    = HAP_COLS$yellow,
    "Hospital beds per 100k" = "#E41A1C",
    "Net migration"          = "#984EA3",
    "Health expenditure"     = "#FF7F00",
    "Distance to city"       = "#A65628",
    "No contribution"        = LANCET$grey_light
  )

  # Globalna paleta i levels ze wszystkich 4 map naraz
  all_vars_global <- purrr::map_dfr(
    c("between", "within"),
    function(cn) {
      vars_use <- get_active_vars(cn)
      res$effect_maps %>%
        dplyr::filter(component == cn, map_type == "delta",
                      variable %in% vars_use) %>%
        dplyr::mutate(var_label = label_var(variable)) %>%
        dplyr::select(var_label, value)
    }
  )

  # Jedna, identyczna dziedzina skali dla wszystkich czterech paneli.
  # Kolejność opiera się najpierw na zdefiniowanej palecie, a następnie
  # dopisuje ewentualne nowe etykiety występujące w danych.
  vars_present <- all_vars_global %>%
    dplyr::filter(!is.na(var_label), var_label != "No contribution") %>%
    dplyr::distinct(var_label) %>%
    dplyr::pull(var_label)
  global_vars <- c(
    intersect(names(fixed_pal), vars_present),
    setdiff(vars_present, names(fixed_pal))
  )
  global_levels <- unique(c(global_vars, "No contribution"))

  global_pal <- fixed_pal[global_levels]
  names(global_pal) <- global_levels
  missing_col <- is.na(global_pal)
  if (any(missing_col)) {
    global_pal[missing_col] <- scales::hue_pal()(sum(missing_col))
  }

  # Helper: mapa dominant contributor dla jednego znaku
  make_dom_map <- function(comp_name, sign, panel_title, show_legend = FALSE) {
    vars_use <- get_active_vars(comp_name)

    df <- res$effect_maps %>%
      dplyr::filter(component == comp_name,
                    map_type  == "delta",
                    variable  %in% vars_use) %>%
      dplyr::mutate(var_label = label_var(variable))

    if (nrow(df) == 0) return(ggplot() + theme_void())

    if (sign == "positive") {
      df_signed <- df %>% dplyr::filter(value > 0)
    } else {
      df_signed <- df %>% dplyr::filter(value < 0) %>%
        dplyr::mutate(value = abs(value))
    }

    best <- df_signed %>%
      dplyr::group_by(code_pow) %>%
      dplyr::slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()

    all_pows <- res$effect_maps %>%
      dplyr::filter(component == comp_name, map_type == "delta") %>%
      dplyr::select(code_pow) %>%
      dplyr::distinct()

    missing_pows <- dplyr::anti_join(all_pows, best, by = "code_pow") %>%
      dplyr::mutate(var_label = "No contribution", value = 0)

    best_full <- dplyr::bind_rows(
      best %>% dplyr::select(code_pow, var_label, value),
      missing_pows
    ) %>%
      dplyr::mutate(var_label = factor(var_label,
                                        levels = global_levels))

    df_map <- best_full %>%
      dplyr::left_join(powiaty_geom, by = "code_pow") %>%
      sf::st_as_sf()

    country_border <- tryCatch(
      sf::st_union(df_map) %>% sf::st_simplify(dTolerance = 500),
      error = function(e) NULL
    )

    p <- ggplot(df_map) +
      geom_sf(aes(fill = var_label), colour = "white",
              linewidth = 0.05) +
      scale_fill_manual(
        values = global_pal,
        drop   = FALSE,          # zachowaj wszystkie levels = jedna legenda
        name   = "Largest contributor",
        guide  = guide_legend(
            nrow           = 2,
          byrow          = TRUE,
          keywidth       = unit(0.4, "cm"),
          keyheight      = unit(0.3, "cm")
        )
      )

    if (!is.null(country_border)) {
      p <- p + geom_sf(data = country_border, fill = NA,
                       colour = "grey25", linewidth = 0.4,
                       inherit.aes = FALSE)
    }

    p + labs(title = panel_title) +
      theme_health_place_map() +
      theme(
        legend.position = if (show_legend) "bottom" else "none",
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 7.5, face = "bold", hjust = 0),
        plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()
      )
  }

  # 2x2: gorny rzad = positive, dolny = negative
  # Tylko jedna mapa ma legendę — patchwork::guide_area() ją zbiera
  p_b_pos <- make_dom_map("between", "positive",
                           "Between — largest positive contribution",
                           show_legend = TRUE)
  p_w_pos <- make_dom_map("within",  "positive",
                           "Within — largest positive contribution",
                           show_legend = FALSE)
  p_b_neg <- make_dom_map("between", "negative",
                           "Between — largest negative contribution",
                           show_legend = FALSE)
  p_w_neg <- make_dom_map("within",  "negative",
                           "Within — largest negative contribution",
                           show_legend = FALSE)

  # Jedna wspólna legenda na dole — osobny obszar legendy.
  p_out <- (p_b_pos | p_w_pos) /
    (p_b_neg | p_w_neg) /
    patchwork::guide_area() +
    plot_layout(heights = c(1, 1, 0.22), guides = "collect") &
    theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7.5, face = "bold", hjust = 0),
      legend.key.size  = unit(0.4, "cm")
    )



  save_hp_plot(file.path(FIG_DIR, filename), p_out,
               width  = 14.0,
               height = 12.0,
               dpi    = 300)
  cat("Zapisano:", filename, "\n")
  invisible(p_out)
}


make_gi_figure <- function(res, comp_name, filename) {
  dat <- res$gi_all %>% dplyr::filter(component == comp_name)
  if (nrow(dat) == 0) {
    warning("Brak gi_all dla: ", res$period, " / ", comp_name)
    return(invisible(NULL))
  }

  vars_u <- dat %>%
    dplyr::distinct(variable, var_label) %>%
    dplyr::arrange(var_label)

  lim_gi <- max(abs(dat$g_i), na.rm = TRUE)
  if (!is.finite(lim_gi) || lim_gi == 0) lim_gi <- 1

  plots <- purrr::map(seq_len(nrow(vars_u)), function(i) {
    vclean <- vars_u$variable[i]
    vlabel <- vars_u$var_label[i]
    df_v <- dat %>% dplyr::filter(variable == vclean)
    make_map(
      df_v, "g_i",
      title = vlabel,
      legend_name = "gᵢ",
      palette = "RdBu",
      symmetric = TRUE,
      limits = c(-lim_gi, lim_gi),
      show_legend = (i == 1L)
    )
  })

  n_cols <- min(4L, length(plots))
  n_rows <- ceiling(length(plots) / n_cols)
  period_tag <- if (identical(res$period, "pre_covid")) "Pre-COVID" else "Full period"
  component_tag <- if (comp_name == "between") "Between" else "Within"

  p_out <- wrap_plots(plots, ncol = n_cols, guides = "collect") /
    patchwork::guide_area() +
    plot_layout(heights = c(1, 0.14), guides = "collect") +
    plot_annotation(
      subtitle = NULL,
      caption = NULL,
      theme = theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()
      )
    ) &
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )

  save_hp_plot(file.path(FIG_DIR, filename), p_out,
               width = 4 * n_cols,
               height = max(5, 3.6 * n_rows + 1.2),
               dpi = 300)
  cat("Zapisano:", filename, "\n")
  invisible(p_out)
}

# ============================================================
# MARGINAL FIGURE — mapy d(ASMR)/d(v) per zmienna (suplement)
# ============================================================
make_marginal_figure <- function(res, comp_name,
                                  filename,
                                  active_only    = TRUE,
                                  zero_thresh    = 1.0,
                                  add_annotation = FALSE) {

  terms_comp <- if (comp_name == "between")
    res$between_terms else res$within_terms
  col <- if (comp_name == "between") LANCET$navy else LANCET$red

  if (active_only && !is.null(res$curves_all)) {
    active_vars <- res$curves_all %>%
      dplyr::filter(component == comp_name) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(max_abs = max(abs(effect_median), na.rm=TRUE),
                       .groups = "drop") %>%
      dplyr::filter(max_abs > zero_thresh) %>%
      dplyr::pull(variable)
    terms_comp <- terms_comp[clean_var(terms_comp) %in% active_vars]
  }

  n_vars <- length(terms_comp)
  if (n_vars == 0) { message("Brak termow: ", comp_name); return(NULL) }

  map_df_all <- res$effect_maps %>%
    dplyr::filter(component == comp_name, map_type == "deriv",
                  variable %in% clean_var(terms_comp))

  has_neg <- any(map_df_all$value < 0, na.rm = TRUE)
  lim_max <- max(abs(map_df_all$value), na.rm = TRUE)
  if (!is.finite(lim_max) || lim_max == 0) lim_max <- 1
  lims_map <- if (has_neg) c(-lim_max, lim_max) else c(0, lim_max)
  pal_map  <- if (has_neg) "RdBu" else "Reds"

  # Etykieta komponentu musi istnieć przed budową listy wykresów.
  comp_label <- if (comp_name == "between") "Between-county" else "Within-county"

  plots <- lapply(terms_comp, function(var) {
    vclean <- clean_var(var)
    vlabel <- label_var(vclean)

    df_m <- res$effect_maps %>%
      dplyr::filter(variable == vclean,
                    component == comp_name,
                    map_type  == "deriv")

    if (nrow(df_m) == 0) return(ggplot() + theme_void())

    make_map(df_m, "value",
             title       = vlabel,
             legend_name = "d(ASMR)/d(v)",
             palette     = pal_map,
             symmetric   = has_neg,
             limits      = lims_map)
  })

  n_cols <- min(3L, n_vars)
  n_rows <- ceiling(n_vars / n_cols)

  p_grid <- wrap_plots(plots, ncol = n_cols, byrow = TRUE)


  h <- 0.4 + n_rows * 3.2
  save_hp_plot(file.path(FIG_DIR, filename), p_grid,
               width = 7.09, height = h, dpi = 300)
  cat("Zapisano marginal:", filename, "\n")
  invisible(p_grid)
}


# ============================================================
# SEKCJA WYWOLAN — generuje wszystkie figury
# ============================================================

cat("\n=== Wczytywanie wyników ===\n")
res_pre  <- load_results("pre_covid")
res_full <- load_results("full")

cat("\n=== Main figures (pre-COVID) ===\n")

# === make_covid_comparison ===

make_covid_comparison <- function(res_pre, res_full,
                                   filename) {

  # APS z obu okresów
  aps_compare <- bind_rows(
    res_pre$aps_all %>%
      mutate(estimator = paste0("Pre-COVID\n",
                                 res_pre$year_label)),
    res_full$aps_all %>%
      mutate(estimator = paste0("Full period\n",
                                 res_full$year_label))
  ) %>%
    mutate(
      estimator = factor(estimator,
                          levels = unique(estimator)),
      effect    = aps,
      ci_low    = aps_q25,
      ci_high   = aps_q75,
      direction = case_when(
        aps_q75 < 0 ~ "lower ASMR",
        aps_q25 > 0 ~ "higher ASMR",
        TRUE        ~ "uncertain"
      )
    )

  # Kolejność wg pre-COVID between
  var_order <- aps_compare %>%
    filter(component == "between",
           str_detect(estimator, "Pre-COVID")) %>%
    arrange(aps) %>%
    pull(var_label) %>% unique()

  est_cols <- c(
    setNames(LANCET$navy,
             levels(aps_compare$estimator)[1]),
    setNames(LANCET$red,
             levels(aps_compare$estimator)[2])
  )

  p_covid <- aps_compare %>%
    mutate(var_label = factor(var_label,
                               levels=var_order)) %>%
    ggplot(aes(x=aps, y=var_label,
               colour=estimator, shape=estimator)) +
    geom_vline(xintercept=0, linetype="dashed",
               colour=LANCET$grey_mid, linewidth=0.35) +
    geom_linerange(
      aes(xmin=ci_low, xmax=ci_high),
      linewidth=0.85,
      position=position_dodge(width=0.6)
    ) +
    geom_point(
      size=2.5,
      position=position_dodge(width=0.6)
    ) +
    scale_colour_manual(values=est_cols, name=NULL) +
    scale_shape_manual(
      values=c(16,17), name=NULL) +
    scale_y_discrete(limits=var_order, drop=FALSE) +
    facet_wrap(~component, scales="free_x") +
    labs(
      title = "Spatial GAM average partial slopes: pre-COVID vs full period",
      subtitle = "Point = average partial slope; interval = IQR across counties.",
      x = "\u0394 ASMR per 100k at +1 SD",
      y = NULL
    ) +
    theme_lancet(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  n_vars <- length(var_order)
  save_hp_plot(file.path(FIG_DIR, filename), p_covid,
               width=14, height=max(6, n_vars*0.5+3))
  cat("Zapisano:", filename, "\n")
}


# ============================================================
# 7. MAPY OPISOWE ZMIENNYCH
# ============================================================


# === make_descriptive_variable_maps ===

make_descriptive_variable_maps <- function(res, filename_means, filename_changes) {
  dat_desc <- res$dat_all
  if (nrow(dat_desc) == 0) return(invisible(NULL))

  all_vars <- unique(c(MODEL_VARS, ONLY_BETWEEN_VARS))
  all_vars <- intersect(all_vars, names(dat_desc))

  means_df <- dat_desc %>%
    group_by(code_pow) %>%
    summarise(across(all_of(all_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(all_vars), names_to = "variable", values_to = "value") %>%
    mutate(var_label = label_var(variable))

  changes_df <- dat_desc %>%
    arrange(code, year) %>%
    group_by(code, code_pow) %>%
    summarise(
      across(all_of(all_vars), ~ {
        x <- .x[is.finite(.x)]
        if (length(x) < 2) return(NA_real_)
        mean(diff(x), na.rm = TRUE)
      }),
      .groups = "drop"
    ) %>%
    group_by(code_pow) %>%
    summarise(across(all_of(all_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(all_vars), names_to = "variable", values_to = "value") %>%
    mutate(var_label = label_var(variable))

  make_desc_grid <- function(df_long, title, subtitle, filename, legend_prefix) {
    vars_u <- df_long %>% distinct(variable, var_label) %>% arrange(var_label)
    plots <- map(seq_len(nrow(vars_u)), function(i) {
      vv <- vars_u$variable[i]
      lbl <- vars_u$var_label[i]
      df_v <- df_long %>% filter(variable == vv)
      if (all(is.na(df_v$value))) return(ggplot() + theme_void())
      make_map(df_v, "value",
               title = lbl,
               legend_name = legend_prefix,
               palette = if (all(df_v$value >= 0, na.rm = TRUE)) "Reds" else "RdBu",
               symmetric = any(df_v$value < 0, na.rm = TRUE),
               limits = NULL)
    })
    n_cols <- 4
    n_rows <- ceiling(nrow(vars_u) / n_cols)
    p_out <- wrap_plots(plots, ncol = n_cols) +
      plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = theme(
          plot.title = element_text(size = 11, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 8, colour = LANCET$grey_mid, hjust = 0)
        )
      )
    save_hp_plot(file.path(FIG_DIR, filename), p_out,
                 width = 16, height = max(7, 3.8 * n_rows))
    cat("Zapisano:", filename, "\n")
    invisible(p_out)
  }

  make_desc_grid(
    means_df,
    subtitle = "Mean raw value for each county across the analysed years. Scales vary across variables.",
    filename = filename_means,
    legend_prefix = "Mean"
  )
  make_desc_grid(
    changes_df,
    subtitle = "Mean year-to-year raw change within each county. Scales vary across variables.",
    filename = filename_changes,
    legend_prefix = "Annual change"
  )
}

# ============================================================
# 8. MAPY ŁĄCZNEJ PREDYKCJI ASMR/ASDR DLA POWIATÓW
# ============================================================


# === make_predicted_powiat_maps ===

make_predicted_powiat_maps <- function(res_pre, res_full,
                                        filename,
                                        value_col = "predicted_total",
                                        map_title = NULL,
                                        legend_name = NULL) {
  if (nrow(res_pre$pred_powiat) == 0 || nrow(res_full$pred_powiat) == 0) {
    warning(
      "Brak pred_powiat w wynikach GAM. Uruchom ponownie 02_stvc_gam.R ",
      "dla pre_covid i full, a potem 03_figures.R."
    )
    return(invisible(NULL))
  }

  dat <- bind_rows(
    res_pre$pred_powiat  %>% mutate(period_label = paste0("Pre-COVID | ", res_pre$year_label)),
    res_full$pred_powiat %>% mutate(period_label = paste0("Full period | ", res_full$year_label))
  )

  if (!value_col %in% names(dat)) {
    warning("Brak kolumny ", value_col, " w pred_powiat.")
    return(invisible(NULL))
  }

  vals <- dat[[value_col]]
  lims <- range(vals, na.rm = TRUE)
  if (!is.finite(lims[1]) || !is.finite(lims[2]) || lims[1] == lims[2]) {
    lims <- c(0, max(vals, na.rm = TRUE))
    if (!is.finite(lims[2]) || lims[2] == 0) lims <- c(0, 1)
  }

  outcome_nm <- outcome_short(unique(dat$outcome)[1] %||% OUTCOME)
  ttl <- map_title %||% paste0("Predicted total ", outcome_nm, " by county")
  leg <- legend_name %||% paste0("Predicted total\n", outcome_nm, " / 100k")

  p_pre <- make_map(
    dat %>% filter(period == "pre_covid"),
    value_col,
    subtitle    = NULL,
    legend_name = leg,
    palette     = "Reds",
    symmetric   = FALSE,
    limits      = lims
  )

  p_full <- make_map(
    dat %>% filter(period == "full"),
    value_col,
    subtitle    = NULL,
    legend_name = leg,
    palette     = "Reds",
    symmetric   = FALSE,
    limits      = lims
  )

  p_out <- (p_pre | p_full) +
    plot_annotation(
      title = ttl,
      subtitle = paste0(
        "Model-based fitted values aggregated to powiat level. ",
        "Common colour scale across periods."
      ),
      theme = theme(
        plot.title    = element_text(size = 11, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 8, colour = LANCET$grey_mid, hjust = 0)
      )
    )

  save_hp_plot(file.path(FIG_DIR, filename), p_out,
               width = 12, height = 6.5)
  cat("Zapisano:", filename, "\n")
  invisible(p_out)
}


# === make_observed_vs_predicted_powiat ===

make_observed_vs_predicted_powiat <- function(res_pre, res_full, filename) {
  if (nrow(res_pre$pred_powiat) == 0 || nrow(res_full$pred_powiat) == 0) {
    warning("Brak pred_powiat; pomijam observed-vs-predicted.")
    return(invisible(NULL))
  }

  dat <- bind_rows(
    res_pre$pred_powiat  %>% mutate(period_label = paste0("Pre-COVID | ", res_pre$year_label)),
    res_full$pred_powiat %>% mutate(period_label = paste0("Full period | ", res_full$year_label))
  )

  outcome_nm <- outcome_short(unique(dat$outcome)[1] %||% OUTCOME)
  lims <- range(c(dat$observed_total, dat$predicted_total), na.rm = TRUE)
  if (!is.finite(lims[1]) || !is.finite(lims[2]) || lims[1] == lims[2]) lims <- c(0, 1)

  p <- ggplot(dat, aes(x = observed_total, y = predicted_total)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = LANCET$grey_mid, linewidth = 0.35) +
    geom_point(alpha = 0.75, size = 1.8, colour = LANCET$navy) +
    facet_wrap(~period_label) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(
      title = paste0("Observed vs predicted total ", outcome_nm, " by county"),
      subtitle = "Each point is one powiat; dashed line = perfect fit.",
      x = paste0("Observed total ", outcome_nm, " / 100k"),
      y = paste0("Predicted total ", outcome_nm, " / 100k")
    ) +
    theme_lancet(base_size = 10) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"))

  save_hp_plot(file.path(FIG_DIR, filename), p,
               width = 10, height = 5.5)
  cat("Zapisano:", filename, "\n")
  invisible(p)
}

# ============================================================
# 9. GENERUJ WSZYSTKIE WYKRESY
# ============================================================


make_grid_figure(res_pre, "between",
                 "Figure_1.png",
                 active_only    = TRUE,
                 add_annotation = FALSE)

make_grid_figure(res_pre, "within",
                 "Figure_2.png",
                 active_only    = TRUE,
                 add_annotation = FALSE)

make_forest_fig(res_pre,
                filename       = "Figure_3.png",
                add_annotation = FALSE)

make_dominant_contributor_figure(res_pre,
                                  filename       = "Figure_4.png",
                                  active_only    = TRUE,
                                  add_annotation = FALSE)

cat("\n=== Supplementary marginal maps ===\n")


cat("\n=== COVID comparison (suplement) ===\n")
make_covid_comparison(res_pre, res_full,
                      "Figure_S2.png")

cat("\n=== Supplementary (full period) ===\n")
make_grid_figure(res_full, "between",
                 "Figure_S5.png",
                 active_only    = TRUE,
                 add_annotation = TRUE)
make_grid_figure(res_full, "within",
                 "Figure_S6.png",
                 active_only    = TRUE,
                 add_annotation = TRUE)

# === make_covid_comparison ===

cat("\n=== Supplementary: descriptive maps ===\n")
make_descriptive_variable_maps(
  res_pre,
  filename_means   = "Figure_S1a.png",
  filename_changes = "Figure_S1b.png"
)
make_descriptive_variable_maps(
  res_full,
  filename_means   = "Figure_S1c.png",
  filename_changes = "Figure_S1d.png"
)

cat("\n=== Supplementary: SVC surfaces ===\n")
make_gi_figure(res_pre, "between", "Figure_S9.png")
make_gi_figure(res_pre, "within",  "Figure_S10.png")

cat("\n=== Supplementary: predicted maps ===\n")
make_predicted_powiat_maps(
  res_pre, res_full,
  filename = "Figure_S13.png"
)
make_observed_vs_predicted_powiat(
  res_pre, res_full,
  filename = "Figure_S15.png"
)

cat("\nWszystkie figury wygenerowane.\n")
