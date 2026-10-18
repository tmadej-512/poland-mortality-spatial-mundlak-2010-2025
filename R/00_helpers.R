# ============================================================
# 00_helpers.R — wspólne funkcje danych/modeli
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

drop_geom <- function(x) {
  if (inherits(x, "sf")) sf::st_drop_geometry(x) else x
}

to_numeric_safe <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

make_code_pow <- function(code) {
  s <- gsub("[^0-9]", "", as.character(code))
  s <- ifelse(nchar(s) == 6L, paste0("0", s), s)
  ifelse(nchar(s) >= 4L, substring(s, 1L, 4L), s)
}

resolve_outcome <- function(d, preferred = OUTCOME) {
  if (preferred %in% names(d)) return(preferred)
  hit <- intersect(OUTCOME_CANDIDATES, names(d))
  if (length(hit) == 0L) {
    stop("Nie znaleziono outcome. Ustaw OUTCOME w 00_config.R. Dostępne kolumny: ",
         paste(names(d), collapse = ", "))
  }
  message("OUTCOME='", preferred, "' nie istnieje; używam: ", hit[1])
  hit[1]
}

normalize_panel_columns <- function(d) {
  # Kanoniczne nazwy wewnętrzne pipeline'u pozostają: code/name/year/code_pow/x/y.
  # Źródłowy panel może używać nazw GUS/TERYT, np. JPT_KOD_JE, powiat, rok.
  if (inherits(d, "sf")) {
    d0 <- d
  } else {
    d0 <- tibble::as_tibble(d)
  }

  find_col <- function(target, candidates, required = TRUE) {
    candidates <- unique(c(target, candidates))

    # Najpierw exact match.
    hit <- candidates[candidates %in% names(d0)]
    if (length(hit) > 0L) return(hit[1])

    # Potem case-insensitive match (bez zgadywania semantyki).
    nms_low <- tolower(names(d0))
    for (cand in candidates) {
      idx <- which(nms_low == tolower(cand))
      if (length(idx)) return(names(d0)[idx[1]])
    }

    if (required) {
      stop(
        "Brak kolumny kanonicznej '", target, "'. Szukano aliasów: ",
        paste(candidates, collapse = ", "),
        "\nDostępne kolumny: ", paste(names(d0), collapse = ", ")
      )
    }
    NULL
  }

  # Rzeczywiste nazwy w finalny_panel_przestrzenny_379_powiatow.rds:
  # JPT_KOD_JE, powiat, rok, centroid_x_2180, centroid_y_2180.
  id_src <- find_col(
    ID_COL,
    c("JPT_KOD_JE", "code", "kod", "id_teryt", "teryt", "TERYT", "jpt_kod_je")
  )
  year_src <- find_col(
    YEAR_COL,
    c("rok", "year", "Rok", "YEAR")
  )
  name_src <- find_col(
    NAME_COL,
    c("powiat", "name", "nazwa", "Nazwa", "JPT_NAZWA_", "jpt_nazwa_"),
    required = FALSE
  )
  pow_src <- find_col(
    POW_CODE_COL,
    c("JPT_KOD_JE", "code_pow", "pow_code", "TERYT", "teryt"),
    required = FALSE
  )
  x_src <- find_col(
    "x",
    c("centroid_x_2180", "x_2180", "centroid_x", "X"),
    required = FALSE
  )
  y_src <- find_col(
    "y",
    c("centroid_y_2180", "y_2180", "centroid_y", "Y"),
    required = FALSE
  )

  # KOPIUJEMY do kanonicznych kolumn zamiast zmieniać nazwy źródłowe.
  # To pozwala użyć JPT_KOD_JE jednocześnie jako ID i kodu do geometrii.
  d0[[ID_COL]] <- as.character(d0[[id_src]])
  d0[[YEAR_COL]] <- suppressWarnings(as.integer(d0[[year_src]]))

  if (!is.null(name_src)) {
    d0[[NAME_COL]] <- as.character(d0[[name_src]])
  } else {
    d0[[NAME_COL]] <- d0[[ID_COL]]
  }

  if (!is.null(pow_src)) {
    d0[[POW_CODE_COL]] <- make_code_pow(d0[[pow_src]])
  } else {
    d0[[POW_CODE_COL]] <- make_code_pow(d0[[ID_COL]])
  }
  d0[[POW_CODE_COL]] <- as.character(d0[[POW_CODE_COL]])

  # Panel ma centroidy EPSG:2180, więc GAM może działać bez ponownego
  # wyliczania centroidów z shapefile. Shapefile nadal jest potrzebny do map polygonowych.
  if (!is.null(x_src)) d0[["x"]] <- to_numeric_safe(d0[[x_src]])
  if (!is.null(y_src)) d0[["y"]] <- to_numeric_safe(d0[[y_src]])

  attr(d0, "column_mapping") <- c(
    code = id_src,
    name = if (is.null(name_src)) "<generated from code>" else name_src,
    year = year_src,
    code_pow = if (is.null(pow_src)) "<generated from code>" else pow_src,
    x = if (is.null(x_src)) "<not found>" else x_src,
    y = if (is.null(y_src)) "<not found>" else y_src
  )

  d0
}

validate_raw_panel <- function(d, outcome_col = OUTCOME) {
  req <- c(ID_COL, NAME_COL, YEAR_COL, POW_CODE_COL, outcome_col, ALL_RAW_VARS)
  miss <- setdiff(req, names(d))
  if (length(miss) > 0L) stop("Brak wymaganych kolumn: ", paste(miss, collapse = ", "))

  yrs <- sort(unique(d[[YEAR_COL]]))
  if (min(yrs, na.rm = TRUE) > 2010L || max(yrs, na.rm = TRUE) < 2025L) {
    warning("Zakres lat panelu: ", min(yrs, na.rm = TRUE), "–", max(yrs, na.rm = TRUE),
            ". Oczekiwano co najmniej 2010–2025.")
  }

  n_county <- dplyr::n_distinct(d[[ID_COL]])
  if (n_county != 379L) warning("Liczba powiatów = ", n_county, "; oczekiwano 379.")

  dup <- d %>%
    drop_geom() %>%
    dplyr::count(.data[[ID_COL]], .data[[YEAR_COL]], name = "n") %>%
    dplyr::filter(n > 1L)
  if (nrow(dup) > 0L) stop("Panel ma duplikaty county-year. Pierwsze rekordy:\n",
                            paste(utils::capture.output(print(utils::head(dup, 10))), collapse = "\n"))
  invisible(TRUE)
}

make_missingness_audit <- function(d, outcome_col = OUTCOME) {
  vars <- intersect(c(outcome_col, ALL_RAW_VARS), names(d))
  d %>%
    drop_geom() %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(vars),
                                   ~ sum(!is.finite(to_numeric_safe(.x))),
                                   .names = "{.col}")) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "n_missing") %>%
    dplyr::mutate(pct_missing = 100 * n_missing / nrow(d))
}

prepare_spatial_objects <- function(dat, powiaty_spdf = NULL) {
  d <- dat

  if (inherits(d, "sf")) {
    geom_by_pow <- d %>%
      dplyr::mutate(!!POW_CODE_COL := as.character(.data[[POW_CODE_COL]])) %>%
      dplyr::group_by(.data[[POW_CODE_COL]]) %>%
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
    powiaty_spdf <- geom_by_pow
    d <- d %>% sf::st_drop_geometry()
  }

  if (!is.null(powiaty_spdf)) {
    pow <- sf::st_as_sf(powiaty_spdf)
    code_candidates <- c("JPT_KOD_JE", POW_CODE_COL, "code_pow", "TERYT", "teryt")
    code_col <- intersect(code_candidates, names(pow))[1]
    if (is.na(code_col) || length(code_col) == 0L) {
      stop("Nie znaleziono kodu powiatu w shapefile. Oczekiwano JPT_KOD_JE/code_pow.")
    }
    pow[[code_col]] <- as.character(pow[[code_col]])
    pow[[code_col]] <- ifelse(nchar(pow[[code_col]]) == 3L,
                              paste0("0", pow[[code_col]]), pow[[code_col]])
    pow <- pow %>%
      dplyr::group_by(.data[[code_col]]) %>%
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
    names(pow)[names(pow) == code_col] <- POW_CODE_COL

    # Centroidy w układzie metrycznym. Jeśli CRS jest geograficzny, transformuj do EPSG:2180.
    pow_xy <- pow
    if (isTRUE(sf::st_is_longlat(pow_xy))) {
      pow_xy <- sf::st_transform(pow_xy, 2180)
    }
    cent <- suppressWarnings(sf::st_centroid(pow_xy))
    xy <- sf::st_coordinates(cent)
    xy_tbl <- tibble::tibble(
      code_pow = as.character(pow_xy[[POW_CODE_COL]]),
      x = xy[, 1], y = xy[, 2]
    )

    d <- d %>%
      dplyr::select(-dplyr::any_of(c("x", "y"))) %>%
      dplyr::left_join(xy_tbl, by = setNames("code_pow", POW_CODE_COL))

    return(list(dat = d, powiaty_spdf = pow))
  }

  if (!all(c("x", "y") %in% names(d))) {
    stop("Brak x/y i brak shapefile. Model przestrzenny wymaga centroidów powiatów.")
  }
  list(dat = d, powiaty_spdf = NULL)
}

fit_scaler <- function(x) {
  x <- to_numeric_safe(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(center = 0, scale = 1))
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) s <- 1
  c(center = mean(x), scale = s)
}

apply_scaler <- function(x, pars) {
  x <- to_numeric_safe(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  out[ok] <- (x[ok] - pars[["center"]]) / pars[["scale"]]
  out
}

# Buduje surowe komponenty Mundlaka, ale jeszcze ich nie standaryzuje.
# mean_years kontroluje lata użyte do county means; lag_source_years określa,
# z jakich lat wolno pobierać ekspozycję t-1.
build_mundlak_components <- function(dat_raw,
                                     period_id = PERIOD,
                                     mean_years = NULL,
                                     lag_source_years = NULL) {
  spec <- PERIOD_SPECS[[period_id]]
  period_years <- seq.int(spec$year_min, spec$year_max)
  if (is.null(mean_years)) mean_years <- period_years
  if (is.null(lag_source_years)) lag_source_years <- mean_years

  d <- dat_raw %>%
    dplyr::filter(.data[[YEAR_COL]] %in% period_years) %>%
    make_period_controls() %>%
    dplyr::arrange(.data[[ID_COL]], .data[[YEAR_COL]])

  train_raw <- d %>% dplyr::filter(.data[[YEAR_COL]] %in% mean_years)
  if (nrow(train_raw) == 0L) stop("Brak danych do county means dla ", period_id)

  means <- train_raw %>%
    drop_geom() %>%
    dplyr::group_by(.data[[ID_COL]]) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(ALL_RAW_VARS),
        ~ if (all(!is.finite(to_numeric_safe(.x)))) NA_real_ else mean(to_numeric_safe(.x), na.rm = TRUE),
        .names = "{.col}_between"
      ),
      .groups = "drop"
    )

  d <- d %>%
    dplyr::left_join(means, by = ID_COL) %>%
    dplyr::group_by(.data[[ID_COL]]) %>%
    dplyr::arrange(.data[[YEAR_COL]], .by_group = TRUE)

  for (v in MODEL_VARS) {
    bnm <- paste0(v, "_between")
    wnm <- paste0(v, "_within")
    lnm <- paste0(v, "_within_lag", WITHIN_LAG_N)
    d <- d %>%
      dplyr::mutate(!!wnm := to_numeric_safe(.data[[v]]) - to_numeric_safe(.data[[bnm]])) %>%
      dplyr::mutate(
        !!lnm := {
          prev_year <- dplyr::lag(.data[[YEAR_COL]], n = WITHIN_LAG_N)
          prev_val  <- dplyr::lag(.data[[wnm]], n = WITHIN_LAG_N)
          ok_lag <- is.finite(prev_year) &
            (.data[[YEAR_COL]] - prev_year == WITHIN_LAG_N) &
            (prev_year %in% lag_source_years)
          ifelse(ok_lag, prev_val, NA_real_)
        }
      )
  }
  d <- d %>% dplyr::ungroup()

  list(
    data = d,
    period_years = period_years,
    mean_years = mean_years,
    lag_source_years = lag_source_years,
    between_raw = paste0(ALL_RAW_VARS, "_between"),
    within_raw = paste0(MODEL_VARS, "_within_lag", WITHIN_LAG_N)
  )
}

fit_component_scalers <- function(component_obj, scale_years = NULL) {
  d <- component_obj$data
  if (is.null(scale_years)) scale_years <- component_obj$mean_years
  tr <- d %>% dplyr::filter(.data[[YEAR_COL]] %in% scale_years)
  out <- list()

  # Between: po jednym county mean na powiat, żeby liczba lat nie ważyła średniej.
  for (nm in component_obj$between_raw) {
    tmp <- tr %>% drop_geom() %>%
      dplyr::select(dplyr::all_of(c(ID_COL, nm))) %>%
      dplyr::distinct()
    out[[nm]] <- fit_scaler(tmp[[nm]])
  }
  # Within: rozkład poprawnych lagowanych odchyleń w obserwacjach treningowych.
  for (nm in component_obj$within_raw) out[[nm]] <- fit_scaler(tr[[nm]])
  out
}

fit_reference_scalers <- function(dat_raw,
                                  reference_period_id = SCALING_REFERENCE_PERIOD,
                                  reference_years = NULL) {
  rspec <- PERIOD_SPECS[[reference_period_id]]
  if (is.null(reference_years)) reference_years <- seq.int(rspec$year_min, rspec$year_max)
  comp <- build_mundlak_components(
    dat_raw,
    period_id = reference_period_id,
    mean_years = reference_years,
    lag_source_years = reference_years
  )
  fit_component_scalers(comp, scale_years = reference_years)
}

scaler_table <- function(scaler_list, source = NA_character_) {
  tibble::tibble(
    raw_term = names(scaler_list),
    center = vapply(scaler_list, function(z) unname(z[["center"]]), numeric(1)),
    scale = vapply(scaler_list, function(z) unname(z[["scale"]]), numeric(1)),
    scaling_source = source
  )
}

build_mundlak_panel <- function(dat_raw,
                                period_id = PERIOD,
                                training_years = NULL,
                                lag_source_years = NULL,
                                scaler_list = NULL,
                                scaling_mode = SCALING_MODE) {
  spec <- PERIOD_SPECS[[period_id]]
  period_years <- seq.int(spec$year_min, spec$year_max)
  if (is.null(training_years)) training_years <- period_years
  if (is.null(lag_source_years)) lag_source_years <- training_years

  # County means są okresowo/foldowo specyficzne; to jest właściwy Mundlak decomposition.
  comp <- build_mundlak_components(
    dat_raw,
    period_id = period_id,
    mean_years = training_years,
    lag_source_years = lag_source_years
  )
  d <- comp$data

  scaler_source <- "supplied"
  if (is.null(scaler_list)) {
    if (identical(scaling_mode, "m1_reference")) {
      scaler_list <- fit_reference_scalers(dat_raw, SCALING_REFERENCE_PERIOD)
      scaler_source <- paste0("reference:", SCALING_REFERENCE_PERIOD)
    } else if (identical(scaling_mode, "training_only")) {
      scaler_list <- fit_component_scalers(comp, scale_years = training_years)
      scaler_source <- paste0("training_only:", paste(range(training_years), collapse = "-"))
    } else {
      stop("Nieznany scaling_mode: ", scaling_mode)
    }
  }

  raw_terms <- c(comp$between_raw, comp$within_raw)
  miss_scalers <- setdiff(raw_terms, names(scaler_list))
  if (length(miss_scalers)) stop("Brak scalerów dla: ", paste(miss_scalers, collapse = ", "))
  for (nm in raw_terms) d[[paste0(nm, "_z")]] <- apply_scaler(d[[nm]], scaler_list[[nm]])

  year_center <- mean(training_years)
  d$year_c <- as.numeric(d[[YEAR_COL]] - year_center)

  between_terms <- paste0(comp$between_raw, "_z")
  within_terms  <- paste0(comp$within_raw, "_z")
  all_terms <- c(between_terms, within_terms)

  list(
    data = d,
    scaler = scaler_list,
    scaler_table = scaler_table(scaler_list, scaler_source),
    scaler_source = scaler_source,
    scaling_mode = scaling_mode,
    year_center = year_center,
    period_years = period_years,
    training_years = training_years,
    lag_source_years = lag_source_years,
    between_terms = between_terms,
    within_terms = within_terms,
    all_terms = all_terms
  )
}

model_complete_cases <- function(d, outcome_col, terms, controls = CONTROL_VARS,
                                 require_xy = TRUE) {
  req <- c(outcome_col, "year_c", terms, controls)
  if (require_xy) req <- c(req, "x", "y")
  missing_req <- setdiff(req, names(d))
  if (length(missing_req)) stop("Brak wymaganych kolumn modelowych: ", paste(missing_req, collapse = ", "))
  chk <- dplyr::select(d, dplyr::all_of(req))
  d[stats::complete.cases(chk), ]
}

smooth1d_term <- function(v, spec = GAM_SPEC) {
  sprintf("s(%s,bs='%s',k=%d)", v, spec$bs_global, as.integer(spec$k_global))
}

svc_term <- function(v, spec = GAM_SPEC) {
  sprintf("s(x,y,by=%s,bs='%s',k=%d)", v, spec$bs_svc, as.integer(spec$k_svc))
}

build_gam_formula <- function(between_terms,
                              within_terms,
                              controls = CONTROL_VARS,
                              spec = GAM_SPEC) {
  arch <- spec$architecture
  allowed <- c("A_global_only", "B_between_svc", "C_all_svc")
  if (!arch %in% allowed) stop("Nieznana architektura GAM: ", arch)

  ctrl_year <- sprintf("s(year_c,bs='%s',k=%d)", spec$bs_year, as.integer(spec$k_year))
  baseline <- sprintf("s(x,y,bs='%s',k=%d)", spec$bs_baseline, as.integer(spec$k_baseline))
  globals <- c(vapply(between_terms, smooth1d_term, character(1), spec = spec),
               vapply(within_terms, smooth1d_term, character(1), spec = spec))

  svc_between <- if (arch %in% c("B_between_svc", "C_all_svc"))
    vapply(between_terms, svc_term, character(1), spec = spec) else character(0)
  svc_within <- if (identical(arch, "C_all_svc"))
    vapply(within_terms, svc_term, character(1), spec = spec) else character(0)

  rhs <- c(ctrl_year, controls, baseline, globals, svc_between, svc_within)
  stats::as.formula(paste(OUTCOME, "~", paste(rhs, collapse = " + ")))
}

make_reference_counties <- function(dat_model, all_terms, controls = CONTROL_VARS) {
  base <- dat_model %>%
    drop_geom() %>%
    dplyr::group_by(.data[[ID_COL]], .data[[POW_CODE_COL]], .data[[NAME_COL]]) %>%
    dplyr::summarise(x = dplyr::first(x), y = dplyr::first(y), .groups = "drop")
  base$year_c <- 0
  for (cc in controls) base[[cc]] <- 0
  for (tt in all_terms) base[[tt]] <- 0
  base
}

safe_range <- function(x, fallback = c(-1, 1)) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(fallback)
  r <- range(x)
  if (r[1] == r[2]) r <- r + c(-0.5, 0.5)
  r
}

strict_fold_data <- function(dat_raw, held_year, period_id = PERIOD) {
  spec <- PERIOD_SPECS[[period_id]]
  yrs <- seq.int(spec$year_min, spec$year_max)
  train_years <- setdiff(yrs, held_year)
  prep <- build_mundlak_panel(
    dat_raw,
    period_id = period_id,
    training_years = train_years,
    lag_source_years = train_years,
    scaling_mode = CV_SCALING_MODE
  )

  train <- prep$data %>%
    dplyr::filter(.data[[YEAR_COL]] %in% train_years) %>%
    model_complete_cases(OUTCOME, prep$all_terms, spec$controls, require_xy = TRUE)
  test <- prep$data %>%
    dplyr::filter(.data[[YEAR_COL]] == held_year) %>%
    model_complete_cases(OUTCOME, prep$all_terms, spec$controls, require_xy = TRUE)

  list(train = train, test = test, prep = prep)
}
